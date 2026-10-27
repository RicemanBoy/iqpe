"""Full iterative QPE (iQPE) on two Steane-code logical qubits, simulated at
circuit level with batched state-vector trajectories (one trajectory per shot).

Logical qubits
  C  iQPE ancilla (control): fresh |+_L> every iteration, measured in X
  S  system register: |1_L>, eigenstate of U = P(2 pi phi); prepared once and
     kept over all iterations (errors on S carry over)

Controlled powers (as in Dustin's `cu`): CP(2a) = P(a)_C P(a)_S CX P(-a)_S CX,
with logical P(a) from {Z, S, T, S_dag, T_dag}, a = k pi/4.  U^(2^m) for
phi = j/8 needs a = j 2^m pi/8, exact in Clifford+T iff j 2^m is even, i.e.
all three iterations are exact for even j (j = 2, 4, 6); odd j would need
sqrt(T) in the last iteration (not implemented here).

Physical conventions (same as phase4_qiskit_hpc/steane_ramsey_qiskit.py):
every operation is a noisy gate, no Pauli frame.  SD6 noise: Pauli error
with total probability p after every 1-qubit gate (uniform X, Y, Z) and CX
(uniform over 15), X flip p after every reset and before every measurement;
basis changes are noisy H gates; corrections are noisy gates; no idle noise.

Building blocks per Steane block (data qubits dq, shared syndrome ancilla A
and flag F): FT |+_L> / |1_L> preparation (gauge round + regular round),
Chao-Reichardt flag QEC after every logical gate, transversal Z/X/H/S/S_dag,
T/T_dag by teleportation (FT magic state |H_L> via flagged H_L measurement +
flagged syndrome round, repeat until accepted, -> |A_L>; transversal CX
data -> magic, transversal Z measurement of the magic block, S_L^m),
transversal CX between blocks + QEC on both, transversal X readout.

Register: data blocks 0..7n-1, A = 7n, F = 7n+1 (n = 1 for Ramsey, 2 for
full iQPE) -> 2^9 or 2^16 amplitudes per shot; the magic block is prepared
in its own 9-qubit register and merged (as qubits 7n..7n+6, ancilla and
flag dropped) only for the teleportation CX and measurement.

  python iqpe_full.py --selftest
  python iqpe_full.py --ftcheck
  python iqpe_full.py --run full --j 2 --p 1e-3 --shots 512
"""

import argparse
import json
import time
from pathlib import Path

import numpy as np

ITERS = 3
SUPPORTS = [(3, 4, 5, 6), (1, 2, 5, 6), (0, 2, 4, 6)]
ORDER = ((3, 4, 5, 6), (1, 2, 6, 5), (0, 6, 2, 4))
ENC = [(2, 4), (2, 5)] + [(piv, t) for piv, sup in ((0, (0, 2, 4, 6)), (1, (1, 2, 5, 6)), (3, (3, 4, 5, 6)))
                          for t in sup if t != piv]
H_OUTCOME = 0
S_PHYS, SDG_PHYS = "sdg", "s"         # transversal S_dag = logical S
TO_A = ("s", "h")                      # transversal gates |H_L> -> |A_L>

s2 = 1 / np.sqrt(2)


def _ry(t):
    return np.array([[np.cos(t / 2), -np.sin(t / 2)], [np.sin(t / 2), np.cos(t / 2)]], dtype=complex)


GATES = {
    "x": np.array([[0, 1], [1, 0]], dtype=complex),
    "y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "z": np.diag([1, -1]).astype(complex),
    "h": np.array([[1, 1], [1, -1]], dtype=complex) * s2,
    "s": np.diag([1, 1j]),
    "sdg": np.diag([1, -1j]),
    "t": np.diag([1, np.exp(1j * np.pi / 4)]),
    "tdg": np.diag([1, np.exp(-1j * np.pi / 4)]),
    "ry+": _ry(np.pi / 4),
    "ry-": _ry(-np.pi / 4),
}
PAULI = [None, GATES["x"], GATES["y"], GATES["z"]]


# ------------------------------------------------------------ fault injection
class Faults:
    """Shot g gets exactly one fault at its loc[g]-th noise location (counted
    per shot), Pauli index pauli[g] (mod 3 for 1-qubit channels; flips fixed)."""

    def __init__(self, loc, pauli):
        self.loc, self.pauli = np.asarray(loc), np.asarray(pauli)
        self.cnt = np.zeros(len(self.loc), int)
        self.kinds = []             # location kinds seen by shot 0 (reference run)


FAULT = None


# ------------------------------------------------------------ simulator
class Sim:
    """Batched noisy state vectors (B, 2^n); ids = global shot index per row."""

    def __init__(self, st, p, rng, ids):
        self.st = st
        self.n = int(round(np.log2(st.shape[1])))
        self.p, self.rng, self.ids = p, rng, np.asarray(ids)

    # ---- exact kernels (rows: None = all, else index array)
    def _get(self, rows):
        return self.st if rows is None else self.st[rows]

    def _put(self, rows, st):
        if rows is not None:
            self.st[rows] = st

    def u1(self, U, q, rows=None):
        st = self._get(rows)
        v = st.reshape(st.shape[0], -1, 2, 1 << q)
        v0, v1 = v[:, :, 0, :], v[:, :, 1, :]
        if U[0, 1] == 0 and U[1, 0] == 0:                 # diagonal
            if U[0, 0] != 1:
                v0 *= U[0, 0]
            if U[1, 1] != 1:
                v1 *= U[1, 1]
        elif U[0, 0] == 0 and U[1, 1] == 0:               # X, Y
            a0 = v0.copy()
            v0[...] = v1
            if U[0, 1] != 1:
                v0 *= U[0, 1]
            v1[...] = a0
            if U[1, 0] != 1:
                v1 *= U[1, 0]
        elif np.isreal(U).all():                           # H, Ry: real 2x2, in place
            (a, b), (c, d) = U.real
            a0 = v0.copy()
            v0 *= a
            v0 += b * v1
            v1 *= d
            v1 += c * a0
        else:
            U = U.astype(st.dtype)
            a0 = v0.copy()
            v0[...] = U[0, 0] * a0 + U[0, 1] * v1
            v1[...] = U[1, 0] * a0 + U[1, 1] * v1
        self._put(rows, st)

    def ucx(self, c, t, rows=None):
        st = self._get(rows)
        hi, lo = max(c, t), min(c, t)
        v = st.reshape(st.shape[0], -1, 2, 1 << (hi - lo - 1), 2, 1 << lo)
        if c > t:
            a = v[:, :, 1, :, 0, :].copy()
            v[:, :, 1, :, 0, :] = v[:, :, 1, :, 1, :]
            v[:, :, 1, :, 1, :] = a
        else:
            a = v[:, :, 0, :, 1, :].copy()
            v[:, :, 0, :, 1, :] = v[:, :, 1, :, 1, :]
            v[:, :, 1, :, 1, :] = a
        self._put(rows, st)

    # ---- noise
    def _hit(self, rows, kind):
        B = self.st.shape[0]
        m = np.ones(B, bool)
        if rows is not None:
            m = np.zeros(B, bool)
            m[rows] = True
        if FAULT is not None:
            g = self.ids[m]
            FAULT.cnt[g] += 1
            if (g == 0).any():
                FAULT.kinds.append(kind)
            m[m] = FAULT.loc[g] == FAULT.cnt[g]
            return np.flatnonzero(m)
        if self.p <= 0:
            return np.zeros(0, int)
        return np.flatnonzero(m & (self.rng.random(B) < self.p))

    def _kind(self, r, n):
        if FAULT is not None:
            return 1 + FAULT.pauli[self.ids[r]] % n
        return self.rng.integers(1, n + 1, r.size)

    def dep1(self, q, rows=None):
        r = self._hit(rows, 1)
        if r.size:
            k = self._kind(r, 3)
            for v in (1, 2, 3):
                if (k == v).any():
                    self.u1(PAULI[v], q, r[k == v])

    def dep2(self, q1, q2, rows=None):
        r = self._hit(rows, 2)
        if r.size:
            k = self._kind(r, 15)
            for v in np.unique(k):
                rr = r[k == v]
                a, b = v // 4, v % 4
                if a:
                    self.u1(PAULI[a], q1, rr)
                if b:
                    self.u1(PAULI[b], q2, rr)

    def flip(self, q, rows=None):
        r = self._hit(rows, "X")
        if r.size:
            self.u1(PAULI[1], q, r)

    # ---- noisy operations
    def g1(self, name, q, rows=None):
        self.u1(GATES[name], q, rows)
        self.dep1(q, rows)

    def cx(self, c, t, rows=None):
        self.ucx(c, t, rows)
        self.dep2(c, t, rows)

    def reset(self, q, rows=None):
        """q is |0> here (fresh or after measure); X flip with probability p."""
        self.flip(q, rows)

    def reset_plus(self, q, rows=None):
        self.reset(q, rows)
        self.g1("h", q, rows)

    def measure(self, q, rows=None):
        """Z measurement (flip before), Born-sampled; qubit left in |0>.
        Rows not measured get -1."""
        B = self.st.shape[0]
        rr = np.arange(B) if rows is None else np.asarray(rows)
        out = np.full(B, -1)
        if rr.size == 0:
            return out
        self.flip(q, rr)
        st = self.st[rr] if rows is not None else self.st
        v = st.reshape(st.shape[0], -1, 2, 1 << q)
        v1 = v[:, :, 1, :]
        p1 = np.einsum("bij,bij->b", v1.real, v1.real) + np.einsum("bij,bij->b", v1.imag, v1.imag)
        p1 = np.clip(p1, 0.0, 1.0)                         # states are kept normalised
        m = (self.rng.random(rr.size) < p1).astype(int)
        one = m == 1
        if one.any():
            v[one, :, 0, :] = v[one, :, 1, :]
        v1[...] = 0
        scale = 1.0 / np.sqrt(np.where(one, p1, 1.0 - p1))
        v[:, :, 0, :] *= scale.astype(st.real.dtype)[:, None, None]
        if rows is not None:
            self.st[rr] = st
        out[rr] = m
        return out

    def measure_x(self, q, rows=None):
        self.g1("h", q, rows)
        return self.measure(q, rows)


# ------------------------------------------------------------ decoding
def syn_tuple(qubits):
    return tuple(int(sum(q in s for q in qubits) % 2) for s in SUPPORTS)


def hook(i):
    return ORDER[i][2:]


def decode_bits(bits):
    syn = np.stack([bits[:, list(s)].sum(1) % 2 for s in SUPPORTS], 1)
    return (bits.sum(1) + syn.any(1)) % 2


class Layout:
    def __init__(self, nblocks):
        self.nb = nblocks
        self.blocks = [list(range(7 * b, 7 * b + 7)) for b in range(nblocks)]
        self.A, self.F = 7 * nblocks, 7 * nblocks + 1
        self.n = 7 * nblocks + 2


DTYPE = np.complex128          # --single: complex64 (half the memory traffic)


def new_state(B, n):
    st = np.zeros((B, 1 << n), dtype=DTYPE)
    st[:, 0] = 1
    return st


# ------------------------------------------------------------ QEC
def measure_stab(sim, L, dq, kind, i, flagged, rows=None):
    A, F = L.A, L.F
    d0, d1, d2, d3 = (dq[q] for q in ORDER[i])
    if kind == "X":
        sim.reset_plus(A, rows)
        if flagged:
            sim.reset(F, rows)
        pairs = [(A, d0)] + ([(A, F)] if flagged else []) + [(A, d1), (A, d2)] + \
                ([(A, F)] if flagged else []) + [(A, d3)]
    else:
        sim.reset(A, rows)
        if flagged:
            sim.reset_plus(F, rows)
        pairs = [(d0, A)] + ([(F, A)] if flagged else []) + [(d1, A), (d2, A)] + \
                ([(F, A)] if flagged else []) + [(d3, A)]
    for c, t in pairs:
        sim.cx(c, t, rows)
    if kind == "X":
        s = sim.measure_x(A, rows)
        f = sim.measure(F, rows) if flagged else None
    else:
        s = sim.measure(A, rows)
        f = sim.measure_x(F, rows) if flagged else None
    return s, f


def syndrome_round(sim, L, dq, flagged, rows=None):
    syn = {"X": [], "Z": []}
    fl = {"X": [], "Z": []}
    for kind in ("X", "Z"):
        for i in range(3):
            s, f = measure_stab(sim, L, dq, kind, i, flagged, rows)
            syn[kind].append(s)
            fl[kind].append(f)
    syn = {k: np.stack(v, 1) for k, v in syn.items()}
    fl = {k: np.stack(v, 1) for k, v in fl.items()} if flagged else None
    return syn, fl


def correct(sim, dq, rows, syn, flags, gate):
    """Noisy correction gates on `rows` (syn, flags: (R, 3) for those rows):
    flagged hook if its syndrome matches (last flag wins), else Hamming."""
    hv = 4 * syn[:, 0] + 2 * syn[:, 1] + syn[:, 2]
    hit = np.zeros((len(rows), 7), bool)
    nz = hv > 0
    hit[nz, hv[nz] - 1] = True
    if flags is not None:
        for i in range(3):
            st_ = syn_tuple(hook(i))
            m = (flags[:, i] == 1) & (hv == 4 * st_[0] + 2 * st_[1] + st_[2])
            hit[m] = False
            for q in hook(i):
                hit[m, q] = True
    for q in range(7):
        r = rows[hit[:, q]]
        if r.size:
            sim.g1(gate, dq[q], r)


def qec(sim, L, dq, gauge=None):
    """Chao-Reichardt round on block dq.  gauge='Z' (after |+>^7): the Z-check
    outcomes are gauge, fixed by X corrections; gauge='X' (after |1>^7): the
    X-check outcomes are gauge, fixed by Z corrections.  Gauge outcomes do
    not trigger the second round."""
    B = sim.st.shape[0]
    allr = np.arange(B)
    syn, fl = syndrome_round(sim, L, dq, True)
    if gauge == "Z":
        correct(sim, dq, allr, syn["Z"], None, "x")
        trig = syn["X"].any(1) | fl["X"].any(1) | fl["Z"].any(1)
    elif gauge == "X":
        correct(sim, dq, allr, syn["X"], None, "z")
        trig = syn["Z"].any(1) | fl["X"].any(1) | fl["Z"].any(1)
    else:
        trig = syn["X"].any(1) | syn["Z"].any(1) | fl["X"].any(1) | fl["Z"].any(1)
    rows = np.flatnonzero(trig)
    if rows.size == 0:
        return
    synB, _ = syndrome_round(sim, L, dq, False, rows)
    correct(sim, dq, rows, synB["X"][rows], fl["Z"][rows], "z")
    correct(sim, dq, rows, synB["Z"][rows], fl["X"][rows], "x")


# ------------------------------------------------------------ logical ops
def prep_plus(sim, L, dq):
    for q in dq:
        sim.reset_plus(q)
    qec(sim, L, dq, gauge="Z")
    qec(sim, L, dq)


def prep_one(sim, L, dq):
    for q in dq:
        sim.reset(q)
        sim.g1("x", q)
    qec(sim, L, dq, gauge="X")
    qec(sim, L, dq)


def transversal(sim, dq, name, rows=None):
    for q in dq:
        sim.g1(name, q, rows)


def magic_prep(B, p, rng, ids):
    """FT |A_L>: (B, 128) states of the 7 magic data qubits."""
    L = Layout(1)
    dq, A, F = L.blocks[0], L.A, L.F
    out = np.zeros((B, 512), dtype=DTYPE)
    todo = np.arange(B)
    while todo.size:
        sim = Sim(new_state(todo.size, 9), p, rng, ids[todo])
        for q in dq:
            sim.reset(q)
        sim.g1("ry+", dq[2])
        for q in (0, 1, 3):
            sim.g1("h", dq[q])
        for c, t in ENC:
            sim.cx(dq[c], dq[t])
        sim.reset_plus(A)
        sim.reset(F)
        for j in range(7):
            sim.g1("ry+", dq[j])            # controlled-H = ry+ . CX . ry- (checked in selftest)
            sim.cx(A, dq[j])
            sim.g1("ry-", dq[j])
            if j in (0, 5):
                sim.cx(A, F)
        ok = (sim.measure_x(A) == H_OUTCOME) & (sim.measure(F) == 0)
        syn, fl = syndrome_round(sim, L, dq, True)
        ok &= ~(syn["X"].any(1) | syn["Z"].any(1) | fl["X"].any(1) | fl["Z"].any(1))
        acc = np.flatnonzero(ok)
        out[todo[acc]] = sim.st[acc]
        todo = todo[~ok]
    sim = Sim(out, p, rng, ids)
    for g in TO_A:
        transversal(sim, dq, g)
    return sim.st.reshape(B, 4, 128)[:, 0, :]


_XOR = np.arange(128)[:, None] ^ np.arange(128)[None, :]      # y ^ z


def gate_t(sim, L, dq, dagger=False):
    """Teleported T (T_dag) on block dq of `sim`, then QEC.

    Transversal CX data -> magic and Z measurement of the magic block are
    done exactly without forming the joint state: with the 7 block bits z of
    a data basis state, outcome y has P(y) = sum_z p_data(z) |M(y^z)|^2 and
    leaves D'(x) = D(x) M(y ^ z(x)).  Noise of the 7 CX (and the measurement
    flips) commutes to the end: the data part acts on the data afterwards,
    an X or Y on a magic qubit flips its recorded bit, a Z there is void."""
    B = sim.st.shape[0]
    M = magic_prep(B, sim.p, sim.rng, sim.ids).astype(sim.st.dtype)
    b0 = dq[0]
    assert list(dq) == list(range(b0, b0 + 7))
    v = sim.st.reshape(B, -1, 128, 1 << b0)               # axis 2 = the 7 block bits
    pz = np.sum(np.abs(v) ** 2, axis=(1, 3))              # (B, 128) data marginal
    py = np.einsum("bz,byz->by", pz, np.abs(M[:, _XOR]) ** 2)
    cum = np.cumsum(py, 1)
    y = np.minimum((sim.rng.random(B)[:, None] * cum[:, -1:] > cum).sum(1), 127)
    v *= M[np.arange(B)[:, None], y[:, None] ^ np.arange(128)[None, :]][:, None, :, None]
    sim.st /= np.sqrt(np.sum(np.abs(sim.st) ** 2, axis=1, keepdims=True))
    flips = np.zeros((B, 7), int)
    data_err = []
    for i in range(7):                                    # noise of the 7 CX
        r = sim._hit(None, 2)
        if r.size:
            k = sim._kind(r, 15)
            a, bm = k // 4, k % 4
            data_err.append((i, r, a))
            flips[r[(bm == 1) | (bm == 2)], i] ^= 1
    for i in range(7):                                    # measurement flips
        r = sim._hit(None, "X")
        flips[r, i] ^= 1
    for i, r, a in data_err:
        for P in (1, 2, 3):
            if (a == P).any():
                sim.u1(PAULI[P], dq[i], r[a == P])
    bits = ((y[:, None] >> np.arange(7)[None, :]) & 1) ^ flips
    m = decode_bits(bits)
    r1 = np.flatnonzero(m == 1)
    if r1.size:
        transversal(sim, dq, S_PHYS, r1)
    if dagger:
        transversal(sim, dq, SDG_PHYS)
    qec(sim, L, dq)


def logical1(sim, L, dq, g):
    if g in ("t", "tdg"):
        gate_t(sim, L, dq, dagger=(g == "tdg"))
        return
    transversal(sim, dq, {"z": "z", "x": "x", "h": "h", "s": S_PHYS, "sdg": SDG_PHYS}[g])
    qec(sim, L, dq)


def logical_cx(sim, L, dc, dt):
    for i in range(7):
        sim.cx(dc[i], dt[i])
    qec(sim, L, dc)
    qec(sim, L, dt)


def readout(sim, dq, basis="X"):
    if basis == "X":
        bits = np.stack([sim.measure_x(q) for q in dq], 1)
    else:
        bits = np.stack([sim.measure(q) for q in dq], 1)
    return decode_bits(bits)


def phase_gates(k):
    """Logical gates for P(k pi/4), k mod 8."""
    k %= 8
    return {0: [], 1: ["t"], 2: ["s"], 3: ["s", "t"], 4: ["z"],
            5: ["sdg", "tdg"], 6: ["sdg"], 7: ["tdg"]}[k]


def cphase(sim, L, dc, ds, k, apply1=None):
    """CP(k pi/2) = P(a)_C P(a)_S CX P(-a)_S CX with a = k pi/4 (only a mod pi
    matters, so k mod 4)."""
    apply1 = apply1 or (lambda dq, g: logical1(sim, L, dq, g))
    k %= 4
    if k == 0:
        return
    for g in phase_gates(k):
        apply1(dc, g)
    for g in phase_gates(k):
        apply1(ds, g)
    logical_cx(sim, L, dc, ds)
    for g in phase_gates(-k):
        apply1(ds, g)
    logical_cx(sim, L, dc, ds)


# ------------------------------------------------------------ protocols
FEEDBACK = {0: "sdg", 1: "tdg"}                       # by age of the earlier bit


def exact_u_power(j, power):
    m = (j * power) % 8
    return (["z"] if m & 4 else []) + (["s"] if m & 2 else []) + (["t"] if m & 1 else [])


def estimate(bits, phi):
    est = (bits * 2.0 ** -(ITERS - np.arange(ITERS))).sum(1)
    d = np.abs(est - phi) % 1
    return np.minimum(d, 1 - d)


def _feedback_groups(bits, t):
    groups = {}
    for b in range(bits.shape[0]):
        g = tuple(FEEDBACK[l] for l, bb in enumerate(bits[b, :t][::-1]) if bb and l in FEEDBACK)
        groups.setdefault(g, []).append(b)
    return groups


def run_full(j, p, shots, seed, check_system=False):
    """Full iQPE (C + S blocks).  Returns (Delta phi per shot, S flipped per shot)."""
    assert (j * 1) % 2 == 0, "exact only for even j"
    rng = np.random.default_rng(seed)
    L = Layout(2)
    dc, ds = L.blocks
    sim = Sim(new_state(shots, L.n), p, rng, np.arange(shots))
    prep_one(sim, L, ds)
    bits = np.zeros((shots, ITERS), int)
    for t in range(ITERS):
        power = 2 ** (ITERS - 1 - t)
        prep_plus(sim, L, dc)
        cphase(sim, L, dc, ds, (j * power // 2) % 4)
        for g, rows in _feedback_groups(bits, t).items():
            rows = np.array(rows)
            sub = Sim(sim.st[rows], p, rng, sim.ids[rows])
            for gg in g:
                logical1(sub, L, dc, gg)
            bits[rows, t] = readout(sub, dc)
            sim.st[rows] = sub.st
    flipped = None
    if check_system:
        flipped = readout(sim, ds, "Z") != 1
    return estimate(bits, j / 8), flipped


SPLIT_STATS = {"splits": 0, "non_product": 0, "max_residual": 0.0}


def _merge(simC, simS, p, rng):
    """Two 9-qubit block registers (ancilla, flag |0>) -> one 16-qubit
    register with C on qubits 0..6, S on 7..13."""
    B = simC.st.shape[0]
    joint = np.zeros((B, 1 << 16), dtype=simC.st.dtype)
    joint[:, :1 << 14] = np.einsum("bs,bc->bsc", simS.st[:, :128], simC.st[:, :128]).reshape(B, 1 << 14)
    return Sim(joint, p, rng, simC.ids)


def _split(joint, simC, simS, tol=1e-4):
    """Factor the 16-qubit state back into C and S (exact when the state is
    a product; checked per shot, the best rank-1 approximation is used and
    counted otherwise)."""
    B = joint.st.shape[0]
    Mx = joint.st[:, :1 << 14].reshape(B, 128, 128).astype(np.complex128)     # [s, c]
    U, sv, Vh = np.linalg.svd(Mx)
    # entanglement measure: weight beyond the leading singular value, relative
    # to it (insensitive to the norm drift of complex64 states)
    res = np.sqrt((sv[:, 1:] ** 2).sum(1)) / sv[:, 0]
    SPLIT_STATS["splits"] += B
    SPLIT_STATS["non_product"] += int((res > tol).sum())
    SPLIT_STATS["max_residual"] = max(SPLIT_STATS["max_residual"], float(res.max()))
    simS.st[:] = 0
    simC.st[:] = 0
    simS.st[:, :128] = (U[:, :, 0] * sv[:, :1]).astype(simS.st.dtype)
    simC.st[:, :128] = Vh[:, 0, :].astype(simC.st.dtype)
    simS.st /= np.sqrt(np.sum(np.abs(simS.st) ** 2, axis=1, keepdims=True))
    simC.st /= np.sqrt(np.sum(np.abs(simC.st) ** 2, axis=1, keepdims=True))


def cphase_split(simC, simS, k):
    """CP(k pi/2) on separate C and S registers; merged only between the CX."""
    k %= 4
    if k == 0:
        return
    L1, L2 = Layout(1), Layout(2)
    d1 = L1.blocks[0]
    dc, ds = L2.blocks
    for g in phase_gates(k):
        logical1(simC, L1, d1, g)
    for g in phase_gates(k):
        logical1(simS, L1, d1, g)
    joint = _merge(simC, simS, simC.p, simC.rng)
    logical_cx(joint, L2, dc, ds)
    for g in phase_gates(-k):
        logical1(joint, L2, ds, g)
    logical_cx(joint, L2, dc, ds)
    _split(joint, simC, simS)


def run_full_split(j, p, shots, seed, check_system=False):
    """Same protocol and same random stream as run_full, but C and S live in
    separate 9-qubit registers except between the two CX of each controlled
    phase (where they are merged into 16 qubits and factored afterwards)."""
    assert j % 2 == 0, "exact only for even j"
    rng = np.random.default_rng(seed)
    L1, L2 = Layout(1), Layout(2)
    d1 = L1.blocks[0]
    dc, ds = L2.blocks
    ids = np.arange(shots)
    simS = Sim(new_state(shots, 9), p, rng, ids)
    simC = Sim(new_state(shots, 9), p, rng, ids)
    prep_one(simS, L1, d1)
    bits = np.zeros((shots, ITERS), int)
    for t in range(ITERS):
        power = 2 ** (ITERS - 1 - t)
        prep_plus(simC, L1, d1)
        cphase_split(simC, simS, (j * power // 2) % 4)
        for g, rows in _feedback_groups(bits, t).items():
            rows = np.array(rows)
            sub = Sim(simC.st[rows], p, rng, ids[rows])
            for gg in g:
                logical1(sub, L1, d1, gg)
            bits[rows, t] = readout(sub, d1)
            simC.st[rows] = sub.st
    flipped = readout(simS, d1, "Z") != 1 if check_system else None
    return estimate(bits, j / 8), flipped


def run_ramsey(j, p, shots, seed):
    """Ramsey iQPE on one block (as phase 4), for cross-checks."""
    rng = np.random.default_rng(seed)
    L = Layout(1)
    dq = L.blocks[0]
    bits = np.zeros((shots, ITERS), int)
    states = new_state(shots, L.n)
    ids = np.arange(shots)
    for t in range(ITERS):
        power = 2 ** (ITERS - 1 - t)
        for g, rows in _feedback_groups(bits, t).items():
            rows = np.array(rows)
            sub = Sim(states[rows], p, rng, ids[rows])
            prep_plus(sub, L, dq)
            for gg in exact_u_power(j, power) + list(g):
                logical1(sub, L, dq, gg)
            bits[rows, t] = readout(sub, dq)
            states[rows] = sub.st
    return estimate(bits, j / 8)


def run_bare_full(j, p, shots, seed):
    """Full iQPE on two physical qubits (ancilla 0, system 1), same noise and
    the same CP decomposition; P(a) is a single noisy phase gate (exact for
    any j)."""
    rng = np.random.default_rng(seed)
    sim = Sim(new_state(shots, 2), p, rng, np.arange(shots))
    sim.reset(1)
    sim.g1("x", 1)
    bits = np.zeros((shots, ITERS), int)
    for t in range(ITERS):
        power = 2 ** (ITERS - 1 - t)
        a = np.pi * j * power / 8
        Pa, Pm = np.diag([1, np.exp(1j * a)]), np.diag([1, np.exp(-1j * a)])
        sim.reset_plus(0)
        if (j * power) % 8:
            for q, U in ((0, Pa), (1, Pa)):
                sim.u1(U, q); sim.dep1(q)
            sim.cx(0, 1)
            sim.u1(Pm, 1); sim.dep1(1)
            sim.cx(0, 1)
        for b in range(t):
            age = t - 1 - b
            if age in FEEDBACK:
                rows = np.flatnonzero(bits[:, b] == 1)
                if rows.size:
                    sim.g1(FEEDBACK[age], 0, rows)
        bits[:, t] = sim.measure_x(0)
    return estimate(bits, j / 8)


# ------------------------------------------------------------ tests
def selftest():
    import itertools
    CX = np.eye(4)[[0, 1, 3, 2]]
    I2 = np.eye(2)
    CH = np.eye(4, dtype=complex)
    CH[2:, 2:] = GATES["h"]
    # circuit order ry+ (target), CX, ry- (target); matrix = ry- CX ry+
    M = np.kron(I2, GATES["ry-"]) @ CX @ np.kron(I2, GATES["ry+"])
    print("CH decomposition:", "ok" if np.allclose(M, CH) else "FAIL")
    rng = np.random.default_rng(1)
    L = Layout(2)
    dc, ds = L.blocks

    def one_iter(gates_c=(), k=None, n=400):
        sim = Sim(new_state(n, L.n), 0.0, rng, np.arange(n))
        prep_one(sim, L, ds)
        prep_plus(sim, L, dc)
        if k is not None:
            cphase(sim, L, dc, ds, k)
        for g in gates_c:
            logical1(sim, L, dc, g)
        c = readout(sim, dc)
        s = readout(sim, ds, "Z")
        return c.mean(), s.mean()

    for label, kw, ideal in (("C: T", dict(gates_c=["t"]), np.sin(np.pi / 8) ** 2),
                             ("CP(pi)   [a=pi/2]", dict(k=2), 1.0),
                             ("CP(pi/2) [a=pi/4]", dict(k=1), 0.5),
                             ("CP(pi/2) then S_dag on C", dict(k=1, gates_c=["sdg"]), 0.0),
                             ("CP(3pi/2) then S on C", dict(k=3, gates_c=["s"]), 0.0)):
        t0 = time.time()
        pc, ps = one_iter(**kw)
        sd = np.sqrt(max(ideal * (1 - ideal), 1 / 400) / 400)
        ok = abs(pc - ideal) <= 4 * sd + 1e-12 and ps == 1.0
        print(f"  {label:28s} P(C=1) = {pc:.3f} (ideal {ideal:.3f}), S in |1_L>: {ps:.3f}  "
              f"{'ok' if ok else 'FAIL'}  [{time.time() - t0:.0f}s]", flush=True)
    for j in (0, 2, 4, 6):
        t0 = time.time()
        d, fl = run_full(j, 0.0, 32, seed=j, check_system=True)
        print(f"  full iQPE j={j}: max Delta phi = {d.max()}, S flipped: {fl.sum()}  "
              f"{'ok' if d.max() == 0 and fl.sum() == 0 else 'FAIL'}  [{time.time() - t0:.0f}s]", flush=True)
    worst = max(run_ramsey(j, 0.0, 32, seed=j).max() for j in range(8))
    print(f"  Ramsey (1 block), all j: max Delta phi = {worst}  {'ok' if worst == 0 else 'FAIL'}")
    for j in range(8):
        d = run_bare_full(j, 0.0, 64, seed=j)
        assert d.max() == 0, j
    print("  bare full iQPE, all j: exact  ok")


FT_SEQS = [("CP(pi)", 2, [], 1), ("CP(pi/2), S_dag", 1, ["sdg"], 0), ("CP(3pi/2), S", 3, ["s"], 0)]


def ftcheck(reps=2, batch=256, seqs=None):
    """Exact single-fault check of one full-iQPE iteration: S prepared in
    |1_L>, C in |+_L>, CP, gates on C, X readout of C and Z readout of S.
    Every location of the (noiseless) reference path gets every Pauli; a
    fault is bad if it changes the C outcome or flips S."""
    global FAULT
    L = Layout(2)
    dc, ds = L.blocks
    rng = np.random.default_rng(5)

    def circuit(loc, pauli):
        global FAULT
        FAULT = Faults(loc, pauli)
        n = len(loc)
        L1 = Layout(1)
        d1 = L1.blocks[0]
        simS = Sim(new_state(n, 9), 0.0, rng, np.arange(n))
        simC = Sim(new_state(n, 9), 0.0, rng, np.arange(n))
        prep_one(simS, L1, d1)
        prep_plus(simC, L1, d1)
        cphase_split(simC, simS, k)
        for g in gates:
            logical1(simC, L1, d1, g)
        c = readout(simC, d1)
        s = readout(simS, d1, "Z")
        kinds, cnt = FAULT.kinds, FAULT.cnt
        FAULT = None
        return c, s, kinds, cnt

    for label, k, gates, ideal in (seqs or FT_SEQS):
        t0 = time.time()
        c, s, kinds, cnt = circuit(np.full(2, -1), np.zeros(2, int))
        assert (c == ideal).all() and (s == 1).all()
        npauli = {1: 3, 2: 15, "X": 1}
        jobs = [(kk, v) for kk, kind in enumerate(kinds, 1) for v in range(npauli[kind]) for _ in range(reps)]
        bad = {}
        for st in range(0, len(jobs), batch):
            chunk = jobs[st:st + batch]
            c, s, _, _ = circuit(np.array([x[0] for x in chunk]), np.array([x[1] for x in chunk]))
            wrong = (c != ideal) | (s != 1)
            for (kk, v), w in zip(chunk, wrong):
                if w:
                    bad[kk] = bad.get(kk, 0) + 1.0 / (npauli[kinds[kk - 1]] * reps)
        print(f"  {label:20s} locations {len(kinds):4d}, faults {len(jobs) // reps:5d}: bad locations {len(bad):3d}, "
              f"c1 = {sum(bad.values()):.3f}" + (f", first bad {sorted(bad)[:10]}" if bad else "")
              + f"  [{time.time() - t0:.0f}s]", flush=True)


# ------------------------------------------------------------ CLI
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--ftcheck", action="store_true")
    ap.add_argument("--run", choices=("full", "ramsey", "bare"))
    ap.add_argument("--j", type=int, default=2)
    ap.add_argument("--p", type=float, default=1e-3)
    ap.add_argument("--shots", type=int, default=256)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default="results")
    ap.add_argument("--single", action="store_true", help="complex64 state vectors")
    a = ap.parse_args()
    global DTYPE
    if a.single:
        DTYPE = np.complex64
    if a.selftest:
        selftest()
        return
    if a.ftcheck:
        ftcheck()
        return
    t0 = time.time()
    fl = None
    if a.run == "full":
        d, fl = run_full(a.j, a.p, a.shots, a.seed, check_system=True)
    elif a.run == "ramsey":
        d = run_ramsey(a.j, a.p, a.shots, a.seed)
    else:
        d = run_bare_full(a.j, a.p, a.shots, a.seed)
    r = {"run": a.run, "j": a.j, "phi": a.j / 8, "p": a.p, "shots": a.shots, "seed": a.seed,
         "sum_dphi": float(d.sum()), "n_wrong": int((d > 0).sum()), "mean_dphi": float(d.mean()),
         "system_flipped": None if fl is None else int(fl.sum()), "seconds": time.time() - t0}
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    json.dump(r, open(out / f"{a.run}_j{a.j}_p{a.p:.1e}_s{a.seed}.json", "w"), indent=1)
    print(f"{a.run} j={a.j} p={a.p:.1e}: mean Delta phi = {r['mean_dphi']:.3e} ({r['n_wrong']} wrong of {a.shots})"
          + ("" if fl is None else f", S flipped {r['system_flipped']}") + f"  [{r['seconds']:.0f}s]")


if __name__ == "__main__":
    main()
