#!/usr/bin/env python3
"""
Rewrite single qubit Clifford+T sequences as Toffoli+Hadamard circuits.

Input  : a text file, one Solovay-Kitaev decomposition per row, gates separated by
         commas, from the alphabet {h, t, tdg, s, sdg, z, x}.
Output : the same file layout, but every row is now a sequence over
         {h0,h1,h2, x0,x1,x2, z0,z1,z2, toff, toff(0,2), toff(1,2), cnot, cnot(c,t)}.
         Read those rows back with parse_output_row() or load_qiskit_circuits() -- a
         plain split(',') would cut token names like toff(0,2) in half.

Why this works
--------------
The Toffoli-Hadamard gate set represents exactly the orthogonal matrices of the form
M / sqrt(2)^k with M integral (Amy, Glaudell, Li, Ross, arXiv:2305.11305).  Realifying a
T gate gives entries 1 and 1/sqrt(2) at the same time, and sqrt(2)^k * 1 is integral only
for even k while sqrt(2)^k / sqrt(2) is integral only for odd k, so realify(T) is not of
that form -- T has no direct Toffoli+H implementation.

The way around it is a catalytic embedding (arXiv:2305.07720).  Let

    Lambda = companion matrix of x^4 + 1  (a 4x4 integral, orthogonal signed 4-cycle),

so Lambda^4 = -I and omega = exp(i pi / 4) is an eigenvalue.  Mapping

    Phi(U) = sum_j U_j  (x)  Lambda^j        for  U = sum_j U_j omega^j

sends every Clifford+T gate to a real matrix of the required M/sqrt(2)^k form, acting on
the data qubit plus a 4 dimensional (2 qubit) catalyst register.  With the catalyst in the
omega eigenvector of Lambda the circuit acts as U on the data and hands the catalyst back
unchanged.  Phi is a ring homomorphism, so a whole word is compiled letter by letter and
concatenated -- no general exact synthesis step is needed here.

Register layout (3 qubits)
--------------------------
    qubit 0 : DATA -- this is the qubit the approximated Rz rotation acts on
    qubit 1 : catalyst high bit b1
    qubit 2 : catalyst low bit  b0          catalyst index j = 2*b1 + b0

Catalyst state to initialise qubits 1,2 with:

    |c_R> = (1/sqrt(2), 1/2, 0, -1/2)  in the basis |b1 b0> = |00>,|01>,|10>,|11>
          = (1/sqrt(2)) |00> + (1/2) |01> - (1/2) |11>

A complex intermediate amplitude psi is carried as Re(psi) (x) |c_R> - Im(psi) (x) |c_I>
with |c_I> = (0, -1/2, -1/sqrt(2), -1/2); measuring the data qubit reproduces |psi_x|^2.
|c_R> is provably not exactly preparable over Toffoli+H (otherwise T would be), so it is
approximated once and reused by every T in every row -- see the module docstring notes in
`catalyst_state`.

Gate naming
-----------
    toff        controls on qubits 0 and 1, target qubit 2   (the unmarked default)
    toff(0,2)   controls on qubits 0 and 2, target qubit 1
    toff(1,2)   controls on qubits 1 and 2, target qubit 0
    cnot        control qubit 0, target qubit 1              (the unmarked default)
    cnot(c,t)   control qubit c, target qubit t
    h1          Hadamard on qubit 1, likewise x<q> and z<q>

Ordering convention
-------------------
A row is read left to right in circuit order: the leftmost gate is the first one applied,
so the row g1,g2,...,gn means U = Gn @ ... @ G2 @ G1.  The output rows use the same
convention as the input rows.  Compilation is order preserving -- each input gate is
replaced in place by its own circuit ordered expansion -- so input and output files are
read the same way.

Usage
-----
    python3 h+toff_decomp.py unitary15.txt
    python3 h+toff_decomp.py unitary15.txt -o unitary15_toffh.txt --angles SK_angles.txt
    python3 h+toff_decomp.py --self-test

and from python, to get the compiled rows back as qiskit circuits:

    import h_toff_decomp as d                    # or importlib, the file name has a '+'
    circuits = d.load_qiskit_circuits('unitary15_toffh.txt')
    qc = d.to_qiskit(d.parse_output_row(row), qc=my_circuit, qubits=(4, 0, 1))
"""

from __future__ import annotations

import argparse
import itertools
import os
import re
import sys

import numpy as np

# =============================================================================
# CATALYST REQUIREMENT -- the compiled circuits are only correct if this holds
# =============================================================================
# The output circuits act on 3 qubits: qubit 0 carries the data, qubits 1 and 2 are the
# catalyst register (catalyst index j = 2*b1 + b0).  Before running any compiled row,
# qubits 1 and 2 MUST be initialised in
#
#     |c_R> = (1/sqrt(2), 1/2, 0, -1/2)              in the basis |b1 b0>
#           = (1/sqrt(2)) |00> + (1/2) |01> - (1/2) |11>
#
# the real part of the omega = exp(i pi/4) eigenvector of LAMBDA, unit normalised.  With
# any other catalyst state the circuits do NOT implement the intended Rz rotations.
#
# One catalyst serves every T of every row.  Each compiled gate hands it back unchanged,
# so it is prepared once for the whole computation, never per gate.  It is provably not
# exactly preparable over Toffoli+H -- if it were, T itself would be exactly
# representable, contradicting the parity argument in the module docstring -- so it is
# the single approximated resource in this pipeline.  See catalyst_state() for the
# preparation recipe (round to sqrt(2)^k, put a four square decomposition of the deficit
# on one extra ancilla, exactly synthesise) and its 2^(-k/4) error scaling.
#
# Complex intermediate amplitudes ride on the catalyst: a data amplitude psi is carried
# as Re(psi) (x) |c_R> - Im(psi) (x) |c_I>, with |c_I> = (0, -1/2, -1/sqrt(2), -1/2).
# Measuring qubit 0 in the computational basis then reproduces |psi_x|^2.
# =============================================================================

# --------------------------------------------------------------------------- constants

N_QUBITS = 3
DIM = 2 ** N_QUBITS

# companion matrix of x^4 + 1 : |j> -> |j+1 mod 4>, with a minus sign on the wrap around
LAMBDA = np.array([[0, 0, 0, -1],
                   [1, 0, 0, 0],
                   [0, 1, 0, 0],
                   [0, 0, 1, 0]], dtype=np.int64)

_I2 = np.eye(2, dtype=np.int64)
_H2 = np.array([[1, 1], [1, -1]], dtype=np.int64)      # Hadamard * sqrt(2)
_X2 = np.array([[0, 1], [1, 0]], dtype=np.int64)
_Z2 = np.array([[1, 0], [0, -1]], dtype=np.int64)


# ------------------------------------------------------- exact arithmetic over Z[1/sqrt2]
# A matrix is stored as a pair (M, k) meaning M / sqrt(2)**k with M integral.  Every
# circuit below stays in this form, so all compilation checks are exact integer identities.

def reduce_pair(M: np.ndarray, k: int):
    """Cancel common factors of 2, keeping the representation canonical."""
    while k >= 2 and np.all(M % 2 == 0):
        M = M // 2
        k -= 2
    return M, k


def mat_mul(a, b):
    """Product of two (M, k) pairs."""
    return reduce_pair(a[0] @ b[0], a[1] + b[1])


def to_float(a) -> np.ndarray:
    return a[0] / np.sqrt(2.0) ** a[1]


# --------------------------------------------------------------------------- gate tables

def _single(q: int, m: np.ndarray) -> np.ndarray:
    """A one qubit integer matrix placed on qubit q of N_QUBITS."""
    ops = [_I2] * N_QUBITS
    ops[q] = m
    out = ops[0]
    for o in ops[1:]:
        out = np.kron(out, o)
    return out


def _perm(rule) -> np.ndarray:
    """Signed permutation from a rule (d,b1,b0) -> (d',b1',b0',sign)."""
    M = np.zeros((DIM, DIM), dtype=np.int64)
    for d, b1, b0 in itertools.product((0, 1), repeat=N_QUBITS):
        nd, n1, n0, s = rule(d, b1, b0)
        M[4 * nd + 2 * n1 + n0, 4 * d + 2 * b1 + b0] = s
    return M


def _cnot_perm(c: int, t: int) -> np.ndarray:
    """CNOT with control qubit c and target qubit t."""
    def rule(*bits):
        b = list(bits)
        b[t] ^= b[c]
        return (b[0], b[1], b[2], 1)
    return _perm(rule)


#: every output token, as an exact (M, k) pair
TOKENS = {
    'toff':      (_perm(lambda d, b1, b0: (d, b1, b0 ^ (d & b1), 1)), 0),   # ctrl 0,1 -> 2
    'toff(0,2)': (_perm(lambda d, b1, b0: (d, b1 ^ (d & b0), b0, 1)), 0),   # ctrl 0,2 -> 1
    'toff(1,2)': (_perm(lambda d, b1, b0: (d ^ (b1 & b0), b1, b0, 1)), 0),  # ctrl 1,2 -> 0
}
for _c, _t in itertools.permutations(range(N_QUBITS), 2):
    TOKENS['cnot(%d,%d)' % (_c, _t)] = (_cnot_perm(_c, _t), 0)
TOKENS['cnot'] = TOKENS['cnot(0,1)']                    # unmarked default
for _q in range(N_QUBITS):
    TOKENS['h%d' % _q] = (_single(_q, _H2), 1)
    TOKENS['x%d' % _q] = (_single(_q, _X2), 0)
    TOKENS['z%d' % _q] = (_single(_q, _Z2), 0)

# Every token is an involution, so inverting a token list is just reversing it.
assert all(np.array_equal(mat_mul(v, v)[0], np.eye(DIM, dtype=np.int64))
           for v in TOKENS.values())


# ----------------------------------------------------------------- compilation of letters
# All lists below are in circuit order: leftmost token applied first.
#
#   Phi(T) applies Lambda to the catalyst controlled on the data qubit, i.e. a controlled
#   increment of the 2 bit catalyst register plus the -1 of the wrap around:
#       CCZ(0,1,2)      the wrap sign, fires on catalyst state |11>
#       toff(0,2)       carry:      b1 ^= b0   (controlled on the data qubit)
#       cnot(0,2)       low bit:    b0 ^= 1    (controlled on the data qubit)
#   Phi(S) = Phi(T)^2 adds 2 instead, so it only touches b1:
#       CZ(0,1)         the wrap sign, fires on b1 = 1
#       cnot(0,1)       b1 ^= 1    (controlled on the data qubit)
#
# CCZ and CZ are a Toffoli / CNOT conjugated by Hadamards on the target.

_CCZ = ['h2', 'toff', 'h2']
_CZ_01 = ['h1', 'cnot(0,1)', 'h1']

COMPILE = {
    'h':   ['h0'],
    'x':   ['x0'],
    'z':   ['z0'],
    't':   _CCZ + ['toff(0,2)', 'cnot(0,2)'],
    's':   _CZ_01 + ['cnot(0,1)'],
}
COMPILE['tdg'] = COMPILE['t'][::-1]     # all tokens are involutions
COMPILE['sdg'] = COMPILE['s'][::-1]

#: accepted spellings in the input files
ALIASES = {'tdag': 'tdg', 't_dg': 'tdg', 'tinv': 'tdg',
           'sdag': 'sdg', 's_dg': 'sdg', 'sinv': 'sdg',
           'i': None, 'id': None}


# ---------------------------------------------------------------- reference embedding Phi

def phi_letter(letter: str):
    """Phi(gate) as an exact (M, k) pair, straight from the algebraic definition."""
    if letter == 'h':
        return (np.kron(_H2, np.eye(4, dtype=np.int64)), 1)
    if letter == 'x':
        return (np.kron(_X2, np.eye(4, dtype=np.int64)), 0)
    powers = {'z': 4, 't': 1, 'tdg': -1, 's': 2, 'sdg': -2}
    if letter not in powers:
        raise KeyError(letter)
    m = powers[letter] % 8
    L = np.linalg.matrix_power(LAMBDA, m).astype(np.int64)
    top = np.zeros((DIM, DIM), dtype=np.int64)
    top[:4, :4] = np.eye(4, dtype=np.int64)
    top[4:, 4:] = L
    return (top, 0)


def unitary_letter(letter: str) -> np.ndarray:
    """The original 2x2 complex gate."""
    w = np.exp(1j * np.pi / 4)
    table = {
        'h': np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2),
        'x': np.array([[0, 1], [1, 0]], dtype=complex),
        'z': np.diag([1, -1]).astype(complex),
        't': np.diag([1, w]), 'tdg': np.diag([1, w.conjugate()]),
        's': np.diag([1, 1j]), 'sdg': np.diag([1, -1j]),
    }
    return table[letter]


# ----------------------------------------------------------------------------- catalyst

def catalyst_state():
    """
    Return (|c_R>, |c_I>), the real and imaginary parts of Lambda's omega eigenvector,
    scaled to unit norm.  Initialise qubits 1,2 in |c_R>.

    |c_R> is not of the form v/sqrt(2)^k, so it cannot be prepared exactly over Toffoli+H
    -- if it could, T itself would be exactly representable.  It is prepared once, to
    precision eps, and then reused: pick k, set v = round(sqrt(2)^k * c_R) on the catalyst
    amplitudes and put a four-square decomposition of D = 2^k - ||v||^2 on one extra
    ancilla (Lagrange guarantees D is a sum of four squares, so the norm matches exactly),
    then exactly synthesise any orthogonal M/sqrt(2)^k with that first column.  The error
    falls as 2^(-k/4).
    """
    ev, evec = np.linalg.eig(LAMBDA.astype(float))
    j = int(np.argmin(np.abs(ev - np.exp(1j * np.pi / 4))))
    c = evec[:, j]
    c = c / np.linalg.norm(c)
    c = c * np.exp(-1j * np.angle(c[0]))          # fix the phase so c[0] is real positive
    return np.sqrt(2) * c.real, np.sqrt(2) * c.imag


def encode(psi: np.ndarray, cR: np.ndarray, cI: np.ndarray) -> np.ndarray:
    """Real 8 dimensional encoding of a complex data state psi."""
    return np.kron(psi.real, cR) - np.kron(psi.imag, cI)


# ------------------------------------------------------------------------------ compiling

#: matches one output token, e.g. h2, z0, toff, toff(0,2), cnot(0,2)
_TOKEN_RE = re.compile(r'[a-z]+\d*(?:\(\s*\d+\s*,\s*\d+\s*\))?')


def parse_output_row(line: str):
    """
    Split one row of an OUTPUT file into tokens.

    Careful: token names such as toff(0,2) and cnot(0,2) contain a comma themselves, so
    line.split(',') does NOT work on the output files -- it would cut 'toff(0,2)' into
    'toff(0' and '2)'.  Use this function (or the same regular expression) when reading
    the generated circuits back in.
    """
    toks = [t.replace(' ', '') for t in _TOKEN_RE.findall(line.strip())]
    unknown = [t for t in toks if t not in TOKENS]
    if unknown:
        raise ValueError('unknown token(s): %s' % ', '.join(unknown))
    return toks


def parse_row(line: str):
    """Split one input row into gate names, tolerating quotes and whitespace."""
    out = []
    for tok in line.strip().split(','):
        g = tok.strip().strip('"').strip("'").lower()
        if not g:
            continue
        g = ALIASES.get(g, g)
        if g is None:
            continue
        if g not in COMPILE:
            raise ValueError('unknown gate %r' % tok)
        out.append(g)
    return out


def compile_row(letters):
    """Clifford+T letters -> Toffoli+H tokens."""
    return [tok for g in letters for tok in COMPILE[g]]


def circuit_matrix(tokens):
    """
    Exact (M, k) matrix of a token list in circuit order: the leftmost token is applied
    first, so it ends up as the rightmost matrix factor.
    """
    acc = (np.eye(DIM, dtype=np.int64), 0)
    for t in tokens:
        acc = mat_mul(TOKENS[t], acc)
    return acc


def phi_word(letters):
    """Exact (M, k) matrix of Phi applied to the whole input word (circuit order)."""
    acc = (np.eye(DIM, dtype=np.int64), 0)
    for g in letters:
        acc = mat_mul(phi_letter(g), acc)
    return acc


def word_unitary(letters) -> np.ndarray:
    """The 2x2 complex unitary of the input word (circuit order)."""
    U = np.eye(2, dtype=complex)
    for g in letters:
        U = unitary_letter(g) @ U
    return U


# -------------------------------------------------------------------------------- qiskit
# Reading the generated files back into qiskit.  QUBIT 0 IS THE DATA QUBIT -- the
# approximated Rz acts on it.  Qubits 1 and 2 are the catalyst register and have to be
# prepared in |c_R> beforehand (see the CATALYST REQUIREMENT block at the top); they must
# stay coherent across the whole computation, since they carry the imaginary part of every
# amplitude, so never reset or measure them between rows.

def token_qargs(tok: str):
    """
    Translate one output token into (gate name, qubit indices).

    Returns ('h'|'x'|'z', (q,)), ('cx', (control, target)) or
    ('ccx', (control, control, target)), with qubit 0 the data qubit.
    """
    if tok not in TOKENS:
        raise ValueError('unknown token %r' % tok)
    if tok.startswith('toff'):
        ctrl = (0, 1) if tok == 'toff' else tuple(int(c) for c in tok[5:-1].split(','))
        target = ({0, 1, 2} - set(ctrl)).pop()
        return 'ccx', (ctrl[0], ctrl[1], target)
    if tok.startswith('cnot'):
        c, t = (0, 1) if tok == 'cnot' else tuple(int(c) for c in tok[5:-1].split(','))
        return 'cx', (c, t)
    return tok[0], (int(tok[1:]),)


def to_qiskit(tokens, qc=None, qubits=None):
    """
    Append one compiled row to a qiskit QuantumCircuit, in circuit order.

    tokens : list of output tokens, e.g. from parse_output_row()
    qc     : circuit to append to; a fresh 3 qubit circuit is made if omitted
    qubits : which physical qubits to use as (data, b1, b0); defaults to (0, 1, 2).
             Pass e.g. qubits=(4, 0, 1) to put the data on qubit 4 while keeping the
             shared catalyst on qubits 0 and 1.

    Uses only qc.h, qc.x, qc.z, qc.cx and qc.ccx.
    """
    from qiskit import QuantumCircuit

    if qc is None:
        qc = QuantumCircuit(N_QUBITS)
    q = list(range(N_QUBITS)) if qubits is None else list(qubits)
    for tok in tokens:
        gate, args = token_qargs(tok)
        if gate == 'h':
            qc.h(q[args[0]])
        elif gate == 'x':
            qc.x(q[args[0]])
        elif gate == 'z':
            qc.z(q[args[0]])
        elif gate == 'cx':
            qc.cx(q[args[0]], q[args[1]])
        elif gate == 'ccx':
            qc.ccx(q[args[0]], q[args[1]], q[args[2]])
        else:                                            # pragma: no cover
            raise ValueError('no qiskit gate for token %r' % tok)
    return qc


def load_qiskit_circuits(path, qubits=None):
    """
    Read a generated output file and return one QuantumCircuit per row.

    The rows come back in file order, so row i is the circuit for the i-th angle of the
    corresponding angle file.
    """
    with open(path) as fh:
        return [to_qiskit(parse_output_row(line), qubits=qubits)
                for line in fh if line.strip()]


def qiskit_matches(tokens) -> bool:
    """
    Check a qiskit circuit against the exact integer matrix of the same token list.

    qiskit orders basis states little endian (qubit 0 is the least significant bit) while
    the matrices here use qubit 0 as the most significant bit, hence the reverse_bits().
    """
    from qiskit.quantum_info import Operator

    got = Operator(to_qiskit(tokens).reverse_bits()).data
    return bool(np.allclose(got, to_float(circuit_matrix(tokens))))


# ------------------------------------------------------------------------------ checking

def phase_invariant_distance(A: np.ndarray, B: np.ndarray) -> float:
    """Operator norm distance minimised over the global phase."""
    t = np.trace(A.conjugate().T @ B)
    ph = np.conjugate(t) / abs(t) if abs(t) > 1e-14 else 1.0
    return float(np.linalg.norm(A - ph * B, 2))


def rz(theta: float) -> np.ndarray:
    return np.diag([np.exp(-0.5j * theta), np.exp(0.5j * theta)])


def check_row(letters, tokens, target_angle=None, rng=None):
    """Verify one compiled row and return a dict of diagnostics."""
    rng = rng or np.random.default_rng(0)
    comp = circuit_matrix(tokens)
    ref = phi_word(letters)
    exact = (comp[1] == ref[1]) and np.array_equal(comp[0], ref[0])

    R = to_float(comp)
    cR, cI = catalyst_state()
    U = word_unitary(letters)
    action = 0.0
    for _ in range(4):
        psi = rng.normal(size=2) + 1j * rng.normal(size=2)
        psi = psi / np.linalg.norm(psi)
        action = max(action, float(np.linalg.norm(
            R @ encode(psi, cR, cI) - encode(U @ psi, cR, cI))))

    out = {'exact': exact, 'orthogonal': bool(np.allclose(R @ R.T, np.eye(DIM))),
           'k': comp[1], 'action_error': action,
           'n_in': len(letters), 'n_out': len(tokens),
           'n_toff': sum(t.startswith('toff') for t in tokens)}
    if target_angle is not None:
        # the adjoint files approximate Rz(-theta), so report whichever sign fits
        e_plus = phase_invariant_distance(U, rz(target_angle))
        e_minus = phase_invariant_distance(U, rz(-target_angle))
        out['sk_error'] = min(e_plus, e_minus)
        out['sign'] = '+' if e_plus <= e_minus else '-'
    return out


def self_test() -> bool:
    """Check every letter identity and the catalyst, in exact integer arithmetic."""
    ok = True
    for g in ('h', 'x', 'z', 't', 'tdg', 's', 'sdg'):
        a, b = circuit_matrix(COMPILE[g]), phi_letter(g)
        good = (a[1] == b[1]) and np.array_equal(a[0], b[0])
        ok &= good
        print('  Phi(%-4s) == %-58s %s'
              % (g, ','.join(COMPILE[g])[:56], 'ok' if good else 'MISMATCH'))
    cR, cI = catalyst_state()

    # A deliberately non palindromic word, to pin the circuit ordering convention: the
    # compiled circuit must implement Gn...G1 (leftmost applied first) and must NOT
    # implement the opposite reading G1...Gn.
    word = ['t', 'h', 's', 'h', 'tdg', 'h']
    R = to_float(circuit_matrix(compile_row(word)))
    psi = np.array([0.6, 0.8j])
    U_circuit = word_unitary(word)
    U_reversed = np.eye(2, dtype=complex)
    for g in word:
        U_reversed = U_reversed @ unitary_letter(g)
    d_circuit = np.linalg.norm(R @ encode(psi, cR, cI) - encode(U_circuit @ psi, cR, cI))
    d_reversed = np.linalg.norm(R @ encode(psi, cR, cI) - encode(U_reversed @ psi, cR, cI))

    checks = [
        ('leftmost gate applied first (%s)' % ','.join(word), d_circuit < 1e-12),
        ('  and not the reverse reading (dist %.2f)' % d_reversed, d_reversed > 1e-3),
        ('catalyst |c_R> = (1/sqrt2, 1/2, 0, -1/2)',
         np.allclose(cR, [1 / np.sqrt(2), 0.5, 0.0, -0.5])),
        ('|c_R>, |c_I> orthonormal',
         abs(cR @ cI) < 1e-12 and abs(np.linalg.norm(cR) - 1) < 1e-12
         and abs(np.linalg.norm(cI) - 1) < 1e-12),
        ('Lambda^4 == -I',
         np.array_equal(np.linalg.matrix_power(LAMBDA, 4), -np.eye(4, dtype=np.int64))),
    ]
    for name, good in checks:
        ok &= good
        print('  %-58s %s' % (name, 'ok' if good else 'FAILED'))

    # the qiskit translation, if qiskit is available
    try:
        good = all(qiskit_matches(compile_row([g]))
                   for g in ('h', 'x', 'z', 't', 'tdg', 's', 'sdg'))
        good &= qiskit_matches(compile_row(word))
        ok &= good
        print('  %-58s %s' % ('qiskit circuit == exact matrix (all letters)',
                              'ok' if good else 'FAILED'))
    except ImportError:
        print('  %-58s %s' % ('qiskit not installed, translation not checked', 'skipped'))
    return bool(ok)


# ---------------------------------------------------------------------------------- main

def convert(in_path, out_path=None, angles_path=None, verbose=True):
    """Compile a whole file and print a verification report.  Returns the report rows."""
    if out_path is None:
        stem, ext = os.path.splitext(in_path)
        out_path = stem + '_toffh' + (ext or '.txt')

    with open(in_path) as fh:
        rows = [parse_row(l) for l in fh if l.strip()]

    angles = None
    if angles_path and os.path.exists(angles_path):
        with open(angles_path) as fh:
            vals = [float(l) for l in fh if l.strip()]
        if len(vals) == len(rows):
            angles = [np.pi * v for v in vals]        # file stores the angle in units of pi
        elif verbose:
            print('  (angle file has %d entries for %d rows -- ignored)'
                  % (len(vals), len(rows)))

    compiled = [compile_row(r) for r in rows]
    with open(out_path, 'w') as fh:
        for toks in compiled:
            fh.write(','.join(toks) + '\n')

    with open(out_path) as fh:
        reread = [parse_output_row(l) for l in fh if l.strip()]
    if reread != compiled:
        raise RuntimeError('the written file does not read back identically')

    report = []
    rng = np.random.default_rng(0)
    for i, (letters, toks) in enumerate(zip(rows, compiled)):
        report.append(check_row(letters, toks,
                                angles[i] if angles else None, rng))

    if verbose:
        print('%s -> %s' % (in_path, out_path))
        head = '  row  in  out  toff   k  exact  action err'
        if angles:
            head += '   target        SK err vs Rz'
        print(head)
        for i, r in enumerate(report):
            line = ('  %3d %4d %4d %5d %3d  %-5s  %.2e'
                    % (i, r['n_in'], r['n_out'], r['n_toff'], r['k'],
                       'yes' if r['exact'] else 'NO', r['action_error']))
            if angles:
                line += ('   Rz(%s%.4f pi)  %.3e'
                         % (r['sign'], angles[i] / np.pi, r['sk_error']))
            print(line)
        print('  file reads back exactly : True (use parse_output_row, not split(\',\'))')
        print('  all rows exact          : %s' % all(r['exact'] for r in report))
        print('  all rows orthogonal     : %s' % all(r['orthogonal'] for r in report))
        print('  max data action error   : %.2e'
              % max(r['action_error'] for r in report))
        print('  gate count %d -> %d (%d Toffoli), factor %.1f'
              % (sum(r['n_in'] for r in report), sum(r['n_out'] for r in report),
                 sum(r['n_toff'] for r in report),
                 sum(r['n_out'] for r in report) / max(1, sum(r['n_in'] for r in report))))
    return report


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n')[1],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('input', nargs='*', help='input txt file(s)')
    p.add_argument('-o', '--output', help='output file (only with a single input)')
    p.add_argument('--angles', help='file with the target angles in units of pi')
    p.add_argument('--self-test', action='store_true', help='check the gate identities')
    a = p.parse_args(argv)

    if a.self_test or not a.input:
        print('self test:')
        ok = self_test()
        print('  -> %s' % ('all identities verified' if ok else 'FAILURES'))
        if not a.input:
            return 0 if ok else 1
        print()

    if a.output and len(a.input) > 1:
        p.error('-o cannot be used with several input files')

    for path in a.input:
        convert(path, a.output, a.angles)
        print()
    return 0


if __name__ == '__main__':
    sys.exit(main())
