"""
rzsynth -- approximate synthesis of the REALIFIED phase gate

    P(theta) = diag(1, e^{i theta}),     M_P(theta) = Phi(P(theta)),
    Phi(A + iB) = A (x) I + B (x) J,     J = [[0,-1],[1,0]],

over the Toffoli-Hadamard gate set {H, X, Z, CCX}  (+ CX with allow_cx=True).

QUBIT LAYOUT (same index convention as rysynth: x = 4*q2 + 2*q1 + q0)
-----------------------------------------------------------------------
    q0 = data qubit
    q1 = imaginary-part qubit   (|0> carries Re, |1> carries Im)
    q2 = clean ancilla          (enters |0>, returns |0> up to eps)

On the inputs (q2 = 0, basis states 0..3) the target acts as

    |0> -> |0>                      data=0, Re     } identity, reproduced EXACTLY
    |2> -> |2>                      data=0, Im     }
    |1> ->  cos t |1> + sin t |3>   data=1, Re     } rotation by t = theta,
    |3> -> -sin t |1> + cos t |3>   data=1, Im     } i.e. R_y(2 theta) on q1

so M_P is a controlled-R_y(2 theta), control q0, target q1.

ALGORITHM (rysynth adapted to four columns)
-------------------------------------------
Every circuit over the gate set is M / sqrt(2)^k with M an INTEGER matrix and
M M^T = 2^k I.  We need such a circuit whose four input columns are, up to eps,
the four target columns.

STEP 1 (rounding).  The data=0 columns must equal sqrt(2)^k e_x exactly, which
is an integer vector only for EVEN k -- so only even k are scanned (this costs
at most one extra unit of k).  The data=1 columns are rounded as in rysynth:

    u = p e1 + q e3 + r1 e4 + r2 e5 + r3 e6 + r4 e7
    v = -q e1 + p e3 - r2 e4 + r1 e5 - r4 e6 + r3 e7

with p ~ 2^{k/2} cos(theta), q = isqrt(2^k - p^2), and Lagrange's four squares
absorbing R = 2^k - p^2 - q^2 into the four q2=1 rows.  u.v = 0, |u|^2=|v|^2=2^k
and both are orthogonal to e0, e2 identically.

STEP 2 (exact synthesis of the 8x4 isometry).  Left-multiply W = [w0 w1 w2 w3]
by a basis permutation and one H_j to drop k by 1 (possible iff the rows paired
by bit j have equal parity vectors), until k = 0; then a signed permutation
returns the columns to their places.  Permutations come from an exact Dijkstra
table over all 8! basis permutations, minimising (#CCX, #gates); the move at
each level minimises CCX first.  beam=n replaces the greedy choice by a beam
search over reduction paths.  H count is always exactly k.

NOTE: M_P(theta) is reproduced exactly only for theta a multiple of pi/2.
(Even theta = pi/4 needs entries 1 and 1/sqrt2 in the same matrix, which no
Toffoli-Hadamard circuit has.)  Every other angle is approximated to eps.

Gates: ('H',t) ('X',t) ('Z',t) ('CX',c,t) ('CCX',c1,c2,t), listed in TIME ORDER.
"""

import heapq
import math
import random
from itertools import permutations
from math import isqrt

INPUTS = (0, 1, 2, 3)           # tracked input columns (q2 = 0)
DATA, IMAG, ANCILLA = 0, 1, 2
DEFAULT_BEAM = 16

# ======================================================================
# number theory (unchanged from rysynth)
# ======================================================================

_SMALL_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]


def is_prime(n):
    if n < 2:
        return False
    for p in _SMALL_PRIMES:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in _SMALL_PRIMES:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def two_squares_prime(p):
    """p = a^2 + b^2 for a prime p = 1 mod 4 (Cornacchia)."""
    if p == 2:
        return (1, 1)
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    x = pow(z, (p - 1) // 4, p)
    a, b, lim = p, x, isqrt(p)
    while b > lim:
        a, b = b, a % b
    c = isqrt(p - b * b)
    assert b * b + c * c == p
    return (b, c)


def four_squares(R, rng=None):
    """R = a^2+b^2+c^2+d^2 with a,b,c,d >= 0 (Lagrange)."""
    if R < 0:
        raise ValueError("negative input")
    if R == 0:
        return (0, 0, 0, 0)
    rng = rng or random.Random(0xC0FFEE)
    scale = 1
    while R % 4 == 0:
        R //= 4
        scale *= 2
    if R <= 20000:
        for a in range(isqrt(R), -1, -1):
            ra = R - a * a
            for b in range(isqrt(ra), -1, -1):
                rb = ra - b * b
                for c in range(isqrt(rb), -1, -1):
                    d = isqrt(rb - c * c)
                    if d * d == rb - c * c:
                        return tuple(scale * v for v in (a, b, c, d))
        raise RuntimeError("unreachable")
    for _ in range(200000):
        a = rng.randrange(0, isqrt(R) + 1)
        rem = R - a * a
        b = rng.randrange(0, isqrt(rem) + 1)
        S = rem - b * b
        if S == 0:
            return tuple(scale * v for v in (a, b, 0, 0))
        if S == 2:
            return tuple(scale * v for v in (a, b, 1, 1))
        if S % 4 == 1 and is_prime(S):
            c, d = two_squares_prime(S)
            return tuple(scale * v for v in (a, b, c, d))
    raise RuntimeError("four_squares: exhausted trials")


# ======================================================================
# target and rounding
# ======================================================================

def target_columns(theta):
    """The 8x4 float target: columns = images of |0>,|1>,|2>,|3> (q2 = 0)."""
    c, s = math.cos(theta), math.sin(theta)
    T = [[0.0] * 4 for _ in range(8)]
    T[0][0] = 1.0
    T[2][2] = 1.0
    T[1][1], T[3][1] = c, s
    T[1][3], T[3][3] = -s, c
    return T


def _residue_variant(r, idx):
    """idx-th rearrangement of the residue (same error, different circuit)."""
    perms = list(permutations(range(4)))
    p = perms[idx % 24]
    signs = idx // 24
    out = [r[p[i]] for i in range(4)]
    return [(-v if (signs >> i) & 1 else v) for i, v in enumerate(out)]


def find_columns(theta, eps, kmax=400, arrangement=0):
    """Integer 8x4 isometry W (W^T W = 2^k I, k even) within eps of the target.

    Returns (W, k, bound); bound is the exact sup-norm error of the encoded map.
    """
    c, s = math.cos(theta), math.sin(theta)
    ac, as_ = abs(c), abs(s)
    sc = 1 if c >= 0 else -1
    ss = 1 if s >= 0 else -1
    phi = math.atan2(as_, ac)
    for k in range(0, kmax, 2):                       # even k only
        N = 1 << k
        rad = 1 << (k // 2)
        w = max(4, int(2.0 * eps * rad) + 4)
        pc, pmax = int(ac * rad), isqrt(N)
        best = None
        for off in range(0, w + 1):
            for p in ({pc + off, pc - off} if off else {pc}):
                if p < 0 or p > pmax:
                    continue
                q = isqrt(N - p * p)
                # ||u/2^{k/2} - t||^2 written without cancellation
                nu2 = math.ldexp(float(N - p * p - q * q), -k)
                rho = math.sqrt(max(0.0, 1.0 - nu2))
                d = math.atan2(q, p) - phi
                err = math.sqrt(max(0.0, 2.0 * nu2 / (1.0 + rho)
                                    + 4.0 * rho * math.sin(d / 2) ** 2))
                if best is None or err < best[0]:
                    best = (err, p, q)
            if best is not None and best[0] <= eps:
                break
        if best is None or best[0] > eps:
            continue
        err, p, q = best
        r = _residue_variant(four_squares(N - p * p - q * q), arrangement)
        p, q = sc * p, ss * q
        W = [[0] * 4 for _ in range(8)]
        W[0][0] = rad                                  # |0> -> |0>
        W[2][2] = rad                                  # |2> -> |2>
        u = [0, p, 0, q, r[0], r[1], r[2], r[3]]       # image of |1>
        v = [0, -q, 0, p, -r[1], r[0], -r[3], r[2]]    # image of |3>
        for i in range(8):
            W[i][1], W[i][3] = u[i], v[i]
        return W, k, err
    raise RuntimeError("no solution found; increase kmax")


# ======================================================================
# gates
# ======================================================================

def _bit(x, j):
    return (x >> j) & 1


def _perm_image(g, x):
    """Image of basis state x under a permutation gate (X, CX, CCX)."""
    if g[0] == 'X':
        return x ^ (1 << g[1])
    if g[0] == 'CX':
        return x ^ (1 << g[2]) if _bit(x, g[1]) else x
    if g[0] == 'CCX':
        return x ^ (1 << g[3]) if (_bit(x, g[1]) and _bit(x, g[2])) else x
    raise ValueError(g)


def apply_gates(V, k, gates):
    """Apply a time-ordered gate list on the left of V / sqrt(2)^k (V is 8 x m)."""
    V = [list(r) for r in V]
    for g in gates:
        if g[0] == 'H':
            t = g[1]
            new = [None] * 8
            for x in range(8):
                if _bit(x, t):
                    continue
                y = x | (1 << t)
                new[x] = [a + b for a, b in zip(V[x], V[y])]
                new[y] = [a - b for a, b in zip(V[x], V[y])]
            V, k = new, k + 1
        elif g[0] == 'Z':
            V = [[-v for v in V[x]] if _bit(x, g[1]) else V[x] for x in range(8)]
        else:
            new = [None] * 8
            for x in range(8):
                new[_perm_image(g, x)] = V[x]
            V = new
        while k >= 2 and all(v % 2 == 0 for r in V for v in r):
            V = [[v // 2 for v in r] for r in V]
            k -= 2
    return V, k


def circuit_matrix(gates):
    """Exact (integer 8x8 numerator, lde) of a time-ordered circuit."""
    I = [[1 if i == j else 0 for j in range(8)] for i in range(8)]
    return apply_gates(I, 0, gates)


def gate_counts(gates):
    c = {}
    for g in gates:
        c[g[0]] = c.get(g[0], 0) + 1
    c['total'] = len(gates)
    return c


def _cost(gates):
    return (sum(1 for g in gates if g[0] == 'CCX'), len(gates))


# ======================================================================
# optimal permutation circuits (Dijkstra over S_8)
# ======================================================================

class _PermTable:
    """Cheapest circuit, by (#CCX, #gates), for every basis permutation."""

    def __init__(self, allow_cx):
        gens = [('X', t) for t in range(3)]
        if allow_cx:
            gens += [('CX', c, t) for c in range(3) for t in range(3) if c != t]
        gens += [('CCX', 1, 2, 0), ('CCX', 0, 2, 1), ('CCX', 0, 1, 2)]
        imgs = [tuple(_perm_image(g, x) for x in range(8)) for g in gens]
        start = tuple(range(8))
        dist, parent = {start: (0, 0)}, {start: None}
        heap = [(0, 0, start)]
        while heap:
            cc, tt, p = heapq.heappop(heap)
            if dist[p] != (cc, tt):
                continue
            for g, im in zip(gens, imgs):
                q = tuple(im[p[x]] for x in range(8))      # g applied after p
                nc = (cc + (g[0] == 'CCX'), tt + 1)
                if q not in dist or nc < dist[q]:
                    dist[q] = nc
                    parent[q] = (p, g)
                    heapq.heappush(heap, (nc[0], nc[1], q))
        self.dist, self.parent = dist, parent
        # For each perfect matching of the 8 rows: the cheapest (perm, j) that
        # puts every matched pair onto a bit-j pair of positions.
        opts = {}
        for p, c in dist.items():
            inv = [0] * 8
            for x in range(8):
                inv[p[x]] = x
            for j in range(3):
                key = tuple(sorted(tuple(sorted((inv[y], inv[y | (1 << j)])))
                                   for y in range(8) if not _bit(y, j)))
                opts.setdefault(key, []).append((c, j, p))
        self.options = {key: sorted(v)[:12] for key, v in opts.items()}

    def cost(self, p):
        return self.dist[p]

    def gates(self, p):
        out = []
        while self.parent[p] is not None:
            p, g = self.parent[p]
            out.append(g)
        return out[::-1]


_TABLES = {}


def _table(allow_cx):
    if allow_cx not in _TABLES:
        _TABLES[allow_cx] = _PermTable(bool(allow_cx))
    return _TABLES[allow_cx]


# ======================================================================
# reduction moves
# ======================================================================

def _moves(W, table):
    """All (cost, j, perm) such that perm followed by H_j lowers the lde."""
    lab = [tuple(v & 1 for v in r) for r in W]
    out = []
    for key, lst in table.options.items():
        if all(lab[a] == lab[b] for a, b in key):
            out.extend(lst)
    out.sort()
    return out


def _apply_move(W, k, j, p):
    V = [None] * 8
    for x in range(8):
        V[p[x]] = W[x]
    new = [None] * 8
    for x in range(8):
        if _bit(x, j):
            continue
        y = x | (1 << j)
        new[x] = [(a + b) // 2 for a, b in zip(V[x], V[y])]
        new[y] = [(a - b) // 2 for a, b in zip(V[x], V[y])]
    k -= 1
    while k >= 2 and all(v % 2 == 0 for r in new for v in r):
        new = [[v // 2 for v in r] for r in new]
        k -= 2
    return new, k


# ======================================================================
# final signed permutation
# ======================================================================

def _sign_library():
    """(cost, gates, sign-vector over 8 points) for every diagonal +-1 pattern
    of degree <= 2 (constant, Z's, and CZ's built as CCX Z CCX)."""
    lib = {}
    pairs = [(0, 1, 2), (0, 2, 1), (1, 2, 0)]          # (a, b, third qubit)
    for quad in range(8):
        for lin in range(8):
            for const in (0, 1):
                gates, l = [], lin
                for i, (a, b, t) in enumerate(pairs):
                    if (quad >> i) & 1:
                        gates += [('CCX', a, b, t), ('Z', t), ('CCX', a, b, t)]
                        l ^= 1 << t                       # that block adds Z_t
                zs = [t for t in range(3) if (l >> t) & 1]
                if const and zs:                          # -Z_t = X_t Z_t X_t
                    t0 = zs[0]
                    gates += [('X', t0), ('Z', t0), ('X', t0)]
                    gates += [('Z', t) for t in zs[1:]]
                else:
                    gates += [('Z', t) for t in zs]
                    if const:
                        gates += [('X', 0), ('Z', 0), ('X', 0), ('Z', 0)]
                sv = []
                for x in range(8):
                    e = const ^ bin(lin & x).count('1') % 2
                    for i, (a, b, _) in enumerate(pairs):
                        if (quad >> i) & 1:
                            e ^= _bit(x, a) & _bit(x, b)
                    sv.append(e)
                key = tuple(sv)
                c = _cost(gates)
                if key not in lib or c < lib[key][0]:
                    lib[key] = (c, gates)
    return sorted((c, g, k) for k, (c, g) in lib.items())


_SIGNS = _sign_library()


def _cheapest_sign(points, minus):
    for c, g, sv in _SIGNS:
        if all(sv[x] == m for x, m in zip(points, minus)):
            return c, g
    return None


def _tail(V, table):
    """Cheapest signed permutation sending the k=0 columns to e_0..e_3."""
    src, minus = [], []
    for c in range(4):
        nz = [(i, V[i][c]) for i in range(8) if V[i][c] != 0]
        assert len(nz) == 1 and abs(nz[0][1]) == 1, "not a signed basis vector"
        src.append(nz[0][0])
        minus.append(1 if nz[0][1] < 0 else 0)
    rest_src = [x for x in range(8) if x not in src]
    rest_dst = [x for x in range(8) if x not in INPUTS]
    sign_before = _cheapest_sign(src, minus)           # fix signs, then permute
    sign_after = _cheapest_sign(list(INPUTS), minus)   # permute, then fix signs
    best = None
    for order in permutations(rest_dst):
        p = [0] * 8
        for c in range(4):
            p[src[c]] = INPUTS[c]
        for x, y in zip(rest_src, order):
            p[x] = y
        p = tuple(p)
        pc = table.cost(p)
        for sg, before in ((sign_before, True), (sign_after, False)):
            if sg is None:
                continue
            tot = (pc[0] + sg[0][0], pc[1] + sg[0][1])
            if best is None or tot < best[0]:
                best = (tot, p, sg[1], before)
    tot, p, sgates, before = best
    pg = table.gates(p)
    return tot, (sgates + pg if before else pg + sgates)


# ======================================================================
# exact synthesis of the 8x4 isometry
# ======================================================================

def _key(W):
    return tuple(map(tuple, W))


def _beam_width(beam):
    """None/False/0 -> None (greedy); True -> DEFAULT_BEAM; int n >= 1 -> n."""
    if beam is None or beam is False or beam == 0:
        return None
    if beam is True:
        return DEFAULT_BEAM
    if isinstance(beam, int) and beam >= 1:
        return beam
    raise ValueError("beam must be None/0/False (off), True, or an int >= 1; "
                     "got %r" % (beam,))


def _branch_width(branch):
    if isinstance(branch, int) and not isinstance(branch, bool) and branch >= 1:
        return branch
    raise ValueError("branch must be an int >= 1; got %r" % (branch,))


def synthesize_isometry(W, k, allow_cx=False, beam=None, branch=4, finalists=8):
    """Circuit (time order) whose input columns equal W / sqrt(2)^k exactly.

    Returns None if every explored path got stuck (never observed in tests).
    """
    table = _table(allow_cx)
    width = _beam_width(beam)
    branch = 1 if width is None else _branch_width(branch)
    width = width or 1
    nodes = [((0, 0), W, k, None, None)]
    frontier = {k: {_key(W): 0}}
    done = []
    while frontier:
        kk = max(frontier)
        bucket = frontier.pop(kk)
        ids = sorted(bucket.values(), key=lambda i: (nodes[i][0], i))[:width]
        for nid in ids:
            cost, V, kv = nodes[nid][0], nodes[nid][1], nodes[nid][2]
            if kv == 0:
                tc, tg = _tail(V, table)
                done.append(((cost[0] + tc[0], cost[1] + tc[1]), nid, tg))
                continue
            seen, taken = set(), 0
            for mc, j, p in _moves(V, table):
                V2, k2 = _apply_move(V, kv, j, p)
                key2 = _key(V2)
                if key2 in seen:
                    continue
                seen.add(key2)
                c2 = (cost[0] + mc[0], cost[1] + mc[1] + 1)       # +1 for the H
                b = frontier.setdefault(k2, {})
                if key2 not in b or c2 < nodes[b[key2]][0]:
                    nodes.append((c2, V2, k2, nid, (j, p)))
                    b[key2] = len(nodes) - 1
                taken += 1
                if taken >= branch:
                    break
    if not done:
        return None
    best = None
    for _, nid, tg in sorted(done, key=lambda t: (t[0], t[1]))[:finalists]:
        red, i = [], nid
        while nodes[i][3] is not None:
            j, p = nodes[i][4]
            red = table.gates(p) + [('H', j)] + red
            i = nodes[i][3]
        raw = list(reversed(red + tg))          # every gate is an involution
        opt = optimize(raw)
        if best is None or _cost(opt) < best[0]:
            best = (_cost(opt), raw)
    return best[1]


# ======================================================================
# peephole optimisation (exact)
# ======================================================================

_COMM = {}


def _commute(g, h):
    if g == h:
        return True
    key = (g, h)
    if key not in _COMM:
        A, _ = circuit_matrix([g, h])
        B, _ = circuit_matrix([h, g])
        _COMM[key] = _COMM[(h, g)] = (A == B)
    return _COMM[key]


def optimize(gates):
    """Cancel pairs of identical gates separated only by commuting gates."""
    g = list(gates)
    while True:
        n, dead = len(g), [False] * len(g)
        for i in range(n):
            if dead[i]:
                continue
            for j in range(i + 1, n):
                if dead[j]:
                    continue
                if g[j] == g[i]:
                    dead[i] = dead[j] = True
                    break
                if not _commute(g[i], g[j]):
                    break
        if not any(dead):
            return g
        g = [x for x, d in zip(g, dead) if not d]


# ======================================================================
# top level
# ======================================================================

def synthesize_p(theta, eps=1e-6, allow_cx=False, beam=None, branch=4,
                 peephole=True, max_arrangements=48):
    """Approximate M_P(theta) = Phi(diag(1, e^{i theta})) to precision eps.

    theta    : rotation angle (any real number).
    eps      : sup over normalised inputs |psi> (on q0,q1) of
               || C(|psi>|0>) - (M_P|psi>)|0> ||.
    allow_cx : also use CX gates (fewer CCX, sometimes more total gates).
    beam     : None/False/0 -> greedy reduction (fastest).
               int n >= 1   -> beam search keeping n partial reductions per
                               lde level (fewer CCX, ~n/branch x runtime).
               True         -> beam = DEFAULT_BEAM (16).
    branch   : moves expanded per kept state in beam mode (default 4).
    peephole : run the exact gate-cancellation pass on the result.

    Returns a dict with 'gates' (time order), 'counts', 'k', 'W', 'bound',
    'error' (computed from the final circuit), and the options used.
    """
    theta = math.remainder(theta, 2 * math.pi)
    for arr in range(max_arrangements):
        W, k, bound = find_columns(theta, eps, arrangement=arr)
        raw = synthesize_isometry(W, k, allow_cx=allow_cx, beam=beam, branch=branch)
        if raw is not None:
            break
    else:
        raise RuntimeError("reduction stuck for every residue arrangement")
    gates = optimize(raw) if peephole else raw
    width = _beam_width(beam)
    # rysynth-compatible keys: the two ROUNDED columns (data = 1 block),
    # v0 = image of |1> (data=1, Re), v1 = image of |3> (data=1, Im).
    # The data = 0 columns are exactly 2^{k/2} e_0 and 2^{k/2} e_2 (see 'W').
    v0 = [W[i][1] for i in range(8)]
    v1 = [W[i][3] for i in range(8)]
    return {'gates': gates, 'counts': gate_counts(gates), 'k': k, 'W': W,
            'v0': v0, 'v1': v1,
            'bound': bound, 'error': encoded_error(gates, theta),
            'theta': theta, 'eps': eps, 'allow_cx': bool(allow_cx),
            'beam': width, 'branch': branch, 'arrangement': arr}


# ======================================================================
# verification and output
# ======================================================================

def encoded_error(gates, theta):
    """sup_{|psi>} || C(|psi>|0>) - (M_P(theta)|psi>)|0> ||, from the circuit."""
    import numpy as np
    M, k = circuit_matrix(gates)
    r = math.sqrt(2.0) ** k
    A = np.array([[M[i][c] / r for c in INPUTS] for i in range(8)])
    return float(np.linalg.norm(A - np.array(target_columns(theta)), 2))


def realified_matrix(theta):
    """4x4 M_P(theta) on (q0 = data, q1 = imag), index 2*q1 + q0."""
    T = target_columns(theta)
    return [[T[i][c] for c in range(4)] for i in range(4)]


def to_qasm(gates, nqubits=3):
    """OpenQASM 2.0 text; q[0] = data, q[1] = imaginary part, q[2] = ancilla."""
    lines = ['OPENQASM 2.0;', 'include "qelib1.inc";', 'qreg q[%d];' % nqubits]
    for g in gates:
        lines.append('%s %s;' % (g[0].lower(), ','.join('q[%d]' % q for q in g[1:])))
    return '\n'.join(lines) + '\n'


def print_circuit(gates, limit=None):
    out = []
    for i, g in enumerate(gates):
        if limit and i >= limit:
            out.append('... (+%d more)' % (len(gates) - limit))
            break
        out.append('%s(%s)' % (g[0], ','.join('q%d' % q for q in g[1:])))
    return ' '.join(out)


if __name__ == '__main__':
    theta, eps = 0.7, 1e-4
    for opts in ({}, {'allow_cx': True}, {'beam': True}, {'allow_cx': True, 'beam': True}):
        r = synthesize_p(theta, eps, **opts)
        print('%-32s k=%3d  %s  error=%.2e' % (opts or 'default', r['k'],
              {g: r['counts'][g] for g in sorted(r['counts'])}, r['error']))
    print(print_circuit(synthesize_p(theta, 1e-1)['gates']))
