"""
rysynth -- approximate synthesis of R_y(theta) over the Toffoli-Hadamard gate
set {H, X, Z, CCX} on 1 data qubit + 2 clean ancillas.

WHY THIS SHAPE OF ALGORITHM
---------------------------
Every generator of {H,X,Z,CCX} is (integer matrix)/sqrt(2)^k, and products
multiply numerators, so every circuit matrix is M/sqrt(2)^k with M integer and
M M^T = 2^k I.  (Equivalently U^sigma = (-1)^k U under sqrt2 -> -sqrt2.)  This
is strictly stronger than "entries in Z[1/sqrt2]", which is why embedding a
Clifford+T circuit via i -> J does NOT land in the group: it produces entries
1 and 1/sqrt2 in the same matrix.

We need a circuit C with  C(|psi>|00>) = (R_y(theta)|psi>)|00> + O(eps), i.e.
only C's first two columns are constrained:

    v0 ~ sqrt(2)^k (c, s, 0,0,0,0,0,0),  v1 ~ sqrt(2)^k (-s, c, 0,...,0)
    c = cos(theta/2), s = sin(theta/2),  v0,v1 in Z^8, |vi|^2 = 2^k, v0.v1 = 0.

STEP 1 (rounding).  Put p ~ sqrt(2)^k c, q = isqrt(2^k - p^2) and let
Lagrange's four-square theorem absorb the residue R = 2^k - p^2 - q^2:

    v0 = (p, q, r1, r2, r3, r4, 0, 0),   v1 = (-q, p, -r2, r1, -r4, r3, 0, 0)

(v1 is v0 with each consecutive pair rotated by i, which makes v0.v1 = 0 and
|v1| = |v0| automatic).  A 1-D scan over p reaches k ~ 3 log2(1/eps).
No grid problem, no factoring, no Diophantine equation over Z[omega].

STEP 2 (exact synthesis of the 8x2 isometry).  Left-multiply V = [v0 v1] by
generators until k = 0.  A Hadamard on qubit j lowers the lde by 1 iff the row
parities of V agree within each bit-j pair.  Row parities live in F_2^2, and
|v0|^2, |v1|^2 even together with v0.v1 = 0 mod 2 force all four parity classes
to have even size -- so a perfect matching always exists and a permutation
turns it into a bit-j matching.  The reduction therefore never gets stuck.
(Reducing a full 8x8 matrix does get stuck: a generic octonion left-multi-
plication matrix has all 8 row parities distinct, forming the [8,4,4] Hamming
code, and then no permutation + H can lower the lde.)
"""

import math
import random
from math import isqrt

# ----------------------------------------------------------------- number theory

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
    x = pow(z, (p - 1) // 4, p)               # x^2 = -1 mod p
    assert (x * x + 1) % p == 0
    a, b, lim = p, x, isqrt(p)
    while b > lim:
        a, b = b, a % b
    c = isqrt(p - b * b)
    assert b * b + c * c == p
    return (b, c)


def four_squares(R, rng=None):
    """R = a^2+b^2+c^2+d^2 with a,b,c,d >= 0 integers (Lagrange)."""
    if R < 0:
        raise ValueError("negative input")
    if R == 0:
        return (0, 0, 0, 0)
    rng = rng or random.Random(0xC0FFEE)
    scale = 1
    while R % 4 == 0:
        R //= 4
        scale *= 2
    if R <= 20000:                                       # exhaustive, for tests
        for a in range(isqrt(R), -1, -1):
            ra = R - a * a
            for b in range(isqrt(ra), -1, -1):
                rb = ra - b * b
                for c in range(isqrt(rb), -1, -1):
                    d = isqrt(rb - c * c)
                    if d * d == rb - c * c:
                        return tuple(scale * v for v in (a, b, c, d))
        raise RuntimeError("unreachable")
    for _ in range(200000):                              # Rabin-Shallit style
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


# ------------------------------------------------------------- rounding search

def find_columns(theta, eps, kmax=400):
    """Integer v0, v1 in Z^8, |vi|^2 = 2^k, v0.v1 = 0, column error <= eps."""
    c, s = math.cos(theta / 2.0), math.sin(theta / 2.0)
    ac, as_ = abs(c), abs(s)
    sc = 1 if c >= 0 else -1
    ss = 1 if s >= 0 else -1
    phi = math.atan2(as_, ac)
    for k in range(0, kmax):
        N = 1 << k
        rad = math.sqrt(2.0) ** k
        w = max(4, int(2.0 * eps * rad) + 4)
        pc, pmax = int(ac * rad), isqrt(N)
        best = None
        for off in range(0, w + 1):
            for p in ({pc + off, pc - off} if off else {pc}):
                if p < 0 or p > pmax:
                    continue
                q = isqrt(N - p * p)
                # ||v/sqrt2^k - u||^2 = 2(1-rho) + 2 rho (1-cos d), written so
                # that neither term cancels: nu^2 = R/2^k, 1-rho = nu^2/(1+rho).
                nu2 = math.ldexp(float(N - p * p - q * q), -k)
                rho = math.sqrt(max(0.0, 1.0 - nu2))
                d = math.atan2(q, p) - phi
                err = math.sqrt(max(0.0, 2.0 * nu2 / (1.0 + rho)
                                   + 4.0 * rho * math.sin(d / 2) ** 2))
                if best is None or err < best[0]:
                    best = (err, p, q)
            if best is not None and best[0] <= eps:
                break
        if best is None:
            continue
        err, p, q = best
        if err <= eps:
            r = four_squares(N - p * p - q * q)
            p, q = sc * p, ss * q
            v0 = [p, q, r[0], r[1], r[2], r[3], 0, 0]
            v1 = [-q, p, -r[1], r[0], -r[3], r[2], 0, 0]
            assert sum(x * x for x in v0) == N == sum(x * x for x in v1)
            assert sum(a * b for a, b in zip(v0, v1)) == 0
            return v0, v1, k, err
    raise RuntimeError("no solution found; increase kmax")


# --------------------------------------------------------------------- gates
# ('H',t) | ('X',t) | ('Z',t) | ('CCX',c1,c2,t).  Qubit 0 = data = LSB of index.

def _bit(x, j):
    return (x >> j) & 1


def _swap_gate(u, v):
    """Swap basis states u,v differing in exactly one bit (X-conjugated CCX)."""
    t = (u ^ v).bit_length() - 1
    o = [j for j in range(3) if j != t]
    pre = [('X', j) for j in o if _bit(u, j) == 0]
    return pre + [('CCX', o[0], o[1], t)] + pre


def _transposition(u, v):
    if u == v:
        return []
    if bin(u ^ v).count('1') == 1:
        return _swap_gate(u, v)
    t = (u ^ v).bit_length() - 1
    w = u ^ (1 << t)
    a = _swap_gate(u, w)
    return a + _transposition(w, v) + a        # (u v) = (u w)(w v)(u w)


def _perm_gates(perm):
    """Gates whose matrix P satisfies P[perm[x]][x] = 1, i.e. |x> -> |perm[x]>."""
    p, gates = [0] * 8, []
    for x in range(8):
        p[perm[x]] = x                        # the loop below builds the inverse
    for x in range(8):
        while p[x] != x:
            y = p[x]
            gates = _transposition(x, y) + gates
            p[x], p[y] = p[y], p[x]
    return gates


def _sign_gates(x):
    """Flip the sign of basis state x only (X-conjugated CCZ = X's + H CCX H)."""
    pre = [('X', j) for j in range(3) if _bit(x, j) == 0]
    return pre + [('H', 2), ('CCX', 0, 1, 2), ('H', 2)] + pre


def gate_matrix(g):
    """(integer 8x8 numerator, lde) of a single gate."""
    M = [[0] * 8 for _ in range(8)]
    if g[0] == 'H':
        t = g[1]
        for x in range(8):
            M[x][x] = 1 if _bit(x, t) == 0 else -1
            M[x][x ^ (1 << t)] = 1
        return M, 1
    if g[0] == 'X':
        for x in range(8):
            M[x ^ (1 << g[1])][x] = 1
        return M, 0
    if g[0] == 'Z':
        for x in range(8):
            M[x][x] = -1 if _bit(x, g[1]) else 1
        return M, 0
    if g[0] == 'CCX':
        c1, c2, t = g[1:]
        for x in range(8):
            M[x ^ (1 << t) if (_bit(x, c1) and _bit(x, c2)) else x][x] = 1
        return M, 0
    raise ValueError(g)


def apply_gates(V, k, gates):
    """Apply a time-ordered gate list on the left of V/sqrt(2)^k (V is 8 x m)."""
    m = len(V[0])
    for g in gates:
        G, kg = gate_matrix(g)
        V = [[sum(G[i][t] * V[t][j] for t in range(8)) for j in range(m)]
             for i in range(8)]
        k += kg
        while k >= 2 and all(v % 2 == 0 for r in V for v in r):
            V = [[v // 2 for v in r] for r in V]
            k -= 2
    return V, k


def circuit_matrix(gates):
    I = [[1 if i == j else 0 for j in range(8)] for i in range(8)]
    return apply_gates(I, 0, gates)


# --------------------------------------------------------- exact synthesis

def _matchings(groups):
    """All perfect matchings of the parity classes (list of lists of indices)."""
    def pair_up(g):
        if not g:
            yield []
            return
        a, rest = g[0], g[1:]
        for i in range(len(rest)):
            b = rest[i]
            for tail in pair_up(rest[:i] + rest[i + 1:]):
                yield [(a, b)] + tail
    out = [[]]
    for g in groups:
        out = [acc + m for acc in out for m in pair_up(g)]
        if len(out) > 400:                      # keep the search bounded
            out = out[:400]
    return out


def _choose_move(V):
    """Qubit j and permutation (or None) so that H_j lowers the lde, min cost."""
    par = [tuple(v & 1 for v in row) for row in V]
    cls = {}
    for x, pv in enumerate(par):
        cls.setdefault(pv, []).append(x)
    for g in cls.values():
        if len(g) % 2:
            raise RuntimeError("odd parity class -- cannot happen for isometries")
    for j in range(3):
        if all(par[x] == par[x ^ (1 << j)] for x in range(8)):
            return j, None                                     # free move
    best = None
    from itertools import permutations
    for j in range(3):
        slots = [(x, x ^ (1 << j)) for x in range(8) if not _bit(x, j)]
        for pairs in _matchings(list(cls.values())):
            for order in permutations(range(4)):
                perm = [0] * 8
                for pi, si in enumerate(order):
                    (a, b), (s0, s1) = pairs[pi], slots[si]
                    # orient to keep elements near their current index
                    if abs(a - s0) + abs(b - s1) <= abs(a - s1) + abs(b - s0):
                        perm[a], perm[b] = s0, s1
                    else:
                        perm[a], perm[b] = s1, s0
                cost = len(_perm_gates(perm))
                if best is None or cost < best[0]:
                    best = (cost, j, perm)
                    if cost == 0:
                        return j, None
    return best[1], best[2]


def _hadamard_step(V, k, j):
    new = [None] * 8
    for x in range(8):
        if _bit(x, j):
            continue
        y = x | (1 << j)
        a, b = V[x], V[y]
        assert all((u + v) % 2 == 0 for u, v in zip(a, b)), "parity mismatch"
        new[x] = [(u + v) // 2 for u, v in zip(a, b)]
        new[y] = [(u - v) // 2 for u, v in zip(a, b)]
    return new, k - 1


def synthesize_isometry(v0, v1, k):
    """Circuit (time order) over {H,X,Z,CCX} with columns 0,1 = v0,v1/sqrt2^k."""
    V = [[v0[i], v1[i]] for i in range(8)]
    red = []
    while k > 0:
        if k >= 2 and all(v % 2 == 0 for r in V for v in r):
            V = [[v // 2 for v in r] for r in V]
            k -= 2
            continue
        j, perm = _choose_move(V)
        if perm is not None:
            g = _perm_gates(perm)
            V, k = apply_gates(V, k, g)
            red += g
        V, k = _hadamard_step(V, k, j)
        red += [('H', j)]
    nz = []
    for c in (0, 1):
        z = [(i, V[i][c]) for i in range(8) if V[i][c] != 0]
        assert len(z) == 1 and abs(z[0][1]) == 1, "not a signed basis vector"
        nz.append(z[0])
    (a, sa), (b, sb) = nz
    tail = []
    if a != 0:
        tail += _transposition(a, 0)
        if b == 0:
            b = a
    if b != 1:
        tail += _transposition(b, 1)
    if sa < 0:
        tail += _sign_gates(0)
    if sb < 0:
        tail += _sign_gates(1)
    red += tail
    return list(reversed(red))                # invert: every gate is an involution


# ------------------------------------------------------------------ top level

def synthesize_ry(theta, eps, peephole=True):
    v0, v1, k, bound = find_columns(theta, eps)
    gates = synthesize_isometry(v0, v1, k)
    if peephole:
        gates = optimize(gates)
    return {'gates': gates, 'k': k, 'v0': v0, 'v1': v1, 'bound': bound,
            'counts': gate_counts(gates), 'theta': theta, 'eps': eps}


def gate_counts(gates):
    c = {}
    for g in gates:
        c[g[0]] = c.get(g[0], 0) + 1
    c['total'] = len(gates)
    return c


def clean_ancilla_error(gates, theta):
    """sup_{|psi|=1} || C(|psi>|00>) - (R_y(theta)|psi>)|00> ||, from the circuit."""
    import numpy as np
    M, k = circuit_matrix(gates)
    r = math.sqrt(2.0) ** k
    A = np.array([[M[i][0] / r, M[i][1] / r] for i in range(8)])
    c, s = math.cos(theta / 2), math.sin(theta / 2)
    B = np.zeros((8, 2))
    B[0, 0], B[1, 0], B[0, 1], B[1, 1] = c, s, -s, c
    return float(np.linalg.norm(A - B, 2))


def print_circuit(gates, limit=None):
    out = []
    for i, g in enumerate(gates):
        if limit and i >= limit:
            out.append("... (+%d more)" % (len(gates) - limit))
            break
        out.append("%s(%s)" % (g[0], ",".join("q%d" % q for q in g[1:])))
    return " ".join(out)


# ------------------------------------------------------------------ peephole

def optimize(gates):
    """Cancel adjacent identical self-inverse gates. Exact, matrix-preserving.

    Every gate in the set is an involution, so g g = I whenever two identical
    gates are adjacent, and X gates on distinct qubits commute past each other,
    which lets more pairs meet.
    """
    out = []
    for g in gates:
        if out and out[-1] == g:
            out.pop()
            continue
        if g[0] == 'X':                    # slide X back past commuting gates
            i = len(out) - 1
            while i >= 0 and _commutes(out[i], g):
                if out[i] == g:
                    out.pop(i)
                    break
                i -= 1
            else:
                out.append(g)
            continue
        out.append(g)
    return out


def _commutes(a, b):
    """Sufficient (conservative) commutation test for an X gate b."""
    if a[0] == 'X':
        return True
    if a[0] == 'Z':
        return a[1] != b[1]
    if a[0] == 'H':
        return a[1] != b[1]
    if a[0] == 'CCX':
        return b[1] not in a[1:]
    return False
