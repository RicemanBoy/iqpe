"""
rysynth -- approximate synthesis of R_y(theta) over {H, X, Z, CCX} (optionally
also CX) on 1 data qubit + 2 clean ancillas, minimising the CCX count.

PIPELINE
--------
1. Rounding (find_columns).  Find the smallest lde k and integers p, q with
       || v0/sqrt(2)^k - (cos(theta/2), sin(theta/2), 0, ..., 0) || <= eps,
   absorbing the residue R = 2^k - p^2 - q^2 with Lagrange's four squares:
       v0 = (p, q, r1, r2, r3, r4, 0, 0),  v1 = (-q, p, -r2, r1, -r4, r3, 0, 0).
   This is the ONLY approximation; everything after it is exact, so eps is the
   error of the final circuit on inputs |psi>|00>.

2. Exact synthesis of the 8x2 isometry V = [v0 v1] (synthesize_isometry).
   Repeatedly: permute rows so equal-parity rows become bit-j partners, then one
   Hadamard H_j lowers k by 1.  At k = 0 fix signs and positions, then reverse.

CCX MINIMISATION
----------------
Every CCX in the circuit comes from a basis permutation.  On first use, a
Dijkstra search over the whole permutation group S_8 (40320 elements) computes
for EVERY permutation a circuit that minimises the CCX count first and the
number of other gates (X, and CX if allowed) second.  Then:

  * every permutation emitted is CCX-optimal for its gate set,
  * at each reduction step the algorithm enumerates ALL permutations that make
    a Hadamard legal (every qubit j, every parity matching, every assignment of
    pairs to slots, every orientation) and takes the one with fewest CCX,
  * the final positioning permutation is the CCX-cheapest of all 720 that send
    the two surviving basis states to |0>, |1>,
  * sign corrections use only X and Z gates (zero CCX): they only have to act
    correctly on the two surviving basis vectors, not on all eight.

Optionally, `tries` tries several equivalent sign/order arrangements of the
residue (r1..r4), which leave the error unchanged, and keeps the circuit with
the fewest CCX.

What is and is not guaranteed: each permutation is optimal and each step is the
CCX-minimal choice, i.e. the reduction is greedy-optimal.  It is not a proof
that no shorter circuit for the same approximation exists anywhere.

BEAM SEARCH (optional, off by default)
---------------------------------------
synthesize_ry(..., beam=n) replaces the one-move-per-step greedy reduction by a
beam search over reduction paths that keeps the n cheapest partial reductions
per lde level.  beam=None (the default) gives exactly the previous output.

Allowed gates: ('H',t) ('X',t) ('Z',t) ('CCX',c1,c2,t), and with allow_cx=True
also ('CX',c,t).  Qubit 0 = data = least significant bit of the basis index
x = 4*q2 + 2*q1 + q0; q1, q2 are ancillas prepared in |00>.
"""

import heapq
import itertools
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

def _bit(x, j):
    return (x >> j) & 1


def _perm_of(g):
    """Basis permutation of a permutation gate: tuple P with |x> -> |P[x]>."""
    out = []
    for x in range(8):
        if g[0] == 'X':
            y = x ^ (1 << g[1])
        elif g[0] == 'CX':
            y = x ^ (1 << g[2]) if _bit(x, g[1]) else x
        elif g[0] == 'CCX':
            y = x ^ (1 << g[3]) if (_bit(x, g[1]) and _bit(x, g[2])) else x
        else:
            raise ValueError("not a permutation gate: %r" % (g,))
        out.append(y)
    return tuple(out)


def gate_matrix(g):
    """(integer 8x8 numerator, lde) of a single gate."""
    M = [[0] * 8 for _ in range(8)]
    if g[0] == 'H':
        t = g[1]
        for x in range(8):
            M[x][x] = 1 if _bit(x, t) == 0 else -1
            M[x][x ^ (1 << t)] = 1
        return M, 1
    if g[0] == 'Z':
        for x in range(8):
            M[x][x] = -1 if _bit(x, g[1]) else 1
        return M, 0
    P = _perm_of(g)
    for x in range(8):
        M[P[x]][x] = 1
    return M, 0


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


# ------------------------------------------------- optimal permutation tables

_CCX_WEIGHT = 1000            # minimise CCX first, other gates second
_TABLES = {}


def _gate_set(allow_cx):
    gs = [('X', q) for q in range(3)]
    gs += [('CCX',) + tuple(q for q in range(3) if q != t) + (t,) for t in range(3)]
    if allow_cx:
        gs += [('CX', c, t) for c in range(3) for t in range(3) if c != t]
    return gs


def _tables(allow_cx):
    """Dijkstra over S_8: cost and predecessor for every permutation (cached)."""
    key = bool(allow_cx)
    if key in _TABLES:
        return _TABLES[key]
    gs = _gate_set(allow_cx)
    perms = {g: _perm_of(g) for g in gs}
    w = {g: (_CCX_WEIGHT if g[0] == 'CCX' else 1) for g in gs}
    ident = tuple(range(8))
    dist, prev, pq = {ident: 0}, {ident: None}, [(0, ident)]
    while pq:
        d, A = heapq.heappop(pq)
        if d > dist[A]:
            continue
        for g in gs:
            G = perms[g]
            B = tuple(G[A[x]] for x in range(8))
            nd = d + w[g]
            if nd < dist.get(B, 1 << 62):
                dist[B], prev[B] = nd, (A, g)
                heapq.heappush(pq, (nd, B))
    assert len(dist) == 40320, "permutation group not fully generated"
    _TABLES[key] = (dist, prev)
    return _TABLES[key]


def perm_cost(perm, allow_cx=False):
    """(CCX count, other-gate count) of the optimal circuit for `perm`."""
    d = _tables(allow_cx)[0][tuple(perm)]
    return divmod(d, _CCX_WEIGHT)


def _perm_gates(perm, allow_cx=False):
    """CCX-optimal gates (time order) whose matrix P has P[perm[x]][x] = 1."""
    prev = _tables(allow_cx)[1]
    A, seq = tuple(perm), []
    while prev[A] is not None:
        A, g = prev[A]
        seq.append(g)
    return seq[::-1]


def _transposition(u, v, allow_cx=False):
    """CCX-optimal gates swapping basis states u and v."""
    p = list(range(8))
    p[u], p[v] = v, u
    return _perm_gates(p, allow_cx)


def _sign_gates(x):
    """Legacy: flip the sign of basis state x alone (X-conjugated CCZ, 1 CCX).
    No longer used by the synthesis, which needs zero-CCX sign fixes only."""
    pre = [('X', j) for j in range(3) if _bit(x, j) == 0]
    return pre + [('H', 2), ('CCX', 0, 1, 2), ('H', 2)] + pre


def _sign_fix(a, sa, b, sb):
    """X/Z gates D with D(sa e_a) = e_a and D(sb e_b) = e_b.  Zero CCX.

    D only has to be right on the two vectors e_a, e_b (a != b), not on all
    eight basis states, so Pauli Z's (conjugated by X) suffice."""
    if sa > 0 and sb > 0:
        return []
    d = a ^ b
    j = (d & -d).bit_length() - 1                     # a bit where a, b differ

    def flip(x):                                      # flips x, not the other one
        return [('Z', j)] if _bit(x, j) else [('X', j), ('Z', j), ('X', j)]

    if sa < 0 < sb:
        return flip(a)
    if sb < 0 < sa:
        return flip(b)
    for i in range(3):                                # both negative
        if _bit(a, i) == _bit(b, i):
            return [('Z', i)] if _bit(a, i) else [('X', i), ('Z', i), ('X', i)]
    return flip(a) + flip(b)                          # a, b differ in every bit


# --------------------------------------------------------- exact synthesis

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


def _matchings(groups):
    """All perfect matchings of the parity classes."""
    def pair_up(g):
        if not g:
            yield []
            return
        a, rest = g[0], g[1:]
        for i in range(len(rest)):
            for tail in pair_up(rest[:i] + rest[i + 1:]):
                yield [(a, rest[i])] + tail
    out = [[]]
    for g in groups:
        out = [acc + m for acc in out for m in pair_up(g)]
    return out


def _choose_move(V, allow_cx=False):
    """(j, perm) with the fewest CCX among ALL permutations making H_j legal.
    perm is None when no permutation is needed."""
    par = [tuple(v & 1 for v in row) for row in V]
    for j in range(3):
        if all(par[x] == par[x ^ (1 << j)] for x in range(8)):
            return j, None                                       # free move
    cls = {}
    for x, pv in enumerate(par):
        cls.setdefault(pv, []).append(x)
    for g in cls.values():
        if len(g) % 2:
            raise RuntimeError("odd parity class -- impossible for an isometry")
    dist = _tables(allow_cx)[0]
    best = None
    matchings = _matchings(list(cls.values()))
    for j in range(3):
        slots = [(x, x | (1 << j)) for x in range(8) if not _bit(x, j)]
        for m in matchings:
            for order in itertools.permutations(range(4)):
                for orient in range(16):
                    perm = [0] * 8
                    for i, si in enumerate(order):
                        a, b = m[i]
                        if (orient >> i) & 1:
                            a, b = b, a
                        perm[a], perm[b] = slots[si]
                    c = dist[tuple(perm)]
                    if best is None or c < best[0]:
                        best = (c, j, perm)
    return best[1], best[2]


def synthesize_isometry(v0, v1, k, allow_cx=False, beam=None, branch=4):
    """Circuit (time order) whose columns 0,1 are exactly v0, v1 / sqrt(2)^k.

    beam   : None (default) -> greedy reduction, one CCX-cheapest move per step.
             int n >= 1     -> beam search keeping the n cheapest partial
                               reductions per lde level (slower, fewer CCX).
    branch : moves expanded per kept state in beam mode (ignored when greedy).

    Returns the raw list; synthesize_ry additionally runs the peephole pass."""
    width = _beam_width(beam)
    if width is None:
        return _synthesize_isometry_greedy(v0, v1, k, allow_cx)
    return _synthesize_isometry_beam(v0, v1, k, allow_cx, width,
                                     _branch_width(branch))


def _synthesize_isometry_greedy(v0, v1, k, allow_cx=False):
    """The default reduction (unchanged from the previous version)."""
    V = [[v0[i], v1[i]] for i in range(8)]
    red = []
    while k > 0:
        if k >= 2 and all(v % 2 == 0 for r in V for v in r):
            V = [[v // 2 for v in r] for r in V]
            k -= 2
            continue
        j, perm = _choose_move(V, allow_cx)
        if perm is not None:
            g = _perm_gates(perm, allow_cx)
            V, k = apply_gates(V, k, g)
            red += g
        V, k = _hadamard_step(V, k, j)
        red.append(('H', j))
    # k == 0: V = [sa e_a, sb e_b]
    (a, sa), (b, sb) = [next((i, V[i][c]) for i in range(8) if V[i][c] != 0)
                        for c in (0, 1)]
    assert abs(sa) == 1 and abs(sb) == 1 and a != b, "not a signed basis pair"
    red += _sign_fix(a, sa, b, sb)
    dist = _tables(allow_cx)[0]
    best = None
    others = [x for x in range(8) if x not in (a, b)]
    for rest in itertools.permutations(range(2, 8)):    # all perms a->0, b->1
        perm = [0] * 8
        perm[a], perm[b] = 0, 1
        for x, y in zip(others, rest):
            perm[x] = y
        c = dist[tuple(perm)]
        if best is None or c < best[0]:
            best = (c, perm)
    red += _perm_gates(best[1], allow_cx)
    return list(reversed(red))               # every gate is an involution



# -------------------------------------------------------------- beam search
# Partial reductions are bucketed by their current lde k.  Buckets are processed
# from high k to low; each keeps its `beam` cheapest states (accumulated cost,
# CCX weighted _CCX_WEIGHT, other gates 1) and expands each with its `branch`
# cheapest distinct moves.  Only states with the same k are compared, so no
# estimate of the remaining cost is needed.  beam=1, branch=1 reproduces the
# greedy reduction exactly.

DEFAULT_BEAM = 16


def _beam_width(beam):
    """Normalise the `beam` argument: None/False/0 -> None (greedy)."""
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


def _normalise(V, k):
    while k >= 2 and all(v % 2 == 0 for r in V for v in r):
        V = [[v // 2 for v in r] for r in V]
        k -= 2
    return V, k


def _all_moves(V, allow_cx):
    """Every (cost, j, perm) that makes H_j legal, in greedy's order."""
    par = [tuple(v & 1 for v in row) for row in V]
    out = []
    for j in range(3):                                     # free moves first
        if all(par[x] == par[x ^ (1 << j)] for x in range(8)):
            out.append((0, j, None))
    cls = {}
    for x, pv in enumerate(par):
        cls.setdefault(pv, []).append(x)
    dist = _tables(allow_cx)[0]
    matchings = _matchings(list(cls.values()))
    for j in range(3):
        slots = [(x, x | (1 << j)) for x in range(8) if not _bit(x, j)]
        for m in matchings:
            for order in itertools.permutations(range(4)):
                for orient in range(16):
                    perm = [0] * 8
                    for i, si in enumerate(order):
                        a, b = m[i]
                        if (orient >> i) & 1:
                            a, b = b, a
                        perm[a], perm[b] = slots[si]
                    out.append((dist[tuple(perm)], j, tuple(perm)))
    return out


def _apply_move(V, k, j, perm):
    if perm is not None:
        W = [None] * 8
        for x in range(8):
            W[perm[x]] = V[x]
        V = W
    V, k = _hadamard_step(V, k, j)
    return _normalise(V, k)


def _final_gates(V, allow_cx):
    """Sign fix + CCX-cheapest positioning permutation at k = 0: (cost, gates)."""
    (a, sa), (b, sb) = [next((i, V[i][c]) for i in range(8) if V[i][c] != 0)
                        for c in (0, 1)]
    assert abs(sa) == 1 and abs(sb) == 1 and a != b, "not a signed basis pair"
    g = _sign_fix(a, sa, b, sb)
    dist = _tables(allow_cx)[0]
    others = [x for x in range(8) if x not in (a, b)]
    best = None
    for rest in itertools.permutations(range(2, 8)):
        perm = [0] * 8
        perm[a], perm[b] = 0, 1
        for x, y in zip(others, rest):
            perm[x] = y
        c = dist[tuple(perm)]
        if best is None or c < best[0]:
            best = (c, perm)
    return best[0] + len(g), g + _perm_gates(best[1], allow_cx)


def _synthesize_isometry_beam(v0, v1, k, allow_cx, beam, branch, finalists=8):
    V0, k0 = _normalise([[v0[i], v1[i]] for i in range(8)], k)
    nodes = [(0, V0, k0, None, None)]         # (cost, V, k, parent, (j, perm))
    buckets = {k0: {tuple(map(tuple, V0)): 0}}
    done = []
    for kk in range(k0, -1, -1):
        bucket = buckets.pop(kk, {})
        ids = sorted(bucket.values(), key=lambda i: (nodes[i][0], i))[:beam]
        for nid in ids:
            cost, V = nodes[nid][0], nodes[nid][1]
            if kk == 0:
                tc, tg = _final_gates(V, allow_cx)
                done.append((cost + tc, nid, tg))
                continue
            seen, taken = set(), 0
            for mc, j, perm in sorted(_all_moves(V, allow_cx), key=lambda t: t[0]):
                W, kw = _apply_move(V, kk, j, perm)
                key = tuple(map(tuple, W))
                if key in seen:
                    continue
                seen.add(key)
                c = cost + mc + 1                                   # +1 for H
                b = buckets.setdefault(kw, {})
                if key not in b or c < nodes[b[key]][0]:
                    nodes.append((c, W, kw, nid, (j, perm)))
                    b[key] = len(nodes) - 1
                taken += 1
                if taken >= branch:
                    break
    best = None
    for _, nid, tg in sorted(done, key=lambda t: t[0])[:finalists]:
        red, i = [], nid
        while nodes[i][3] is not None:
            j, perm = nodes[i][4]
            red = ([] if perm is None else _perm_gates(list(perm), allow_cx)) \
                + [('H', j)] + red
            i = nodes[i][3]
        raw = list(reversed(red + tg))          # every gate is an involution
        c = gate_counts(optimize(raw))          # judge by the final circuit
        score = (c.get('CCX', 0), c['total'])
        if best is None or score < best[0]:
            best = (score, raw)
    return best[1]


# ------------------------------------------------------------------ peephole

def _pauli_commutes(a, g):
    """Does gate a commute with the Pauli g = ('X',q) or ('Z',q)?"""
    kind, q = g
    if a[0] in ('X', 'Z'):
        return a[0] == kind or a[1] != q
    if a[0] == 'H':
        return a[1] != q
    if a[0] == 'CX':                       # X on target / Z on control commute
        return q != a[1] if kind == 'X' else q != a[2]
    if a[0] == 'CCX':
        return q not in a[1:3] if kind == 'X' else q != a[3]
    return False


def optimize(gates):
    """Exact peephole pass: cancel identical self-inverse gates, sliding X and Z
    back past gates they commute with so that more pairs meet."""
    out = []
    for g in gates:
        if g[0] in ('X', 'Z'):
            i = len(out) - 1
            while i >= 0 and _pauli_commutes(out[i], g):
                if out[i] == g:
                    out.pop(i)
                    break
                i -= 1
            else:
                out.append(g)
            continue
        if out and out[-1] == g:
            out.pop()
            continue
        out.append(g)
    return out


# ------------------------------------------------------------------ top level

def _variants(v0, n, rng):
    """Equivalent column pairs: same p, q and residue, rearranged tail.
    The error depends only on p, q and R, so every variant has identical error."""
    p, q, tail = v0[0], v0[1], list(v0[2:6])
    out = [list(v0)]
    seen = {tuple(v0)}
    pairs = [(2, 3), (4, 5), (6, 7)]
    tries = 0
    while len(out) < n and tries < 50 * n:
        tries += 1
        t = tail[:]
        rng.shuffle(t)
        t = [x * rng.choice((1, -1)) for x in t]
        empty = rng.randrange(3)
        used = [pr for i, pr in enumerate(pairs) if i != empty]
        w = [p, q, 0, 0, 0, 0, 0, 0]
        (i0, i1), (i2, i3) = used
        w[i0], w[i1], w[i2], w[i3] = t
        if tuple(w) not in seen:
            seen.add(tuple(w))
            out.append(w)
    return out


def _companion(v0):
    """v1 = J v0: rotate every consecutive pair (a, b) -> (-b, a)."""
    v1 = [0] * 8
    for i in range(0, 8, 2):
        v1[i], v1[i + 1] = -v0[i + 1], v0[i]
    return v1


def synthesize_ry(theta, eps, peephole=True, allow_cx=False, tries=1, seed=0,
                  beam=None, branch=4):
    """Approximate R_y(theta) to precision eps, minimising the CCX count.

    allow_cx : also use CX gates (fewer CCX, slightly more total gates).
    tries    : number of equivalent residue arrangements to try (same error,
               different circuits); the one with fewest CCX, then fewest
               gates, is returned.  tries=1 is fastest.
    beam     : None (default) -> greedy exact synthesis, as before.
               int n >= 1     -> beam search over reduction paths keeping n
                                 candidates per lde level; larger n gives fewer
                                 CCX (roughly -3% per 4x beam without CX) at
                                 about n/4 times the greedy runtime.
               True           -> beam = DEFAULT_BEAM (16).
               With allow_cx=True the CCX count is essentially fixed and the
               beam reduces the total gate count instead.
    branch   : moves expanded per kept state in beam mode (default 4).
    """
    width = _beam_width(beam)
    branch = _branch_width(branch)
    v0, v1, k, bound = find_columns(theta, eps)
    best = None
    for w0 in _variants(v0, max(1, tries), random.Random(seed)):
        w1 = _companion(w0)
        assert sum(x * x for x in w0) == 1 << k
        assert sum(x * y for x, y in zip(w0, w1)) == 0
        g = synthesize_isometry(w0, w1, k, allow_cx, beam=width, branch=branch)
        if peephole:
            g = optimize(g)
        c = gate_counts(g)
        score = (c.get('CCX', 0), c['total'])
        if best is None or score < best[0]:
            best = (score, g, w0, w1)
    _, gates, v0, v1 = best
    return {'gates': gates, 'k': k, 'v0': v0, 'v1': v1, 'bound': bound,
            'counts': gate_counts(gates), 'theta': theta, 'eps': eps,
            'allow_cx': bool(allow_cx), 'beam': width, 'branch': branch}


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
