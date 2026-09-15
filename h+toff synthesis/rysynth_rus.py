"""
rysynth_rus -- repeat-until-success synthesis of R_y(theta) over {H, X, Z, CCX}.

Companion to rysynth.py (which it imports for the exact-synthesis engine, the
four-square routine and the peephole optimiser).

IDEA
----
The deterministic algorithm needs  a^2 + b^2 = 2^k  EXACTLY, because the residue
R = 2^k - a^2 - b^2 leaks amplitude into the ancillas and therefore has to obey
R <= eps^2 * 2^k.  That two-sided constraint (right direction AND norm almost
exactly a power of two) is what forces k ~ 3 log2(1/eps).

If the two ancillas are measured, R stops being error and becomes failure
probability.  Only the direction constraint survives,

        a / b  ~  cot(theta/2)   to precision eps,

which needs a, b ~ eps^(-1/2), hence  k ~ log2(1/eps):  three times smaller.

Counting lattice points in the annular sector {angle within eps, radius^2 in
[p*2^k, 2^k]} gives area eps*(1-p)*2^k >= 1, i.e. k = log2(1/eps) + O(1) for any
fixed success probability p bounded away from 1.

FALLBACK STRUCTURE
------------------
With the quaternionic pairing used for v1, EVERY measurement outcome m has
Kraus operator proportional to [[r1,-r2],[r2,r1]], i.e. an exact rotation by a
known angle.  So a "failure" never damages the data: it applies a known exact
R_y(phi_m), and the next round simply retargets  theta - phi_m.  This is
repeat-until-success with a free fallback, so intermediate rounds may be run at
coarse precision -- only the final, successful round must meet eps.

LAYOUT
------
q0 = data, q1 q2 = ancillas, prepared in |00> and measured at the end of each
round.  Basis index x = 4*q2 + 2*q1 + q0, so the success subspace is {0,1}.
"""

import math
import random
from fractions import Fraction
from math import isqrt

import rysynth as _rs

__all__ = ["find_columns_rus", "synthesize_rus", "branches", "run_adaptive",
           "expected_cost", "verify"]


# --------------------------------------------------------------- rounding step

def find_columns_rus(theta, eps, pmin=0.5, kmax=200):
    """Integer columns for the RUS circuit.

    Returns (v0, v1, k, angle_err, p_success) where the success-branch Kraus
    operator is [[a,-b],[b,a]]/sqrt(2)^k, an exact rotation by 2*atan2(b,a)
    after normalisation, fired with probability (a^2+b^2)/2^k.

    Among all (a, b) whose angle is within eps, the one with the LARGEST norm is
    taken, i.e. success probability is maximised subject to the accuracy budget.
    """
    c, s = math.cos(theta / 2.0), math.sin(theta / 2.0)
    ac, as_ = abs(c), abs(s)
    sc = 1 if c >= 0 else -1
    ss = 1 if s >= 0 else -1
    phi = math.atan2(as_, ac)
    t = math.tan(phi) if ac > 1e-15 else None
    for k in range(0, kmax):
        N = 1 << k
        rad = math.sqrt(N)
        best = None
        amax = int(ac * rad) + 2
        amin = max(0, int(ac * math.sqrt(pmin) * rad) - 2)
        for a in range(amin, amax + 1):
            b = isqrt(N - a * a) if t is None else round(a * t)
            for bb in ({b, b + 1, b - 1} if t is not None else {b}):
                if bb < 0:
                    continue
                n = a * a + bb * bb
                if n > N or n < pmin * N:
                    continue
                err = 2.0 * abs(math.atan2(bb, a) - phi)     # error on theta
                if err <= eps and (best is None or n > best[3]):
                    best = (err, a, bb, n)
        if best is not None:
            err, a, b, n = best
            r = _rs.four_squares(N - n)
            A, B = sc * a, ss * b
            v0 = [A, B, r[0], r[1], r[2], r[3], 0, 0]
            v1 = [-B, A, -r[1], r[0], -r[3], r[2], 0, 0]
            assert sum(x * x for x in v0) == N == sum(x * x for x in v1)
            assert sum(x * y for x, y in zip(v0, v1)) == 0
            return v0, v1, k, err, Fraction(n, N)
    raise RuntimeError("no solution up to kmax=%d" % kmax)


# ------------------------------------------------------------------ one round

def synthesize_rus(theta, eps, pmin=0.5, peephole=True):
    """One RUS round. Gates are over {H,X,Z,CCX} on qubits q0,q1,q2."""
    v0, v1, k, err, psucc = find_columns_rus(theta, eps, pmin)
    gates = _rs.synthesize_isometry(v0, v1, k)
    if peephole:
        gates = _rs.optimize(gates)
    return {'gates': gates, 'k': k, 'v0': v0, 'v1': v1,
            'angle_err': err, 'p_success': psucc,
            'counts': _rs.gate_counts(gates), 'theta': theta, 'eps': eps}


def branches(gates):
    """Exact per-outcome data: [(m, probability, is_rotation, angle), ...].

    m is the ancilla outcome as an integer, m = 0 (that is |q1 q2> = |00>) being
    success. Probabilities are exact Fractions and are state-independent.
    """
    M, k = _rs.circuit_matrix(gates)
    out = []
    for m in range(4):
        i, j = 2 * m, 2 * m + 1
        A, B = M[i][0], M[i][1]
        C, D = M[j][0], M[j][1]
        p = Fraction(A * A + C * C, 1 << k)
        if p == 0:
            continue
        out.append((m, p, (A == D and B == -C), 2.0 * math.atan2(C, A)))
    return out


# ------------------------------------------------- adaptive repeat-until-success

def _wrap(x):
    """Fold an angle into (-2pi, 2pi], the natural period of R_y."""
    return (x + 2 * math.pi) % (4 * math.pi) - 2 * math.pi


def run_adaptive(theta, eps, rng=None, pmin=0.5, max_rounds=64, cache=None):
    """Simulate the full adaptive loop.

    Each round synthesises the REMAINING angle, samples a measurement outcome
    from the exact branch distribution, and -- on a non-success outcome --
    subtracts the known rotation that was applied and goes round again.

    Returns a dict with the realised gate cost, the rounds used, and the final
    angle error. Two X gates per round are charged for resetting the ancillas
    to |00> (classically conditioned on the outcome, still inside the gate set).
    """
    rng = rng or random.Random()
    cache = {} if cache is None else cache
    remaining = _wrap(theta)
    total, ccx, rounds, seq = 0, 0, 0, []
    for _ in range(max_rounds):
        key = (round(remaining, 15), eps, pmin)
        if key not in cache:
            cache[key] = synthesize_rus(remaining, eps, pmin)
        res = cache[key]
        br = branches(res['gates'])
        rounds += 1
        total += res['counts']['total']
        ccx += res['counts'].get('CCX', 0)
        seq.append(res)
        u, acc = rng.random(), 0.0
        for m, p, _isrot, ang in br:
            acc += float(p)
            if u <= acc:
                break
        if m == 0:                                    # success
            remaining = _wrap(remaining - ang)
            break
        remaining = _wrap(remaining - ang)            # fallback: known rotation
        total += 2                                    # ancilla reset (X gates)
    return {'rounds': rounds, 'gates': total, 'ccx': ccx,
            'angle_error': abs(remaining), 'rounds_detail': seq}


def expected_cost(theta, eps, pmin=0.5):
    """Expected gate cost of the loop, assuming each round is statistically
    like the first: cost / p_success."""
    res = synthesize_rus(theta, eps, pmin)
    p = float(res['p_success'])
    return {'k': res['k'], 'gates_per_round': res['counts']['total'],
            'ccx_per_round': res['counts'].get('CCX', 0), 'p_success': p,
            'expected_gates': res['counts']['total'] / p,
            'expected_ccx': res['counts'].get('CCX', 0) / p,
            'expected_rounds': 1.0 / p}


# ------------------------------------------------------------------ verification

def verify(res, theta=None, tol=1e-12):
    """Check the round: gate set, exactness, unitarity, branch structure."""
    import numpy as np
    g = res['gates']
    rep = {}
    rep['gate_set_ok'] = all(x[0] in ('H', 'X', 'Z', 'CCX') and
                             all(0 <= q < 3 for q in x[1:]) for x in g)
    M, k = _rs.circuit_matrix(g)
    rep['orthogonal'] = all(
        sum(M[i][m] * M[j][m] for m in range(8)) == ((1 << k) if i == j else 0)
        for i in range(8) for j in range(8))
    rep['columns_exact'] = ([[M[i][0], M[i][1]] for i in range(8)] ==
                            [[res['v0'][i], res['v1'][i]] for i in range(8)]
                            and k == res['k'])
    br = branches(g)
    rep['all_branches_are_rotations'] = all(b[2] for b in br)
    rep['probabilities_sum_to_1'] = sum(b[1] for b in br) == 1
    rep['p_success'] = float(next(b[1] for b in br if b[0] == 0))
    a, b_ = res['v0'][0], res['v0'][1]
    ang = 2.0 * math.atan2(b_, a)
    rep['success_angle'] = ang
    if theta is not None:
        # R_y has period 4*pi (R_y(t+2pi) = -R_y(t), a real global sign this
        # construction tracks exactly), so compare modulo 4*pi.
        d = (ang - theta) % (4 * math.pi)
        rep['angle_error'] = min(d, 4 * math.pi - d)
        rep['within_eps'] = rep['angle_error'] <= res['eps'] * (1 + 1e-9)
    # success branch is an exact rotation after normalising by sqrt(a^2+b^2)
    n = math.sqrt(a * a + b_ * b_)
    K = np.array([[M[0][0], M[0][1]], [M[1][0], M[1][1]]]) / n
    rep['success_branch_unitary'] = bool(
        np.allclose(K @ K.T, np.eye(2), atol=tol))
    return rep


# ------------------------------------------------------------------------ demo

if __name__ == "__main__":
    import statistics as st

    theta, eps = 0.4321, 1e-6
    print("theta = %.6f   eps = %.0e\n" % (theta, eps))

    res = synthesize_rus(theta, eps)
    print("k = %d   gates = %d   %s" %
          (res['k'], res['counts']['total'],
           {x: y for x, y in res['counts'].items() if x != 'total'}))
    print("v0 =", res['v0'])
    print("v1 =", res['v1'])
    print("first 14 gates:", _rs.print_circuit(res['gates'], 14), "\n")

    print("branches (ancilla |q1 q2>):")
    for m, p, isrot, ang in branches(res['gates']):
        print("   %s   p = %-10.6f exact rotation: %-5s  angle = %+.9f"
              % (format(m, '02b'), float(p), isrot, ang))
    print()

    for key, val in verify(res, theta).items():
        print("   %-28s %s" % (key, val))
    print()

    ec = expected_cost(theta, eps)
    print("expected rounds %.2f, expected gates %.0f, expected CCX %.0f"
          % (ec['expected_rounds'], ec['expected_gates'], ec['expected_ccx']))

    rng = random.Random(7)
    cache = {}
    runs = [run_adaptive(theta, eps, rng, cache=cache) for _ in range(300)]
    print("simulated 300 adaptive runs: mean rounds %.2f, mean gates %.0f, "
          "max angle error %.2e"
          % (st.mean(r['rounds'] for r in runs),
             st.mean(r['gates'] for r in runs),
             max(r['angle_error'] for r in runs)))
