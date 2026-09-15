"""Test suite for rysynth. Run: python3 test_rysynth.py"""
import math
import random
import numpy as np
import rysynth as R

PASS, FAIL = [], []


def check(name, cond, info=""):
    (PASS if cond else FAIL).append(name)
    print(("  ok   " if cond else "  FAIL ") + name + (("  " + info) if info else ""))


# ------------------------------------------------------------------ 1. number theory
print("1. number theory")
rng = random.Random(11)
bad = [n for n in range(0, 4000) if sum(x * x for x in R.four_squares(n)) != n]
big = [rng.randrange(10 ** 6, 10 ** 18) for _ in range(400)]
bad += [n for n in big if sum(x * x for x in R.four_squares(n)) != n]
check("four_squares(n) exact for 4000 small + 400 large n", not bad, str(bad[:3]))
primes = [p for p in range(5, 5000) if R.is_prime(p) and p % 4 == 1]
check("Cornacchia two_squares_prime on all p=1 mod 4 below 5000",
      all(sum(x * x for x in R.two_squares_prime(p)) == p for p in primes))
check("is_prime agrees with a sieve below 20000",
      [n for n in range(2, 20000) if R.is_prime(n)] ==
      [n for n in range(2, 20000)
       if all(n % d for d in range(2, int(n ** .5) + 1))])

# --------------------------------------------------- 2. group-membership invariant
print("2. the gate group is (integer matrix)/sqrt(2)^k, not O(Z[1/sqrt2])")
rng = random.Random(5)
ok = True
for _ in range(200):
    gates = []
    for _ in range(rng.randint(1, 12)):
        t = rng.randrange(3)
        gates.append(rng.choice([('H', t), ('X', t), ('Z', t),
                                 ('CCX', 0, 1, 2), ('CCX', 1, 2, 0)]))
    M, k = R.circuit_matrix(gates)
    MMt = [[sum(M[i][m] * M[j][m] for m in range(8)) for j in range(8)]
           for i in range(8)]
    ok &= all(MMt[i][j] == ((1 << k) if i == j else 0)
              for i in range(8) for j in range(8))
    ok &= all(isinstance(v, int) for r in M for v in r)
check("random circuits satisfy M integer and M M^T = 2^k I", ok)
J = np.array([[0, -1], [1, 0]])
Temb = np.block([[np.eye(2), np.zeros((2, 2))],
                 [np.zeros((2, 2)), (np.eye(2) + J) / math.sqrt(2)]])
sig = Temb.copy()                                     # sqrt2 -> -sqrt2
sig[2:, 2:] = (np.eye(2) - J) / math.sqrt(2) * -1 + 2 * np.eye(2) * 0  # recompute
sig = np.block([[np.eye(2), np.zeros((2, 2))],
                [np.zeros((2, 2)), (np.eye(2) + J) / (-math.sqrt(2))]])
check("i->J embedding of T is orthogonal over Z[1/sqrt2] ...",
      np.allclose(Temb @ Temb.T, np.eye(4)))
check("... but violates U^sigma = +-U, so it is NOT a Toffoli+H circuit",
      not (np.allclose(sig, Temb) or np.allclose(sig, -Temb)))

# --------------------------------------------------------- 3. ancilla-free no-go
print("3. no ancilla-free approximation exists (SO(2,Z[1/sqrt2]) is finite)")
angles = set()
for k in range(0, 20):
    N = 1 << k
    for p in range(0, math.isqrt(N) + 1):
        q2 = N - p * p
        q = math.isqrt(q2)
        if q * q == q2:
            angles.add(round(math.degrees(math.atan2(q, p)), 6))
check("2x2 integer M with M M^T = 2^k I only realises multiples of 45 deg",
      angles == {0.0, 45.0, 90.0}, str(sorted(angles)))

# --------------------------------------------------------------- 4. column search
print("4. rounding step")
ok = True
for th in [0.0, 0.1, 1.0, 2.0, math.pi / 2, math.pi, 3.0, 5.7]:
    for eps in [1e-2, 1e-4]:
        v0, v1, k, b = R.find_columns(th, eps)
        ok &= sum(x * x for x in v0) == (1 << k)
        ok &= sum(x * x for x in v1) == (1 << k)
        ok &= sum(a * b for a, b in zip(v0, v1)) == 0
        ok &= b <= eps
check("v0,v1 integral, |vi|^2=2^k, v0.v1=0, predicted error <= eps", ok)
v0, v1, k, b = R.find_columns(math.pi / 2, 1e-9)
check("theta=pi/2 is exact at k=1 (R_y(pi/2) is a Toffoli+H element)",
      k == 1 and b == 0.0, "k=%d" % k)
v0, v1, k, b = R.find_columns(math.pi, 1e-9)
check("theta=pi is exact at k=0", k == 0 and b == 0.0, "k=%d" % k)

# ------------------------------------------------------------ 5. exact synthesis
print("5. exact synthesis reproduces the two columns exactly (integer identity)")


def norm_lde(M, k):
    while k >= 2 and all(v % 2 == 0 for r in M for v in r):
        M = [[v // 2 for v in r] for r in M]
        k -= 2
    return M, k


ok, okset = True, True
rng = random.Random(2)
for _ in range(25):
    th = rng.uniform(0, 2 * math.pi)
    eps = 10 ** -rng.randint(1, 4)
    r = R.synthesize_ry(th, eps)
    M, kc = R.circuit_matrix(r['gates'])
    Vt = [[r['v0'][i], r['v1'][i]] for i in range(8)]
    A, ka = norm_lde([[M[i][0], M[i][1]] for i in range(8)], kc)
    B, kb = norm_lde(Vt, r['k'])
    ok &= (A == B and ka == kb)
    okset &= all(g[0] in ('H', 'X', 'Z', 'CCX') and all(0 <= q < 3 for q in g[1:])
                 for g in r['gates'])
check("circuit columns 0,1 == v0,v1/sqrt2^k, exactly, over 25 random cases", ok)
check("every emitted gate is in {H,X,Z,CCX} on qubits q0,q1,q2", okset)

# --------------------------------------------------------- 6. functional accuracy
print("6. functional accuracy and ancilla cleanliness")
rng = random.Random(4)
worst = 0.0
ok = True
for _ in range(40):
    th = rng.uniform(0, 2 * math.pi)
    eps = 10 ** -rng.randint(1, 5)
    r = R.synthesize_ry(th, eps)
    err = R.clean_ancilla_error(r['gates'], th)
    ok &= err <= eps
    worst = max(worst, err / eps)
check("sup_psi ||C(|psi>|00>) - (R_y|psi>)|00>|| <= eps over 40 random (theta,eps)",
      ok, "worst err/eps = %.3f" % worst)

r = R.synthesize_ry(1.234, 1e-5)
M, k = R.circuit_matrix(r['gates'])
rad = math.sqrt(2.0) ** k
leak = math.sqrt(sum((M[i][0] / rad) ** 2 + (M[i][1] / rad) ** 2 for i in range(2, 8)))
check("ancilla leakage out of |00> is O(eps)", leak <= 2e-5, "leak = %.2e" % leak)

# full state-vector simulation on a random input
psi = np.array([0.6 + 0j, -0.8])
full = np.zeros(8, dtype=complex)
full[0], full[1] = psi[0], psi[1]                      # |psi>|00>, data = LSB
Mn = np.array([[M[i][j] / rad for j in range(8)] for i in range(8)])
out = Mn @ full
th = 1.234
Ry = np.array([[math.cos(th / 2), -math.sin(th / 2)],
               [math.sin(th / 2), math.cos(th / 2)]])
want = np.zeros(8, dtype=complex)
want[0], want[1] = Ry @ psi
check("state-vector simulation matches R_y(1.234) to 1e-5",
      np.linalg.norm(out - want) < 1e-5, "||diff|| = %.2e" % np.linalg.norm(out - want))

# ------------------------------------------------------------------ 7. edge cases
print("7. edge cases")
r0 = R.synthesize_ry(0.0, 1e-6)
check("theta=0 gives the identity (0 or only trivial gates)",
      R.clean_ancilla_error(r0['gates'], 0.0) < 1e-12,
      "%d gates" % r0['counts']['total'])
ok = True
for th in [math.pi / 2, math.pi, 3 * math.pi / 2, 2 * math.pi]:
    rr = R.synthesize_ry(th, 1e-8)
    ok &= R.clean_ancilla_error(rr['gates'], th) < 1e-12 and rr['k'] <= 2
check("Clifford angles (pi/2, pi, 3pi/2, 2pi) come out exact with k<=2", ok)
ok = True
for th in [-1.0, 7.5, 100.0]:                          # out of [0,2pi)
    rr = R.synthesize_ry(th % (4 * math.pi), 1e-4)
    ok &= R.clean_ancilla_error(rr['gates'], th % (4 * math.pi)) <= 1e-4
check("angles outside [0,2pi) handled (reduce mod 4pi first)", ok)

print("\n%d passed, %d failed" % (len(PASS), len(FAIL)))
if FAIL:
    print("failures:", FAIL)
    raise SystemExit(1)
