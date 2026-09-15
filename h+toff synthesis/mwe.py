"""Minimal working example."""
import math
import numpy as np
import rysynth as R

theta, eps = 0.4321, 1e-3

res = R.synthesize_ry(theta, eps)
print("theta = %.4f   eps = %.0e" % (theta, eps))
print("lde k = %d   gates = %d   %s" % (res['k'], res['counts']['total'],
                                        {k: v for k, v in res['counts'].items()
                                         if k != 'total'}))
print("v0 =", res['v0'])
print("v1 =", res['v1'])
print("|v0|^2 = 2^%d :" % res['k'], sum(x * x for x in res['v0']) == 1 << res['k'])
print("\nfirst 12 gates:", R.print_circuit(res['gates'], 12))

# --- verify by direct simulation -------------------------------------------
M, k = R.circuit_matrix(res['gates'])          # exact integer numerator + lde
U = np.array([[M[i][j] / math.sqrt(2) ** k for j in range(8)] for i in range(8)])
psi = np.array([0.6, -0.8])                    # data qubit state (qubit 0 = LSB)
inp = np.zeros(8); inp[0], inp[1] = psi        # |psi> (x) |00>
out = U @ inp
c, s = math.cos(theta / 2), math.sin(theta / 2)
want = np.array([c * psi[0] - s * psi[1], s * psi[0] + c * psi[1]])
print("\nout (ancilla |00> block) =", np.round(out[:2], 8))
print("R_y(theta)|psi>          =", np.round(want, 8))
print("ancilla leakage          = %.2e" % np.linalg.norm(out[2:]))
print("operator error           = %.2e  (<= eps: %s)"
      % (R.clean_ancilla_error(res['gates'], theta),
         R.clean_ancilla_error(res['gates'], theta) <= eps))
