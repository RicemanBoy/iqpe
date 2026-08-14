import qutip as qt
import gates as gate
import numpy as np
import qutip_qip.operations as op
from solovay_kitaev import solovay_kitaev
from utils import load_tree

#U = gate.R([0.21, 0.14, 0.7], 7 * np.pi / 6)

U = qt.tensor(qt.qeye(2), gate.Rz(0.25*np.pi))

tree = load_tree('SolovayKitaev/trees/HCPhase_5.pkl')

U_approx = solovay_kitaev(U, tree, 2)

print(U.full())
print(U_approx.full())
print(qt.tracedist(U, U_approx))