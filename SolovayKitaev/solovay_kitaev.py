import gates as gate
import numpy as np
import qutip as qt
from utils import multiply, gcd, su, phase_features

# Solovay-Kitaev approximation algorithm
def solovay_kitaev(U, tree, n):
    d = tree['basis'][0].shape[0]
    if U.shape[0] != d:
        raise ValueError('U has dimension {0}, but the tree is built from {1} dimensional '
                         'gates'.format(U.shape[0], d))
    if np.asarray(tree['tree'].data).shape[1] != 2 * d ** 2:
        raise ValueError('The tree was built with an older version of create_tree, '
                         'please create it again')

    def skfunc(U, n):
        if n == 0:
            # One query per SU(d) representative, so the lookup ignores global phases
            distance, index = tree['tree'].query(phase_features(su(U)), k=1)
            name = tree['names'][index[np.argmin(distance[:, 0]), 0]]
            if name == '':
                return qt.qeye(U.dims[0])
            else:
                basis = tree['basis']
                V = multiply([basis[int(x)] for x in name])
                return V if V.dims == U.dims else qt.Qobj(V.full(), dims=U.dims)
        else:
            U_next = skfunc(U, n - 1)
            V, W = gcd(U * U_next.dag())
            V_next = skfunc(V, n - 1)
            W_next = skfunc(W, n - 1)
            return V_next * W_next * V_next.dag() * W_next.dag() * U_next

    return skfunc(U, n)