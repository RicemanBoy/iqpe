import numpy as np
import pickle as pkl
import warnings
from scipy.linalg import logm, schur
from sklearn.neighbors import KDTree
import gates as gate
import qutip as qt

# List multiplication
def multiply(arr):
    A = arr[0]
    for i in range(1, len(arr)):
        A = A * arr[i]
    return A


# Determinant
def det(A):
    return np.linalg.det(A.full())


# Diagonalization
def diagonalize(A):
    d, V = np.linalg.eig(A.full())
    return d, qt.Qobj(V)


# Traceless hermitian H with expm(1j * H) == U up to a global phase
def su_log(U):
    M = U.full() if isinstance(U, qt.Qobj) else np.asarray(U)
    d = M.shape[0]
    if np.linalg.norm(M.conj().T @ M - np.eye(d)) < 1e-10:
        # A unitary is diagonalized by its schur decomposition, which is much faster
        # (and more accurate) than the general matrix logarithm
        T, Z = schur(M, output='complex')
        H = (Z * np.angle(np.diag(T))) @ Z.conj().T
    else:
        H = -1j * logm(M)
    H = 0.5 * (H + H.conj().T)
    return H - np.trace(H) / d * np.eye(d)


# expm(1j * A) for a hermitian A
def herm_exp(A):
    w, Q = np.linalg.eigh(A)
    return (Q * np.exp(1j * w)) @ Q.conj().T


# Orthonormal basis of the traceless hermitian matrices (su(d) generators)
_basis_cache = {}


def su_basis(d):
    if d not in _basis_cache:
        T = []
        for j in range(d):
            for k in range(j + 1, d):
                M = np.zeros((d, d), dtype=complex)
                M[j, k] = M[k, j] = 1 / np.sqrt(2)
                T.append(M)
                M = np.zeros((d, d), dtype=complex)
                M[j, k], M[k, j] = 1j / np.sqrt(2), -1j / np.sqrt(2)
                T.append(M)
        for m in range(1, d):
            v = np.zeros(d)
            v[:m], v[m] = 1, -m
            T.append(np.diag(v / np.linalg.norm(v)).astype(complex))
        _basis_cache[d] = T
    return _basis_cache[d]


# Real expansion coefficients of a traceless hermitian matrix in su_basis
def su_coeffs(R, T):
    return np.array([np.real(np.trace(t.conj().T @ R)) for t in T])


def ladder_generators(H):
    """
    Hermitian A, B that solve the commutator equation 1j * (A @ B - B @ A) == H exactly.

    Both generators are built on the eigenbasis of H as a chain of sigma_x / sigma_y
    rotations on neighbouring eigenvectors, with amplitudes fixed by the partial sums of
    the (descending) eigenvalues. The cross terms between neighbouring links cancel, so
    the identity holds for every dimension and ||A|| = ||B|| = O(sqrt(||H||)).

    Parameters
    ----------
    H : Traceless hermitian d x d array

    Returns
    -------
    Tuple (A, B) of hermitian d x d arrays.
    """
    d = H.shape[0]
    w, Q = np.linalg.eigh(H)
    order = np.argsort(w)[::-1]
    w, Q = w[order], Q[:, order]
    # Partial sums of the descending eigenvalues are non-negative for a traceless H
    alpha = np.sqrt(np.clip(np.cumsum(w)[:-1], 0, None) / 2)
    A = np.zeros((d, d), dtype=complex)
    B = np.zeros((d, d), dtype=complex)
    for k in range(d - 1):
        A[k, k + 1] = A[k + 1, k] = alpha[k]
        B[k, k + 1], B[k + 1, k] = 1j * alpha[k], -1j * alpha[k]
    return Q @ A @ Q.conj().T, Q @ B @ Q.conj().T


# Logarithm of the group commutator of expm(1j * A) and expm(1j * B)
def group_commutator_log(A, B):
    V, W = herm_exp(A), herm_exp(B)
    return su_log(V @ W @ V.conj().T @ W.conj().T), V, W


# Group Commutator Decomposition in arbitrary dimension
def gcd_general(U, tol=1e-12, max_iter=50, h=1e-6):
    """
    Balanced group commutator decomposition of a d x d unitary.

    The ladder construction solves the decomposition to leading order in the Lie algebra;
    the remaining higher order BCH terms are removed by a Gauss-Newton iteration on the
    exact group equation V @ W @ V.dag() @ W.dag() == U. The generators keep their
    O(sqrt(||log U||)) size, which is what the Solovay-Kitaev recursion needs.

    Parameters
    ----------
    U : Unitary Qobj of dimension d, close to the identity
    tol : Residual (in the algebra) at which the iteration stops
    max_iter : Maximum number of Gauss-Newton steps
    h : Step size of the finite difference jacobian

    Returns
    -------
    Tuple (V, W) of Qobjs with V * W * V.dag() * W.dag() == U up to a global phase.
    """
    H = su_log(U)
    d = H.shape[0]
    T = su_basis(d)
    n = len(T)

    A, B = ladder_generators(H)
    R, V, W = group_commutator_log(A, B)
    r = su_coeffs(R - H, T)
    res = np.linalg.norm(r)

    for _ in range(max_iter):
        if res < tol:
            break
        # Finite difference jacobian of the group commutator w.r.t. A and B
        J = np.zeros((n, 2 * n))
        for j, t in enumerate(T):
            Rp, _, _ = group_commutator_log(A + h * t, B)
            Rm, _, _ = group_commutator_log(A - h * t, B)
            J[:, j] = su_coeffs((Rp - Rm) / (2 * h), T)
            Rp, _, _ = group_commutator_log(A, B + h * t)
            Rm, _, _ = group_commutator_log(A, B - h * t)
            J[:, n + j] = su_coeffs((Rp - Rm) / (2 * h), T)
        # Minimum norm solution keeps the generators close to the leading order guess
        delta = np.linalg.lstsq(J, -r, rcond=None)[0]
        dA = sum(delta[j] * T[j] for j in range(n))
        dB = sum(delta[n + j] * T[j] for j in range(n))
        # Backtracking, so a step never makes the residual worse
        step = 1.0
        for _ in range(20):
            R_new, V_new, W_new = group_commutator_log(A + step * dA, B + step * dB)
            r_new = su_coeffs(R_new - H, T)
            if np.linalg.norm(r_new) < res:
                A, B, V, W = A + step * dA, B + step * dB, V_new, W_new
                r, res = r_new, np.linalg.norm(r_new)
                break
            step /= 2
        else:
            break

    if res > 1e-8:
        warnings.warn('gcd: group commutator decomposition converged only to {0:.2e}'.format(res))

    dims = U.dims if isinstance(U, qt.Qobj) else None
    return qt.Qobj(V, dims=dims), qt.Qobj(W, dims=dims)


# Group Commutator Decomposition
def gcd(U):
    if U.shape[0] != 2:
        return gcd_general(U)

    # The formulas below need det(U) == 1, physical gates like H or T do not have it
    U = qt.Qobj(su(U))

    theta = 2 * np.arccos(np.real(U.tr() / 2))
    phi = 2 * np.arcsin(np.sqrt(np.sqrt((0.5 - 0.5 * np.cos(theta / 2)))))

    axis, angle = gate.bloch(U)
    V = gate.Rx(phi)
    if axis[2] < 0:
        W = gate.Ry(2 * np.pi - phi)
    else:
        W = gate.Ry(phi)

    _, V1 = diagonalize(U)
    _, V2 = diagonalize(V * W * V.dag() * W.dag())
    S = V1 * V2.dag()
    V_tilde = S * V * S.dag()
    W_tilde = S * W * S.dag()
    return V_tilde, W_tilde


# Representative of U in SU(d), i.e. U with its global phase removed
def su(U):
    M = U.full() if isinstance(U, qt.Qobj) else np.asarray(U)
    return M * np.linalg.det(M) ** (-1 / M.shape[0])


# Feature vector of an SU(d) matrix, so that the euclidean distance of two feature
# vectors is the frobenius distance of the two matrices
def features(M):
    return np.concatenate([M.real.ravel(), M.imag.ravel()])


def phase_features(M):
    """
    Feature vectors of all SU(d) representatives of a unitary.

    An operator only fixes its SU(d) representative up to a d-th root of unity, so a
    nearest neighbour search has to compare against all d of them. Querying the tree with
    every row and keeping the closest hit makes the lookup independent of global phases.

    Parameters
    ----------
    M : SU(d) matrix as returned by su()

    Returns
    -------
    Array of shape (d, 2 * d ** 2).
    """
    d = M.shape[0]
    return np.array([features(np.exp(2j * np.pi * k / d) * M) for k in range(d)])


# Hashable key that is equal exactly for unitaries that agree up to a global phase
def phase_key(M, decimals=9):
    return min((np.round(v, decimals) + 0.0).tobytes() for v in phase_features(M))


# Operator norm distance of two gates, minimized over the global phase
def dist(U, V):
    A = U.full() if isinstance(U, qt.Qobj) else np.asarray(U)
    B = V.full() if isinstance(V, qt.Qobj) else np.asarray(V)
    t = np.trace(A.conj().T @ B)
    phase = np.conj(t) / abs(t) if abs(t) > 1e-14 else 1
    return np.linalg.norm(A - phase * B, 2)


# Creating tree
def create_tree(basis, max_depth=10):
    """
    Nearest neighbour tree of all gate sequences up to a given length.

    Works for a basis of any dimension (2, 4, 8, ... for one, two, three qubit gates).
    The sequences are enumerated breadth first and a sequence is only extended if its
    product is an operator that no shorter sequence produces already, which prunes the
    identities and all other duplicates. Note that the number of surviving sequences still
    grows like len(basis) ** max_depth once the basis stops producing collisions.

    Parameters
    ----------
    basis : List of unitary Qobjs, all of the same dimension
    max_depth : Maximum number of gates in a sequence

    Returns
    -------
    Dict with the KDTree, the sequence names and the basis.
    """
    n = len(basis)
    d = basis[0].shape[0]
    arrays = [su(b) for b in basis]

    identity = np.eye(d, dtype=complex)
    X = [features(identity)]
    names = ['']
    seen = {phase_key(identity)}
    frontier = [('', identity)]

    for i in range(1, max_depth + 1):
        print('Creating tree: {0}/{1}'.format(i, max_depth))
        next_frontier = []
        for name, U in frontier:
            for k in range(n):
                M = su(U @ arrays[k])
                key = phase_key(M)
                if key in seen:
                    continue
                seen.add(key)
                names.append(name + str(k))
                X.append(features(M))
                next_frontier.append((name + str(k), M))
        frontier = next_frontier
        if not frontier:
            break
    print('Creating tree')
    tree = KDTree(np.array(X), metric='euclidean')
    return {'tree': tree, 'names': names, 'basis': basis}


def save_tree(tree, filename):
    with open(filename, 'wb') as f:
        pkl.dump(tree, f, pkl.HIGHEST_PROTOCOL)


def load_tree(filename):
    with open(filename, 'rb') as f:
        tree = pkl.load(f)
    return tree