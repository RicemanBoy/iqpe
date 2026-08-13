import qutip as qt
import qutip_qip.operations as op
import numpy as np

# 1 qubit gates
I = qt.qeye(2)
X = qt.sigmax()
Y = qt.sigmay()
Z = qt.sigmaz()
H = qt.Qobj([[1, 1],
             [1, -1]]) / np.sqrt(2)

S = qt.Qobj([
    [1, 0],
    [0, 1j]
            ])

T = qt.Qobj([
    [1, 0],
    [0, np.exp(1j * np.pi / 4)]
            ])

SQNOT = 0.5 * qt.Qobj([
    [1+1j, 1-1j],
    [1-1j, 1+1j]
])

def det(A):
    return np.linalg.det(A.full())

def Rx(angle):
    return op.rx(angle)

def Ry(angle):
    return op.ry(angle)

def Rz(angle):
    return op.rz(angle)

def R(axis, angle):
    angle = np.remainder(angle, 2 * np.pi)
    if not (type(axis) is np.ndarray):
        axis = np.array(axis)
    axis = axis / np.linalg.norm(axis)
    return qt.Qobj(np.cos(angle / 2) * I - 1j * np.sin(angle / 2) * (axis[0] * X + axis[1] * Y + axis[2] * Z))

def Phase(angle):
    return qt.phasegate(angle)

def H1():
    return qt.tensor(H, I)

def H2():
    return qt.tensor(I, H)

def Real_Rz(angle):
    """
    Creates the real version of Rz with the help of an ancilla.

    Parameters
    ----------
    angle : Angle of the Rz-gate

    Returns
    -------
    Real 4x4 matrix that applies on the data qubit AND ancilla qubit, where the data qubit is the first one.
    """
    Rz_matrix = Rz(angle).full()
    Rz_real = qt.Qobj(np.real(Rz_matrix))
    Rz_imag = qt.Qobj(np.imag(Rz_matrix))
    return qt.tensor(qt.qeye(2), Rz_real) + qt.tensor(-X*Z, Rz_imag)

# 2 qubit gates
CNOT = op.cnot()
CZ = op.csign()
Berkeley = op.berkeley()
SWAP = op.swap()
iSWAP = op.iswap()
SQSWAP = op.sqrtswap()
SQiSWAP = op.sqrtiswap()
CS = op.cs_gate()

def aSWAP(angle):
    return op.swapalpha(angle)

# 3 qubit gates
Fredkin = op.fredkin()
Toffoli = op.toffoli()



# Get unitary axis and angle
def bloch(U):
    if isinstance(U, qt.Qobj):
        U = U.full()
    angle = np.real(2 * np.arccos(np.trace(U) / 2))
    sin = np.sin(angle / 2)
    eps = 1e-10
    if sin < eps:
        axis = [0, 0, 1]
    else:
        nz = np.imag(U[1, 1] - U[0, 0]) / (2 * sin)
        nx = -np.imag(U[1, 0]) / sin
        ny = np.real(U[1, 0]) / sin
        axis = [nx, ny, nz]
    return axis, angle

# Create SU2 operator
def su2(U):
    t = 0.0j + det(U)
    phase = np.sqrt(1 / t)
    return phase * U
