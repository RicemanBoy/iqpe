from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister
import numpy as np
import matplotlib.pyplot as plt
#import bitstring
from qiskit_aer import AerSimulator
from qiskit.quantum_info import Operator

from qiskit import transpile
from qiskit.quantum_info import Statevector
from qiskit.primitives import StatevectorSampler, StatevectorEstimator
from qiskit.quantum_info import SparsePauliOp
from qiskit.circuit.classical import expr

from qiskit_aer.noise import (NoiseModel, QuantumError, ReadoutError,
    pauli_error, depolarizing_error, thermal_relaxation_error)

from qiskit.circuit.library import UnitaryGate

matrix_h = ([[2**(-0.5),2**(-0.5)],[2**(-0.5),-2**(-0.5)]])
h_ideal = UnitaryGate(matrix_h)

matrix_cx = ([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
cx_ideal = UnitaryGate(matrix_cx)       #Erst Target, dann Control Qubit!!

matrix_x = ([[0,1],[1,0]])
x_ideal = UnitaryGate(matrix_x)

matrix_z = ([[1,0],[0,-1]])
z_ideal = UnitaryGate(matrix_z)

matrix_t = ([[1,0],[0,np.exp(1j*(np.pi/4))]])
t_ideal = UnitaryGate(matrix_t)

matrix_tdg = ([[1,0],[0,np.exp(-1j*(np.pi/4))]])
tdg_ideal = UnitaryGate(matrix_tdg)

def parity_values(n):
    return [
        value
        for value in range(2**n)
        if value.bit_count() % 2 == 1
    ]

def gates(qc:QuantumCircuit):
    hmm = dict(qc.count_ops())
    hmm["reset"] = 0
    hmm["measure"] = 00
    hmm["swap"] = 0
    return sum(hmm.values())

def convert(bin: str):                  #konvertiert den bitstring in decimal, e.g. 0110 = 0.375
    k = list(bin)
    a = [int(i) for i in k]
    n = 0
    for i in range(len(a)):
        if a[i] == 1:
            n += 1/2**(i+1)
    return n

def avg7_ramsey(code: str, iter: int, noise: float, qec = False, k = 1, bias = 0, post = False, path = ""):       #only exact angles!  
    n = 15
    angle = np.linspace(0,1,n+2)
    angle = np.delete(angle, [n+1])
    angle = np.delete(angle, [0])

    a, b = [], []
    with open("{}unitary{}_improved_with_rs.txt".format(path, n), "r") as file:
        for line in file:
            a.append(list(map(str, line.strip().split(","))))
    with open("{}adjunitary{}_improved_with_rs.txt".format(path, n), "r") as file:
        for line in file:
            b.append(list(map(str, line.strip().split(","))))
    
    y = 0
    y_list, bruh1 = [], []
    for m in range(k):
        for o in range(7):
            gatecount = 0
            bitstring = ""
            rots = []
            for t in range(iter):
                rots = [k*0.5 for k in rots]
                counter = 0
                while True:
                    if code == "rotsurf":
                        self = RotSurf9q(1, noise=noise)

                    self.err = qec
                    self.h(pos=0)
                    #############################
                    for j in range(2**(iter-t-1)):
                        self.cu_ramsey(a[2*o+1])
                    ###############################
                    for l in rots:
                        if l == 0.25:
                            self.sdg(pos=0)
                        if l == 0.125:
                            self.tdg(pos=0)
                    self.h(pos=0)
                    if self.err:
                        self.flagFTec(0)

                    self.readout(pos=0, shots=1, p = noise)
                    gatecount += gates(self.qc)
            
                    if self.zeros == 1:
                        bitstring += "0"
                        break
                    if self.ones == 1:
                        bitstring += "1"
                        rots.append(0.5)
                        break
                    if not post:
                        if np.random.rand() < 0.5:
                            bitstring += "0"
                            break
                        else:
                            bitstring += "1"
                            rots.append(0.5)
                            break
                    counter += 1
                    print("Angle {}, {}%% errorrate, Iteration {}: {} Repetition".format(2*o+1, noise*100, t, counter))
                    del self
            bitstring = bitstring[::-1]
            hmm = convert(bitstring)
            diff = min(np.abs(hmm-angle[2*o+1]), 1-np.abs(hmm-angle[2*o+1]))
            y += diff
            print("Performance {}for angle {} ({} gates): ".format("(QEC) " if qec else "", 2*o+1, gatecount), diff)
            bruh1.append(diff), y_list.append(diff)
    y = y/(7*k)
    arg = 0
    for i in range(len(bruh1)):
        arg += (y-bruh1[i])**2
    sigma = ((1/(k*n))*arg)**0.5
    sigma = sigma/((k*n)**0.5)

    return y_list

# Rotated surface [[9,1,3]] code, data qubits on a 3x3 grid    0 1 2 / 3 4 5 / 6 7 8
# A logical Hadamard rotates the lattice by 90 degrees, (row,col) -> (2-col,row), which is exactly the
# permutation readout() already applies. That rotation maps STABS_A_9Q[i] onto STABS_B_9Q[i] and swaps
# the two stabilizer types, so the supports below stay fixed and only the measured Pauli type follows
# self.hadamards[pos]:
#   hadamards even -> set A is measured as X-Stabilizers (detects Z errors), set B as Z-Stabilizers
#   hadamards odd  -> set A is measured as Z-Stabilizers (detects X errors), set B as X-Stabilizers
STABS_A_9Q = ((0,1), (1,2,4,5), (3,4,6,7), (7,8))       #order matches qecc[3],qecc[2],qecc[1],qecc[0]
STABS_B_9Q = ((3,6), (0,1,3,4), (4,5,7,8), (2,5))       #order matches qecc[7],qecc[6],qecc[5],qecc[4]

# E_min(s) per stabilizer set. Key = the 4 bit sub-syndrome of that set, read MSB..LSB, i.e. in the
# order the stabilizers are listed above. Value = a minimum weight error support producing it.
CORR_A_9Q = {
    "0000": (),    "0001": (8,),   "0010": (3,),   "0011": (7,),
    "0100": (2,),  "0101": (2,8),  "0110": (4,),   "0111": (2,7),
    "1000": (0,),  "1001": (0,8),  "1010": (0,3),  "1011": (0,7),
    "1100": (1,),  "1101": (1,8),  "1110": (0,4),  "1111": (1,7),
}
CORR_B_9Q = {
    "0000": (),    "0001": (2,),   "0010": (7,),   "0011": (5,),
    "0100": (0,),  "0101": (0,2),  "0110": (4,),   "0111": (0,5),
    "1000": (6,),  "1001": (2,6),  "1010": (3,4),  "1011": (5,6),
    "1100": (3,),  "1101": (2,3),  "1110": (3,7),  "1111": (3,5),
}

# Flag error sets E(g_i) of the 1-flag circuits in flagsyndrome(), needed for case 3 of the Flag 1-FTEC
# Protocol. Key = the flag register bit, which is aligned with the stabilizer index: flags[0..3] are the
# set A circuits, flags[4..7] the set B circuits. Only the weight 4 stabilizers carry a flag ancilla, the
# weight 2 ones are fault tolerant on their own. Value = the 8 bit qecc slice of the following unflagged
# syndrome measurement mapped onto the qubits the flagged fault hit.
# Propagating every single fault through the flag circuits gives E(g) = {I, {d1}, {d3,d4}, {d4}} for a
# CNOT order [flag,d1,d2,d3,flag,d4]. For set A all of those already agree with E_min(s) up to a
# stabilizer, so they fall through to CORR_A_9Q/CORR_B_9Q as in the "otherwise apply E_min(s)" of the
# protocol. For set B the weight 2 element does not: E_min would return a weight 1 error and leave a
# logical Z (resp. X) behind, so those two get an entry.
CORR_FLAG_9Q = {
    1: {},                                  #C(A[1]) = {1,2,4,5}
    2: {},                                  #C(A[2]) = {3,4,6,7}
    5: {"00000100": (0,1)},                 #C(B[1]) = {0,1,3,4}, E_min would give (2,) -> logical
    6: {"00000010": (7,8)},                 #C(B[2]) = {4,5,7,8}, E_min would give (3,) -> logical
}

class RotSurf9q:
    def __init__(self, n: int, magic = 0, noise = 0):
        self.n = n
        self.magic = magic

        self.zeros = 0
        self.ones = 0
        self.preselected = 0
        self.post = 0

        self.err = False
        self.postselection = True
        self.classical_ec = False
        self.qec_counter = 0

        self.hadamards = [0 for i in range(n+magic)]

        self.noise_model = self.__noise_model__(noise, 0)

        qr = QuantumRegister(9*(n+magic)+4, "q")
        cbit = ClassicalRegister(4,"c")
        self.cbit = cbit
        self.qc = QuantumCircuit(qr,cbit)
        for i in range(9*n):
            self.qc.id(i)
        for i in range(n):
            self.qc.h(9*i+1)
            self.qc.h(9*i+3)
            self.qc.h(9*i+5)
            self.qc.h(9*i+7)

            self.qc.cx(9*i+1,9*i)
            self.qc.cx(9*i+5,9*i+4)
            self.qc.cx(9*i+7,9*i+8)

            self.qc.cx(9*i+5,9*i+2)

            self.qc.cx(9*i+3,9*i+4)
            self.qc.cx(9*i+2,9*i+1)

            self.qc.cx(9*i+3,9*i+6)

            self.qc.cx(9*i+6,9*i+7)

        self.qecc = ClassicalRegister(8)
        self.qc.add_register(self.qecc)

    def __noise_model__(self, p: float, bias: float):
        p_x, p_z = 0, 0
        if bias > 0:
            p_x += (bias/(1+bias))*p
            p_z += p - p_x
        elif bias < 0:
            p_z += (np.abs(bias)/(1+np.abs(bias)))*p
            p_x += p - p_z
        else:
            p_x += p/2
            p_z += p/2
        noise_model = NoiseModel()
        p_error = pauli_error([["X",p_x],["I",1-p],["Z",p_z]])
        p_error_2 = pauli_error([["XI",p_x/2],["IX",p_x/2],["II",1-p],["ZI",p_z/2],["IZ",p_z/2]])
        p_error_3 = pauli_error([["XII",p_x/3],["IXI",p_x/3],["IIX",p_x/3],["III",1-p],["ZII",p_z/3],["IZI",p_z/3],["IIZ",p_z/3]])
        noise_model.add_all_qubit_quantum_error(p_error, ['x', "z", 'h', "s", "sdg", "t", "tdg", 'id',"rx"])  # Apply to single-qubit gates
        noise_model.add_all_qubit_quantum_error(p_error_2, ['cx'])  # Apply to 2-qubit gates
        noise_model.add_all_qubit_quantum_error(p_error_3, ['ccx'])  # Apply to 3-qubit gates
        return noise_model

    def testing_magic(self, target: int):
        self.h(pos=target)
        self.sdg_cheat(pos=target)
        self.h(pos=target)

        pos = target + 1        #pos = welcher qubit als magic state verwendet wird, target = qubit auf den das gate teleportiert werden soll
        for i in range(9):
            self.qc.reset(9*pos+i)
        
        self.qc.h(9*pos+1)
        self.qc.h(9*pos+3)

        # self.qc.h(9*pos+4)
        # self.qc.t(9*pos+4)

        self.qc.ry(np.pi/4,9*pos+4)

        self.qc.h(9*pos+5)
        self.qc.h(9*pos+7)

        self.qc.cx(9*pos+4,9*pos+0)
        # self.qc.cx(9*pos+4,9*pos+2)
        # self.qc.cx(9*pos+4,9*pos+6)
        self.qc.cx(9*pos+4,9*pos+8)

        self.qc.cx(9*pos+1,9*pos+0)
        self.qc.cx(9*pos+7,9*pos+8)

        self.qc.cx(9*pos+5,9*pos+2)
        self.qc.cx(9*pos+5,9*pos+4)

        self.qc.cx(9*pos+3,9*pos+4)
        self.qc.cx(9*pos+3,9*pos+6)

        self.qc.cx(9*pos+5,9*pos+1)
        self.qc.cx(9*pos+3,9*pos+7)

        ################ Controlled Hadamard injection#################################
        anc = self.qc.num_qubits - 1
        ancc = anc - 1
        flag = ancc - 1
        self.qc.reset(anc), self.qc.reset(ancc)

        self.qc.h(ancc)

        if self.hadamards[pos]%2==0:
            self.qc.cx(3+9*pos, anc)
            self.qc.cx(4+9*pos, anc)
            self.qc.cx(5+9*pos, anc)
        else:
            self.qc.cx(1+9*pos, anc)
            self.qc.cx(4+9*pos, anc)
            self.qc.cx(7+9*pos, anc)

        self.qc.h(anc)

        if self.hadamards[pos]%2==0:
            self.qc.cx(anc, 1+9*pos)
            self.qc.cx(anc, 4+9*pos)
            self.qc.cx(anc, 7+9*pos)
        else:
            self.qc.cx(anc, 3+9*pos)
            self.qc.cx(anc, 4+9*pos)
            self.qc.cx(anc, 5+9*pos)

        self.qc.h(anc)

        self.qc.x(ancc)
        self.qc.ry(-np.pi/4, anc)
        self.qc.cz(ancc, anc)
        self.qc.ry(np.pi/4, anc)
        self.qc.x(ancc)

        self.qc.measure(anc,0)         #anc
        if self.hadamards[pos]%2==0:
            with self.qc.if_test((0,0)): #anc
                self.qc.cx(ancc, 9*pos+1)
                self.qc.cx(ancc, 9*pos+4)
                self.qc.cx(ancc, 9*pos+7)
        else:
            with self.qc.if_test((0,0)):
                self.qc.cx(ancc, 9*pos+3)
                self.qc.cx(ancc, 9*pos+4)
                self.qc.cx(ancc, 9*pos+5)

        if self.hadamards[pos]%2==0:
            with self.qc.if_test((0,1)): #anc
                self.qc.x(9*pos+1)
                self.qc.x(9*pos+4)
                self.qc.x(9*pos+7)

                self.qc.z(9*pos+3)
                self.qc.z(9*pos+4)
                self.qc.z(9*pos+5)

                self.qc.cx(ancc, 9*pos+1)
                self.qc.cx(ancc, 9*pos+4)
                self.qc.cx(ancc, 9*pos+7)

                self.qc.z(ancc)
        else:
            with self.qc.if_test((0,1)):
                self.qc.x(9*pos+3)
                self.qc.x(9*pos+4)
                self.qc.x(9*pos+5)

                self.qc.z(9*pos+1)
                self.qc.z(9*pos+4)
                self.qc.z(9*pos+7)

                self.qc.cx(ancc, 9*pos+3)
                self.qc.cx(ancc, 9*pos+4)
                self.qc.cx(ancc, 9*pos+5)

                self.qc.z(ancc)
        
        self.qc.h(ancc)
        self.qc.measure(ancc, 0)
        ################ qec cylce #################################
        #empty for now
        ############### Controlled Y-gate###########################
        self.sdg_cheat(pos=target)
        for i in range(9):
            self.qc.cx(i+9*pos,i+9*target)
        self.s_cheat(pos=target)
        ###### Measure logical state of the magic state for state injection##########
        self.qc.reset(anc)
        self.sdg_cheat(pos=pos)
        self.h(pos=pos)
        for i in range(9):
            self.qc.cx(i+pos*9, anc)
        self.qc.measure(anc,0)
        ###############Apply conditioned Ry(pi/2) onto the Target####################
        for i in range(9):
            with self.qc.if_test((0,1)):
                self.qc.h(i+9*target)             

        self.qc.reset(ancc)                      #Need this to track hadamard rots on ancilla because of conditional hadamard above
        if self.hadamards[target]%2 != 0:
            self.qc.x(ancc)
        with self.qc.if_test((0,1)):
            self.qc.x(ancc)
        self.qc.measure(ancc, 1)

        with self.qc.if_test((1,0)):
            with self.qc.if_test((0,1)):
                self.qc.x(1+9*target)
                self.qc.x(4+9*target)
                self.qc.x(7+9*target)

        with self.qc.if_test((1,1)):
            with self.qc.if_test((0,1)):
                self.qc.x(3+9*target)
                self.qc.x(4+9*target)
                self.qc.x(5+9*target)
        self.h(pos=target)
        self.s_cheat(pos=target)
        self.h(pos=target)
        if self.err:
            self.qec(pos=target)
    
    def id(self, pos: int):
        for i in range(9):
            self.qc.id(9*pos + i)

    def x(self, pos: int):
        if self.hadamards[pos]%2==0:
            self.qc.x(9*pos+1)
            self.qc.x(9*pos+4)
            self.qc.x(9*pos+7)
        else:
            self.qc.x(9*pos+3)
            self.qc.x(9*pos+4)
            self.qc.x(9*pos+5)
    
    def z(self, pos: int):
        if self.hadamards[pos]%2==0:
            self.qc.z(9*pos+3)
            self.qc.z(9*pos+4)
            self.qc.z(9*pos+5)
        else:
            self.qc.z(9*pos+1)
            self.qc.z(9*pos+4)
            self.qc.z(9*pos+7)

    def h(self, pos: int):
        for i in range(9):
            self.qc.h(9*pos+i)
        self.hadamards[pos] += 1

    def cnot(self, control: int, target: int):
        if control == 0:
            if self.hadamards[0]%2==1 and self.hadamards[1]%2==0:
                self.qc.cx(0,9+6)
                self.qc.cx(1,9+3)
                self.qc.cx(2,9+0)
                self.qc.cx(3,9+7)
                self.qc.cx(4,9+4)
                self.qc.cx(5,9+1)
                self.qc.cx(6,9+8)
                self.qc.cx(7,9+5)
                self.qc.cx(8,9+2)
            elif self.hadamards[0]%2==0 and self.hadamards[1]%2==1:
                self.qc.cx(0,9+2)
                self.qc.cx(1,9+5)
                self.qc.cx(2,9+8)
                self.qc.cx(3,9+1)
                self.qc.cx(4,9+4)
                self.qc.cx(5,9+7)
                self.qc.cx(6,9+0)
                self.qc.cx(7,9+3)
                self.qc.cx(8,9+6)
            else:
                for i in range(9):
                    self.qc.cx(i,9+i)
        elif control == 1:
            if self.hadamards[0]%2==0 and self.hadamards[1]%2==1:
                self.qc.cx(9+0,6)
                self.qc.cx(9+1,3)
                self.qc.cx(9+2,0)
                self.qc.cx(9+3,7)
                self.qc.cx(9+4,4)
                self.qc.cx(9+5,1)
                self.qc.cx(9+6,8)
                self.qc.cx(9+7,5)
                self.qc.cx(9+8,2)
            elif self.hadamards[0]%2==1 and self.hadamards[1]%2==0:
                self.qc.cx(9+0,2)
                self.qc.cx(9+1,5)
                self.qc.cx(9+2,8)
                self.qc.cx(9+3,1)
                self.qc.cx(9+4,4)
                self.qc.cx(9+5,7)
                self.qc.cx(9+6,0)
                self.qc.cx(9+7,3)
                self.qc.cx(9+8,6)
            else: 
                for i in range(9):
                    self.qc.cx(9+i,i)

    def cz(self):
        self.h(pos=1)
        self.cnot(control = 0, target=1)
        self.h(pos=1)

    def s(self, pos: int):
        anc = self.qc.num_qubits - 1
        self.qc.reset(anc)

        self.qc.h(anc)
        #self.qc.append(h_ideal,[anc])
        self.qc.s(anc)

        if self.hadamards[pos]%2 == 0:
            self.qc.cx(3+9*pos, anc)
            self.qc.cx(4+9*pos, anc)
            self.qc.cx(5+9*pos, anc)        
        else:
            self.qc.cx(1+9*pos, anc)
            self.qc.cx(4+9*pos, anc)
            self.qc.cx(7+9*pos, anc)  

        self.qc.measure(anc, 0)

        if self.hadamards[pos]%2 == 0:
            with self.qc.if_test((0,1)):
                self.qc.z(3+9*pos)
                self.qc.z(4+9*pos)
                self.qc.z(5+9*pos)
        else:
            with self.qc.if_test((0,1)):
                self.qc.z(1+9*pos)
                self.qc.z(4+9*pos)
                self.qc.z(7+9*pos)
   
    def s_cheat(self, pos: int):
        if self.hadamards[pos]%2==0:
            self.qc.cx(9*pos+3, 9*pos+5)
            self.qc.cx(9*pos+4, 9*pos+5)
            self.qc.s(9*pos+5)
            self.qc.cx(9*pos+3, 9*pos+5)
            self.qc.cx(9*pos+4, 9*pos+5)
        else:
            self.qc.cx(9*pos+1, 9*pos+7)
            self.qc.cx(9*pos+4, 9*pos+7)
            self.qc.s(9*pos+7)
            self.qc.cx(9*pos+1, 9*pos+7)
            self.qc.cx(9*pos+4, 9*pos+7)

    def ft_s(self, pos: int):
        anc = self.qc.num_qubits - 1
        ancc = anc - 1
        self.qc.reset(anc), self.qc.reset(ancc)
        magic = int((self.qc.num_qubits - 2)/9)-1
        #################Magic state initialization#################
        for j in range(9):
            self.qc.reset(9*magic+j)
        self.qc.h(9*magic+1)
        self.qc.h(9*magic+3)
        self.qc.h(9*magic+5)
        self.qc.h(9*magic+7)

        self.qc.cx(9*magic+1,9*magic)
        self.qc.cx(9*magic+5,9*magic+4)
        self.qc.cx(9*magic+7,9*magic+8)

        self.qc.cx(9*magic+5,9*magic+2)

        self.qc.cx(9*magic+3,9*magic+4)
        self.qc.cx(9*magic+2,9*magic+1)

        self.qc.cx(9*magic+3,9*magic+6)

        self.qc.cx(9*magic+6,9*magic+7)
        
        self.h(magic)
        self.s_cheat(magic)
        ################# Entanglement with Ancilla, check +-1 eigenstate of pauli y #################
        self.qc.h(anc)

        self.qc.cz(anc, 9*magic+1)
        #self.qc.cx(anc, ancc)
        self.qc.cz(anc, 9*magic+4)
        self.qc.cz(anc, 9*magic+7)

        self.qc.cx(anc, 9*magic+3)
        self.qc.cx(anc, 9*magic+4)
        #self.qc.cx(anc, ancc)
        self.qc.cx(anc, 9*magic+5)

        self.qc.h(anc)

        self.qc.measure(anc, 0)             #für noisefree case
        # self.qc.measure(ancc, 0)
        ################# Entanglement with target logical qubit #################
        if self.hadamards[pos]%2 == 0:
            self.qc.cx(0+9*pos,9*magic+2)
            self.qc.cx(1+9*pos,9*magic+5)
            self.qc.cx(2+9*pos,9*magic+8)
            self.qc.cx(3+9*pos,9*magic+1)
            self.qc.cx(4+9*pos,9*magic+4)
            self.qc.cx(5+9*pos,9*magic+7)
            self.qc.cx(6+9*pos,9*magic+0)
            self.qc.cx(7+9*pos,9*magic+3)
            self.qc.cx(8+9*pos,9*magic+6)
        else:    
            self.qc.cx(0+9*pos,9*magic+0)
            self.qc.cx(1+9*pos,9*magic+1)
            self.qc.cx(2+9*pos,9*magic+2)
            self.qc.cx(3+9*pos,9*magic+3)
            self.qc.cx(4+9*pos,9*magic+4)
            self.qc.cx(5+9*pos,9*magic+5)
            self.qc.cx(6+9*pos,9*magic+6)
            self.qc.cx(7+9*pos,9*magic+7)
            self.qc.cx(8+9*pos,9*magic+8)

        ################ Readout of Magic State #################
        self.qc.reset(anc)
        for j in range(9):
            self.qc.cx(9*magic+j, anc)
        ################# Z-rotation based on readout of magic state #################
        self.qc.measure(anc, 0)
        if self.hadamards[pos]%2 == 0:
            with self.qc.if_test((0,1)):
                self.qc.z(9*pos+3)
                self.qc.z(9*pos+4)
                self.qc.z(9*pos+5)
        else:
            with self.qc.if_test((0,1)):
                self.qc.z(9*pos+1)
                self.qc.z(9*pos+4)
                self.qc.z(9*pos+7)
        self.qc.reset(anc), self.qc.reset(ancc)

    def sdg(self, pos: int):
        anc = self.qc.num_qubits - 1
        self.qc.reset(anc)

        self.qc.h(anc)
        self.qc.sdg(anc)

        if self.hadamards[pos]%2 == 0:
            self.qc.cx(3+9*pos, anc)
            self.qc.cx(4+9*pos, anc)
            self.qc.cx(5+9*pos, anc)        
        else:
            self.qc.cx(1+9*pos, anc)
            self.qc.cx(4+9*pos, anc)
            self.qc.cx(7+9*pos, anc)  

        self.qc.measure(anc, 0)

        if self.hadamards[pos]%2 == 0:
            with self.qc.if_test((0,1)):
                self.qc.z(3+9*pos)
                self.qc.z(4+9*pos)
                self.qc.z(5+9*pos)
        else:
            with self.qc.if_test((0,1)):
                self.qc.z(1+9*pos)
                self.qc.z(4+9*pos)
                self.qc.z(7+9*pos)
    
    def sdg_cheat(self, pos: int):
        if self.hadamards[pos]%2==0:
            self.qc.cx(9*pos+3, 9*pos+5)
            self.qc.cx(9*pos+4, 9*pos+5)
            self.qc.sdg(9*pos+5)
            self.qc.cx(9*pos+3, 9*pos+5)
            self.qc.cx(9*pos+4, 9*pos+5)
        else:
            self.qc.cx(9*pos+1, 9*pos+7)
            self.qc.cx(9*pos+4, 9*pos+7)
            self.qc.sdg(9*pos+7)
            self.qc.cx(9*pos+1, 9*pos+7)
            self.qc.cx(9*pos+4, 9*pos+7)

    def t(self, pos: int):
        anc = self.qc.num_qubits - 1
        self.qc.reset(anc)

        self.qc.h(anc)
        #self.qc.append(h_ideal,[anc])
        self.qc.t(anc)

        if self.hadamards[pos]%2 == 0:
            self.qc.cx(3+9*pos, anc)
            self.qc.cx(4+9*pos, anc)
            self.qc.cx(5+9*pos, anc)        
        else:
            self.qc.cx(1+9*pos, anc)
            self.qc.cx(4+9*pos, anc)
            self.qc.cx(7+9*pos, anc)  

        self.qc.measure(anc, 0)
        # if z_stab:
        #     z_qec_ideal(qc, had=had, pos=pos)

        if self.hadamards[pos]%2 == 0:
            with self.qc.if_test((0,1)):
                self.qc.reset(anc)
                self.qc.h(anc)
                self.qc.s(anc)
                self.qc.cx(3+9*pos, anc)
                self.qc.cx(4+9*pos, anc)
                self.qc.cx(5+9*pos, anc) 
                self.qc.measure(anc, 0)
                with self.qc.if_test((0,1)):
                    self.qc.z(3+9*pos)
                    self.qc.z(4+9*pos)
                    self.qc.z(5+9*pos)
        else:
            with self.qc.if_test((0,1)):
                self.qc.reset(anc)
                self.qc.h(anc)
                self.qc.s(anc)
                self.qc.cx(1+9*pos, anc)
                self.qc.cx(4+9*pos, anc)
                self.qc.cx(7+9*pos, anc) 
                self.qc.measure(anc, 0)
                with self.qc.if_test((0,1)):
                    self.qc.z(1+9*pos)
                    self.qc.z(4+9*pos)
                    self.qc.z(7+9*pos)

    def t_cheat(self, pos: int):
            if self.hadamards[pos]%2==0:
                self.qc.cx(9*pos+3, 9*pos+5)
                self.qc.cx(9*pos+4, 9*pos+5)
                self.qc.t(9*pos+5)
                self.qc.cx(9*pos+3, 9*pos+5)
                self.qc.cx(9*pos+4, 9*pos+5)
            else:
                self.qc.cx(9*pos+1, 9*pos+7)
                self.qc.cx(9*pos+4, 9*pos+7)
                self.qc.t(9*pos+7)
                self.qc.cx(9*pos+1, 9*pos+7)
                self.qc.cx(9*pos+4, 9*pos+7)

    def rz_cheat(self, angle: float, pos: int):
        if self.hadamards[pos]%2==0:
            self.qc.cx(9*pos+3, 9*pos+5)
            self.qc.cx(9*pos+4, 9*pos+5)
            self.qc.rz(angle, 9*pos+5)
            self.qc.cx(9*pos+3, 9*pos+5)
            self.qc.cx(9*pos+4, 9*pos+5)
        else:
            self.qc.cx(9*pos+1, 9*pos+7)
            self.qc.cx(9*pos+4, 9*pos+7)
            self.qc.rz(angle, 9*pos+7)
            self.qc.cx(9*pos+1, 9*pos+7)
            self.qc.cx(9*pos+4, 9*pos+7)

    def rx_cheat(self, angle: float, pos: int):
        if self.hadamards[pos]%2==0:
            self.qc.cx(9*pos+7, 9*pos+1)
            self.qc.cx(9*pos+7, 9*pos+4)
            self.qc.rx(angle, 9*pos+7)
            self.qc.cx(9*pos+7, 9*pos+1)
            self.qc.cx(9*pos+7, 9*pos+4)
        else:
            self.qc.cx(9*pos+5, 9*pos+3)
            self.qc.cx(9*pos+5, 9*pos+4)
            self.qc.rx(angle, 9*pos+5)
            self.qc.cx(9*pos+5, 9*pos+3)
            self.qc.cx(9*pos+5, 9*pos+4)

    def ry_cheat(self, angle: float, pos: int):
        if self.hadamards[pos]%2==0:
            self.qc.h(9*pos+3), self.qc.h(9*pos+4), self.qc.h(9*pos+5)
            self.qc.s(9*pos+3), self.qc.s(9*pos+4), self.qc.s(9*pos+5)
            self.qc.h(9*pos+3), self.qc.h(9*pos+4), self.qc.h(9*pos+5)
            self.qc.cx(9*pos+3, 9*pos+5)
            self.qc.cx(9*pos+4, 9*pos+5)
            self.qc.rz(angle, 9*pos+5)
            self.qc.cx(9*pos+3, 9*pos+5)
            self.qc.cx(9*pos+4, 9*pos+5)
            self.qc.h(9*pos+3), self.qc.h(9*pos+4), self.qc.h(9*pos+5)
            self.qc.sdg(9*pos+3), self.qc.sdg(9*pos+4), self.qc.sdg(9*pos+5)
            self.qc.h(9*pos+3), self.qc.h(9*pos+4), self.qc.h(9*pos+5)
        else:
            self.qc.h(9*pos+1), self.qc.h(9*pos+4), self.qc.h(9*pos+7)
            self.qc.s(9*pos+1), self.qc.s(9*pos+4), self.qc.s(9*pos+7)
            self.qc.h(9*pos+1), self.qc.h(9*pos+4), self.qc.h(9*pos+7)
            self.qc.cx(9*pos+1, 9*pos+7)
            self.qc.cx(9*pos+4, 9*pos+7)
            self.qc.rz(angle, 9*pos+7)
            self.qc.cx(9*pos+1, 9*pos+7)
            self.qc.cx(9*pos+4, 9*pos+7)
            self.qc.h(9*pos+1), self.qc.h(9*pos+4), self.qc.h(9*pos+7)
            self.qc.sdg(9*pos+1), self.qc.sdg(9*pos+4), self.qc.sdg(9*pos+7)
            self.qc.h(9*pos+1), self.qc.h(9*pos+4), self.qc.h(9*pos+7)

    def tdg(self, pos: int):
        anc = self.qc.num_qubits - 1
        self.qc.reset(anc)

        self.qc.h(anc)
        #self.qc.append(h_ideal,[anc])
        self.qc.tdg(anc)

        if self.hadamards[pos]%2 == 0:
            self.qc.cx(3+9*pos, anc)
            self.qc.cx(4+9*pos, anc)
            self.qc.cx(5+9*pos, anc)        
        else:
            self.qc.cx(1+9*pos, anc)
            self.qc.cx(4+9*pos, anc)
            self.qc.cx(7+9*pos, anc)  

        self.qc.measure(anc, 0)
        # if z_stab:
        #     z_qec_ideal(qc, had=had, pos=pos)

        if self.hadamards[pos]%2 == 0:
            with self.qc.if_test((0,1)):
                self.qc.reset(anc)
                self.qc.h(anc)
                self.qc.sdg(anc)
                self.qc.cx(3+9*pos, anc)
                self.qc.cx(4+9*pos, anc)
                self.qc.cx(5+9*pos, anc) 
                self.qc.measure(anc, 0)
                with self.qc.if_test((0,1)):
                    self.qc.z(3+9*pos)
                    self.qc.z(4+9*pos)
                    self.qc.z(5+9*pos)
        else:
            with self.qc.if_test((0,1)):
                self.qc.reset(anc)
                self.qc.h(anc)
                self.qc.sdg(anc)
                self.qc.cx(1+9*pos, anc)
                self.qc.cx(4+9*pos, anc)
                self.qc.cx(7+9*pos, anc) 
                self.qc.measure(anc, 0)
                with self.qc.if_test((0,1)):
                    self.qc.z(1+9*pos)
                    self.qc.z(4+9*pos)
                    self.qc.z(7+9*pos)

    def tdg_cheat(self, pos: int):
        if self.hadamards[pos]%2==0:
            self.qc.cx(9*pos+3, 9*pos+5)
            self.qc.cx(9*pos+4, 9*pos+5)
            self.qc.tdg(9*pos+5)
            self.qc.cx(9*pos+3, 9*pos+5)
            self.qc.cx(9*pos+4, 9*pos+5)
        else:
            self.qc.cx(9*pos+1, 9*pos+7)
            self.qc.cx(9*pos+4, 9*pos+7)
            self.qc.tdg(9*pos+7)
            self.qc.cx(9*pos+1, 9*pos+7)
            self.qc.cx(9*pos+4, 9*pos+7)

    def cs(self):
        self.t(pos=0)
        self.t(pos=1)
        self.cnot(control=0, target=1)
        self.tdg(pos=1)
        self.cnot(control=0, target=1)

    def u2(self, pos: int, gate: list):
        for i in gate:
            if i == "s":
                self.s(pos=pos)
            if i == "sdg":
                self.sdg(pos=pos)
            if i == "t":
                self.t(pos=pos)
            if i == "tdg":
                self.tdg(pos=pos)
            if i == "h":
                self.h(pos=pos)
            if i == "z":
                self.z(pos=pos)
            if i == "x":
                self.x(pos=pos)

    def cu_ramsey(self, Ugates: list):
        self.u2(0, Ugates)
        if self.err:
            self.flagFTec(pos=0)
        self.u2(0, Ugates)

    def qec(self, pos: int):
        self.qec_counter += 1
        anc = self.qc.num_qubits - 1
        if self.hadamards[pos]%2==1:
            #X3 X6 Stabilizer:
            self.qc.reset(anc)
            self.qc.h(anc)
            self.qc.cx(anc, 3+9*pos)
            self.qc.cx(anc, 6+9*pos)
            self.qc.h(anc)
            # self.qc.id(anc)
            self.qc.measure(anc,0)

            #X0 X1 X3 X4 Stabilizer:
            self.qc.reset(anc)
            self.qc.h(anc)
            self.qc.cx(anc, 1+9*pos)
            self.qc.cx(anc, 4+9*pos)
            self.qc.cx(anc, 0+9*pos)
            self.qc.cx(anc, 3+9*pos)
            self.qc.h(anc)
            # self.qc.id(anc)
            self.qc.measure(anc,1)

            #X4 X5 X7 X8 Stabilizer:
            self.qc.reset(anc)
            self.qc.h(anc)
            self.qc.cx(anc, 4+9*pos)
            self.qc.cx(anc, 7+9*pos)
            self.qc.cx(anc, 5+9*pos)
            self.qc.cx(anc, 8+9*pos)
            self.qc.h(anc)
            # self.qc.id(anc)
            self.qc.measure(anc,2)

            #X2 X5 Stabilizer:
            self.qc.reset(anc)
            self.qc.h(anc)
            self.qc.cx(anc, 2+9*pos)
            self.qc.cx(anc, 5+9*pos)
            self.qc.h(anc)
            # self.qc.id(anc)
            self.qc.measure(anc,3)

            with self.qc.if_test((0,1)):             #6
                with self.qc.if_test((1,0)):
                    self.qc.z(6+9*pos)

            with self.qc.if_test((0,1)):             #3
                with self.qc.if_test((1,1)):
                    self.qc.z(3+9*pos)

            with self.qc.if_test((3,1)):             #2
                with self.qc.if_test((2,0)):
                    self.qc.z(2+9*pos)
            
            with self.qc.if_test((3,1)):             #5
                with self.qc.if_test((2,1)):
                    self.qc.z(5+9*pos)
            
            with self.qc.if_test((1,1)):             #4
                with self.qc.if_test((2,1)):
                    self.qc.z(4+9*pos)

            with self.qc.if_test((0,0)):             #0 und 1
                with self.qc.if_test((1,1)):
                    with self.qc.if_test((2,0)):
                        self.qc.z(0+9*pos)
            
            with self.qc.if_test((1,0)):             #7 und 8
                with self.qc.if_test((2,1)):
                    with self.qc.if_test((3,0)):
                        self.qc.z(7+9*pos)

        ###########################################################################################################

            #Z0 Z1 Stabilizer:
            self.qc.reset(anc)
            self.qc.cx(0+9*pos, anc)
            self.qc.cx(1+9*pos, anc)
            # self.qc.id(anc)
            self.qc.measure(anc,0)

            #Z1 Z2 Z4 Z5 Stabilizer:
            self.qc.reset(anc)
            self.qc.cx(4+9*pos, anc)
            self.qc.cx(5+9*pos, anc)
            self.qc.cx(1+9*pos, anc)
            self.qc.cx(2+9*pos, anc)
            # self.qc.id(anc)
            self.qc.measure(anc,1)
        
            #Z3 Z4 Z6 Z7 Stabilizer:
            self.qc.reset(anc)
            self.qc.cx(3+9*pos, anc)
            self.qc.cx(4+9*pos, anc)
            self.qc.cx(6+9*pos, anc)
            self.qc.cx(7+9*pos, anc)
            # self.qc.id(anc)
            self.qc.measure(anc,2)

            #Z7 Z8 Stabilizer:
            self.qc.reset(anc)
            self.qc.cx(7+9*pos, anc)
            self.qc.cx(8+9*pos, anc)
            # self.qc.id(anc)
            self.qc.measure(anc,3)
            
            with self.qc.if_test((0,1)):             #0
                with self.qc.if_test((1,0)):
                    self.qc.x(0+9*pos)

            with self.qc.if_test((0,1)):             #1
                with self.qc.if_test((1,1)):
                    self.qc.x(1+9*pos)
            
            with self.qc.if_test((3,1)):             #8
                with self.qc.if_test((2,0)):
                    self.qc.x(8+9*pos)
            
            with self.qc.if_test((3,1)):             #7
                with self.qc.if_test((2,1)):
                    self.qc.x(7+9*pos)
            
            with self.qc.if_test((1,1)):             #4
                with self.qc.if_test((2,1)):
                    self.qc.x(4+9*pos)

            with self.qc.if_test((0,0)):             #2 und 5
                with self.qc.if_test((1,1)):
                    with self.qc.if_test((2,0)):
                        self.qc.x(2+9*pos)
            
            with self.qc.if_test((1,0)):             #3 und 6
                with self.qc.if_test((2,1)):
                    with self.qc.if_test((3,0)):
                        self.qc.x(3+9*pos)

        else:
            #X0 X1 Stabilizer:
            self.qc.reset(anc)
            self.qc.h(anc)
            self.qc.cx(anc, 0+9*pos)
            self.qc.cx(anc, 1+9*pos)
            self.qc.h(anc)
            # self.qc.id(anc)
            self.qc.measure(anc,0)
            
            #X1 X2 X4 X5 Stabilizer:
            self.qc.reset(anc)
            self.qc.h(anc)
            self.qc.cx(anc, 4+9*pos)
            self.qc.cx(anc, 5+9*pos)
            self.qc.cx(anc, 1+9*pos)
            self.qc.cx(anc, 2+9*pos)
            self.qc.h(anc)
            # self.qc.id(anc)
            self.qc.measure(anc,1)

            #X3 X4 X6 X7 Stabilizer:
            self.qc.reset(anc)
            self.qc.h(anc)
            self.qc.cx(anc, 3+9*pos)
            self.qc.cx(anc, 4+9*pos)
            self.qc.cx(anc, 6+9*pos)
            self.qc.cx(anc, 7+9*pos)
            self.qc.h(anc)
            # self.qc.id(anc)
            self.qc.measure(anc,2)

            #X7 X8 Stabilizer:
            self.qc.reset(anc)
            self.qc.h(anc)
            self.qc.cx(anc, 7+9*pos)
            self.qc.cx(anc, 8+9*pos)
            self.qc.h(anc)
            # self.qc.id(anc)
            self.qc.measure(anc,3)

            with self.qc.if_test((0,1)):             #0
                with self.qc.if_test((1,0)):    
                    self.qc.z(0+9*pos)

            with self.qc.if_test((0,1)):             #1
                with self.qc.if_test((1,1)):
                    self.qc.z(1+9*pos)

            with self.qc.if_test((1,1)):             #4
                with self.qc.if_test((2,1)):
                    self.qc.z(4+9*pos)

            with self.qc.if_test((2,1)):             #7
                with self.qc.if_test((3,1)):
                    self.qc.z(7+9*pos)

            with self.qc.if_test((2,0)):             #8
                with self.qc.if_test((3,1)):
                    self.qc.z(8+9*pos)

            with self.qc.if_test((0,0)):             #2 und 5
                with self.qc.if_test((1,1)):        
                    with self.qc.if_test((2,0)):    
                        self.qc.z(2+9*pos)

            with self.qc.if_test((1,0)):             #3 und 6
                with self.qc.if_test((2,1)):
                    with self.qc.if_test((3,0)):
                        self.qc.z(3+9*pos)

        ###########################################################################################################

            #Z3 Z6 Stabilizer:
            self.qc.reset(anc)
            self.qc.cx(3+9*pos, anc)
            self.qc.cx(6+9*pos, anc)
            # self.qc.id(anc)
            self.qc.measure(anc,0)

            #Z0 Z1 Z3 Z4 Stabilizer:
            self.qc.reset(anc)
            self.qc.cx(4+9*pos, anc)
            self.qc.cx(1+9*pos, anc)
            self.qc.cx(0+9*pos, anc)
            self.qc.cx(3+9*pos, anc)
            # self.qc.id(anc)
            self.qc.measure(anc,1)
        
            #Z4 Z5 Z7 Z8 Stabilizer:
            self.qc.reset(anc)
            self.qc.cx(4+9*pos, anc)
            self.qc.cx(7+9*pos, anc)
            self.qc.cx(5+9*pos, anc)
            self.qc.cx(8+9*pos, anc)
            # self.qc.id(anc)
            self.qc.measure(anc,2)

            #Z2 Z5 Stabilizer:
            self.qc.reset(anc)
            self.qc.cx(2+9*pos, anc)
            self.qc.cx(5+9*pos, anc)
            # self.qc.id(anc)
            self.qc.measure(anc,3)
            
            with self.qc.if_test((0,1)):             #6
                with self.qc.if_test((1,0)):
                    self.qc.x(6+9*pos)

            with self.qc.if_test((0,1)):             #3
                with self.qc.if_test((1,1)):
                    self.qc.x(3+9*pos)

            with self.qc.if_test((1,1)):             #4
                with self.qc.if_test((2,1)):
                    self.qc.x(4+9*pos)

            with self.qc.if_test((2,0)):             #2
                with self.qc.if_test((3,1)):
                    self.qc.x(2+9*pos)

            with self.qc.if_test((2,1)):             #5
                with self.qc.if_test((3,1)):
                    self.qc.x(5+9*pos)
            
            with self.qc.if_test((0,0)):             #0 und 1
                with self.qc.if_test((1,1)):
                    with self.qc.if_test((2,0)):
                        self.qc.x(0+9*pos)
            
            with self.qc.if_test((1,0)):             #7 und 8
                with self.qc.if_test((2,1)):
                    with self.qc.if_test((3,0)):
                        self.qc.x(7+9*pos)

    def qec_ideal(self, pos: int):
        self.qec_counter += 1
        anc = self.qc.num_qubits - 1
        if self.hadamards[pos]%2==1:
            #X3 X6 Stabilizer:
            self.qc.reset(anc)
            self.qc.append(h_ideal,[anc])
            self.qc.append(cx_ideal, [3+9*pos, anc])
            self.qc.append(cx_ideal, [6+9*pos, anc])
            self.qc.append(h_ideal,[anc])
            # self.qc.id(anc)
            self.qc.measure(anc,0)

            #X0 X1 X3 X4 Stabilizer:
            self.qc.reset(anc)
            self.qc.append(h_ideal,[anc])
            self.qc.append(cx_ideal, [1+9*pos, anc])
            self.qc.append(cx_ideal, [4+9*pos, anc])
            self.qc.append(cx_ideal, [0+9*pos, anc])
            self.qc.append(cx_ideal, [3+9*pos, anc])
            self.qc.append(h_ideal,[anc])
            # self.qc.id(anc)
            self.qc.measure(anc,1)

            #X4 X5 X7 X8 Stabilizer:
            self.qc.reset(anc)
            self.qc.append(h_ideal,[anc])
            self.qc.append(cx_ideal, [4+9*pos, anc])
            self.qc.append(cx_ideal, [7+9*pos, anc])
            self.qc.append(cx_ideal, [5+9*pos, anc])
            self.qc.append(cx_ideal, [8+9*pos, anc])
            self.qc.append(h_ideal,[anc])
            # self.qc.id(anc)
            self.qc.measure(anc,2)

            #X2 X5 Stabilizer:
            self.qc.reset(anc)
            self.qc.append(h_ideal,[anc])
            self.qc.append(cx_ideal, [2+9*pos, anc])
            self.qc.append(cx_ideal, [5+9*pos, anc])
            self.qc.append(h_ideal,[anc])
            # self.qc.id(anc)
            self.qc.measure(anc,3)

            with self.qc.if_test((0,1)):             #6
                with self.qc.if_test((1,0)):
                    self.qc.append(z_ideal,[6+9*pos])

            with self.qc.if_test((0,1)):             #3
                with self.qc.if_test((1,1)):
                    self.qc.append(z_ideal,[3+9*pos])

            with self.qc.if_test((3,1)):             #2
                with self.qc.if_test((2,0)):
                    self.qc.append(z_ideal,[2+9*pos])
            
            with self.qc.if_test((3,1)):             #5
                with self.qc.if_test((2,1)):
                    self.qc.append(z_ideal,[5+9*pos])
            
            with self.qc.if_test((1,1)):             #4
                with self.qc.if_test((2,1)):
                    self.qc.append(z_ideal,[4+9*pos])

            with self.qc.if_test((0,0)):             #0 und 1
                with self.qc.if_test((1,1)):
                    with self.qc.if_test((2,0)):
                        self.qc.append(z_ideal,[0+9*pos])
            
            with self.qc.if_test((1,0)):             #7 und 8
                with self.qc.if_test((2,1)):
                    with self.qc.if_test((3,0)):
                        self.qc.append(z_ideal,[7+9*pos])

        ###########################################################################################################

            #Z0 Z1 Stabilizer:
            self.qc.reset(anc)
            self.qc.append(cx_ideal, [anc, 0+9*pos])
            self.qc.append(cx_ideal, [anc, 1+9*pos])
            # self.qc.id(anc)
            self.qc.measure(anc,0)

            #Z1 Z2 Z4 Z5 Stabilizer:
            self.qc.reset(anc)
            self.qc.append(cx_ideal, [anc, 4+9*pos])
            self.qc.append(cx_ideal, [anc, 5+9*pos])
            self.qc.append(cx_ideal, [anc, 1+9*pos])
            self.qc.append(cx_ideal, [anc, 2+9*pos])
            # self.qc.id(anc)
            self.qc.measure(anc,1)
        
            #Z3 Z4 Z6 Z7 Stabilizer:
            self.qc.reset(anc)
            self.qc.append(cx_ideal, [anc, 3+9*pos])
            self.qc.append(cx_ideal, [anc, 4+9*pos])
            self.qc.append(cx_ideal, [anc, 6+9*pos])
            self.qc.append(cx_ideal, [anc, 7+9*pos])
            # self.qc.id(anc)
            self.qc.measure(anc,2)

            #Z7 Z8 Stabilizer:
            self.qc.reset(anc)
            self.qc.append(cx_ideal, [anc, 7+9*pos])
            self.qc.append(cx_ideal, [anc, 8+9*pos])
            # self.qc.id(anc)
            self.qc.measure(anc,3)
            
            with self.qc.if_test((0,1)):             #0
                with self.qc.if_test((1,0)):
                    self.qc.append(x_ideal,[0+9*pos])

            with self.qc.if_test((0,1)):             #1
                with self.qc.if_test((1,1)):
                    self.qc.append(x_ideal,[1+9*pos])
            
            with self.qc.if_test((3,1)):             #8
                with self.qc.if_test((2,0)):
                    self.qc.append(x_ideal,[8+9*pos])
            
            with self.qc.if_test((3,1)):             #7
                with self.qc.if_test((2,1)):
                    self.qc.append(x_ideal,[7+9*pos])
            
            with self.qc.if_test((1,1)):             #4
                with self.qc.if_test((2,1)):
                    self.qc.append(x_ideal,[4+9*pos])

            with self.qc.if_test((0,0)):             #2 und 5
                with self.qc.if_test((1,1)):
                    with self.qc.if_test((2,0)):
                        self.qc.append(x_ideal,[2+9*pos])
            
            with self.qc.if_test((1,0)):             #3 und 6
                with self.qc.if_test((2,1)):
                    with self.qc.if_test((3,0)):
                        self.qc.append(x_ideal,[3+9*pos])

        else:
            #X0 X1 Stabilizer:
            self.qc.reset(anc)
            self.qc.append(h_ideal,[anc])
            self.qc.append(cx_ideal, [0+9*pos, anc])
            self.qc.append(cx_ideal, [1+9*pos, anc])
            self.qc.append(h_ideal,[anc])
            # self.qc.id(anc)
            self.qc.measure(anc,0)
            
            #X1 X2 X4 X5 Stabilizer:
            self.qc.reset(anc)
            self.qc.append(h_ideal,[anc])
            self.qc.append(cx_ideal, [4+9*pos, anc])
            self.qc.append(cx_ideal, [5+9*pos, anc])
            self.qc.append(cx_ideal, [1+9*pos, anc])
            self.qc.append(cx_ideal, [2+9*pos, anc])
            self.qc.append(h_ideal,[anc])
            # self.qc.id(anc)
            self.qc.measure(anc,1)

            #X3 X4 X6 X7 Stabilizer:
            self.qc.reset(anc)
            self.qc.append(h_ideal,[anc])
            self.qc.append(cx_ideal, [3+9*pos, anc])
            self.qc.append(cx_ideal, [4+9*pos, anc])
            self.qc.append(cx_ideal, [6+9*pos, anc])
            self.qc.append(cx_ideal, [7+9*pos, anc])
            self.qc.append(h_ideal,[anc])
            # self.qc.id(anc)
            self.qc.measure(anc,2)

            #X7 X8 Stabilizer:
            self.qc.reset(anc)
            self.qc.append(h_ideal,[anc])
            self.qc.append(cx_ideal, [7+9*pos, anc])
            self.qc.append(cx_ideal, [8+9*pos, anc])
            self.qc.append(h_ideal,[anc])
            # self.qc.id(anc)
            self.qc.measure(anc,3)

            with self.qc.if_test((0,1)):             #0
                with self.qc.if_test((1,0)):    
                    self.qc.append(z_ideal,[0+9*pos])

            with self.qc.if_test((0,1)):             #1
                with self.qc.if_test((1,1)):
                    self.qc.append(z_ideal,[1+9*pos])

            with self.qc.if_test((1,1)):             #4
                with self.qc.if_test((2,1)):
                    self.qc.append(z_ideal,[4+9*pos])

            with self.qc.if_test((2,1)):             #7
                with self.qc.if_test((3,1)):
                    self.qc.append(z_ideal,[7+9*pos])

            with self.qc.if_test((2,0)):             #8
                with self.qc.if_test((3,1)):
                    self.qc.append(z_ideal,[8+9*pos])

            with self.qc.if_test((0,0)):             #2 und 5
                with self.qc.if_test((1,1)):        
                    with self.qc.if_test((2,0)):    
                        self.qc.append(z_ideal,[2+9*pos])

            with self.qc.if_test((1,0)):             #3 und 6
                with self.qc.if_test((2,1)):
                    with self.qc.if_test((3,0)):
                        self.qc.append(z_ideal,[3+9*pos])

        ###########################################################################################################

            #Z3 Z6 Stabilizer:
            self.qc.reset(anc)
            self.qc.append(cx_ideal, [anc, 3+9*pos])
            self.qc.append(cx_ideal, [anc, 6+9*pos])
            # self.qc.id(anc)
            self.qc.measure(anc,0)

            #Z0 Z1 Z3 Z4 Stabilizer:
            self.qc.reset(anc)
            self.qc.append(cx_ideal, [anc, 4+9*pos])
            self.qc.append(cx_ideal, [anc, 1+9*pos])
            self.qc.append(cx_ideal, [anc, 0+9*pos])
            self.qc.append(cx_ideal, [anc, 3+9*pos])
            # self.qc.id(anc)
            self.qc.measure(anc,1)
        
            #Z4 Z5 Z7 Z8 Stabilizer:
            self.qc.reset(anc)
            self.qc.append(cx_ideal, [anc, 4+9*pos])
            self.qc.append(cx_ideal, [anc, 7+9*pos])
            self.qc.append(cx_ideal, [anc, 5+9*pos])
            self.qc.append(cx_ideal, [anc, 8+9*pos])
            # self.qc.id(anc)
            self.qc.measure(anc,2)

            #Z2 Z5 Stabilizer:
            self.qc.reset(anc)
            self.qc.append(cx_ideal, [anc, 2+9*pos])
            self.qc.append(cx_ideal, [anc, 5+9*pos])
            # self.qc.id(anc)
            self.qc.measure(anc,3)
            
            with self.qc.if_test((0,1)):             #6
                with self.qc.if_test((1,0)):
                    self.qc.append(x_ideal,[6+9*pos])

            with self.qc.if_test((0,1)):             #3
                with self.qc.if_test((1,1)):
                    self.qc.append(x_ideal,[3+9*pos])

            with self.qc.if_test((1,1)):             #4
                with self.qc.if_test((2,1)):
                    self.qc.append(x_ideal,[4+9*pos])

            with self.qc.if_test((2,0)):             #2
                with self.qc.if_test((3,1)):
                    self.qc.append(x_ideal,[2+9*pos])

            with self.qc.if_test((2,1)):             #5
                with self.qc.if_test((3,1)):
                    self.qc.append(x_ideal,[5+9*pos])
            
            with self.qc.if_test((0,0)):             #0 und 1
                with self.qc.if_test((1,1)):
                    with self.qc.if_test((2,0)):
                        self.qc.append(x_ideal,[0+9*pos])
            
            with self.qc.if_test((1,0)):             #7 und 8
                with self.qc.if_test((2,1)):
                    with self.qc.if_test((3,0)):
                        self.qc.append(x_ideal,[7+9*pos])

    def qec_zstab(self, pos: int):
        self.qec_counter += 1
        anc = self.qc.num_qubits - 1
        if self.hadamards[pos]%2==1:
            #Z0 Z1 Stabilizer:
            self.qc.reset(anc)
            self.qc.cx(0+9*pos, anc)
            self.qc.cx(1+9*pos, anc)
            self.qc.id(anc)
            self.qc.measure(anc,0)

            #Z1 Z2 Z4 Z5 Stabilizer:
            self.qc.reset(anc)
            self.qc.cx(4+9*pos, anc)
            self.qc.cx(5+9*pos, anc)
            self.qc.cx(1+9*pos, anc)
            self.qc.cx(2+9*pos, anc)
            self.qc.id(anc)
            self.qc.measure(anc,1)
        
            #Z3 Z4 Z6 Z7 Stabilizer:
            self.qc.reset(anc)
            self.qc.cx(3+9*pos, anc)
            self.qc.cx(4+9*pos, anc)
            self.qc.cx(6+9*pos, anc)
            self.qc.cx(7+9*pos, anc)
            self.qc.id(anc)
            self.qc.measure(anc,2)

            #Z7 Z8 Stabilizer:
            self.qc.reset(anc)
            self.qc.cx(7+9*pos, anc)
            self.qc.cx(8+9*pos, anc)
            self.qc.id(anc)
            self.qc.measure(anc,3)
            
            with self.qc.if_test((0,1)):             #0
                with self.qc.if_test((1,0)):
                    self.qc.x(0+9*pos)

            with self.qc.if_test((0,1)):             #1
                with self.qc.if_test((1,1)):
                    self.qc.x(1+9*pos)
            
            with self.qc.if_test((3,1)):             #8
                with self.qc.if_test((2,0)):
                    self.qc.x(8+9*pos)
            
            with self.qc.if_test((3,1)):             #7
                with self.qc.if_test((2,1)):
                    self.qc.x(7+9*pos)
            
            with self.qc.if_test((1,1)):             #4
                with self.qc.if_test((2,1)):
                    self.qc.x(4+9*pos)

            with self.qc.if_test((0,0)):             #2 und 5
                with self.qc.if_test((1,1)):
                    with self.qc.if_test((2,0)):
                        self.qc.x(2+9*pos)
            
            with self.qc.if_test((1,0)):             #3 und 6
                with self.qc.if_test((2,1)):
                    with self.qc.if_test((3,0)):
                        self.qc.x(3+9*pos)

        else:
            #Z3 Z6 Stabilizer:
            self.qc.reset(anc)
            self.qc.cx(3+9*pos, anc)
            self.qc.cx(6+9*pos, anc)
            self.qc.id(anc)
            self.qc.measure(anc,0)

            #Z0 Z1 Z3 Z4 Stabilizer:
            self.qc.reset(anc)
            self.qc.cx(4+9*pos, anc)
            self.qc.cx(1+9*pos, anc)
            self.qc.cx(0+9*pos, anc)
            self.qc.cx(3+9*pos, anc)
            self.qc.id(anc)
            self.qc.measure(anc,1)
        
            #Z4 Z5 Z7 Z8 Stabilizer:
            self.qc.reset(anc)
            self.qc.cx(4+9*pos, anc)
            self.qc.cx(7+9*pos, anc)
            self.qc.cx(5+9*pos, anc)
            self.qc.cx(8+9*pos, anc)
            self.qc.id(anc)
            self.qc.measure(anc,2)

            #Z2 Z5 Stabilizer:
            self.qc.reset(anc)
            self.qc.cx(2+9*pos, anc)
            self.qc.cx(5+9*pos, anc)
            self.qc.id(anc)
            self.qc.measure(anc,3)
            
            with self.qc.if_test((0,1)):             #6
                with self.qc.if_test((1,0)):
                    self.qc.x(6+9*pos)

            with self.qc.if_test((0,1)):             #3
                with self.qc.if_test((1,1)):
                    self.qc.x(3+9*pos)

            with self.qc.if_test((1,1)):             #4
                with self.qc.if_test((2,1)):
                    self.qc.x(4+9*pos)

            with self.qc.if_test((2,0)):             #2
                with self.qc.if_test((3,1)):
                    self.qc.x(2+9*pos)

            with self.qc.if_test((2,1)):             #5
                with self.qc.if_test((3,1)):
                    self.qc.x(5+9*pos)
            
            with self.qc.if_test((0,0)):             #0 und 1
                with self.qc.if_test((1,1)):
                    with self.qc.if_test((2,0)):
                        self.qc.x(0+9*pos)
            
            with self.qc.if_test((1,0)):             #7 und 8
                with self.qc.if_test((2,1)):
                    with self.qc.if_test((3,0)):
                        self.qc.x(7+9*pos)

    def qec_flags(self, pos: int):
        flags = ClassicalRegister(8)
        self.qc.add_register(flags)
        anc = self.qc.num_qubits - 1
        ancc = anc - 1
        if self.hadamards[pos]%2==1:
            #X3 X6 Stabilizer:
            self.qc.reset(anc)
            self.qc.id(anc)
            self.qc.h(anc)
            self.qc.cx(anc, 3+9*pos)
            self.qc.cx(anc, 6+9*pos)
            self.qc.h(anc)
            self.qc.id(anc)
            self.qc.measure(anc,0)

            #X0 X1 X3 X4 Stabilizer:
            self.qc.reset(anc), self.qc.reset(ancc)
            self.qc.id(anc), self.qc.id(ancc)
            self.qc.h(anc)
            self.qc.cx(anc, ancc)
            self.qc.cx(anc, 0+9*pos)
            self.qc.cx(anc, 1+9*pos)
            self.qc.cx(anc, 3+9*pos)
            self.qc.cx(anc, ancc)
            self.qc.cx(anc, 4+9*pos)
            self.qc.h(anc)
            self.qc.id(anc), self.qc.id(ancc)
            self.qc.measure(anc,1), self.qc.measure(ancc, flags[1])

            #X4 X5 X7 X8 Stabilizer:
            self.qc.reset(anc), self.qc.reset(ancc)
            self.qc.id(anc), self.qc.id(ancc)
            self.qc.h(anc)
            self.qc.cx(anc, ancc)
            self.qc.cx(anc, 4+9*pos)
            self.qc.cx(anc, 5+9*pos)
            self.qc.cx(anc, 7+9*pos)
            self.qc.cx(anc, ancc)
            self.qc.cx(anc, 8+9*pos)
            self.qc.h(anc)
            self.qc.id(anc), self.qc.id(ancc)
            self.qc.measure(anc,2), self.qc.measure(ancc, flags[2])

            #X2 X5 Stabilizer:
            self.qc.reset(anc)
            self.qc.id(anc)
            self.qc.h(anc)
            self.qc.cx(anc, 2+9*pos)
            self.qc.cx(anc, 5+9*pos)
            self.qc.h(anc)
            self.qc.id(anc)
            self.qc.measure(anc,3)

            with self.qc.if_test((0,1)):             #6
                with self.qc.if_test((1,0)):
                    self.qc.z(6+9*pos)

            with self.qc.if_test((0,1)):             #3
                with self.qc.if_test((1,1)):
                    self.qc.z(3+9*pos)

            with self.qc.if_test((3,1)):             #2
                with self.qc.if_test((2,0)):
                    self.qc.z(2+9*pos)
            
            with self.qc.if_test((3,1)):             #5
                with self.qc.if_test((2,1)):
                    self.qc.z(5+9*pos)
            
            with self.qc.if_test((1,1)):             #4
                with self.qc.if_test((2,1)):
                    self.qc.z(4+9*pos)

            with self.qc.if_test((0,0)):             #0 und 1
                with self.qc.if_test((1,1)):
                    with self.qc.if_test((2,0)):
                        self.qc.z(0+9*pos)
            
            with self.qc.if_test((1,0)):             #7 und 8
                with self.qc.if_test((2,1)):
                    with self.qc.if_test((3,0)):
                        self.qc.z(7+9*pos)

        ###########################################################################################################

            #Z0 Z1 Stabilizer:
            self.qc.reset(anc)
            self.qc.id(anc)
            self.qc.cx(0+9*pos, anc)
            self.qc.cx(1+9*pos, anc)
            self.qc.id(anc)
            self.qc.measure(anc,0)

            #Z1 Z2 Z4 Z5 Stabilizer:
            self.qc.reset(anc), self.qc.reset(ancc)
            self.qc.id(anc), self.qc.id(ancc)
            self.qc.h(ancc)
            self.qc.cx(ancc, anc)
            self.qc.cx(1+9*pos, anc)
            self.qc.cx(2+9*pos, anc)
            self.qc.cx(4+9*pos, anc)
            self.qc.cx(ancc, anc)
            self.qc.cx(5+9*pos, anc)
            self.qc.h(ancc)
            self.qc.id(anc), self.qc.id(ancc)
            self.qc.measure(anc,1), self.qc.measure(ancc, flags[5])
        
            #Z3 Z4 Z6 Z7 Stabilizer:
            self.qc.reset(anc), self.qc.reset(ancc)
            self.qc.id(anc), self.qc.id(ancc)
            self.qc.h(ancc)
            self.qc.cx(ancc, anc)
            self.qc.cx(3+9*pos, anc)
            self.qc.cx(4+9*pos, anc)
            self.qc.cx(6+9*pos, anc)
            self.qc.cx(ancc, anc)
            self.qc.cx(7+9*pos, anc)
            self.qc.h(ancc)
            self.qc.id(anc), self.qc.id(ancc)
            self.qc.measure(anc,2), self.qc.measure(ancc, flags[6])

            #Z7 Z8 Stabilizer:
            self.qc.reset(anc)
            self.qc.id(anc)
            self.qc.cx(7+9*pos, anc)
            self.qc.cx(8+9*pos, anc)
            self.qc.id(anc)
            self.qc.measure(anc,3)
            
            with self.qc.if_test((0,1)):             #0
                with self.qc.if_test((1,0)):
                    self.qc.x(0+9*pos)

            with self.qc.if_test((0,1)):             #1
                with self.qc.if_test((1,1)):
                    self.qc.x(1+9*pos)
            
            with self.qc.if_test((3,1)):             #8
                with self.qc.if_test((2,0)):
                    self.qc.x(8+9*pos)
            
            with self.qc.if_test((3,1)):             #7
                with self.qc.if_test((2,1)):
                    self.qc.x(7+9*pos)
            
            with self.qc.if_test((1,1)):             #4
                with self.qc.if_test((2,1)):
                    self.qc.x(4+9*pos)

            with self.qc.if_test((0,0)):             #2 und 5
                with self.qc.if_test((1,1)):
                    with self.qc.if_test((2,0)):
                        self.qc.x(2+9*pos)
            
            with self.qc.if_test((1,0)):             #3 und 6
                with self.qc.if_test((2,1)):
                    with self.qc.if_test((3,0)):
                        self.qc.x(3+9*pos)

        else:
            #X0 X1 Stabilizer:
            self.qc.reset(anc)
            self.qc.id(anc)
            self.qc.h(anc)
            self.qc.cx(anc, 0+9*pos)
            self.qc.cx(anc, 1+9*pos)
            self.qc.h(anc)
            self.qc.id(anc)
            self.qc.measure(anc,0)
            
            #X1 X2 X4 X5 Stabilizer:
            self.qc.reset(anc), self.qc.reset(ancc)
            self.qc.id(anc), self.qc.id(ancc)
            self.qc.h(anc)
            self.qc.cx(anc, ancc)
            self.qc.cx(anc, 1+9*pos)
            self.qc.cx(anc, 2+9*pos)
            self.qc.cx(anc, 4+9*pos)
            self.qc.cx(anc, ancc)
            self.qc.cx(anc, 5+9*pos)
            self.qc.h(anc)
            self.qc.id(anc), self.qc.id(ancc)
            self.qc.measure(anc,1), self.qc.measure(ancc, flags[1])

            #X3 X4 X6 X7 Stabilizer:
            self.qc.reset(anc), self.qc.reset(ancc)
            self.qc.id(anc), self.qc.id(ancc)
            self.qc.h(anc)
            self.qc.cx(anc, ancc)
            self.qc.cx(anc, 3+9*pos)
            self.qc.cx(anc, 4+9*pos)
            self.qc.cx(anc, 6+9*pos)
            self.qc.cx(anc, ancc)
            self.qc.cx(anc, 7+9*pos)
            self.qc.h(anc)
            self.qc.id(anc), self.qc.id(ancc)
            self.qc.measure(anc,2), self.qc.measure(ancc, flags[2])

            #X7 X8 Stabilizer:
            self.qc.reset(anc)
            self.qc.id(anc)
            self.qc.h(anc)
            self.qc.cx(anc, 7+9*pos)
            self.qc.cx(anc, 8+9*pos)
            self.qc.h(anc)
            self.qc.id(anc)
            self.qc.measure(anc,3)

            with self.qc.if_test((0,1)):             #0
                with self.qc.if_test((1,0)):    
                    self.qc.z(0+9*pos)

            with self.qc.if_test((0,1)):             #1
                with self.qc.if_test((1,1)):
                    self.qc.z(1+9*pos)

            with self.qc.if_test((1,1)):             #4
                with self.qc.if_test((2,1)):
                    self.qc.z(4+9*pos)

            with self.qc.if_test((2,1)):             #7
                with self.qc.if_test((3,1)):
                    self.qc.z(7+9*pos)

            with self.qc.if_test((2,0)):             #8
                with self.qc.if_test((3,1)):
                    self.qc.z(8+9*pos)

            with self.qc.if_test((0,0)):             #2 und 5
                with self.qc.if_test((1,1)):        
                    with self.qc.if_test((2,0)):    
                        self.qc.z(2+9*pos)

            with self.qc.if_test((1,0)):             #3 und 6
                with self.qc.if_test((2,1)):
                    with self.qc.if_test((3,0)):
                        self.qc.z(3+9*pos)

        ###########################################################################################################

            #Z3 Z6 Stabilizer:
            self.qc.reset(anc)
            self.qc.id(anc)
            self.qc.cx(3+9*pos, anc)
            self.qc.cx(6+9*pos, anc)
            self.qc.id(anc)
            self.qc.measure(anc,0)

            #Z0 Z1 Z3 Z4 Stabilizer:
            self.qc.reset(anc), self.qc.reset(ancc)
            self.qc.id(anc), self.qc.id(ancc)
            self.qc.h(ancc)
            self.qc.cx(ancc, anc)
            self.qc.cx(0+9*pos, anc)
            self.qc.cx(1+9*pos, anc)
            self.qc.cx(3+9*pos, anc)
            self.qc.cx(ancc, anc)
            self.qc.cx(4+9*pos, anc)
            self.qc.h(ancc)
            self.qc.id(anc), self.qc.id(ancc)
            self.qc.measure(anc,1), self.qc.measure(ancc, flags[5])
        
            #Z4 Z5 Z7 Z8 Stabilizer:
            self.qc.reset(anc), self.qc.reset(ancc)
            self.qc.id(anc), self.qc.id(ancc)
            self.qc.h(ancc)
            self.qc.cx(ancc, anc)
            self.qc.cx(4+9*pos, anc)
            self.qc.cx(5+9*pos, anc)
            self.qc.cx(7+9*pos, anc)
            self.qc.cx(ancc, anc)
            self.qc.cx(8+9*pos, anc)
            self.qc.h(ancc)
            self.qc.id(anc), self.qc.id(ancc)
            self.qc.measure(anc,2), self.qc.measure(ancc, flags[6])

            #Z2 Z5 Stabilizer:
            self.qc.reset(anc)
            self.qc.id(anc)
            self.qc.cx(2+9*pos, anc)
            self.qc.cx(5+9*pos, anc)
            self.qc.id(anc)
            self.qc.measure(anc,3)
            
            with self.qc.if_test((0,1)):             #6
                with self.qc.if_test((1,0)):
                    self.qc.x(6+9*pos)

            with self.qc.if_test((0,1)):             #3
                with self.qc.if_test((1,1)):
                    self.qc.x(3+9*pos)

            with self.qc.if_test((1,1)):             #4
                with self.qc.if_test((2,1)):
                    self.qc.x(4+9*pos)

            with self.qc.if_test((2,0)):             #2
                with self.qc.if_test((3,1)):
                    self.qc.x(2+9*pos)

            with self.qc.if_test((2,1)):             #5
                with self.qc.if_test((3,1)):
                    self.qc.x(5+9*pos)
            
            with self.qc.if_test((0,0)):             #0 und 1
                with self.qc.if_test((1,1)):
                    with self.qc.if_test((2,0)):
                        self.qc.x(0+9*pos)
            
            with self.qc.if_test((1,0)):             #7 und 8
                with self.qc.if_test((2,1)):
                    with self.qc.if_test((3,0)):
                        self.qc.x(7+9*pos)

################################# Neues FTQEC Procol based on arxiv:1708.02246, page 5, "Flag 1-FTEC Protocol" ##########################
    def _measure_stab(self, pos: int, stab: tuple, xtype: bool, bit, flags = None, flagbit = None):
        #one stabilizer measurement. xtype: the ancilla is the control (X-Stabilizer), otherwise the data
        #qubits are the controls and the ancilla is the target (Z-Stabilizer). A flag ancilla is only
        #inserted for the weight 4 stabilizers, the weight 2 ones are fault tolerant on their own: a single
        #ancilla fault there propagates to at most one data qubit or to the whole stabilizer.
        anc = self.qc.num_qubits - 1
        ancc = anc - 1
        flag = flags is not None and len(stab) == 4

        self.qc.reset(anc)
        self.qc.id(anc)
        if flag:
            self.qc.reset(ancc)
            self.qc.id(ancc)

        #the flag ancilla brackets the first three data CNOTs, same ordering as in qec_flags()
        order = ["f", stab[0], stab[1], stab[2], "f", stab[3]] if flag else list(stab)

        if xtype:
            self.qc.h(anc)
            for q in order:
                self.qc.cx(anc, ancc if q == "f" else q+9*pos)
            self.qc.h(anc)
        else:
            if flag:
                self.qc.h(ancc)
            for q in order:
                self.qc.cx(ancc if q == "f" else q+9*pos, anc)
            if flag:
                self.qc.h(ancc)

        self.qc.id(anc)
        self.qc.measure(anc, bit)
        if flag:
            self.qc.id(ancc)
            self.qc.measure(ancc, flags[flagbit])

    def _round(self, pos: int, flags = None):
        #one full syndrome extraction round. self.hadamards decides which set is measured as X- and which
        #as Z-Stabilizers, the supports themselves do not move under H_L.
        self.qec_counter += 1
        xtype = self.hadamards[pos]%2 == 0
        for i, stab in enumerate(STABS_A_9Q):
            self._measure_stab(pos, stab, xtype, self.qecc[3-i], flags, i)
        for i, stab in enumerate(STABS_B_9Q):
            self._measure_stab(pos, stab, not xtype, self.qecc[7-i], flags, 4+i)

    def flagsyndrome(self, pos: int):
        flags = ClassicalRegister(8)                #flags[0..3] -> set A, flags[4..7] -> set B
        self.qc.add_register(flags)
        self._round(pos=pos, flags=flags)
        return flags

    def reset_qec(self, result):
        psi_full = result.data(0)["psi"]
        qr2 = QuantumRegister(9*(self.n+self.magic)+4, "q")
        self.qc = QuantumCircuit(qr2, self.cbit)        #same register as in __init__, otherwise the
                                                        #classical layout shifts and readout() misslices
        self.qc.add_register(self.qecc)
        self.qc.set_statevector(psi_full)
        del psi_full

    def _creg_bits(self, bitstring: str, creg):
        #the counts bitstring lists the classical registers in reverse order of when they were added,
        #so the position of a register depends on how many registers the circuit currently has
        start = 0
        for reg in reversed(self.qc.cregs):
            if reg == creg:                     #qiskit rebuilds the register objects, so compare by name/size
                return bitstring[start:start+reg.size]
            start += reg.size
        return None

    def _flag_round(self, pos: int):
        #one round of syndrome extraction with the 1-flag circuits --> (flags, syndrome, bitstring, result)
        flags = self.flagsyndrome(pos=pos)
        self.qc.save_statevector(label="psi")

        sim = AerSimulator(method="statevector", noise_model = self.noise_model)
        result = sim.run(self.qc, shots=1).result()
        counts = result.get_counts()                #one shot --> one bitstring

        bitstring = list(counts.keys())
        bitstring = [i.replace(" ","") for i in bitstring][0]

        return self._creg_bits(bitstring, flags), self._creg_bits(bitstring, self.qecc), bitstring, result

    def syndrome(self, pos: int):
        #unflagged syndrome measurement, measure only. The outcome is read out classically (so it can be
        #fed into a lookup table) instead of being corrected in the circuit like in qec().
        self._round(pos=pos)

        self.qc.save_statevector(label="psi")
        sim = AerSimulator(method="statevector", noise_model = self.noise_model)
        result = sim.run(self.qc, shots=1).result()
        bitstring = list(result.get_counts().keys())[0].replace(" ","")
        syndrome = self._creg_bits(bitstring, self.qecc)
        self.reset_qec(result)

        return syndrome

    def correct_9q(self, syndrome: str, pos: int):
        #E_min(s). syndrome = the 8 bit qecc slice:  qecc[7]qecc[6]qecc[5]qecc[4] | qecc[3]qecc[2]qecc[1]qecc[0]
        #                                           = set B syndrome              | set A syndrome
        #While the code sits in its original orientation set A is measured as X-Stabilizers, so its
        #syndrome asks for Z corrections and set B for X corrections. After an odd number of H_L the two
        #swap, the supports stay where they are.
        even = self.hadamards[pos]%2 == 0
        for q in CORR_A_9Q[syndrome[4:8]]:
            (self.qc.z if even else self.qc.x)(q+9*pos)
        for q in CORR_B_9Q[syndrome[0:4]]:
            (self.qc.x if even else self.qc.z)(q+9*pos)

    def flagcorrect(self, flags: str, pos: int):
        #3.Case of the protocol: a circuit C(g_i) flagged, so the syndrome is measured again with the
        #non-flag circuits. If the outcome s matches an element E of the flag error set E(g_i) that
        #element is applied, otherwise E_min(s).
        flag_i = len(flags)-1 - flags.index("1")        #the flags slice is printed high bit first
        syndrome = self.syndrome(pos=pos)

        qubits = CORR_FLAG_9Q[flag_i].get(syndrome)
        if qubits is None:
            self.correct_9q(syndrome, pos=pos)
            return
        even = self.hadamards[pos]%2 == 0
        for q in qubits:
            if flag_i < 4:                              #set A, measured as X-Stabilizers while even, so
                (self.qc.x if even else self.qc.z)(q+9*pos)     #the ancilla propagates X errors
            else:                                       #set B, measured as Z-Stabilizers while even, so
                (self.qc.z if even else self.qc.x)(q+9*pos)     #the ancilla propagates Z errors

    def flagFTec(self, pos: int):
        #based on https://arxiv.org/pdf/1708.02246 , protocol on page 5, "Flag 1-FTEC Protocol"  --> no need of postselection, but we need statevector simulation for classical decoding
        #The protocol repeats the flagged syndrome measurement until case 1, 2 or 3 applies. After the
        #second round one of the three always applies, so two flagged rounds are enough.

        flags_r1, syndrome_r1, bitstring, result = self._flag_round(pos=pos)

        if flags_r1.count("1") != 0:            #3.Case, falls direkt geflagged wird
            self.reset_qec(result)
            self.flagcorrect(flags_r1, pos=pos)
            return

        self.reset_qec(result)

        flags_r2, syndrome_r2, bitstring, result = self._flag_round(pos=pos)

        self.reset_qec(result)

        #3.Case: a circuit C(g_i) flagged --> non-flag syndrome measurement + flag error set E(g_i)
        if flags_r2.count("1") != 0:
            self.flagcorrect(flags_r2, pos=pos)
        #1.Case: no flags and the syndrome was repeated twice in a row --> E_min(s)
        elif syndrome_r1 == syndrome_r2:
            self.correct_9q(syndrome_r2, pos=pos)
        #2.Case: no flags but the two syndromes differ --> non-flag syndrome measurement + E_min(s)
        else:
            self.qec(pos=pos)

#######################################################################################################################################

    def readout(self, pos: int, shots: int, p = 0):
        code0 = ['000110101', '110110110', '110110101', '110000000', '000110110', '101101101', '011101101', '011011000', '011011011', '110000011', '000000000', '011101110', '101011011', '101101110', '000000011', '101011000']
        code1 = ['010100111', '010010001', '111111111', '001001010', '111001010', '001111111', '100010010', '111111100', '100100100', '100010001', '001001001', '010010010', '100100111', '111001001', '001111100', '010100100']

        
        read = ClassicalRegister(9)
        self.qc.add_register(read)
        # for i in range(9):
        #     self.qc.id(i+9*pos)
        if self.hadamards[pos]%2 == 0:
            for i in range(9):
                self.qc.measure(i+9*pos, read[8-i])
        else:
            self.qc.measure(0+9*pos, read[8-6])
            self.qc.measure(1+9*pos, read[8-3])
            self.qc.measure(2+9*pos, read[8-0])
            self.qc.measure(3+9*pos, read[8-7])
            self.qc.measure(4+9*pos, read[8-4])
            self.qc.measure(5+9*pos, read[8-1])
            self.qc.measure(6+9*pos, read[8-8])
            self.qc.measure(7+9*pos, read[8-5])
            self.qc.measure(8+9*pos, read[8-2])

        p_error = pauli_error([["X",p/2],["I",1-p],["Z",p/2]])
        p_error_2 = pauli_error([["XI",p/4],["IX",p/4],["II",1-p],["ZI",p/4],["IZ",p/4]])

        noise_model = NoiseModel()
        noise_model.add_all_qubit_quantum_error(p_error, ['x', "z", 'h', "id", "s", "sdg", "t", "tdg"])  # Apply to single-qubit gates
        noise_model.add_all_qubit_quantum_error(p_error_2, ['cx'])  # Apply to 2-qubit gates

        sim = AerSimulator()
        job = sim.run(self.qc, noise_model = noise_model, shots=shots)
        result = job.result()
        counts = result.get_counts()

        bitstring = list(counts.keys())
        # print(bitstring)
        bitstring = [i.replace(" ","") for i in bitstring]
        hmm = list(counts.values())

        bits = [i[:9] for i in bitstring]
        # print(bits)
        flags = [i[9:len(bitstring[0])-12] for i in bitstring]      #4 cbits + 8 qecc bits at the end

        if self.classical_ec:
            bits = self.__classical_error_correction__(bits)

        if self.postselection:
            for i in range(len(bits)):
                for j in code0:
                    if j == bits[i]:
                        bits[i] = 0
                        break
                if bits[i] != 0:
                    for j in code1:
                        if j == bits[i]:
                            bits[i] = 1
                            break
                if bits[i] != 1 and bits[i] != 0:
                    # print("AHA")
                    bits[i] = "post"
        else:
            for i in range(len(bits)):
                for j in code0:
                    if j == bits[i]:
                        bits[i] = 0
                        break
                if bits[i] != 0:
                    for j in code1:
                        if j == bits[i]:
                            bits[i] = 1
                            break
                if bits[i] != 1 and bits[i] != 0:
                    if np.random.rand() < 0.5:
                        bits[i] = 0
                    else:
                        bits[i] = 1

        for i in range(len(flags)):
            if flags[i].count("1") != 0:
                # print("AHAAA")
                bits[i] = "post"

        ones = 0
        zeros = 0
        err = 0

        for i in range(len(bits)):
            if bits[i] == 0:
                zeros += hmm[i]
            if bits[i] == 1:
                ones += hmm[i]
            if bits[i] == "post":
                err += hmm[i]
        
        ones = (ones/shots)
        zeros = (zeros/shots)
        err = (err/shots)

        self.ones = ones
        self.zeros = zeros
        self.post = err
        #return zeros, ones, err

    def magic_readout(self, pos: int, shots: int, p = 0):
        code0 = ['000110101', '110110110', '110110101', '110000000', '000110110', '101101101', '011101101', '011011000', '011011011', '110000011', '000000000', '011101110', '101011011', '101101110', '000000011', '101011000']
        code1 = ['010100111', '010010001', '111111111', '001001010', '111001010', '001111111', '100010010', '111111100', '100100100', '100010001', '001001001', '010010010', '100100111', '111001001', '001111100', '010100100']

        mapping = {0: 2, 1: 5, 2: 8, 3: 1, 4: 4, 5: 7, 6: 0, 7: 3, 8: 6}

        # def apply_mapping(s, mapping):
        #     result = [''] * len(s)
        #     for i, new_pos in mapping.items():
        #         result[new_pos] = s[i]
        #     return ''.join(result)

        # code0_h = [apply_mapping(s, mapping) for s in code0]
        # code1_h = [apply_mapping(s, mapping) for s in code1]
        
        read = ClassicalRegister(9)
        self.qc.add_register(read)
        # for i in range(9):
        #     self.qc.id(i+9*pos)
        if self.hadamards[pos]%2 == 0:
            for i in range(9):
                self.qc.measure(i+9*pos, read[8-i])
        else:
            self.qc.measure(0+9*pos, read[8-6])
            self.qc.measure(1+9*pos, read[8-3])
            self.qc.measure(2+9*pos, read[8-0])
            self.qc.measure(3+9*pos, read[8-7])
            self.qc.measure(4+9*pos, read[8-4])
            self.qc.measure(5+9*pos, read[8-1])
            self.qc.measure(6+9*pos, read[8-8])
            self.qc.measure(7+9*pos, read[8-5])
            self.qc.measure(8+9*pos, read[8-2])

        p_error = pauli_error([["X",p/2],["I",1-p],["Z",p/2]])
        p_error_2 = pauli_error([["XI",p/4],["IX",p/4],["II",1-p],["ZI",p/4],["IZ",p/4]])

        noise_model = NoiseModel()
        noise_model.add_all_qubit_quantum_error(p_error, ['x', "z", 'h', "id", "s", "sdg", "t", "tdg"])  # Apply to single-qubit gates
        noise_model.add_all_qubit_quantum_error(p_error_2, ['cx'])  # Apply to 2-qubit gates

        sim = AerSimulator()
        new_qc = transpile(self.qc, optimization_level=2)
        job = sim.run(new_qc, noise_model = noise_model, shots=shots)
        result = job.result()
        counts = result.get_counts()


        bitstring = list(counts.keys())
        print(bitstring)
        bitstring = [i.replace(" ","") for i in bitstring]
        hmm = list(counts.values())

        bits = [i[:9] for i in bitstring]
        # print(bits)
        flags = [i[9:len(bitstring[0])-12] for i in bitstring]      #4 cbits + 8 qecc bits at the end

        magic = [i[len(bitstring[0])-1] for i in bitstring]

        print(magic)

        if self.classical_ec:
            bits = self.__classical_error_correction__(bits)
 
        if self.postselection:
            for i in range(len(bits)):
                for j in code0:
                    if j == bits[i]:
                        bits[i] = 0
                        break
                if bits[i] != 0:
                    for j in code1:
                        if j == bits[i]:
                            bits[i] = 1
                            break
                if bits[i] != 1 and bits[i] != 0:
                    print(magic[i])
                    bits[i] = "post"
        else:
            for i in range(len(bits)):
                for j in code0:
                    if j == bits[i]:
                        bits[i] = 0
                        break
                if bits[i] != 0:
                    for j in code1:
                        if j == bits[i]:
                            bits[i] = 1
                            break
                if bits[i] != 1 and bits[i] != 0:
                    if np.random.rand() < 0.5:
                        bits[i] = 0
                    else:
                        bits[i] = 1

        for i in range(len(flags)):
            if flags[i].count("1") != 0:
                bits[i] = "post"

        ones = 0
        zeros = 0
        err = 0

        for i in range(len(bits)):
            if bits[i] == 0:
                zeros += hmm[i]
            if bits[i] == 1:
                ones += hmm[i]
            if bits[i] == "post":
                err += hmm[i]
        
        ones = (ones/shots)
        zeros = (zeros/shots)
        err = (err/shots)

        self.ones = ones
        self.zeros = zeros
        self.post = err
