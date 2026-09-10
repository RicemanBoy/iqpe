from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister
import numpy as np
import matplotlib.pyplot as plt

from qiskit_aer import AerSimulator

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
                    if code == "steane17":
                        self = Steane17q(1)

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
                        self.qec(0)

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

############################## [[17,1,5]] 4.8.8 Color Code ##############################
# The 8 faces of the square-octagon lattice (https://arxiv.org/pdf/2601.13313, Fig. 3b).
# The code is self-dual, so the X- and Z-type generators sit on the same faces.
# 36 CNOTs per basis. Qubits 0, 1, 12 are the three corners of the triangle,
# the bottom row of the lattice (0,1,2,4,6) is the logical operator used in x()/z().
STABS_17Q = [(5, 9, 11, 12),                # square
             (8, 10, 13, 15),               # square
             (3, 7, 14, 16),                # square
             (2, 6, 14, 16),                # octagon, cut by the boundary
             (4, 5, 6, 7, 8, 10, 11, 16),   # octagon
             (1, 4, 10, 13),                # octagon, cut by the boundary
             (8, 9, 11, 15),                # octagon, cut by the boundary
             (0, 2, 3, 14)]                  # octagon, cut by the boundary

# syndrome (bit i = face i) -> minimum weight correction. Covers every single- and
# twoqubit error; the 140 syndromes that only weight>=3 errors produce are left alone,
# so those runs stay outside the code space and get postselected in readout().
CORR_17Q = {
    0x01: (12,), 0x02: (1, 13), 0x03: (8, 11), 0x04: (0, 3),
    0x05: (5, 7), 0x08: (0, 2), 0x09: (5, 6), 0x0C: (0, 14),
    0x0D: (5, 16), 0x10: (1, 4), 0x11: (5,), 0x12: (1, 10),
    0x13: (8, 9), 0x14: (7,), 0x15: (7, 12), 0x18: (6,),
    0x19: (6, 12), 0x1C: (16,), 0x1D: (12, 16), 0x20: (1,),
    0x21: (1, 12), 0x22: (13,), 0x23: (5, 10), 0x24: (4, 7),
    0x26: (7, 10), 0x28: (4, 6), 0x2A: (6, 10), 0x2C: (4, 16),
    0x2E: (10, 16), 0x30: (4,), 0x31: (1, 5), 0x32: (10,),
    0x33: (5, 13), 0x34: (1, 7), 0x36: (7, 13), 0x38: (1, 6),
    0x3A: (6, 13), 0x3C: (1, 16), 0x3E: (13, 16), 0x40: (5, 11),
    0x41: (9,), 0x42: (15,), 0x43: (5, 8), 0x45: (7, 11),
    0x46: (7, 8), 0x49: (6, 11), 0x4A: (6, 8), 0x4D: (11, 16),
    0x4E: (8, 16), 0x50: (5, 9), 0x51: (11,), 0x52: (8,),
    0x53: (5, 15), 0x55: (7, 9), 0x56: (7, 15), 0x59: (6, 9),
    0x5A: (6, 15), 0x5D: (9, 16), 0x5E: (15, 16), 0x60: (8, 10),
    0x61: (1, 9), 0x62: (1, 15), 0x63: (9, 13), 0x70: (8, 13),
    0x71: (1, 11), 0x72: (1, 8), 0x73: (9, 10), 0x80: (0,),
    0x81: (0, 12), 0x84: (3,), 0x85: (3, 12), 0x88: (2,),
    0x89: (2, 12), 0x8C: (14,), 0x8D: (12, 14), 0x90: (2, 6),
    0x91: (0, 5), 0x94: (0, 7), 0x95: (3, 5), 0x98: (0, 6),
    0x99: (2, 5), 0x9C: (0, 16), 0x9D: (5, 14), 0xA0: (0, 1),
    0xA2: (0, 13), 0xA4: (1, 3), 0xA6: (3, 13), 0xA8: (1, 2),
    0xAA: (2, 13), 0xAC: (1, 14), 0xAE: (13, 14), 0xB0: (0, 4),
    0xB2: (0, 10), 0xB4: (3, 4), 0xB6: (3, 10), 0xB8: (2, 4),
    0xBA: (2, 10), 0xBC: (4, 14), 0xBE: (10, 14), 0xC1: (0, 9),
    0xC2: (0, 15), 0xC5: (3, 9), 0xC6: (3, 15), 0xC9: (2, 9),
    0xCA: (2, 15), 0xCD: (9, 14), 0xCE: (14, 15), 0xD1: (0, 11),
    0xD2: (0, 8), 0xD5: (3, 11), 0xD6: (3, 8), 0xD9: (2, 11),
    0xDA: (2, 8), 0xDD: (11, 14), 0xDE: (8, 14),
}

class Steane17q:
    def __init__(self, n = 1):
        self.n = n

        self.zeros = 0
        self.ones = 0
        self.preselected = 0
        self.post = 0
        self.err = False
        self.qec_counter = 0
        self.magiccounter = 0

        self.classical_ec = False
        self.postselection = True

        qr = QuantumRegister(17*n+2,"q")
        cbits = ClassicalRegister(1, "c")
        
        self.qc = QuantumCircuit(qr, cbits)
        
        anc = self.qc.num_qubits - 1
        ### Preparation from this paper: https://arxiv.org/pdf/2601.13313 : 23 CNOTs + 8 Hadamards ######
        for i in range(n):
            self.qc.h(0+17*i)
            self.qc.h(1+17*i)
            self.qc.h(2+17*i)
            self.qc.h(3+17*i)
            self.qc.h(4+17*i)
            self.qc.h(9+17*i)
            self.qc.h(10+17*i)
            self.qc.h(12+17*i)

            self.qc.cx(2+17*i, 16+17*i)
            self.qc.cx(3+17*i, 14+17*i)
            self.qc.cx(4+17*i, 7+17*i)
            self.qc.cx(10+17*i, 15+17*i)

            self.qc.cx(2+17*i, 6+17*i)
            self.qc.cx(4+17*i, 5+17*i)
            self.qc.cx(7+17*i, 8+17*i)
            self.qc.cx(10+17*i, 13+17*i)
            self.qc.cx(9+17*i, 15+17*i)

            self.qc.cx(0+17*i, 2+17*i)
            self.qc.cx(3+17*i, 7+17*i)
            self.qc.cx(5+17*i, 6+17*i)
            self.qc.cx(12+17*i, 9+17*i)
            self.qc.cx(8+17*i, 11+17*i)
            self.qc.cx(1+17*i, 4+17*i)

            self.qc.cx(0+17*i, 3+17*i)
            self.qc.cx(2+17*i, 14+17*i)
            self.qc.cx(7+17*i, 16+17*i)
            self.qc.cx(12+17*i, 5+17*i)
            self.qc.cx(9+17*i, 11+17*i)
            self.qc.cx(15+17*i, 8+17*i)
            self.qc.cx(4+17*i, 10+17*i)
            self.qc.cx(1+17*i, 13+17*i)

        ### Preparation from this paper: https://arxiv.org/pdf/2402.17761 : 25 CNOTs + 8 Hadamards ######
        # for i in range(n):
        #     self.qc.h(6+17*i)
        #     self.qc.h(1+17*i)
        #     self.qc.h(2+17*i)
        #     self.qc.h(16+17*i)
        #     self.qc.h(4+17*i)
        #     self.qc.h(11+17*i)
        #     self.qc.h(12+17*i)
        #     self.qc.h(14+17*i)

        #     self.qc.cx(1+17*i, 13+17*i)
        #     self.qc.cx(12+17*i, 0+17*i)
        #     self.qc.cx(0+17*i, 3+17*i)
        #     self.qc.cx(4+17*i, 5+17*i)
        #     self.qc.cx(6+17*i, 9+17*i)
        #     self.qc.cx(11+17*i, 10+17*i)

        #     self.qc.cx(4+17*i, 1+17*i)
        #     self.qc.cx(16+17*i, 7+17*i)
        #     self.qc.cx(14+17*i, 6+17*i)
        #     self.qc.cx(9+17*i, 8+17*i)

        #     self.qc.cx(2+17*i, 13+17*i)
        #     self.qc.cx(0+17*i, 15+17*i)

        #     self.qc.cx(16+17*i, 15+17*i)
        #     self.qc.cx(13+17*i, 12+17*i)
        #     self.qc.cx(3+17*i, 8+17*i)

        #     self.qc.cx(0+17*i, 5+17*i)
        #     self.qc.cx(14+17*i, 16+17*i)

        #     self.qc.cx(8+17*i, 1+17*i)
        #     self.qc.cx(15+17*i, 11+17*i)
        #     self.qc.cx(1+17*i, 3+17*i)
        #     self.qc.cx(8+17*i, 6+17*i)
        #     self.qc.cx(14+17*i, 11+17*i)

        #     self.qc.cx(2+17*i, 0+17*i)
        #     self.qc.cx(10+17*i, 6+17*i)

        #     self.qc.cx(10+17*i, 7+17*i)
            
        self.qecc = ClassicalRegister(16)
        #self.qc.add_register(self.qecc)

    def x(self, pos = 0):
        self.qc.x(0+17*pos)
        self.qc.x(1+17*pos)
        self.qc.x(2+17*pos)
        self.qc.x(4+17*pos)
        self.qc.x(6+17*pos)

    def z(self, pos = 0):
        self.qc.z(0+17*pos)
        self.qc.z(1+17*pos)
        self.qc.z(2+17*pos)
        self.qc.z(4+17*pos)
        self.qc.z(6+17*pos)

    def h(self, pos = 0):
        for i in range(17):
            self.qc.h(i+17*pos)

    def s(self, pos = 0):
        for i in range(17):
            self.qc.s(i+17*pos)

    def sdg(self, pos = 0):
        for i in range(17):
            self.qc.sdg(i+17*pos)

    def t_anc(self, pos = 0):
        anc = self.qc.num_qubits - 1
        self.qc.reset(anc)
        self.qc.h(anc)
        self.qc.t(anc)

        self.qc.cx(0+17*pos, anc)
        self.qc.cx(1+17*pos, anc)
        self.qc.cx(2+17*pos, anc)
        self.qc.cx(4+17*pos, anc)
        self.qc.cx(6+17*pos, anc)

        self.qc.measure(anc, 0)

        with self.qc.if_test((0,1)):
            for i in range(17):
                self.qc.s(i+17*pos)

    def t(self, pos = 0):
        for i in range(17):
            self.qc.h(i+17*pos)
            self.qc.h(i+17*pos)
            self.qc.s(i+17*pos)

        anc = self.qc.num_qubits - 1
        self.qc.reset(anc)
        self.qc.append(h_ideal,[anc])
        self.qc.append(t_ideal,[anc])

        self.qc.append(cx_ideal, [anc, 0+17*pos])
        self.qc.append(cx_ideal, [anc, 1+17*pos])
        self.qc.append(cx_ideal, [anc, 2+17*pos])
        self.qc.append(cx_ideal, [anc, 4+17*pos])
        self.qc.append(cx_ideal, [anc, 6+17*pos])

        self.qc.measure(anc, 0)

        with self.qc.if_test((0,1)):
            for i in range(17):
                self.qc.s(i+17*pos)

        for i in range(17):
            self.qc.sdg(i+17*pos)
            self.qc.h(i+17*pos)
            self.qc.h(i+17*pos)

        if self.err:
            self.qec(pos=pos)

    def tdg_anc(self, pos = 0):
            anc = self.qc.num_qubits - 1
            self.qc.reset(anc)
            self.qc.h(anc), self.qc.tdg(anc)
    
            self.qc.cx(0+17*pos, anc)
            self.qc.cx(1+17*pos, anc)
            self.qc.cx(2+17*pos, anc)
            self.qc.cx(4+17*pos, anc)
            self.qc.cx(6+17*pos, anc)
    
            self.qc.measure(anc, 0)
    
            for i in range(17):
                with self.qc.if_test((0,1)):
                    self.qc.sdg(i+17*pos)

    def tdg(self, pos = 0):
        for i in range(17):
            self.qc.h(i+17*pos)
            self.qc.h(i+17*pos)
            self.qc.s(i+17*pos)

        anc = self.qc.num_qubits - 1
        self.qc.reset(anc)
        self.qc.append(h_ideal,[anc])
        self.qc.append(tdg_ideal,[anc])

        self.qc.append(cx_ideal, [anc, 0+17*pos])
        self.qc.append(cx_ideal, [anc, 1+17*pos])
        self.qc.append(cx_ideal, [anc, 2+17*pos])
        self.qc.append(cx_ideal, [anc, 4+17*pos])
        self.qc.append(cx_ideal, [anc, 6+17*pos])

        self.qc.measure(anc, 0)

        with self.qc.if_test((0,1)):
            for i in range(17):
                self.qc.sdg(i+17*pos)

        for i in range(17):
            self.qc.sdg(i+17*pos)
            self.qc.h(i+17*pos)
            self.qc.h(i+17*pos)

        if self.err:
            self.qec(pos=pos)
        
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

    def cu_ramsey(self, gate: list):
        self.u2(0, gate=gate)
        if self.err:
            self.qec(0)
        self.u2(0, gate=gate)

    def qec(self, pos = 0):
        self.qec_counter += 1
        if self.qecc not in self.qc.cregs:              # qecc is not added in __init__
            self.qc.add_register(self.qecc)
        anc = self.qc.num_qubits - 1
        self.qc.reset(anc)
        ##################################Z-Stabilizers##########################################
        for i, face in enumerate(STABS_17Q):
            for q in face:
                self.qc.cx(q+17*pos, anc)
            self.qc.id(anc)
            self.qc.measure(anc, self.qecc[i])
            self.qc.reset(anc)
            self.qc.id(anc)
        ##################################X-Stabilizers##############################################
        for i, face in enumerate(STABS_17Q):
            self.qc.h(anc)
            for q in face:
                self.qc.cx(anc, q+17*pos)
            self.qc.h(anc)
            self.qc.id(anc)
            self.qc.measure(anc, self.qecc[8+i])
            self.qc.reset(anc)
            self.qc.id(anc)
        ##################################Bitflip Error correction##############################################
        zsyn = expr.bit_and(self.qecc, 0x00FF)
        for syndrome, qubits in CORR_17Q.items():
            with self.qc.if_test(expr.equal(zsyn, syndrome)):
                for q in qubits:
                    self.qc.x(q+17*pos)
        ##################################Phaseflip Error correction##############################################
        xsyn = expr.bit_and(self.qecc, 0xFF00)
        for syndrome, qubits in CORR_17Q.items():
            with self.qc.if_test(expr.equal(xsyn, syndrome << 8)):
                for q in qubits:
                    self.qc.z(q+17*pos)

    def qec_ideal(self, pos = 0):
        self.qec_counter += 1
        if self.qecc not in self.qc.cregs:              # qecc is not added in __init__
            self.qc.add_register(self.qecc)
        anc = self.qc.num_qubits - 1
        self.qc.reset(anc)
        ##################################Z-Stabilizers##########################################
        for i, face in enumerate(STABS_17Q):
            for q in face:
                self.qc.append(cx_ideal, [anc, q+17*pos])
            self.qc.measure(anc, self.qecc[i])
            self.qc.reset(anc)
        ##################################X-Stabilizers##############################################
        for i, face in enumerate(STABS_17Q):
            self.qc.append(h_ideal, [anc])
            for q in face:
                self.qc.append(cx_ideal, [q+17*pos, anc])
            self.qc.append(h_ideal, [anc])
            self.qc.measure(anc, self.qecc[8+i])
            self.qc.reset(anc)
        ##################################Bitflip Error correction##############################################
        zsyn = expr.bit_and(self.qecc, 0x00FF)
        for syndrome, qubits in CORR_17Q.items():
            with self.qc.if_test(expr.equal(zsyn, syndrome)):
                for q in qubits:
                    self.qc.append(x_ideal, [q+17*pos])
        ##################################Phaseflip Error correction##############################################
        xsyn = expr.bit_and(self.qecc, 0xFF00)
        for syndrome, qubits in CORR_17Q.items():
            with self.qc.if_test(expr.equal(xsyn, syndrome << 8)):
                for q in qubits:
                    self.qc.append(z_ideal, [q+17*pos])

    def readout(self, pos: int, shots: int, p = 0):
        p_error = pauli_error([["X",p/2],["I",1-p],["Z",p/2]])
        p_error_2 = pauli_error([["XI",p/4],["IX",p/4],["II",1-p],["ZI",p/4],["IZ",p/4]])

        noise_model = NoiseModel()
        noise_model.add_all_qubit_quantum_error(p_error, ['x', "z", 'h', "s", "sdg", "id", "t", "tdg"])  # Apply to single-qubit gates
        noise_model.add_all_qubit_quantum_error(p_error_2, ['cx'])  # Apply to 2-qubit gates
        
        read = ClassicalRegister(17)
        self.qc.add_register(read)

        for i in range(17):
            self.qc.id(i+7*pos)
            self.qc.measure(i+17*pos,read[16-i])

        # self.qc = transpile(self.qc, optimization_level=1)

        code0 = ['11101001100000011', '00111100101100001', '00110011011101000', '11011010010100001', '01010110111000100', '00101001001110110', '00000000110100010', '11111000010100100', '00000000000000000', '01011101110110111', '10001000010011111', '01100101010001110', '01101110000010101', '00101101000101110', '10010010000000001', '11001111101011110', '01001000111101010', '00010101010110101', '01001100000010000', '01001000010100000', '11010001110011000', '10101110110001000', '01101010001001101', '00111000010011011', '11101101101011011', '01111111000010000', '11111100011111100', '00110011000000000', '10001000100111101', '00000100111111010', '11101001001001001', '11111100101011110', '10000111111111110', '11000000101110101', '01100001110011100', '11111100000010100', '10101010010011010', '11111100110110110', '01000111001100011', '00100110100010111', '00110011101001010', '10110100001011100', '10000011000000100', '01100001101110100', '00111100011000011', '10001000001110111', '10101110101100000', '01011001100000111', '00000000101001010', '01110100111000001', '10011001010011010', '01000111100101001', '11001011010100100', '00111100000101011', '10000011011101100', '11110011000111111', '00111000111010001', '10110000110100110', '11001011111101110', '01111111101011010', '00001111000101011', '10111011100111101', '11100110010001010', '10011101101100000', '00110111100010010', '01000011011010011', '11111000100000110', '00011010100111100', '00110111001011000', '10001100000101111', '10100001000000001', '01110000000111011', '00010101100010111', '11100010110011000', '11011110110110011', '11001111011111100', '01010010000111110', '11010101111000000', '01100101100101100', '10101010001110010', '11100110100101000', '11000100100101101', '01010110100101100', '10100101001011001', '00101001111010100', '10000011101001110', '01000111111000001', '01100001000111110', '10111011010011111', '01110100010001011', '11100010000111010', '00111100110001001', '01011001010100101', '00110111111111010', '10111111101100101', '01101110011111101', '00000000011101000', '10011101011000010', '00100010000000101', '01001000001001000', '01110000101110001', '10110100010110100', '10011101110001000', '00010001000000101', '11000000110011101', '00101101110001100', '11110111010001111', '01111011100000010', '00011010001110110', '10000011110100110', '00011110011000110', '00001011001110011', '01101010111101111', '10101110011000010', '00101101101100100', '00010001011101101', '00100110001011101', '00001011010011011', '01011001111101111', '01000011110011001', '01010010011010110', '01000011000111011', '01110000011010011', '00101001010011110', '00001011111010001', '00011010010011110', '01101110101011111', '11010101100101000', '00010101111111111', '10010010110100011', '10100101100010011', '10000111100010110', '10001100101100101', '00000100001011000', '11000100111000101', '10001100011000111', '01011101011111101', '10101110000101010', '11010101010001010', '11000100001100111', '01111111110110010', '00000100010110000', '10100101111111011', '10011001111010000', '11010101001100010', '00110011110100010', '10111111110001101', '10101010100111000', '11110011011010111', '01001100011111000', '10100001101001011', '01100101111000100', '00011110110001100', '00010101001011101', '10001000111010101', '10000111010110100', '10111011111010101', '01101110110110111', '10110000000000100', '00001111101100001', '10010010101001011', '01100101001100110', '00110111010110000', '11001111110110110', '01010110001100110', '00101001100111100', '11100010101110000', '11011010100000011', '11001011100000110', '00100010110100111', '00111000100111001', '00001111011000011', '10010010011101001', '00010001101001111', '00101101011000110', '01011101101011111', '11011110101011011', '10010110001011001', '00010001110100111', '11110011110011101', '10011001100111000', '10110100111111110', '11010001000111010', '11111000111101110', '01110000110011001', '00011010111010100', '11011010001001001', '01000111010001011', '11101001111101011', '01001100110110010', '10011101000101010', '11101101011111001', '10111011001110111', '01111011001001000', '01011101000010101', '10001100110001101', '10100001110100011', '01010010110011100', '11110111111000101', '01110100001100011', '11011010111101011', '10101010111010000', '01101010100000111', '10100101010110001', '00100110111111111', '00001011100111001', '10010110100010011', '01101010010100101', '11100010011010010', '11011110011111001', '01011001001001101', '11010001011010010', '01111011111101010', '11100110111000000', '11101001010100001', '11110111001100111', '10111111000101111', '01000011101110001', '01010010101110100', '00001111110001001', '11011110000010001', '11010001101110000', '11000000000111111', '01110100100101001', '00011110000101110', '00100010011101101', '00100110010110101', '10100001011101001', '11101101110110011', '01111111011111000', '11000000011010111', '10000111001011100', '01001100101011010', '11000100010001111', '01111011010100000', '00000100100010010', '00111000001110011', '10110000101001110', '11110111100101101', '10011001001110010', '11110011101110101', '11001111000010100', '10110000011101100', '00011110101100100', '11101101000010001', '01001000100000010', '11100110001100010', '10110100100010110', '10111111011000111', '01010110010001110', '11001011001001100', '00100010101001111', '01100001011010110', '10010110111111011', '10010110010110001', '11111000001001100']
        code1 = ['01100010111010101', '01111000101001011', '00110100011111001', '00110000010100001', '00100101101011110', '10100010100000010', '01100010010011111', '11011101111111010', '11100001111010001', '00000011111101011', '00101010000111111', '01011110100010110', '10100110110110010', '00001000101110000', '10011010101110001', '01011110111111110', '00011101001100111', '11010110110001001', '10110111101011111', '00100101000010100', '10011010000111011', '10000100000010101', '10011010110011001', '01110111011000010', '11010010100111001', '01011110001011100', '01110111110001000', '10010101101011010', '01001111100010011', '00100001100000110', '00101010110011101', '10110011100000111', '01000000111010000', '01001011101001011', '01011010110100110', '10111100111000100', '00100101110110110', '00110100101011011', '01111100100010011', '01111000000000001', '01111100111111011', '10101101010001011', '11000011010011110', '10101101100101001', '00010010010100100', '01100010001110111', '10100010001001000', '11101010101001010', '00101110111000101', '10000000111101111', '00000011100000011', '10001011101110100', '10001111100101100', '01101101100010110', '11100101000101011', '01010001001110111', '00001000110011000', '00001000011010010', '10000000010100101', '00011101100101101', '00000111000010001', '01001111001011001', '11111011101001111', '01001011110100011', '10010001010100000', '01101001110100110', '10001011000111110', '01000000010011010', '00101110010001111', '00011001011010111', '10101101001100011', '01110011111010000', '10011110001100011', '01001011000000001', '10101001101110001', '01110111101100000', '10011110111000001', '01010101000101111', '11111011000000101', '11110100101100100', '01011110010110100', '10110111110110111', '00000011010100001', '11110100000101110', '11001100100010111', '11010010001110011', '10101101111000001', '11010110101100001', '00011101111000101', '00000111110110011', '00100001010100100', '01011010011101100', '10110111000010101', '00110100110110011', '11000111110001100', '01010101110001101', '11100001100111001', '10000100011111101', '11010010010011011', '10011110100101001', '11011101001011000', '01000100110001000', '11111011011101101', '11111111010110101', '11101010000000000', '10010101110110010', '11001000101001111', '00111011101110000', '11110100011000110', '01101001011101100', '10110111011111101', '01111000011101001', '00101110001100111', '10101001110011001', '00010010111101110', '10100110000010000', '00010110101011110', '00111111010001010', '01000100101100000', '10010001100000010', '11001000110100111', '11110000100111100', '11000011111010100', '11111111001011101', '01110011001110010', '10000000100000111', '10001011011010110', '01101101010110100', '10010101011111000', '11100001010011011', '00111111100101000', '00111011011010010', '00000011001001001', '11000011100111100', '00011101010001111', '11011101100010010', '00111111111000000', '01110011100111000', '01100110110001101', '11101010011101000', '00001100001100010', '10110011010100101', '10011010011010011', '00110100000010001', '11011101010110000', '00011001101110101', '01010101011000111', '00000111101011011', '01000100011000010', '10001011110011100', '11010010111010001', '01101101001011100', '00111111001100010', '01010001010011111', '11001000000000101', '10111000110011100', '11001100001011101', '11000011001110110', '01100010100111101', '00010010100000110', '00001100100101000', '01011010000000100', '11101110001011000', '01000000001110010', '00001000000111010', '11101110010110000', '11111111100010111', '11100101110001001', '11111111111111111', '00010010001001100', '10111100001100110', '00110000001001001', '11011001110100010', '01100110011000111', '11010110011000011', '00011001110011101', '10100010111101010', '11111011110100111', '10010001111101010', '11100101101100001', '11101110100010010', '11110000010011110', '10100110011111000', '01111000110100011', '01101001101001110', '10100110101011010', '00010110000010100', '11110000001110110', '01001011011101001', '10000100110110111', '10010001001001000', '11001000011101101', '00110000100000011', '00011001000111111', '00000111011111001', '01100110101100101', '00100101011111100', '10111000011010110', '10110011111101111', '11100101011000011', '10111100010001110', '01111100001011001', '00111011000111010', '10001111010001110', '11011001000000000', '11001100010110101', '01110111000101010', '00001100111000000', '10011110010001011', '10111100100101100', '00001100010001010', '11011001011101000', '00100001111101110', '11100001001110011', '01111100010110001', '10000100101011111', '10010101000010000', '01001111010110001', '10001111111000100', '11110100110001100', '01101001000000100', '01010101101100101', '10110011001001101', '10101001011010011', '00010110110110110', '01100110000101111', '11000111101100100', '11011001101001010', '01010001100111101', '00101110100101101', '01000100000101010', '10101001000111011', '00100001001001100', '10111000101110100', '10001111001100110', '11000111000101110', '01010001111010101', '10000000001001101', '11101010110100010', '01001111111111011', '00110000111101011', '00101010011010111', '00101010101110101', '11001100111111111', '11110000111010100', '10111000000111110', '00111011110011000', '01101101111111110', '11101110111111010', '11010110000101011', '01110011010011010', '01011010101001110', '11000111011000110', '00010110011111100', '01000000100111000', '10100010010100000']
                
        sim = AerSimulator()
        job = sim.run(self.qc, shots=shots, noise_model=noise_model)
        result = job.result()
        counts = result.get_counts()

        bitstring = list(counts.keys())
        bitstring = [i.replace(" ","") for i in bitstring]
        hmm = list(counts.values())

        bits = [i[:17] for i in bitstring]

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
            if bits[i] != 1 and bits[i] != 0 and bits[i] != "pre":
                if self.postselection:
                    bits[i] = "post"
                else:
                    if np.random.rand() < 0.5:
                        bits[i] = 0
                    else:
                        bits[i] = 1
    
        for i in range(len(bits)):
            if bits[i] == 0:
                self.zeros += hmm[i]
            if bits[i] == 1:
                self.ones += hmm[i]
            if bits[i] == "post":
                self.post += hmm[i]
            if bits[i] == "pre":
                self.preselected += hmm[i]

        self.ones = (self.ones/shots)
        self.zeros = (self.zeros/shots)
        self.post = (self.post/shots)
        self.preselected = (self.preselected/shots)
        return counts
