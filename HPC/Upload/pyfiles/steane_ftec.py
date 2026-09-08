from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister
import numpy as np
import matplotlib.pyplot as plt

from qiskit_aer import AerSimulator

from qiskit_aer.noise import (NoiseModel, pauli_error)

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
    hmm["measure"] = 0
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
    for _ in range(k):
        for o in range(7):
            gatecount = 0
            bitstring = ""
            rots = []
            for t in range(iter):
                rots = [k*0.5 for k in rots]
                counter = 0
                while True:
                    if code == "steane":
                        self = Steane7q(1, 1, noise=noise)
                    else:
                        print("Error hahaha")
                        return

                    self.err = qec
                    self.h(pos=0)
                    #############################
                    for _ in range(2**(iter-t-1)):
                        self.cu_ramsey(a[2*o+1])
                    ###############################
                    for l in rots:
                        if l == 0.25:
                            self.sdg(pos=0)
                        if l == 0.125:
                            self.tdg(pos=0)
                    self.h(pos=0)
                    if self.err:
                        self.qec_ft(0)

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
                    print("Angle {}, {}% errorrate, Iteration {}: {} Repetition".format(2*o+1, noise*100, t, counter))
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

# Steane [[7,1,3]] lookup table. Key = 3 bit sub-syndrome of the qecc register, read MSB..LSB.
# Z-Stabilizers  qecc[2]qecc[1]qecc[0]  -> X correction
# X-Stabilizers  qecc[5]qecc[4]qecc[3]  -> Z correction
CORR_7Q = {
    "000": (),    "001": (3,), "010": (1,), "011": (5,),
    "100": (0,),  "101": (4,), "110": (2,), "111": (6,),
}

STABS_7Q = ((0,2,4,6), (1,2,5,6), (3,4,5,6))    #order matches qecc[2],qecc[1],qecc[0] and qecc[5],qecc[4],qecc[3]

# Flag error sets E(g_i) of the 1-flag circuits in flagsyndrome(), needed for case 3 of the Flag 1-FTEC
# Protocol. Key = the 6 bit qecc slice of the following unflagged syndrome measurement, value = the qubits
# the flagged fault hit. flags[0..2] are the Z-Stabilizer circuits (Z errors on the data), flags[3..5] the
# X-Stabilizer circuits (X errors). Entries that agree with E_min(s) up to a stabilizer are left out and
# fall through to CORR_7Q, as in the "otherwise apply E_min(s)" of the protocol.
CORR_FLAG_7Q = {
    0: {"010000": (4,6),   "100000": (2,4,6), "111000": (6,)},      #C(Z{0,2,4,6})
    1: {"010000": (2,5,6), "100000": (5,6),   "111000": (6,)},      #C(Z{1,2,5,6})
    2: {"001000": (4,5,6), "100000": (5,6),   "111000": (6,)},      #C(Z{3,4,5,6})
    3: {"000010": (4,6),   "000100": (2,4,6), "000111": (6,)},      #C(X{0,2,4,6})
    4: {"000010": (2,5,6), "000100": (5,6),   "000111": (6,)},      #C(X{1,2,5,6})
    5: {"000001": (4,5,6), "000100": (5,6),   "000111": (6,)},      #C(X{3,4,5,6})
}

class Steane7q:
    def __init__(self, n: int, magic = 1, noise = 0):
        self.n = n

        self.zeros = 0
        self.ones = 0
        self.preselected = 0
        self.post = 0
        self.err = False
        self.qec_counter = 0
        self.magiccounter = 0

        self.preselection_flag = False      #change to true if the initialization flag striked!

        self.classical_ec = False
        self.postselection = True

        self.noise_model = self.__noise_model__(noise,0)

        qr = QuantumRegister(7*(n+magic)+2,"q")
        cbits = ClassicalRegister(3, "c")
        self.cbits = cbits
        
        self.qc = QuantumCircuit(qr, cbits)
        
        anc = self.qc.num_qubits - 1

        for i in range(7*n):
            self.qc.id(i)

        for i in range(n):
            self.qc.h(1+7*i)
            self.qc.h(2+7*i)
            self.qc.h(3+7*i)

            self.qc.cx(1+7*i,0+7*i)
            self.qc.cx(3+7*i,5+7*i)

            self.qc.cx(2+7*i,6+7*i)

            self.qc.cx(1+7*i,4+7*i)

            self.qc.cx(2+7*i,0+7*i)
            self.qc.cx(3+7*i,6+7*i)

            self.qc.cx(1+7*i,5+7*i)

            self.qc.cx(6+7*i,4+7*i)

            self.qc.cx(0+7*i,anc)
            self.qc.cx(5+7*i,anc)
            self.qc.cx(6+7*i,anc)

            self.qc.id(anc)
            self.qc.measure(anc,cbits[i+1])

        self.qecc = ClassicalRegister(6)
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
    
    def id(self, pos: int):
        for i in range(7):
            self.qc.id(i+7*pos)

    def x(self, pos: int):
        self.qc.x(0+7*pos)
        self.qc.x(1+7*pos)
        self.qc.x(2+7*pos)

    def z(self, pos: int):
        self.qc.z(0+7*pos)
        self.qc.z(1+7*pos)
        self.qc.z(2+7*pos)

    def h(self, pos: int):
        for i in range(7):
            self.qc.h(i+7*pos)

    def s(self, pos: int):
        for i in range(7):
            self.qc.sdg(i+7*pos)

    def cz(self, control: int, target: int):
        self.h(pos=control)
        self.cnot(control=target, target=control)
        self.h(pos=control)

    def cnot(self, control: int, target: int):
        for i in range(7):
            self.qc.cx(i+7*control, i+7*target)

    def sdg(self, pos: int):
        for i in range(7):
            self.qc.s(i+7*pos)

    def t(self, pos: int):
        self.magiccounter += 1
        self.h(pos=pos)
        self.sdg(pos=pos)
        self.h(pos=pos)
        # state_inj = ClassicalRegister(1)
        # self.qc.add_register(state_inj)

        anc = self.qc.num_qubits - 1
        ancc = anc - 1

        for i in range(7):
            self.qc.reset(i+7*self.n)

        self.qc.append(h_ideal,[0+7*self.n])
        self.qc.append(h_ideal,[1+7*self.n])
        self.qc.ry(np.pi/4,2+7*self.n)
        self.qc.append(h_ideal,[3+7*self.n])

        self.qc.append(cx_ideal, [4+7*self.n, 2+7*self.n])
        self.qc.append(cx_ideal, [6+7*self.n, 0+7*self.n])
        
        self.qc.append(cx_ideal, [5+7*self.n, 3+7*self.n])

        self.qc.append(cx_ideal, [5+7*self.n, 2+7*self.n])

        self.qc.append(cx_ideal, [4+7*self.n, 0+7*self.n])
        self.qc.append(cx_ideal, [6+7*self.n, 1+7*self.n])

        self.qc.append(cx_ideal, [2+7*self.n, 0+7*self.n])

        self.qc.append(cx_ideal, [5+7*self.n, 1+7*self.n])

        self.qc.append(cx_ideal, [2+7*self.n, 1+7*self.n])
        self.qc.append(cx_ideal, [4+7*self.n, 3+7*self.n])
        self.qc.append(cx_ideal, [6+7*self.n, 3+7*self.n])
        #################################Controlled Hadamards##########################################
        self.qc.reset(ancc)
        self.qc.append(h_ideal,[ancc])
        for i in range(7):
            self.qc.ry(-np.pi/4,6-i+self.n*7)
            self.qc.cz(ancc,6-i+self.n*7)
            self.qc.ry(np.pi/4,6-i+self.n*7)
        self.qc.append(h_ideal,[ancc])
        # self.qc.measure(ancc, state_inj[0])
        ########################Controlled-Y Gate####################################################
        self.sdg(pos=pos)
        for i in range(7):
            self.qc.cx(i+7*self.n,i+7*pos)
        self.s(pos=pos)
        self.qc.reset(anc-1)
        #############################Measure logical state of the magic state for state injection#############################
        self.sdg(pos=self.n)
        self.h(pos=self.n)
        for i in range(7):
            self.qc.cx(i+self.n*7, anc-1)
        self.qc.measure(anc-1,0)
        #################################Apply conditioned Ry(pi/2) onto the Target###########################
        for i in range(7):
            with self.qc.if_test((0,1)):
                self.qc.h(i+7*pos)
        for i in range(3):
            with self.qc.if_test((0,1)):
                self.qc.x(i+7*pos)
        self.h(pos=pos)
        self.s(pos=pos)
        self.h(pos=pos)
        if self.err:
            self.qec_ft(pos=pos)

    def t_cheat(self, pos: int):
        self.qc.cx(0+7*pos, 2+7*pos)
        self.qc.cx(1+7*pos, 2+7*pos)
        self.qc.t(2+7*pos)
        self.qc.cx(0+7*pos, 2+7*pos)
        self.qc.cx(1+7*pos, 2+7*pos)

    def tdg_cheat(self, pos: int):
        self.qc.cx(0+7*pos, 2+7*pos)
        self.qc.cx(1+7*pos, 2+7*pos)
        self.qc.tdg(2+7*pos)
        self.qc.cx(0+7*pos, 2+7*pos)
        self.qc.cx(1+7*pos, 2+7*pos)

    def tdg(self, pos: int):
        self.magiccounter += 1
        self.h(pos=pos)
        self.sdg(pos=pos)
        self.h(pos=pos)
        # state_inj = ClassicalRegister(1)
        # self.qc.add_register(state_inj)

        anc = self.qc.num_qubits - 1
        ancc = anc - 1

        for i in range(7):
            self.qc.reset(i+7*self.n)

        self.qc.append(h_ideal,[0+7*self.n])
        self.qc.append(h_ideal,[1+7*self.n])
        self.qc.ry(np.pi/4,2+7*self.n)
        self.qc.append(h_ideal,[3+7*self.n])

        self.qc.append(cx_ideal, [4+7*self.n, 2+7*self.n])
        self.qc.append(cx_ideal, [6+7*self.n, 0+7*self.n])
        
        self.qc.append(cx_ideal, [5+7*self.n, 3+7*self.n])

        self.qc.append(cx_ideal, [5+7*self.n, 2+7*self.n])

        self.qc.append(cx_ideal, [4+7*self.n, 0+7*self.n])
        self.qc.append(cx_ideal, [6+7*self.n, 1+7*self.n])

        self.qc.append(cx_ideal, [2+7*self.n, 0+7*self.n])

        self.qc.append(cx_ideal, [5+7*self.n, 1+7*self.n])

        self.qc.append(cx_ideal, [2+7*self.n, 1+7*self.n])
        self.qc.append(cx_ideal, [4+7*self.n, 3+7*self.n])
        self.qc.append(cx_ideal, [6+7*self.n, 3+7*self.n])
        #################################Controlled Hadamards##########################################
        self.qc.reset(ancc)
        self.qc.append(h_ideal,[ancc])
        for i in range(7):
            self.qc.ry(-np.pi/4,6-i+self.n*7)
            self.qc.cz(ancc,6-i+self.n*7)
            self.qc.ry(np.pi/4,6-i+self.n*7)
        self.qc.append(h_ideal,[ancc])
        # self.qc.measure(ancc, state_inj[0])
        ########################Controlled-Y Gate####################################################
        self.sdg(pos=pos)
        for i in range(7):
            self.qc.cx(i+7*self.n,i+7*pos)
        self.s(pos=pos)
        self.qc.reset(anc-1)
        #############################Measure logical state of the magic state for state injection#############################
        self.sdg(pos=self.n)
        self.h(pos=self.n)
        for i in range(7):
            self.qc.cx(i+self.n*7, anc-1)
        self.qc.measure(anc-1,0)
        #################################Apply conditioned Ry(pi/2) onto the Target###########################
        for i in range(3):
            with self.qc.if_test((0,0)):
                self.qc.x(i+7*pos)
        for i in range(7):
            with self.qc.if_test((0,0)):
                self.qc.h(i+7*pos)
        self.h(pos=pos)
        self.s(pos=pos)
        self.h(pos=pos)
        if self.err:
            self.qec_ft(pos=pos)

    def cs(self, control: int, target: int):
        self.t(pos=control)
        self.t(pos=target)
        self.cnot(control=control, target=target)
        self.tdg(pos=target)
        self.cnot(control=control, target=target)

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
        self.u2(0, gate=gate)

    def qec_ft(self, pos: int):
        self.qec_counter += 1
        flags = ClassicalRegister(6)
        self.qc.add_register(flags)
        anc = self.qc.num_qubits - 1
        ancc = anc - 1
        self.qc.reset(anc), self.qc.reset(ancc)
        ##################################Z-Stabilizers##########################################
        self.qc.h(ancc)
        self.qc.cx(0+7*pos, anc)
        self.qc.cx(ancc,anc)
        self.qc.cx(2+7*pos, anc)
        self.qc.cx(4+7*pos, anc)
        self.qc.cx(ancc,anc)
        self.qc.cx(6+7*pos, anc)

        self.qc.id(anc), self.qc.h(ancc), self.qc.id(ancc)
        self.qc.measure(anc, self.qecc[2]), self.qc.measure(ancc, flags[0])
        self.qc.reset(anc), self.qc.reset(ancc)
        self.qc.id(anc), self.qc.id(ancc)

        self.qc.h(ancc)
        self.qc.cx(1+7*pos, anc)
        self.qc.cx(ancc, anc)
        self.qc.cx(2+7*pos, anc)
        self.qc.cx(5+7*pos, anc)
        self.qc.cx(ancc, anc)
        self.qc.cx(6+7*pos, anc)

        self.qc.id(anc), self.qc.h(ancc), self.qc.id(ancc)
        self.qc.measure(anc, self.qecc[1]), self.qc.measure(ancc, flags[1])
        self.qc.reset(anc), self.qc.reset(ancc)
        self.qc.id(anc), self.qc.id(ancc)

        self.qc.h(ancc)
        self.qc.cx(3+7*pos, anc)
        self.qc.cx(ancc, anc)
        self.qc.cx(4+7*pos, anc)
        self.qc.cx(5+7*pos, anc)
        self.qc.cx(ancc, anc)
        self.qc.cx(6+7*pos, anc)

        self.qc.id(anc), self.qc.h(ancc), self.qc.id(ancc)
        self.qc.measure(anc, self.qecc[0]), self.qc.measure(ancc, flags[2])
        self.qc.reset(anc), self.qc.reset(ancc)
        self.qc.id(anc), self.qc.id(ancc)
        ##################################X-Stabilizers##############################################
        self.qc.h(anc)
        self.qc.cx(anc, 0+7*pos)
        self.qc.cx(anc, ancc)
        self.qc.cx(anc, 2+7*pos)
        self.qc.cx(anc, 4+7*pos)
        self.qc.cx(anc, ancc)
        self.qc.cx(anc, 6+7*pos)
        self.qc.h(anc)

        self.qc.id(anc), self.qc.id(ancc)
        self.qc.measure(anc, self.qecc[5]), self.qc.measure(ancc, flags[3])
        self.qc.reset(anc), self.qc.reset(ancc)
        self.qc.id(anc), self.qc.id(ancc)

        self.qc.h(anc)
        self.qc.cx(anc, 1+7*pos)
        self.qc.cx(anc, ancc)
        self.qc.cx(anc, 2+7*pos)
        self.qc.cx(anc, 5+7*pos)
        self.qc.cx(anc, ancc)
        self.qc.cx(anc, 6+7*pos)
        self.qc.h(anc)

        self.qc.id(anc), self.qc.id(ancc)
        self.qc.measure(anc, self.qecc[4]), self.qc.measure(ancc, flags[4])
        self.qc.reset(anc), self.qc.reset(ancc)
        self.qc.id(anc), self.qc.id(ancc)

        self.qc.h(anc)
        self.qc.cx(anc, 3+7*pos)
        self.qc.cx(anc, ancc)
        self.qc.cx(anc, 4+7*pos)
        self.qc.cx(anc, 5+7*pos)
        self.qc.cx(anc, ancc)
        self.qc.cx(anc, 6+7*pos)
        self.qc.h(anc)

        self.qc.id(anc), self.qc.id(ancc)
        self.qc.measure(anc, self.qecc[3]), self.qc.measure(ancc, flags[5])
        self.qc.reset(anc), self.qc.reset(ancc)
        ##################################Bitflip Error correction##############################################
        
        with self.qc.if_test((self.qecc[0],0)):             #qbit 0
            with self.qc.if_test((self.qecc[1],0)):
                with self.qc.if_test((self.qecc[2],1)):
                    self.qc.x(0+7*pos)

        with self.qc.if_test((self.qecc[0],0)):             #qbit 1
            with self.qc.if_test((self.qecc[1],1)):
                with self.qc.if_test((self.qecc[2],0)):
                    self.qc.x(1+7*pos)
        
        with self.qc.if_test((self.qecc[0],0)):             #qbit 2
            with self.qc.if_test((self.qecc[1],1)):
                with self.qc.if_test((self.qecc[2],1)):
                    self.qc.x(2+7*pos)
        
        with self.qc.if_test((self.qecc[0],1)):             #qbit 3
            with self.qc.if_test((self.qecc[1],0)):
                with self.qc.if_test((self.qecc[2],0)):
                    self.qc.x(3+7*pos)
        
        with self.qc.if_test((self.qecc[0],1)):             #qbit 4
            with self.qc.if_test((self.qecc[1],0)):
                with self.qc.if_test((self.qecc[2],1)):
                    self.qc.x(4+7*pos)
        
        with self.qc.if_test((self.qecc[0],1)):             #qbit 5
            with self.qc.if_test((self.qecc[1],1)):
                with self.qc.if_test((self.qecc[2],0)):
                    self.qc.x(5+7*pos)
        
        with self.qc.if_test((self.qecc[0],1)):             #qbit 6
            with self.qc.if_test((self.qecc[1],1)):
                with self.qc.if_test((self.qecc[2],1)):
                    self.qc.x(6+7*pos)

        ##################################Phaseflip Error correction##############################################
        
        with self.qc.if_test((self.qecc[3],0)):             #qbit 0
            with self.qc.if_test((self.qecc[4],0)):
                with self.qc.if_test((self.qecc[5],1)):
                    self.qc.z(0+7*pos)

        with self.qc.if_test((self.qecc[3],0)):             #qbit 1
            with self.qc.if_test((self.qecc[4],1)):
                with self.qc.if_test((self.qecc[5],0)):
                    self.qc.z(1+7*pos)
        
        with self.qc.if_test((self.qecc[3],0)):             #qbit 2
            with self.qc.if_test((self.qecc[4],1)):
                with self.qc.if_test((self.qecc[5],1)):
                    self.qc.z(2+7*pos)
        
        with self.qc.if_test((self.qecc[3],1)):             #qbit 3
            with self.qc.if_test((self.qecc[4],0)):
                with self.qc.if_test((self.qecc[5],0)):
                    self.qc.z(3+7*pos)
        
        with self.qc.if_test((self.qecc[3],1)):             #qbit 4
            with self.qc.if_test((self.qecc[4],0)):
                with self.qc.if_test((self.qecc[5],1)):
                    self.qc.z(4+7*pos)
        
        with self.qc.if_test((self.qecc[3],1)):             #qbit 5
            with self.qc.if_test((self.qecc[4],1)):
                with self.qc.if_test((self.qecc[5],0)):
                    self.qc.z(5+7*pos)
        
        with self.qc.if_test((self.qecc[3],1)):             #qbit 6
            with self.qc.if_test((self.qecc[4],1)):
                with self.qc.if_test((self.qecc[5],1)):
                    self.qc.z(6+7*pos)

################################# Neues FTQEC Procol based on arxiv:1708.02246, page 5, "Flag 1-FTEC Protocol" ##########################
    def flagsyndrome(self, pos: int):
        self.qec_counter += 1
        flags = ClassicalRegister(6)
        self.qc.add_register(flags)
        anc = self.qc.num_qubits - 1
        ancc = anc - 1
        self.qc.reset(anc), self.qc.reset(ancc)
        ##################################Z-Stabilizers##########################################
        self.qc.h(ancc)
        self.qc.cx(0+7*pos, anc)
        self.qc.cx(ancc,anc)
        self.qc.cx(2+7*pos, anc)
        self.qc.cx(4+7*pos, anc)
        self.qc.cx(ancc,anc)
        self.qc.cx(6+7*pos, anc)

        self.qc.id(anc), self.qc.h(ancc), self.qc.id(ancc)
        self.qc.measure(anc, self.qecc[2]), self.qc.measure(ancc, flags[0])
        self.qc.reset(anc), self.qc.reset(ancc)
        self.qc.id(anc), self.qc.id(ancc)

        self.qc.h(ancc)
        self.qc.cx(1+7*pos, anc)
        self.qc.cx(ancc, anc)
        self.qc.cx(2+7*pos, anc)
        self.qc.cx(5+7*pos, anc)
        self.qc.cx(ancc, anc)
        self.qc.cx(6+7*pos, anc)

        self.qc.id(anc), self.qc.h(ancc), self.qc.id(ancc)
        self.qc.measure(anc, self.qecc[1]), self.qc.measure(ancc, flags[1])
        self.qc.reset(anc), self.qc.reset(ancc)
        self.qc.id(anc), self.qc.id(ancc)

        self.qc.h(ancc)
        self.qc.cx(3+7*pos, anc)
        self.qc.cx(ancc, anc)
        self.qc.cx(4+7*pos, anc)
        self.qc.cx(5+7*pos, anc)
        self.qc.cx(ancc, anc)
        self.qc.cx(6+7*pos, anc)

        self.qc.id(anc), self.qc.h(ancc), self.qc.id(ancc)
        self.qc.measure(anc, self.qecc[0]), self.qc.measure(ancc, flags[2])
        self.qc.reset(anc), self.qc.reset(ancc)
        self.qc.id(anc), self.qc.id(ancc)
        ##################################X-Stabilizers##############################################
        self.qc.h(anc)
        self.qc.cx(anc, 0+7*pos)
        self.qc.cx(anc, ancc)
        self.qc.cx(anc, 2+7*pos)
        self.qc.cx(anc, 4+7*pos)
        self.qc.cx(anc, ancc)
        self.qc.cx(anc, 6+7*pos)
        self.qc.h(anc)

        self.qc.id(anc), self.qc.id(ancc)
        self.qc.measure(anc, self.qecc[5]), self.qc.measure(ancc, flags[3])
        self.qc.reset(anc), self.qc.reset(ancc)
        self.qc.id(anc), self.qc.id(ancc)

        self.qc.h(anc)
        self.qc.cx(anc, 1+7*pos)
        self.qc.cx(anc, ancc)
        self.qc.cx(anc, 2+7*pos)
        self.qc.cx(anc, 5+7*pos)
        self.qc.cx(anc, ancc)
        self.qc.cx(anc, 6+7*pos)
        self.qc.h(anc)

        self.qc.id(anc), self.qc.id(ancc)
        self.qc.measure(anc, self.qecc[4]), self.qc.measure(ancc, flags[4])
        self.qc.reset(anc), self.qc.reset(ancc)
        self.qc.id(anc), self.qc.id(ancc)

        self.qc.h(anc)
        self.qc.cx(anc, 3+7*pos)
        self.qc.cx(anc, ancc)
        self.qc.cx(anc, 4+7*pos)
        self.qc.cx(anc, 5+7*pos)
        self.qc.cx(anc, ancc)
        self.qc.cx(anc, 6+7*pos)
        self.qc.h(anc)

        self.qc.id(anc), self.qc.id(ancc)
        self.qc.measure(anc, self.qecc[3]), self.qc.measure(ancc, flags[5])
        self.qc.reset(anc), self.qc.reset(ancc)

        return flags

    def reset_qec(self, result):
        psi_full = result.data(0)["psi"]
        qr2 = QuantumRegister(7*(self.n+1)+2, "q")
        self.qc = QuantumCircuit(qr2, self.cbits)       #same register as in __init__, otherwise the
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
        self.qec_counter += 1
        anc = self.qc.num_qubits - 1
        self.qc.reset(anc)
        ##################################Z-Stabilizers##########################################
        for bit, stab in zip((2,1,0), STABS_7Q):
            for q in stab:
                self.qc.cx(q+7*pos, anc)
            self.qc.id(anc)
            self.qc.measure(anc, self.qecc[bit])
            self.qc.reset(anc)
            self.qc.id(anc)
        ##################################X-Stabilizers##############################################
        for bit, stab in zip((5,4,3), STABS_7Q):
            self.qc.h(anc)
            for q in stab:
                self.qc.cx(anc, q+7*pos)
            self.qc.h(anc)
            self.qc.id(anc)
            self.qc.measure(anc, self.qecc[bit])
            self.qc.reset(anc)
            self.qc.id(anc)

        self.qc.save_statevector(label="psi")
        sim = AerSimulator(method="statevector", noise_model = self.noise_model)
        result = sim.run(self.qc, shots=1).result()
        bitstring = list(result.get_counts().keys())[0].replace(" ","")
        syndrome = self._creg_bits(bitstring, self.qecc)
        self.reset_qec(result)

        return syndrome

    def correct_7q(self, syndrome: str, pos: int):
        #E_min(s). syndrome = the 6 bit qecc slice:  qecc[5]qecc[4]qecc[3] | qecc[2]qecc[1]qecc[0]
        for q in CORR_7Q[syndrome[3:6]]:                #Z-Stabilizers --> X errors
            self.qc.x(q+7*pos)
        for q in CORR_7Q[syndrome[0:3]]:                #X-Stabilizers --> Z errors
            self.qc.z(q+7*pos)

    def flagcorrect(self, flags: str, pos: int):
        #3.Case of the protocol: a circuit C(g_i) flagged, so the syndrome is measured again with the
        #non-flag circuits. If the outcome s matches an element E of the flag error set E(g_i) that
        #element is applied, otherwise E_min(s).
        flag_i = len(flags)-1 - flags.index("1")        #the flags slice is printed high bit first
        syndrome = self.syndrome(pos=pos)

        qubits = CORR_FLAG_7Q[flag_i].get(syndrome)
        if qubits is None:
            self.correct_7q(syndrome, pos=pos)
            return
        for q in qubits:
            if flag_i < 3:
                self.qc.z(q+7*pos)                      #Z-Stabilizer circuits flag on Z errors
            else:
                self.qc.x(q+7*pos)                      #X-Stabilizer circuits flag on X errors

    def flagFTec(self, pos: int):
        #based on https://arxiv.org/pdf/1708.02246 , protocol on page 5, "Flag 1-FTEC Protocol"  --> no need of postselection, but we need statevector simulation for classical decoding
        #The protocol repeats the flagged syndrome measurement until case 1, 2 or 3 applies. After the
        #second round one of the three always applies, so two flagged rounds are enough.

        flags_r1, syndrome_r1, bitstring, result = self._flag_round(pos=pos)

        if flags_r1.count("1") != 0:            #3.Case, falls direkt geflagged wird
            self.reset_qec(result)
            self.flagcorrect(flags_r1, pos=pos)
            return

        preselection = self._creg_bits(bitstring, self.cbits)
        if preselection is not None:            #cbits, one bit for each logical qubit
            #cbits[0] is NOT an initialization flag, t()/tdg() measure the state injection into it, so it is
            #1 in about half of the shots. The flags sit in cbits[1..n], and the slice is printed MSB first.
            if preselection[:-1].count("1") != 0:
                self.preselection_flag = True               
                self.reset_qec(result)                      #trotzdem zuruecksetzen, sonst bleibt das save_statevector("psi") im Circuit
                return 

        self.reset_qec(result)

        flags_r2, syndrome_r2, bitstring, result = self._flag_round(pos=pos)

        self.reset_qec(result)

        #3.Case: a circuit C(g_i) flagged --> non-flag syndrome measurement + flag error set E(g_i)
        if flags_r2.count("1") != 0:
            self.flagcorrect(flags_r2, pos=pos)
        #1.Case: no flags and the syndrome was repeated twice in a row --> E_min(s)
        elif syndrome_r1 == syndrome_r2:
            self.correct_7q(syndrome_r2, pos=pos)
        #2.Case: no flags but the two syndromes differ --> non-flag syndrome measurement + E_min(s)
        else:
            self.qec(pos=pos)

#######################################################################################################################################

    def qec(self, pos: int):
        self.qec_counter += 1
        anc = self.qc.num_qubits - 1
        self.qc.reset(anc)
        ##################################Z-Stabilizers##########################################
        self.qc.cx(0+7*pos, anc)
        self.qc.cx(2+7*pos, anc)
        self.qc.cx(4+7*pos, anc)
        self.qc.cx(6+7*pos, anc)

        self.qc.id(anc)
        self.qc.measure(anc, self.qecc[2])
        self.qc.reset(anc)
        self.qc.id(anc)

        self.qc.cx(1+7*pos, anc)
        self.qc.cx(2+7*pos, anc)
        self.qc.cx(5+7*pos, anc)
        self.qc.cx(6+7*pos, anc)

        self.qc.id(anc)
        self.qc.measure(anc, self.qecc[1])
        self.qc.reset(anc)
        self.qc.id(anc)

        self.qc.cx(3+7*pos, anc)
        self.qc.cx(4+7*pos, anc)
        self.qc.cx(5+7*pos, anc)
        self.qc.cx(6+7*pos, anc)

        self.qc.id(anc)
        self.qc.measure(anc, self.qecc[0])
        self.qc.reset(anc)
        self.qc.id(anc)
        ##################################X-Stabilizers##############################################
        self.qc.h(anc)
        self.qc.cx(anc, 0+7*pos)
        self.qc.cx(anc, 2+7*pos)
        self.qc.cx(anc, 4+7*pos)
        self.qc.cx(anc, 6+7*pos)
        self.qc.h(anc)

        self.qc.id(anc)
        self.qc.measure(anc, self.qecc[5])
        self.qc.reset(anc)
        self.qc.id(anc)

        self.qc.h(anc)
        self.qc.cx(anc, 1+7*pos)
        self.qc.cx(anc, 2+7*pos)
        self.qc.cx(anc, 5+7*pos)
        self.qc.cx(anc, 6+7*pos)
        self.qc.h(anc)

        self.qc.id(anc)
        self.qc.measure(anc, self.qecc[4])
        self.qc.reset(anc)
        self.qc.id(anc)

        self.qc.h(anc)
        self.qc.cx(anc, 3+7*pos)
        self.qc.cx(anc, 4+7*pos)
        self.qc.cx(anc, 5+7*pos)
        self.qc.cx(anc, 6+7*pos)
        self.qc.h(anc)

        self.qc.id(anc)
        self.qc.measure(anc, self.qecc[3])
        self.qc.reset(anc)
        ##################################Bitflip Error correction##############################################
        
        with self.qc.if_test((self.qecc[0],0)):             #qbit 0
            with self.qc.if_test((self.qecc[1],0)):
                with self.qc.if_test((self.qecc[2],1)):
                    self.qc.x(0+7*pos)

        with self.qc.if_test((self.qecc[0],0)):             #qbit 1
            with self.qc.if_test((self.qecc[1],1)):
                with self.qc.if_test((self.qecc[2],0)):
                    self.qc.x(1+7*pos)
        
        with self.qc.if_test((self.qecc[0],0)):             #qbit 2
            with self.qc.if_test((self.qecc[1],1)):
                with self.qc.if_test((self.qecc[2],1)):
                    self.qc.x(2+7*pos)
        
        with self.qc.if_test((self.qecc[0],1)):             #qbit 3
            with self.qc.if_test((self.qecc[1],0)):
                with self.qc.if_test((self.qecc[2],0)):
                    self.qc.x(3+7*pos)
        
        with self.qc.if_test((self.qecc[0],1)):             #qbit 4
            with self.qc.if_test((self.qecc[1],0)):
                with self.qc.if_test((self.qecc[2],1)):
                    self.qc.x(4+7*pos)
        
        with self.qc.if_test((self.qecc[0],1)):             #qbit 5
            with self.qc.if_test((self.qecc[1],1)):
                with self.qc.if_test((self.qecc[2],0)):
                    self.qc.x(5+7*pos)
        
        with self.qc.if_test((self.qecc[0],1)):             #qbit 6
            with self.qc.if_test((self.qecc[1],1)):
                with self.qc.if_test((self.qecc[2],1)):
                    self.qc.x(6+7*pos)

        ##################################Phaseflip Error correction##############################################
        
        with self.qc.if_test((self.qecc[3],0)):             #qbit 0
            with self.qc.if_test((self.qecc[4],0)):
                with self.qc.if_test((self.qecc[5],1)):
                    self.qc.z(0+7*pos)

        with self.qc.if_test((self.qecc[3],0)):             #qbit 1
            with self.qc.if_test((self.qecc[4],1)):
                with self.qc.if_test((self.qecc[5],0)):
                    self.qc.z(1+7*pos)
        
        with self.qc.if_test((self.qecc[3],0)):             #qbit 2
            with self.qc.if_test((self.qecc[4],1)):
                with self.qc.if_test((self.qecc[5],1)):
                    self.qc.z(2+7*pos)
        
        with self.qc.if_test((self.qecc[3],1)):             #qbit 3
            with self.qc.if_test((self.qecc[4],0)):
                with self.qc.if_test((self.qecc[5],0)):
                    self.qc.z(3+7*pos)
        
        with self.qc.if_test((self.qecc[3],1)):             #qbit 4
            with self.qc.if_test((self.qecc[4],0)):
                with self.qc.if_test((self.qecc[5],1)):
                    self.qc.z(4+7*pos)
        
        with self.qc.if_test((self.qecc[3],1)):             #qbit 5
            with self.qc.if_test((self.qecc[4],1)):
                with self.qc.if_test((self.qecc[5],0)):
                    self.qc.z(5+7*pos)
        
        with self.qc.if_test((self.qecc[3],1)):             #qbit 6
            with self.qc.if_test((self.qecc[4],1)):
                with self.qc.if_test((self.qecc[5],1)):
                    self.qc.z(6+7*pos)

    def qec_ideal(self, pos: int):
        self.qec_counter += 1
        anc = self.qc.num_qubits - 1
        self.qc.reset(anc)
        ##################################Z-Stabilizers##########################################
        self.qc.append(cx_ideal, [anc, 0+7*pos])
        self.qc.append(cx_ideal, [anc, 2+7*pos])
        self.qc.append(cx_ideal, [anc, 4+7*pos])
        self.qc.append(cx_ideal, [anc, 6+7*pos])

        self.qc.measure(anc, self.qecc[2])
        self.qc.reset(anc)

        self.qc.append(cx_ideal, [anc, 1+7*pos])
        self.qc.append(cx_ideal, [anc, 2+7*pos])
        self.qc.append(cx_ideal, [anc, 5+7*pos])
        self.qc.append(cx_ideal, [anc, 6+7*pos])

        self.qc.measure(anc, self.qecc[1])
        self.qc.reset(anc)

        self.qc.append(cx_ideal, [anc, 3+7*pos])
        self.qc.append(cx_ideal, [anc, 4+7*pos])
        self.qc.append(cx_ideal, [anc, 5+7*pos])
        self.qc.append(cx_ideal, [anc, 6+7*pos])

        self.qc.measure(anc, self.qecc[0])
        self.qc.reset(anc)
        ##################################X-Stabilizers##############################################
        self.qc.append(h_ideal, [anc])
        self.qc.append(cx_ideal, [0+7*pos, anc])
        self.qc.append(cx_ideal, [2+7*pos, anc])
        self.qc.append(cx_ideal, [4+7*pos, anc])
        self.qc.append(cx_ideal, [6+7*pos, anc])
        self.qc.append(h_ideal, [anc])

        self.qc.measure(anc, self.qecc[5])
        self.qc.reset(anc)

        self.qc.append(h_ideal, [anc])
        self.qc.append(cx_ideal, [1+7*pos, anc])
        self.qc.append(cx_ideal, [2+7*pos, anc])
        self.qc.append(cx_ideal, [5+7*pos, anc])
        self.qc.append(cx_ideal, [6+7*pos, anc])
        self.qc.append(h_ideal, [anc])

        self.qc.measure(anc, self.qecc[4])
        self.qc.reset(anc)

        self.qc.append(h_ideal, [anc])
        self.qc.append(cx_ideal, [3+7*pos, anc])
        self.qc.append(cx_ideal, [4+7*pos, anc])
        self.qc.append(cx_ideal, [5+7*pos, anc])
        self.qc.append(cx_ideal, [6+7*pos, anc])
        self.qc.append(h_ideal, [anc])

        self.qc.measure(anc, self.qecc[3])
        self.qc.reset(anc)
        ##################################Bitflip Error correction##############################################
        
        with self.qc.if_test((self.qecc[0],0)):             #qbit 0
            with self.qc.if_test((self.qecc[1],0)):
                with self.qc.if_test((self.qecc[2],1)):
                    self.qc.append(x_ideal, [0+7*pos])

        with self.qc.if_test((self.qecc[0],0)):             #qbit 1
            with self.qc.if_test((self.qecc[1],1)):
                with self.qc.if_test((self.qecc[2],0)):
                    self.qc.append(x_ideal, [1+7*pos])
        
        with self.qc.if_test((self.qecc[0],0)):             #qbit 2
            with self.qc.if_test((self.qecc[1],1)):
                with self.qc.if_test((self.qecc[2],1)):
                    self.qc.append(x_ideal, [2+7*pos])
        
        with self.qc.if_test((self.qecc[0],1)):             #qbit 3
            with self.qc.if_test((self.qecc[1],0)):
                with self.qc.if_test((self.qecc[2],0)):
                    self.qc.append(x_ideal, [3+7*pos])
        
        with self.qc.if_test((self.qecc[0],1)):             #qbit 4
            with self.qc.if_test((self.qecc[1],0)):
                with self.qc.if_test((self.qecc[2],1)):
                    self.qc.append(x_ideal, [4+7*pos])
        
        with self.qc.if_test((self.qecc[0],1)):             #qbit 5
            with self.qc.if_test((self.qecc[1],1)):
                with self.qc.if_test((self.qecc[2],0)):
                    self.qc.append(x_ideal, [5+7*pos])
        
        with self.qc.if_test((self.qecc[0],1)):             #qbit 6
            with self.qc.if_test((self.qecc[1],1)):
                with self.qc.if_test((self.qecc[2],1)):
                    self.qc.append(x_ideal, [6+7*pos])

        ##################################Phaseflip Error correction##############################################
        
        with self.qc.if_test((self.qecc[3],0)):             #qbit 0
            with self.qc.if_test((self.qecc[4],0)):
                with self.qc.if_test((self.qecc[5],1)):
                    self.qc.append(z_ideal, [0+7*pos])

        with self.qc.if_test((self.qecc[3],0)):             #qbit 1
            with self.qc.if_test((self.qecc[4],1)):
                with self.qc.if_test((self.qecc[5],0)):
                    self.qc.append(z_ideal, [1+7*pos])
        
        with self.qc.if_test((self.qecc[3],0)):             #qbit 2
            with self.qc.if_test((self.qecc[4],1)):
                with self.qc.if_test((self.qecc[5],1)):
                    self.qc.append(z_ideal, [2+7*pos])
        
        with self.qc.if_test((self.qecc[3],1)):             #qbit 3
            with self.qc.if_test((self.qecc[4],0)):
                with self.qc.if_test((self.qecc[5],0)):
                    self.qc.append(z_ideal, [3+7*pos])
        
        with self.qc.if_test((self.qecc[3],1)):             #qbit 4
            with self.qc.if_test((self.qecc[4],0)):
                with self.qc.if_test((self.qecc[5],1)):
                    self.qc.append(z_ideal, [4+7*pos])
        
        with self.qc.if_test((self.qecc[3],1)):             #qbit 5
            with self.qc.if_test((self.qecc[4],1)):
                with self.qc.if_test((self.qecc[5],0)):
                    self.qc.append(z_ideal, [5+7*pos])
        
        with self.qc.if_test((self.qecc[3],1)):             #qbit 6
            with self.qc.if_test((self.qecc[4],1)):
                with self.qc.if_test((self.qecc[5],1)):
                    self.qc.append(z_ideal, [6+7*pos])

    def readout(self, pos: int, shots: int, p = 0):
        p_error = pauli_error([["X",p/2],["I",1-p],["Z",p/2]])
        p_error_2 = pauli_error([["XI",p/4],["IX",p/4],["II",1-p],["ZI",p/4],["IZ",p/4]])

        noise_model = NoiseModel()
        noise_model.add_all_qubit_quantum_error(p_error, ['x', "z", 'h', "s", "sdg", "id", "t", "tdg"])  # Apply to single-qubit gates
        noise_model.add_all_qubit_quantum_error(p_error_2, ['cx'])  # Apply to 2-qubit gates

        read = ClassicalRegister(7)
        self.qc.add_register(read)

        for i in range(7):
            self.qc.id(i+7*pos)
            self.qc.measure(i+7*pos,read[6-i])

        sim = AerSimulator()
        job = sim.run(self.qc, shots=shots, noise_model=noise_model)
        result = job.result()
        counts = result.get_counts()

        #print(counts)

        bitstring = list(counts.keys())
        bitstring = [i.replace(" ","") for i in bitstring]


        hmm = list(counts.values())

        allcbits = len(bitstring[0])                
        pre, preselected = [i[allcbits-3:allcbits-1] for i in bitstring], 0                 #Flags during intialization
        bits = [i[:7] for i in bitstring]                                                   #Bits that make up the logical qubits
        postprocess = [i[7:-9] for i in bitstring]                                 #Flags during qec to make it fault tolerant, if at least one strikes, need to discard shot

        # print(bitstring)
        # print(bits)
        # print("FLAGS: ", postprocess)

        for i in range(len(pre)):
            if pre[i].count("1") != 0:
                bits[i] = "pre"
                # print("AHAAAA wieso flag???")

        test_0 = ["0000000","1010101","0110011","1100110","0001111","1011010","0111100","1101001"]
        test_1 = ["1111111","0101010","1001100","0011001","1110000","0100101","1000011","0010110"]

        # if self.classical_ec:
        #     bits = self.__classical_error_correction__(bits)

        for i in range(len(bits)):
            for j in test_0:
                if j == bits[i]:
                    bits[i] = 0
                    break
            if bits[i] != 0:
                for j in test_1:
                    if j == bits[i]:
                        bits[i] = 1
                        break
            if bits[i] != 1 and bits[i] != 0 and bits[i] != "pre":
                # print("Wrong bitstring: ", bits[i])
                if self.postselection:
                    bits[i] = "post"
                    # print("Postselected!")
                else:
                    if np.random.rand() < 0.5:
                        bits[i] = 0
                    else:
                        bits[i] = 1

        if not self.postselection:                  #wenn postselect aus, dann soll er alle shots, die durch preselection raus sind, coinflip machen
            for i in range(len(bits)):
                if bits[i] == "pre":
                    if np.random.rand() < 0.5:
                        bits[i] = 0
                    else:
                        bits[i] = 1

        for i in range(len(postprocess)):                           #Postprocess von Flags bei QEC!
            if postprocess[i].count("1") != 0:
                if bits[i] != "pre" and bits[i] != "post":
                    bits[i] = "post"

        #print(bits)
        ones = 0
        zeros = 0
        post = 0
        #magic = 0

        for i in range(len(bits)):
            if bits[i] == 0:
                zeros += hmm[i]
            if bits[i] == 1:
                ones += hmm[i]
            if bits[i] == "post":
                post += hmm[i]
            if bits[i] == "pre":
                preselected += hmm[i]

        if shots == 1 and self.preselection_flag:                  #only valid for shots = 1 and we use the 1-FTEC-protocol!!!
            ones, zeros, post = 0, 0, 0
            preselected = 1
        
        ones = (ones/shots)
        zeros = (zeros/shots)
        post = (post/shots)
        preselected = (preselected/shots)

        self.ones = ones
        self.zeros = zeros
        self.post = post
        self.preselected = preselected

