from HPC.Upload.functions import code_goto
from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister
from qiskit.visualization import plot_histogram
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
#import bitstring
from qiskit_aer import AerSimulator
from qiskit.transpiler.passes.synthesis import SolovayKitaev
from qiskit.synthesis import generate_basic_approximations
from qiskit.quantum_info import Operator

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

def convert(bin: str):                  #konvertiert den bitstring in decimal, e.g. 0110 = 0.375
    k = list(bin)
    a = [int(i) for i in k]
    n = 0
    for i in range(len(a)):
        if a[i] == 1:
            n += 1/2**(i+1)
    return n

class Steane7q:
    def __init__(self, n: int, magic = 1):
        self.n = n

        self.zeros = 0
        self.ones = 0
        self.preselected = 0
        self.post = 0

        self.postselection = True

        qr = QuantumRegister(7*(n+magic)+2,"q")
        cbits = ClassicalRegister(3, "c")
        
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
        self.qc.s(0+7*pos), self.qc.s(1+7*pos), self.qc.s(3+7*pos), self.qc.s(6+7*pos)
        self.qc.sdg(2+7*pos), self.qc.sdg(4+7*pos), self.qc.sdg(5+7*pos)

    def cz(self, control: int, target: int):
        self.h(pos=control)
        self.cnot(control=target, target=control)
        self.h(pos=control)

    def cnot(self, control: int, target: int):
        for i in range(7):
            self.qc.cx(i+7*control, i+7*target)

    def sdg(self, pos: int):
        self.qc.sdg(0+7*pos), self.qc.sdg(1+7*pos), self.qc.sdg(3+7*pos), self.qc.sdg(6+7*pos)
        self.qc.s(2+7*pos), self.qc.s(4+7*pos), self.qc.s(5+7*pos)

    def t(self, pos: int):
        self.h(pos=pos)
        self.sdg(pos=pos)
        self.h(pos=pos)
        state_inj = ClassicalRegister(1)
        self.qc.add_register(state_inj)

        anc = self.qc.num_qubits - 1
        ancc = anc - 1

        for i in range(7):
            self.qc.reset(i+7*2)

        self.qc.append(h_ideal,[0+7*2])
        self.qc.append(h_ideal,[1+7*2])
        self.qc.ry(np.pi/4,2+7*2)
        self.qc.append(h_ideal,[3+7*2])

        self.qc.append(cx_ideal, [4+7*2, 2+7*2])
        self.qc.append(cx_ideal, [6+7*2, 0+7*2])
        
        self.qc.append(cx_ideal, [5+7*2, 3+7*2])

        self.qc.append(cx_ideal, [5+7*2, 2+7*2])

        self.qc.append(cx_ideal, [4+7*2, 0+7*2])
        self.qc.append(cx_ideal, [6+7*2, 1+7*2])

        self.qc.append(cx_ideal, [2+7*2, 0+7*2])

        self.qc.append(cx_ideal, [5+7*2, 1+7*2])

        self.qc.append(cx_ideal, [2+7*2, 1+7*2])
        self.qc.append(cx_ideal, [4+7*2, 3+7*2])
        self.qc.append(cx_ideal, [6+7*2, 3+7*2])
        #################################Controlled Hadamards##########################################
        self.qc.reset(ancc)
        self.qc.append(h_ideal,[ancc])
        for i in range(7):
            self.qc.ry(-np.pi/4,6-i+2*7)
            self.qc.cz(ancc,6-i+2*7)
            self.qc.ry(np.pi/4,6-i+2*7)
        self.qc.append(h_ideal,[ancc])
        self.qc.measure(ancc, state_inj[0])
        ########################Controlled-Y Gate####################################################
        self.sdg(pos=pos)
        for i in range(7):
            self.qc.cx(i+7*2,i+7*pos)
        self.s(pos=pos)
        self.qc.reset(anc-1)
        #############################Measure logical state of the magic state for state injection#############################
        self.sdg(pos=2)
        self.h(pos=2)
        for i in range(7):
            self.qc.cx(i+2*7, anc-1)
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

    def tdg(self, pos: int):
        self.h(pos=pos)
        self.sdg(pos=pos)
        self.h(pos=pos)
        state_inj = ClassicalRegister(1)
        self.qc.add_register(state_inj)

        anc = self.qc.num_qubits - 1
        ancc = anc - 1

        for i in range(7):
            self.qc.reset(i+7*2)

        self.qc.append(h_ideal,[0+7*2])
        self.qc.append(h_ideal,[1+7*2])
        self.qc.ry(np.pi/4,2+7*2)
        self.qc.append(h_ideal,[3+7*2])

        self.qc.append(cx_ideal, [4+7*2, 2+7*2])
        self.qc.append(cx_ideal, [6+7*2, 0+7*2])
        
        self.qc.append(cx_ideal, [5+7*2, 3+7*2])

        self.qc.append(cx_ideal, [5+7*2, 2+7*2])

        self.qc.append(cx_ideal, [4+7*2, 0+7*2])
        self.qc.append(cx_ideal, [6+7*2, 1+7*2])

        self.qc.append(cx_ideal, [2+7*2, 0+7*2])

        self.qc.append(cx_ideal, [5+7*2, 1+7*2])

        self.qc.append(cx_ideal, [2+7*2, 1+7*2])
        self.qc.append(cx_ideal, [4+7*2, 3+7*2])
        self.qc.append(cx_ideal, [6+7*2, 3+7*2])
        #################################Controlled Hadamards##########################################
        self.qc.reset(ancc)
        self.qc.append(h_ideal,[ancc])
        for i in range(7):
            self.qc.ry(-np.pi/4,6-i+2*7)
            self.qc.cz(ancc,6-i+2*7)
            self.qc.ry(np.pi/4,6-i+2*7)
        self.qc.append(h_ideal,[ancc])
        self.qc.measure(ancc, state_inj[0])
        ########################Controlled-Y Gate####################################################
        self.sdg(pos=pos)
        for i in range(7):
            self.qc.cx(i+7*2,i+7*pos)
        self.s(pos=pos)
        self.qc.reset(anc-1)
        #############################Measure logical state of the magic state for state injection#############################
        self.sdg(pos=2)
        self.h(pos=2)
        for i in range(7):
            self.qc.cx(i+2*7, anc-1)
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

    def cu(self, gate: list, adjgate: list):
        self.u2(0, gate=gate)
        self.u2(1, gate=gate)
        self.qc.cnot(control=0, target=1)
        self.u2(1, gate=adjgate)
        self.qc.cnot(control=0, target=1)

    def qec_ft(self, pos: int, ft: bool):
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
            with self.c.if_test((self.qecc[1],1)):
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
        return True

    def err_mitigation(self, setting: bool):
        self.postselection = setting
    
    def avg15_coin(self, iter: int, noise: float, err = False, k = 1, path = ""):       #each iteration own circuit
        n = 15
        angle = np.linspace(0,1,n+2)
        angle = np.delete(angle, [n+1])
        angle = np.delete(angle, [0])

        a, b = [], []
        with open("{}unitary{}.txt".format(path, n), "r") as file:
            for line in file:
                a.append(list(map(str, line.strip().split(","))))
        with open("{}adjunitary{}.txt".format(path, n), "r") as file:
            for line in file:
                b.append(list(map(str, line.strip().split(","))))
        
        y = 0
        bruh1 = []
        for m in range(k):
            for o in range(n):
                bitstring = ""
                rots = []
                for t in range(iter):
                    rots = [k*0.5 for k in rots]

                    self.x(pos=1)
                    self.h(pos=0)
                    #############################
                    for j in range(2**(iter-t-1)):
                        self.cu(a[o], b[o])
                    ###############################
                    for l in rots:
                        if l == 0.25:
                            self.sdg(pos=0)
                        if l == 0.125:
                            self.tdg(pos=0, err=err)
                    self.h(pos=0)
                    if err:
                        self.qec_ft(pos=0)
                    self.readout(pos=0, shots=1, noise=noise)
            
                    if self.zeros == 1:
                        bitstring += "0"
                    elif self.ones == 1:
                        bitstring += "1"
                        rots.append(0.5)
                    else:
                        if np.random.rand() < 0.5:
                            bitstring += "0"
                        else:
                            bitstring += "1"
                            rots.append(0.5)
                bitstring = bitstring[::-1]
                hmm = convert(bitstring)
                diff = np.abs(hmm-angle[o])
                y += diff
                bruh1.append(diff)
        y = y/(n*k)
        arg = 0
        for i in range(len(bruh1)):
            arg += (y-bruh1[i])**2
        sigma = ((1/(k*n))*arg)**0.5
        sigma = sigma/((k*n)**0.5)

        return y, sigma

    def readout(self, pos: int, shots: int, noise = 0):
        p = noise
        p_error = pauli_error([["X",p/2],["I",1-p],["Z",p/2]])
        p_error_2 = pauli_error([["XI",p/4],["IX",p/4],["II",1-p],["ZI",p/4],["IZ",p/4]])

        noise_model = NoiseModel()
        noise_model.add_all_qubit_quantum_error(p_error, ['x', "z", 'h', "s", "sdg", "id"])  # Apply to single-qubit gates
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

        #print(counts)

        bitstring = list(counts.keys())
        bitstring = [i.replace(" ","") for i in bitstring]


        hmm = list(counts.values())

        allcbits = len(bitstring[0])                
        pre, preselected = [i[allcbits-3:allcbits-1] for i in bitstring], 0                 #Flags during intialization
        bits = [i[:7] for i in bitstring]                                                   #Bits that make up the logical qubits
        postprocess = [i[7:allcbits-10] for i in bitstring]                                 #Flags during qec to make it fault tolerant, if at least one strikes, need to discard shot

        #print(bits)
        #print(postprocess)

        for i in range(len(pre)):
            if pre[i].count("1") != 0:
                bits[i] = "pre"
                preselected += hmm[i]

        test_0 = ["0000000","1010101","0110011","1100110","0001111","1011010","0111100","1101001"]
        test_1 = ["1111111","0101010","1001100","0011001","1110000","0100101","1000011","0010110"]

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
                if self.postselection:
                    bits[i] = "post"
                else:
                    if np.random.rand() < 0.5:
                        bits[i] = 0
                    else:
                        bits[i] = 1

        for i in range(len(postprocess)):
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
            # if bits[i] == "magic":
            #     magic += hmm[i]
        
        ones = (ones/shots)
        zeros = (zeros/shots)
        post = (post/shots)
        preselected = (preselected/shots)

        self.ones = ones
        self.zeros = zeros
        self.post = post
        self.preselected = preselected
        #magic = (magic/shots)

        # print("0: ", zeros*100, "%")
        # print("1: ", ones*100, "%")
        # print("Preselection discarded: ", (preselected/shots)*100, "%")
        # print("Postselection discarded: ", (post/shots)*100, "%")
        #return zeros, ones, preselected, post

class RotSurf9q:
    def __init__(self, n: int):
        self.n = n

        self.zeros = 0
        self.ones = 0
        self.preselected = 0
        self.post = 0

        self.postselection = True

        self.hadamards = [0,0]

        qr = QuantumRegister(9*n+1, "q")
        cbit = ClassicalRegister(9,"c")
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
    
    def sdg(self, pos: int):
        anc = self.qc.num_qubits - 1
        self.qc.reset(anc)

        self.qc.h(anc)
        #self.qc.append(h_ideal,[anc])
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

        self.qc.measure(anc, 4)
        # if z_stab:
        #     z_qec_ideal(qc, had=had, pos=pos)

        if self.hadamards[pos]%2 == 0:
            with self.qc.if_test((4,1)):
                self.qc.reset(anc)
                self.qc.append(h_ideal,[anc])
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
            with self.qc.if_test((4,1)):
                self.qc.reset(anc)
                self.qc.append(h_ideal,[anc])
                self.qc.s(anc)
                self.qc.cx(1+9*pos, anc)
                self.qc.cx(4+9*pos, anc)
                self.qc.cx(7+9*pos, anc) 
                self.qc.measure(anc, 0)
                with self.qc.if_test((0,1)):
                    self.qc.z(1+9*pos)
                    self.qc.z(4+9*pos)
                    self.qc.z(7+9*pos)

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

        self.qc.measure(anc, 4)
        # if z_stab:
        #     z_qec_ideal(qc, had=had, pos=pos)

        if self.hadamards[pos]%2 == 0:
            with self.qc.if_test((4,1)):
                self.qc.reset(anc)
                self.qc.append(h_ideal,[anc])
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
            with self.qc.if_test((4,1)):
                self.qc.reset(anc)
                self.qc.append(h_ideal,[anc])
                self.qc.sdg(anc)
                self.qc.cx(1+9*pos, anc)
                self.qc.cx(4+9*pos, anc)
                self.qc.cx(7+9*pos, anc) 
                self.qc.measure(anc, 0)
                with self.qc.if_test((0,1)):
                    self.qc.z(1+9*pos)
                    self.qc.z(4+9*pos)
                    self.qc.z(7+9*pos)

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

    def CU_L(self, Ugates: list, adjUgates: list, err = False):
        self.u2(0, Ugates)
        if err:
            self.qec(pos=0)
        self.u2(1, Ugates)
        if err:
            self.qec(pos=1)
        self.cnot(control=0, target=1)
        self.u2(1, adjUgates)
        if err:
            self.qec(pos=0)
        self.cnot(control=0, target=1)

    def avg15coin(self, iter: int, n:int, noise: float, err = False, k = 1, path = ""):       #each iteration own circuit
        angle = np.linspace(0,1,n+2)
        angle = np.delete(angle, [n+1])
        angle = np.delete(angle, [0])

        a, b = [], []
        with open("{}unitary{}.txt".format(path, n), "r") as file:
            for line in file:
                a.append(list(map(str, line.strip().split(","))))
        with open("{}adjunitary{}.txt".format(path, n), "r") as file:
            for line in file:
                b.append(list(map(str, line.strip().split(","))))
        
        y = 0
        bruh1 = []
        for m in range(k):
            for o in range(n):
                bitstring = ""
                rots = []
                for t in range(iter):
                    rots = [k*0.5 for k in rots]

                    self.x(pos=1)
                    self.h(pos=0)
                    #############################
                    for j in range(2**(iter-t-1)):
                        self.cu(a[o], b[o], err=err)
                    ###############################
                    for l in rots:
                        if l == 0.25:
                            self.sdg(pos=0)
                        if l == 0.125:
                            self.tdg(pos=0)
                    self.h(pos=0)
                    if err:
                        self.qec(pos = 0)
                    self.readout(pos=0, shots=1, noise=noise)
            
                    if self.zeros == 1:
                        bitstring += "0"
                    elif self.ones == 1:
                        bitstring += "1"
                        rots.append(0.5)
                    else:
                        if np.random.rand() < 0.5:
                            bitstring += "0"
                        else:
                            bitstring += "1"
                            rots.append(0.5)
                bitstring = bitstring[::-1]
                hmm = convert(bitstring)
                diff = np.abs(hmm-angle[o])
                y += diff
                bruh1.append(diff)
        y = y/(n*k)
        arg = 0
        for i in range(len(bruh1)):
            arg += (y-bruh1[i])**2
        sigma = ((1/(k*n))*arg)**0.5
        sigma = sigma/((k*n)**0.5)

        return y, sigma

    def qec(self, pos: int):
        anc = self.qc.num_qubits - 1
        if self.hadamards[pos]%2==1:
            #X3 X6 Stabilizer:
            self.qc.reset(anc)
            self.qc.h(anc)
            self.qc.cx(anc, 3+9*pos)
            self.qc.cx(anc, 6+9*pos)
            self.qc.h(anc)
            self.qc.id(anc)
            self.qc.measure(anc,0)

            #X0 X1 X3 X4 Stabilizer:
            self.qc.reset(anc)
            self.qc.h(anc)
            self.qc.cx(anc, 0+9*pos)
            self.qc.cx(anc, 1+9*pos)
            self.qc.cx(anc, 3+9*pos)
            self.qc.cx(anc, 4+9*pos)
            self.qc.h(anc)
            self.qc.id(anc)
            self.qc.measure(anc,1)

            #X4 X5 X7 X8 Stabilizer:
            self.qc.reset(anc)
            self.qc.h(anc)
            self.qc.cx(anc, 4+9*pos)
            self.qc.cx(anc, 5+9*pos)
            self.qc.cx(anc, 7+9*pos)
            self.qc.cx(anc, 8+9*pos)
            self.qc.h(anc)
            self.qc.id(anc)
            self.qc.measure(anc,2)

            #X2 X5 Stabilizer:
            self.qc.reset(anc)
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
            self.qc.cx(0+9*pos, anc)
            self.qc.cx(1+9*pos, anc)
            self.qc.id(anc)
            self.qc.measure(anc,0)

            #Z1 Z2 Z4 Z5 Stabilizer:
            self.qc.reset(anc)
            self.qc.cx(1+9*pos, anc)
            self.qc.cx(2+9*pos, anc)
            self.qc.cx(4+9*pos, anc)
            self.qc.cx(5+9*pos, anc)
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
            #X0 X1 Stabilizer:
            self.qc.reset(anc)
            self.qc.h(anc)
            self.qc.cx(anc, 0+9*pos)
            self.qc.cx(anc, 1+9*pos)
            self.qc.h(anc)
            self.qc.id(anc)
            self.qc.measure(anc,0)
            
            #X1 X2 X4 X5 Stabilizer:
            self.qc.reset(anc)
            self.qc.h(anc)
            self.qc.cx(anc, 1+9*pos)
            self.qc.cx(anc, 2+9*pos)
            self.qc.cx(anc, 4+9*pos)
            self.qc.cx(anc, 5+9*pos)
            self.qc.h(anc)
            self.qc.id(anc)
            self.qc.measure(anc,1)

            #X3 X4 X6 X7 Stabilizer:
            self.qc.reset(anc)
            self.qc.h(anc)
            self.qc.cx(anc, 3+9*pos)
            self.qc.cx(anc, 4+9*pos)
            self.qc.cx(anc, 6+9*pos)
            self.qc.cx(anc, 7+9*pos)
            self.qc.h(anc)
            self.qc.id(anc)
            self.qc.measure(anc,2)

            #X7 X8 Stabilizer:
            self.qc.reset(anc)
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
            self.qc.cx(3+9*pos, anc)
            self.qc.cx(6+9*pos, anc)
            self.qc.id(anc)
            self.qc.measure(anc,0)

            #Z0 Z1 Z3 Z4 Stabilizer:
            self.qc.reset(anc)
            self.qc.cx(0+9*pos, anc)
            self.qc.cx(1+9*pos, anc)
            self.qc.cx(3+9*pos, anc)
            self.qc.cx(4+9*pos, anc)
            self.qc.id(anc)
            self.qc.measure(anc,1)
        
            #Z4 Z5 Z7 Z8 Stabilizer:
            self.qc.reset(anc)
            self.qc.cx(4+9*pos, anc)
            self.qc.cx(5+9*pos, anc)
            self.qc.cx(7+9*pos, anc)
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

    def err_mitigation(self, setting: bool):
        self.postselection = setting

    def readout(self, pos: int, shots: int, noise = 0):
        code0 = ['000110101', '110110110', '110110101', '110000000', '000110110', '101101101', '011101101', '011011000', '011011011', '110000011', '000000000', '011101110', '101011011', '101101110', '000000011', '101011000']
        code1 = ['010100111', '010010001', '111111111', '001001010', '111001010', '001111111', '100010010', '111111100', '100100100', '100010001', '001001001', '010010010', '100100111', '111001001', '001111100', '010100100']

        # for i in range(9):
        #     qc.id(i+9*pos)
        if self.hadamards[pos]%2 == 0:
            for i in range(9):
                self.qc.measure(i+9*pos, 8-i)
        else:
            self.qc.measure(0+9*pos, 8-6)
            self.qc.measure(1+9*pos, 8-3)
            self.qc.measure(2+9*pos, 8-0)
            self.qc.measure(3+9*pos, 8-7)
            self.qc.measure(4+9*pos, 8-4)
            self.qc.measure(5+9*pos, 8-1)
            self.qc.measure(6+9*pos, 8-8)
            self.qc.measure(7+9*pos, 8-5)
            self.qc.measure(8+9*pos, 8-2)

        p = noise
        p_error = pauli_error([["X",p/2],["I",1-p],["Z",p/2]])
        p_error_2 = pauli_error([["XI",p/4],["IX",p/4],["II",1-p],["ZI",p/4],["IZ",p/4]])

        noise_model = NoiseModel()
        noise_model.add_all_qubit_quantum_error(p_error, ['x', "z", 'h', "id", "s", "sdg", "t", "tdg"])  # Apply to single-qubit gates
        noise_model.add_all_qubit_quantum_error(p_error_2, ['cx'])  # Apply to 2-qubit gates

        sim = AerSimulator()
        job = sim.run(self.qc, noise_model = noise_model, shots=shots)
        result = job.result()
        counts = result.get_counts()

        #print(counts)

        bits = list(counts.keys())
        #print("Anzahl der versch. Bitstring: ", len(bits))
        hmm = list(counts.values())

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
                    bits[i] = 2
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

        ones = 0
        zeros = 0
        err = 0

        for i in range(len(bits)):
            if bits[i] == 0:
                zeros += hmm[i]
            if bits[i] == 1:
                ones += hmm[i]
            if bits[i] == 2:
                err += hmm[i]
        
        ones = (ones/shots)
        zeros = (zeros/shots)
        err = (err/shots)

        self.ones = ones
        self.zeros = zeros
        self.post = err
        #return zeros, ones, err
