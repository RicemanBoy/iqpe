from pdb import pm

from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister, qasm3, qasm2, qpy
from qiskit.visualization import plot_histogram
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
#import bitstring
from qiskit_aer import AerSimulator
from qiskit.transpiler.passes.synthesis import SolovayKitaev
from qiskit.synthesis import generate_basic_approximations
from qiskit.quantum_info import Operator

from qiskit import transpile
from qiskit.transpiler import PassManager
from qiskit.transpiler.passes import BasisTranslator
from qiskit.circuit.equivalence_library import SessionEquivalenceLibrary

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

def gates(qc:QuantumCircuit):
    hmm = dict(qc.count_ops())
    hmm["reset"] = 0
    hmm["measure"] = 0
    hmm["if_else"] = 0
    #hmm["h"] = 0
    return print("Amount of gates in this circuit: ", sum(hmm.values()))

def convert(bin: str):                  #konvertiert den bitstring in decimal, e.g. 0110 = 0.375
    k = list(bin)
    a = [int(i) for i in k]
    n = 0
    for i in range(len(a)):
        if a[i] == 1:
            n += 1/2**(i+1)
    return n

def inv_covert(num: int, precision: int):
    max = sum(2**i for i in range(precision))
    assert num <= max, "Precision is not sufficient for this integer!"
    bitstring = ""
    for i in range(precision):
        if num >= 2**(precision-i-1):
            num = num - 2**(precision-i-1)
            bitstring += "1"
        else:
            bitstring += "0"
    
    return bitstring

def closest_bitstring(num: float, depth: int):
    bitstring = ""
    
    for i in range(depth):
        if num >= 2**(-i-1):
            bitstring += "1"
            num = num - 2**(-i-1)
        else:
            bitstring += "0"
    
    return bitstring

def avg15_coin(code: str, iter: int, noise: float, qec = False, k = 1, path = ""):       #each iteration own circuit
    assert code in ["steane", "rotsurf9", "rotsurf16"], "Code not supported! Choose between 'steane', 'rotsurf9' and 'rotsurf16'."
    
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
                if code == "steane":
                    self = Steane7q(2)
                elif code == "rotsurf9":
                    self = RotSurf9q(2)
                elif code == "rotsurf16":
                    self = RotSurf16q(2)
                self.err = qec
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
                        self.tdg(pos=0)
                self.h(pos=0)
                if self.err:
                    self.qec(pos = 0)
                # print("Unoptimized: ")
                # gates(self.qc)
                self.qc = transpile(self.qc, optimization_level=1)
                # print("Optimized: ")
                # gates(self.qc)
                # print(self.qec_counter)
                self.readout(pos=0, shots=1, p=noise)
        
                if self.zeros == 1:
                    bitstring += "0"
                elif self.ones == 1:
                    bitstring += "1"
                    rots.append(0.5)
                else:
                    print("Error")
                    if np.random.rand() < 0.5:
                        bitstring += "0"
                    else:
                        bitstring += "1"
                        rots.append(0.5)
                del self
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

def avg14_coin_success(code: str, iter: int, noise: float, qec = False, k = 1, path = ""):       #each iteration own circuit
    n = 14
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
                if code == "steane":
                    self = Steane7q(2)
                else:
                    self = RotSurf9q(2)
                # else:
                #     self = RotSurf16q(2)
                self.err = qec
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
                        self.tdg(pos=0)
                self.h(pos=0)
                if self.err:
                    self.qec(pos = 0)
                # print("Unoptimized: ")
                # gates(self.qc)
                #new_qc = transpile(self.qc, optimization_level=1)
                self.qc = transpile(self.qc, optimization_level=1)
                # print("Optimized: ")
                # gates(new_qc)

                self.readout(pos=0, shots=1, p=noise)

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
                del self
            bitstring = bitstring[::-1]
            angle_bit = closest_bitstring(angle[o], iter)
            print("bitstring measured: ", bitstring)
            print("Expected bitstring: ", angle_bit)
            if bitstring == angle_bit:
                y += 1

    y = y/(k*n)
    return y

def avg_15_coin_circ(thisangle: int, iter: int, rots: list, qec = False, path = "", name=""):       #each iteration own circuit
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
    
    bitstring = ""
    self = RotSurf16q(2)
    self.err = qec
    rots = [k*0.5 for k in rots]

    self.x(pos=1)
    self.h(pos=0)
    #############################
    for j in range(2**(2-iter)):
        self.cu(a[thisangle], b[thisangle])
    ###############################
    for l in rots:
        if l == 0.25:
            self.sdg(pos=0)
        if l == 0.125:
            self.tdg(pos=0)
    self.h(pos=0)
    # if self.err:
    #     self.qec(pos = 0)
    # for i in range(16):
    #         self.qc.measure(i, i)

    # print(dict(self.qc.count_ops()))
    # gates(self.qc)

    new_qc = transpile(self.qc, optimization_level=1)
    # text_circuit = new_qc.draw(output="text")
    # with open('Circuits for Timo/{}.txt'.format(name), 'w') as f:
    #     f.write(str(text_circuit))


    # qasm_str = qasm3.dumps(new_qc)
    # with open("Circuits for Timo/{}.qasm".format(name), "w") as f:
    #     f.write(qasm_str)


    with open("Circuits for Timo/{}.qpy".format(name), "wb") as file:
        qpy.dump(new_qc, file)

    # print(dict(new_qc.count_ops()))
    # gates(new_qc)

def avg15(code: str, iter: int, noise: float, qec = False, k = 1, path=""):       #each iteration own circuit
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
                while True:
                    if code == "steane":
                        self = Steane7q(2, magic=0)
                    elif code == "rotsurf9":
                        self = RotSurf9q(2)
                    else:
                        self = RotSurf16q(2)
                    self.err = qec
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
                            self.nFTtdg(pos=0)
                    self.h(pos=0)
                    if self.err:
                        self.qec(pos = 0)
                    self.qc = transpile(self.qc, optimization_level=1)
                    self.readout(pos=0, shots=1, p=noise)

                    if self.zeros == 1:
                        bitstring += "0"
                        break
                    if self.ones == 1:
                        bitstring += "1"
                        rots.append(0.5)
                        break
            bitstring = bitstring[::-1]
            print(bitstring)
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

# def decoding(counts: dict):


class Steane7q:
    def __init__(self, n: int, magic = 1):
        self.n = n

        self.zeros = 0
        self.ones = 0
        self.preselected = 0
        self.post = 0
        self.err = False
        self.qec_counter = 0

        self.classical_ec = False
        self.postselection = False

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

    def t_ghz(self, pos: int):
        anc = self.qc.num_qubits - 1
        self.qc.reset(anc)
        self.qc.reset(anc-1)
        self.qc.reset(anc-2)
        self.qc.reset(anc-3)

        self.qc.h(anc)
        self.qc.t(anc)
        self.qc.cx(anc, anc-1)
        self.qc.cx(anc, anc-2)
        
        self.qc.cx(0+7*pos, anc)
        self.qc.cx(1+7*pos, anc-1)
        self.qc.cx(2+7*pos, anc-2)          

        self.qc.cx(anc, anc-3)
        self.qc.cx(anc-1, anc-3)
        self.qc.cx(anc-2, anc-3)

        flags = ClassicalRegister(1)
        self.qc.add_register(flags)
        self.qc.measure(anc-3, flags[0])

        with self.qc.if_test((flags[0],1)):
            self.qc.s(0+7*pos)
            self.qc.s(1+7*pos)
            self.qc.s(3+7*pos)
            self.qc.s(6+7*pos)
            self.qc.sdg(2+7*pos)
            self.qc.sdg(4+7*pos)
            self.qc.sdg(5+7*pos)
        
        self.qc.cx(1+7*pos, anc-1)

        self.qc.cx(anc-1, anc)
        self.qc.cx(anc-1, anc-2)

        self.qc.cx(anc, anc-1)
        self.qc.cx(anc-2, anc-1)

        self.qc.cx(0+7*pos, anc)
        self.qc.cx(1+7*pos, anc-1)
        self.qc.cx(2+7*pos, anc-2)

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
        # state_inj = ClassicalRegister(1)
        # self.qc.add_register(state_inj)

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
        # self.qc.measure(ancc, state_inj[0])
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
        if self.err:
            self.qec(pos=pos)

    def nFTt(self,pos: int):
        anc = self.qc.num_qubits - 1
        self.qc.reset(anc)

        self.qc.h(anc)
        self.qc.t(anc)

        self.qc.cx(0+7*pos, anc)
        self.qc.cx(1+7*pos, anc)
        self.qc.cx(2+7*pos, anc)        

        self.qc.measure(anc, 0)

        with self.qc.if_test((0,1)):
            self.qc.reset(anc)
            self.qc.h(anc)
            self.qc.s(anc)
            self.qc.cx(0+7*pos, anc)
            self.qc.cx(1+7*pos, anc)
            self.qc.cx(2+7*pos, anc) 
            self.qc.measure(anc, 0)
            with self.qc.if_test((0,1)):
                self.qc.z(0+7*pos)
                self.qc.z(1+7*pos)
                self.qc.z(2+7*pos)

    def rz_cheat(self, angle: float, pos: int):
        self.qc.cx(0+7*pos, 2+7*pos)
        self.qc.cx(1+7*pos, 2+7*pos)
        self.qc.rz(angle, 2+7*pos)
        self.qc.cx(0+7*pos, 2+7*pos)
        self.qc.cx(1+7*pos, 2+7*pos)
    
    def rx_cheat(self, angle: float, pos: int):
        self.qc.cx(2+7*pos, 1+7*pos)
        self.qc.cx(2+7*pos, 0+7*pos)
        self.qc.rx(angle, 2+7*pos)
        self.qc.cx(2+7*pos, 0+7*pos)
        self.qc.cx(2+7*pos, 1+7*pos)

    def ry_cheat(self, angle: float, pos: int):
        self.qc.h(0+7*pos), self.qc.h(1+7*pos), self.qc.h(2+7*pos)
        self.qc.s(0+7*pos), self.qc.s(1+7*pos), self.qc.s(2+7*pos)
        self.qc.h(0+7*pos), self.qc.h(1+7*pos), self.qc.h(2+7*pos)
        self.qc.cx(0+7*pos, 2+7*pos)
        self.qc.cx(1+7*pos, 2+7*pos)
        self.qc.rz(angle, 2+7*pos)
        self.qc.cx(0+7*pos, 2+7*pos)
        self.qc.cx(1+7*pos, 2+7*pos)
        self.qc.h(0+7*pos), self.qc.h(1+7*pos), self.qc.h(2+7*pos)
        self.qc.sdg(0+7*pos), self.qc.sdg(1+7*pos), self.qc.sdg(2+7*pos)
        self.qc.h(0+7*pos), self.qc.h(1+7*pos), self.qc.h(2+7*pos)

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
        self.h(pos=pos)
        self.sdg(pos=pos)
        self.h(pos=pos)
        # state_inj = ClassicalRegister(1)
        # self.qc.add_register(state_inj)

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
        # self.qc.measure(ancc, state_inj[0])
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
            self.qec(pos=pos)

    def nFTtdg(self,pos: int):
        anc = self.qc.num_qubits - 1
        self.qc.reset(anc)

        self.qc.h(anc)
        self.qc.tdg(anc)

        self.qc.cx(0+7*pos, anc)
        self.qc.cx(1+7*pos, anc)
        self.qc.cx(2+7*pos, anc)        

        self.qc.measure(anc, 0)

        with self.qc.if_test((0,1)):
            self.qc.reset(anc)
            self.qc.h(anc)
            self.qc.sdg(anc)
            self.qc.cx(0+7*pos, anc)
            self.qc.cx(1+7*pos, anc)
            self.qc.cx(2+7*pos, anc) 
            self.qc.measure(anc, 0)
            with self.qc.if_test((0,1)):
                self.qc.z(0+7*pos)
                self.qc.z(1+7*pos)
                self.qc.z(2+7*pos)

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
                #self.nFTt(pos=pos)
                self.t_cheat(pos=pos)
                # if self.err and self.qec_counter%8==0:
                #     self.qec(pos = pos)
            if i == "tdg":
                #self.nFTtdg(pos=pos)
                self.tdg_cheat(pos=pos)
                # if self.err and self.qec_counter%8==0:
                #     self.qec(pos = pos)
            if i == "h":
                self.h(pos=pos)
            if i == "z":
                self.z(pos=pos)

    def cu(self, gate: list, adjgate: list):
        self.u2(0, gate=gate)
        if self.err:
            self.qec(pos = 0)
        self.u2(1, gate=gate)
        if self.err:
            self.qec(pos = 1)
        self.cnot(control=0, target=1)
        self.u2(1, gate=adjgate)
        # if self.err:
        #     self.qec(pos = 1)
        self.cnot(control=0, target=1)

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

    def err_mitigation(self, setting: bool):
        self.postselection = setting
    
    def __classical_error_correction__(self, bits: list):
        code0 = ["0000000","1010101","0110011","1100110","0001111","1011010","0111100","1101001"]
        code1 = ["1111111","0101010","1001100","0011001","1110000","0100101","1000011","0010110"]
        hmm = 0
        for i,val in enumerate(bits):
            for j in code0:
                diff = 0
                for a,b in zip(val, j):
                    if a!=b:
                        diff += 1
                if diff == 0:
                    hmm += 1
                    break
                if diff == 1:
                    bits[i] = j
                    break
            for j in code1:
                diff = 0
                for a,b in zip(val, j):
                    if a!=b:
                        diff += 1
                if diff == 0:
                    hmm += 1
                    break
                if diff == 1:
                    bits[i] = j
                    break
        return bits

    def readout(self, pos: int, shots: int, p = 0):
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

        self.qc = transpile(self.qc, optimization_level=1)

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
        postprocess = [i[7:allcbits-10] for i in bitstring]                                 #Flags during qec to make it fault tolerant, if at least one strikes, need to discard shot

        #print(bits)
        # print(postprocess)

        for i in range(len(pre)):
            if pre[i].count("1") != 0:
                bits[i] = "pre"
                # print("AHAAAA wieso flag???")
                preselected += hmm[i]

        test_0 = ["0000000","1010101","0110011","1100110","0001111","1011010","0111100","1101001"]
        test_1 = ["1111111","0101010","1001100","0011001","1110000","0100101","1000011","0010110"]

        if self.classical_ec:
            bits = self.__classical_error_correction__(bits)

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
                else:
                    if np.random.rand() < 0.5:
                        bits[i] = 0
                    else:
                        bits[i] = 1

        # for i in range(len(postprocess)):
        #     if postprocess[i].count("1") != 0:
        #         if bits[i] != "pre" and bits[i] != "post":
        #             bits[i] = "post"

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
    def __init__(self, n: int, magic = 0):
        self.n = n

        self.zeros = 0
        self.ones = 0
        self.preselected = 0
        self.post = 0

        self.err = False
        self.postselection = False
        self.classical_ec = False
        self.qec_counter = 0

        self.hadamards = [0 for i in range(n+magic)]

        qr = QuantumRegister(9*(n+magic)+2, "q")
        cbit = ClassicalRegister(4,"c")
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
    
    def s_timo(self, pos: int):
        S_alt = np.diag([1, 1j, 1j, 1, 1j, 1, 1, 1j])
        #threshold = 1e-10
        #T_alt[np.abs(T_alt) < threshold] = np.nan
        s_timo = UnitaryGate(S_alt, label="s_timo")

        if self.hadamards[pos]%2 == 0:  
            self.qc.append(s_timo, [3+9*pos, 4+9*pos, 5+9*pos])            #ich glaube man muss hier die reihenfolge reversen, ist ein 50/50
        else:
            self.qc.append(s_timo, [1+9*pos, 4+9*pos, 7+9*pos])

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

    def sdg_timo(self, pos: int):
        Sdg_alt = np.diag([1, -1j, -1j, 1, -1j, 1, 1, -1j])
        #threshold = 1e-10
        #T_alt[np.abs(T_alt) < threshold] = np.nan
        sdg_timo = UnitaryGate(Sdg_alt, label="s_timo")

        if self.hadamards[pos]%2 == 0:  
            self.qc.append(sdg_timo, [3+9*pos, 4+9*pos, 5+9*pos])            #ich glaube man muss hier die reihenfolge reversen, ist ein 50/50
        else:
            self.qc.append(sdg_timo, [1+9*pos, 4+9*pos, 7+9*pos])

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

    def t_timo(self, pos: int):
        T_alt = np.diag([1, (1+1j)/np.sqrt(2), (1+1j)/np.sqrt(2), 1, (1+1j)/np.sqrt(2), 1, 1, (1+1j)/np.sqrt(2)])
        #threshold = 1e-10
        #T_alt[np.abs(T_alt) < threshold] = np.nan
        t_timo = UnitaryGate(T_alt, label="t_timo")

        if self.hadamards[pos]%2 == 0:  
            self.qc.append(t_timo, [3+9*pos, 4+9*pos, 5+9*pos])            #ich glaube man muss hier die reihenfolge reversen, ist ein 50/50
        else:
            self.qc.append(t_timo, [1+9*pos, 4+9*pos, 7+9*pos])

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

    def tdg_timo(self, pos: int):
        T_alt = np.diag([1, (1+1j)/np.sqrt(2), (1+1j)/np.sqrt(2), 1, (1+1j)/np.sqrt(2), 1, 1, (1+1j)/np.sqrt(2)])
        Tdg_alt = np.conjugate(T_alt)
        #threshold = 1e-10
        #T_alt[np.abs(T_alt) < threshold] = np.nan
        tdg_timo = UnitaryGate(Tdg_alt, label="tdg_timo")

        if self.hadamards[pos]%2 == 0:  
            self.qc.append(tdg_timo, [3+9*pos, 4+9*pos, 5+9*pos])            #ich glaube man muss hier die reihenfolge reversen, ist ein 50/50
        else:
            self.qc.append(tdg_timo, [1+9*pos, 4+9*pos, 7+9*pos])

    def cs(self):
        self.t(pos=0)
        self.t(pos=1)
        self.cnot(control=0, target=1)
        self.tdg(pos=1)
        self.cnot(control=0, target=1)

    def u2(self, pos: int, gate: list):
        for i in gate:
            if i == "s":
                # if self.err == True and self.qec_counter%5==0:
                #     self.qec_zstab(pos=pos)
                self.s_cheat(pos=pos)
            if i == "sdg":
                #  if self.err == True and self.qec_counter%5==0:
                #     self.qec_zstab(pos=pos)
                 self.sdg_cheat(pos=pos)
            if i == "t":
                #  if self.err == True and self.qec_counter%5==0:
                #     self.qec_zstab(pos=pos)
                 self.t_cheat(pos=pos)
            if i == "tdg":
                #  if self.err == True and self.qec_counter%5==0:
                #     self.qec_zstab(pos=pos)
                 self.tdg_cheat(pos=pos)
            if i == "h":
                 self.h(pos=pos)
            if i == "z":
                 self.z(pos=pos)

    def cu(self, Ugates: list, adjUgates: list):
        self.u2(0, Ugates)
        # if self.err:
        #     self.qec(pos=0)
        self.u2(1, Ugates)
        # if self.err:
        #     self.qec(pos=1)
        self.cnot(control=0, target=1)
        self.u2(1, adjUgates)
        # if self.err:
        #     self.qec(pos=1)
        self.cnot(control=0, target=1)

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
            self.qc.id(anc)
            self.qc.measure(anc,0)

            #X0 X1 X3 X4 Stabilizer:
            self.qc.reset(anc)
            self.qc.h(anc)
            self.qc.cx(anc, 1+9*pos)
            self.qc.cx(anc, 4+9*pos)
            self.qc.cx(anc, 0+9*pos)
            self.qc.cx(anc, 3+9*pos)
            self.qc.h(anc)
            self.qc.id(anc)
            self.qc.measure(anc,1)

            #X4 X5 X7 X8 Stabilizer:
            self.qc.reset(anc)
            self.qc.h(anc)
            self.qc.cx(anc, 4+9*pos)
            self.qc.cx(anc, 7+9*pos)
            self.qc.cx(anc, 5+9*pos)
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
            self.qc.cx(anc, 4+9*pos)
            self.qc.cx(anc, 5+9*pos)
            self.qc.cx(anc, 1+9*pos)
            self.qc.cx(anc, 2+9*pos)
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

    def err_mitigation(self, setting: bool):
        self.postselection = setting

    def __classical_error_correction__(self, bits: list):
        code0 = ['000110101', '110110110', '110110101', '110000000', '000110110', '101101101', '011101101', '011011000', '011011011', '110000011', '000000000', '011101110', '101011011', '101101110', '000000011', '101011000']
        code1 = ['010100111', '010010001', '111111111', '001001010', '111001010', '001111111', '100010010', '111111100', '100100100', '100010001', '001001001', '010010010', '100100111', '111001001', '001111100', '010100100']
        hmm = 0
        for i,val in enumerate(bits):
            for j in code0:
                diff = 0
                for a,b in zip(val, j):
                    if a!=b:
                        diff += 1
                if diff == 0:
                    hmm += 1
                    break
                if diff == 1:
                    bits[i] = j
                    break
            for j in code1:
                diff = 0
                for a,b in zip(val, j):
                    if a!=b:
                        diff += 1
                if diff == 0:
                    hmm += 1
                    break
                if diff == 1:
                    bits[i] = j
                    break
        return bits

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
        flags = [i[9:len(bitstring[0])-4] for i in bitstring]

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

class RotSurf16q:
    def __init__(self, n: int):
        self.n = n

        self.zeros = 0
        self.ones = 0
        self.preselected = 0
        self.post = 0

        self.err = False
        self.postselection = False
        self.classical_ec = False

        self.hadamards = [0,0]

        qr = QuantumRegister(16*n, "q")
        cbit = ClassicalRegister(0,"c")
        self.qc = QuantumCircuit(qr, cbit)
        # for i in range(16*n):
        #     self.qc.id(i)
            #1st step
        for i in range(n):
            self.qc.h(0+16*i)
            self.qc.h(2+16*i)
            self.qc.h(5+16*i)
            self.qc.h(7+16*i)
            self.qc.h(8+16*i)
            self.qc.h(10+16*i)
            self.qc.h(13+16*i)
            self.qc.h(15+16*i)

            self.qc.cx(0+16*i, 1+16*i)
            self.qc.cx(2+16*i, 3+16*i)
            self.qc.cx(13+16*i, 12+16*i)
            self.qc.cx(15+16*i, 14+16*i)

            self.qc.cx(5+16*i, 1+16*i)
            self.qc.cx(5+16*i, 6+16*i)

            self.qc.cx(7+16*i, 6+16*i)
            self.qc.cx(7+16*i, 11+16*i)

            self.qc.cx(8+16*i, 4+16*i)
            self.qc.cx(8+16*i, 9+16*i)

            self.qc.cx(10+16*i, 9+16*i)
            self.qc.cx(10+16*i, 14+16*i)
            #### 2nd step
            self.qc.cx(5+16*i, 2+16*i)
            self.qc.cx(10+16*i, 13+16*i)
            #### 3rd step
            self.qc.cx(7+16*i, 10+16*i)
            self.qc.cx(8+16*i, 5+16*i)


            # self.qc.h(1+16*i)
            # self.qc.h(2+16*i)
            # self.qc.h(5+16*i)
            # self.qc.h(7+16*i)
            # self.qc.h(8+16*i)
            # self.qc.h(10+16*i)
            # self.qc.h(13+16*i)
            # self.qc.h(14+16*i)

            # self.qc.cx(1+16*i, 0+16*i)
            # self.qc.cx(2+16*i, 3+16*i)
            # self.qc.cx(13+16*i, 12+16*i)
            # self.qc.cx(14+16*i, 15+16*i)

            # self.qc.cx(7+16*i, 11+16*i)
            # self.qc.cx(7+16*i, 6+16*i)

            # self.qc.cx(8+16*i, 4+16*i)
            # self.qc.cx(8+16*i, 9+16*i)
            # #### 2nd step
            # self.qc.cx(5+16*i, 1+16*i)
            # self.qc.cx(5+16*i, 2+16*i)
            # self.qc.cx(5+16*i, 6+16*i)
            # self.qc.cx(10+16*i, 9+16*i)
            # self.qc.cx(10+16*i, 13+16*i)
            # self.qc.cx(10+16*i, 14+16*i)
            # #### 3rd step
            # self.qc.cx(7+16*i, 10+16*i)
            # self.qc.cx(8+16*i, 5+16*i)

    def x(self, pos: int):
        if self.hadamards[pos]%2==0:
            self.qc.x(0+16*pos)
            self.qc.x(4+16*pos)
            self.qc.x(8+16*pos)
            self.qc.x(12+16*pos)
        else:
            self.qc.x(0+16*pos)
            self.qc.x(1+16*pos)
            self.qc.x(2+16*pos)
            self.qc.x(3+16*pos)

    def z(self, pos: int):
        if self.hadamards[pos]%2==1:
            self.qc.z(0+16*pos)
            self.qc.z(4+16*pos)
            self.qc.z(8+16*pos)
            self.qc.z(12+16*pos)
        else:
            self.qc.z(0+16*pos)
            self.qc.z(1+16*pos)
            self.qc.z(2+16*pos)
            self.qc.z(3+16*pos)

    def h(self, pos: int):
        for i in range(16):
            self.qc.h(i+16*pos)
        self.hadamards[pos] += 1

    def cnot(self, control: int, target: int):
        if control == 0:
            if self.hadamards[0]%2==1 and self.hadamards[1]%2==0:
                for i in range(16):
                    self.qc.h(i+16)
                    self.qc.cz(i, i+16)
                    self.qc.h(i+16)
            elif self.hadamards[0]%2==0 and self.hadamards[1]%2==1:
                for i in range(16):
                    self.qc.h(i+16*target)
                    self.qc.cz(i, i+16)
                    self.qc.h(i+16)
            else:
                for i in range(16):
                    self.qc.cx(i,16+i)
        elif control == 1:
            if self.hadamards[0]%2==0 and self.hadamards[1]%2==1:
                for i in range(16):
                    self.qc.h(i)
                    self.qc.cz(i, i+16)
                    self.qc.h(i)
            elif self.hadamards[0]%2==1 and self.hadamards[1]%2==0:
                for i in range(16):
                    self.qc.h(i)
                    self.qc.cz(i, i+16)
                    self.qc.h(i)
            else: 
                for i in range(16):
                    self.qc.cx(16+i,i)

    # def s(self, pos: int):
    #     anc = self.qc.num_qubits - 1
    #     self.qc.reset(anc)

    #     self.qc.h(anc)
    #     #self.qc.append(h_ideal,[anc])
    #     self.qc.s(anc)

    #     if self.hadamards[pos]%2 == 0:
    #         self.qc.cx(0+16*pos, anc)
    #         self.qc.cx(1+16*pos, anc)
    #         self.qc.cx(2+16*pos, anc)
    #         self.qc.cx(3+16*pos, anc)        
    #     else:
    #         self.qc.cx(0+16*pos, anc)
    #         self.qc.cx(4+16*pos, anc)
    #         self.qc.cx(8+16*pos, anc)
    #         self.qc.cx(12+16*pos, anc)    

    #     self.qc.measure(anc, 0)

    #     if self.hadamards[pos]%2 == 0:
    #         with self.qc.if_test((0,1)):
    #             self.qc.z(0+16*pos)
    #             self.qc.z(1+16*pos)
    #             self.qc.z(2+16*pos)
    #             self.qc.z(3+16*pos)        
    #     else:
    #         with self.qc.if_test((0,1)):
    #             self.qc.z(0+16*pos)
    #             self.qc.z(4+16*pos)
    #             self.qc.z(8+16*pos)
    #             self.qc.z(12+16*pos) 

    def s(self, pos: int):
        S_alt = np.diag([1, 1j, 1j, 1,
            1j, 1, 1, 1j,
            1j, 1, 1, 1j,
            1, 1j, 1j, 1])
        #threshold = 1e-10
        #T_alt[np.abs(T_alt) < threshold] = np.nan
        S_timo = UnitaryGate(S_alt, label="t_timo")

        if self.hadamards[pos]%2 == 0:  
            self.qc.append(S_timo, [0+16*pos, 1+16*pos, 2+16*pos, 3+16*pos])            #ich glaube man muss hier die reihenfolge reversen, ist ein 50/50
        else:
            self.qc.append(S_timo, [0+16*pos, 4+16*pos, 8+16*pos, 12+16*pos])

    # def sdg(self, pos: int):
    #     anc = self.qc.num_qubits - 1
    #     self.qc.reset(anc)

    #     self.qc.h(anc)
    #     #self.qc.append(h_ideal,[anc])
    #     self.qc.sdg(anc)

    #     if self.hadamards[pos]%2 == 0:
    #         self.qc.cx(0+16*pos, anc)
    #         self.qc.cx(1+16*pos, anc)
    #         self.qc.cx(2+16*pos, anc)
    #         self.qc.cx(3+16*pos, anc)        
    #     else:
    #         self.qc.cx(0+16*pos, anc)
    #         self.qc.cx(4+16*pos, anc)
    #         self.qc.cx(8+16*pos, anc)
    #         self.qc.cx(12+16*pos, anc)    

    #     self.qc.measure(anc, 0)

    #     if self.hadamards[pos]%2 == 0:
    #         with self.qc.if_test((0,1)):
    #             self.qc.z(0+16*pos)
    #             self.qc.z(1+16*pos)
    #             self.qc.z(2+16*pos)
    #             self.qc.z(3+16*pos)        
    #     else:
    #         with self.qc.if_test((0,1)):
    #             self.qc.z(0+16*pos)
    #             self.qc.z(4+16*pos)
    #             self.qc.z(8+16*pos)
    #             self.qc.z(12+16*pos) 

    def sdg(self, pos: int):
        Sdg_alt = np.diag([1, -1j, -1j, 1,
            -1j, 1, 1, -1j,
            -1j, 1, 1, -1j,
            1, -1j, -1j, 1])
        #threshold = 1e-10
        #T_alt[np.abs(T_alt) < threshold] = np.nan
        S_timo = UnitaryGate(Sdg_alt, label="t_timo")

        if self.hadamards[pos]%2 == 0:  
            self.qc.append(S_timo, [0+16*pos, 1+16*pos, 2+16*pos, 3+16*pos])            #ich glaube man muss hier die reihenfolge reversen, ist ein 50/50
        else:
            self.qc.append(S_timo, [0+16*pos, 4+16*pos, 8+16*pos, 12+16*pos])

    # def t(self, pos: int):
    #     anc = self.qc.num_qubits - 1
    #     self.qc.reset(anc)

    #     self.qc.h(anc)
    #     #self.qc.append(h_ideal,[anc])
    #     self.qc.t(anc)

    #     if self.hadamards[pos]%2 == 0:
    #         self.qc.cx(0+16*pos, anc)
    #         self.qc.cx(1+16*pos, anc)
    #         self.qc.cx(2+16*pos, anc)
    #         self.qc.cx(3+16*pos, anc)       
    #     else:
    #         self.qc.cx(0+16*pos, anc)
    #         self.qc.cx(4+16*pos, anc)
    #         self.qc.cx(8+16*pos, anc)
    #         self.qc.cx(12+16*pos, anc)     

    #     self.qc.measure(anc, 0)
    #     # if z_stab:
    #     #     z_qec_ideal(qc, had=had, pos=pos)
    #     self.qc.reset(anc)
    #     if self.hadamards[pos]%2 == 0:
    #         with self.qc.if_test((0,1)):
    #             self.qc.h(anc)
    #             self.qc.s(anc)
    #             self.qc.cx(0+16*pos, anc)
    #             self.qc.cx(1+16*pos, anc)
    #             self.qc.cx(2+16*pos, anc)
    #             self.qc.cx(3+16*pos, anc)  
    #             self.qc.measure(anc, 0)
    #             with self.qc.if_test((0,1)):
    #                 self.qc.z(0+16*pos)
    #                 self.qc.z(1+16*pos)
    #                 self.qc.z(2+16*pos)
    #                 self.qc.z(3+16*pos)  
    #     else:
    #         with self.qc.if_test((0,1)):
    #             self.qc.h(anc)
    #             self.qc.s(anc)
    #             self.qc.cx(0+16*pos, anc)
    #             self.qc.cx(4+16*pos, anc)
    #             self.qc.cx(8+16*pos, anc)
    #             self.qc.cx(12+16*pos, anc)   
    #             self.qc.measure(anc, 0)
    #             with self.qc.if_test((0,1)):
    #                 self.qc.z(0+16*pos)
    #                 self.qc.z(4+16*pos)
    #                 self.qc.z(8+16*pos)
    #                 self.qc.z(12+16*pos)   
    
    def t(self, pos: int):
        T_alt = np.diag([1, (1+1j)/np.sqrt(2), (1+1j)/np.sqrt(2), 1,
            (1+1j)/np.sqrt(2), 1, 1, (1+1j)/np.sqrt(2),
            (1+1j)/np.sqrt(2), 1, 1, (1+1j)/np.sqrt(2),
            1, (1+1j)/np.sqrt(2), (1+1j)/np.sqrt(2), 1])
        #threshold = 1e-10
        #T_alt[np.abs(T_alt) < threshold] = np.nan
        t_timo = UnitaryGate(T_alt, label="t_timo")

        if self.hadamards[pos]%2 == 0:  
            self.qc.append(t_timo, [0+16*pos, 1+16*pos, 2+16*pos, 3+16*pos])            #ich glaube man muss hier die reihenfolge reversen, ist ein 50/50
        else:
            self.qc.append(t_timo, [0+16*pos, 4+16*pos, 8+16*pos, 12+16*pos])

    # def tdg(self, pos: int):
    #     anc = self.qc.num_qubits - 1
    #     self.qc.reset(anc)

    #     self.qc.h(anc)
    #     #self.qc.append(h_ideal,[anc])
    #     self.qc.tdg(anc)

    #     if self.hadamards[pos]%2 == 0:
    #         self.qc.cx(0+16*pos, anc)
    #         self.qc.cx(1+16*pos, anc)
    #         self.qc.cx(2+16*pos, anc)
    #         self.qc.cx(3+16*pos, anc)       
    #     else:
    #         self.qc.cx(0+16*pos, anc)
    #         self.qc.cx(4+16*pos, anc)
    #         self.qc.cx(8+16*pos, anc)
    #         self.qc.cx(12+16*pos, anc)     

    #     self.qc.measure(anc, 0)
    #     # if z_stab:
    #     #     z_qec_ideal(qc, had=had, pos=pos)

    #     if self.hadamards[pos]%2 == 0:
    #         with self.qc.if_test((0,1)):
    #             self.qc.reset(anc)
    #             self.qc.h(anc)
    #             self.qc.sdg(anc)
    #             self.qc.cx(0+16*pos, anc)
    #             self.qc.cx(1+16*pos, anc)
    #             self.qc.cx(2+16*pos, anc)
    #             self.qc.cx(3+16*pos, anc)  
    #             self.qc.measure(anc, 0)
    #             with self.qc.if_test((0,1)):
    #                 self.qc.z(0+16*pos)
    #                 self.qc.z(1+16*pos)
    #                 self.qc.z(2+16*pos)
    #                 self.qc.z(3+16*pos)  
    #     else:
    #         with self.qc.if_test((0,1)):
    #             self.qc.reset(anc)
    #             self.qc.h(anc)
    #             self.qc.sdg(anc)
    #             self.qc.cx(0+16*pos, anc)
    #             self.qc.cx(4+16*pos, anc)
    #             self.qc.cx(8+16*pos, anc)
    #             self.qc.cx(12+16*pos, anc)   
    #             self.qc.measure(anc, 0)
    #             with self.qc.if_test((0,1)):
    #                 self.qc.z(0+16*pos)
    #                 self.qc.z(4+16*pos)
    #                 self.qc.z(8+16*pos)
    #                 self.qc.z(12+16*pos)   

    def tdg(self, pos: int):
        T_alt = np.diag([1, (1+1j)/np.sqrt(2), (1+1j)/np.sqrt(2), 1,
            (1+1j)/np.sqrt(2), 1, 1, (1+1j)/np.sqrt(2),
            (1+1j)/np.sqrt(2), 1, 1, (1+1j)/np.sqrt(2),
            1, (1+1j)/np.sqrt(2), (1+1j)/np.sqrt(2), 1])
        # threshold = 1e-10
        # T_alt[np.abs(T_alt) < threshold] = np.nan
        Tdg_alt = np.conjugate(T_alt)
        tdg_timo = UnitaryGate(Tdg_alt, label="t_timo")

        if self.hadamards[pos]%2 == 0:
            self.qc.append(tdg_timo, [0+16*pos, 1+16*pos, 2+16*pos, 3+16*pos])
        else:
            self.qc.append(tdg_timo, [0+16*pos, 4+16*pos, 8+16*pos, 12+16*pos])

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

    def cu(self, Ugates: list, adjUgates: list):
        self.u2(0, Ugates)
        # if self.err:
        #     self.qec(pos=0)
        self.u2(1, Ugates)
        # if self.err:
        #     self.qec(pos=1)
        self.cnot(control=0, target=1)
        self.u2(1, adjUgates)
        # if self.err:
        #     self.qec(pos=1)
        self.cnot(control=0, target=1)

    def qec(self, pos: int):
        anc = self.qc.num_qubits - 1
        if self.hadamards[pos]%2==1:
            #Z0 Z1 Stabilizer:
            self.qc.reset(anc)
            self.qc.cx(0+16*pos, anc)
            self.qc.cx(1+16*pos, anc)
            self.qc.id(anc)
            self.qc.measure(anc,0)

            #Z2 Z3 Stabilizer:
            self.qc.reset(anc)
            self.qc.cx(2+16*pos, anc)
            self.qc.cx(3+16*pos, anc)
            self.qc.id(anc)
            self.qc.measure(anc,1)

            #Z1 Z2 Z5 Z6 Stabilizer:
            self.qc.reset(anc)
            self.qc.cx(5+16*pos, anc)
            self.qc.cx(6+16*pos, anc)
            self.qc.cx(1+16*pos, anc)
            self.qc.cx(2+16*pos, anc)
            self.qc.id(anc)
            self.qc.measure(anc,2)

            #Z4 Z5 Z8 Z9 Stabilizer:
            self.qc.reset(anc)
            self.qc.cx(4+16*pos, anc)
            self.qc.cx(9+16*pos, anc)
            self.qc.cx(8+16*pos, anc)
            self.qc.cx(5+16*pos, anc)
            self.qc.id(anc)
            self.qc.measure(anc,3)

            #Z6 Z7 Z10 Z11 Stabilizer:
            self.qc.reset(anc)
            self.qc.cx(6+16*pos, anc)
            self.qc.cx(11+16*pos, anc)
            self.qc.cx(10+16*pos, anc)
            self.qc.cx(7+16*pos, anc)
            self.qc.id(anc)
            self.qc.measure(anc,4)


            #Z9 Z10 Z13 Z14 Stabilizer:
            self.qc.reset(anc)
            self.qc.cx(9+16*pos, anc)
            self.qc.cx(10+16*pos, anc)
            self.qc.cx(13+16*pos, anc)
            self.qc.cx(14+16*pos, anc)
            self.qc.id(anc)
            self.qc.measure(anc,5)


            #Z12 Z13 Stabilizer:
            self.qc.reset(anc)
            self.qc.cx(12+16*pos, anc)
            self.qc.cx(13+16*pos, anc)
            self.qc.id(anc)
            self.qc.measure(anc,6)


            #Z14 Z15 Stabilizer:
            self.qc.reset(anc)
            self.qc.cx(14+16*pos, anc)
            self.qc.cx(15+16*pos, anc)
            self.qc.id(anc)
            self.qc.measure(anc,7)


            with self.qc.if_test((0,1)):             #0
                with self.qc.if_test((2,0)):    
                    self.qc.x(0+16*pos)

            with self.qc.if_test((0,1)):             #1
                with self.qc.if_test((2,1)):
                    with self.qc.if_test((1,0)):
                        with self.qc.if_test((3,0)):
                            with self.qc.if_test((4,0)):
                                self.qc.x(1+16*pos)

            with self.qc.if_test((1,1)):             #2
                with self.qc.if_test((2,1)):
                    with self.qc.if_test((0,0)):
                        with self.qc.if_test((3,0)):
                            with self.qc.if_test((4,0)):
                                self.qc.x(2+16*pos)

            with self.qc.if_test((1,1)):             #3
                with self.qc.if_test((2,0)):
                    self.qc.x(3+16*pos)

            with self.qc.if_test((3,1)):             #4 und 8
                with self.qc.if_test((2,0)):
                    with self.qc.if_test((5,0)):
                        self.qc.x(4+16*pos)

            with self.qc.if_test((2,1)):             #5
                with self.qc.if_test((3,1)):
                    with self.qc.if_test((0,0)):
                        with self.qc.if_test((5,0)):
                            with self.qc.if_test((4,0)):
                                with self.qc.if_test((1,0)):
                                    self.qc.x(5+16*pos)
            
            with self.qc.if_test((2,1)):             #6
                with self.qc.if_test((4,1)):
                    with self.qc.if_test((0,0)):
                        with self.qc.if_test((1,0)):
                            with self.qc.if_test((5,0)):
                                with self.qc.if_test((3,0)):
                                    self.qc.x(6+16*pos)

            with self.qc.if_test((4,1)):             #7 und 11
                with self.qc.if_test((2,0)):
                    with self.qc.if_test((5,0)):
                        self.qc.x(7+16*pos)
            
            with self.qc.if_test((5,1)):             #9
                with self.qc.if_test((3,1)):
                    with self.qc.if_test((6,0)): 
                        with self.qc.if_test((7,0)):
                            with self.qc.if_test((2,0)):
                                with self.qc.if_test((4,0)):
                                    self.qc.x(9+16*pos)

            with self.qc.if_test((5,1)):             #10
                with self.qc.if_test((4,1)):
                    with self.qc.if_test((6,0)):
                        with self.qc.if_test((2,0)):
                            with self.qc.if_test((7,0)):
                                with self.qc.if_test((3,0)):
                                    self.qc.x(10+16*pos)
            
            with self.qc.if_test((6,1)):             #12
                with self.qc.if_test((5,0)):
                    self.qc.x(12+16*pos)

            with self.qc.if_test((5,1)):             #13
                with self.qc.if_test((6,1)):
                    with self.qc.if_test((3,0)):
                        with self.qc.if_test((4,0)):
                            with self.qc.if_test((7,0)):
                                self.qc.x(13+16*pos)

            with self.qc.if_test((7,1)):             #14
                with self.qc.if_test((5,1)):
                    with self.qc.if_test((3,0)):
                        with self.qc.if_test((4,0)):
                            with self.qc.if_test((6,0)):
                                self.qc.x(14+16*pos)

            with self.qc.if_test((7,1)):             #15
                with self.qc.if_test((5,0)):
                    self.qc.x(15+16*pos)

            with self.qc.if_test((2,1)):             #special case for hook error on X5 X6 X9 X10 Stabilizer
                with self.qc.if_test((3,1)):
                    with self.qc.if_test((4,1)):
                        with self.qc.if_test((5,1)):
                            with self.qc.if_test((0,0)):
                                with self.qc.if_test((1,0)):
                                    with self.qc.if_test((6,0)):
                                        with self.qc.if_test((7,0)):
                                            self.qc.x(5+16*pos), self.qc.x(10+16*pos)
        ###########################################################################################################

            #X4 X8 Stabilizer:
            self.qc.reset(anc)
            self.qc.h(anc)
            self.qc.cx(anc, 4+16*pos)
            self.qc.cx(anc, 8+16*pos)
            self.qc.h(anc)
            self.qc.id(anc)
            self.qc.measure(anc,0)

            #X0 X1 X4 X5 Stabilizer:
            self.qc.reset(anc)
            self.qc.h(anc)
            self.qc.cx(anc, 5+16*pos)
            self.qc.cx(anc, 1+16*pos)
            self.qc.cx(anc, 4+16*pos)
            self.qc.cx(anc, 0+16*pos)
            self.qc.h(anc)
            self.qc.id(anc)
            self.qc.measure(anc,1)
        
            #X2 X3 X6 X7 Stabilizer:
            self.qc.reset(anc)
            self.qc.h(anc)
            self.qc.cx(anc, 2+16*pos)
            self.qc.cx(anc, 6+16*pos)
            self.qc.cx(anc, 3+16*pos)
            self.qc.cx(anc, 7+16*pos)
            self.qc.h(anc)
            self.qc.id(anc)
            self.qc.measure(anc,2)

            #X5 X6 X9 X10 Stabilizer:
            self.qc.reset(anc)
            self.qc.h(anc)
            self.qc.cx(anc, 5+16*pos)
            self.qc.cx(anc, 10+16*pos)
            self.qc.cx(anc, 9+16*pos)
            self.qc.cx(anc, 6+16*pos)
            self.qc.h(anc)
            self.qc.id(anc)
            self.qc.measure(anc,3)

            #X8 X9 X12 X13 Stabilizer:
            self.qc.reset(anc)
            self.qc.h(anc)
            self.qc.cx(anc, 13+16*pos)
            self.qc.cx(anc, 9+16*pos)
            self.qc.cx(anc, 12+16*pos)
            self.qc.cx(anc, 8+16*pos)
            self.qc.h(anc)
            self.qc.id(anc)
            self.qc.measure(anc,4)

            #X10 X11 X14 X15 Stabilizer:
            self.qc.reset(anc)
            self.qc.h(anc)
            self.qc.cx(anc, 10+16*pos)
            self.qc.cx(anc, 14+16*pos)
            self.qc.cx(anc, 11+16*pos)
            self.qc.cx(anc, 15+16*pos)
            self.qc.h(anc)
            self.qc.id(anc)
            self.qc.measure(anc,5)

            #X7 X11 Stabilizer:
            self.qc.reset(anc)
            self.qc.h(anc)
            self.qc.cx(anc, 7+16*pos)
            self.qc.cx(anc, 11+16*pos)
            self.qc.h(anc)
            self.qc.id(anc)
            self.qc.measure(anc,6)
            
            with self.qc.if_test((1,1)):             #0 und 1
                with self.qc.if_test((0,0)):
                    with self.qc.if_test((3,0)):
                        self.qc.z(0+16*pos)

            with self.qc.if_test((2,1)):             #2 und 3
                with self.qc.if_test((3,0)):
                    with self.qc.if_test((6,0)):
                        self.qc.z(2+16*pos)

            with self.qc.if_test((0,1)):             #4
                with self.qc.if_test((1,1)):
                    with self.qc.if_test((4,0)):
                        with self.qc.if_test((3,0)):
                            self.qc.z(4+16*pos)

            with self.qc.if_test((1,1)):             #5
                with self.qc.if_test((3,1)):
                    with self.qc.if_test((0,0)):
                        with self.qc.if_test((2,0)):
                            with self.qc.if_test((4,0)):
                                with self.qc.if_test((5,0)):
                                    self.qc.z(5+16*pos)

            with self.qc.if_test((2,1)):             #6
                with self.qc.if_test((3,1)):
                    with self.qc.if_test((1,0)):
                        with self.qc.if_test((6,0)):
                            with self.qc.if_test((4,0)):
                                with self.qc.if_test((5,0)):
                                    self.qc.z(6+16*pos)

            with self.qc.if_test((2,1)):             #7
                with self.qc.if_test((6,1)):
                    with self.qc.if_test((3,0)):
                        with self.qc.if_test((5,0)):
                            self.qc.z(7+16*pos)

            with self.qc.if_test((0,1)):             #8
                with self.qc.if_test((4,1)):
                    with self.qc.if_test((1,0)):
                        with self.qc.if_test((3,0)):
                            self.qc.z(7+16*pos)
            
            with self.qc.if_test((3,1)):             #9
                with self.qc.if_test((4,1)):
                    with self.qc.if_test((0,0)):
                        with self.qc.if_test((1,0)):
                            with self.qc.if_test((2,0)):
                                with self.qc.if_test((5,0)):
                                    self.qc.z(9+16*pos)

            with self.qc.if_test((3,1)):             #10
                with self.qc.if_test((5,1)):
                    with self.qc.if_test((4,0)):
                        with self.qc.if_test((1,0)):
                            with self.qc.if_test((2,0)):
                                with self.qc.if_test((6,0)):
                                    self.qc.z(10+16*pos)
            
            with self.qc.if_test((5,1)):             #11
                with self.qc.if_test((6,1)):
                    with self.qc.if_test((2,0)):
                        with self.qc.if_test((3,0)):
                            self.qc.z(11+16*pos)
            
            with self.qc.if_test((4,1)):             #12 und 13
                with self.qc.if_test((0,0)):
                    with self.qc.if_test((3,0)):
                        self.qc.z(12+16*pos)

            with self.qc.if_test((5,1)):             #14 und 15
                with self.qc.if_test((3,0)):
                    with self.qc.if_test((6,0)):
                        self.qc.z(14+16*pos)

            with self.qc.if_test((2,1)):             #special case for hook error on Z6 Z7 Z10 Z11 Stabilizer
                with self.qc.if_test((3,1)):
                    with self.qc.if_test((6,1)):
                        with self.qc.if_test((5,1)):
                            with self.qc.if_test((4,0)):
                                with self.qc.if_test((1,0)):
                                    self.qc.z(6+16*pos), self.qc.z(11+16*pos)

            with self.qc.if_test((0,1)):             #special case for hook error on Z4 Z5 Z8 Z9 Stabilizer
                with self.qc.if_test((3,1)):
                    with self.qc.if_test((1,1)):
                        with self.qc.if_test((4,1)):
                            with self.qc.if_test((2,0)):
                                with self.qc.if_test((5,0)):
                                    self.qc.z(4+16*pos), self.qc.z(9+16*pos)

        else:
            #X0 X1 Stabilizer:
            self.qc.reset(anc)
            self.qc.h(anc)
            self.qc.cx(anc, 0+16*pos)
            self.qc.cx(anc, 1+16*pos)
            self.qc.h(anc)
            self.qc.id(anc)
            self.qc.measure(anc,0)

            #X2 X3 Stabilizer:
            self.qc.reset(anc)
            self.qc.h(anc)
            self.qc.cx(anc, 2+16*pos)
            self.qc.cx(anc, 3+16*pos)
            self.qc.h(anc)
            self.qc.id(anc)
            self.qc.measure(anc,1)

            #X1 X2 X5 X6 Stabilizer:
            self.qc.reset(anc)
            self.qc.h(anc)
            self.qc.cx(anc, 5+16*pos)
            self.qc.cx(anc, 6+16*pos)
            self.qc.cx(anc, 1+16*pos)
            self.qc.cx(anc, 2+16*pos)
            self.qc.h(anc)
            self.qc.id(anc)
            self.qc.measure(anc,2)


            #X4 X5 X8 X9 Stabilizer:
            self.qc.reset(anc)
            self.qc.h(anc)
            self.qc.cx(anc, 8+16*pos)
            self.qc.cx(anc, 5+16*pos)
            self.qc.cx(anc, 9+16*pos)
            self.qc.cx(anc, 4+16*pos)
            self.qc.h(anc)
            self.qc.id(anc)
            self.qc.measure(anc,3)

            #X6 X7 X10 X11 Stabilizer:
            self.qc.reset(anc)
            self.qc.h(anc)
            self.qc.cx(anc, 6+16*pos)
            self.qc.cx(anc, 11+16*pos)
            self.qc.cx(anc, 7+16*pos)
            self.qc.cx(anc, 10+16*pos)
            self.qc.h(anc)
            self.qc.id(anc)
            self.qc.measure(anc,4)

            #X9 X10 X13 X14 Stabilizer:
            self.qc.reset(anc)
            self.qc.h(anc)
            self.qc.cx(anc, 9+16*pos)
            self.qc.cx(anc, 10+16*pos)
            self.qc.cx(anc, 13+16*pos)
            self.qc.cx(anc, 14+16*pos)
            self.qc.h(anc)
            self.qc.id(anc)
            self.qc.measure(anc,5)

            #X12 X13 Stabilizer:
            self.qc.reset(anc)
            self.qc.h(anc)
            self.qc.cx(anc, 12+16*pos)
            self.qc.cx(anc, 13+16*pos)
            self.qc.h(anc)
            self.qc.id(anc)
            self.qc.measure(anc,6)

            #X14 X15 Stabilizer:
            self.qc.reset(anc)
            self.qc.h(anc)
            self.qc.cx(anc, 14+16*pos)
            self.qc.cx(anc, 15+16*pos)
            self.qc.h(anc)
            self.qc.id(anc)
            self.qc.measure(anc,7)

            with self.qc.if_test((0,1)):             #0
                with self.qc.if_test((2,0)):    
                    self.qc.z(0+16*pos)

            with self.qc.if_test((0,1)):             #1
                with self.qc.if_test((2,1)):
                    with self.qc.if_test((1,0)):
                        with self.qc.if_test((3,0)):
                            with self.qc.if_test((4,0)):
                                self.qc.z(1+16*pos)

            with self.qc.if_test((1,1)):             #2
                with self.qc.if_test((2,1)):
                    with self.qc.if_test((0,0)):
                        with self.qc.if_test((3,0)):
                            with self.qc.if_test((4,0)):
                                self.qc.z(2+16*pos)

            with self.qc.if_test((1,1)):             #3
                with self.qc.if_test((2,0)):
                    self.qc.z(3+16*pos)

            with self.qc.if_test((3,1)):             #4 und 8
                with self.qc.if_test((2,0)):
                    with self.qc.if_test((5,0)):
                        self.qc.z(4+16*pos)

            with self.qc.if_test((2,1)):             #5
                with self.qc.if_test((3,1)):
                    with self.qc.if_test((0,0)):
                        with self.qc.if_test((5,0)):
                            with self.qc.if_test((4,0)):
                                with self.qc.if_test((1,0)):
                                    self.qc.z(5+16*pos)
            
            with self.qc.if_test((2,1)):             #6
                with self.qc.if_test((4,1)):
                    with self.qc.if_test((0,0)):
                        with self.qc.if_test((1,0)):
                            with self.qc.if_test((5,0)):
                                with self.qc.if_test((3,0)):
                                    self.qc.z(6+16*pos)

            with self.qc.if_test((4,1)):             #7 und 11
                with self.qc.if_test((2,0)):
                    with self.qc.if_test((5,0)):
                        self.qc.z(7+16*pos)
            
            with self.qc.if_test((5,1)):             #9
                with self.qc.if_test((3,1)):
                    with self.qc.if_test((6,0)): 
                        with self.qc.if_test((7,0)):
                            with self.qc.if_test((2,0)):
                                with self.qc.if_test((4,0)):
                                    self.qc.z(9+16*pos)

            with self.qc.if_test((5,1)):             #10
                with self.qc.if_test((4,1)):
                    with self.qc.if_test((6,0)):
                        with self.qc.if_test((2,0)):
                            with self.qc.if_test((7,0)):
                                with self.qc.if_test((3,0)):
                                    self.qc.z(10+16*pos)
            
            with self.qc.if_test((6,1)):             #12
                with self.qc.if_test((5,0)):
                    self.qc.z(12+16*pos)

            with self.qc.if_test((5,1)):             #13
                with self.qc.if_test((6,1)):
                    with self.qc.if_test((3,0)):
                        with self.qc.if_test((4,0)):
                            with self.qc.if_test((7,0)):
                                self.qc.z(13+16*pos)

            with self.qc.if_test((7,1)):             #14
                with self.qc.if_test((5,1)):
                    with self.qc.if_test((3,0)):
                        with self.qc.if_test((4,0)):
                            with self.qc.if_test((6,0)):
                                self.qc.z(14+16*pos)

            with self.qc.if_test((7,1)):             #15
                with self.qc.if_test((5,0)):
                    self.qc.z(15+16*pos)

            with self.qc.if_test((2,1)):             #special case for hook error on Z5 Z6 Z9 Z10 Stabilizer
                with self.qc.if_test((3,1)):
                    with self.qc.if_test((4,1)):
                        with self.qc.if_test((5,1)):
                            with self.qc.if_test((0,0)):
                                with self.qc.if_test((1,0)):
                                    with self.qc.if_test((6,0)):
                                        with self.qc.if_test((7,0)):
                                            self.qc.z(5+16*pos), self.qc.z(10+16*pos)

        ###########################################################################################################

            #Z4 Z8 Stabilizer:
            self.qc.reset(anc)
            self.qc.cx(4+16*pos, anc)
            self.qc.cx(8+16*pos, anc)
            self.qc.id(anc)
            self.qc.measure(anc,0)

            #Z0 Z1 Z4 Z5 Stabilizer:
            self.qc.reset(anc)
            self.qc.cx(1+16*pos, anc)
            self.qc.cx(5+16*pos, anc)
            self.qc.cx(0+16*pos, anc)
            self.qc.cx(4+16*pos, anc)
            self.qc.id(anc)
            self.qc.measure(anc,1)
        
            #Z2 Z3 Z6 Z7 Stabilizer:
            self.qc.reset(anc)
            self.qc.cx(2+16*pos, anc)
            self.qc.cx(6+16*pos, anc)
            self.qc.cx(3+16*pos, anc)
            self.qc.cx(7+16*pos, anc)
            self.qc.id(anc)
            self.qc.measure(anc,2)

            #Z5 Z6 Z9 Z10 Stabilizer:
            self.qc.reset(anc)
            self.qc.cx(5+16*pos, anc)
            self.qc.cx(10+16*pos, anc)
            self.qc.cx(9+16*pos, anc)
            self.qc.cx(6+16*pos, anc)
            self.qc.id(anc)
            self.qc.measure(anc,3)

            #Z8 Z9 Z12 Z13 Stabilizer:
            self.qc.reset(anc)
            self.qc.cx(13+16*pos, anc)
            self.qc.cx(9+16*pos, anc)
            self.qc.cx(12+16*pos, anc)
            self.qc.cx(8+16*pos, anc)
            self.qc.id(anc)
            self.qc.measure(anc,4)

            #Z10 Z11 Z14 Z15 Stabilizer:
            self.qc.reset(anc)
            self.qc.cx(10+16*pos, anc)
            self.qc.cx(14+16*pos, anc)
            self.qc.cx(11+16*pos, anc)
            self.qc.cx(15+16*pos, anc)
            self.qc.id(anc)
            self.qc.measure(anc,5)

            #Z7 Z11 Stabilizer:
            self.qc.reset(anc)
            self.qc.cx(7+16*pos, anc)
            self.qc.cx(11+16*pos, anc)
            self.qc.id(anc)
            self.qc.measure(anc,6)
            
            with self.qc.if_test((1,1)):             #0 und 1
                with self.qc.if_test((0,0)):
                    with self.qc.if_test((3,0)):
                        self.qc.x(0+16*pos)

            with self.qc.if_test((2,1)):             #2 und 3
                with self.qc.if_test((3,0)):
                    with self.qc.if_test((6,0)):
                        self.qc.x(2+16*pos)

            with self.qc.if_test((0,1)):             #4
                with self.qc.if_test((1,1)):
                    with self.qc.if_test((4,0)):
                        with self.qc.if_test((3,0)):
                            self.qc.x(4+16*pos)

            with self.qc.if_test((1,1)):             #5
                with self.qc.if_test((3,1)):
                    with self.qc.if_test((0,0)):
                        with self.qc.if_test((2,0)):
                            with self.qc.if_test((4,0)):
                                with self.qc.if_test((5,0)):
                                    self.qc.x(5+16*pos)

            with self.qc.if_test((2,1)):             #6
                with self.qc.if_test((3,1)):
                    with self.qc.if_test((1,0)):
                        with self.qc.if_test((6,0)):
                            with self.qc.if_test((4,0)):
                                with self.qc.if_test((5,0)):
                                    self.qc.x(6+16*pos)

            with self.qc.if_test((2,1)):             #7
                with self.qc.if_test((6,1)):
                    with self.qc.if_test((3,0)):
                        with self.qc.if_test((5,0)):
                            self.qc.x(7+16*pos)

            with self.qc.if_test((0,1)):             #8
                with self.qc.if_test((4,1)):
                    with self.qc.if_test((1,0)):
                        with self.qc.if_test((3,0)):
                            self.qc.x(7+16*pos)
            
            with self.qc.if_test((3,1)):             #9
                with self.qc.if_test((4,1)):
                    with self.qc.if_test((0,0)):
                        with self.qc.if_test((1,0)):
                            with self.qc.if_test((2,0)):
                                with self.qc.if_test((5,0)):
                                    self.qc.x(9+16*pos)

            with self.qc.if_test((3,1)):             #10
                with self.qc.if_test((5,1)):
                    with self.qc.if_test((4,0)):
                        with self.qc.if_test((1,0)):
                            with self.qc.if_test((2,0)):
                                with self.qc.if_test((6,0)):
                                    self.qc.x(10+16*pos)
            
            with self.qc.if_test((5,1)):             #11
                with self.qc.if_test((6,1)):
                    with self.qc.if_test((2,0)):
                        with self.qc.if_test((3,0)):
                            self.qc.x(11+16*pos)
            
            with self.qc.if_test((4,1)):             #12 und 13
                with self.qc.if_test((0,0)):
                    with self.qc.if_test((3,0)):
                        self.qc.x(12+16*pos)

            with self.qc.if_test((5,1)):             #14 und 15
                with self.qc.if_test((3,0)):
                    with self.qc.if_test((6,0)):
                        self.qc.x(14+16*pos)

            with self.qc.if_test((2,1)):             #special case for hook error on X6 X7 X10 X11 Stabilizer
                with self.qc.if_test((3,1)):
                    with self.qc.if_test((6,1)):
                        with self.qc.if_test((5,1)):
                            with self.qc.if_test((4,0)):
                                with self.qc.if_test((1,0)):
                                    self.qc.x(6+16*pos), self.qc.x(11+16*pos)

            with self.qc.if_test((0,1)):             #special case for hook error on X4 X5 X8 X9 Stabilizer
                with self.qc.if_test((3,1)):
                    with self.qc.if_test((1,1)):
                        with self.qc.if_test((4,1)):
                            with self.qc.if_test((2,0)):
                                with self.qc.if_test((5,0)):
                                    self.qc.x(4+16*pos), self.qc.x(9+16*pos)

    def __classical_error_correction__(self, bits: list):
        code0 = ['0110010100110000', '1111000001100110', '1111001101010101', '0011110011000011', '0011001100110000', '0011001100111111', '1001011001101001', '0011111110010101', '0101011001100110', '0000000000001111', '1010010100111100', '1010010101010101', '0011001101011001', '1100111110011010', '1100111110010101', '1100001101010110', '1111110010100110', '0101010100111100', '1010100110011010', '1010011001100101', '0011000001100101', '1111110011001100', '1100111110011001', '1100000001100101', '0000111111111111', '0011000000000011', '1001100110010101', '0110010101010110', '1010010100110011', '1010101011001100', '1001011001100101', '0000001101010110', '0101101011000000', '0011001101011010', '0000110010100110', '1001101010100101', '0000000000000011', '0011111111110000', '1010010101010110', '0110011001100110', '1111111111111100', '1100001100110011', '1010101011000000', '0000000001100101', '0110100110011001', '1001011001100110', '0011001100111100', '0101100110011001', '1100110010100101', '0011110010101001', '0101010100110011', '0110100111111100', '1100001101011001', '1100110010101001', '0000001100111111', '1010101010100101', '1001010100110000', '1100111111110011', '0110101011000000', '0101100110010101', '0101011001100101', '0000111110011010', '0000001100110000', '1010101010101001', '0000001101011010', '1001101010101001', '1111110011000011', '1001101011000000', '0101011001101010', '0101100110010110', '1111110010101010', '0110010100111100', '1100000000000000', '1001010101011010', '0110100110010101', '1010101011000011', '0101101011001100', '0110101011001111', '0110011000001100', '0110101010100110', '0110100110010110', '1100000000001111', '0101101010100101', '0110011000001111', '0000111110010110', '0011000000001111', '1010011000000000', '0110100111110011', '1100001101011010', '1001011001101010', '1010100111111100', '0011111111111111', '1010010101011001', '0000001101011001', '1111000000000000', '0101010101011001', '1100111110010110', '1100000000001100', '1001101010100110', '0011111110010110', '0000110010101010', '0101010100111111', '0011111110011010', '1010011000001100', '1001010100111100', '0000000001101010', '0000111111110011', '0110100111110000', '1001010101011001', '1100110010101010', '0110011001101010', '0101100111110011', '1010100111110011', '1001011000000011', '0101011000001100', '0000001100110011', '1001100111110011', '0000110011001100', '0011110011001111', '0011111110011001', '0101101010101010', '1111111110011001', '0011001100110011', '1001011000001111', '1010101010100110', '1100110011001100', '0101101010100110', '1010011000000011', '1111001100111100', '0101101011001111', '0101100111111100', '0101100110011010', '0110010100111111', '1010100110010110', '0000001100111100', '0110101010101001', '1001010101010110', '1010011001101001', '1001100111110000', '1001100110011010', '0000000000000000', '1111001101011010', '0101011000000000', '1111001101010110', '1100110011000011', '0101100111111111', '0000111111111100', '0110011000000000', '1010100111111111', '0011110011001100', '0110011001100101', '1010010100110000', '0101101011000011', '1100000001100110', '1010010101011010', '0110101011000011', '1100110010100110', '0110100111111111', '1111000001100101', '1010101011001111', '0000001101010101', '0110010101011010', '1111110010100101', '0000110011001111', '0110100110011010', '0110101010100101', '1111001100110011', '0101011000001111', '1111111110010101', '1001101011000011', '1100001101010101', '0011111111110011', '1010100111110000', '1001010101010101', '1001100110010110', '0101010100110000', '1111001100110000', '1111110010101001', '1001100111111111', '1001011000000000', '1001010100111111', '1010011001100110', '0011110010100101', '0110010100110011', '1010011001101010', '0000111110010101', '1100000001101001', '0101011000000011', '1100111111111111', '1111110011000000', '0000000001101001', '1001101011001100', '0110010101010101', '0101010101010101', '1111000001101010', '1111111111110011', '1111000000001111', '1010010100111111', '0011000001100110', '0101101010101001', '0011000000000000', '1111001100111111', '0110010101011001', '0000000000001100', '1001101011001111', '0011110010100110', '0101100111110000', '0011110010101010', '0000110011000011', '1010011000001111', '0101010101010110', '0101011001101001', '1111000001101001', '1001010100110011', '1100001100111100', '0011000000001100', '0011110011000000', '1100000001101010', '0011001101010110', '0000111110011001', '1100110011000000', '1100000000000011', '1100111111111100', '1111001101011001', '0011000001101010', '1100111111110000', '1001101010101010', '1111111110010110', '0011111111111100', '1111110011001111', '1010100110011001', '1111111111111111', '0110011000000011', '1111111111110000', '1100110011001111', '1111000000000011', '0110011001101001', '1100001100111111', '0011000001101001', '1001100111111100', '0000110010101001', '1001100110011001', '0101010101011010', '1100001100110000', '0011001101010101', '0000110011000000', '1001011000001100', '1010100110010101', '1010101010101010', '0110101010101010', '0000111111110000', '0110101011001100', '0000000001100110', '1111111110011010', '0000110010100101', '1111000000001100']
        code1 = ['0111100011101110', '0111100010001011', '0100100011101110', '0010000100011110', '1101111011100010', '0100011101110111', '0010110111010001', '0111100011101101', '1110000101111011', '1101110110111011', '0001111010001000', '0001110111011101', '1101001001000111', '0010111010000111', '1101000101110100', '1000010000101110', '0111010000101110', '0001110110111011', '1101111010000111', '1011010000101101', '0100011101111011', '1000100011101110', '0111101110111011', '1101001000101101', '1101000101111000', '1101110111010010', '0010000101110111', '0010001000100010', '1000101111011101', '0001111011101110', '0010000101110100', '1011100011101110', '1110110111010010', '0001110110110111', '1110111011101101', '1110111011100001', '0010001000101110', '1101111010000100', '0111011100011110', '1101001001001011', '1110001000100010', '1101111011101101', '1101110111010001', '1000011101111000', '0111011100010010', '0111010000100010', '0100101110111000', '0111010000100001', '1011101110111000', '0001111010000111', '1011100011100001', '1101111010001011', '0001001000101110', '0111011101110111', '1011011100010001', '1110000100011110', '1011101110110111', '1110110110111011', '1000101110111000', '0100011100010001', '0010111010001000', '0111011101111000', '0010110110110111', '0010111010000100', '1000101110111011', '0100011100011110', '1000100010001011', '0111100010000111', '1110000101111000', '1101000100010010', '1000011100010001', '0010001001001011', '1110111010001000', '1101111011101110', '1101000101110111', '0100010001000100', '0010110110110100', '1000011100011101', '1011101111010010', '1011101111010001', '1000100011101101', '1110000101110100', '1011011101111000', '1101000101111011', '1011100010001000', '0111011100011101', '0100100010000111', '0010111011100010', '0111100010000100', '1110000101110111', '1110111011101110', '1110001000100001', '0001000101111000', '1101110111011110', '1110111011100010', '1000010000100010', '0010000100010010', '1000100011100010', '0100011101111000', '1011011100011110', '0010111010001011', '1011010001001011', '1011011100011101', '0001110110110100', '1110111010001011', '1011010001000100', '0100011100010010', '1000100010000100', '0111010001000111', '1101000100011110', '1011101111011110', '0010000101111011', '1000100010001000', '1110001000101101', '0001000100010001', '1000101110110111', '0010111011101110', '0100101111011101', '0001111010000100', '1000101111011110', '0100101110111011', '0001000100011110', '1110111010000111', '1110001001000111', '1110110110111000', '1101001000100010', '1011010000101110', '1110000100011101', '0111101110110100', '0100010000101110', '0001001001001011', '0111011101111011', '1000010001000100', '1110111010000100', '1110110111011101', '0111010000101101', '1101001001000100', '1000010000100001', '1011010001001000', '0010001000101101', '0001001000101101', '1000010000101101', '0001001001000100', '1000101111010010', '0100101111010001', '1110001000101110', '0100101111010010', '0100100010001011', '0111101110110111', '0111100011100010', '0010110111011110', '1000011101110111', '0001000101110100', '1000101111010001', '0111010001001011', '0100100011100010', '1000011100011110', '0100011101110100', '0001000100011101', '1101000100011101', '1101110110111000', '0001110111010001', '1110110111010001', '0010110110111011', '1011100011100010', '1000011101111011', '0100011100011101', '0100010000100010', '0001111011100010', '0010000100010001', '0001001000100001', '0001001001000111', '0100100011100001', '1000100011100001', '1110000100010010', '1011010000100001', '0111101111010010', '0100010000101101', '0111010001000100', '0111011101110100', '1000100010000111', '1101001000100001', '0100101111011110', '0010110111010010', '0010111011101101', '0001000101111011', '0100010001000111', '0100010001001000', '0111011100010001', '0001001000100010', '0010001001001000', '1011101111011101', '1101111011100001', '0111101111011101', '0100101110110100', '1101001001001000', '1110001001001011', '0111101110111000', '0010001001000111', '0010110111011101', '1011010001000111', '1000010001001000', '1110110111011110', '0001000101110111', '0010001000100001', '1011011101110111', '0100010000100001', '0010111011100001', '1110001001000100', '1011010000100010', '1110000100010001', '1101110111011101', '0111101111011110', '1011100010000100', '1011011101111011', '0100010001001011', '0111100010001000', '1011011100010010', '1101000100010001', '0100100011101101', '0111100011100001', '0001111011101101', '0111101111010001', '0001110110111000', '1000010001000111', '0001110111011110', '0010001001000100', '0100101110110111', '0010110110111000', '1000101110110100', '1000010001001011', '1101111010001000', '1011101110110100', '0001110111010010', '1101110110110100', '1000011101110100', '0001111011100001', '0001001001001000', '1000011100010010', '0001111010001011', '1011100011101101', '1101001000101110', '0100100010000100', '1011101110111011', '0100100010001000', '1011011101110100', '0010000100011101', '1110110110110100', '1011100010000111', '0010000101111000', '1110001001001000', '1110110110110111', '0111010001001000', '1011100010001011', '1101110110110111', '0001000100010010']
        hmm = 0
        for i,val in enumerate(bits):
            for j in code0:
                diff = 0
                for a,b in zip(val, j):
                    if a!=b:
                        diff += 1
                if diff == 0:
                    hmm += 1
                    break
                if diff == 1:
                    bits[i] = j
                    break
            for j in code1:
                diff = 0
                for a,b in zip(val, j):
                    if a!=b:
                        diff += 1
                if diff == 0:
                    hmm += 1
                    break
                if diff == 1:
                    bits[i] = j
                    break
        return bits

    def save(self):
        text_circuit = self.qc.draw(output="text")
        with open('circuit.txt', 'w') as f:
            f.write(str(text_circuit))

    def readout(self, pos: int, shots: int, p = 0):
        code0 = ['0110010100110000', '1111000001100110', '1111001101010101', '0011110011000011', '0011001100110000', '0011001100111111', '1001011001101001', '0011111110010101', '0101011001100110', '0000000000001111', '1010010100111100', '1010010101010101', '0011001101011001', '1100111110011010', '1100111110010101', '1100001101010110', '1111110010100110', '0101010100111100', '1010100110011010', '1010011001100101', '0011000001100101', '1111110011001100', '1100111110011001', '1100000001100101', '0000111111111111', '0011000000000011', '1001100110010101', '0110010101010110', '1010010100110011', '1010101011001100', '1001011001100101', '0000001101010110', '0101101011000000', '0011001101011010', '0000110010100110', '1001101010100101', '0000000000000011', '0011111111110000', '1010010101010110', '0110011001100110', '1111111111111100', '1100001100110011', '1010101011000000', '0000000001100101', '0110100110011001', '1001011001100110', '0011001100111100', '0101100110011001', '1100110010100101', '0011110010101001', '0101010100110011', '0110100111111100', '1100001101011001', '1100110010101001', '0000001100111111', '1010101010100101', '1001010100110000', '1100111111110011', '0110101011000000', '0101100110010101', '0101011001100101', '0000111110011010', '0000001100110000', '1010101010101001', '0000001101011010', '1001101010101001', '1111110011000011', '1001101011000000', '0101011001101010', '0101100110010110', '1111110010101010', '0110010100111100', '1100000000000000', '1001010101011010', '0110100110010101', '1010101011000011', '0101101011001100', '0110101011001111', '0110011000001100', '0110101010100110', '0110100110010110', '1100000000001111', '0101101010100101', '0110011000001111', '0000111110010110', '0011000000001111', '1010011000000000', '0110100111110011', '1100001101011010', '1001011001101010', '1010100111111100', '0011111111111111', '1010010101011001', '0000001101011001', '1111000000000000', '0101010101011001', '1100111110010110', '1100000000001100', '1001101010100110', '0011111110010110', '0000110010101010', '0101010100111111', '0011111110011010', '1010011000001100', '1001010100111100', '0000000001101010', '0000111111110011', '0110100111110000', '1001010101011001', '1100110010101010', '0110011001101010', '0101100111110011', '1010100111110011', '1001011000000011', '0101011000001100', '0000001100110011', '1001100111110011', '0000110011001100', '0011110011001111', '0011111110011001', '0101101010101010', '1111111110011001', '0011001100110011', '1001011000001111', '1010101010100110', '1100110011001100', '0101101010100110', '1010011000000011', '1111001100111100', '0101101011001111', '0101100111111100', '0101100110011010', '0110010100111111', '1010100110010110', '0000001100111100', '0110101010101001', '1001010101010110', '1010011001101001', '1001100111110000', '1001100110011010', '0000000000000000', '1111001101011010', '0101011000000000', '1111001101010110', '1100110011000011', '0101100111111111', '0000111111111100', '0110011000000000', '1010100111111111', '0011110011001100', '0110011001100101', '1010010100110000', '0101101011000011', '1100000001100110', '1010010101011010', '0110101011000011', '1100110010100110', '0110100111111111', '1111000001100101', '1010101011001111', '0000001101010101', '0110010101011010', '1111110010100101', '0000110011001111', '0110100110011010', '0110101010100101', '1111001100110011', '0101011000001111', '1111111110010101', '1001101011000011', '1100001101010101', '0011111111110011', '1010100111110000', '1001010101010101', '1001100110010110', '0101010100110000', '1111001100110000', '1111110010101001', '1001100111111111', '1001011000000000', '1001010100111111', '1010011001100110', '0011110010100101', '0110010100110011', '1010011001101010', '0000111110010101', '1100000001101001', '0101011000000011', '1100111111111111', '1111110011000000', '0000000001101001', '1001101011001100', '0110010101010101', '0101010101010101', '1111000001101010', '1111111111110011', '1111000000001111', '1010010100111111', '0011000001100110', '0101101010101001', '0011000000000000', '1111001100111111', '0110010101011001', '0000000000001100', '1001101011001111', '0011110010100110', '0101100111110000', '0011110010101010', '0000110011000011', '1010011000001111', '0101010101010110', '0101011001101001', '1111000001101001', '1001010100110011', '1100001100111100', '0011000000001100', '0011110011000000', '1100000001101010', '0011001101010110', '0000111110011001', '1100110011000000', '1100000000000011', '1100111111111100', '1111001101011001', '0011000001101010', '1100111111110000', '1001101010101010', '1111111110010110', '0011111111111100', '1111110011001111', '1010100110011001', '1111111111111111', '0110011000000011', '1111111111110000', '1100110011001111', '1111000000000011', '0110011001101001', '1100001100111111', '0011000001101001', '1001100111111100', '0000110010101001', '1001100110011001', '0101010101011010', '1100001100110000', '0011001101010101', '0000110011000000', '1001011000001100', '1010100110010101', '1010101010101010', '0110101010101010', '0000111111110000', '0110101011001100', '0000000001100110', '1111111110011010', '0000110010100101', '1111000000001100']
        code1 = ['0111100011101110', '0111100010001011', '0100100011101110', '0010000100011110', '1101111011100010', '0100011101110111', '0010110111010001', '0111100011101101', '1110000101111011', '1101110110111011', '0001111010001000', '0001110111011101', '1101001001000111', '0010111010000111', '1101000101110100', '1000010000101110', '0111010000101110', '0001110110111011', '1101111010000111', '1011010000101101', '0100011101111011', '1000100011101110', '0111101110111011', '1101001000101101', '1101000101111000', '1101110111010010', '0010000101110111', '0010001000100010', '1000101111011101', '0001111011101110', '0010000101110100', '1011100011101110', '1110110111010010', '0001110110110111', '1110111011101101', '1110111011100001', '0010001000101110', '1101111010000100', '0111011100011110', '1101001001001011', '1110001000100010', '1101111011101101', '1101110111010001', '1000011101111000', '0111011100010010', '0111010000100010', '0100101110111000', '0111010000100001', '1011101110111000', '0001111010000111', '1011100011100001', '1101111010001011', '0001001000101110', '0111011101110111', '1011011100010001', '1110000100011110', '1011101110110111', '1110110110111011', '1000101110111000', '0100011100010001', '0010111010001000', '0111011101111000', '0010110110110111', '0010111010000100', '1000101110111011', '0100011100011110', '1000100010001011', '0111100010000111', '1110000101111000', '1101000100010010', '1000011100010001', '0010001001001011', '1110111010001000', '1101111011101110', '1101000101110111', '0100010001000100', '0010110110110100', '1000011100011101', '1011101111010010', '1011101111010001', '1000100011101101', '1110000101110100', '1011011101111000', '1101000101111011', '1011100010001000', '0111011100011101', '0100100010000111', '0010111011100010', '0111100010000100', '1110000101110111', '1110111011101110', '1110001000100001', '0001000101111000', '1101110111011110', '1110111011100010', '1000010000100010', '0010000100010010', '1000100011100010', '0100011101111000', '1011011100011110', '0010111010001011', '1011010001001011', '1011011100011101', '0001110110110100', '1110111010001011', '1011010001000100', '0100011100010010', '1000100010000100', '0111010001000111', '1101000100011110', '1011101111011110', '0010000101111011', '1000100010001000', '1110001000101101', '0001000100010001', '1000101110110111', '0010111011101110', '0100101111011101', '0001111010000100', '1000101111011110', '0100101110111011', '0001000100011110', '1110111010000111', '1110001001000111', '1110110110111000', '1101001000100010', '1011010000101110', '1110000100011101', '0111101110110100', '0100010000101110', '0001001001001011', '0111011101111011', '1000010001000100', '1110111010000100', '1110110111011101', '0111010000101101', '1101001001000100', '1000010000100001', '1011010001001000', '0010001000101101', '0001001000101101', '1000010000101101', '0001001001000100', '1000101111010010', '0100101111010001', '1110001000101110', '0100101111010010', '0100100010001011', '0111101110110111', '0111100011100010', '0010110111011110', '1000011101110111', '0001000101110100', '1000101111010001', '0111010001001011', '0100100011100010', '1000011100011110', '0100011101110100', '0001000100011101', '1101000100011101', '1101110110111000', '0001110111010001', '1110110111010001', '0010110110111011', '1011100011100010', '1000011101111011', '0100011100011101', '0100010000100010', '0001111011100010', '0010000100010001', '0001001000100001', '0001001001000111', '0100100011100001', '1000100011100001', '1110000100010010', '1011010000100001', '0111101111010010', '0100010000101101', '0111010001000100', '0111011101110100', '1000100010000111', '1101001000100001', '0100101111011110', '0010110111010010', '0010111011101101', '0001000101111011', '0100010001000111', '0100010001001000', '0111011100010001', '0001001000100010', '0010001001001000', '1011101111011101', '1101111011100001', '0111101111011101', '0100101110110100', '1101001001001000', '1110001001001011', '0111101110111000', '0010001001000111', '0010110111011101', '1011010001000111', '1000010001001000', '1110110111011110', '0001000101110111', '0010001000100001', '1011011101110111', '0100010000100001', '0010111011100001', '1110001001000100', '1011010000100010', '1110000100010001', '1101110111011101', '0111101111011110', '1011100010000100', '1011011101111011', '0100010001001011', '0111100010001000', '1011011100010010', '1101000100010001', '0100100011101101', '0111100011100001', '0001111011101101', '0111101111010001', '0001110110111000', '1000010001000111', '0001110111011110', '0010001001000100', '0100101110110111', '0010110110111000', '1000101110110100', '1000010001001011', '1101111010001000', '1011101110110100', '0001110111010010', '1101110110110100', '1000011101110100', '0001111011100001', '0001001001001000', '1000011100010010', '0001111010001011', '1011100011101101', '1101001000101110', '0100100010000100', '1011101110111011', '0100100010001000', '1011011101110100', '0010000100011101', '1110110110110100', '1011100010000111', '0010000101111000', '1110001001001000', '1110110110110111', '0111010001001000', '1011100010001011', '1101110110110111', '0001000100010010']

        p_error = pauli_error([["X",p/2],["I",1-p],["Z",p/2]])
        p_error_2 = pauli_error([["XI",p/4],["IX",p/4],["II",1-p],["ZI",p/4],["IZ",p/4]])

        noise_model = NoiseModel()
        noise_model.add_all_qubit_quantum_error(p_error, ['x', "z", 'h', "s", "sdg", "id"])  # Apply to single-qubit gates
        noise_model.add_all_qubit_quantum_error(p_error_2, ['cx'])  # Apply to 2-qubit gates

        for i in range(16):
            self.qc.measure(i+16*pos, i)

        sim = AerSimulator()
        job = sim.run(self.qc, shots=shots, noise_model=noise_model)

        result = job.result()
        counts = result.get_counts()

        bitstring = list(counts.keys())

        bits = [i.replace(" ","") for i in bitstring]
        hmm = list(counts.values())

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