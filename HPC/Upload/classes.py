from pdb import pm
from turtle import pos

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
from qiskit.circuit.classical import expr

from itertools import product

matrix_h = ([[2**(-0.5),2**(-0.5)],[2**(-0.5),-2**(-0.5)]])
h_ideal = UnitaryGate(matrix_h)

matrix_cx = ([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
cx_ideal = UnitaryGate(matrix_cx)       #Erst Target, dann Control Qubit!!

matrix_x = ([[0,1],[1,0]])
x_ideal = UnitaryGate(matrix_x)

matrix_z = ([[1,0],[0,-1]])
z_ideal = UnitaryGate(matrix_z)


def majority_values(n):
    threshold = n // 2 + 1
    return [
        value
        for value in range(2**n)
        if value.bit_count() >= threshold
    ]

def gates(qc:QuantumCircuit):
    hmm = dict(qc.count_ops())
    hmm["reset"] = 0
    hmm["measure"] = 0
    hmm["if_else"] = 0
    hmm["id"] = 0
    hmm["swap"] = 0
    return print("Amount of gates in this circuit: ", sum(hmm.values()))
    #return print(hmm)

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
    assert code in ["steane", "rotsurf9", "rotsurf16", "repcode"], "Code not supported! Choose between 'steane', 'rotsurf9', 'rotsurf16' and 'repcode'."
    
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
                elif code == "repcode":
                    self = RepCode(3, 2)
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
                        self.tdg_switch(pos=0)
                self.h(pos=0)
                if self.err:
                    self.qec_ft(pos = 0)
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
    self = RotSurf9q(2)
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
            self.sdg_timo(pos=0)
        if l == 0.125:
            self.tdg_timo(pos=0)
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
                rots = [k*0.5 for k in rots]
                counter = 0                     #counter wie oft was wiederholt werden muss
                while True:
                    if code == "steane":
                        self = Steane7q(2, magic=1)
                    elif code == "rotsurf9":
                        self = RotSurf9q(2)
                    elif code == "rotsurf16":
                        self = RotSurf16q(2)
                    self.postselection = True
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
                            self.tdg(pos=0)
                    self.h(pos=0)
                    if self.err:
                        self.qec_ft(pos = 0)
                    self.qc = transpile(self.qc, optimization_level=1)
                    self.readout(pos=0, shots=1, p=noise)
                    if self.zeros == 1:
                        bitstring += "0"
                        break
                    if self.ones == 1:
                        bitstring += "1"
                        rots.append(0.5)
                        break
                    counter += 1
                    print("Angle {}, {}%% error, Iteration {}: {} Repetition".format(o, noise*100, t, counter))
            bitstring = bitstring[::-1]
            # print(bitstring)
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

def avg15_repcode(code: str, distance: int, iter: int, noise: float, qec = False, post = False, k = 1, bias = 0, path = ""):       #for bitflip protected rep code  
    assert code == "x" or code == "z", "Error: Only accept \"x\" or \"z\" as repetition codes!"
    n = 15
    angle = np.linspace(0,1,n+2)
    angle = np.delete(angle, [n+1])
    angle = np.delete(angle, [0])

    a, b = [], []
    with open("{}unitary{}_repz.txt".format(path, n), "r") as file:
        for line in file:
            a.append(list(map(str, line.strip().split(","))))
    with open("{}adjunitary{}_repz.txt".format(path, n), "r") as file:
        for line in file:
            b.append(list(map(str, line.strip().split(","))))
    
    y = 0
    y_list, bruh1 = [], []
    for m in range(k):
        for o in range(n):
            bitstring = ""
            rots = []
            for t in range(iter):
                rots = [k*0.5 for k in rots]
                counter = 0
                while True:
                    if code == "z":
                        self = RepCode_z(distance, 2)
                    elif code == "x":
                        self = RepCode(distance, 2)
                    self.err = qec
                    self.postselection = post
                    
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
                    # print("Unoptimized: ")
                    # gates(self.qc)
                    # self.qc = transpile(self.qc, optimization_level=1)
                    # print("Optimized: ")
                    gates(self.qc)
                    print("QEC counter: {}".format(self.qec_counter))
                    # if self.err:
                    #     self.qec_ideal(pos=0)
                    self.readout(pos=0, shots=1, p=noise, bias=bias)
            
                    if self.zeros == 1:
                        bitstring += "0"
                        break
                    if self.ones == 1:
                        bitstring += "1"
                        rots.append(0.5)
                        break
                    counter += 1
                    print("Angle {}, {}%% error, Iteration {}: {} Repetition".format(o, noise*100, t, counter))
                    del self
            bitstring = bitstring[::-1]
            hmm = convert(bitstring)
            diff = np.abs(hmm-angle[o])
            y += diff
            print("Performance {}for angle {}: ".format("(QEC) " if qec else "", o), diff)
            bruh1.append(diff), y_list.append(diff)
    y = y/(n*k)
    arg = 0
    for i in range(len(bruh1)):
        arg += (y-bruh1[i])**2
    sigma = ((1/(k*n))*arg)**0.5
    sigma = sigma/((k*n)**0.5)

    return y, sigma, y_list



class Steane7q:
    def __init__(self, n: int, magic = 1):
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
        for i in range(7):
            self.qc.sdg(i+7*pos)

        # self.qc.s(0+7*pos), self.qc.s(1+7*pos), self.qc.s(3+7*pos), self.qc.s(6+7*pos)
        # self.qc.sdg(2+7*pos), self.qc.sdg(4+7*pos), self.qc.sdg(5+7*pos)

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
        self.magiccounter += 1
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
        if self.err and self.magiccounter%2==0:
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

    def t_switch(self, pos: int):
        anc = self.qc.num_qubits - 1
        n = int((self.qc.num_qubits-2)/7) - 1
        self.qc.reset(anc)

        for i in range(3):
            self.qc.reset(i+7*n)
        
        for i in range(3):
            self.qc.h(i+7*n)

        self.qc.cx(5+7*pos, anc)
        self.qc.cx(6+7*pos, anc)
        self.qc.cx(0+7*n, anc)

        self.qc.measure(anc, self.qecc[0])
        self.qc.reset(anc)

        self.qc.cx(4+7*pos, anc)
        self.qc.cx(6+7*pos, anc)
        self.qc.cx(1+7*n, anc)

        self.qc.measure(anc, self.qecc[2])
        self.qc.reset(anc)

        self.qc.cx(2+7*pos, anc)
        self.qc.cx(6+7*pos, anc)
        self.qc.cx(2+7*n, anc)

        self.qc.measure(anc, self.qecc[1])

        #B_Z^4 = 0, B_Z^5 = 1, B_Z^6 = 2

        with self.qc.if_test((self.qecc[0],1)):
            with self.qc.if_test((self.qecc[1],0)):
                with self.qc.if_test((self.qecc[2],0)):
                    self.qc.x(0+7*pos)
                    self.qc.x(2+7*pos)
                    self.qc.x(4+7*pos)
                    self.qc.x(6+7*pos)
        
        with self.qc.if_test((self.qecc[0],0)):
            with self.qc.if_test((self.qecc[1],1)):
                with self.qc.if_test((self.qecc[2],0)):
                    self.qc.x(3+7*pos)
                    self.qc.x(4+7*pos)
                    self.qc.x(5+7*pos)
                    self.qc.x(6+7*pos)
        
        with self.qc.if_test((self.qecc[0],0)):
            with self.qc.if_test((self.qecc[1],0)):
                with self.qc.if_test((self.qecc[2],1)):
                    self.qc.x(1+7*pos)
                    self.qc.x(2+7*pos)
                    self.qc.x(5+7*pos)
                    self.qc.x(6+7*pos)
        
        with self.qc.if_test((self.qecc[0],1)):
            with self.qc.if_test((self.qecc[1],1)):
                with self.qc.if_test((self.qecc[2],0)):
                    self.qc.x(0+7*pos)
                    self.qc.x(2+7*pos)
                    self.qc.x(3+7*pos)
                    self.qc.x(5+7*pos)

        with self.qc.if_test((self.qecc[0],1)):
            with self.qc.if_test((self.qecc[1],0)):
                with self.qc.if_test((self.qecc[2],1)):
                    self.qc.x(0+7*pos)
                    self.qc.x(1+7*pos)
                    self.qc.x(4+7*pos)
                    self.qc.x(5+7*pos)
        
        with self.qc.if_test((self.qecc[0],0)):
            with self.qc.if_test((self.qecc[1],1)):
                with self.qc.if_test((self.qecc[2],1)):
                    self.qc.x(1+7*pos)
                    self.qc.x(2+7*pos)
                    self.qc.x(3+7*pos)
                    self.qc.x(4+7*pos)
        
        with self.qc.if_test((self.qecc[0],1)):
            with self.qc.if_test((self.qecc[1],1)):
                with self.qc.if_test((self.qecc[2],1)):
                    self.qc.x(0+7*pos)
                    self.qc.x(1+7*pos)
                    self.qc.x(3+7*pos)
                    self.qc.x(6+7*pos)
        
        self.qc.t(0+7*pos)              #The three edge qubits and the central one has T, the three qubits on the side have tdg, so in total 4 T, 3 Tdg, and one ccz on the bulk qubits
        self.qc.t(3+7*pos)
        self.qc.tdg(4+7*pos)
        self.qc.t(6+7*pos)
        self.qc.t(1+7*pos)
        self.qc.tdg(2+7*pos)
        self.qc.tdg(5+7*pos)

        self.qc.ccz(0+7*n, 1+7*n, 2+7*n)

        flags = ClassicalRegister(3)
        self.qc.add_register(flags)
        ancc = anc - 1

        self.qc.reset(anc), self.qc.reset(ancc)
        self.qc.h(anc)
        self.qc.cx(anc, 0+7*pos)
        self.qc.cx(anc, ancc)
        self.qc.cx(anc, 2+7*pos)
        self.qc.cx(anc, 4+7*pos)
        self.qc.cx(anc, ancc)
        self.qc.cx(anc, 6+7*pos)
        self.qc.h(anc)

        self.qc.measure(anc, self.qecc[0]), self.qc.measure(ancc, flags[0])
        self.qc.reset(anc), self.qc.reset(ancc)

        self.qc.h(anc)
        self.qc.cx(anc, 1+7*pos)
        self.qc.cx(anc, ancc)
        self.qc.cx(anc, 2+7*pos)
        self.qc.cx(anc, 5+7*pos)
        self.qc.cx(anc, ancc)
        self.qc.cx(anc, 6+7*pos)
        self.qc.h(anc)

        self.qc.measure(anc, self.qecc[2]), self.qc.measure(ancc, flags[1])
        self.qc.reset(anc), self.qc.reset(ancc)

        self.qc.h(anc)
        self.qc.cx(anc, 3+7*pos)
        self.qc.cx(anc, ancc)
        self.qc.cx(anc, 4+7*pos)
        self.qc.cx(anc, 5+7*pos)
        self.qc.cx(anc, ancc)
        self.qc.cx(anc, 6+7*pos)
        self.qc.h(anc)

        self.qc.measure(anc, self.qecc[1]), self.qc.measure(ancc, flags[2])

        #A_X^1 = 0 , A_X^2 = 1 , A_X^3 = 2

        with self.qc.if_test((self.qecc[0],1)):
            with self.qc.if_test((self.qecc[1],0)):
                with self.qc.if_test((self.qecc[2],0)):
                    self.qc.z(5+7*pos)
                    self.qc.z(6+7*pos)
                    self.qc.z(0+7*n)
        
        with self.qc.if_test((self.qecc[0],0)):
            with self.qc.if_test((self.qecc[1],1)):
                with self.qc.if_test((self.qecc[2],0)):
                    self.qc.z(2+7*pos)
                    self.qc.z(6+7*pos)
                    self.qc.z(2+7*n)
        
        with self.qc.if_test((self.qecc[0],0)):
            with self.qc.if_test((self.qecc[1],0)):
                with self.qc.if_test((self.qecc[2],1)):
                    self.qc.z(4+7*pos)
                    self.qc.z(6+7*pos)
                    self.qc.z(1+7*n)
        
        with self.qc.if_test((self.qecc[0],1)):
            with self.qc.if_test((self.qecc[1],1)):
                with self.qc.if_test((self.qecc[2],0)):
                    self.qc.z(2+7*pos)
                    self.qc.z(5+7*pos)
                    self.qc.z(0+7*n)
                    self.qc.z(1+7*n)
        
        with self.qc.if_test((self.qecc[0],1)):
            with self.qc.if_test((self.qecc[1],0)):
                with self.qc.if_test((self.qecc[2],1)):
                    self.qc.z(4+7*pos)
                    self.qc.z(5+7*pos)
                    self.qc.z(0+7*n)
                    self.qc.z(1+7*n)
        
        with self.qc.if_test((self.qecc[0],0)):
            with self.qc.if_test((self.qecc[1],1)):
                with self.qc.if_test((self.qecc[2],1)):
                    self.qc.z(2+7*pos)
                    self.qc.z(4+7*pos)
                    self.qc.z(1+7*n)
                    self.qc.z(2+7*n)
        
        with self.qc.if_test((self.qecc[0],1)):
            with self.qc.if_test((self.qecc[1],1)):
                with self.qc.if_test((self.qecc[2],1)):
                    self.qc.z(2+7*pos)
                    self.qc.z(4+7*pos)
                    self.qc.z(5+7*pos)
                    self.qc.z(6+7*pos)
                    self.qc.z(0+7*n)
                    self.qc.z(1+7*n)
                    self.qc.z(2+7*n)
    
    def tdg_switch(self, pos: int):
        anc = self.qc.num_qubits - 1
        n = int((self.qc.num_qubits-2)/7) - 1
        self.qc.reset(anc)

        for i in range(3):
            self.qc.reset(i+7*n)
        
        for i in range(3):
            self.qc.h(i+7*n)

        self.qc.cx(5+7*pos, anc)
        self.qc.cx(6+7*pos, anc)
        self.qc.cx(0+7*n, anc)

        self.qc.measure(anc, self.qecc[0])
        self.qc.reset(anc)

        self.qc.cx(4+7*pos, anc)
        self.qc.cx(6+7*pos, anc)
        self.qc.cx(1+7*n, anc)

        self.qc.measure(anc, self.qecc[2])
        self.qc.reset(anc)

        self.qc.cx(2+7*pos, anc)
        self.qc.cx(6+7*pos, anc)
        self.qc.cx(2+7*n, anc)

        self.qc.measure(anc, self.qecc[1])

        #B_Z^4 = 0, B_Z^5 = 1, B_Z^6 = 2

        with self.qc.if_test((self.qecc[0],1)):
            with self.qc.if_test((self.qecc[1],0)):
                with self.qc.if_test((self.qecc[2],0)):
                    self.qc.x(0+7*pos)
                    self.qc.x(2+7*pos)
                    self.qc.x(4+7*pos)
                    self.qc.x(6+7*pos)
        
        with self.qc.if_test((self.qecc[0],0)):
            with self.qc.if_test((self.qecc[1],1)):
                with self.qc.if_test((self.qecc[2],0)):
                    self.qc.x(3+7*pos)
                    self.qc.x(4+7*pos)
                    self.qc.x(5+7*pos)
                    self.qc.x(6+7*pos)
        
        with self.qc.if_test((self.qecc[0],0)):
            with self.qc.if_test((self.qecc[1],0)):
                with self.qc.if_test((self.qecc[2],1)):
                    self.qc.x(1+7*pos)
                    self.qc.x(2+7*pos)
                    self.qc.x(5+7*pos)
                    self.qc.x(6+7*pos)
        
        with self.qc.if_test((self.qecc[0],1)):
            with self.qc.if_test((self.qecc[1],1)):
                with self.qc.if_test((self.qecc[2],0)):
                    self.qc.x(0+7*pos)
                    self.qc.x(2+7*pos)
                    self.qc.x(3+7*pos)
                    self.qc.x(5+7*pos)

        with self.qc.if_test((self.qecc[0],1)):
            with self.qc.if_test((self.qecc[1],0)):
                with self.qc.if_test((self.qecc[2],1)):
                    self.qc.x(0+7*pos)
                    self.qc.x(1+7*pos)
                    self.qc.x(4+7*pos)
                    self.qc.x(5+7*pos)
        
        with self.qc.if_test((self.qecc[0],0)):
            with self.qc.if_test((self.qecc[1],1)):
                with self.qc.if_test((self.qecc[2],1)):
                    self.qc.x(1+7*pos)
                    self.qc.x(2+7*pos)
                    self.qc.x(3+7*pos)
                    self.qc.x(4+7*pos)
        
        with self.qc.if_test((self.qecc[0],1)):
            with self.qc.if_test((self.qecc[1],1)):
                with self.qc.if_test((self.qecc[2],1)):
                    self.qc.x(0+7*pos)
                    self.qc.x(1+7*pos)
                    self.qc.x(3+7*pos)
                    self.qc.x(6+7*pos)
        
        self.qc.tdg(0+7*pos)              #The three edge qubits and the central one has T, the three qubits on the side have tdg, so in total 4 T, 3 Tdg, and one ccz on the bulk qubits
        self.qc.tdg(3+7*pos)
        self.qc.t(4+7*pos)
        self.qc.tdg(6+7*pos)
        self.qc.tdg(1+7*pos)
        self.qc.t(2+7*pos)
        self.qc.t(5+7*pos)

        self.qc.ccz(0+7*n, 1+7*n, 2+7*n)

        flags = ClassicalRegister(3)
        self.qc.add_register(flags)
        ancc = anc - 1

        self.qc.reset(anc), self.qc.reset(ancc)
        self.qc.h(anc)
        self.qc.cx(anc, 0+7*pos)
        self.qc.cx(anc, ancc)
        self.qc.cx(anc, 2+7*pos)
        self.qc.cx(anc, 4+7*pos)
        self.qc.cx(anc, ancc)
        self.qc.cx(anc, 6+7*pos)
        self.qc.h(anc)

        self.qc.measure(anc, self.qecc[0]), self.qc.measure(ancc, flags[0])
        self.qc.reset(anc), self.qc.reset(ancc)

        self.qc.h(anc)
        self.qc.cx(anc, 1+7*pos)
        self.qc.cx(anc, ancc)
        self.qc.cx(anc, 2+7*pos)
        self.qc.cx(anc, 5+7*pos)
        self.qc.cx(anc, ancc)
        self.qc.cx(anc, 6+7*pos)
        self.qc.h(anc)

        self.qc.measure(anc, self.qecc[2]), self.qc.measure(ancc, flags[1])
        self.qc.reset(anc), self.qc.reset(ancc)

        self.qc.h(anc)
        self.qc.cx(anc, 3+7*pos)
        self.qc.cx(anc, ancc)
        self.qc.cx(anc, 4+7*pos)
        self.qc.cx(anc, 5+7*pos)
        self.qc.cx(anc, ancc)
        self.qc.cx(anc, 6+7*pos)
        self.qc.h(anc)

        self.qc.measure(anc, self.qecc[1]), self.qc.measure(ancc, flags[2])

        #A_X^1 = 0 , A_X^2 = 1 , A_X^3 = 2

        with self.qc.if_test((self.qecc[0],1)):
            with self.qc.if_test((self.qecc[1],0)):
                with self.qc.if_test((self.qecc[2],0)):
                    self.qc.z(5+7*pos)
                    self.qc.z(6+7*pos)
                    self.qc.z(0+7*n)
        
        with self.qc.if_test((self.qecc[0],0)):
            with self.qc.if_test((self.qecc[1],1)):
                with self.qc.if_test((self.qecc[2],0)):
                    self.qc.z(2+7*pos)
                    self.qc.z(6+7*pos)
                    self.qc.z(2+7*n)
        
        with self.qc.if_test((self.qecc[0],0)):
            with self.qc.if_test((self.qecc[1],0)):
                with self.qc.if_test((self.qecc[2],1)):
                    self.qc.z(4+7*pos)
                    self.qc.z(6+7*pos)
                    self.qc.z(1+7*n)
        
        with self.qc.if_test((self.qecc[0],1)):
            with self.qc.if_test((self.qecc[1],1)):
                with self.qc.if_test((self.qecc[2],0)):
                    self.qc.z(2+7*pos)
                    self.qc.z(5+7*pos)
                    self.qc.z(0+7*n)
                    self.qc.z(1+7*n)
        
        with self.qc.if_test((self.qecc[0],1)):
            with self.qc.if_test((self.qecc[1],0)):
                with self.qc.if_test((self.qecc[2],1)):
                    self.qc.z(4+7*pos)
                    self.qc.z(5+7*pos)
                    self.qc.z(0+7*n)
                    self.qc.z(1+7*n)
        
        with self.qc.if_test((self.qecc[0],0)):
            with self.qc.if_test((self.qecc[1],1)):
                with self.qc.if_test((self.qecc[2],1)):
                    self.qc.z(2+7*pos)
                    self.qc.z(4+7*pos)
                    self.qc.z(1+7*n)
                    self.qc.z(2+7*n)
        
        with self.qc.if_test((self.qecc[0],1)):
            with self.qc.if_test((self.qecc[1],1)):
                with self.qc.if_test((self.qecc[2],1)):
                    self.qc.z(2+7*pos)
                    self.qc.z(4+7*pos)
                    self.qc.z(5+7*pos)
                    self.qc.z(6+7*pos)
                    self.qc.z(0+7*n)
                    self.qc.z(1+7*n)
                    self.qc.z(2+7*n)

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
        self.magiccounter += 1
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
        if self.err and self.magiccounter%2==0:
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
                self.t(pos=pos)
                #self.t_cheat(pos=pos)
                # if self.err and self.magiccounter%2==0:
                #     self.qec_ft(pos = pos)
            if i == "tdg":
                self.tdg(pos=pos)
                #self.tdg_cheat(pos=pos)
                # if self.err and self.magiccounter%2==0:
                #     self.qec_ft(pos = pos)
            if i == "h":
                self.h(pos=pos)
            if i == "z":
                self.z(pos=pos)

    def cu(self, gate: list, adjgate: list):
        self.u2(0, gate=gate)
        # if self.err:
        #     self.qec(pos = 0)
        self.u2(1, gate=gate)
        if self.err:
            self.qec(pos = 0)
            self.qec(pos = 1)
        self.cnot(control=0, target=1)
        self.u2(1, gate=adjgate)
        # if self.err:
        #     self.qec(pos = 1)
        if self.err:
            self.qec(pos = 0)
            self.qec(pos = 1)
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
                    print("Postselected!")
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

class RepCode:      #Bitflip protected repetition code
    def __init__(self, n: int, logical_q: int):
        self.ones = 0
        self.zeros = 0
        self.post = 0
        self.n = n          # number of physical qubits per logical qubit
        self.qec_counter = 0
        self.postselection = False
        self.logicalq = logical_q
        self.err = False

        qr = QuantumRegister(n*(logical_q+2)+1, "q")
        cbit = ClassicalRegister(0, "c")
        self.qc = QuantumCircuit(qr, cbit)

        self.qecc = ClassicalRegister(n-1)
        self.qc.add_register(self.qecc)

        for i in range(n*logical_q):
            self.qc.id(i)

    def x(self, pos: int):
        for i in range(self.n):
            self.qc.x(self.n*pos + i)

    def z(self, pos: int):
        self.qc.z(self.n*pos)
    
    def h(self, pos: int):
        for i in range(self.n - 1):
            self.qc.cx(self.n*pos, self.n*pos + i + 1)
        self.qc.h(self.n*pos)
        for i in range(self.n - 1):
            self.qc.cx(self.n*pos, self.n*pos + i + 1)
    
    def h_toff(self, pos: int):
        for i in range(2):
            for j in range(self.n):
                self.qc.reset(self.n*(i+self.logicalq)+j)
        self.qc.h(pos=self.n*(self.logicalq))
        self.qc.cx(pos=self.n*(self.logicalq), control=self.n*(self.logicalq)+1)
        self.qc.cx(pos=self.n*(self.logicalq), control=self.n*(self.logicalq)+2)

        self.qc.h(pos=self.n*(self.logicalq+1))
        self.qc.cx(pos=self.n*(self.logicalq+1), control=self.n*(self.logicalq+1)+1)
        self.qc.cx(pos=self.n*(self.logicalq+1), control=self.n*(self.logicalq+1)+2)
        self.z(pos=self.n*(self.logicalq+1))
        
        self.toff(control1=self.logicalq, control2=pos, targ=self.logicalq+1)

        for i in range(self.n):                     #measure X_L
            self.qc.h(self.n*pos + i)
            self.qc.id(self.n*pos + i)
            self.qc.measure(self.n*pos + i, self.qecc[i])

        maj = majority_values(self.n)               #do majority vote to ensure FT
        for value in maj:
            with self.qc.if_test((self.qecc, value)):
                self.qc.x(self.n*self.logicalq)

        for i in range(self.n):                     #swap logical qubits such that the target qubit is at the same spot as before for convenience
            self.qc.swap(self.n*pos+i, self.n*self.logicalq+i)

    def s(self, pos: int):
        self.qc.s(self.n*pos)
    
    def sdg(self, pos: int):
        self.qc.sdg(self.n*pos)
    
    def t(self, pos: int):
        self.qc.t(self.n*pos)
    
    def tdg(self, pos: int):
        self.qc.tdg(self.n*pos)

    def toff(self, control1: int, control2: int, targ: int):
        for i in range(self.n):
            for j in range(self.n):
                self.qc.ccx(self.n*control1 + i, self.n*control2 + j, self.n*targ + j)
            self.qec(pos=targ)

    def rz(self, pos: int, angle: float):
        self.qc.rz(angle, self.n*pos)

    def cnot(self, control: int, target: int):
        for i in range(self.n):
            self.qc.cx(self.n*control + i, self.n*target + i)
    
    def u2(self, pos: int, gate: list):
        for i in gate:
            if i == "s":
                self.s(pos=pos)
            if i == "sdg":
                self.sdg(pos=pos)
            if i == "t":
                self.t(pos=pos)
                if self.err and self.qec_counter%2==0:
                    self.qec(pos = pos)
            if i == "tdg":
                self.tdg(pos=pos)
                #self.tdg_cheat(pos=pos)
                if self.err and self.qec_counter%2==0:
                    self.qec(pos = pos)
            if i == "h":
                self.h(pos=pos)
            if i == "z":
                self.z(pos=pos)

    def cu(self, gate: list, adjgate: list):
        self.u2(0, gate=gate)
        # if self.err:
        #     self.qec(pos = 0)
        self.u2(1, gate=gate)
        # if self.err:
        #     self.qec(pos = 1)
        self.cnot(control=0, target=1)
        self.u2(1, gate=adjgate)
        # if self.err:
        #     self.qec(pos = 1)
        self.cnot(control=0, target=1)

    def qec(self, pos: int):
        anc = self.qc.num_qubits - 1

        for i in range(self.n-1):
            self.qc.reset(anc)
            self.qc.cx(self.n*pos + i, anc)
            self.qc.cx(self.n*pos + i + 1, anc)
            self.qc.id(anc)
            self.qc.measure(anc, self.qecc[i])

        with self.qc.if_test((self.qecc[0], 1)):
            with self.qc.if_test((self.qecc[1], 0)):
                self.qc.x(self.n*pos)
        
        with self.qc.if_test((self.qecc[self.n-2], 1)):
            with self.qc.if_test((self.qecc[self.n-3], 0)):
                self.qc.x(self.n*pos+self.n-1)

        for i in range(self.n-1-1):       
            with self.qc.if_test((self.qecc[i], 1)):
                with self.qc.if_test((self.qecc[i+1], 1)):
                    self.qc.x(self.n*pos + i + 1)
        
        self.qec_counter += 1

    def readout(self, pos: int, shots: int, p: float, bias = 0):
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
        noise_model.add_all_qubit_quantum_error(p_error, ['x', "z", 'h', "s", "sdg", "t", "tdg", 'id'])  # Apply to single-qubit gates
        noise_model.add_all_qubit_quantum_error(p_error_2, ['cx'])  # Apply to 2-qubit gates

        read = ClassicalRegister(self.n)
        self.qc.add_register(read)

        for i in range(self.n):
            self.qc.id(self.n*pos + i)
            self.qc.measure(self.n*pos + i, read[self.n-1-i])

        test_0, test_1 = ["0"*self.n], ["1"*self.n]         #logical qubits

        sim = AerSimulator()
        job = sim.run(self.qc, shots=shots, noise_model=noise_model)
        result = job.result()
        counts = result.get_counts()

        bitstring = list(counts.keys())
        bits = [i.replace(" ","") for i in bitstring]
        counter = list(counts.values())

        bits = [i[:self.n] for i in bits]

        zero, one, post = 0, 0, 0

        if self.postselection:
            for i in range(len(bits)):
                for j in test_0:
                    if j == bits[i]:
                        bits[i] = 0
                        zero += counter[i]
                        break
                if bits[i] != 0:
                    for j in test_1:
                        if j == bits[i]:
                            bits[i] = 1
                            one += counter[i]
                            break
                if bits[i] != 1 and bits[i] != 0:
                    # print("Wrong bitstring: ", bits[i])
                    bits[i] = "post"
                    post += counter[i]
        else:
            for i, val in enumerate(bits):        #ohne Postselection, aka readout mit majority vote
                counter0, counter1 = 0, 0
                for j in val:
                    if j == "0":
                        counter0 += 1
                    else:
                        counter1 += 1
                if counter0 > counter1:
                    zero += counter[i]
                elif counter1 > counter0:
                    one += counter[i]
                else:
                    for k in range(counter[i]):
                        if np.random.rand() < 0.5:
                            zero += 1
                        else:
                            one += 1
        
        self.ones += one/shots
        self.zeros += zero/shots
        self.post += post/shots

class RepCode_z:      #Phaseflip protected repetition code
    def __init__(self, n: int, logical_q: int):              
        self.ones = 0
        self.post = 0
        self.zeros = 0
        self.n = n          # number of physical qubits per logical qubit
        self.qec_counter = 0
        self.magiccounter = 0
        self.logicalq = logical_q
        self.err = False
        self.postselection = False                  #useless hier, PS geht nicht bei dem Code, da jeder Bitstring einem logischen Zustand entspricht

        qr = QuantumRegister(n*(logical_q+2)+2, "q")
        cbit = ClassicalRegister(0, "c")
        self.qc = QuantumCircuit(qr, cbit)

        self.qecc = ClassicalRegister(n)
        self.qc.add_register(self.qecc)

        for i in range(n*logical_q):
            self.qc.id(i)
            self.qc.h(i)
        
        for i in range(logical_q):
            self.h(pos=i)

    def z(self, pos: int):
        for i in range(self.n):
            self.qc.z(self.n*pos + i)

    def x(self, pos: int):
        self.qc.x(self.n*pos)
    
    def h_nft(self, pos: int):
        for i in range(self.n - 1):
            self.qc.cx(self.n*pos + i + 1, self.n*pos)
        self.qc.h(self.n*pos)
        for i in range(self.n - 1):
            self.qc.cx(self.n*pos + i + 1, self.n*pos)
    
    def h(self, pos: int):
        self.magiccounter += 1
        for i in range(2):
            for j in range(self.n):
                self.qc.reset(self.n*(i+self.logicalq)+j)
        for i in range(self.n):
            self.qc.h(self.n*(self.logicalq)+i)             #prep +_L
            self.qc.h(self.n*(self.logicalq+1)+i)
            self.qc.z(self.n*(self.logicalq+1)+i)           #prep -_L
        # if self.err:
        #     self.qec(pos = self.logicalq+1)   #qec bei -_L
        
        self.toff(control1=self.logicalq, control2=pos, targ=self.logicalq+1)           #computationally very though

        for i in range(self.n):                     #measure X_L
            self.qc.h(self.n*pos + i)
            self.qc.id(self.n*pos + i)
            self.qc.measure(self.n*pos + i, self.qecc[i])

        maj = majority_values(self.n)               #do majority vote to ensure FT, somewhat of an QEC step in itself
        for value in maj:
            with self.qc.if_test((self.qecc, value)):
                self.qc.x(self.n*self.logicalq)

        for i in range(self.n):                     #swap logical qubits such that the target qubit is at the same spot as before for convenience
            self.qc.swap(self.n*pos+i, self.n*self.logicalq+i)

    def rx(self, pos: int, angle: float):
        i = np.random.randint(0,2)
        self.qc.rx(angle, self.n*pos+i)
    
    def sqrt_x(self, pos: int):
        # self.qc.h(self.n*pos)
        # self.qc.s(self.n*pos)
        # self.qc.h(self.n*pos)
        i = np.random.randint(0,2)
        self.qc.rx(np.pi/2, self.n*pos+i)
    
    def sqrt_xdg(self, pos: int):
        # self.qc.h(self.n*pos)
        # self.qc.sdg(self.n*pos)
        # self.qc.h(self.n*pos)
        i = np.random.randint(0,2)
        self.qc.rx(-np.pi/2, self.n*pos+i)
    
    def sqrt2_x(self, pos: int):
        # self.qc.h(self.n*pos)
        # self.qc.t(self.n*pos)
        # self.qc.h(self.n*pos)
        i = np.random.randint(0,2)
        self.qc.rx(np.pi/4, self.n*pos+i)
    
    def sqrt2_xdg(self, pos: int):
        # self.qc.h(self.n*pos)
        # self.qc.tdg(self.n*pos)
        # self.qc.h(self.n*pos)
        i = np.random.randint(0,2)
        self.qc.rx(-np.pi/4, self.n*pos+i)

    def s(self, pos: int):
        self.h(pos=pos)
        self.sqrt_x(pos=pos)
        self.h(pos=pos)

    def sdg(self, pos: int):
        self.h(pos=pos)
        self.sqrt_xdg(pos=pos)
        self.h(pos=pos)
    
    def t(self, pos: int):
        self.h(pos=pos)
        self.sqrt2_x(pos=pos)
        self.h(pos=pos)

    def tdg(self, pos: int):
        self.h(pos=pos)
        self.sqrt2_xdg(pos=pos)
        self.h(pos=pos)

    def toff(self, control1: int, control2: int, targ: int):
        for i in range(self.n):
            for j in range(self.n):
                self.qc.ccx(self.n*control1 + i, self.n*control2 + j, self.n*targ + j)
            if self.err:
                self.qec(pos=targ)               #needed for FT
                # self.qec_counter -= 1

    def cnot(self, control: int, target: int):
        for i in range(self.n):
            self.qc.cx(self.n*control + i, self.n*target + i)
    
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
            if i == "t_x":
                self.sqrt2_x(pos=pos)
            if i == "tdg_x":
                self.sqrt2_xdg(pos=pos)
            if i == "s_x":
                self.sqrt_x(pos=pos)
            if i == "sdg_x":
                self.sqrt_xdg(pos=pos)

    def cu(self, gate: list, adjgate: list):
        self.u2(0, gate=gate)
        # if self.err:
        #     self.qec(pos = 0)
        self.u2(1, gate=gate)
        # if self.err:
        #     self.qec(pos = 1)
        self.cnot(control=0, target=1)
        self.u2(1, gate=adjgate)
        # if self.err:
        #     self.qec(pos = 1)
        self.cnot(control=0, target=1)

    def qec_nft(self, pos: int):
        anc = self.qc.num_qubits - 1
        self.qec_counter += 1

        # for i in range(self.n-1):
        #     self.qc.reset(anc)
        #     self.qc.h(anc)
        #     self.qc.cx(anc, self.n*pos + i)
        #     self.qc.cx(anc, self.n*pos + i + 1)
        #     self.qc.h(anc)
        #     self.qc.id(anc)
        #     self.qc.measure(anc, self.qecc[i])

        # with self.qc.if_test((self.qecc[0], 1)):                #first
        #     with self.qc.if_test((self.qecc[1], 0)):               #second
        #         self.qc.z(self.n*pos)
        
        # with self.qc.if_test((self.qecc[self.n-2], 1)):                 #last
        #     with self.qc.if_test((self.qecc[self.n-3], 0)):                #one before last
        #         self.qc.z(self.n*pos+self.n-1)

        # for i in range(self.n-2):       
        #     with self.qc.if_test((self.qecc[i], 1)):
        #         with self.qc.if_test((self.qecc[i+1], 1)):
        #             self.qc.z(self.n*pos + i + 1)

        self.qc.reset(anc)
        self.qc.h(anc)
        self.qc.cx(anc, 3*pos + 0)
        self.qc.cx(anc, 3*pos + 1)
        self.qc.h(anc)
        self.qc.id(anc)
        self.qc.measure(anc, self.qecc[0])

        self.qc.reset(anc)
        self.qc.h(anc)
        self.qc.cx(anc, 3*pos + 1)
        self.qc.cx(anc, 3*pos + 2)
        self.qc.h(anc)
        self.qc.id(anc)
        self.qc.measure(anc, self.qecc[1])

        with self.qc.if_test((self.qecc[0], 1)):                
            with self.qc.if_test((self.qecc[1], 0)):               
                self.qc.z(3*pos)

        with self.qc.if_test((self.qecc[0], 1)):                
            with self.qc.if_test((self.qecc[1], 1)):               
                self.qc.z(3*pos + 1)
        
        with self.qc.if_test((self.qecc[0], 0)):                
            with self.qc.if_test((self.qecc[1], 1)):               
                self.qc.z(3*pos + 2)

    def qec_ft(self, pos: int):
        anc = self.qc.num_qubits - 1
        ancc = anc - 1
        self.qec_counter += 1

        # for i in range(self.n-1):
        #     self.qc.reset(anc)
        #     self.qc.h(anc)
        #     self.qc.cx(anc, self.n*pos + i)
        #     self.qc.cx(anc, self.n*pos + i + 1)
        #     self.qc.h(anc)
        #     self.qc.id(anc)
        #     self.qc.measure(anc, self.qecc[i])

        # with self.qc.if_test((self.qecc[0], 1)):                #first
        #     with self.qc.if_test((self.qecc[1], 0)):               #second
        #         self.qc.z(self.n*pos)
        
        # with self.qc.if_test((self.qecc[self.n-2], 1)):                 #last
        #     with self.qc.if_test((self.qecc[self.n-3], 0)):                #one before last
        #         self.qc.z(self.n*pos+self.n-1)

        # for i in range(self.n-2):       
        #     with self.qc.if_test((self.qecc[i], 1)):
        #         with self.qc.if_test((self.qecc[i+1], 1)):
        #             self.qc.z(self.n*pos + i + 1)

        self.qc.reset(anc), self.qc.reset(ancc)
        self.qc.h(anc), self.qc.h(ancc)
        self.qc.cx(anc, 3*pos + 0)
        self.qc.cx(ancc, 3*pos + 1)
        self.qc.cx(anc, ancc)
        self.qc.h(anc)
        self.qc.id(anc)
        self.qc.measure(anc, self.qecc[0])

        self.qc.reset(anc), self.qc.reset(ancc)
        self.qc.h(anc), self.qc.h(ancc)
        self.qc.cx(anc, 3*pos + 1)
        self.qc.cx(ancc, 3*pos + 2)
        self.qc.cx(anc, ancc)
        self.qc.h(anc)
        self.qc.id(anc)
        self.qc.measure(anc, self.qecc[1])

        with self.qc.if_test((self.qecc[0], 1)):                
            with self.qc.if_test((self.qecc[1], 0)):               
                self.qc.z(3*pos)

        with self.qc.if_test((self.qecc[0], 1)):                
            with self.qc.if_test((self.qecc[1], 1)):               
                self.qc.z(3*pos + 1)
        
        with self.qc.if_test((self.qecc[0], 0)):                
            with self.qc.if_test((self.qecc[1], 1)):               
                self.qc.z(3*pos + 2)

    def qec(self, pos: int):
        anc = self.qc.num_qubits - 1
        self.qec_counter += 1

        # for i in range(self.n-1):           #general single qubit correction for arbitrary n
        #     self.qc.reset(anc)
        #     self.qc.append(h_ideal, [anc])
        #     self.qc.append(cx_ideal, [self.n*pos + i, anc])
        #     self.qc.append(cx_ideal, [self.n*pos + i + 1, anc])
        #     self.qc.append(h_ideal, [anc])
        #     self.qc.measure(anc, self.qecc[i])

        # with self.qc.if_test((self.qecc[0], 1)):                #first
        #     with self.qc.if_test((self.qecc[1], 0)):               #second
        #         self.qc.append(z_ideal, [self.n*pos])
        
        # with self.qc.if_test((self.qecc[self.n-2], 1)):                 #last
        #     with self.qc.if_test((self.qecc[self.n-3], 0)):                #one before last
        #         self.qc.append(z_ideal, [self.n*pos+self.n-1])

        # for i in range(self.n-2):       
        #     with self.qc.if_test((self.qecc[i], 1)):
        #         with self.qc.if_test((self.qecc[i+1], 1)):
        #             self.qc.append(z_ideal, [self.n*pos + i + 1])
        

        #specific case for n = 3
        self.qc.reset(anc)      
        self.qc.append(h_ideal, [anc])
        self.qc.append(cx_ideal, [3*pos+0, anc])
        self.qc.append(cx_ideal, [3*pos+1, anc])
        self.qc.append(h_ideal, [anc])
        self.qc.measure(anc, self.qecc[0])

        self.qc.reset(anc)      
        self.qc.append(h_ideal, [anc])
        self.qc.append(cx_ideal, [3*pos+1, anc])
        self.qc.append(cx_ideal, [3*pos+2, anc])
        self.qc.append(h_ideal, [anc])
        self.qc.measure(anc, self.qecc[1])

        with self.qc.if_test((self.qecc[0], 1)):                #first
            with self.qc.if_test((self.qecc[1], 0)):               #second
                self.qc.append(z_ideal, [3*pos+0])
        
        with self.qc.if_test((self.qecc[0], 1)):                #first
            with self.qc.if_test((self.qecc[1], 1)):               #second
                self.qc.append(z_ideal, [3*pos+1])
        
        with self.qc.if_test((self.qecc[0], 0)):                #first
            with self.qc.if_test((self.qecc[1], 1)):               #second
                self.qc.append(z_ideal, [3*pos+2])

    def readout(self, pos: int, shots: int, p: float, bias = 0):        
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
        noise_model.add_all_qubit_quantum_error(p_error, ['x', "z", 'h', "s", "sdg", "t", "tdg", 'id'])  # Apply to single-qubit gates
        noise_model.add_all_qubit_quantum_error(p_error_2, ['cx'])  # Apply to 2-qubit gates
        noise_model.add_all_qubit_quantum_error(p_error_3, ['ccx'])  # Apply to 3-qubit gates

        count0, count1 = [], []                 #alle statevectors für 0_L und 1_L
        for i in range(2**self.n):
            bit = inv_covert(i,self.n)
            if bit.count("1")%2 == 0:
                count0.append(bit)
            else:
                count1.append(bit)

        read = ClassicalRegister(self.n)
        self.qc.add_register(read)

        for i in range(self.n):
            self.qc.id(self.n*pos + i)
            self.qc.measure(self.n*pos + i, read[self.n-1-i])

        sim = AerSimulator()
        job = sim.run(self.qc, shots=shots, noise_model=noise_model)
        result = job.result()
        counts = result.get_counts()

        bitstring = list(counts.keys())
        bits = [i.replace(" ","") for i in bitstring]
        counter = list(counts.values())

        allcbits = len(bitstring[0])

        bits = [i[:self.n] for i in bits]

        for i in range(len(bits)):
            for j in count0:
                if j == bits[i]:
                    bits[i] = 0
                    break
            if bits[i] != 0:
                for j in count1:
                    if j == bits[i]:
                        bits[i] = 1
                        break
            if bits[i] != 1 and bits[i] != 0:
                # print("Wrong bitstring: ", bits[i])
                if self.postselection:
                    bits[i] = "post"
                    print("Postselected!")
                else:
                    if np.random.rand() < 0.5:
                        bits[i] = 0
                    else:
                        bits[i] = 1

        zero, one, err = 0, 0, 0

        for i, val in enumerate(bits):
            if val == 0:
                zero += counter[i]
            elif val == 1:
                one += counter[i]
            else:
                err += counter[i]

        self.ones += one/shots
        self.zeros += zero/shots
        self.post += err/shots
