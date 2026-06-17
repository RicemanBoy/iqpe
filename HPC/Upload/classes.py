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
    with open("{}unitary{}.txt".format(path, n), "r") as file:
        for line in file:
            a.append(list(map(str, line.strip().split(","))))
    with open("{}adjunitary{}.txt".format(path, n), "r") as file:
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
                    self.qc = transpile(self.qc, optimization_level=1)
                    # print("Optimized: ")
                    gates(self.qc)
                    print("QEC counter: {}".format(self.qec_counter))
                    # if self.err:
                    #     self.qec(pos=0)
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
            print("Performance for angle {}: ".format(o), diff)
            bruh1.append(diff), y_list.append(diff)
    y = y/(n*k)
    arg = 0
    for i in range(len(bruh1)):
        arg += (y-bruh1[i])**2
    sigma = ((1/(k*n))*arg)**0.5
    sigma = sigma/((k*n)**0.5)

    return y, sigma, y_list


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
    def __init__(self, n: int, logical_q: int):                 #number of physical qubits per logical qubit
        self.ones = 0
        self.post = 0
        self.zeros = 0
        self.n = n          # number of physical qubits per logical qubit
        self.qec_counter = 0
        self.logicalq = logical_q
        self.err = False
        self.postselection = False                  #useless hier, PS geht nicht bei dem Code, da jeder Bitstring einem logischen Zustand entspricht

        qr = QuantumRegister(n*(logical_q+2)+1, "q")
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
        for i in range(2):
            for j in range(self.n):
                self.qc.reset(self.n*(i+self.logicalq)+j)
        for i in range(self.n):
            self.qc.h(self.n*(self.logicalq)+i)             #prep +_L
            self.qc.h(self.n*(self.logicalq+1)+i)
            self.qc.z(self.n*(self.logicalq+1)+i)           #prep -_L
        
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

    def rx(self, pos: int, angle: float):
        self.qc.rx(angle, self.n*pos)
    
    def sqrt_x(self, pos: int):
        # self.qc.h(self.n*pos)
        # self.qc.s(self.n*pos)
        # self.qc.h(self.n*pos)
        self.qc.rx(np.pi/2, self.n*pos)
    
    def sqrt_xdg(self, pos: int):
        # self.qc.h(self.n*pos)
        # self.qc.sdg(self.n*pos)
        # self.qc.h(self.n*pos)
        self.qc.rx(-np.pi/2, self.n*pos)
    
    def sqrt2_x(self, pos: int):
        # self.qc.h(self.n*pos)
        # self.qc.t(self.n*pos)
        # self.qc.h(self.n*pos)
        self.qc.rx(np.pi/4, self.n*pos)
    
    def sqrt2_xdg(self, pos: int):
        # self.qc.h(self.n*pos)
        # self.qc.tdg(self.n*pos)
        # self.qc.h(self.n*pos)
        self.qc.rx(-np.pi/4, self.n*pos)

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
            self.qec_ideal(pos=targ)
            self.qec_counter -= 1

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
                if self.err and self.qec_counter%1==0:
                    self.qec_ideal(pos = pos)
                    # print("QEC applied, counter:  {}".format(self.qec_counter))
            if i == "tdg":
                self.tdg(pos=pos)
                if self.err and self.qec_counter%1==0:
                    self.qec_ideal(pos = pos)
                    # print("QEC applied, counter:  {}".format(self.qec_counter))
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
            self.qc.h(anc)
            self.qc.cx(anc, self.n*pos + i)
            self.qc.cx(anc, self.n*pos + i + 1)
            self.qc.h(anc)
            self.qc.id(anc)
            self.qc.measure(anc, self.qecc[i])

        with self.qc.if_test((self.qecc[0], 1)):                #first
            with self.qc.if_test((self.qecc[1], 0)):               #second
                self.qc.z(self.n*pos)
        
        with self.qc.if_test((self.qecc[self.n-2], 1)):                 #last
            with self.qc.if_test((self.qecc[self.n-3], 0)):                #one before last
                self.qc.z(self.n*pos+self.n-1)

        for i in range(self.n-2):       
            with self.qc.if_test((self.qecc[i], 1)):
                with self.qc.if_test((self.qecc[i+1], 1)):
                    self.qc.z(self.n*pos + i + 1)
        
        self.qec_counter += 1

    def qec_ideal(self, pos: int):
        anc = self.qc.num_qubits - 1

        for i in range(self.n-1):
            self.qc.reset(anc)
            self.qc.append(h_ideal, [anc])
            self.qc.append(cx_ideal, [self.n*pos + i, anc])
            self.qc.append(cx_ideal, [self.n*pos + i + 1, anc])
            self.qc.append(h_ideal, [anc])
            self.qc.measure(anc, self.qecc[i])

        with self.qc.if_test((self.qecc[0], 1)):                #first
            with self.qc.if_test((self.qecc[1], 0)):               #second
                self.qc.append(z_ideal, [self.n*pos])
        
        with self.qc.if_test((self.qecc[self.n-2], 1)):                 #last
            with self.qc.if_test((self.qecc[self.n-3], 0)):                #one before last
                self.qc.append(z_ideal, [self.n*pos+self.n-1])

        for i in range(self.n-2):       
            with self.qc.if_test((self.qecc[i], 1)):
                with self.qc.if_test((self.qecc[i+1], 1)):
                    self.qc.append(z_ideal, [self.n*pos + i + 1])
        
        self.qec_counter += 1

    def readout(self, pos: int, shots: int, p: float, bias = 0):            #ICH MUSS READOUT MACHEN
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
