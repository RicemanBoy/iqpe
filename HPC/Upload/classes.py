from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister, qasm3, qasm2, qpy
from qiskit.visualization import plot_histogram
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from functools import wraps
from dataclasses import dataclass
#import bitstring
from qiskit_aer import AerSimulator
from qiskit.transpiler.passes.synthesis import SolovayKitaev
from qiskit.synthesis import generate_basic_approximations
from qiskit.quantum_info import Operator

from qiskit import transpile
from qiskit.quantum_info import Statevector
from qiskit.primitives import StatevectorSampler, StatevectorEstimator
from qiskit.quantum_info import SparsePauliOp

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

@dataclass
class Command:
    method: str
    args: tuple
    kwargs: dict

#Damit kann ich speichern welche Gates im QC aufgerufen wurde für Statevectorsimulation
def record(func):
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        if self._recording and self._call_depth == 0:
            self.history.append(Command(func.__name__, args, kwargs))

        self._call_depth += 1
        try:
            return func(self, *args, **kwargs)
        finally:
            self._call_depth -= 1

    return wrapper

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
                    #gates(self.qc)
                    #print("QEC counter: {}".format(self.qec_counter))
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
            #diff = np.abs(hmm-angle[o])
            diff = min(np.abs(hmm-angle[o]), 1-np.abs(hmm-angle[o]))
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

def avg7_repcode(code: str, distance: int, iter: int, noise: float, qec = False, post = False, k = 1, bias = 0, path = ""):       #only exact angles!  
    assert code == "x" or code == "z", "Error: Only accept \"x\" or \"z\" as repetition codes!"
    n = 15
    angle = np.linspace(0,1,n+2)
    angle = np.delete(angle, [n+1])
    angle = np.delete(angle, [0])

    a, b = [], []
    if code == "z":
        with open("{}unitary{}_repz.txt".format(path, n), "r") as file:
            for line in file:
                a.append(list(map(str, line.strip().split(","))))
        with open("{}adjunitary{}_repz.txt".format(path, n), "r") as file:
            for line in file:
                b.append(list(map(str, line.strip().split(","))))
    else:
        with open("{}unitary{}.txt".format(path, n), "r") as file:
            for line in file:
                a.append(list(map(str, line.strip().split(","))))
        with open("{}adjunitary{}.txt".format(path, n), "r") as file:
            for line in file:
                b.append(list(map(str, line.strip().split(","))))
    
    y = 0
    y_list, bruh1 = [], []
    for m in range(k):
        for o in range(7):
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
                        self.cu(a[2*o+1], b[2*o+1])
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
                    #gates(self.qc)
                    #print("QEC counter: {}".format(self.qec_counter))
                    # if self.err:
                    #     self.qec_ideal(pos=0)
                    self.readout(pos=0, shots=1, p=noise, bias=bias)
                    print("T/Tdg counter for iteration {}: {}".format(t, self.magiccounter))
            
                    if self.zeros == 1:
                        bitstring += "0"
                        break
                    if self.ones == 1:
                        bitstring += "1"
                        rots.append(0.5)
                        break
                    counter += 1
                    print("Angle {}, {}%% error, Iteration {}: {} Repetition".format(2*o+1, noise*100, t, counter))
                    del self
            bitstring = bitstring[::-1]
            hmm = convert(bitstring)
            diff = min(np.abs(hmm-angle[2*o+1]), 1-np.abs(hmm-angle[2*o+1]))
            y += diff
            print("Performance {}for angle {}: ".format("(QEC) " if qec else "", 2*o+1), diff)
            bruh1.append(diff), y_list.append(diff)
    y = y/(n*k)
    arg = 0
    for i in range(len(bruh1)):
        arg += (y-bruh1[i])**2
    sigma = ((1/(k*n))*arg)**0.5
    sigma = sigma/((k*n)**0.5)

    return y, sigma, y_list

def avg7_repcode_timo(code: str, distance: int, iter: int, noise: float, qec = False, post = False, k = 1, bias = 0, path = ""):       #only exact angles and take avg of N shots of each angle with pos/neg performance possible  
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
        for o in range(7):
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
                        self.cu(a[2*o+1], b[2*o+1])
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
                    #gates(self.qc)
                    #print("QEC counter: {}".format(self.qec_counter))
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
                    print("Angle {}, {}%% error, Iteration {}: {} Repetition".format(2*o+1, noise*100, t, counter))
                    del self
            bitstring = bitstring[::-1]
            hmm = convert(bitstring)
            diff = 0
            if np.abs(hmm-angle[2*o+1]) <= 1-np.abs(hmm-angle[2*o+1]):
                diff += hmm-angle[2*o+1]
            else:
                diff += 1-np.abs(hmm-angle[2*o+1])
            y += diff
            print("Performance {}for angle {}: ".format("(QEC) " if qec else "", 2*o+1), diff)
            bruh1.append(diff), y_list.append(hmm)              #y_list enthält die gemessenen Winkel, nicht die Performance!
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
        noise_model.add_all_qubit_quantum_error(p_error, ['x', "z", 'h', "s", "sdg", "id", "t", "tdg"])  # Apply to single-qubit gates
        noise_model.add_all_qubit_quantum_error(p_error_2, ['cx'])  # Apply to 2-qubit gates

        read = ClassicalRegister(7)
        self.qc.add_register(read)

        for i in range(7):
            # self.qc.id(i+7*pos)
            self.qc.measure(i+7*pos,read[6-i])

        # self.qc = transpile(self.qc, optimization_level=1)

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

        print(bitstring)
        print(bits)
        print("FLAGS: ", postprocess)

        for i in range(len(pre)):
            if pre[i].count("1") != 0:
                bits[i] = "pre"
                # print("AHAAAA wieso flag???")

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
        self.postselection = True
        self.classical_ec = False
        self.qec_counter = 0

        self.hadamards = [0 for i in range(n+magic)]

        qr = QuantumRegister(9*(n+magic)+4, "q")
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
    
    def cnot_test(self, control: int, target: int):
        if self.hadamards[0]%2==1 and self.hadamards[1]%2==0:
            self.qc.h(target)
            self.cz_test()
            self.qc.h(target)
        elif self.hadamards[0]%2==0 and self.hadamards[1]%2==1:
            self.qc.h(target)
            self.cz_test()
            self.qc.h(target)
        else:
            for i in range(9):
                self.qc.cx(9*control+i,9*target+i)

    def cz(self):
        self.h(pos=1)
        self.cnot(control = 0, target=1)
        self.h(pos=1)

    def cz_test(self):
        if self.hadamards[0]%2==1 and self.hadamards[1]%2==0:
            self.qc.cz(0,9+6)
            self.qc.cz(1,9+3)
            self.qc.cz(2,9+0)
            self.qc.cz(3,9+7)
            self.qc.cz(4,9+4)
            self.qc.cz(5,9+1)
            self.qc.cz(6,9+8)
            self.qc.cz(7,9+5)
            self.qc.cz(8,9+2)
        elif self.hadamards[0]%2==0 and self.hadamards[1]%2==1:
            self.qc.cz(0,9+2)
            self.qc.cz(1,9+5)
            self.qc.cz(2,9+8)
            self.qc.cz(3,9+1)
            self.qc.cz(4,9+4)
            self.qc.cz(5,9+7)
            self.qc.cz(6,9+0)
            self.qc.cz(7,9+3)
            self.qc.cz(8,9+6)
        else:
            print("Same orientation!")
            for i in range(9):
                self.qc.cz(i,9+i)

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
                self.s_timo(pos=pos)
            if i == "sdg":
                #  if self.err == True and self.qec_counter%5==0:
                #     self.qec_zstab(pos=pos)
                 self.sdg_timo(pos=pos)
            if i == "t":
                 if self.err == True and self.qec_counter%8==0:
                    self.qec(pos=pos)
                 self.t_timo(pos=pos)
            if i == "tdg":
                 if self.err == True and self.qec_counter%8==0:
                    self.qec(pos=pos)
                 self.tdg_timo(pos=pos)
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
        flags = [i[9:len(bitstring[0])-4] for i in bitstring]

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
        self.magiccounter = 0

        qr = QuantumRegister(n*(logical_q+2)+3, "q")
        cbit = ClassicalRegister(0, "c")
        self.qc = QuantumCircuit(qr, cbit)

        self.qecc = ClassicalRegister(n)
        self.qc.add_register(self.qecc)

        # for i in range(n*logical_q):
        #     self.qc.id(i)

    def x(self, pos: int):
        for i in range(self.n):
            self.qc.x(self.n*pos + i)

    def z(self, pos: int):
        self.qc.z(self.n*pos)
    
    def h_nft(self, pos: int):
        for i in range(self.n - 1):
            self.qc.cx(self.n*pos, self.n*pos + i + 1)
        self.qc.h(self.n*pos)
        for i in range(self.n - 1):
            self.qc.cx(self.n*pos, self.n*pos + i + 1)
    
    def h(self, pos: int):
        self.magiccounter += 1
        for i in range(2):
            for j in range(self.n):
                self.qc.reset(self.n*(i+self.logicalq)+j)

        self.qc.h(self.n*(self.logicalq)+0)             #prep +_L
        self.qc.cx(self.n*(self.logicalq)+0, self.n*(self.logicalq)+1)
        self.qc.cx(self.n*(self.logicalq)+0, self.n*(self.logicalq)+2)

        self.qc.h(self.n*(self.logicalq+1)+0)           #prep -_L
        self.qc.cx(self.n*(self.logicalq+1)+0, self.n*(self.logicalq+1)+1)
        self.qc.cx(self.n*(self.logicalq+1)+0, self.n*(self.logicalq+1)+2)
        self.qc.z(self.n*(self.logicalq+1)+0)
        # if self.err:
        #     self.qec(pos = self.logicalq+1)   #qec bei -_L
        
        self.toff(control1=self.logicalq, control2=pos, targ=self.logicalq+1)           #computationally very though

        for i in range(self.n):                     #measure X_L
            self.qc.h(self.n*pos + i)
            # self.qc.id(self.n*pos + i)
            self.qc.measure(self.n*pos + i, self.qecc[i])

        maj = majority_values(self.n)               #do majority vote to ensure FT, somewhat of an QEC step in itself
        for i in range(self.n):
            for value in maj:
                with self.qc.if_test((self.qecc, value)):
                    self.qc.x(self.n*self.logicalq+i)

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
            if self.err:
                self.qec(pos=targ)               #needed for FT
                # self.qec_counter -= 1

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
        self.zeros = 0
        self.n = n          # number of physical qubits per logical qubit
        self.qec_counter = 0
        self.magiccounter = 0
        self.logicalq = logical_q
        self.err = False
        self.noise_model = self.__noise_model__(0,0)

        ### Speichert welche Gates benutzt wurden, wichtig für Statevector simulation ####
        self.history = []
        self._call_depth = 0
        self._recording = False
        ##################################################################################

        qr = QuantumRegister(n*(logical_q+2)+3, "q")
        self.qc = QuantumCircuit(qr)

        self.qecc = ClassicalRegister(n)
        self.qc.add_register(self.qecc)

        for i in range(n*logical_q):
            # self.qc.id(i)
            self.qc.h(i)
        
        for i in range(logical_q):
            self.h(i)

        self._recording = True

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
   
    @record
    def z(self, pos: int):
        for i in range(self.n):
            self.qc.z(self.n*pos + i)

    @record
    def id(self, pos: int):
        for i in range(self.n):
            self.qc.id(self.n*pos + i)

    @record
    def x(self, pos: int):
        self.qc.x(self.n*pos)
    
    def h_nft(self, pos: int):
        for i in range(self.n - 1):
            self.qc.cx(self.n*pos + i + 1, self.n*pos)
        self.qc.h(self.n*pos)
        for i in range(self.n - 1):
            self.qc.cx(self.n*pos + i + 1, self.n*pos)

    @record
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
            # self.qc.id(self.n*pos + i)
            self.qc.measure(self.n*pos + i, self.qecc[i])

        maj = majority_values(self.n)               #do majority vote to ensure FT, somewhat of an QEC step in itself
        for value in maj:
            with self.qc.if_test((self.qecc, value)):
                self.qc.x(self.n*self.logicalq)

        for i in range(self.n):                     #swap logical qubits such that the target qubit is at the same spot as before for convenience
            self.qc.swap(self.n*pos+i, self.n*self.logicalq+i)

    @record
    def rx(self, pos: int, angle: float):
        i = np.random.randint(0,2)
        self.qc.rx(angle, self.n*pos+i)

    @record
    def sqrt_x(self, pos: int):
        # self.qc.h(self.n*pos)
        # self.qc.s(self.n*pos)
        # self.qc.h(self.n*pos)
        i = np.random.randint(0,2)
        self.qc.rx(np.pi/2, self.n*pos+i)

    @record
    def sqrt_xdg(self, pos: int):
        # self.qc.h(self.n*pos)
        # self.qc.sdg(self.n*pos)
        # self.qc.h(self.n*pos)
        i = np.random.randint(0,2)
        self.qc.rx(-np.pi/2, self.n*pos+i)

    @record
    def sqrt2_x(self, pos: int):
        self.magiccounter += 1
        # self.qc.h(self.n*pos)
        # self.qc.t(self.n*pos)
        # self.qc.h(self.n*pos)
        i = np.random.randint(0,2)
        self.qc.rx(np.pi/4, self.n*pos+i)

    @record
    def sqrt2_xdg(self, pos: int):
        self.magiccounter += 1
        # self.qc.h(self.n*pos)
        # self.qc.tdg(self.n*pos)
        # self.qc.h(self.n*pos)
        i = np.random.randint(0,2)
        self.qc.rx(-np.pi/4, self.n*pos+i)

    @record
    def s(self, pos: int):
        self.h(pos=pos)
        self.sqrt_x(pos=pos)
        self.h(pos=pos)

    @record
    def sdg(self, pos: int):
        self.h(pos=pos)
        self.sqrt_xdg(pos=pos)
        self.h(pos=pos)

    @record
    def t(self, pos: int):
        self.h(pos=pos)
        self.sqrt2_x(pos=pos)
        self.h(pos=pos)

    @record    
    def tdg(self, pos: int):
        self.h(pos=pos)
        self.sqrt2_xdg(pos=pos)
        self.h(pos=pos)

    @record
    def toff(self, control1: int, control2: int, targ: int):
        for i in range(self.n):
            for j in range(self.n):
                self.qc.ccx(self.n*control1 + i, self.n*control2 + j, self.n*targ + j)
            if self.err:
                self.qec_statevector(pos=targ)               #needed for FT
                # self.qec_counter -= 1

    @record
    def cnot(self, control: int, target: int):
        for i in range(self.n):
            self.qc.cx(self.n*control + i, self.n*target + i)

    @record
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

    @record
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

    @record
    def qec_statevector(self, pos: int):
        anc = self.qc.num_qubits - 1
        self.qec_counter += 1

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

        self.qc.save_statevector(label="psi")
        sim = AerSimulator(method="statevector", noise_model = self.noise_model)
        result = sim.run(self.qc, shots=1).result()
        psi_full = result.data(0)["psi"]

        bit = list(result.get_counts().keys())[0]           #Das ist das gesamte klassische Register, z.B. 0100...0101 , Achtung: Ganz rechts ist der erste Bit und ganz links der letzte Bit!

        qr2 = QuantumRegister(self.n*(self.logicalq+2)+3, "q")         #Qcirq+ClassicalRegister für H_L/QEC wiederherstellen und weitermachen
        self.qc = QuantumCircuit(qr2)
        self.qc.add_register(self.qecc)
        self.qc.set_statevector(psi_full)

        del psi_full

        #Nur die letzten 2 Bits nehmen, da wir nur 2 Stabilizer haben

        if bit[1:] == "01":                                     #S_0    strike
            self.qc.z(3*pos)
        elif bit[1:] == "11":                                   #S_0 S_1        strike
            self.qc.z(3*pos + 1)
        elif bit[1:] == "10":                                   #S_1        strike
            self.qc.z(3*pos + 2)

    def qec(self, pos: int):
        if self.n == 5:
            self.qec5(pos=pos)
            return
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

    def qec5(self, pos: int):
        anc = self.qc.num_qubits - 1
        self.qec_counter += 1

        self.qc.reset(anc)
        self.qc.h(anc)
        self.qc.cx(anc, 5*pos + 0)
        self.qc.cx(anc, 5*pos + 1)
        self.qc.h(anc)
        self.qc.measure(anc, self.qecc[0])

        self.qc.reset(anc)
        self.qc.h(anc)
        self.qc.cx(anc, 5*pos + 1)
        self.qc.cx(anc, 5*pos + 2)
        self.qc.h(anc)
        self.qc.measure(anc, self.qecc[1])

        self.qc.reset(anc)
        self.qc.h(anc)
        self.qc.cx(anc, 5*pos + 2)
        self.qc.cx(anc, 5*pos + 3)
        self.qc.h(anc)
        self.qc.measure(anc, self.qecc[2])

        self.qc.reset(anc)
        self.qc.h(anc)
        self.qc.cx(anc, 5*pos + 3)
        self.qc.cx(anc, 5*pos + 4)
        self.qc.h(anc)
        self.qc.measure(anc, self.qecc[3])

        ### Single-qubit errors #####

        with self.qc.if_test((self.qecc[0], 1)):                
            with self.qc.if_test((self.qecc[1], 0)):       
                with self.qc.if_test((self.qecc[2], 0)):  
                    with self.qc.if_test((self.qecc[3], 0)):        
                        self.qc.z(5*pos)

        with self.qc.if_test((self.qecc[0], 1)):                
            with self.qc.if_test((self.qecc[1], 1)):   
                with self.qc.if_test((self.qecc[2], 0)):    
                    with self.qc.if_test((self.qecc[3], 0)):  
                        self.qc.z(5*pos + 1)
        
        with self.qc.if_test((self.qecc[0], 0)):                
            with self.qc.if_test((self.qecc[1], 1)):               
                with self.qc.if_test((self.qecc[2], 1)):  
                    with self.qc.if_test((self.qecc[3], 0)):  
                        self.qc.z(5*pos + 2)

        with self.qc.if_test((self.qecc[0], 0)):  
            with self.qc.if_test((self.qecc[1], 0)):                
                with self.qc.if_test((self.qecc[2], 1)):     
                    with self.qc.if_test((self.qecc[3], 1)):             
                        self.qc.z(5*pos + 3)

        with self.qc.if_test((self.qecc[0], 0)):  
            with self.qc.if_test((self.qecc[1], 0)):   
                with self.qc.if_test((self.qecc[2], 0)):                
                    with self.qc.if_test((self.qecc[3], 1)):              
                        self.qc.z(5*pos + 4)

        ### Two-qubit errors #####

        with self.qc.if_test((self.qecc[0], 0)):                
            with self.qc.if_test((self.qecc[1], 1)):  
                with self.qc.if_test((self.qecc[2], 0)):  
                    with self.qc.if_test((self.qecc[3], 0)):            
                        self.qc.z(5*pos), self.qc.z(5*pos + 1)

        with self.qc.if_test((self.qecc[0], 1)):                
            with self.qc.if_test((self.qecc[1], 1)):  
                with self.qc.if_test((self.qecc[2], 1)):              
                    with self.qc.if_test((self.qecc[3], 0)): 
                        self.qc.z(5*pos), self.qc.z(5*pos + 2)

        with self.qc.if_test((self.qecc[0], 1)):                
            with self.qc.if_test((self.qecc[1], 0)):  
                with self.qc.if_test((self.qecc[2], 1)):              
                    with self.qc.if_test((self.qecc[3], 1)): 
                        self.qc.z(5*pos), self.qc.z(5*pos + 3)

        with self.qc.if_test((self.qecc[0], 1)):                
            with self.qc.if_test((self.qecc[1], 0)):  
                with self.qc.if_test((self.qecc[2], 0)):              
                    with self.qc.if_test((self.qecc[3], 1)): 
                        self.qc.z(5*pos), self.qc.z(5*pos + 4)
        
        with self.qc.if_test((self.qecc[0], 1)):                
            with self.qc.if_test((self.qecc[1], 0)):  
                with self.qc.if_test((self.qecc[2], 1)):              
                    with self.qc.if_test((self.qecc[3], 0)): 
                        self.qc.z(5*pos + 1), self.qc.z(5*pos + 2)
        
        with self.qc.if_test((self.qecc[0], 1)):                
            with self.qc.if_test((self.qecc[1], 1)):  
                with self.qc.if_test((self.qecc[2], 1)):              
                    with self.qc.if_test((self.qecc[3], 1)): 
                        self.qc.z(5*pos + 1), self.qc.z(5*pos + 3)
        
        with self.qc.if_test((self.qecc[0], 1)):                
            with self.qc.if_test((self.qecc[1], 1)):  
                with self.qc.if_test((self.qecc[2], 0)):              
                    with self.qc.if_test((self.qecc[3], 1)): 
                        self.qc.z(5*pos + 1), self.qc.z(5*pos + 4)

        with self.qc.if_test((self.qecc[0], 0)):                
            with self.qc.if_test((self.qecc[1], 1)):  
                with self.qc.if_test((self.qecc[2], 0)):              
                    with self.qc.if_test((self.qecc[3], 1)): 
                        self.qc.z(5*pos + 2), self.qc.z(5*pos + 3)

        with self.qc.if_test((self.qecc[0], 0)):                
            with self.qc.if_test((self.qecc[1], 1)):  
                with self.qc.if_test((self.qecc[2], 1)):              
                    with self.qc.if_test((self.qecc[3], 1)): 
                        self.qc.z(5*pos + 2), self.qc.z(5*pos + 4)
        
        with self.qc.if_test((self.qecc[0], 0)):                
            with self.qc.if_test((self.qecc[1], 0)):  
                with self.qc.if_test((self.qecc[2], 1)):              
                    with self.qc.if_test((self.qecc[3], 0)): 
                        self.qc.z(5*pos + 3), self.qc.z(5*pos + 4)
        
    def qec_ft(self, pos: int):
        anc = self.qc.num_qubits - 1
        ancc = anc - 1
        anccc = ancc - 1
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

        self.qc.reset(anc), self.qc.reset(ancc)                 #works for d=3
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

        # self.qc.reset(anc), self.qc.reset(ancc), self.qc.reset(anccc)                 
        # self.qc.h(anc), self.qc.h(ancc), self.qc.h(anccc)
        # self.qc.cx(anc, 3*pos + 0)
        # self.qc.cx(ancc, 3*pos + 1)
        # self.qc.cx(anccc, 3*pos + 2)
        # self.qc.cx(ancc, anc), self.qc.cx(ancc, anccc)
        # self.qc.h(anc), self.qc.h(anccc)
        # self.qc.id(anc), self.qc.id(anccc)
        # self.qc.measure(anc, self.qecc[0])
        # self.qc.measure(anccc, self.qecc[1])

        # with self.qc.if_test((self.qecc[0], 1)):                
        #     with self.qc.if_test((self.qecc[1], 0)):               
        #         self.qc.z(3*pos)

        # with self.qc.if_test((self.qecc[0], 1)):                
        #     with self.qc.if_test((self.qecc[1], 1)):               
        #         self.qc.z(3*pos + 1)
        
        # with self.qc.if_test((self.qecc[0], 0)):                
        #     with self.qc.if_test((self.qecc[1], 1)):               
        #         self.qc.z(3*pos + 2)

    def qec_ideal(self, pos: int):
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

    def readout_old(self, pos: int, shots: int):      #funktioniert auch, wenn man nur 1 Shot machen und dann extern drüber for loop macht lol

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

        sim = AerSimulator(method = "statevector", noise_model=self.noise_model)
        result = sim.run(self.qc, shots=shots).result()
        counts = result.get_counts()
            

        bitstring = list(counts.keys())
        bits = [i.replace(" ","") for i in bitstring]
        counter = list(counts.values())

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

        zero, one = 0, 0

        for i, val in enumerate(bits):
            if val == 0:
                zero += counter[i]
            elif val == 1:
                one += counter[i]

        self.ones += one/shots
        self.zeros += zero/shots

    def readout(self, pos: int, shots: int):        
            self.qc = None

            count0, count1 = [], []                 #alle statevectors für 0_L und 1_L
            for i in range(2**self.n):
                bit = inv_covert(i,self.n)
                if bit.count("1")%2 == 0:
                    count0.append(bit)
                else:
                    count1.append(bit)

            zero, one, err = 0, 0, 0
    
            for o in range(shots):
                circuit = RepCode_z(self.n, self.logicalq)
                for command in self.history:
                    method = getattr(circuit, command.method)
                    method(*command.args, **command.kwargs)
    
                read = ClassicalRegister(self.n)
                circuit.qc.add_register(read)
                
                for i in range(self.n):
                    circuit.qc.id(self.n*pos + i)
                    circuit.qc.measure(self.n*pos + i, read[self.n-1-i])

                sim = AerSimulator(method = "statevector", noise_model=self.noise_model)
                result = sim.run(circuit.qc, shots=1).result()
                counts = result.get_counts()

                bitstring = list(counts.keys())
                bits = [i.replace(" ","") for i in bitstring]
                counter = list(counts.values())
                
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
    
                for i, val in enumerate(bits):
                    if val == 0:
                        zero += counter[i]
                    elif val == 1:
                        one += counter[i]
    
            self.ones += one/shots
            self.zeros += zero/shots
