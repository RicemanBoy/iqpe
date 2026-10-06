from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister, qasm3, qasm2, qpy
from qiskit.visualization import plot_histogram
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
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

# try:
#     from .rep_code_decoder import decode_qec_block, data_readout_to_syndrome
# except ImportError:
#     from rep_code_decoder import decode_qec_block, data_readout_to_syndrome



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

def avg7_repcode(code: str, distance: int, iter: int, noise: float, qec = False, k = 1, bias = 0, path = ""):       #only exact angles!  
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
                    self.noise_model = self.__noise_model__(noise, bias)
                    self.err = qec
                    
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
                    self.readout(pos=0, shots=1)
                    # print("T/Tdg counter for iteration {}: {}".format(t, self.magiccounter))
            
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
    y = y/(7*k)
    arg = 0
    for i in range(len(bruh1)):
        arg += (y-bruh1[i])**2
    sigma = ((1/(k*n))*arg)**0.5
    sigma = sigma/((k*n)**0.5)

    return y_list

def avg7_ramsey(code: str, distance: int, iter: int, noise: float, qec = False, k = 1, bias = 0, path = ""):       #only exact angles!  
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
                        self = RepCode_z(distance, 1)
                    elif code == "x":
                        self = RepCode(distance, 1)
                    self.noise_model = self.__noise_model__(noise, bias)
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

                    self.readout(pos=0, shots=1)
            
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
    y = y/(7*k)
    arg = 0
    for i in range(len(bruh1)):
        arg += (y-bruh1[i])**2
    sigma = ((1/(k*n))*arg)**0.5
    sigma = sigma/((k*n)**0.5)

    return y_list

def avg7_ramsey_htoff(distance: int, iter: int, noise: float, qec = False, k = 1, bias = 0, path = ""):       #only exact angles!  
    n = 15
    angle = np.linspace(0,1,n+2)
    angle = np.delete(angle, [n+1])
    angle = np.delete(angle, [0])
    
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
                    self = RepCode_z(distance, 3)
                    self.noise_model = self.__noise_model__(noise, bias)
                    self.err = qec
                    
                    #############################
                    for _ in range(2**(iter-t-1)):
                        self.U_htoff(index=2*o+1)
                        self.U_htoff(index=2*o+1)
                    ###############################
                    for l in rots:
                        if l == 0.25:
                            self.adj_ry_pi_2()
                        if l == 0.125:
                            self.adj_ry_pi_4()

                    self.readout(pos=0, shots=1)
            
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
    y = y/(7*k)
    arg = 0
    for i in range(len(bruh1)):
        arg += (y-bruh1[i])**2
    sigma = ((1/(k*n))*arg)**0.5
    sigma = sigma/((k*n)**0.5)

    return y_list


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
   
    def z(self, pos: int):
        for i in range(self.n):
            self.qc.z(self.n*pos + i)

    def id(self, pos: int):
        for i in range(self.n):
            self.qc.id(self.n*pos + i)

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
            # self.qc.id(self.n*pos + i)
        for i in range(self.n):
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
        i = np.random.randint(0,self.n-1)
        self.qc.rx(np.pi/2, self.n*pos+i)

    def sqrt_xdg(self, pos: int):
        i = np.random.randint(0,self.n-1)
        self.qc.rx(-np.pi/2, self.n*pos+i)

    def sqrt2_x(self, pos: int):
        self.magiccounter += 1
        i = np.random.randint(0,self.n-1)
        self.qc.rx(np.pi/4, self.n*pos+i)

    def sqrt2_xdg(self, pos: int):
        self.magiccounter += 1
        i = np.random.randint(0,self.n-1)
        self.qc.rx(-np.pi/4, self.n*pos+i)

    def adj_ry_pi_2(self):
        #implements ry(-pi/2). more or less accurate
        [[('X', 0), ('H', 0)]]
        self.x(0)
        self.h(0)

    def adj_ry_pi_4(self):
        #implements ry(-pi/4) with 1e-2 error, error is defined over the algorithm
        #lists = [('X', 2), ('CCX', 0, 2, 1), ('X', 2), ('X', 1), ('X', 0), ('CCX', 0, 1, 2), ('X', 1), ('X', 0), ('H', 1), ('CCX', 1, 2, 0), ('H', 2), ('CCX', 0, 1, 2), ('CCX', 0, 2, 1), ('X', 1), ('CCX', 1, 2, 0), ('X', 1), ('CCX', 0, 2, 1), ('CCX', 0, 1, 2), ('CCX', 0, 2, 1), ('CCX', 0, 1, 2), ('H', 2), ('X', 2), ('CCX', 0, 2, 1), ('X', 2), ('X', 1), ('CCX', 1, 2, 0), ('X', 1), ('H', 2), ('X', 2), ('CCX', 0, 2, 1), ('X', 2), ('X', 1), ('CCX', 1, 2, 0), ('X', 1), ('H', 2), ('CCX', 0, 1, 2), ('CCX', 0, 2, 1), ('X', 1), ('CCX', 1, 2, 0), ('X', 1), ('CCX', 0, 2, 1), ('CCX', 0, 1, 2), ('CCX', 0, 2, 1), ('H', 2), ('CCX', 0, 1, 2), ('CCX', 0, 2, 1), ('X', 1), ('CCX', 1, 2, 0), ('X', 1), ('CCX', 0, 2, 1), ('CCX', 0, 1, 2), ('CCX', 0, 2, 1), ('CCX', 0, 1, 2), ('H', 2), ('CCX', 0, 1, 2), ('CCX', 0, 2, 1), ('X', 1), ('CCX', 1, 2, 0), ('X', 1), ('CCX', 0, 2, 1), ('CCX', 0, 1, 2), ('CCX', 0, 2, 1), ('CCX', 0, 1, 2), ('H', 2), ('CCX', 0, 1, 2), ('CCX', 0, 2, 1), ('X', 1), ('CCX', 1, 2, 0), ('X', 1), ('CCX', 0, 2, 1), ('CCX', 0, 1, 2), ('CCX', 1, 2, 0), ('CCX', 0, 1, 2), ('H', 0), ('CCX', 0, 1, 2), ('CCX', 1, 2, 0), ('X', 1), ('CCX', 1, 2, 0), ('X', 1), ('CCX', 0, 2, 1), ('CCX', 0, 1, 2), ('H', 2), ('CCX', 0, 1, 2), ('CCX', 0, 2, 1), ('X', 1), ('CCX', 1, 2, 0), ('X', 1), ('CCX', 0, 1, 2), ('H', 2), ('X', 2), ('CCX', 0, 2, 1), ('X', 2), ('CCX', 1, 2, 0), ('X', 0), ('CCX', 0, 1, 2), ('X', 0), ('X', 1), ('CCX', 1, 2, 0), ('X', 1), ('H', 0), ('CCX', 0, 1, 2), ('CCX', 0, 2, 1), ('CCX', 1, 2, 0), ('CCX', 0, 1, 2), ('H', 2), ('X', 0), ('CCX', 0, 2, 1), ('X', 0), ('X', 1), ('CCX', 1, 2, 0), ('X', 1), ('H', 2), ('CCX', 0, 1, 2), ('CCX', 0, 2, 1), ('X', 1), ('CCX', 1, 2, 0), ('X', 1), ('CCX', 0, 2, 1), ('CCX', 0, 1, 2), ('CCX', 0, 2, 1), ('H', 0), ('CCX', 0, 1, 2), ('CCX', 0, 2, 1), ('CCX', 1, 2, 0), ('CCX', 0, 1, 2), ('H', 1), ('CCX', 0, 1, 2), ('CCX', 0, 2, 1), ('X', 1), ('CCX', 1, 2, 0), ('X', 1), ('CCX', 0, 2, 1), ('CCX', 0, 1, 2), ('CCX', 1, 2, 0), ('H', 2), ('X', 2), ('X', 1), ('CCX', 1, 2, 0), ('X', 2), ('X', 1), ('CCX', 0, 1, 2), ('CCX', 0, 2, 1), ('CCX', 0, 1, 2), ('X', 0), ('CCX', 0, 2, 1), ('X', 0), ('H', 1), ('CCX', 0, 1, 2), ('CCX', 0, 2, 1), ('X', 1), ('CCX', 1, 2, 0), ('X', 1), ('CCX', 0, 2, 1), ('CCX', 0, 1, 2), ('CCX', 0, 2, 1), ('CCX', 1, 2, 0), ('CCX', 0, 2, 1), ('H', 0), ('CCX', 0, 1, 2), ('CCX', 0, 2, 1), ('X', 1), ('CCX', 1, 2, 0), ('X', 1), ('CCX', 0, 2, 1), ('CCX', 0, 1, 2), ('CCX', 0, 2, 1), ('CCX', 1, 2, 0), ('CCX', 0, 2, 1), ('H', 1), ('X', 2), ('X', 1), ('CCX', 1, 2, 0), ('X', 2), ('X', 1), ('CCX', 0, 1, 2), ('CCX', 1, 2, 0), ('CCX', 0, 2, 1), ('CCX', 0, 1, 2), ('H', 2), ('CCX', 0, 1, 2), ('CCX', 1, 2, 0), ('CCX', 0, 1, 2)]
        lists = [('X', 1), ('CX', 0, 2), ('CX', 0, 1), ('X', 0), ('Z', 0), ('X', 0), ('H', 0), ('CX', 0, 2), ('CCX', 1, 2, 0), ('CX', 0, 1), ('H', 1), ('CCX', 0, 1, 2), ('CX', 2, 1), ('CCX', 1, 2, 0), ('CX', 0, 2), ('H', 0), ('CX', 0, 1), ('CCX', 0, 1, 2), ('CCX', 1, 2, 0), ('H', 1), ('CCX', 1, 2, 0), ('CCX', 0, 1, 2), ('CX', 1, 0), ('H', 1), ('CCX', 1, 2, 0), ('CCX', 0, 1, 2), ('CX', 0, 1), ('H', 0), ('CCX', 0, 1, 2), ('CCX', 0, 2, 1), ('CX', 0, 2), ('CX', 1, 0), ('H', 0), ('CCX', 0, 1, 2), ('CCX', 0, 2, 1), ('CX', 0, 2), ('CX', 1, 0), ('CX', 2, 1), ('H', 0), ('CCX', 0, 1, 2), ('CCX', 0, 2, 1), ('CX', 0, 1), ('H', 0), ('CCX', 0, 2, 1), ('CX', 0, 1), ('CCX', 0, 1, 2), ('H', 0), ('CCX', 0, 1, 2), ('CCX', 0, 2, 1), ('CX', 0, 1), ('H', 0), ('CCX', 0, 1, 2), ('CCX', 0, 2, 1), ('CX', 0, 1), ('H', 0), ('CX', 0, 1), ('CCX', 0, 1, 2), ('CX', 2, 0), ('CCX', 0, 2, 1), ('H', 2), ('CCX', 0, 2, 1), ('CCX', 1, 2, 0), ('CX', 1, 2), ('H', 1), ('CCX', 1, 2, 0), ('CCX', 0, 1, 2), ('CX', 0, 1), ('H', 0), ('CCX', 0, 2, 1), ('CCX', 0, 1, 2), ('CX', 1, 0), ('H', 1), ('CCX', 0, 1, 2), ('CCX', 1, 2, 0), ('CX', 2, 1), ('H', 2), ('CCX', 1, 2, 0), ('CCX', 0, 2, 1), ('CX', 0, 2), ('H', 0), ('CCX', 0, 1, 2), ('CX', 2, 1), ('CX', 0, 2), ('CCX', 1, 2, 0), ('H', 1), ('CCX', 1, 2, 0), ('CX', 0, 1), ('CCX', 0, 1, 2), ('CX', 2, 0), ('H', 0), ('CX', 0, 1), ('CCX', 0, 1, 2), ('CX', 2, 0), ('CCX', 0, 2, 1), ('H', 2), ('CCX', 1, 2, 0), ('CX', 0, 2)]


        for i in lists:
            if i[0] == "X":
                self.x(pos=i[1])
            elif i[0] == "H":
                self.h(pos=i[1])
            elif i[0] == "Z":
                self.z(pos=i[1])
            elif i[0] == "CCX":
                self.toff(i[1],i[2],i[3])
            elif i[0] == "CX":
                self.cnot(i[1],i[2])

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
  
    def U_htoff(self, index=0):             #latest solution, this implements sqrt(U), so we have to apply this twice!
        with open("htoffcx15_1e-2.txt", "r") as f:
            htoff_decomp = [
                [tuple(item.split("_")[0:1] + [int(x) for x in item.split("_")[1:]])
                for item in line.strip().split(",")]
                for line in f
            ]
        for i in htoff_decomp[index]:
            if i[0] == "X":
                self.x(pos=i[1])

            elif i[0] == "H":
                self.h(pos=i[1])

            elif i[0] == "Z":
                self.z(pos=i[1])

            elif i[0] == "CCX":
                self.toff(i[1],i[2],i[3])

            elif i[0] == "CX":
                self.cnot(i[1],i[2])
            else:
                print("MISSING GATE")

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

    def cu_ramsey(self, gate: list):
        self.u2(0, gate=gate)
        self.u2(0, gate=gate)

    def cu_ramsey_htoff(self, gate: list):
            self.u2_htoff(0, gate=gate)
            self.u2_htoff(0, gate=gate)

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

    def qec5(self, pos: int):           #qec for d=5 repcode
        anc = self.qc.num_qubits - 1
        self.qec_counter += 1

        self.qc.reset(anc)
        self.qc.h(anc)
        self.qc.cx(anc, 5*pos + 0)
        self.qc.cx(anc, 5*pos + 1)
        self.qc.h(anc)
        self.qc.id(anc)
        self.qc.measure(anc, self.qecc[0])

        self.qc.reset(anc)
        self.qc.h(anc)
        self.qc.cx(anc, 5*pos + 1)
        self.qc.cx(anc, 5*pos + 2)
        self.qc.h(anc)
        self.qc.id(anc)
        self.qc.measure(anc, self.qecc[1])

        self.qc.reset(anc)
        self.qc.h(anc)
        self.qc.cx(anc, 5*pos + 2)
        self.qc.cx(anc, 5*pos + 3)
        self.qc.h(anc)
        self.qc.id(anc)
        self.qc.measure(anc, self.qecc[2])

        self.qc.reset(anc)
        self.qc.h(anc)
        self.qc.cx(anc, 5*pos + 3)
        self.qc.cx(anc, 5*pos + 4)
        self.qc.h(anc)
        self.qc.id(anc)
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

    def qec5_ideal(self, pos: int):         #ideal qec for d=5 repcode
        anc = self.qc.num_qubits - 1
        self.qec_counter += 1

        self.qc.reset(anc)
        self.qc.append(h_ideal, [anc])
        self.qc.append(cx_ideal, [5*pos + 0, anc])
        self.qc.append(cx_ideal, [5*pos + 1, anc])
        self.qc.append(h_ideal, [anc])
        self.qc.measure(anc, self.qecc[0])

        self.qc.reset(anc)
        self.qc.append(h_ideal, [anc])
        self.qc.append(cx_ideal, [5*pos + 1, anc])
        self.qc.append(cx_ideal, [5*pos + 2, anc])
        self.qc.append(h_ideal, [anc])
        self.qc.measure(anc, self.qecc[1])

        self.qc.reset(anc)
        self.qc.append(h_ideal, [anc])
        self.qc.append(cx_ideal, [5*pos + 2, anc])
        self.qc.append(cx_ideal, [5*pos + 3, anc])
        self.qc.append(h_ideal, [anc])
        self.qc.measure(anc, self.qecc[2])

        self.qc.reset(anc)
        self.qc.append(h_ideal, [anc])
        self.qc.append(cx_ideal, [5*pos + 3, anc])
        self.qc.append(cx_ideal, [5*pos + 4, anc])
        self.qc.append(h_ideal, [anc])
        self.qc.measure(anc, self.qecc[3])

        ### Single-qubit errors #####

        with self.qc.if_test((self.qecc[0], 1)):                
            with self.qc.if_test((self.qecc[1], 0)):       
                with self.qc.if_test((self.qecc[2], 0)):  
                    with self.qc.if_test((self.qecc[3], 0)):        
                        self.qc.append(z_ideal, [5*pos+0])

        with self.qc.if_test((self.qecc[0], 1)):                
            with self.qc.if_test((self.qecc[1], 1)):   
                with self.qc.if_test((self.qecc[2], 0)):    
                    with self.qc.if_test((self.qecc[3], 0)):  
                        self.qc.append(z_ideal, [5*pos+1])
        
        with self.qc.if_test((self.qecc[0], 0)):                
            with self.qc.if_test((self.qecc[1], 1)):               
                with self.qc.if_test((self.qecc[2], 1)):  
                    with self.qc.if_test((self.qecc[3], 0)):  
                        self.qc.append(z_ideal, [5*pos+2])

        with self.qc.if_test((self.qecc[0], 0)):  
            with self.qc.if_test((self.qecc[1], 0)):                
                with self.qc.if_test((self.qecc[2], 1)):     
                    with self.qc.if_test((self.qecc[3], 1)):             
                        self.qc.append(z_ideal, [5*pos+3])

        with self.qc.if_test((self.qecc[0], 0)):  
            with self.qc.if_test((self.qecc[1], 0)):   
                with self.qc.if_test((self.qecc[2], 0)):                
                    with self.qc.if_test((self.qecc[3], 1)):              
                        self.qc.append(z_ideal, [5*pos+4])

        ### Two-qubit errors  --> 5*4/2 = 10 Possiblities #####

        with self.qc.if_test((self.qecc[0], 0)):                
            with self.qc.if_test((self.qecc[1], 1)):  
                with self.qc.if_test((self.qecc[2], 0)):  
                    with self.qc.if_test((self.qecc[3], 0)):            
                        self.qc.append(z_ideal, [5*pos+0]), self.qc.append(z_ideal, [5*pos+1])

        with self.qc.if_test((self.qecc[0], 1)):                
            with self.qc.if_test((self.qecc[1], 1)):  
                with self.qc.if_test((self.qecc[2], 1)):              
                    with self.qc.if_test((self.qecc[3], 0)): 
                        self.qc.append(z_ideal, [5*pos+0]), self.qc.append(z_ideal, [5*pos+2])

        with self.qc.if_test((self.qecc[0], 1)):                
            with self.qc.if_test((self.qecc[1], 0)):  
                with self.qc.if_test((self.qecc[2], 1)):              
                    with self.qc.if_test((self.qecc[3], 1)): 
                        self.qc.append(z_ideal, [5*pos+0]), self.qc.append(z_ideal, [5*pos+3])

        with self.qc.if_test((self.qecc[0], 1)):                
            with self.qc.if_test((self.qecc[1], 0)):  
                with self.qc.if_test((self.qecc[2], 0)):              
                    with self.qc.if_test((self.qecc[3], 1)): 
                        self.qc.append(z_ideal, [5*pos+0]), self.qc.append(z_ideal, [5*pos+4])
        
        with self.qc.if_test((self.qecc[0], 1)):                
            with self.qc.if_test((self.qecc[1], 0)):  
                with self.qc.if_test((self.qecc[2], 1)):              
                    with self.qc.if_test((self.qecc[3], 0)): 
                        self.qc.append(z_ideal, [5*pos+1]), self.qc.append(z_ideal, [5*pos+2])
        
        with self.qc.if_test((self.qecc[0], 1)):                
            with self.qc.if_test((self.qecc[1], 1)):  
                with self.qc.if_test((self.qecc[2], 1)):              
                    with self.qc.if_test((self.qecc[3], 1)): 
                        self.qc.append(z_ideal, [5*pos+1]), self.qc.append(z_ideal, [5*pos+3])
        
        with self.qc.if_test((self.qecc[0], 1)):                
            with self.qc.if_test((self.qecc[1], 1)):  
                with self.qc.if_test((self.qecc[2], 0)):              
                    with self.qc.if_test((self.qecc[3], 1)): 
                        self.qc.append(z_ideal, [5*pos+1]), self.qc.append(z_ideal, [5*pos+4])

        with self.qc.if_test((self.qecc[0], 0)):                
            with self.qc.if_test((self.qecc[1], 1)):  
                with self.qc.if_test((self.qecc[2], 0)):              
                    with self.qc.if_test((self.qecc[3], 1)): 
                        self.qc.append(z_ideal, [5*pos+2]), self.qc.append(z_ideal, [5*pos+3])

        with self.qc.if_test((self.qecc[0], 0)):                
            with self.qc.if_test((self.qecc[1], 1)):  
                with self.qc.if_test((self.qecc[2], 1)):              
                    with self.qc.if_test((self.qecc[3], 1)): 
                        self.qc.append(z_ideal, [5*pos+2]), self.qc.append(z_ideal, [5*pos+4])
        
        with self.qc.if_test((self.qecc[0], 0)):                
            with self.qc.if_test((self.qecc[1], 0)):  
                with self.qc.if_test((self.qecc[2], 1)):              
                    with self.qc.if_test((self.qecc[3], 0)): 
                        self.qc.append(z_ideal, [5*pos+3]), self.qc.append(z_ideal, [5*pos+4])

    def syndromes(self, pos: int):                 #ONLY syndrome extraction for d=3
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
        
    def syndromes5(self, pos: int):                 #ONLY syndrome extraction for d=5
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

    def qec_block(self, pos: int):                     #n Blocks of syndrome extraction, followed by decoding and correction

        spacetime_syndromes = []

        for i in range(self.n):
            self.syndromes(pos=pos)
            self.qc.save_statevector(label="psi")
            sim = AerSimulator(method="statevector", noise_model = self.noise_model)
            result = sim.run(self.qc, shots=1).result()
            psi_full = result.data(0)["psi"]

            bit = list(result.get_counts().keys())[0]           #Das ist das gesamte klassische Register, z.B. 0100...0101 , Achtung: Ganz rechts ist der erste Bit und ganz links der letzte Bit!

            spacetime_syndromes.append(bit[1:])

            qr2 = QuantumRegister(self.n*(self.logicalq+2)+3, "q")         #Qcirq+ClassicalRegister für H_L/QEC wiederherstellen und weitermachen
            self.qc = QuantumCircuit(qr2)
            self.qc.add_register(self.qecc)
            self.qc.set_statevector(psi_full)

            del psi_full

        # print(spacetime_syndromes)

        correction = decode_qec_block(spacetime_syndromes)

        # print("Syndromes: ", spacetime_syndromes, "  |   Korrektur hier: ", correction)

        for i in correction:
            self.qc.z(self.n*pos + i)

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

    def readout(self, pos: int, shots: int):      #funktioniert, wenn man nur 1 Shot machen und dann extern drüber for loop macht lol
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

        sim = AerSimulator(noise_model=self.noise_model)
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

