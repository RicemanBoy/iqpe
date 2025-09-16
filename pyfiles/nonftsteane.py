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

#https://www.nature.com/articles/srep19578 , https://arxiv.org/pdf/1106.2190 Figure 5b)
def code_goto(n=2):             #encodes |00>_L
    qr = QuantumRegister(7*n+2,"q")
    cbits = ClassicalRegister(8, "c")      #Cbit 1-6: Ancillas, Cbit 7-8: Preselection, nicht mehr überschreiben bis Readout!
    qc = QuantumCircuit(qr, cbits)
    
    anc = qc.num_qubits - 1

    for i in range(7*(n-1)):                        #start noise
        qc.id(i)

    for i in range(n-1):
        qc.h(1+7*i)
        qc.h(2+7*i)
        qc.h(3+7*i)

        qc.cx(1+7*i,0+7*i)
        qc.cx(3+7*i,5+7*i)

        qc.cx(2+7*i,6+7*i)

        qc.cx(1+7*i,4+7*i)

        qc.cx(2+7*i,0+7*i)
        qc.cx(3+7*i,6+7*i)

        qc.cx(1+7*i,5+7*i)

        qc.cx(6+7*i,4+7*i)

        qc.cx(0+7*i,anc)
        qc.cx(5+7*i,anc)
        qc.cx(6+7*i,anc)

        qc.id(anc)
        qc.measure(anc,cbits[8-i-1])      
    return qc

def gates(qc:QuantumCircuit):
    hmm = dict(qc.count_ops())
    hmm["reset"] = 0
    hmm["measure"] = 0
    hmm["if_else"] = 0
    return print("Amount of gates in this circuit: ", sum(hmm.values()))

def X_L(qc: QuantumCircuit, pos: int):
    qc.x(0+7*pos)
    qc.x(1+7*pos)
    qc.x(2+7*pos)

def Z_L(qc: QuantumCircuit, pos: int):
    qc.z(0+7*pos)
    qc.z(1+7*pos)
    qc.z(2+7*pos)
    
def H_L(qc: QuantumCircuit, pos: int):
    for i in range(7):
        qc.h(i+7*pos)

def S_L(qc: QuantumCircuit, pos: int):
    # for i in range(7):
    #     qc.z(i+7*pos)
    #     qc.s(i+7*pos)
    qc.s(0+7*pos), qc.s(1+7*pos), qc.s(3+7*pos), qc.s(6+7*pos)
    qc.sdg(2+7*pos), qc.sdg(4+7*pos), qc.sdg(5+7*pos)

def CZ_L(qc: QuantumCircuit):
    H_L(qc, 0)
    CNOT_L(qc, 1)
    H_L(qc, 0)

def adj_S_L(qc: QuantumCircuit, pos: int):
    # for i in range(7):
    #     qc.z(i+7*pos)
    #     qc.sdg(i+7*pos)
    qc.sdg(0+7*pos), qc.sdg(1+7*pos), qc.sdg(3+7*pos), qc.sdg(6+7*pos)
    qc.s(2+7*pos), qc.s(4+7*pos), qc.s(5+7*pos)

def CNOT_L(qc: QuantumCircuit, control: int):
    if control == 0:
        for i in range(7):
            qc.cx(i, i+7)
    else:
        for i in range(7):
            qc.cx(i+7,i)

def CS_L(qc: QuantumCircuit, control: int, target: int):
    T_L(qc, 0)
    T_L(qc, 1)
    CNOT_L(qc, control=control)
    adj_T_L(qc, pos = target)
    CNOT_L(qc, control=control)

def T_L(qc: QuantumCircuit, pos: int):
    anc = qc.num_qubits - 1
    qc.reset(anc)

    # qc.append(h_ideal,[anc])
    # qc.tdg(anc)

    qc.id(anc)
    qc.h(anc)
    qc.t(anc)

    qc.cx(7*pos+0,anc)
    qc.cx(7*pos+1,anc)
    qc.cx(7*pos+2,anc)

    qc.measure(anc,0)

    with qc.if_test((0,1)):
        qc.s(0+7*pos), qc.s(1+7*pos), qc.s(3+7*pos), qc.s(6+7*pos)
        qc.sdg(2+7*pos), qc.sdg(4+7*pos), qc.sdg(5+7*pos)
    
def adj_T_L(qc: QuantumCircuit, pos: int):
    anc = qc.num_qubits - 1
    qc.reset(anc)

    qc.id(anc)
    qc.h(anc)
    qc.tdg(anc)

    qc.cx(7*pos+0,anc)
    qc.cx(7*pos+1,anc)
    qc.cx(7*pos+2,anc)

    qc.measure(anc,0)

    with qc.if_test((0,1)):
        qc.sdg(0+7*pos), qc.sdg(1+7*pos), qc.sdg(3+7*pos), qc.sdg(6+7*pos)
        qc.s(2+7*pos), qc.s(4+7*pos), qc.s(5+7*pos)

circ = QuantumCircuit(1)
circ.rz(np.pi/8, 0)
basis = ["t", "tdg", "z", "h"]
approx = generate_basic_approximations(basis, depth=3)
skd = SolovayKitaev(recursion_degree=2, basic_approximations=approx)
rootT = skd(circ)

def root_T_L(qc: QuantumCircuit, pos: int, err = False):
    instruction = rootT.data
    counter = 0
    for i in instruction:
        if counter%3 == 0:
            if err:
                qec(qc, pos=pos)
        if i.name == "t":
            T_L(qc, pos=pos)
            counter += 1
        if i.name == "tdg":
            adj_T_L(qc, pos=pos)
            counter += 1
        if i.name == "h":
            H_L(qc, pos=pos)
        if i.name == "s":
            S_L(qc, pos=pos)
        if i.name == "sdg":
            adj_S_L(qc, pos=pos)
        if i.name == "z":
            Z_L(qc, pos)

circ = QuantumCircuit(1)
circ.rz(-np.pi/8, 0)
basis = ["t", "tdg", "z", "h"]
approx = generate_basic_approximations(basis, depth=3)
skd = SolovayKitaev(recursion_degree=2, basic_approximations=approx)
adj_rootT = skd(circ)

def adj_root_T_L(qc: QuantumCircuit, pos: int, err=False):
    instruction = adj_rootT.data
    counter = 0
    for i in instruction:
        if counter%3 == 0:
            if err:
                qec(qc, pos=pos)
        if i.name == "t":
            T_L(qc, pos=pos)
            counter += 1
        if i.name == "tdg":
            adj_T_L(qc, pos=pos)
            counter += 1
        if i.name == "h":
            H_L(qc, pos=pos)

def CT_L(qc: QuantumCircuit, err=False):
    root_T_L(qc, 0, err=err)
    root_T_L(qc, 1, err=err)
    CNOT_L(qc, 0)
    adj_root_T_L(qc, 1, err=err)
    CNOT_L(qc, 0)
#######################################################################
def convert(bin: str):                  #konvertiert den bitstring in decimal, e.g. 0110 = 0.375
    k = list(bin)
    a = [int(i) for i in k]
    n = 0
    for i in range(len(a)):
        if a[i] == 1:
            n += 1/2**(i+1)
    return n

def U2(qc: QuantumCircuit, pos: int, gate: list):
    for i in gate:
        if i == "s":
            S_L(qc, pos=pos)
        if i == "sdg":
            adj_S_L(qc, pos=pos)
        if i == "t":
            T_L(qc, pos=pos)
        if i == "tdg":
            adj_T_L(qc, pos=pos)
        if i == "h":
            H_L(qc, pos=pos)
        if i == "z":
            Z_L(qc, pos=pos)

def CU_L(qc: QuantumCircuit, Ugates: list, adjUgates: list):
    U2(qc, 0, Ugates)
    U2(qc, 1, Ugates)
    CNOT_L(qc, control=0)
    U2(qc, 1, adjUgates)
    CNOT_L(qc, control=0)

def avg15(iter: int, n:int, argh: float, err = False, k = 1):       #each iteration own circuit
    angle = np.linspace(0,1,n+2)
    angle = np.delete(angle, [n+1])
    angle = np.delete(angle, [0])

    a, b = [], []
    with open("text/unitary{}.txt".format(n), "r") as file:
        for line in file:
            a.append(list(map(str, line.strip().split(","))))
    with open("text/adjunitary{}.txt".format(n), "r") as file:
        for line in file:
            b.append(list(map(str, line.strip().split(","))))
    
    y = 0
    bruh1 = []
    global hads
    for m in range(k):
        for o in range(n):
            bitstring = ""
            rots = []
            for t in range(iter):
                rots = [k*0.5 for k in rots]
                while True:
                    qc = code_goto()

                    X_L(qc,1)
                    H_L(qc,0)
                    #############################
                    for j in range(2**(iter-t-1)):
                        CU_L(qc, a[o], b[o], err=err)
                    ###############################
                    for l in rots:
                        if l == 0.25:
                            adj_S_L(qc, pos=0)
                        if l == 0.125:
                            adj_T_L(qc, pos=0)
                    H_L(qc, pos=0)
                    if err:
                        qec(qc, pos=0)
                    zeros, ones, _,_ = readout(qc, pos=0, shots=1, noise=argh)
            
                    if zeros == 1:
                        bitstring += "0"
                        break
                    if ones == 1:
                        bitstring += "1"
                        rots.append(0.5)
                        break
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
#################################################################

def readout(qc: QuantumCircuit, pos: int, shots: int, noise = 0):
    p = noise
    p_error = pauli_error([["X",p/2],["I",1-p],["Z",p/2]])
    p_error_2 = pauli_error([["XI",p/4],["IX",p/4],["II",1-p],["ZI",p/4],["IZ",p/4]])

    noise_model = NoiseModel()
    noise_model.add_all_qubit_quantum_error(p_error, ['x', "z", 'h', "s", "sdg", "id"])  # Apply to single-qubit gates
    noise_model.add_all_qubit_quantum_error(p_error_2, ['cx'])  # Apply to 2-qubit gates

    read = ClassicalRegister(7)
    qc.add_register(read)

    for i in range(7):
        qc.id(i+7*pos)
        qc.measure(i+7*pos,read[6-i])

    sim = AerSimulator()
    
    job = sim.run(qc, shots=shots, noise_model=noise_model)

    result = job.result()
    counts = result.get_counts()

    #print(counts)

    print(counts)

    bitstring = list(counts.keys())
    bitstring = [i.replace(" ","") for i in bitstring]


    hmm = list(counts.values())

    allcbits = len(bitstring[0])                
    pre, preselected = [i[allcbits-8:allcbits-6] for i in bitstring], 0
    bits = [i[:7] for i in bitstring]

    print(pre)

    #print(bits)
    #print(postprocess)

    for i in range(len(pre)):
        if pre[i].count("1") != 0:
            bits[i] = "pre"
            preselected += hmm[i]

    test_0 = ["0000000","1010101","0110011","1100110","0001111","1011010","0111100","1101001"]
    test_1 = ["1111111","0101010","1001100","0011001","1110000","0100101","1000011","0010110"]

    ############
    #Klassische Error-correction hier noch rein programmieren!
    ############

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
    #magic = (magic/shots)

    # print("0: ", zeros*100, "%")
    # print("1: ", ones*100, "%")
    # print("Preselection discarded: ", (preselected/shots)*100, "%")
    # print("Postselection discarded: ", (post/shots)*100, "%")
    return zeros, ones, preselected, post#,magic

def qec(qc: QuantumCircuit, pos: int):
    anc = qc.num_qubits - 1
    qc.reset(anc)
    qc.id(anc)
    ##################################Z-Stabilizers##########################################
    qc.cx(0+7*pos, anc)
    qc.cx(2+7*pos, anc)
    qc.cx(4+7*pos, anc)
    qc.cx(6+7*pos, anc)

    qc.id(anc)
    qc.measure(anc, 0)
    qc.reset(anc)
    qc.id(anc)

    qc.cx(1+7*pos, anc)
    qc.cx(2+7*pos, anc)
    qc.cx(5+7*pos, anc)
    qc.cx(6+7*pos, anc)

    qc.id(anc)
    qc.measure(anc, 1)
    qc.reset(anc)
    qc.id(anc)

    qc.cx(3+7*pos, anc)
    qc.cx(4+7*pos, anc)
    qc.cx(5+7*pos, anc)
    qc.cx(6+7*pos, anc)

    qc.id(anc)
    qc.measure(anc, 2)
    qc.reset(anc)
    qc.id(anc)
    ##################################X-Stabilizers##############################################
    qc.h(anc)
    qc.cx(anc, 0+7*pos)
    qc.cx(anc, 2+7*pos)
    qc.cx(anc, 4+7*pos)
    qc.cx(anc, 6+7*pos)
    qc.h(anc)

    qc.id(anc)
    qc.measure(anc, 3)
    qc.reset(anc)
    qc.id(anc)

    qc.h(anc)
    qc.cx(anc, 1+7*pos)
    qc.cx(anc, 2+7*pos)
    qc.cx(anc, 5+7*pos)
    qc.cx(anc, 6+7*pos)
    qc.h(anc)

    qc.id(anc)
    qc.measure(anc, 4)
    qc.reset(anc)
    qc.id(anc)

    qc.h(anc)
    qc.cx(anc, 3+7*pos)
    qc.cx(anc, 4+7*pos)
    qc.cx(anc, 5+7*pos)
    qc.cx(anc, 6+7*pos)
    qc.h(anc)

    qc.id(anc)
    qc.measure(anc, 5)
    qc.reset(anc)
    ##################################Bitflip Error correction##############################################
    
    with qc.if_test((2,0)):             #qbit 0
        with qc.if_test((1,0)):
            with qc.if_test((0,1)):
                qc.x(0+7*pos)

    with qc.if_test((2,0)):             #qbit 1
        with qc.if_test((1,1)):
            with qc.if_test((0,0)):
                qc.x(1+7*pos)
    
    with qc.if_test((2,0)):             #qbit 2
        with qc.if_test((1,1)):
            with qc.if_test((0,1)):
                qc.x(2+7*pos)
    
    with qc.if_test((2,1)):             #qbit 3
        with qc.if_test((1,0)):
            with qc.if_test((0,0)):
                qc.x(3+7*pos)
    
    with qc.if_test((2,1)):             #qbit 4
        with qc.if_test((1,0)):
            with qc.if_test((0,1)):
                qc.x(4+7*pos)
    
    with qc.if_test((2,1)):             #qbit 5
        with qc.if_test((1,1)):
            with qc.if_test((0,0)):
                qc.x(5+7*pos)
    
    with qc.if_test((2,1)):             #qbit 6
        with qc.if_test((1,1)):
            with qc.if_test((0,1)):
                qc.x(6+7*pos)

    ##################################Phaseflip Error correction##############################################
    
    with qc.if_test((5,0)):             #qbit 0
        with qc.if_test((4,0)):
            with qc.if_test((3,1)):
                qc.z(0+7*pos)

    with qc.if_test((5,0)):             #qbit 1
        with qc.if_test((4,1)):
            with qc.if_test((3,0)):
                qc.z(1+7*pos)
    
    with qc.if_test((5,0)):             #qbit 2
        with qc.if_test((4,1)):
            with qc.if_test((3,1)):
                qc.z(2+7*pos)
    
    with qc.if_test((5,1)):             #qbit 3
        with qc.if_test((4,0)):
            with qc.if_test((3,0)):
                qc.z(3+7*pos)
    
    with qc.if_test((5,1)):             #qbit 4
        with qc.if_test((4,0)):
            with qc.if_test((3,1)):
                qc.z(4+7*pos)
    
    with qc.if_test((5,1)):             #qbit 5
        with qc.if_test((4,1)):
            with qc.if_test((3,0)):
                qc.z(5+7*pos)
    
    with qc.if_test((5,1)):             #qbit 6
        with qc.if_test((4,1)):
            with qc.if_test((3,1)):
                qc.z(6+7*pos)

