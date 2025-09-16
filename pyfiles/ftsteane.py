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
def code_goto(n=3):             #encodes |00>_L
    qr = QuantumRegister(7*n+2,"q")
    cbits = ClassicalRegister(3, "c")
    qc = QuantumCircuit(qr, cbits)
    
    anc = qc.num_qubits - 1

    for i in range(7*(n-1)):                        #start noise
        qc.id(i)

    for i in range(n-1):
        # qc.h(0+7*i)                           #funktioniert auch, aber der goto code geht auch, bin bloß dumm
        # qc.h(1+7*i)
        # qc.h(3+7*i)

        # qc.cx(0+7*i,2+7*i)
        # qc.cx(3+7*i,4+7*i)

        # qc.cx(1+7*i,5+7*i)

        # qc.cx(0+7*i,6+7*i)

        # qc.cx(1+7*i,2+7*i)
        # qc.cx(3+7*i,5+7*i)

        # qc.cx(0+7*i,4+7*i)

        # qc.cx(5+7*i,6+7*i)

        # qc.cx(2+7*i,anc)
        # qc.cx(4+7*i,anc)
        # qc.cx(5+7*i,anc)

        #qc.h(7)                            #directly implements a t gate onto the steane code
        #qc.t(7)

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
        qc.measure(anc,cbits[i+1])      
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

# def Ty_L(qc: QuantumCircuit, pos: int):
#     state_inj = ClassicalRegister(14)
#     qc.add_register(state_inj)

#     anc = qc.num_qubits - 1
#     ancc = anc - 1

#     for i in range(7):
#         qc.reset(i+7*2)

#     for i in range(7):                        #start noise
#         qc.id(i+7*2)

#     qc.h(0+7*2)
#     qc.h(1+7*2)
#     qc.ry(np.pi/4,2+7*2)
#     qc.h(3+7*2)

#     qc.cx(2+7*2,4+7*2)
#     qc.cx(0+7*2,6+7*2)

#     qc.cx(3+7*2,5+7*2)

#     qc.cx(2+7*2,5+7*2)

#     qc.cx(0+7*2,4+7*2)
#     qc.cx(1+7*2,6+7*2)

#     qc.cx(0+7*2,2+7*2)

#     qc.cx(1+7*2,5+7*2)

#     qc.cx(1+7*2,2+7*2)
#     qc.cx(3+7*2,4+7*2)
#     qc.cx(3+7*2,6+7*2)
#     #################################Controlled Hadamards##########################################
#     qc.reset(anc), qc.reset(ancc)
#     qc.h(ancc)
#     for i in range(7):
#         #qc.ch(anc-1,6-i+2*7)
#         qc.ry(-np.pi/4,6-i+2*7)
#         qc.cz(ancc,6-i+2*7)
#         qc.ry(np.pi/4,6-i+2*7)
#         if i == 0:
#             qc.cx(ancc,anc)
#         if i == 5:
#             qc.cx(ancc,anc)
#     qc.h(ancc)
#     qc.measure(ancc, state_inj[0])
#     qc.measure(anc, state_inj[1])
#     ##########################################QEC Block#######################################
#     qc.reset(anc), qc.reset(ancc)
#     ##################################Z-Stabilizers##########################################
#     qc.id(anc), qc.id(ancc)
#     qc.h(ancc)
#     qc.cx(0+7*2, anc)
#     qc.cx(2+7*2, anc)
#     qc.cx(4+7*2, anc)
#     qc.cx(6+7*2, anc)
#     qc.h(ancc)

#     qc.id(anc), qc.id(ancc)
#     qc.measure(anc, state_inj[2]), qc.measure(ancc, state_inj[8])
#     qc.reset(anc), qc.reset(ancc)

#     qc.id(anc), qc.id(ancc)
#     qc.h(ancc)
#     qc.cx(ancc, anc)
#     qc.cx(1+7*2, anc)
#     qc.cx(2+7*2, anc)
#     qc.cx(ancc, anc)
#     qc.cx(5+7*2, anc)
#     qc.cx(6+7*2, anc)
#     qc.h(ancc)

#     qc.id(anc), qc.id(ancc)
#     qc.measure(anc, state_inj[3]), qc.measure(ancc, state_inj[9])
#     qc.reset(anc), qc.reset(ancc)

#     qc.id(anc), qc.id(ancc)
#     qc.h(ancc)
#     qc.cx(ancc, anc)
#     qc.cx(3+7*2, anc)
#     qc.cx(4+7*2, anc)
#     qc.cx(ancc, anc)
#     qc.cx(5+7*2, anc)
#     qc.cx(6+7*2, anc)
#     qc.h(ancc)

#     qc.id(anc), qc.id(ancc)
#     qc.measure(anc, state_inj[4]), qc.measure(ancc, state_inj[10])
#     qc.reset(anc), qc.reset(ancc)
#     ##################################X-Stabilizers##############################################
#     qc.id(anc), qc.id(ancc)
#     qc.h(anc)
#     qc.cx(anc, 0+7*2)
#     qc.cx(anc, 2+7*2)
#     qc.cx(anc, 4+7*2)
#     qc.cx(anc, 6+7*2)
#     qc.h(anc)

#     qc.id(anc), qc.id(ancc)
#     qc.measure(anc, state_inj[5]), qc.measure(ancc, state_inj[11])
#     qc.reset(anc), qc.reset(ancc)

#     qc.id(anc), qc.id(ancc)
#     qc.h(anc)
#     qc.cx(anc, ancc)
#     qc.cx(anc, 1+7*2)
#     qc.cx(anc, 2+7*2)
#     qc.cx(anc, ancc)
#     qc.cx(anc, 5+7*2)
#     qc.cx(anc, 6+7*2)
#     qc.h(anc), qc.reset(ancc)

#     qc.id(anc), qc.id(ancc)
#     qc.measure(anc, state_inj[6]), qc.measure(ancc, state_inj[12])
#     qc.reset(anc)

#     qc.id(anc), qc.id(ancc)
#     qc.h(anc)
#     qc.cx(anc, ancc)
#     qc.cx(anc, 3+7*2)
#     qc.cx(anc, 4+7*2)
#     qc.cx(anc, ancc)
#     qc.cx(anc, 5+7*2)
#     qc.cx(anc, 6+7*2)
#     qc.h(anc)

#     qc.id(anc), qc.id(ancc)
#     qc.measure(anc, state_inj[7]), qc.measure(ancc, state_inj[13])
    
#     ########################Controlled-Y Gate####################################################
#     adj_S_L(qc, pos)
#     for i in range(7):
#         qc.cx(i+7*2,i+7*pos)
#     S_L(qc, pos)
#     # read = ClassicalRegister(7)
#     # qc.add_register(read)
#     qc.reset(anc-1)
#     #############################Measure logical state for state injection#############################
#     adj_S_L(qc, pos=2)
#     H_L(qc, pos=2)
#     for i in range(7):
#         qc.cx(i+2*7, anc-1)
#     qc.measure(anc-1,0)
#     #################################Apply conditioned Ry(pi/2) onto the Target###########################
#     for i in range(7):
#         with qc.if_test((0,1)):
#             qc.h(i+7*pos)
#     for i in range(3):
#         with qc.if_test((0,1)):
#             qc.x(i+7*pos)

# def adj_Ty_L(qc: QuantumCircuit, pos: int):
#     state_inj = ClassicalRegister(14)
#     qc.add_register(state_inj)

#     anc = qc.num_qubits - 1
#     ancc = anc - 1

#     for i in range(7):
#         qc.reset(i+7*2)

#     for i in range(7):                        #start noise
#         qc.id(i+7*2)

#     qc.h(0+7*2)
#     qc.h(1+7*2)
#     qc.ry(np.pi/4,2+7*2)
#     qc.h(3+7*2)

#     qc.cx(2+7*2,4+7*2)
#     qc.cx(0+7*2,6+7*2)

#     qc.cx(3+7*2,5+7*2)

#     qc.cx(2+7*2,5+7*2)

#     qc.cx(0+7*2,4+7*2)
#     qc.cx(1+7*2,6+7*2)

#     qc.cx(0+7*2,2+7*2)

#     qc.cx(1+7*2,5+7*2)

#     qc.cx(1+7*2,2+7*2)
#     qc.cx(3+7*2,4+7*2)
#     qc.cx(3+7*2,6+7*2)
#     #################################Controlled Hadamards##########################################
#     qc.reset(anc), qc.reset(anc-1)
#     qc.h(ancc)
#     for i in range(7):
#         #qc.ch(anc-1,6-i+2*7)
#         qc.ry(-np.pi/4,6-i+2*7)
#         qc.cz(ancc,6-i+2*7)
#         qc.ry(np.pi/4,6-i+2*7)
#         if i == 0:
#             qc.cx(ancc,anc)
#         if i == 5:
#             qc.cx(ancc,anc)
#     qc.h(ancc)
#     qc.measure(ancc, state_inj[0])
#     qc.measure(anc, state_inj[1])
#     ##########################################QEC Block#######################################
#     qc.reset(anc), qc.reset(ancc)
#     ##################################Z-Stabilizers##########################################
#     qc.id(anc), qc.id(ancc)
#     qc.h(ancc)
#     qc.cx(0+7*2, anc)
#     qc.cx(2+7*2, anc)
#     qc.cx(4+7*2, anc)
#     qc.cx(6+7*2, anc)
#     qc.h(ancc)

#     qc.id(anc), qc.id(ancc)
#     qc.measure(anc, state_inj[2]), qc.measure(ancc, state_inj[8])
#     qc.reset(anc), qc.reset(ancc)

#     qc.id(anc), qc.id(ancc)
#     qc.h(ancc)
#     qc.cx(ancc, anc)
#     qc.cx(1+7*2, anc)
#     qc.cx(2+7*2, anc)
#     qc.cx(ancc, anc)
#     qc.cx(5+7*2, anc)
#     qc.cx(6+7*2, anc)
#     qc.h(ancc)

#     qc.id(anc), qc.id(ancc)
#     qc.measure(anc, state_inj[3]), qc.measure(ancc, state_inj[9])
#     qc.reset(anc), qc.reset(ancc)

#     qc.id(anc), qc.id(ancc)
#     qc.h(ancc)
#     qc.cx(ancc, anc)
#     qc.cx(3+7*2, anc)
#     qc.cx(4+7*2, anc)
#     qc.cx(ancc, anc)
#     qc.cx(5+7*2, anc)
#     qc.cx(6+7*2, anc)
#     qc.h(ancc)

#     qc.id(anc), qc.id(ancc)
#     qc.measure(anc, state_inj[4]), qc.measure(ancc, state_inj[10])
#     qc.reset(anc), qc.reset(ancc)
#     ##################################X-Stabilizers##############################################
#     qc.id(anc), qc.id(ancc)
#     qc.h(anc)
#     qc.cx(anc, 0+7*2)
#     qc.cx(anc, 2+7*2)
#     qc.cx(anc, 4+7*2)
#     qc.cx(anc, 6+7*2)
#     qc.h(anc)

#     qc.id(anc), qc.id(ancc)
#     qc.measure(anc, state_inj[5]), qc.measure(ancc, state_inj[11])
#     qc.reset(anc), qc.reset(ancc)

#     qc.id(anc), qc.id(ancc)
#     qc.h(anc)
#     qc.cx(anc, ancc)
#     qc.cx(anc, 1+7*2)
#     qc.cx(anc, 2+7*2)
#     qc.cx(anc, ancc)
#     qc.cx(anc, 5+7*2)
#     qc.cx(anc, 6+7*2)
#     qc.h(anc), qc.reset(ancc)

#     qc.id(anc), qc.id(ancc)
#     qc.measure(anc, state_inj[6]), qc.measure(ancc, state_inj[12])
#     qc.reset(anc)

#     qc.id(anc), qc.id(ancc)
#     qc.h(anc)
#     qc.cx(anc, ancc)
#     qc.cx(anc, 3+7*2)
#     qc.cx(anc, 4+7*2)
#     qc.cx(anc, ancc)
#     qc.cx(anc, 5+7*2)
#     qc.cx(anc, 6+7*2)
#     qc.h(anc)

#     qc.id(anc), qc.id(ancc)
#     qc.measure(anc, state_inj[7]), qc.measure(ancc, state_inj[13])
    
#     ########################Controlled-Y Gate####################################################
#     adj_S_L(qc, pos)
#     for i in range(7):
#         qc.cx(i+7*2,i+7*pos)
#     S_L(qc, pos)
#     # read = ClassicalRegister(7)
#     # qc.add_register(read)
#     qc.reset(anc-1)
#     #############################Measure logical state for state injection#############################
#     adj_S_L(qc, pos=2)
#     H_L(qc, pos=2)
#     for i in range(7):
#         qc.cx(i+2*7, anc-1)
#     qc.measure(anc-1,0)
#     #################################Apply conditioned Ry(-pi/2) onto the Target###########################
#     for i in range(3):
#         with qc.if_test((0,0)):
#             qc.x(i+7*pos)
#     for i in range(7):
#         with qc.if_test((0,0)):
#             qc.h(i+7*pos)

def Ty_ec_L(qc: QuantumCircuit, pos: int):
    state_inj = ClassicalRegister(8)
    qc.add_register(state_inj)
    flags = ClassicalRegister(6)
    qc.add_register(flags)

    anc = qc.num_qubits - 1
    ancc = anc - 1

    for i in range(7):
        qc.reset(i+7*2)

    for i in range(7):                        #start noise
        qc.id(i+7*2)

    qc.h(0+7*2)
    qc.h(1+7*2)
    qc.ry(np.pi/4,2+7*2)
    qc.h(3+7*2)

    qc.cx(2+7*2,4+7*2)
    qc.cx(0+7*2,6+7*2)

    qc.cx(3+7*2,5+7*2)

    qc.cx(2+7*2,5+7*2)

    qc.cx(0+7*2,4+7*2)
    qc.cx(1+7*2,6+7*2)

    qc.cx(0+7*2,2+7*2)

    qc.cx(1+7*2,5+7*2)

    qc.cx(1+7*2,2+7*2)
    qc.cx(3+7*2,4+7*2)
    qc.cx(3+7*2,6+7*2)
    #################################Controlled Hadamards##########################################
    qc.reset(anc), qc.reset(anc-1)
    qc.h(anc-1)
    for i in range(7):
        #qc.ch(anc-1,6-i+2*7)
        qc.ry(-np.pi/4,6-i+2*7)
        qc.cz(anc-1,6-i+2*7)
        qc.ry(np.pi/4,6-i+2*7)
        if i == 0:
            qc.cx(anc-1,anc)
        if i == 5:
            qc.cx(anc-1,anc)
    qc.h(anc-1)
    qc.measure(anc-1, state_inj[0])
    qc.measure(anc, state_inj[1])
    ##########################################QEC Block#######################################
    qc.reset(anc), qc.reset(ancc)
    ##################################Z-Stabilizers##########################################
    qc.id(anc), qc.id(ancc)
    qc.h(ancc)
    qc.cx(0+7*2, anc)
    qc.cx(ancc, anc)
    qc.cx(2+7*2, anc)
    qc.cx(4+7*2, anc)
    qc.cx(ancc, anc)
    qc.cx(6+7*2, anc)

    qc.id(anc), qc.h(ancc), qc.id(ancc)
    qc.measure(anc, state_inj[4]), qc.measure(ancc, flags[0])
    qc.reset(anc), qc.reset(ancc)

    qc.id(anc), qc.id(ancc)
    qc.h(ancc)
    qc.cx(1+7*2, anc)
    qc.cx(ancc, anc)
    qc.cx(2+7*2, anc)
    qc.cx(5+7*2, anc)
    qc.cx(ancc, anc)
    qc.cx(6+7*2, anc)

    qc.id(anc), qc.h(ancc), qc.id(ancc)
    qc.measure(anc, state_inj[3]), qc.measure(ancc, flags[1])
    qc.reset(anc), qc.reset(ancc)

    qc.id(anc), qc.id(ancc)
    qc.h(ancc)
    qc.cx(3+7*2, anc)
    qc.cx(ancc, anc)
    qc.cx(4+7*2, anc)
    qc.cx(5+7*2, anc)
    qc.cx(ancc, anc)
    qc.cx(6+7*2, anc)

    qc.id(anc), qc.h(ancc), qc.id(ancc)
    qc.measure(anc, state_inj[2]), qc.measure(ancc, flags[2])
    qc.reset(anc), qc.reset(ancc)
    ##################################X-Stabilizers##############################################
    qc.id(anc), qc.id(ancc)
    qc.h(anc)
    qc.cx(anc, 0+7*2)
    qc.cx(anc, ancc)
    qc.cx(anc, 2+7*2)
    qc.cx(anc, 4+7*2)
    qc.cx(anc, ancc)
    qc.cx(anc, 6+7*2)
    qc.h(anc)

    qc.id(anc), qc.id(ancc)
    qc.measure(anc, state_inj[7]), qc.measure(ancc, flags[3])
    qc.reset(anc), qc.reset(ancc)

    qc.id(anc), qc.id(ancc)
    qc.h(anc)
    qc.cx(anc, 1+7*2)
    qc.cx(anc, ancc)
    qc.cx(anc, 2+7*2)
    qc.cx(anc, 5+7*2)
    qc.cx(anc, ancc)
    qc.cx(anc, 6+7*2)
    qc.h(anc)

    qc.id(anc), qc.id(ancc)
    qc.measure(anc, state_inj[6]), qc.measure(ancc, flags[4])
    qc.reset(anc), qc.reset(ancc)

    qc.id(anc), qc.id(ancc)
    qc.h(anc)
    qc.cx(anc, 3+7*2)
    qc.cx(anc, ancc)
    qc.cx(anc, 4+7*2)
    qc.cx(anc, 5+7*2)
    qc.cx(anc, ancc)
    qc.cx(anc, 6+7*2)
    qc.h(anc)

    qc.id(anc), qc.id(ancc)
    qc.measure(anc, state_inj[5]), qc.measure(ancc, flags[5])
    ###############################QEC-Block####################################################
    
    with qc.if_test((state_inj[2],0)):             #qbit 0
        with qc.if_test((state_inj[3],0)):
            with qc.if_test((state_inj[4],1)):
                qc.x(0+7*pos)

    with qc.if_test((state_inj[2],0)):             #qbit 1
        with qc.if_test((state_inj[3],1)):
            with qc.if_test((state_inj[4],0)):
                qc.x(1+7*pos)
    
    with qc.if_test((state_inj[2],0)):             #qbit 2
        with qc.if_test((state_inj[3],1)):
            with qc.if_test((state_inj[4],1)):
                qc.x(2+7*pos)
    
    with qc.if_test((state_inj[2],1)):             #qbit 3
        with qc.if_test((state_inj[3],0)):
            with qc.if_test((state_inj[4],0)):
                qc.x(3+7*pos)
    
    with qc.if_test((state_inj[2],1)):             #qbit 4
        with qc.if_test((state_inj[3],0)):
            with qc.if_test((state_inj[4],1)):
                qc.x(4+7*pos)
    
    with qc.if_test((state_inj[2],1)):             #qbit 5
        with qc.if_test((state_inj[3],1)):
            with qc.if_test((state_inj[4],0)):
                qc.x(5+7*pos)
    
    with qc.if_test((state_inj[2],1)):             #qbit 6
        with qc.if_test((state_inj[3],1)):
            with qc.if_test((state_inj[4],1)):
                qc.x(6+7*pos)

    with qc.if_test((state_inj[5],0)):             #qbit 0
        with qc.if_test((state_inj[6],0)):
            with qc.if_test((state_inj[7],1)):
                qc.z(0+7*pos)

    with qc.if_test((state_inj[5],0)):             #qbit 1
        with qc.if_test((state_inj[6],1)):
            with qc.if_test((state_inj[7],0)):
                qc.z(1+7*pos)
    
    with qc.if_test((state_inj[5],0)):             #qbit 2
        with qc.if_test((state_inj[6],1)):
            with qc.if_test((state_inj[7],1)):
                qc.z(2+7*pos)
    
    with qc.if_test((state_inj[5],1)):             #qbit 3
        with qc.if_test((state_inj[6],0)):
            with qc.if_test((state_inj[7],0)):
                qc.z(3+7*pos)
    
    with qc.if_test((state_inj[5],1)):             #qbit 4
        with qc.if_test((state_inj[6],0)):
            with qc.if_test((state_inj[7],1)):
                qc.z(4+7*pos)
    
    with qc.if_test((state_inj[5],1)):             #qbit 5
        with qc.if_test((state_inj[6],1)):
            with qc.if_test((state_inj[7],0)):
                qc.z(5+7*pos)
    
    with qc.if_test((state_inj[5],1)):             #qbit 6
        with qc.if_test((state_inj[6],1)):
            with qc.if_test((state_inj[7],1)):
                qc.z(6+7*pos)
    qc.reset(anc)
    for i in range(8):
        qc.measure(anc, state_inj[i])
    ########################Controlled-Y Gate####################################################
    adj_S_L(qc, pos)
    for i in range(7):
        qc.cx(i+7*2,i+7*pos)
    S_L(qc, pos)
    # read = ClassicalRegister(7)
    # qc.add_register(read)
    qc.reset(anc-1)
    #############################Measure logical state for state injection#############################
    adj_S_L(qc, pos=2)
    H_L(qc, pos=2)
    for i in range(7):
        qc.cx(i+2*7, anc-1)
    qc.measure(anc-1,0)
    #################################Apply conditioned Ry(pi/2) onto the Target###########################
    for i in range(7):
        with qc.if_test((0,1)):
            qc.h(i+7*pos)
    for i in range(3):
        with qc.if_test((0,1)):
            qc.x(i+7*pos)

def T_L(qc: QuantumCircuit, pos: int, qecc, err = False , ecc = False):
    H_L(qc, pos=pos)
    adj_S_L(qc, pos=pos)
    H_L(qc, pos=pos)
    if ecc:
        Ty_ec_L(qc, pos=pos)
    else:
        Ty_L(qc, pos=pos)
    H_L(qc, pos=pos)
    S_L(qc, pos=pos)
    H_L(qc, pos=pos)
    if err:
        qec_ft(qc, qecc=qecc, pos=pos)

def adj_T_L(qc: QuantumCircuit, pos: int, qecc, err = False, ecc = False):
    H_L(qc, pos=pos)
    adj_S_L(qc, pos=pos)
    H_L(qc, pos=pos)
    if ecc:
        adj_Ty_ec_L(qc, pos=pos)
    else:
        adj_Ty_L(qc, pos=pos)
    H_L(qc, pos=pos)
    S_L(qc, pos=pos)
    H_L(qc, pos=pos)
    if err:
        qec_ft(qc, qecc=qecc, pos=pos)

def adj_Ty_ec_L(qc: QuantumCircuit, pos: int):
    state_inj = ClassicalRegister(8)
    qc.add_register(state_inj)
    flags = ClassicalRegister(6)
    qc.add_register(flags)

    anc = qc.num_qubits - 1
    ancc = anc - 1

    for i in range(7):
        qc.reset(i+7*2)

    for i in range(7):                        #start noise
        qc.id(i+7*2)

    qc.h(0+7*2)
    qc.h(1+7*2)
    qc.ry(np.pi/4,2+7*2)
    qc.h(3+7*2)

    qc.cx(2+7*2,4+7*2)
    qc.cx(0+7*2,6+7*2)

    qc.cx(3+7*2,5+7*2)

    qc.cx(2+7*2,5+7*2)

    qc.cx(0+7*2,4+7*2)
    qc.cx(1+7*2,6+7*2)

    qc.cx(0+7*2,2+7*2)

    qc.cx(1+7*2,5+7*2)

    qc.cx(1+7*2,2+7*2)
    qc.cx(3+7*2,4+7*2)
    qc.cx(3+7*2,6+7*2)
    #################################Controlled Hadamards##########################################
    qc.reset(anc), qc.reset(anc-1)
    qc.h(anc-1)
    for i in range(7):
        #qc.ch(anc-1,6-i+2*7)
        qc.ry(-np.pi/4,6-i+2*7)
        qc.cz(anc-1,6-i+2*7)
        qc.ry(np.pi/4,6-i+2*7)
        if i == 0:
            qc.cx(anc-1,anc)
        if i == 5:
            qc.cx(anc-1,anc)
    qc.h(anc-1)
    qc.measure(anc-1, state_inj[0])
    qc.measure(anc, state_inj[1])
    ##########################################QEC Block#######################################
    qc.reset(anc), qc.reset(ancc)
    ##################################Z-Stabilizers##########################################
    qc.id(anc), qc.id(ancc)
    qc.h(ancc)
    qc.cx(0+7*2, anc)
    qc.cx(ancc, anc)
    qc.cx(2+7*2, anc)
    qc.cx(4+7*2, anc)
    qc.cx(ancc, anc)
    qc.cx(6+7*2, anc)

    qc.id(anc), qc.h(ancc), qc.id(ancc)
    qc.measure(anc, state_inj[4]), qc.measure(ancc, flags[0])
    qc.reset(anc), qc.reset(ancc)

    qc.id(anc), qc.id(ancc)
    qc.h(ancc)
    qc.cx(1+7*2, anc)
    qc.cx(ancc, anc)
    qc.cx(2+7*2, anc)
    qc.cx(5+7*2, anc)
    qc.cx(ancc, anc)
    qc.cx(6+7*2, anc)

    qc.id(anc), qc.h(ancc), qc.id(ancc)
    qc.measure(anc, state_inj[3]), qc.measure(ancc, flags[1])
    qc.reset(anc), qc.reset(ancc)

    qc.id(anc), qc.id(ancc)
    qc.h(ancc)
    qc.cx(3+7*2, anc)
    qc.cx(ancc, anc)
    qc.cx(4+7*2, anc)
    qc.cx(5+7*2, anc)
    qc.cx(ancc, anc)
    qc.cx(6+7*2, anc)

    qc.id(anc), qc.h(ancc), qc.id(ancc)
    qc.measure(anc, state_inj[2]), qc.measure(ancc, flags[2])
    qc.reset(anc), qc.reset(ancc)
    ##################################X-Stabilizers##############################################
    qc.id(anc), qc.id(ancc)
    qc.h(anc)
    qc.cx(anc, 0+7*2)
    qc.cx(anc, ancc)
    qc.cx(anc, 2+7*2)
    qc.cx(anc, 4+7*2)
    qc.cx(anc, ancc)
    qc.cx(anc, 6+7*2)
    qc.h(anc)

    qc.id(anc), qc.id(ancc)
    qc.measure(anc, state_inj[7]), qc.measure(ancc, flags[3])
    qc.reset(anc), qc.reset(ancc)

    qc.id(anc), qc.id(ancc)
    qc.h(anc)
    qc.cx(anc, 1+7*2)
    qc.cx(anc, ancc)
    qc.cx(anc, 2+7*2)
    qc.cx(anc, 5+7*2)
    qc.cx(anc, ancc)
    qc.cx(anc, 6+7*2)
    qc.h(anc)

    qc.id(anc), qc.id(ancc)
    qc.measure(anc, state_inj[6]), qc.measure(ancc, flags[4])
    qc.reset(anc), qc.reset(ancc)

    qc.id(anc), qc.id(ancc)
    qc.h(anc)
    qc.cx(anc, 3+7*2)
    qc.cx(anc, ancc)
    qc.cx(anc, 4+7*2)
    qc.cx(anc, 5+7*2)
    qc.cx(anc, ancc)
    qc.cx(anc, 6+7*2)
    qc.h(anc)

    qc.id(anc), qc.id(ancc)
    qc.measure(anc, state_inj[5]), qc.measure(ancc, flags[5])
    ###############################QEC-Block####################################################
    
    with qc.if_test((state_inj[2],0)):             #qbit 0
        with qc.if_test((state_inj[3],0)):
            with qc.if_test((state_inj[4],1)):
                qc.x(0+7*pos)

    with qc.if_test((state_inj[2],0)):             #qbit 1
        with qc.if_test((state_inj[3],1)):
            with qc.if_test((state_inj[4],0)):
                qc.x(1+7*pos)
    
    with qc.if_test((state_inj[2],0)):             #qbit 2
        with qc.if_test((state_inj[3],1)):
            with qc.if_test((state_inj[4],1)):
                qc.x(2+7*pos)
    
    with qc.if_test((state_inj[2],1)):             #qbit 3
        with qc.if_test((state_inj[3],0)):
            with qc.if_test((state_inj[4],0)):
                qc.x(3+7*pos)
    
    with qc.if_test((state_inj[2],1)):             #qbit 4
        with qc.if_test((state_inj[3],0)):
            with qc.if_test((state_inj[4],1)):
                qc.x(4+7*pos)
    
    with qc.if_test((state_inj[2],1)):             #qbit 5
        with qc.if_test((state_inj[3],1)):
            with qc.if_test((state_inj[4],0)):
                qc.x(5+7*pos)
    
    with qc.if_test((state_inj[2],1)):             #qbit 6
        with qc.if_test((state_inj[3],1)):
            with qc.if_test((state_inj[4],1)):
                qc.x(6+7*pos)

    with qc.if_test((state_inj[5],0)):             #qbit 0
        with qc.if_test((state_inj[6],0)):
            with qc.if_test((state_inj[7],1)):
                qc.z(0+7*pos)

    with qc.if_test((state_inj[5],0)):             #qbit 1
        with qc.if_test((state_inj[6],1)):
            with qc.if_test((state_inj[7],0)):
                qc.z(1+7*pos)
    
    with qc.if_test((state_inj[5],0)):             #qbit 2
        with qc.if_test((state_inj[6],1)):
            with qc.if_test((state_inj[7],1)):
                qc.z(2+7*pos)
    
    with qc.if_test((state_inj[5],1)):             #qbit 3
        with qc.if_test((state_inj[6],0)):
            with qc.if_test((state_inj[7],0)):
                qc.z(3+7*pos)
    
    with qc.if_test((state_inj[5],1)):             #qbit 4
        with qc.if_test((state_inj[6],0)):
            with qc.if_test((state_inj[7],1)):
                qc.z(4+7*pos)
    
    with qc.if_test((state_inj[5],1)):             #qbit 5
        with qc.if_test((state_inj[6],1)):
            with qc.if_test((state_inj[7],0)):
                qc.z(5+7*pos)
    
    with qc.if_test((state_inj[5],1)):             #qbit 6
        with qc.if_test((state_inj[6],1)):
            with qc.if_test((state_inj[7],1)):
                qc.z(6+7*pos)
    qc.reset(anc)
    for i in range(8):
        qc.measure(anc, state_inj[i])
    ########################Controlled-Y Gate####################################################
    adj_S_L(qc, pos)
    for i in range(7):
        qc.cx(i+7*2,i+7*pos)
    S_L(qc, pos)
    # read = ClassicalRegister(7)
    # qc.add_register(read)
    qc.reset(anc-1)
    #############################Measure logical state for state injection#############################
    adj_S_L(qc, pos=2)
    H_L(qc, pos=2)
    for i in range(7):
        qc.cx(i+2*7, anc-1)
    qc.measure(anc-1,0)
    #################################Apply conditioned Ry(pi/2) onto the Target###########################
    for i in range(3):
        with qc.if_test((0,0)):
            qc.x(i+7*pos)
    for i in range(7):
        with qc.if_test((0,0)):
            qc.h(i+7*pos)
##############################################################
#Ideal magic state, das spart viele shots und ist ja auch das was man in der Praxis macht gell
def Ty_L(qc: QuantumCircuit, pos: int):
    state_inj = ClassicalRegister(1)
    qc.add_register(state_inj)

    anc = qc.num_qubits - 1
    ancc = anc - 1

    for i in range(7):
        qc.reset(i+7*2)

    qc.append(h_ideal,[0+7*2])
    qc.append(h_ideal,[1+7*2])
    qc.ry(np.pi/4,2+7*2)
    qc.append(h_ideal,[3+7*2])

    qc.append(cx_ideal, [4+7*2, 2+7*2])
    qc.append(cx_ideal, [6+7*2, 0+7*2])
    
    qc.append(cx_ideal, [5+7*2, 3+7*2])

    qc.append(cx_ideal, [5+7*2, 2+7*2])

    qc.append(cx_ideal, [4+7*2, 0+7*2])
    qc.append(cx_ideal, [6+7*2, 1+7*2])

    qc.append(cx_ideal, [2+7*2, 0+7*2])

    qc.append(cx_ideal, [5+7*2, 1+7*2])

    qc.append(cx_ideal, [2+7*2, 1+7*2])
    qc.append(cx_ideal, [4+7*2, 3+7*2])
    qc.append(cx_ideal, [6+7*2, 3+7*2])
    #################################Controlled Hadamards##########################################
    qc.reset(ancc)
    qc.append(h_ideal,[ancc])
    for i in range(7):
        #qc.ch(anc-1,6-i+2*7)
        qc.ry(-np.pi/4,6-i+2*7)
        qc.cz(ancc,6-i+2*7)
        qc.ry(np.pi/4,6-i+2*7)
    qc.append(h_ideal,[ancc])
    qc.measure(ancc, state_inj[0])
    ########################Controlled-Y Gate####################################################
    adj_S_L(qc, pos)
    for i in range(7):
        qc.cx(i+7*2,i+7*pos)
    S_L(qc, pos)
    # read = ClassicalRegister(7)
    # qc.add_register(read)
    qc.reset(anc-1)
    #############################Measure logical state for state injection#############################
    adj_S_L(qc, pos=2)
    H_L(qc, pos=2)
    for i in range(7):
        qc.cx(i+2*7, anc-1)
    qc.measure(anc-1,0)
    #################################Apply conditioned Ry(pi/2) onto the Target###########################
    for i in range(7):
        with qc.if_test((0,1)):
            qc.h(i+7*pos)
    for i in range(3):
        with qc.if_test((0,1)):
            qc.x(i+7*pos)

def adj_Ty_L(qc: QuantumCircuit, pos: int):
    state_inj = ClassicalRegister(1)
    qc.add_register(state_inj)

    anc = qc.num_qubits - 1
    ancc = anc - 1

    for i in range(7):
        qc.reset(i+7*2)

    qc.append(h_ideal,[0+7*2])
    qc.append(h_ideal,[1+7*2])
    qc.ry(np.pi/4,2+7*2)
    qc.append(h_ideal,[3+7*2])

    qc.append(cx_ideal, [4+7*2, 2+7*2])
    qc.append(cx_ideal, [6+7*2, 0+7*2])
    
    qc.append(cx_ideal, [5+7*2, 3+7*2])

    qc.append(cx_ideal, [5+7*2, 2+7*2])

    qc.append(cx_ideal, [4+7*2, 0+7*2])
    qc.append(cx_ideal, [6+7*2, 1+7*2])

    qc.append(cx_ideal, [2+7*2, 0+7*2])

    qc.append(cx_ideal, [5+7*2, 1+7*2])

    qc.append(cx_ideal, [2+7*2, 1+7*2])
    qc.append(cx_ideal, [4+7*2, 3+7*2])
    qc.append(cx_ideal, [6+7*2, 3+7*2])
    #################################Controlled Hadamards##########################################
    qc.reset(ancc)
    qc.h(ancc)
    for i in range(7):
        #qc.ch(anc-1,6-i+2*7)
        qc.ry(-np.pi/4,6-i+2*7)
        qc.cz(ancc,6-i+2*7)
        qc.ry(np.pi/4,6-i+2*7)
    qc.h(ancc)
    qc.measure(ancc, state_inj[0])
    ########################Controlled-Y Gate####################################################
    adj_S_L(qc, pos)
    for i in range(7):
        qc.cx(i+7*2,i+7*pos)
    S_L(qc, pos)
    # read = ClassicalRegister(7)
    # qc.add_register(read)
    qc.reset(anc-1)
    #############################Measure logical state for state injection#############################
    adj_S_L(qc, pos=2)
    H_L(qc, pos=2)
    for i in range(7):
        qc.cx(i+2*7, anc-1)
    qc.measure(anc-1,0)
    #################################Apply conditioned Ry(-pi/2) onto the Target###########################
    for i in range(3):
        with qc.if_test((0,0)):
            qc.x(i+7*pos)
    for i in range(7):
        with qc.if_test((0,0)):
            qc.h(i+7*pos)
##############################################################
circ = QuantumCircuit(1)
circ.rz(np.pi/8, 0)
basis = ["t", "tdg", "z", "h"]
approx = generate_basic_approximations(basis, depth=3)
skd = SolovayKitaev(recursion_degree=2, basic_approximations=approx)
rootT = skd(circ)

def root_T_L(qc: QuantumCircuit, pos: int, qecc, err = False, ecc = False):
    instruction = rootT.data
    counter = 0
    for i in instruction:
        if counter%3 == 0:
            if err:
                qec_ft(qc, qecc=qecc, pos=pos)
        if i.name == "t":
            T_L(qc, pos=pos, qecc=qecc, err=err, ecc=ecc)
            counter += 1
        if i.name == "tdg":
            adj_T_L(qc, pos=pos, qecc=qecc, err=err, ecc=ecc)
            counter += 1
        if i.name == "h":
            H_L(qc, pos=pos)

circ = QuantumCircuit(1)
circ.rz(-np.pi/8, 0)
basis = ["t", "tdg", "z", "h"]
approx = generate_basic_approximations(basis, depth=3)
skd = SolovayKitaev(recursion_degree=2, basic_approximations=approx)
adj_rootT = skd(circ)

def adj_root_T_L(qc: QuantumCircuit, pos: int, qecc, err=False, ecc = False):
    instruction = adj_rootT.data
    counter = 0
    for i in instruction:
        if counter%3 == 0:
            if err:
                qec_ft(qc, qecc=qecc, pos=pos)
        if i.name == "t":
            T_L(qc, pos=pos, qecc=qecc, err=err, ecc=ecc)
            counter += 1
        if i.name == "tdg":
            adj_T_L(qc, pos=pos, qecc=qecc, err=err, ecc=ecc)
            counter += 1
        if i.name == "h":
            H_L(qc, pos=pos)

def CT_L(qc: QuantumCircuit, qecc, err=False, ecc = False):
    root_T_L(qc, 0, qecc = qecc, err=err, ecc=ecc)
    root_T_L(qc, 1, qecc = qecc, err=err, ecc=ecc)
    CNOT_L(qc, 0)
    adj_root_T_L(qc, 1, qecc = qecc, err=err, ecc=ecc)
    CNOT_L(qc, 0)
#######################################################################

# def U2(qc: QuantumCircuit, pos: int, qecc, ecc=False, err=False):
#     gate = ['s', 'h', 'tdg', 'h', 'tdg', 'h', 't', 'h', 't', 'h', 't', 'h', 't', 'h', 'tdg', 'h', 'sdg', 'h', 'tdg', 'h', 't', 'h', 't', 'h', 't', 'h', 't', 'h', 'tdg', 'h', 'tdg', 'h', 'sdg', 'h', 'tdg', 'h', 't', 'h', 't', 'h', 't', 'h', 't', 'h', 's', 'h', 't', 'h', 'tdg', 'h', 'sdg', 'h', 'sdg', 'tdg', 'h', 'tdg', 'h', 'tdg', 'h', 't', 'h', 't', 'h', 'tdg', 'h', 't', 'h', 't', 'h', 't', 'h', 'tdg', 'h', 'tdg', 'h', 't', 'h', 's', 'h', 's', 'h', 't', 'h', 'tdg', 'h', 'sdg', 'h', 'tdg', 'h', 'tdg', 'h', 'tdg', 'h', 'tdg', 'h', 't', 'h', 't', 'h', 'tdg', 'h', 't', 'h', 't', 'h', 'tdg', 'h', 'tdg', 'h', 'tdg', 'h', 't', 'h', 'tdg', 'h', 'tdg', 'h', 't', 'h', 't', 'h', 't']
#     for i in gate:
#         if i == "s":
#             S_L(qc, pos=pos)
#         if i == "sdg":
#             adj_S_L(qc, pos=pos)
#         if i == "t":
#             T_L(qc, pos=pos, qecc=qecc, err=err, ecc=ecc)
#             if err:
#                 qec_ft(qc, qecc, pos=pos)
#         if i == "tdg":
#             adj_T_L(qc, pos=pos, qecc=qecc, err=err, ecc=ecc)
#             if err:
#                 qec_ft(qc, qecc, pos=pos)
#         if i == "h":
#             H_L(qc, pos=pos)

# def adjU2(qc: QuantumCircuit, pos: int, qecc, ecc=False, err=False):
#     gate = ['s', 'h', 'tdg', 'h', 'tdg', 'h', 't', 'h', 't', 'h', 't', 'h', 't', 'h', 'tdg', 'h', 'sdg', 'h', 'tdg', 'h', 't', 'h', 't', 'h', 't', 'h', 't', 'h', 'tdg', 'h', 'tdg', 'h', 'sdg', 'h', 'tdg', 'h', 't', 'h', 't', 'h', 't', 'h', 't', 'h', 's', 'h', 't', 'h', 'tdg', 'h', 'sdg', 'h', 'sdg', 'tdg', 'h', 'tdg', 'h', 'tdg', 'h', 't', 'h', 't', 'h', 'tdg', 'h', 't', 'h', 't', 'h', 't', 'h', 'tdg', 'h', 'tdg', 'h', 't', 'h', 's', 'h', 's', 'h', 't', 'h', 'tdg', 'h', 'sdg', 'h', 'tdg', 'h', 'tdg', 'h', 'tdg', 'h', 'tdg', 'h', 't', 'h', 't', 'h', 'tdg', 'h', 't', 'h', 't', 'h', 'tdg', 'h', 'tdg', 'h', 'tdg', 'h', 't', 'h', 'tdg', 'h', 'tdg', 'h', 't', 'h', 't', 'h', 't']
#     gate.reverse()
#     for i in gate:
#         if i == "sdg":
#             S_L(qc, pos=pos)
#         if i == "s":
#             adj_S_L(qc, pos=pos)
#         if i == "tdg":
#             T_L(qc, pos=pos, qecc=qecc, err=err, ecc=ecc)
#             if err:
#                 qec_ft(qc, qecc, pos=pos)
#         if i == "t":
#             adj_T_L(qc, pos=pos, qecc=qecc, err=err, ecc=ecc)
#             if err:
#                 qec_ft(qc, qecc, pos=pos)
#         if i == "h":
#             H_L(qc, pos=pos)

# def CU_L(qc: QuantumCircuit, qecc, ecc=False, err = False):
#     U2(qc, pos=0, qecc=qecc, ecc=False, err=err)
#     U2(qc, pos=1, qecc=qecc, ecc=False, err=err)
#     CNOT_L(qc, 0)
#     adjU2(qc, pos=1, qecc=qecc, ecc=False, err=err)
#     CNOT_L(qc, 0)

################################################################
def convert(bin: str):                  #konvertiert den bitstring in decimal, e.g. 0110 = 0.375
    k = list(bin)
    a = [int(i) for i in k]
    n = 0
    for i in range(len(a)):
        if a[i] == 1:
            n += 1/2**(i+1)
    return n

def U2(qc: QuantumCircuit, pos: int, gate: list, qecc, err=False, ecc=False):
    for i in gate:
        if i == "s":
            S_L(qc, pos=pos)
        if i == "sdg":
            adj_S_L(qc, pos=pos)
        if i == "t":
            T_L(qc, pos=pos, qecc=qecc, err=err , ecc=ecc)
        if i == "tdg":
            adj_T_L(qc, pos=pos, qecc=qecc, err=err , ecc=ecc)
        if i == "h":
            H_L(qc, pos=pos)
        if i == "z":
            Z_L(qc, pos=pos)

def CU_L(qc: QuantumCircuit, Ugates: list, adjUgates: list, qecc, err=False, ecc=False):
    U2(qc, 0, Ugates, qecc, err=err, ecc=ecc)
    U2(qc, 1, Ugates, qecc, err=err, ecc=ecc)
    CNOT_L(qc, control=0)
    U2(qc, 1, adjUgates, qecc, err=err, ecc=ecc)
    CNOT_L(qc, control=0)

def Leon(iter: int, n:int, argh: float, err = False, ecc = False, k = 1):       #each iteration own circuit
    angle = np.linspace(0,1,n+2)
    angle = np.delete(angle, [n+1])
    angle = np.delete(angle, [0])

    a, b = [], []
    with open("Unitary/unitary{}.txt".format(n), "r") as file:
        for line in file:
            a.append(list(map(str, line.strip().split(","))))
    with open("Unitary/adjunitary{}.txt".format(n), "r") as file:
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
                    qecc = ClassicalRegister(6)
                    qc.add_register(qecc)

                    X_L(qc,1)
                    H_L(qc,0)
                    #############################
                    for j in range(2**(iter-t-1)):
                        CU_L(qc, a[o], b[o], qecc, err=err, ecc=ecc)
                    ###############################
                    for l in rots:
                        if l == 0.25:
                            adj_S_L(qc, pos=0)
                        if l == 0.125:
                            adj_T_L(qc, pos=0, qecc=qecc, err=err, ecc=ecc)
                    H_L(qc, pos=0)
                    if err:
                        qec_ft(qc, qecc=qecc, pos=0)
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

    #print(counts)

    bitstring = list(counts.keys())
    bitstring = [i.replace(" ","") for i in bitstring]


    hmm = list(counts.values())

    allcbits = len(bitstring[0])                
    pre, preselected = [i[allcbits-3:allcbits-1] for i in bitstring], 0
    bits = [i[:7] for i in bitstring]
    postprocess = [i[7:allcbits-10] for i in bitstring]

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
            bits[i] = "post"

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
    #magic = (magic/shots)

    # print("0: ", zeros*100, "%")
    # print("1: ", ones*100, "%")
    # print("Preselection discarded: ", (preselected/shots)*100, "%")
    # print("Postselection discarded: ", (post/shots)*100, "%")
    return zeros, ones, preselected, post#,magic

def qec(qc: QuantumCircuit, qecc, pos: int):
    anc = qc.num_qubits - 1
    qc.reset(anc)
    ##################################Z-Stabilizers##########################################
    qc.cx(0+7*pos, anc)
    qc.cx(2+7*pos, anc)
    qc.cx(4+7*pos, anc)
    qc.cx(6+7*pos, anc)

    qc.id(anc)
    qc.measure(anc, qecc[2])
    qc.reset(anc)
    qc.id(anc)

    qc.cx(1+7*pos, anc)
    qc.cx(2+7*pos, anc)
    qc.cx(5+7*pos, anc)
    qc.cx(6+7*pos, anc)

    qc.id(anc)
    qc.measure(anc, qecc[1])
    qc.reset(anc)
    qc.id(anc)

    qc.cx(3+7*pos, anc)
    qc.cx(4+7*pos, anc)
    qc.cx(5+7*pos, anc)
    qc.cx(6+7*pos, anc)

    qc.id(anc)
    qc.measure(anc, qecc[0])
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
    qc.measure(anc, qecc[5])
    qc.reset(anc)
    qc.id(anc)

    qc.h(anc)
    qc.cx(anc, 1+7*pos)
    qc.cx(anc, 2+7*pos)
    qc.cx(anc, 5+7*pos)
    qc.cx(anc, 6+7*pos)
    qc.h(anc)

    qc.id(anc)
    qc.measure(anc, qecc[4])
    qc.reset(anc)
    qc.id(anc)

    qc.h(anc)
    qc.cx(anc, 3+7*pos)
    qc.cx(anc, 4+7*pos)
    qc.cx(anc, 5+7*pos)
    qc.cx(anc, 6+7*pos)
    qc.h(anc)

    qc.id(anc)
    qc.measure(anc, qecc[3])
    qc.reset(anc)
    ##################################Bitflip Error correction##############################################
    
    with qc.if_test((qecc[0],0)):             #qbit 0
        with qc.if_test((qecc[1],0)):
            with qc.if_test((qecc[2],1)):
                qc.x(0+7*pos)

    with qc.if_test((qecc[0],0)):             #qbit 1
        with qc.if_test((qecc[1],1)):
            with qc.if_test((qecc[2],0)):
                qc.x(1+7*pos)
    
    with qc.if_test((qecc[0],0)):             #qbit 2
        with qc.if_test((qecc[1],1)):
            with qc.if_test((qecc[2],1)):
                qc.x(2+7*pos)
    
    with qc.if_test((qecc[0],1)):             #qbit 3
        with qc.if_test((qecc[1],0)):
            with qc.if_test((qecc[2],0)):
                qc.x(3+7*pos)
    
    with qc.if_test((qecc[0],1)):             #qbit 4
        with qc.if_test((qecc[1],0)):
            with qc.if_test((qecc[2],1)):
                qc.x(4+7*pos)
    
    with qc.if_test((qecc[0],1)):             #qbit 5
        with qc.if_test((qecc[1],1)):
            with qc.if_test((qecc[2],0)):
                qc.x(5+7*pos)
    
    with qc.if_test((qecc[0],1)):             #qbit 6
        with qc.if_test((qecc[1],1)):
            with qc.if_test((qecc[2],1)):
                qc.x(6+7*pos)

    ##################################Phaseflip Error correction##############################################
    
    with qc.if_test((qecc[3],0)):             #qbit 0
        with qc.if_test((qecc[4],0)):
            with qc.if_test((qecc[5],1)):
                qc.z(0+7*pos)

    with qc.if_test((qecc[3],0)):             #qbit 1
        with qc.if_test((qecc[4],1)):
            with qc.if_test((qecc[5],0)):
                qc.z(1+7*pos)
    
    with qc.if_test((qecc[3],0)):             #qbit 2
        with qc.if_test((qecc[4],1)):
            with qc.if_test((qecc[5],1)):
                qc.z(2+7*pos)
    
    with qc.if_test((qecc[3],1)):             #qbit 3
        with qc.if_test((qecc[4],0)):
            with qc.if_test((qecc[5],0)):
                qc.z(3+7*pos)
    
    with qc.if_test((qecc[3],1)):             #qbit 4
        with qc.if_test((qecc[4],0)):
            with qc.if_test((qecc[5],1)):
                qc.z(4+7*pos)
    
    with qc.if_test((qecc[3],1)):             #qbit 5
        with qc.if_test((qecc[4],1)):
            with qc.if_test((qecc[5],0)):
                qc.z(5+7*pos)
    
    with qc.if_test((qecc[3],1)):             #qbit 6
        with qc.if_test((qecc[4],1)):
            with qc.if_test((qecc[5],1)):
                qc.z(6+7*pos)

def qec_ft(qc: QuantumCircuit, qecc, pos: int):         #70 gates, 72 depth
    flags = ClassicalRegister(6)
    qc.add_register(flags)
    anc = qc.num_qubits - 1
    ancc = anc - 1
    qc.reset(anc), qc.reset(ancc)
    ##################################Z-Stabilizers##########################################
    qc.h(ancc)
    qc.cx(0+7*pos, anc)
    qc.cx(ancc,anc)
    qc.cx(2+7*pos, anc)
    qc.cx(4+7*pos, anc)
    qc.cx(ancc,anc)
    qc.cx(6+7*pos, anc)

    qc.id(anc), qc.h(ancc), qc.id(ancc)
    qc.measure(anc, qecc[2]), qc.measure(ancc, flags[0])
    qc.reset(anc), qc.reset(ancc)
    qc.id(anc), qc.id(ancc)

    qc.h(ancc)
    qc.cx(1+7*pos, anc)
    qc.cx(ancc, anc)
    qc.cx(2+7*pos, anc)
    qc.cx(5+7*pos, anc)
    qc.cx(ancc, anc)
    qc.cx(6+7*pos, anc)

    qc.id(anc), qc.h(ancc), qc.id(ancc)
    qc.measure(anc, qecc[1]), qc.measure(ancc, flags[1])
    qc.reset(anc), qc.reset(ancc)
    qc.id(anc), qc.id(ancc)

    qc.h(ancc)
    qc.cx(3+7*pos, anc)
    qc.cx(ancc, anc)
    qc.cx(4+7*pos, anc)
    qc.cx(5+7*pos, anc)
    qc.cx(ancc, anc)
    qc.cx(6+7*pos, anc)

    qc.id(anc), qc.h(ancc), qc.id(ancc)
    qc.measure(anc, qecc[0]), qc.measure(ancc, flags[2])
    qc.reset(anc), qc.reset(ancc)
    qc.id(anc), qc.id(ancc)
    ##################################X-Stabilizers##############################################
    qc.h(anc)
    qc.cx(anc, 0+7*pos)
    qc.cx(anc, ancc)
    qc.cx(anc, 2+7*pos)
    qc.cx(anc, 4+7*pos)
    qc.cx(anc, ancc)
    qc.cx(anc, 6+7*pos)
    qc.h(anc)

    qc.id(anc), qc.id(ancc)
    qc.measure(anc, qecc[5]), qc.measure(ancc, flags[3])
    qc.reset(anc), qc.reset(ancc)
    qc.id(anc), qc.id(ancc)

    qc.h(anc)
    qc.cx(anc, 1+7*pos)
    qc.cx(anc, ancc)
    qc.cx(anc, 2+7*pos)
    qc.cx(anc, 5+7*pos)
    qc.cx(anc, ancc)
    qc.cx(anc, 6+7*pos)
    qc.h(anc)

    qc.id(anc), qc.id(ancc)
    qc.measure(anc, qecc[4]), qc.measure(ancc, flags[4])
    qc.reset(anc), qc.reset(ancc)
    qc.id(anc), qc.id(ancc)

    qc.h(anc)
    qc.cx(anc, 3+7*pos)
    qc.cx(anc, ancc)
    qc.cx(anc, 4+7*pos)
    qc.cx(anc, 5+7*pos)
    qc.cx(anc, ancc)
    qc.cx(anc, 6+7*pos)
    qc.h(anc)

    qc.id(anc), qc.id(ancc)
    qc.measure(anc, qecc[3]), qc.measure(ancc, flags[5])
    qc.reset(anc), qc.reset(ancc)
    ##################################Bitflip Error correction##############################################
    
    with qc.if_test((qecc[0],0)):             #qbit 0
        with qc.if_test((qecc[1],0)):
            with qc.if_test((qecc[2],1)):
                qc.x(0+7*pos)

    with qc.if_test((qecc[0],0)):             #qbit 1
        with qc.if_test((qecc[1],1)):
            with qc.if_test((qecc[2],0)):
                qc.x(1+7*pos)
    
    with qc.if_test((qecc[0],0)):             #qbit 2
        with qc.if_test((qecc[1],1)):
            with qc.if_test((qecc[2],1)):
                qc.x(2+7*pos)
    
    with qc.if_test((qecc[0],1)):             #qbit 3
        with qc.if_test((qecc[1],0)):
            with qc.if_test((qecc[2],0)):
                qc.x(3+7*pos)
    
    with qc.if_test((qecc[0],1)):             #qbit 4
        with qc.if_test((qecc[1],0)):
            with qc.if_test((qecc[2],1)):
                qc.x(4+7*pos)
    
    with qc.if_test((qecc[0],1)):             #qbit 5
        with qc.if_test((qecc[1],1)):
            with qc.if_test((qecc[2],0)):
                qc.x(5+7*pos)
    
    with qc.if_test((qecc[0],1)):             #qbit 6
        with qc.if_test((qecc[1],1)):
            with qc.if_test((qecc[2],1)):
                qc.x(6+7*pos)

    ##################################Phaseflip Error correction##############################################
    
    with qc.if_test((qecc[3],0)):             #qbit 0
        with qc.if_test((qecc[4],0)):
            with qc.if_test((qecc[5],1)):
                qc.z(0+7*pos)

    with qc.if_test((qecc[3],0)):             #qbit 1
        with qc.if_test((qecc[4],1)):
            with qc.if_test((qecc[5],0)):
                qc.z(1+7*pos)
    
    with qc.if_test((qecc[3],0)):             #qbit 2
        with qc.if_test((qecc[4],1)):
            with qc.if_test((qecc[5],1)):
                qc.z(2+7*pos)
    
    with qc.if_test((qecc[3],1)):             #qbit 3
        with qc.if_test((qecc[4],0)):
            with qc.if_test((qecc[5],0)):
                qc.z(3+7*pos)
    
    with qc.if_test((qecc[3],1)):             #qbit 4
        with qc.if_test((qecc[4],0)):
            with qc.if_test((qecc[5],1)):
                qc.z(4+7*pos)
    
    with qc.if_test((qecc[3],1)):             #qbit 5
        with qc.if_test((qecc[4],1)):
            with qc.if_test((qecc[5],0)):
                qc.z(5+7*pos)
    
    with qc.if_test((qecc[3],1)):             #qbit 6
        with qc.if_test((qecc[4],1)):
            with qc.if_test((qecc[5],1)):
                qc.z(6+7*pos)
