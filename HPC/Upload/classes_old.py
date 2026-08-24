from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister, qasm3, qasm2, qpy
from qiskit.visualization import plot_histogram
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from functools import wraps, reduce
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

def parity_values(n):
    return [
        value
        for value in range(2**n)
        if value.bit_count() % 2 == 1
    ]

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

###################### Phaselflip code #######################

    def t_phaseflip(self, pos: int):
        for i in range(3):                      #initialize +_L on the phaseflip code
            self.qc.reset(7*self.n + i)
            self.qc.h(7*self.n + i)

        ##### T_L on the phaseflip code #####
        for i in range(2):
            self.qc.cx(self.n*7 + i + 1, self.n*7)
        self.qc.h(self.n*7)
        for i in range(2):
            self.qc.cx(self.n*7 + i + 1, self.n*7)
        
        self.qc.rx(np.pi/4, self.n*7)

        for i in range(2):
            self.qc.cx(self.n*7 + i + 1, self.n*7)
        self.qc.h(self.n*7)
        for i in range(2):
            self.qc.cx(self.n*7 + i + 1, self.n*7)
        ##### T_L on the phaseflip code #####

        for i in range(3):
            self.qc.cx(7*pos + i, self.n*7 + i)

        for i in range(3):
            self.qc.measure(self.n*7 + i, self.qecc[i])

        MASK = 0b000111  # qecc[0], qecc[1], qecc[2]
        masked = expr.bit_and(self.qecc, MASK)
        cond = reduce(expr.logic_or,                #Z_L of the phaseflip code is Z0Z1Z2, so the logical outcome is the parity of the 3 bits, not the majority
                    (expr.equal(masked, v) for v in parity_values(3)))
        with self.qc.if_test(cond):
            for i in range(7):
                self.qc.sdg(7*pos+i)

    def tdg_phaseflip(self, pos: int):
        for i in range(3):                      #initialize +_L on the phaseflip code
            self.qc.reset(7*self.n + i)
            self.qc.h(7*self.n + i)

        ##### T_L on the phaseflip code #####
        for i in range(2):
            self.qc.cx(self.n*7 + i + 1, self.n*7)
        self.qc.h(self.n*7)
        for i in range(2):
            self.qc.cx(self.n*7 + i + 1, self.n*7)
        
        self.qc.rx(-np.pi/4, self.n*7)

        for i in range(2):
            self.qc.cx(self.n*7 + i + 1, self.n*7)
        self.qc.h(self.n*7)
        for i in range(2):
            self.qc.cx(self.n*7 + i + 1, self.n*7)
        ##### T_L on the phaseflip code #####

        for i in range(3):
            self.qc.cx(7*pos + i, self.n*7 + i)

        for i in range(3):
            self.qc.measure(self.n*7 + i, self.qecc[i])

        MASK = 0b000111  # qecc[0], qecc[1], qecc[2]
        masked = expr.bit_and(self.qecc, MASK)
        cond = reduce(expr.logic_or,                #Z_L of the phaseflip code is Z0Z1Z2, so the logical outcome is the parity of the 3 bits, not the majority
                    (expr.equal(masked, v) for v in parity_values(3)))
        with self.qc.if_test(cond):
            for i in range(7):
                self.qc.s(7*pos+i)
     
##################### Phaseflip code #########################

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
            self.qc.id(i+7*pos)
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
        
        ones = (ones/shots)
        zeros = (zeros/shots)
        post = (post/shots)
        preselected = (preselected/shots)

        self.ones = ones
        self.zeros = zeros
        self.post = post
        self.preselected = preselected

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
        cbits = ClassicalRegister(0, "c")
        
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
            self.qc.sdg(i+17*pos)

    def sdg(self, pos = 0):
        for i in range(17):
            self.qc.s(i+17*pos)

    def t(self, pos = 0):
        anc = self.qc.num_qubits - 1
        self.qc.reset(anc)
        self.qc.h(anc), self.qc.t(anc)

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
                    self.qc.s(i+17*pos)
            
    def readout(self, pos: int, shots: int, p = 0):
        read = ClassicalRegister(17)
        self.qc.add_register(read)

        for i in range(17):
            # self.qc.id(i+7*pos)
            self.qc.measure(i+17*pos,read[16-i])

        # self.qc = transpile(self.qc, optimization_level=1)

        code0 = ['11101001100000011', '00111100101100001', '00110011011101000', '11011010010100001', '01010110111000100', '00101001001110110', '00000000110100010', '11111000010100100', '00000000000000000', '01011101110110111', '10001000010011111', '01100101010001110', '01101110000010101', '00101101000101110', '10010010000000001', '11001111101011110', '01001000111101010', '00010101010110101', '01001100000010000', '01001000010100000', '11010001110011000', '10101110110001000', '01101010001001101', '00111000010011011', '11101101101011011', '01111111000010000', '11111100011111100', '00110011000000000', '10001000100111101', '00000100111111010', '11101001001001001', '11111100101011110', '10000111111111110', '11000000101110101', '01100001110011100', '11111100000010100', '10101010010011010', '11111100110110110', '01000111001100011', '00100110100010111', '00110011101001010', '10110100001011100', '10000011000000100', '01100001101110100', '00111100011000011', '10001000001110111', '10101110101100000', '01011001100000111', '00000000101001010', '01110100111000001', '10011001010011010', '01000111100101001', '11001011010100100', '00111100000101011', '10000011011101100', '11110011000111111', '00111000111010001', '10110000110100110', '11001011111101110', '01111111101011010', '00001111000101011', '10111011100111101', '11100110010001010', '10011101101100000', '00110111100010010', '01000011011010011', '11111000100000110', '00011010100111100', '00110111001011000', '10001100000101111', '10100001000000001', '01110000000111011', '00010101100010111', '11100010110011000', '11011110110110011', '11001111011111100', '01010010000111110', '11010101111000000', '01100101100101100', '10101010001110010', '11100110100101000', '11000100100101101', '01010110100101100', '10100101001011001', '00101001111010100', '10000011101001110', '01000111111000001', '01100001000111110', '10111011010011111', '01110100010001011', '11100010000111010', '00111100110001001', '01011001010100101', '00110111111111010', '10111111101100101', '01101110011111101', '00000000011101000', '10011101011000010', '00100010000000101', '01001000001001000', '01110000101110001', '10110100010110100', '10011101110001000', '00010001000000101', '11000000110011101', '00101101110001100', '11110111010001111', '01111011100000010', '00011010001110110', '10000011110100110', '00011110011000110', '00001011001110011', '01101010111101111', '10101110011000010', '00101101101100100', '00010001011101101', '00100110001011101', '00001011010011011', '01011001111101111', '01000011110011001', '01010010011010110', '01000011000111011', '01110000011010011', '00101001010011110', '00001011111010001', '00011010010011110', '01101110101011111', '11010101100101000', '00010101111111111', '10010010110100011', '10100101100010011', '10000111100010110', '10001100101100101', '00000100001011000', '11000100111000101', '10001100011000111', '01011101011111101', '10101110000101010', '11010101010001010', '11000100001100111', '01111111110110010', '00000100010110000', '10100101111111011', '10011001111010000', '11010101001100010', '00110011110100010', '10111111110001101', '10101010100111000', '11110011011010111', '01001100011111000', '10100001101001011', '01100101111000100', '00011110110001100', '00010101001011101', '10001000111010101', '10000111010110100', '10111011111010101', '01101110110110111', '10110000000000100', '00001111101100001', '10010010101001011', '01100101001100110', '00110111010110000', '11001111110110110', '01010110001100110', '00101001100111100', '11100010101110000', '11011010100000011', '11001011100000110', '00100010110100111', '00111000100111001', '00001111011000011', '10010010011101001', '00010001101001111', '00101101011000110', '01011101101011111', '11011110101011011', '10010110001011001', '00010001110100111', '11110011110011101', '10011001100111000', '10110100111111110', '11010001000111010', '11111000111101110', '01110000110011001', '00011010111010100', '11011010001001001', '01000111010001011', '11101001111101011', '01001100110110010', '10011101000101010', '11101101011111001', '10111011001110111', '01111011001001000', '01011101000010101', '10001100110001101', '10100001110100011', '01010010110011100', '11110111111000101', '01110100001100011', '11011010111101011', '10101010111010000', '01101010100000111', '10100101010110001', '00100110111111111', '00001011100111001', '10010110100010011', '01101010010100101', '11100010011010010', '11011110011111001', '01011001001001101', '11010001011010010', '01111011111101010', '11100110111000000', '11101001010100001', '11110111001100111', '10111111000101111', '01000011101110001', '01010010101110100', '00001111110001001', '11011110000010001', '11010001101110000', '11000000000111111', '01110100100101001', '00011110000101110', '00100010011101101', '00100110010110101', '10100001011101001', '11101101110110011', '01111111011111000', '11000000011010111', '10000111001011100', '01001100101011010', '11000100010001111', '01111011010100000', '00000100100010010', '00111000001110011', '10110000101001110', '11110111100101101', '10011001001110010', '11110011101110101', '11001111000010100', '10110000011101100', '00011110101100100', '11101101000010001', '01001000100000010', '11100110001100010', '10110100100010110', '10111111011000111', '01010110010001110', '11001011001001100', '00100010101001111', '01100001011010110', '10010110111111011', '10010110010110001', '11111000001001100']
        code1 = ['01100010111010101', '01111000101001011', '00110100011111001', '00110000010100001', '00100101101011110', '10100010100000010', '01100010010011111', '11011101111111010', '11100001111010001', '00000011111101011', '00101010000111111', '01011110100010110', '10100110110110010', '00001000101110000', '10011010101110001', '01011110111111110', '00011101001100111', '11010110110001001', '10110111101011111', '00100101000010100', '10011010000111011', '10000100000010101', '10011010110011001', '01110111011000010', '11010010100111001', '01011110001011100', '01110111110001000', '10010101101011010', '01001111100010011', '00100001100000110', '00101010110011101', '10110011100000111', '01000000111010000', '01001011101001011', '01011010110100110', '10111100111000100', '00100101110110110', '00110100101011011', '01111100100010011', '01111000000000001', '01111100111111011', '10101101010001011', '11000011010011110', '10101101100101001', '00010010010100100', '01100010001110111', '10100010001001000', '11101010101001010', '00101110111000101', '10000000111101111', '00000011100000011', '10001011101110100', '10001111100101100', '01101101100010110', '11100101000101011', '01010001001110111', '00001000110011000', '00001000011010010', '10000000010100101', '00011101100101101', '00000111000010001', '01001111001011001', '11111011101001111', '01001011110100011', '10010001010100000', '01101001110100110', '10001011000111110', '01000000010011010', '00101110010001111', '00011001011010111', '10101101001100011', '01110011111010000', '10011110001100011', '01001011000000001', '10101001101110001', '01110111101100000', '10011110111000001', '01010101000101111', '11111011000000101', '11110100101100100', '01011110010110100', '10110111110110111', '00000011010100001', '11110100000101110', '11001100100010111', '11010010001110011', '10101101111000001', '11010110101100001', '00011101111000101', '00000111110110011', '00100001010100100', '01011010011101100', '10110111000010101', '00110100110110011', '11000111110001100', '01010101110001101', '11100001100111001', '10000100011111101', '11010010010011011', '10011110100101001', '11011101001011000', '01000100110001000', '11111011011101101', '11111111010110101', '11101010000000000', '10010101110110010', '11001000101001111', '00111011101110000', '11110100011000110', '01101001011101100', '10110111011111101', '01111000011101001', '00101110001100111', '10101001110011001', '00010010111101110', '10100110000010000', '00010110101011110', '00111111010001010', '01000100101100000', '10010001100000010', '11001000110100111', '11110000100111100', '11000011111010100', '11111111001011101', '01110011001110010', '10000000100000111', '10001011011010110', '01101101010110100', '10010101011111000', '11100001010011011', '00111111100101000', '00111011011010010', '00000011001001001', '11000011100111100', '00011101010001111', '11011101100010010', '00111111111000000', '01110011100111000', '01100110110001101', '11101010011101000', '00001100001100010', '10110011010100101', '10011010011010011', '00110100000010001', '11011101010110000', '00011001101110101', '01010101011000111', '00000111101011011', '01000100011000010', '10001011110011100', '11010010111010001', '01101101001011100', '00111111001100010', '01010001010011111', '11001000000000101', '10111000110011100', '11001100001011101', '11000011001110110', '01100010100111101', '00010010100000110', '00001100100101000', '01011010000000100', '11101110001011000', '01000000001110010', '00001000000111010', '11101110010110000', '11111111100010111', '11100101110001001', '11111111111111111', '00010010001001100', '10111100001100110', '00110000001001001', '11011001110100010', '01100110011000111', '11010110011000011', '00011001110011101', '10100010111101010', '11111011110100111', '10010001111101010', '11100101101100001', '11101110100010010', '11110000010011110', '10100110011111000', '01111000110100011', '01101001101001110', '10100110101011010', '00010110000010100', '11110000001110110', '01001011011101001', '10000100110110111', '10010001001001000', '11001000011101101', '00110000100000011', '00011001000111111', '00000111011111001', '01100110101100101', '00100101011111100', '10111000011010110', '10110011111101111', '11100101011000011', '10111100010001110', '01111100001011001', '00111011000111010', '10001111010001110', '11011001000000000', '11001100010110101', '01110111000101010', '00001100111000000', '10011110010001011', '10111100100101100', '00001100010001010', '11011001011101000', '00100001111101110', '11100001001110011', '01111100010110001', '10000100101011111', '10010101000010000', '01001111010110001', '10001111111000100', '11110100110001100', '01101001000000100', '01010101101100101', '10110011001001101', '10101001011010011', '00010110110110110', '01100110000101111', '11000111101100100', '11011001101001010', '01010001100111101', '00101110100101101', '01000100000101010', '10101001000111011', '00100001001001100', '10111000101110100', '10001111001100110', '11000111000101110', '01010001111010101', '10000000001001101', '11101010110100010', '01001111111111011', '00110000111101011', '00101010011010111', '00101010101110101', '11001100111111111', '11110000111010100', '10111000000111110', '00111011110011000', '01101101111111110', '11101110111111010', '11010110000101011', '01110011010011010', '01011010101001110', '11000111011000110', '00010110011111100', '01000000100111000', '10100010010100000']
                
        sim = AerSimulator()
        job = sim.run(self.qc, shots=shots)
        result = job.result()
        counts = result.get_counts()

        bitstring = list(counts.keys())
        bitstring = [i.replace(" ","") for i in bitstring]
        hmm = list(counts.values())

        print("Total shots:", sum(hmm))

        bits = [i[:17] for i in bitstring]

        print(bits)

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
            shots = 1
            self.ones = (self.ones/shots)
            self.zeros = (self.zeros/shots)
            self.post = (self.post/shots)
            self.preselected = (self.preselected/shots)
        return counts
