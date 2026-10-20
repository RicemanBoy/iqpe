from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister, qasm3, qasm2, qpy
from qiskit.visualization import plot_histogram
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from dataclasses import dataclass
import time
from contextlib import ExitStack
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

def phaseflip_lookup_table(n: int):         #syndrome -> data qubits that get a Z, for the phaseflip repcode with stabilizers X_i X_{i+1}
    table = {}
    for value in range(1, 2**(n-1)):        #syndrome 00..0 needs no correction, so it is not in the table
        syndrome = tuple((value >> i) & 1 for i in range(n-1))       #syndrome[i] = qecc[i] = outcome of X_i X_{i+1}
        error = [0]
        for s in syndrome:                  #Z on qubit i flips X_{i-1}X_i and X_i X_{i+1}  -->  e_{i+1} = e_i XOR s_i
            error.append(error[-1] ^ s)
        if sum(error) > n//2:               #two error patterns per syndrome (e and its complement), take the lighter one
            error = [1 - e for e in error]
        table[syndrome] = tuple(i for i, e in enumerate(error) if e)
    return table

LUT_PHASEFLIP_3 = phaseflip_lookup_table(3)
LUT_PHASEFLIP_5 = phaseflip_lookup_table(5)

def phaseflip_mwpm_order(n: int, order = None):            #order in which the checks X_i X_{i+1} are measured each round, default 0, 1, ..., n-2
    order = tuple(range(n - 1)) if order is None else tuple(int(i) for i in order)
    if sorted(order) != list(range(n - 1)):
        raise ValueError(f"order has to be a permutation of 0..{n-2}, got {order}")
    return order

def phaseflip_mwpm_model(n: int, rounds: int, p_x: float, p_z: float, order = None):      #stim copy of the syndrome rounds in qec_mwpm, with the noise of __noise_model__
    import stim
    order = phaseflip_mwpm_order(n, order)
    c = stim.Circuit()
    anc = n
    p1 = [p_x, 0, p_z]                                          #X, Y, Z after every single-qubit gate
    p2 = [p_x/2, 0, p_z/2, p_x/2] + [0]*7 + [p_z/2, 0, 0, 0]    #IX, IY, IZ, XI, ..., ZI after every cx
    c.append("RX", range(n))
    c.append("Z_ERROR", range(n), p_z)                          #errors the block inherits from the gates before it
    for r in range(rounds):
        for i in order:
            c.append("R", [anc])
            c.append("H", [anc]); c.append("PAULI_CHANNEL_1", [anc], p1)
            for q in (i, i + 1):
                c.append("CX", [anc, q]); c.append("PAULI_CHANNEL_2", [anc, q], p2)
            c.append("H", [anc]); c.append("PAULI_CHANNEL_1", [anc], p1)
            c.append("I", [anc]); c.append("PAULI_CHANNEL_1", [anc], p1)
            c.append("M", [anc])
        for i in range(n - 1):                                  #detectors sorted by check (not by measurement order), so detector r*(n-1)+i <-> syndromes[r][i]
            now = stim.target_rec(-(n - 1) + order.index(i))
            if r == 0:                                          #first round is compared to the code space, later rounds to the round before
                c.append("DETECTOR", [now])
            else:
                c.append("DETECTOR", [now, stim.target_rec(-2*(n - 1) + order.index(i))])
    c.append("MX", range(n))                                    #not in the real circuit, only tells the decoder which data qubits carry a Z
    for q in range(n):
        c.append("OBSERVABLE_INCLUDE", [stim.target_rec(-n + q)], q)
    return c

_MWPM_CACHE = {}

def phaseflip_mwpm_decoder(n: int, rounds: int, p_x: float, p_z: float, order = None):    #matching graphs from the circuit-level error model, built once per setting
    import pymatching
    order = phaseflip_mwpm_order(n, order)
    key = (n, rounds, p_x, p_z, order)
    if key not in _MWPM_CACHE:
        model = phaseflip_mwpm_model(n, rounds, p_x, max(p_z, 1e-6), order)  #p_z = 0 would give an empty graph
        dem = model.detector_error_model(decompose_errors=True, approximate_disjoint_errors=True).flattened()

        edges = {}                  #detectors -> [probability, {data qubits with a Z: probability}]
        for ins in dem:
            if ins.type != "error":
                continue
            p = ins.args_copy()[0]
            dets, obs = set(), set()
            for t in ins.targets_copy() + [None]:
                if t is None or t.is_separator():                       #end of one graphlike piece of the fault
                    if dets:                                            #undetectable pieces can't be matched anyway
                        e = edges.setdefault(frozenset(dets), [0.0, {}])
                        e[0] = e[0]*(1 - p) + p*(1 - e[0])
                        e[1][frozenset(obs)] = e[1].get(frozenset(obs), 0) + p
                    dets, obs = set(), set()
                elif t.is_relative_detector_id():
                    dets ^= {t.val}
                elif t.is_logical_observable_id():
                    obs ^= {t.val}

        #one matcher per data qubit with the same graph -> same matching, each tells if its qubit gets a Z
        #(pymatching 2.2.1 + numpy 2 returns wrong vectors for more than one observable)
        #parallel faults (only in the last round: readout error vs. Z on a data qubit) --> most likely one
        matchers = []
        for q in range(n):
            m = pymatching.Matching()
            for dets, (p, obs_probs) in edges.items():
                obs = max(obs_probs, key=obs_probs.get)
                d = sorted(dets)
                fid = {0} if q in obs else set()
                w = np.log((1 - p)/p)
                if len(d) == 1:
                    m.add_boundary_edge(d[0], fault_ids=fid, weight=w, error_probability=p)
                else:
                    m.add_edge(d[0], d[1], fault_ids=fid, weight=w, error_probability=p)
            matchers.append(m)
        _MWPM_CACHE[key] = matchers
    return _MWPM_CACHE[key]

def phaseflip_mwpm_correction(syndromes, n: int, p_x: float, p_z: float, order = None):   #syndromes[r][i] = outcome of X_i X_{i+1} in round r  -->  data qubits that get a Z
    syn = np.asarray(syndromes, dtype=np.uint8)
    det = syn.copy()
    det[1:] ^= syn[:-1]                                         #detection events: changes between rounds
    det = det.reshape(-1)
    return [q for q, m in enumerate(phaseflip_mwpm_decoder(n, len(syn), p_x, p_z, order)) if m.decode(det)[0]]

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

def avg7_ramsey_htoff(distance: int, iter: int, noise: float, qec = False, k = 1, bias = 0, path = "", short=False):       #only exact angles!  
    """
    Simulates the new iQPE Circuit.

    Parameters
    ----------
    distance : int
        Distance of the phaseflip code. Choose either 3 or 5, above will likely not run due to memory reuqirements.
    iter : int
        Amount of iterations of iQPE. Always stick to <3, since the Unitarys are not build for more precision.
    noise : float
        Noiserate of the simulation. Stochastic Pauli Noise is used.
    qec : bool
        De-/activates QEC for the simulation.
    k : int
        Shots
    bias : float
        Noisebias of the qubits. Default set to zero, i.e. no bias. Positive bias --> p_x > p_z, negative bias --> p_z > p_x
    short : bool
        True implements U of each angle, False implements sqrt(U) twice. The False option is an artifact of the original iQPE Circuit where implementing of sqrt(U) was necessary due to the decomp. of controlled-U.

    Returns
    -------
    list
        Performance for each angle, i.e. the distance to the desired angle.
    """
    
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
            start = time.perf_counter()
            for t in range(iter):
                rots = [k*0.5 for k in rots]
                self = RepCode_z(distance, 3)
                self.noise_model = self.__noise_model__(noise, bias)
                self.err = qec
                
                #############################
                for _ in range(2**(iter-t-1)):
                    self.U_htoff(index=2*o+1, short=short)
                ###############################
                for l in rots:
                    if l == 0.25:
                        self.adj_ry_pi_2()
                    if l == 0.125:
                        self.adj_ry_pi_4()

                self.readout(pos=0, shots=1)
        
                if self.zeros == 1:
                    bitstring += "0"
                if self.ones == 1:
                    bitstring += "1"
                    rots.append(0.5)
                del self
            bitstring = bitstring[::-1]
            hmm = convert(bitstring)
            diff = min(np.abs(hmm-angle[2*o+1]), 1-np.abs(hmm-angle[2*o+1]))
            y += diff
            elapsed = time.perf_counter() - start
            print("Performance {}for angle {}: {} ({:.2f}s)".format("(QEC) " if qec else "", 2*o+1, diff, elapsed))
            bruh1.append(diff), y_list.append(diff)
    y = y/(7*k)
    arg = 0
    for i in range(len(bruh1)):
        arg += (y-bruh1[i])**2
    sigma = ((1/(k*n))*arg)**0.5
    sigma = sigma/((k*n)**0.5)

    return y_list

def avg7_RZ_ramsey_htoff(distance: int, iter: int, noise: float, qec = False, k = 1, bias = 0, short = False):       #only exact angles!  
    """
    Simulates the new iQPE Circuit.

    Parameters
    ----------
    distance : int
        Distance of the phaseflip code. Choose either 3 or 5, above will likely not run due to memory reuqirements.
    iter : int
        Amount of iterations of iQPE. Always stick to <3, since the Unitarys are not build for more precision.
    noise : float
        Noiserate of the simulation. Stochastic Pauli Noise is used.
    qec : bool
        De-/activates QEC for the simulation.
    k : int
        Shots
    bias : float
        Noisebias of the qubits. Default set to zero, i.e. no bias. Positive bias --> p_x > p_z, negative bias --> p_z > p_x
    short : bool
        True implements U of each angle, False implements sqrt(U) twice. The False option is an artifact of the original iQPE Circuit where implementing of sqrt(U) was necessary due to the decomp. of controlled-U.

    Returns
    -------
    list
        Performance for each angle, i.e. the distance to the desired angle.
    """
    
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
            start = time.perf_counter()
            for t in range(iter):
                rots = [k*0.5 for k in rots]
                self = RepCode_z(distance, 3)
                self.noise_model = self.__noise_model__(noise, bias)
                self.err = qec
                self.h(pos=0)
                #############################
                for _ in range(2**(iter-t-1)):
                    self.U_RZ_htoff(index=2*o+1, short=short)
                ###############################
                for l in rots:
                    if l == 0.25:
                        self.adj_rz_pi_2()
                    if l == 0.125:
                        self.adj_rz_pi_4()
                self.h(pos=0)
                self.readout(pos=0, shots=1)
        
                if self.zeros == 1:
                    bitstring += "0"
                if self.ones == 1:
                    bitstring += "1"
                    rots.append(0.5)
                del self
            bitstring = bitstring[::-1]
            hmm = convert(bitstring)
            diff = min(np.abs(hmm-angle[2*o+1]), 1-np.abs(hmm-angle[2*o+1]))
            y += diff
            elapsed = time.perf_counter() - start
            print("Performance {}for angle {}: {} ({:.2f}s)".format("(QEC) " if qec else "", 2*o+1, diff, elapsed))
            bruh1.append(diff), y_list.append(diff)
    y = y/(7*k)
    arg = 0
    for i in range(len(bruh1)):
        arg += (y-bruh1[i])**2
    sigma = ((1/(k*n))*arg)**0.5
    sigma = sigma/((k*n)**0.5)

    return y_list

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

        qr = QuantumRegister(n*(logical_q+2)+1, "q")
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

        for i in range(self.n):                     #swap logical qubits such that the target qubit is at the same spot as before for convenience, therefore no noise on these one
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

    def adj_rz_pi_2(self):
        #implements rz(-pi/2). more or less accurate
        lists = [('CX', 0, 1), ('CCX', 0, 1, 2), ('Z', 2), ('CCX', 0, 1, 2)]
        self.cnot(0,1)
        self.toff(0,1,2)
        self.z(2)
        self.toff(0,1,2)

    def adj_rz_pi_4(self):
        #implements rz(-pi/4) with 1e-2 error, error is defined over the algorithm
        lists = [('CCX', 0, 1, 2), ('CX', 0, 2), ('H', 1), ('CCX', 0, 1, 2), ('CX', 2, 1), ('CCX', 1, 2, 0), ('H', 1), ('CX', 2, 1), ('CCX', 0, 1, 2), ('CCX', 1, 2, 0), ('H', 1), ('CCX', 1, 2, 0), ('CCX', 0, 1, 2), ('H', 1), ('CCX', 0, 1, 2), ('CCX', 1, 2, 0), ('H', 1), ('CCX', 0, 1, 2), ('CCX', 1, 2, 0), ('H', 1), ('CCX', 1, 2, 0), ('CCX', 0, 1, 2), ('CX', 0, 1), ('CX', 2, 0), ('H', 1), ('CCX', 0, 1, 2), ('CCX', 1, 2, 0), ('CX', 2, 1), ('H', 1), ('CCX', 0, 1, 2), ('CCX', 1, 2, 0), ('H', 1), ('CCX', 0, 1, 2), ('CCX', 1, 2, 0), ('CX', 0, 1), ('H', 1), ('CX', 0, 1), ('CCX', 0, 1, 2), ('CCX', 1, 2, 0), ('H', 1), ('CCX', 1, 2, 0), ('CCX', 0, 1, 2), ('H', 1), ('CCX', 0, 1, 2), ('CCX', 1, 2, 0), ('H', 1), ('CCX', 0, 1, 2), ('CCX', 1, 2, 0), ('H', 1), ('CCX', 1, 2, 0), ('CCX', 0, 1, 2), ('H', 1), ('CCX', 0, 1, 2), ('CCX', 1, 2, 0), ('H', 1), ('CCX', 0, 1, 2), ('CCX', 1, 2, 0), ('H', 1), ('CCX', 1, 2, 0), ('CCX', 0, 1, 2), ('CX', 2, 1), ('H', 1), ('CCX', 0, 1, 2), ('CCX', 1, 2, 0), ('CX', 0, 1), ('CX', 2, 0), ('H', 1), ('CCX', 1, 2, 0), ('CCX', 0, 1, 2), ('H', 1), ('CCX', 1, 2, 0), ('CCX', 0, 1, 2), ('H', 1), ('CCX', 1, 2, 0), ('CCX', 0, 1, 2), ('H', 1), ('CCX', 0, 1, 2), ('CCX', 1, 2, 0), ('H', 1), ('CCX', 1, 2, 0), ('CCX', 0, 1, 2), ('CX', 0, 1), ('H', 1), ('CCX', 0, 1, 2), ('CX', 2, 1), ('CCX', 1, 2, 0), ('H', 1), ('CX', 0, 1), ('CCX', 0, 1, 2), ('CCX', 1, 2, 0), ('H', 1), ('CCX', 0, 1, 2), ('CCX', 1, 2, 0), ('CX', 0, 1)]

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
                self.qec_shor(pos=targ)               #needed for FT
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
  
    def U_RZ_htoff(self, index=0, short=False):             #latest solution, this implements sqrt(U), so we have to apply this twice!
        if not short:
            with open("rz_htoffcxbeam64_15_1e-2.txt", "r") as f:
                htoff_decomp = [
                    [tuple(item.split("_")[0:1] + [int(x) for x in item.split("_")[1:]])
                    for item in line.strip().split(",")]
                    for line in f
                ]
            for _ in range(2):
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
        else:
            with open("rz_htoffcxbeam64_15_U2_1e-2.txt", "r") as f:
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

    def U_htoff(self, index=0, short=False):             #latest solution, this implements sqrt(U), so we have to apply this twice!
        if not short:
            with open("htoffcxbeam64_15_1e-2.txt", "r") as f:
                htoff_decomp = [
                    [tuple(item.split("_")[0:1] + [int(x) for x in item.split("_")[1:]])
                    for item in line.strip().split(",")]
                    for line in f
                ]
            for _ in range(2):
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
        else:
            with open("htoffcxbeam64_15_U2_1e-2.txt", "r") as f:
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

    def qec_lut(self, pos: int):            #qec for arbitrary d, correction via lookup table
        anc = self.qc.num_qubits - 1
        self.qec_counter += 1

        for i in range(self.n - 1):         #measure X_i X_{i+1} into qecc[i]
            self.qc.reset(anc)
            self.qc.h(anc)
            self.qc.cx(anc, self.n*pos + i)
            self.qc.cx(anc, self.n*pos + i + 1)
            self.qc.h(anc)
            self.qc.id(anc)
            self.qc.measure(anc, self.qecc[i])

        for syndrome, flips in phaseflip_lookup_table(self.n).items():
            with ExitStack() as stack:      #nested if_tests on qecc[0..n-2] only, qecc[n-1] can still hold the H_L measurement
                for i, s in enumerate(syndrome):
                    stack.enter_context(self.qc.if_test((self.qecc[i], s)))
                for q in flips:
                    self.qc.z(self.n*pos + q)

    def _noise_params(self):                #(p_x, p_z) of the single-qubit gates, read back from self.noise_model
        for err in self.noise_model.to_dict()["errors"]:
            if "h" in err["operations"]:
                probs = {ins[0]["name"]: pr for ins, pr in zip(err["instructions"], err["probabilities"])}
                return float(probs.get("x", 0)), float(probs.get("z", 0))
        return 0.0, 0.0

    def qec_mwpm(self, pos: int, rounds = None, order = None):     #qec for arbitrary d: `rounds` syndrome rounds (default d), decoded together with MWPM
        rounds = self.n if rounds is None else rounds                   #order: sequence of checks per round, e.g. (0, 2, 1, 3) for d=5 (fewer fatal fault pairs than the default 0, 1, 2, 3)
        order = phaseflip_mwpm_order(self.n, order)
        anc = self.qc.num_qubits - 1
        self.qec_counter += 1

        syn_reg = ClassicalRegister(rounds*(self.n - 1), "mwpm")
        self.qc.add_register(syn_reg)
        for r in range(rounds):
            for i in order:                 #measure X_i X_{i+1} into syn_reg[r*(n-1) + i]
                self.qc.reset(anc)
                self.qc.h(anc)
                self.qc.cx(anc, self.n*pos + i)
                self.qc.cx(anc, self.n*pos + i + 1)
                self.qc.h(anc)
                self.qc.id(anc)
                self.qc.measure(anc, syn_reg[r*(self.n - 1) + i])

        self.qc.save_statevector(label="psi")           #simulate this trajectory up to here (1 shot), decode the syndromes in python
        sim = AerSimulator(method="statevector", noise_model = self.noise_model)
        result = sim.run(self.qc, shots=1).result()
        psi_full = result.data(0)["psi"]

        bits = list(result.get_counts().keys())[0].split()[0][::-1]       #syn_reg was added last, so it is printed first; [::-1] --> bits[k] = syn_reg[k]
        syndromes = [[int(bits[r*(self.n - 1) + i]) for i in range(self.n - 1)] for r in range(rounds)]

        qr2 = QuantumRegister(self.qc.num_qubits, "q")                     #continue from the sampled state
        self.qc = QuantumCircuit(qr2)
        self.qc.add_register(self.qecc)
        self.qc.set_statevector(psi_full)

        del psi_full

        for q in phaseflip_mwpm_correction(syndromes, self.n, *self._noise_params(), order):
            self.qc.z(self.n*pos + q)

    def _sample_syndrome_rounds(self, pos: int, rounds: int, order):     #append `rounds` syndrome rounds, simulate this trajectory (1 shot), continue from the sampled state
        anc = self.qc.num_qubits - 1

        syn_reg = ClassicalRegister(rounds*(self.n - 1), "syn")
        self.qc.add_register(syn_reg)
        for r in range(rounds):
            for i in order:                 #measure X_i X_{i+1} into syn_reg[r*(n-1) + i]
                self.qc.reset(anc)
                self.qc.h(anc)
                self.qc.cx(anc, self.n*pos + i)
                self.qc.cx(anc, self.n*pos + i + 1)
                self.qc.h(anc)
                self.qc.id(anc)
                self.qc.measure(anc, syn_reg[r*(self.n - 1) + i])

        self.qc.save_statevector(label="psi")
        sim = AerSimulator(method="statevector", noise_model = self.noise_model)
        result = sim.run(self.qc, shots=1).result()
        psi_full = result.data(0)["psi"]

        bits = list(result.get_counts().keys())[0].split()[0][::-1]       #syn_reg was added last, so it is printed first; [::-1] --> bits[k] = syn_reg[k]
        syndromes = [tuple(int(bits[r*(self.n - 1) + i]) for i in range(self.n - 1)) for r in range(rounds)]

        qr2 = QuantumRegister(self.qc.num_qubits, "q")
        self.qc = QuantumCircuit(qr2)
        self.qc.add_register(self.qecc)
        self.qc.set_statevector(psi_full)

        del psi_full
        return syndromes

    def qec_shor(self, pos: int, max_rounds = None, order = None):     #FT qec (Shor): repeat syndrome rounds until the last t+1 agree, then lookup-table correction
        t = (self.n - 1)//2                                             #with <= t faults this stops after at most (t+1)^2 rounds (default cap)
        max_rounds = (t + 1)**2 if max_rounds is None else max_rounds
        order = phaseflip_mwpm_order(self.n, order)
        self.qec_counter += 1

        history = self._sample_syndrome_rounds(pos, min(t + 1, max_rounds), order)
        while len(set(history[-(t + 1):])) > 1 and len(history) < max_rounds:
            history += self._sample_syndrome_rounds(pos, 1, order)
        self.shor_rounds = len(history)                                 #rounds used by the last qec_shor call

        if any(history[-1]):                                            #cap reached without agreement --> last syndrome
            for q in phaseflip_lookup_table(self.n)[history[-1]]:
                self.qc.z(self.n*pos + q)

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

