from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister
import itertools
import numpy as np
import matplotlib.pyplot as plt

from qiskit_aer import AerSimulator

from qiskit_aer.noise import (NoiseModel, pauli_error)

from qiskit.circuit.library import UnitaryGate

matrix_h = ([[2**(-0.5),2**(-0.5)],[2**(-0.5),-2**(-0.5)]])
h_ideal = UnitaryGate(matrix_h)

matrix_cx = ([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
cx_ideal = UnitaryGate(matrix_cx)       #Erst Target, dann Control Qubit!!

matrix_x = ([[0,1],[1,0]])
x_ideal = UnitaryGate(matrix_x)

matrix_z = ([[1,0],[0,-1]])
z_ideal = UnitaryGate(matrix_z)

matrix_t = ([[1,0],[0,np.exp(1j*(np.pi/4))]])
t_ideal = UnitaryGate(matrix_t)

matrix_tdg = ([[1,0],[0,np.exp(-1j*(np.pi/4))]])
tdg_ideal = UnitaryGate(matrix_tdg)

def parity_values(n):
    return [
        value
        for value in range(2**n)
        if value.bit_count() % 2 == 1
    ]

def gates(qc:QuantumCircuit):
    hmm = dict(qc.count_ops())
    hmm["reset"] = 0
    hmm["measure"] = 0
    hmm["swap"] = 0
    return sum(hmm.values())

def convert(bin: str):                  #konvertiert den bitstring in decimal, e.g. 0110 = 0.375
    k = list(bin)
    a = [int(i) for i in k]
    n = 0
    for i in range(len(a)):
        if a[i] == 1:
            n += 1/2**(i+1)
    return n

def avg7_ramsey(code: str, iter: int, noise: float, qec = False, k = 1, bias = 0, post = False, path = "", ftec = False):       #only exact angles!
    #ftec = True runs the Flag 2-FTEC Protocol instead of the plain qec() round
    n = 15
    angle = np.linspace(0,1,n+2)
    angle = np.delete(angle, [n+1])
    angle = np.delete(angle, [0])

    a, b = [], []
    with open("{}unitary{}_improved_with_rs.txt".format(path, n), "r") as file:
        for line in file:
            a.append(list(map(str, line.strip().split(","))))
    with open("{}adjunitary{}_improved_with_rs.txt".format(path, n), "r") as file:
        for line in file:
            b.append(list(map(str, line.strip().split(","))))

    y = 0
    y_list, bruh1 = [], []
    for _ in range(k):
        for o in range(7):
            gatecount = 0
            bitstring = ""
            rots = []
            for t in range(iter):
                rots = [k*0.5 for k in rots]
                counter = 0
                while True:
                    if code == "bigsteane":
                        self = BigSteane17q(1, 0, noise=noise)
                    else:
                        print("Error hahaha")
                        return

                    self.err = qec
                    self.ftec = ftec
                    self.h(pos=0)
                    #############################
                    for _ in range(2**(iter-t-1)):
                        self.cu_ramsey(a[2*o+1])
                    ###############################
                    for l in rots:
                        if l == 0.25:
                            self.sdg(pos=0)
                        if l == 0.125:
                            self.tdg(pos=0)
                    self.h(pos=0)
                    if self.err:
                        self.ec(0)

                    self.readout(pos=0, shots=1, p = noise)
                    gatecount += self.gatecount

                    if self.zeros == 1:
                        bitstring += "0"
                        break
                    if self.ones == 1:
                        bitstring += "1"
                        rots.append(0.5)
                        break
                    if not post:
                        if np.random.rand() < 0.5:
                            bitstring += "0"
                            break
                        else:
                            bitstring += "1"
                            rots.append(0.5)
                            break
                    counter += 1
                    print("Angle {}, {}% errorrate, Iteration {}: {} Repetition".format(2*o+1, noise*100, t, counter))
                    del self
            bitstring = bitstring[::-1]
            hmm = convert(bitstring)
            diff = min(np.abs(hmm-angle[2*o+1]), 1-np.abs(hmm-angle[2*o+1]))
            y += diff
            print("Performance {}for angle {} ({} gates): ".format("(QEC) " if qec else "", 2*o+1, gatecount), diff)
            bruh1.append(diff), y_list.append(diff)
    y = y/(7*k)
    arg = 0
    for i in range(len(bruh1)):
        arg += (y-bruh1[i])**2
    sigma = ((1/(k*n))*arg)**0.5
    sigma = sigma/((k*n)**0.5)

    return y_list

############################## [[17,1,5]] 4.8.8 Color Code ##############################
# The 8 faces of the square-octagon lattice (https://arxiv.org/pdf/2601.13313, Fig. 3b).
# The code is self-dual, so the X- and Z-type generators sit on the same faces. It is the
# same [[17,1,5]] color code that is listed in Table 7 of arXiv:1708.02246 (seven weight-4
# and one weight-8 generator, 36 CNOTs per basis). Qubits 0, 1, 12 are the three corners of
# the triangle, the bottom row of the lattice (0,1,2,4,6) is the logical operator used in
# x()/z().
STABS_17Q = [(5, 9, 11, 12),                # square
             (8, 10, 13, 15),               # square
             (3, 7, 14, 16),                # square
             (2, 6, 14, 16),                # octagon, cut by the boundary
             (4, 5, 6, 7, 8, 10, 11, 16),   # octagon
             (1, 4, 10, 13),                # octagon, cut by the boundary
             (8, 9, 11, 15),                # octagon, cut by the boundary
             (0, 2, 3, 14)]                 # octagon, cut by the boundary

LOGICAL_17Q = (0, 1, 2, 4, 6)

# The order in which the data qubits of a face are coupled to the measurement qubit. For the
# weight-4 faces any order works, but the weight-8 face needs a specific one: with the
# straightforward order the flag error sets of the octagon contain two logically inequivalent
# errors with the same syndrome, so the flag 2-FTEC condition (Definition 10) is violated. Of
# the 8! orders exactly one family works, ORDER_17Q[4] is its lexicographically first member
# (verified by enumerating all errors, as the paper does for this code in Section 2.3).
ORDER_17Q = {i: list(face) for i, face in enumerate(STABS_17Q)}
ORDER_17Q[4] = [4, 5, 6, 8, 7, 10, 16, 11]

# 2-flag circuits (Definition 6). The list is the order of the CNOTs on the measurement qubit,
# ("d", k) is a CNOT_dm coupling the k'th data qubit of ORDER_17Q, ("f", j) a CNOT_fm coupling
# the j'th flag qubit. The measurement qubit is always the target (Section 3.3).
FLAG_SEQ_17Q = {
    # weight-4: Fig. 2(b), which is Fig. 11(a) with w = 4 and needs a single flag qubit. The
    # paper notes it is even a 4-flag circuit, so it doubles as a 2-flag circuit.
    4: [("d", 0), ("f", 0), ("d", 1), ("d", 2), ("f", 0), ("d", 3)],
    # weight-8: the general 2-flag construction of Section 3.3 with w/2-1 = 3 flag qubits.
    #   1. a CNOT_fm pair between the first and the second last CNOT_dm         -> flag 0
    #   2. a CNOT_fm pair between the second and the last CNOT_dm               -> flag 1
    #   3. the first CNOT_fm of the remaining pair two CNOT_dm after the second
    #      CNOT_fm, its partner three CNOT_dm later                             -> flag 2
    8: [("d", 0), ("f", 0), ("d", 1), ("f", 1), ("d", 2), ("d", 3), ("f", 2),
        ("d", 4), ("d", 5), ("f", 0), ("d", 6), ("f", 1), ("f", 2), ("d", 7)],
}
NFLAG_17Q = {w: 1 + max(k for kind, k in seq if kind == "f") for w, seq in FLAG_SEQ_17Q.items()}

# Each 2-flag circuit gets its own slice of the flags register, the weight-8 face needs three
# bits. flags[0..9] are the Z-Stabilizer circuits, flags[10..19] the X-Stabilizer ones.
FLAGBIT_17Q, NFLAGBITS_17Q = {}, 0
for _xtype in (False, True):
    for _i, _face in enumerate(STABS_17Q):
        FLAGBIT_17Q[(_i, _xtype)] = NFLAGBITS_17Q
        NFLAGBITS_17Q += NFLAG_17Q[len(_face)]

# n_max = (t^2+3t+2)/2 rounds are enough for t = 2 ("the protocol is repeated at most 6 times")
NMAX_17Q = 6


def _mask(qubits):
    m = 0
    for q in qubits:
        m |= 1 << q
    return m

FACE_17Q = [_mask(face) for face in STABS_17Q]
LMASK_17Q = _mask(LOGICAL_17Q)

# the 256 elements of the stabilizer group, as bit masks over the 17 data qubits
STABGROUP_17Q = {0}
for _face in FACE_17Q:
    STABGROUP_17Q |= {s ^ _face for s in STABGROUP_17Q}

# readout(): the computational basis support of the two logical states
CODE0_17Q = frozenset("".join("1" if (s >> q) & 1 else "0" for q in range(17)) for s in STABGROUP_17Q)
CODE1_17Q = frozenset("".join("1" if ((s ^ LMASK_17Q) >> q) & 1 else "0" for q in range(17)) for s in STABGROUP_17Q)


def _sector_syndrome(err: int):
    #syndrome of a single sector: bit i = parity of the error on face i
    return sum((bin(err & FACE_17Q[i]).count("1") & 1) << i for i in range(8))

def _min_weight_table():
    #E_min(s) for one sector. The code is self-dual, so the same table decodes the X- and the
    #Z-sector. Errors of weight >= 3 are beyond the correctable radius, for those any
    #representative of the syndrome is fine: it only has to project back into the code space,
    #which is the second criterion of Definition 3.
    table = {}
    for w in range(18):
        for qubits in itertools.combinations(range(17), w):
            table.setdefault(_sector_syndrome(_mask(qubits)), qubits)
        if len(table) == 256:
            break
    return table

CORR_17Q = _min_weight_table()

# every error of weight <= 1, the set E_1 that multiplies the flag error sets in Eq. 4
WEIGHT1_17Q = [(0, 0)]
for _q in range(17):
    WEIGHT1_17Q += [(1 << _q, 0), (0, 1 << _q), (1 << _q, 1 << _q)]


def syndrome_of(err):
    #s(E) of a Pauli E = (xmask, zmask): the Z-Stabilizers detect the X part, the
    #X-Stabilizers the Z part.
    return _sector_syndrome(err[0]), _sector_syndrome(err[1])

def equivalent(e1, e2):
    #E ~ E' : the two differ by an element of the stabilizer group
    return (e1[0] ^ e2[0]) in STABGROUP_17Q and (e1[1] ^ e2[1]) in STABGROUP_17Q


############################ Flag error sets of Definition 9 / Eq. 9 ############################
# E_m(g_i, f_i) is built by enumerating every set of exactly m CNOT faults in the 2-flag circuit
# C(g_i) and keeping those runs that flagged. Section 3.3: "any set of v faults (including those
# at idle, preparation or measurement locations) will have the same output Pauli operator and
# flag measurement results as some set of at most v faults on CNOT gates", so enumerating CNOT
# faults is enough. A fault is one of the 15 non-trivial two qubit Paulis after the gate.

_PAULIS2 = [p for p in itertools.product((0, 1), repeat=4) if any(p)]

def _propagate(seq, w, faults):
    #Pauli frame of the Z-type circuit. Data/flag qubits are the controls, the measurement qubit
    #is the target, so X errors run control -> target and Z errors target -> control.
    err = {}
    for gi, (kind, k) in enumerate(seq):
        c = (kind, k)
        xc, zc = err.get(c, (0, 0))
        xm, zm = err.get("m", (0, 0))
        err[c] = (xc, zc ^ zm)
        err["m"] = (xm ^ xc, zm)
        if gi in faults:
            fxc, fzc, fxm, fzm = faults[gi]
            a, b = err[c]; err[c] = (a ^ fxc, b ^ fzc)
            a, b = err["m"]; err["m"] = (a ^ fxm, b ^ fzm)
    #flag qubits sit in |+> and are read out in X, so a Z error on them flags
    flags = tuple(err.get(("f", j), (0, 0))[1] for j in range(NFLAG_17Q[w]))
    xmask = zmask = 0
    for k in range(w):
        x, z = err.get(("d", k), (0, 0))
        xmask |= x << k
        zmask |= z << k
    #the measurement qubit is prepared in |0> and read out in Z, so an X error on it flips the
    #outcome. Under the H-conjugation that turns the circuit into the X-Stabilizer version the
    #same component flips the X basis readout.
    return flags, xmask, zmask, err.get("m", (0, 0))[0]

def _local_flag_sets(w, nfaults):
    #{flag pattern -> set of local data errors}, reduced modulo the stabilizer being measured
    seq = FLAG_SEQ_17Q[w]
    full, out = (1 << w) - 1, {}
    for gates_ in itertools.combinations(range(len(seq)), nfaults):
        for combo in itertools.product(_PAULIS2, repeat=nfaults):
            flags, xmask, zmask, _ = _propagate(seq, w, dict(zip(gates_, combo)))
            if not any(flags):
                continue
            err = min((xmask, zmask), (xmask, zmask ^ full),
                      key=lambda e: bin(e[0]).count("1") + bin(e[1]).count("1"))
            out.setdefault(flags, set()).add(err)
    return out

_LOCAL_SETS = {(w, m): _local_flag_sets(w, m) for w in FLAG_SEQ_17Q for m in (1, 2)}

def _lift(i, xtype, xmask, zmask):
    #local masks over the qubits of ORDER_17Q[i] -> masks over the 17 data qubits. The
    #X-Stabilizer circuit is the H-conjugate of the Z-Stabilizer one, so its errors are the
    #ones of the Z-type circuit with X and Z exchanged.
    if xtype:
        xmask, zmask = zmask, xmask
    gx = gz = 0
    for k, q in enumerate(ORDER_17Q[i]):
        if (xmask >> k) & 1:
            gx |= 1 << q
        if (zmask >> k) & 1:
            gz |= 1 << q
    return gx, gz

# FLAGSETS_17Q[(i, xtype, m)][flag pattern] = E_m(g_i, f_i) as global Paulis
FLAGSETS_17Q = {}
for _i, _face in enumerate(STABS_17Q):
    for _xtype in (False, True):
        for _m in (1, 2):
            FLAGSETS_17Q[(_i, _xtype, _m)] = {
                f: frozenset(_lift(_i, _xtype, x, z) for x, z in s)
                for f, s in _LOCAL_SETS[(len(_face), _m)].items()
            }


def flag_candidates(flagged, m: int):
    #E_m(g_i1,..,g_ik, f_i1,..,f_ik) of Eq. 9: the errors of exactly m faults spread over the k
    #flagged circuits, which produced the observed flag patterns. Every flagged circuit needs at
    #least one fault, so k <= m. flagged is the list of flag events, one entry per circuit that
    #flagged: if the same generator flags in two different rounds those are two circuits, the two
    #faults sat in two different runs of C(g_i), so the error comes from E_1(g_i,f_1) x E_1(g_i,f_2)
    #and not from E_2(g_i,f).
    items = list(flagged)
    if len(items) == m:                                     #one fault per flagged circuit
        sets = [FLAGSETS_17Q[(i, xtype, 1)].get(pattern, frozenset())
                for (i, xtype), pattern in items]
    elif len(items) == 1 and m == 2:                        #both faults inside the same circuit
        (i, xtype), pattern = items[0]
        sets = [FLAGSETS_17Q[(i, xtype, 2)].get(pattern, frozenset())]
    else:
        return []
    out = [(0, 0)]
    for s in sets:
        out = [(a[0] ^ b[0], a[1] ^ b[1]) for a in out for b in s]
    return out

def correction_set(flagged, m: int, syndrome):
    #E~_2^m(g_i1,..,g_ik, s, f_i1,..,f_ik) of Eq. 4/9: an element of
    #E_m(g_i1,..,g_ik) x E_(2-m) whose syndrome is s, or None if there is no such element.
    rest = WEIGHT1_17Q if m == 1 else [(0, 0)]
    for base in flag_candidates(flagged, m):
        for extra in rest:
            err = (base[0] ^ extra[0], base[1] ^ extra[1])
            if syndrome_of(err) == syndrome:
                return err
    return None

def correction(flagged = (), syndrome = (0, 0), orders = ()):
    #any correction from the union of the E~_2^m(..., s) for m in orders, falling back to
    #E_min(s) ("{E_min(s)} if above set empty")
    for m in orders:
        err = correction_set(flagged, m, syndrome)
        if err is not None:
            return err
    return _mask(CORR_17Q[syndrome[0]]), _mask(CORR_17Q[syndrome[1]])


def flag2_ftec_protocol(flag_round, nonflag_round, correct, clean_run = False):
    #"Flag 2-FTEC Protocol" of arxiv:1708.02246, Section 2.2 (pages 7 and 8), together with the
    #update rules for n_diff and n_same given right above it.
    #   flag_round()    -> ([((generator, xtype), flag pattern), ..], syndrome) of a round with
    #                      the 2-flag circuits
    #   nonflag_round() -> syndrome of a round with the non-flag circuits
    #   correct(flagged, syndrome, orders) applies a correction from the union of the
    #                      E~_2^m(flagged, syndrome) for m in orders, E_min(syndrome) if orders = ()
    #
    #clean_run tightens the 1.Case, see the comment at the threshold below. It is off by default,
    #so that the protocol is the one printed in the paper.
    n_diff, n_same = 0, 0
    diff_rose = False               #did n_diff increase in the previous round?
    run_charged = False             #did the first round of the current run of equal syndromes
                                    #increase n_diff?
    last_syn = None                 #syndrome of the previous round, None if that round flagged
    flagged = []                    #the circuits C(g_i) that flagged so far, with their flag patterns

    for _ in range(NMAX_17Q):
        new, syn = flag_round()

        ##################################### update rules #####################################
        if new:
            n_same = 0                                  #2. rule: a flag resets n_same
            flagged.extend(new)
            last_syn, diff_rose = None, False           #a flagged round records no syndrome
            run_charged = False
        else:
            if last_syn is None:
                run_charged = False                     #a fresh run after a flagged round
            elif syn == last_syn:
                n_same += 1                             #3. rule
                diff_rose = False
            else:
                #n_same counts how often the syndrome repeated in a row, so a syndrome that
                #differs from the previous one starts a new run. Without this the 1.Case would
                #fire on a syndrome that was seen only once.
                n_same = 0
                if not diff_rose:
                    n_diff += 1                         #1. rule
                    diff_rose = True
                    run_charged = True
                else:
                    diff_rose = False
                    run_charged = False
            last_syn = syn

        ###################################### corrections #####################################
        #3.Case: two circuits C(g_i) and C(g_j) flagged --> non-flag syndrome measurement and any
        #correction from E~_2^2(g_i, g_j, s)
        if len(flagged) >= 2:
            return correct(flagged, nonflag_round(), (2,))
        #4.Case: a circuit C(g_i) flagged and n_diff = 1 --> non-flag syndrome measurement and any
        #correction from E~_2^1(g_i, s)
        if flagged and n_diff == 1:
            return correct(flagged, nonflag_round(), (1,))
        #5.Case: a circuit C(g_i) flagged, n_diff = 0 and the same syndrome came up twice in a row
        #--> keep that syndrome and correct with E~_2^1(g_i, s) u E~_2^2(g_i, s)
        if flagged and n_diff == 0 and n_same == 1:
            return correct(flagged, syn, (1, 2))
        if not flagged:
            #2.Case: no flags and n_diff = 2 --> non-flag syndrome measurement and E_min(s)
            if n_diff == 2:
                return correct([], nonflag_round(), ())
            #1.Case: no flags and the same syndrome repeated 3 - n_diff times in a row --> E_min(s)
            #The proof of this case (Appendix A) needs one of those rounds to have been fault
            #free. That does not follow when the fault n_diff was charged for sits in the first
            #round of the run itself: with the sequence s_0, s_1, s_1 and n_diff = 1 one fault can
            #produce s_1 in round 1 and a second fault can reproduce the same wrong syndrome in
            #round 2, so neither round is fault free. Two faults then leave a weight-3 residual
            #(0.005% of random 2-fault configurations, see the test suite). clean_run demands one
            #extra repetition in that situation, which restores the guarantee.
            if n_same == 2 - n_diff + (1 if clean_run and run_charged else 0):
                return correct([], syn, ())

    #with more than two faults none of the five cases has to trigger within n_max rounds, so fall
    #back to the 2.Case, which still projects the state back into the code space
    return correct([], nonflag_round(), ())


class BigSteane17q:
    def __init__(self, n: int, magic = 0, noise = 0):
        self.n = n

        self.zeros = 0
        self.ones = 0
        self.preselected = 0
        self.post = 0
        self.err = False
        self.qec_counter = 0
        self.magiccounter = 0
        self.gatecount = 0

        self.preselection_flag = False      #the |0_L> preparation below is unverified, so there
                                            #is nothing to preselect on. Kept for readout().

        self.classical_ec = False
        self.postselection = True

        self.ftec = False                   #ec() then runs the plain qec() round, set it to True
                                            #to run the Flag 2-FTEC Protocol instead
        self.clean_run = False              #see the 1.Case in flag2_ftec_protocol(). False runs
                                            #the protocol exactly as it is printed in the paper.

        self.noise_model = self.__noise_model__(noise, 0)

        #one measurement qubit plus the three flag qubits of the weight-8 2-flag circuit
        self.nanc = 1 + max(NFLAG_17Q.values())
        qr = QuantumRegister(17*n + self.nanc, "q")
        cbits = ClassicalRegister(1, "c")
        self.cbits = cbits

        self.qc = QuantumCircuit(qr, cbits)

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

        self.qecc = ClassicalRegister(16)       #qecc[0..7] Z-Stabilizers, qecc[8..15] X-Stabilizers
        self.qc.add_register(self.qecc)

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

    def id(self, pos: int):
        for i in range(17):
            self.qc.id(i+17*pos)

    def x(self, pos: int):
        for q in LOGICAL_17Q:
            self.qc.x(q+17*pos)

    def z(self, pos: int):
        for q in LOGICAL_17Q:
            self.qc.z(q+17*pos)

    def h(self, pos: int):
        for i in range(17):
            self.qc.h(i+17*pos)

    def s(self, pos: int):
        #every stabilizer element has weight 0 mod 4 and every representative of the logical
        #operator weight 1 mod 4, so the transversal S acts as diag(1, i) = S on the logical
        #qubit. Unlike the [[7,1,3]] code in steane_ftec.py, s and sdg are therefore not swapped.
        for i in range(17):
            self.qc.s(i+17*pos)

    def sdg(self, pos: int):
        for i in range(17):
            self.qc.sdg(i+17*pos)

    def cnot(self, control: int, target: int):
        for i in range(17):
            self.qc.cx(i+17*control, i+17*target)

    def cz(self, control: int, target: int):
        self.h(pos=control)
        self.cnot(control=target, target=control)
        self.h(pos=control)

    def t(self, pos: int):
        self.h(pos)
        self.h(pos)
        self.s(pos)
        #state injection with a single ancilla, as in the [[17,1,5]] class of oldclasses.py. The
        #17 qubit code has no room for a magic state block in a statevector simulation, so this
        #gate is not fault tolerant.
        self.magiccounter += 1
        anc = self.qc.num_qubits - 1
        self.qc.reset(anc)
        self.qc.h(anc)
        self.qc.t(anc)

        for q in LOGICAL_17Q:
            self.qc.cx(q+17*pos, anc)

        self.qc.measure(anc, 0)

        with self.qc.if_test((0,1)):
            for i in range(17):
                self.qc.s(i+17*pos)

        self.sdg(pos)
        self.h(pos)
        self.h(pos)
        if self.err:
            self.ec(pos)

    def tdg(self, pos: int):
        self.h(pos)
        self.h(pos)
        self.s(pos)

        self.magiccounter += 1
        anc = self.qc.num_qubits - 1
        self.qc.reset(anc)
        self.qc.h(anc)
        self.qc.tdg(anc)

        for q in LOGICAL_17Q:
            self.qc.cx(q+17*pos, anc)

        self.qc.measure(anc, 0)

        with self.qc.if_test((0,1)):
            for i in range(17):
                self.qc.sdg(i+17*pos)

        self.sdg(pos)
        self.h(pos)
        self.h(pos)
        if self.err:
            self.ec(pos)

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
            if i == "x":
                self.x(pos=pos)

    def cu_ramsey(self, gate: list):
        self.u2(0, gate=gate)
        if self.err:
            self.ec(pos=0)
        self.u2(0, gate=gate)

################################# FTQEC Protocol based on arxiv:1708.02246, pages 7 and 8, "Flag 2-FTEC Protocol" ##########################
    def _measure_stab(self, pos: int, i: int, xtype: bool, flags = None):
        #one stabilizer measurement. xtype: the measurement qubit is the control (X-Stabilizer),
        #otherwise the data qubits are the controls and the measurement qubit is the target
        #(Z-Stabilizer). With flags != None the CNOT_fm gates of FLAG_SEQ_17Q are inserted, which
        #turns the circuit into the 2-flag circuit C(g_i), without them it is the plain
        #("non-flag") circuit the protocol falls back to.
        face = ORDER_17Q[i]
        w = len(face)
        anc = self.qc.num_qubits - 1
        fq = [anc-1-j for j in range(NFLAG_17Q[w])]
        seq = FLAG_SEQ_17Q[w] if flags is not None else [("d", k) for k in range(w)]
        used = sorted({k for kind, k in seq if kind == "f"})

        self.qc.reset(anc)
        self.qc.id(anc)
        for j in used:
            self.qc.reset(fq[j])
            self.qc.id(fq[j])

        if xtype:
            #H-conjugate of the Z-type circuit: the measurement qubit is prepared in |+> and read
            #out in X, the flag qubits are the targets and stay in |0>/Z.
            self.qc.h(anc)
            for kind, k in seq:
                self.qc.cx(anc, face[k]+17*pos if kind == "d" else fq[k])
            self.qc.h(anc)
        else:
            for j in used:
                self.qc.h(fq[j])                    #flag qubits are prepared in |+>
            for kind, k in seq:
                self.qc.cx(face[k]+17*pos if kind == "d" else fq[k], anc)
            for j in used:
                self.qc.h(fq[j])                    #and read out in the X basis

        self.qc.id(anc)
        self.qc.measure(anc, self.qecc[i + (8 if xtype else 0)])
        for j in used:
            self.qc.id(fq[j])
            self.qc.measure(fq[j], flags[FLAGBIT_17Q[(i, xtype)] + j])
        self.qc.reset(anc)
        for j in used:
            self.qc.reset(fq[j])

    def _round(self, pos: int, flags = None):
        #one full syndrome extraction round over all 16 generators
        self.qec_counter += 1
        for xtype in (False, True):
            for i in range(8):
                self._measure_stab(pos, i, xtype, flags)

    def flagsyndrome(self, pos: int):
        flags = ClassicalRegister(NFLAGBITS_17Q)
        self.qc.add_register(flags)
        self._round(pos=pos, flags=flags)
        return flags

    def reset_qec(self, result):
        psi_full = result.data(0)["psi"]
        qr2 = QuantumRegister(17*self.n + self.nanc, "q")
        self.qc = QuantumCircuit(qr2, self.cbits)       #same register as in __init__, otherwise the
                                                        #classical layout shifts and readout() misslices
        self.qc.add_register(self.qecc)
        self.qc.set_statevector(psi_full)
        del psi_full

    def _creg_bits(self, bitstring: str, creg):
        #the counts bitstring lists the classical registers in reverse order of when they were added,
        #so the position of a register depends on how many registers the circuit currently has
        start = 0
        for reg in reversed(self.qc.cregs):
            if reg == creg:                     #qiskit rebuilds the register objects, so compare by name/size
                return bitstring[start:start+reg.size]
            start += reg.size
        return None

    def _flag_round(self, pos: int):
        #one round of syndrome extraction with the 2-flag circuits --> (flags, syndrome, bitstring, result)
        flags = self.flagsyndrome(pos=pos)
        self.qc.save_statevector(label="psi")
        self.gatecount += gates(self.qc)
        sim = AerSimulator(method="statevector", noise_model = self.noise_model)
        result = sim.run(self.qc, shots=1).result()
        counts = result.get_counts()                #one shot --> one bitstring

        bitstring = list(counts.keys())
        bitstring = [i.replace(" ","") for i in bitstring][0]

        return self._creg_bits(bitstring, flags), self._creg_bits(bitstring, self.qecc), bitstring, result

    def syndrome(self, pos: int):
        #syndrome measurement with the non-flag circuits. The outcome is read out classically (so it
        #can be fed into the correction sets) instead of being corrected in the circuit.
        self._round(pos=pos)

        self.qc.save_statevector(label="psi")
        self.gatecount += gates(self.qc)
        sim = AerSimulator(method="statevector", noise_model = self.noise_model)
        result = sim.run(self.qc, shots=1).result()
        bitstring = list(result.get_counts().keys())[0].replace(" ","")
        syndrome = self._creg_bits(bitstring, self.qecc)
        self.reset_qec(result)

        return self._split_syndrome(syndrome)

    def _split_syndrome(self, qecc_bits: str):
        #the qecc slice is printed qecc[15]..qecc[0], so int() puts qecc[i] on bit i
        value = int(qecc_bits, 2)
        return value & 0xFF, (value >> 8) & 0xFF        #(Z-Stabilizers, X-Stabilizers)

    def _parse_flags(self, flag_bits: str):
        #the flag events of this round: [((generator, xtype), flag pattern), ..]
        value = int(flag_bits, 2)
        out = []
        for (i, xtype), base in FLAGBIT_17Q.items():
            pattern = tuple((value >> (base+j)) & 1 for j in range(NFLAG_17Q[len(STABS_17Q[i])]))
            if any(pattern):
                out.append(((i, xtype), pattern))
        return out

    def _apply(self, err, pos: int):
        for q in range(17):
            if (err[0] >> q) & 1:
                self.qc.x(q+17*pos)
            if (err[1] >> q) & 1:
                self.qc.z(q+17*pos)

    def correct_17q(self, syndrome, pos: int):
        #E_min(s): the Z-Stabilizer syndrome asks for X corrections, the X-Stabilizer syndrome for
        #Z corrections. The code is self-dual, so one table serves both sectors.
        self._apply(correction([], syndrome), pos=pos)

    def flagcorrect(self, flagged, pos: int, syndrome, orders = ()):
        #apply any correction from the union of the E~_2^m(..., s) for m in orders, and fall back
        #to E_min(s) if all of them are empty
        self._apply(correction(flagged, syndrome, orders), pos=pos)

    def flagFTec(self, pos: int):
        #based on https://arxiv.org/pdf/1708.02246 , protocol on pages 7 and 8, "Flag 2-FTEC
        #Protocol" --> no need of postselection, but we need statevector simulation for classical
        #decoding. Every generator is measured with a 2-flag circuit, the weight-8 octagon needs
        #three flag qubits, the weight-4 faces one.
        def flag_round():
            flag_bits, qecc_bits, _, result = self._flag_round(pos=pos)
            self.reset_qec(result)
            return self._parse_flags(flag_bits), self._split_syndrome(qecc_bits)

        def nonflag_round():
            return self.syndrome(pos=pos)

        def correct(flagged, syndrome, orders):
            self.flagcorrect(flagged, pos, syndrome, orders)

        flag2_ftec_protocol(flag_round, nonflag_round, correct, clean_run=self.clean_run)

#######################################################################################################################################

    def ec(self, pos: int):
        #one round of error correction, either the Flag 2-FTEC Protocol or the plain qec() round
        if self.ftec:
            self.flagFTec(pos=pos)
        else:
            self.qec(pos=pos)

    def qec(self, pos: int):
        #plain (non fault tolerant) error correction: one round with the non-flag circuits followed
        #by E_min(s). Decoded classically like flagFTec, a lookup table over all 256 syndromes would
        #need 512 if_test blocks in the circuit.
        self.correct_17q(self.syndrome(pos=pos), pos=pos)

    def readout(self, pos: int, shots: int, p = 0):
        p_error = pauli_error([["X",p/2],["I",1-p],["Z",p/2]])
        p_error_2 = pauli_error([["XI",p/4],["IX",p/4],["II",1-p],["ZI",p/4],["IZ",p/4]])

        noise_model = NoiseModel()
        noise_model.add_all_qubit_quantum_error(p_error, ['x', "z", 'h', "s", "sdg", "id", "t", "tdg"])  # Apply to single-qubit gates
        noise_model.add_all_qubit_quantum_error(p_error_2, ['cx'])  # Apply to 2-qubit gates

        read = ClassicalRegister(17)
        self.qc.add_register(read)

        for i in range(17):
            self.qc.id(i+17*pos)
            self.qc.measure(i+17*pos, read[16-i])
        self.gatecount += gates(self.qc)
        sim = AerSimulator()
        job = sim.run(self.qc, shots=shots, noise_model=noise_model)
        result = job.result()
        counts = result.get_counts()

        bitstring = list(counts.keys())
        bitstring = [i.replace(" ","") for i in bitstring]

        hmm = list(counts.values())

        bits = [i[:17] for i in bitstring]                  #the 17 data qubits, bits[k] is qubit k

        for i in range(len(bits)):
            if bits[i] in CODE0_17Q:
                bits[i] = 0
            elif bits[i] in CODE1_17Q:
                bits[i] = 1
            else:
                if self.postselection:
                    bits[i] = "post"
                else:
                    if np.random.rand() < 0.5:
                        bits[i] = 0
                    else:
                        bits[i] = 1

        ones = 0
        zeros = 0
        post = 0
        preselected = 0

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


################################################ Tests ################################################
# Run with:  python bigsteane_ftec.py           (the Pauli frame checks, no qiskit simulation needed)
#            python bigsteane_ftec.py circuit   (additionally run the qiskit circuit end to end)
#
# The checks below verify, in this order:
#   1. every stabilizer measurement circuit is a 2-flag circuit          (Definition 6)
#   2. the code is a self-dual [[17,1,5]] CSS code and E_min is complete
#   3. the flag 2-FTEC condition holds for the chosen circuits           (Definition 10)
#   4. the protocol itself is fault tolerant                             (Definition 3, criterion 1)
# Point 4 drives the real flag2_ftec_protocol() and the real correction() with a Pauli frame model
# of the syndrome extraction rounds, which makes exhaustive fault injection cheap.

def _wt(xmask, zmask):
    return bin(xmask | zmask).count("1")

def _check_tflag_circuits(t = 2):
    #Definition 6: for any set of v <= t faults in C(P) resulting in an error E with
    #min(wt(E), wt(EP)) > v, the circuit must flag.
    bad = 0
    for w, seq in sorted(FLAG_SEQ_17Q.items()):
        full = (1 << w) - 1
        clean = _propagate(seq, w, {})
        assert clean == ((0,)*NFLAG_17Q[w], 0, 0, 0), "fault free run of the w=%d circuit is wrong" % w
        for v in range(1, t+1):
            for gates_ in itertools.combinations(range(len(seq)), v):
                for combo in itertools.product(_PAULIS2, repeat=v):
                    flags, xm, zm, _ = _propagate(seq, w, dict(zip(gates_, combo)))
                    if any(flags):
                        continue
                    if min(_wt(xm, zm), _wt(xm, zm ^ full)) > v:
                        bad += 1
        print("  w=%d circuit, %d flag qubit(s), %d CNOTs: %d-flag circuit"
              % (w, NFLAG_17Q[w], len(seq), t), "OK" if not bad else "VIOLATED (%d)" % bad)
    return bad

def _check_code():
    bad = 0
    #self-dual CSS: any two faces overlap in an even number of qubits
    for a in FACE_17Q:
        for b in FACE_17Q:
            if bin(a & b).count("1") % 2:
                bad += 1
    print("  generators commute:", "OK" if not bad else "VIOLATED")
    print("  weights:", [len(f) for f in STABS_17Q], "-> %d CNOTs per basis" % sum(len(f) for f in STABS_17Q))
    print("  stabilizer group size:", len(STABGROUP_17Q))
    d = min(bin(LMASK_17Q ^ s).count("1") for s in STABGROUP_17Q)
    print("  distance:", d, "OK" if d == 5 else "EXPECTED 5")
    bad += (d != 5) + (LMASK_17Q in STABGROUP_17Q)
    #E_min must be defined for every syndrome and must correct every error of weight <= 2 exactly
    print("  E_min covers all 256 syndromes:", "OK" if len(CORR_17Q) == 256 else "INCOMPLETE")
    bad += (len(CORR_17Q) != 256)
    miss = 0
    for w in (0, 1, 2):
        for qs in itertools.combinations(range(17), w):
            e = _mask(qs)
            if (e ^ _mask(CORR_17Q[_sector_syndrome(e)])) not in STABGROUP_17Q:
                miss += 1
    print("  E_min corrects every error of weight <= 2:", "OK" if not miss else "VIOLATED (%d)" % miss)
    print("  logical states disjoint and complete:",
          "OK" if len(CODE0_17Q) == len(CODE1_17Q) == 256 and not (CODE0_17Q & CODE1_17Q) else "BROKEN")
    return bad + miss

def _distinct_or_equivalent(errors):
    #E, E' in the set => s(E) != s(E') or E ~ E'
    seen = {}
    for e in errors:
        s = syndrome_of(e)
        ref = seen.get(s)
        if ref is None:
            seen[s] = e
        elif not equivalent(ref, e):
            return (ref, e, s)
    return None

def _check_flag2_condition():
    #Definition 10, with the flag error sets resolved by the flag pattern (Eq. 9 of Section 3.3).
    #The straightforward qubit order for the octagon fails condition 2, which is why ORDER_17Q[4]
    #is permuted.
    gens = [(i, xtype) for i in range(8) for xtype in (False, True)]
    bad1 = 0
    for a, b in itertools.combinations_with_replacement(gens, 2):
        #a == b covers the case that the same circuit flags in two different rounds
        for fa, sa in FLAGSETS_17Q[(a[0], a[1], 1)].items():
            for fb, sb in FLAGSETS_17Q[(b[0], b[1], 1)].items():
                if _distinct_or_equivalent((x[0] ^ y[0], x[1] ^ y[1]) for x in sa for y in sb):
                    bad1 += 1
    print("  condition 1, E_2(g_i, g_j):", "OK" if not bad1 else "VIOLATED (%d pairs)" % bad1)

    bad2 = []
    for g in gens:
        for f in set(FLAGSETS_17Q[(g[0], g[1], 1)]) | set(FLAGSETS_17Q[(g[0], g[1], 2)]):
            errors = list(FLAGSETS_17Q[(g[0], g[1], 2)].get(f, ()))
            for e in FLAGSETS_17Q[(g[0], g[1], 1)].get(f, ()):
                errors += [(e[0] ^ b[0], e[1] ^ b[1]) for b in WEIGHT1_17Q]
            if _distinct_or_equivalent(errors):
                bad2.append((g, f))
    print("  condition 2, E_2(g_i) u (E_1(g_i) x E_1):",
          "OK" if not bad2 else "VIOLATED %s" % bad2[:4])
    return bad1 + len(bad2)


#---------------------------------- Pauli frame model of the protocol ----------------------------------
def _coset_min_table():
    #minimum weight in a coset of the stabilizer group, keyed by (sector syndrome, logical parity)
    out = {}
    for w in range(10):
        for qs in itertools.combinations(range(17), w):
            m = _mask(qs)
            out.setdefault((_sector_syndrome(m), bin(m & LMASK_17Q).count("1") & 1), w)
        if len(out) == 512:
            break
    return out

def _residual_weight(residual, budget):
    #min over the stabilizer group of wt(residual * S). The sum of the two sector minima is an
    #upper bound (exact unless the X and the Z part share a qubit), so the expensive exact search
    #only runs when that bound alone would fail the budget.
    table = _COSET_MIN
    bound = (table[(_sector_syndrome(residual[0]), bin(residual[0] & LMASK_17Q).count("1") & 1)]
             + table[(_sector_syndrome(residual[1]), bin(residual[1] & LMASK_17Q).count("1") & 1)])
    if bound <= budget:
        return bound
    best = 18
    for sx in STABGROUP_17Q:
        ax = residual[0] ^ sx
        if bin(ax).count("1") >= best:
            continue
        for sz in STABGROUP_17Q:
            best = min(best, bin(ax | (residual[1] ^ sz)).count("1"))
    return best

def _sim_round(state, rnd, withflags, plan):
    #one syndrome extraction round over all 16 generators in the Pauli frame. state = [x, z] masks
    #of the data block, updated in place. plan = [(round, (generator, xtype), gate, pauli), ..]
    flagged, s_z, s_x = [], 0, 0
    for xtype in (False, True):
        for i in range(8):
            face, w = ORDER_17Q[i], len(ORDER_17Q[i])
            seq = FLAG_SEQ_17Q[w] if withflags else [("d", k) for k in range(w)]
            faults = {g: p for (r, cid, g, p) in plan
                      if r == rnd and cid == (i, xtype) and g < len(seq)}
            flags, lx, lz, mflip = _propagate(seq, w, faults)

            #the outcome is the parity of the data error on the face, flipped by a fault that
            #reached the measurement qubit
            bit = (bin(state[1 if xtype else 0] & FACE_17Q[i]).count("1") & 1) ^ mflip
            if xtype:
                s_x |= bit << i
            else:
                s_z |= bit << i

            #what the circuit deposited on the data block, visible to every later generator
            gx, gz = _lift(i, xtype, lx, lz)
            state[0] ^= gx
            state[1] ^= gz

            if withflags and any(flags):
                flagged.append(((i, xtype), flags))
    return flagged, (s_z, s_x)

def _run_protocol(input_err, plan, clean_run = False):
    #-> residual error after the protocol
    state = [input_err[0], input_err[1]]
    counter = [0]

    def flag_round():
        counter[0] += 1
        return _sim_round(state, counter[0]-1, True, plan)

    def nonflag_round():
        counter[0] += 1
        return _sim_round(state, counter[0]-1, False, plan)[1]

    def correct(flagged, syndrome, orders):
        err = correction(flagged, syndrome, orders)
        state[0] ^= err[0]
        state[1] ^= err[1]

    flag2_ftec_protocol(flag_round, nonflag_round, correct, clean_run=clean_run)
    return tuple(state)

_LOCATIONS = [(cid, g) for cid in FLAGBIT_17Q for g in range(max(len(s) for s in FLAG_SEQ_17Q.values()))]

def _report(name, fails, total):
    print("  %-4s %-46s %7d runs, %d violations"
          % ("OK" if not fails else "FAIL", name, total, len(fails)))
    for f in fails[:3]:
        print("         ", f)
    return len(fails)

def _check_fault_tolerance(trials = 20000, clean_run = False, seed = 0):
    #Definition 3, criterion 1: with w_in input errors and v faults, w_in + v <= 2, the output must
    #differ from the input codeword by an error of weight at most w_in + v.
    import random
    rng = random.Random(seed)
    total = 0

    fails, n = [], 0
    for xw in range(3):
        for zw in range(3-xw):
            for xq in itertools.combinations(range(17), xw):
                for zq in itertools.combinations(range(17), zw):
                    e = (_mask(xq), _mask(zq))
                    n += 1
                    if _residual_weight(_run_protocol(e, [], clean_run), xw+zw) > xw+zw:
                        fails.append(e)
    total += _report("input error of weight <= 2, no faults", fails, n)

    fails, n = [], 0
    for rnd in range(NMAX_17Q+1):
        for cid, g in _LOCATIONS:
            for p in _PAULIS2:
                n += 1
                if _residual_weight(_run_protocol((0, 0), [(rnd, cid, g, p)], clean_run), 1) > 1:
                    fails.append((rnd, cid, g, p))
    total += _report("1 fault, exhaustive", fails, n)

    weight1 = [e for e in WEIGHT1_17Q if e != (0, 0)]
    fails = []
    for _ in range(trials):
        plan = [(rng.randrange(NMAX_17Q+1), *rng.choice(_LOCATIONS), rng.choice(_PAULIS2))]
        e = rng.choice(weight1)
        if _residual_weight(_run_protocol(e, plan, clean_run), 2) > 2:
            fails.append((e, plan))
    total += _report("1 fault + weight-1 input error, random", fails, trials)

    for label, rounds in (("2 faults, random", NMAX_17Q+1), ("2 faults in the first 3 rounds", 3)):
        fails = []
        for _ in range(trials):
            plan = [(rng.randrange(rounds), *rng.choice(_LOCATIONS), rng.choice(_PAULIS2))
                    for _ in range(2)]
            if _residual_weight(_run_protocol((0, 0), plan, clean_run), 2) > 2:
                fails.append(plan)
        total += _report(label, fails, trials)

    #both faults inside the same circuit -> E_2(g_i), and one fault in each of two circuits
    #-> E_2(g_i, g_j). These are the two cases the correction sets are built for.
    for label, same in (("2 faults in the same circuit", True), ("2 faults in two circuits", False)):
        fails = []
        for _ in range(trials):
            if same:
                cid = rng.choice(list(FLAGBIT_17Q))
                seq = FLAG_SEQ_17Q[len(ORDER_17Q[cid[0]])]
                rnd = rng.randrange(3)
                plan = [(rnd, cid, g, rng.choice(_PAULIS2))
                        for g in rng.sample(range(len(seq)), 2)]
            else:
                plan = []
                for cid in rng.sample(list(FLAGBIT_17Q), 2):
                    seq = FLAG_SEQ_17Q[len(ORDER_17Q[cid[0]])]
                    plan.append((rng.randrange(2), cid, rng.randrange(len(seq)),
                                 rng.choice(_PAULIS2)))
            if _residual_weight(_run_protocol((0, 0), plan, clean_run), 2) > 2:
                fails.append(plan)
        total += _report(label, fails, trials)
    return total


def _check_circuit():
    #the same statements, but through the actual qiskit circuit
    bad = 0
    for label, ops in (("|0_L>", []), ("X_L", ["x"]), ("H_L H_L", ["h", "h"]),
                       ("H_L X_L H_L", ["h", "x", "h"]), ("S_L Sdg_L", ["s", "sdg"]),
                       ("H_L S_L S_L H_L", ["h", "s", "s", "h"])):
        q = BigSteane17q(1, 0, noise=0.0)
        for o in ops:
            getattr(q, o)(pos=0)
        q.readout(pos=0, shots=100, p=0)
        expect = 1.0 if label in ("|0_L>", "H_L H_L", "H_L X_L H_L", "S_L Sdg_L") else 0.0
        ok = abs(q.zeros - expect) < 1e-9
        bad += not ok
        print("  %-16s zeros=%.2f ones=%.2f post=%.2f %s"
              % (label, q.zeros, q.ones, q.post, "OK" if ok else "UNEXPECTED"))

    q = BigSteane17q(1, 0, noise=0.0)
    s = q.syndrome(pos=0)
    print("  syndrome of a clean |0_L>:", s, "OK" if s == (0, 0) else "UNEXPECTED")
    bad += s != (0, 0)
    for inj in ((3,), (3, 11), (0, 5, 9)):
        q = BigSteane17q(1, 0, noise=0.0)
        for qb in inj:
            q.qc.x(qb)
        s, exp = q.syndrome(pos=0), _sector_syndrome(_mask(inj))
        print("  X on %-10s -> s_Z=%3d, expected %3d %s"
              % (str(inj), s[0], exp, "OK" if s[0] == exp else "MISMATCH"))
        bad += s[0] != exp

    for inj_x, inj_z in (((), ()), ((4,), ()), ((), (12,)), ((2, 15), ()), ((), (6, 9)), ((7,), (13,))):
        q = BigSteane17q(1, 0, noise=0.0)
        for qb in inj_x:
            q.qc.x(qb)
        for qb in inj_z:
            q.qc.z(qb)
        q.flagFTec(pos=0)
        q.readout(pos=0, shots=1, p=0)
        print("  flagFTec on X%s Z%s -> zeros=%s ones=%s post=%s"
              % (list(inj_x), list(inj_z), q.zeros, q.ones, q.post),
              "OK" if q.zeros == 1 else "UNEXPECTED")
        bad += q.zeros != 1
    return bad


if __name__ == "__main__":
    import sys
    _COSET_MIN = _coset_min_table()

    print("\n[1] 2-flag circuits (Definition 6)")
    bad = _check_tflag_circuits()
    print("\n[2] the [[17,1,5]] code and E_min")
    bad += _check_code()
    print("\n[3] flag 2-FTEC condition (Definition 10)")
    bad += _check_flag2_condition()

    for clean_run in (False, True):
        print("\n[4] fault tolerance of the protocol (Definition 3), clean_run = %s%s"
              % (clean_run, "   <- the protocol as printed in the paper" if not clean_run else ""))
        n = _check_fault_tolerance(trials=20000, clean_run=clean_run)
        if not clean_run:
            print("      (a handful of violations here is expected, see the 1.Case in"
                  " flag2_ftec_protocol)")
        else:
            bad += n

    if "circuit" in sys.argv:
        print("\n[5] the qiskit circuit end to end (slow)")
        bad += _check_circuit()

    print("\n=> %s" % ("all checks passed" if not bad else "%d PROBLEMS" % bad))
