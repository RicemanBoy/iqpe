"""
Spacetime MWPM decoding of a repetition-code QEC block (Qiskit -> PyMatching).

Input : list of stabilizer-readout bitstrings, one per QEC cycle.
        d=5 repetition code -> 4 checks -> bitstrings of length 4.
Output: list of data-qubit indices where an X correction must be applied.
"""

import numpy as np
import pymatching

# Matching graphs are cached: building one per block would dominate runtime.
_MATCHER_CACHE = {}


def _node(t, j, n_checks):
    """Detector index of check j in round t (row-major, matches det.reshape(-1))."""
    return t * n_checks + j


def _build_matcher(n_checks, T, p_data, p_meas):
    """Spacetime matching graph: space edges = data errors, time edges = readout errors."""
    key = (n_checks, T, p_data, p_meas)
    if key in _MATCHER_CACHE:
        return _MATCHER_CACHE[key]

    n_data = n_checks + 1
    w_data = np.log((1.0 - p_data) / p_data)
    w_meas = np.log((1.0 - p_meas) / p_meas)

    m = pymatching.Matching()

    for t in range(T):
        # --- space-like edges: an X error on data qubit q flips checks q-1 and q
        for q in range(n_data):
            left, right = q - 1, q
            if left < 0:                                    # leftmost qubit -> boundary
                m.add_boundary_edge(_node(t, right, n_checks), fault_ids={q},
                                    weight=w_data, error_probability=p_data)
            elif right > n_checks - 1:                      # rightmost qubit -> boundary
                m.add_boundary_edge(_node(t, left, n_checks), fault_ids={q},
                                    weight=w_data, error_probability=p_data)
            else:                                           # bulk qubit
                m.add_edge(_node(t, left, n_checks), _node(t, right, n_checks),
                           fault_ids={q}, weight=w_data, error_probability=p_data)

        # --- time-like edges: a readout error on check j flips det[t][j] and det[t+1][j]
        if t < T - 1:
            for j in range(n_checks):
                m.add_edge(_node(t, j, n_checks), _node(t + 1, j, n_checks),
                           fault_ids=set(), weight=w_meas, error_probability=p_meas)

    _MATCHER_CACHE[key] = m
    return m


def syndromes_to_detection_events(a, boundary=None, qiskit_endian=True):
    """
    Raw stabilizer outcomes -> detection events.

    det[0] = syn[0] XOR boundary      (boundary = expected syndrome after state prep)
    det[t] = syn[t] XOR syn[t-1]
    """
    T = len(a)
    n_checks = len(a[0])
    if any(len(s) != n_checks for s in a):
        raise ValueError("All bitstrings in the block must have the same length.")

    def parse(bstr):
        bits = [int(b) for b in bstr]
        if qiskit_endian:
            bits = bits[::-1]        # Qiskit prints classical bit 0 rightmost
        return np.array(bits, dtype=np.uint8)

    syn = np.stack([parse(s) for s in a])            # (T, n_checks)

    if boundary is None:
        boundary = np.zeros(n_checks, dtype=np.uint8)
    else:
        boundary = parse(boundary) if isinstance(boundary, str) \
            else np.asarray(boundary, dtype=np.uint8)

    det = syn.copy()
    det[0] ^= boundary
    det[1:] ^= syn[:-1]
    return det


def data_readout_to_syndrome(m_bits, qiskit_endian=True):
    """
    Perfect final syndrome from a destructive data-qubit measurement.

    Repetition code: check j is the parity of data qubits j and j+1, so a
    5-bit readout yields a 4-bit syndrome. Append the result to the block as
    one extra item to close the time boundary at the end of the algorithm.
    """
    bits = [int(b) for b in m_bits]
    if qiskit_endian:
        bits = bits[::-1]                 # bits[q] = data qubit q
    m = np.asarray(bits, dtype=np.uint8)
    s = m[:-1] ^ m[1:]                    # length n_data - 1 = n_checks
    if qiskit_endian:
        s = s[::-1]
    return "".join(str(int(b)) for b in s)


def decode_qec_block(a, p=1e-3, p_meas=None, boundary=None,
                     qiskit_endian=True, return_array=False):
    """
    Decode one QEC block and return the net physical correction.

    Parameters
    ----------
    a : list[str]
        One stabilizer-readout bitstring per QEC cycle, e.g. ["0000", "0100", ...].
    p : float
        Data/gate fault probability per round (sets space-edge weights).
    p_meas : float or None
        Readout fault probability (time-edge weights). Defaults to p, which is the
        right choice when all gates and qubits share one error rate.
    boundary : str | array | None
        Expected syndrome at the start of the block. None -> all zeros, i.e. the
        block begins right after a clean logical state preparation.
    return_array : bool
        If True, also return the length-n_data 0/1 correction vector.

    Returns
    -------
    list[int]  (and optionally np.ndarray)
        Data-qubit indices that need an X gate. Empty list = no correction.
    """
    if not a:
        return ([], np.zeros(0, dtype=np.uint8)) if return_array else []

    if p_meas is None:
        p_meas = p

    det = syndromes_to_detection_events(a, boundary=boundary,
                                        qiskit_endian=qiskit_endian)
    T, n_checks = det.shape

    matcher = _build_matcher(n_checks, T, p, p_meas)

    # reshape(-1) flattens row-major -> exactly the t*n_checks + j order used above
    correction = matcher.decode(det.reshape(-1)).astype(np.uint8)

    qubits = np.flatnonzero(correction).tolist()
    return (qubits, correction) if return_array else qubits


# ----------------------------------------------------------------------
if __name__ == "__main__":
    # d=5 repetition code, 4 checks, 5 cycles.

    # (1) clean block -> no correction
    print(decode_qec_block(["0000"] * 5))

    # (2) persistent syndrome on checks 1,2 from cycle 1 on -> data qubit 2 flipped
    #     (Qiskit endianness: "0110" -> bits [0,1,1,0] for checks 0..3)
    print(decode_qec_block(["0000", "0110", "0110", "0110", "0110"]))

    # (3) transient blip on check 0 in cycle 2 only -> readout error, no correction
    print(decode_qec_block(["0000", "0000", "0001", "0000", "0000"]))

    # (4) boundary error: leftmost qubit flipped -> only check 0 fires, persistently
    print(decode_qec_block(["0001", "0001", "0001", "0001", "0001"]))

    # (5) two separate data errors, net correction on both qubits
    print(decode_qec_block(["0000", "0110", "0110", "1100", "1100"]))