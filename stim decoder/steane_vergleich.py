"""
Steane-Code [[7,1,3]] Memory-Experiment mit Stim + Sinter.
Vergleicht zwei Decoder:
  - pymatching  (Minimum-Weight Perfect Matching, Sparse Blossom)
  - tesseract   (Most-Likely-Error-Decoder von Google, kann Hyperkanten)

Installation:
  pip install stim sinter pymatching tesseract-decoder matplotlib
"""

import stim
import sinter
import matplotlib.pyplot as plt
from tesseract_decoder import make_tesseract_sinter_decoders_dict

# Die drei Paritätschecks des [7,4,3]-Hamming-Codes (0-indiziert).
# Der Steane-Code nutzt sie doppelt: einmal als X-, einmal als Z-Stabilisatoren.
CHECKS = [(3, 4, 5, 6), (1, 2, 5, 6), (0, 2, 4, 6)]


def steane_memory_z(rounds: int, p: float) -> stim.Circuit:
    """Phänomenologisches Memory-Z-Experiment für den Steane-Code.

    - 7 Datenqubits, Start in |0...0>
    - pro Runde: Depolarisierungsrauschen (Stärke p) auf allen Datenqubits,
      dann verrauschte Messung (Fehlerrate p) aller 6 Stabilisatoren via MPP
    - Detektoren vergleichen Stabilisator-Ergebnisse aufeinanderfolgender Runden
    - Observable: logisches Z = Z auf allen 7 Qubits
    """
    c = stim.Circuit()
    c.append("R", range(7))

    def measure_stabilizers(basis: str):
        for chk in CHECKS:
            targets = []
            for q in chk:
                if targets:
                    targets.append(stim.target_combiner())
                targets.append(stim.target_x(q) if basis == "X" else stim.target_z(q))
            c.append("MPP", targets, p)  # verrauschte Stabilisatormessung

    for r in range(rounds):
        c.append("DEPOLARIZE1", range(7), p)
        measure_stabilizers("Z")
        measure_stabilizers("X")

        # Z-Stabilisatoren: in Runde 0 deterministisch (+1 wegen |0>^7),
        # danach Vergleich mit der Vorrunde.
        for i in range(3):
            if r == 0:
                c.append("DETECTOR", [stim.target_rec(-6 + i)])
            else:
                c.append("DETECTOR", [stim.target_rec(-6 + i), stim.target_rec(-12 + i)])
        # X-Stabilisatoren: erste Messung ist zufällig -> erst ab Runde 1 vergleichen.
        if r > 0:
            for i in range(3):
                c.append("DETECTOR", [stim.target_rec(-3 + i), stim.target_rec(-9 + i)])

    # Finale Messung aller Datenqubits in der Z-Basis.
    c.append("M", range(7))
    for i, chk in enumerate(CHECKS):
        recs = [stim.target_rec(q - 7) for q in chk]
        recs.append(stim.target_rec(-7 - 6 + i))  # letzte Z-Stabilisator-Messung
        c.append("DETECTOR", recs)
    c.append("OBSERVABLE_INCLUDE", [stim.target_rec(q - 7) for q in range(7)], 0)
    return c


def main():
    tasks = [
        sinter.Task(
            circuit=steane_memory_z(rounds=3, p=p),
            json_metadata={"p": p, "rounds": 3, "d": 3},
        )
        for p in [0.002, 0.005, 0.01, 0.02, 0.05]
    ]

    stats = sinter.collect(
        num_workers=4,
        tasks=tasks,
        decoders=["pymatching", "tesseract"],
        custom_decoders=make_tesseract_sinter_decoders_dict(),
        max_shots=100_000,
        max_errors=1_000,
        print_progress=True,
    )

    # Ergebnisse als Tabelle ausgeben
    print(f"\n{'Decoder':<12} {'p':>7} {'Shots':>8} {'Fehler':>7} {'log. Fehlerrate':>16}")
    for s in sorted(stats, key=lambda s: (s.decoder, s.json_metadata["p"])):
        print(f"{s.decoder:<12} {s.json_metadata['p']:>7} {s.shots:>8} "
              f"{s.errors:>7} {s.errors / s.shots:>15.4%}")

    # Plot: logische Fehlerrate vs. physikalische Fehlerrate
    fig, ax = plt.subplots(figsize=(7, 5))
    sinter.plot_error_rate(
        ax=ax,
        stats=stats,
        x_func=lambda s: s.json_metadata["p"],
        group_func=lambda s: s.decoder,
    )
    ax.loglog()
    ax.plot([1e-3, 1e-1], [1e-3, 1e-1], "k--", alpha=0.4, label="p_log = p (Referenz)")
    ax.set_xlabel("Physikalische Fehlerrate p")
    ax.set_ylabel("Logische Fehlerrate")
    ax.set_title("Steane-Code [[7,1,3]], 3 Runden: PyMatching vs. Tesseract")
    ax.grid(which="both", alpha=0.3)
    ax.legend()
    fig.savefig("stim decoder/steane_pymatching_vs_tesseract.png", dpi=150, bbox_inches="tight")
    print("\nPlot gespeichert: steane_pymatching_vs_tesseract.png")


if __name__ == "__main__":
    main()
