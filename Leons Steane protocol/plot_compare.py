"""Full iQPE (this folder, results/sweep_full.json) next to Ramsey iQPE
(phase4_qiskit_hpc/results/sweep.json): Steane code vs bare qubits.

  python plot_compare.py      # writes results/compare.png
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
INK2, MUTED, SURFACE, GRID = "#52514e", "#8a8984", "#fcfcfb", "#e6e5e1"
BLUE, LIGHT = "#0d366b", "#6b8fc4"
plt.rcParams.update({"figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
                     "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.grid": True, "grid.color": GRID, "axes.spines.top": False, "axes.spines.right": False,
                     "legend.frameon": False, "lines.linewidth": 2, "lines.markersize": 6.5, "font.size": 10})


def series(rows, key, val):
    rr = sorted((r for r in rows if r[key] == val), key=lambda r: r["p"])
    return [r["p"] for r in rr], [r["mean_dphi"] for r in rr], [r["err"] for r in rr]


def breakeven(ps, a, b):
    """log-log interpolation of the crossing of a (logical) and b (bare)."""
    lr = np.log(np.array(a) / np.array(b))
    for i in range(len(ps) - 1):
        if lr[i] < 0 <= lr[i + 1]:
            x0, x1 = np.log(ps[i]), np.log(ps[i + 1])
            return float(np.exp(x0 + (x1 - x0) * (-lr[i]) / (lr[i + 1] - lr[i])))
    return None


def main():
    full = json.load(open(HERE / "results" / "sweep_full.json"))
    ram = json.load(open(HERE.parent / "phase4_qiskit_hpc" / "results" / "sweep.json"))
    fig, ax = plt.subplots(figsize=(6.8, 5.0))
    p, y, e = series(full, "run", "bare")
    ax.errorbar(p, y, yerr=e, ls=":", marker="o", color=MUTED, mec=SURFACE, capsize=3,
                label="full iQPE, 2 bare qubits")
    pf, yf, ef = series(full, "run", "full")
    ax.errorbar(pf, yf, yerr=ef, ls="-", marker="D", color=BLUE, mec=SURFACE, capsize=3,
                label="full iQPE, 2 Steane logical qubits")
    be_full = breakeven(pf, yf, [dict(zip(p, y))[q] for q in pf])
    pr, yr, er = series(ram, "code", "bare")
    ax.errorbar(pr, yr, yerr=er, ls=":", marker="s", color="#b9b8b3", mec=SURFACE, capsize=3, ms=5,
                label="Ramsey iQPE, bare qubit (Qiskit)")
    prs, yrs, ers = series(ram, "code", "steane")
    ax.errorbar(prs, yrs, yerr=ers, ls="--", marker="s", color=LIGHT, mec=SURFACE, capsize=3, ms=5,
                label="Ramsey iQPE, Steane (Qiskit)")
    be_ram = breakeven(prs, yrs, [dict(zip(pr, yr))[q] for q in prs])
    for be, col in ((be_full, BLUE), (be_ram, LIGHT)):
        if be:
            ax.axvline(be, color=col, lw=1, ls="-.", alpha=0.6)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("physical error rate p")
    ax.set_ylabel("mean phase error  Δφ")
    ax.set_title("FT iQPE on the Steane code with real T gates", loc="left", fontsize=11, fontweight="bold")
    ax.legend(fontsize=8, loc="upper left")
    fig.text(0.01, -0.02,
             f"Break-even (dash-dotted): full iQPE p ≈ {be_full:.1e}, Ramsey p ≈ {be_ram:.1e}.  3 iterations; full iQPE: "
             "angles j/8, j = 2, 4, 6,\nancilla + system block, controlled-Z / controlled-S with teleported T gates; "
             "Ramsey: one block, j = 1..7.\nSD6 circuit noise on every gate, reset and measurement, no idle noise; "
             "flag QEC after every logical gate.", fontsize=7, color=INK2, va="top")
    fig.tight_layout()
    fig.savefig(HERE / "results" / "compare.png", dpi=150, bbox_inches="tight")
    print(f"break-even: full iQPE {be_full:.2e}, Ramsey {be_ram:.2e}; wrote results/compare.png")


if __name__ == "__main__":
    main()
