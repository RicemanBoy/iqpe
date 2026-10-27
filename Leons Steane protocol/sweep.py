"""Resumable local sweep for iqpe_full.py: full iQPE on the Steane code
(angles j = 2, 4, 6) and the bare two-qubit iQPE, several p, in parallel.

  python sweep.py --procs 8            # runs until all jobs are done
  python sweep.py --summary            # table + results/sweep_full.png

One JSON per (run, j, p, seed) in results/; existing files are skipped, so
the sweep can be stopped and restarted at any time.  Jobs are ordered so
that every p first gets a balanced set of chunks over all angles.
"""

import argparse
import json
import time
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path

import numpy as np

import iqpe_full as m

HERE = Path(__file__).resolve().parent
OUT = HERE / "results"
ANGLES = (2, 4, 6)
CHUNK = 512
# (p, chunks of CHUNK shots per angle), in priority order
PLAN = [(2e-3, 32), (1e-3, 64), (5e-3, 16), (5e-4, 256), (2e-4, 512), (1e-4, 768)]
BARE_SHOTS = 2_000_000


def fname(run, j, p, seed):
    return OUT / f"{run}_j{j}_p{p:.1e}_s{seed}.json"


def jobs():
    js = []
    for p, chunks in PLAN:
        for j in ANGLES:
            js.append(("bare", j, p, 1, BARE_SHOTS))
        for seed in range(1, chunks + 1):
            for j in ANGLES:
                js.append(("full", j, p, seed, CHUNK))
    return js


def work(job):
    run, j, p, seed, shots = job
    f = fname(run, j, p, seed)
    if f.exists():
        return None
    m.DTYPE = np.complex64
    t0 = time.time()
    s = 1_000_003 * seed + 1009 * j + int(round(p * 1e7))
    fl = None
    m.SPLIT_STATS.update(splits=0, non_product=0, max_residual=0.0)
    if run == "full":
        d, fl = m.run_full_split(j, p, shots, s, check_system=True)
    else:
        d = m.run_bare_full(j, p, shots, s)
    r = {"run": run, "j": j, "phi": j / 8, "p": p, "shots": shots, "seed": seed,
         "sum_dphi": float(d.sum()), "n_wrong": int((d > 0).sum()),
         "system_flipped": None if fl is None else int(fl.sum()), "seconds": time.time() - t0,
         "non_product": m.SPLIT_STATS["non_product"], "max_split_residual": m.SPLIT_STATS["max_residual"]}
    tmp = f.with_suffix(".tmp")
    json.dump(r, open(tmp, "w"))
    tmp.rename(f)
    return r


def summary():
    agg = defaultdict(lambda: defaultdict(lambda: [0, 0.0, 0, 0]))
    for f in list(OUT.glob("full_*.json")) + list(OUT.glob("bare_*.json")):
        r = json.load(open(f))
        a = agg[(r["run"], r["p"])][r["j"]]
        a[0] += r["shots"]; a[1] += r["sum_dphi"]; a[2] += r["n_wrong"]; a[3] += r["system_flipped"] or 0
    rows = []
    for (run, p), per in sorted(agg.items()):
        means = [v[1] / v[0] for v in per.values()]
        nw = sum(v[2] for v in per.values())
        mean = float(np.mean(means))
        rows.append({"run": run, "p": p, "angles": sorted(per), "shots": sum(v[0] for v in per.values()),
                     "n_wrong": nw, "mean_dphi": mean, "err": mean / np.sqrt(max(nw, 1)),
                     "system_flipped": sum(v[3] for v in per.values()),
                     "per_angle": {j: v[1] / v[0] for j, v in sorted(per.items())}})
    json.dump(rows, open(OUT / "sweep_full.json", "w"), indent=1)
    print(f"{'run':5s} {'p':>8s} {'shots':>8s} {'wrong':>6s} {'mean dphi':>20s}  per angle")
    for r in rows:
        print(f"{r['run']:5s} {r['p']:8.1e} {r['shots']:8d} {r['n_wrong']:6d} {r['mean_dphi']:10.3e} +- {r['err']:.1e}  "
              + " ".join(f"j{j}:{v:.2e}" for j, v in r["per_angle"].items())
              + (f"  S flipped at end: {r['system_flipped']}" if r["run"] == "full" else ""))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    for run, lab, kw in (("bare", "bare physical qubits (2-qubit iQPE)", dict(color="#8a8984", marker="o", ls=":")),
                         ("full", "Steane code, FT (2 logical qubits)", dict(color="#0d366b", marker="D", ls="-"))):
        rr = sorted((r for r in rows if r["run"] == run), key=lambda r: r["p"])
        if rr:
            ax.errorbar([r["p"] for r in rr], [r["mean_dphi"] for r in rr], yerr=[r["err"] for r in rr],
                        capsize=3, label=lab, **kw)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("physical error rate p")
    ax.set_ylabel("mean phase error  Δφ  (j = 2, 4, 6)")
    ax.set_title("Full iQPE (ancilla + system), real controlled-S / CZ", loc="left", fontsize=10)
    ax.grid(True, color="#e6e5e1")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "sweep_full.png", dpi=150)
    print(f"wrote {OUT / 'sweep_full.png'}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--procs", type=int, default=8)
    ap.add_argument("--summary", action="store_true")
    ap.add_argument("--task", type=int, default=None,
                    help="HPC array mode: this task takes every ntasks-th job of the plan")
    ap.add_argument("--ntasks", type=int, default=1)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if a.summary:
        summary()
        return
    plan = jobs() if a.task is None else jobs()[a.task::a.ntasks]
    todo = [jb for jb in plan if not fname(*jb[:4]).exists()]
    print(f"{len(todo)} jobs to do" + ("" if a.task is None else f" (task {a.task}/{a.ntasks})"), flush=True)
    t0 = time.time()
    with Pool(a.procs) as pool:
        for i, r in enumerate(pool.imap(work, todo, chunksize=1)):
            if r is not None:
                print(f"[{(time.time() - t0) / 3600:5.2f} h] {i + 1}/{len(todo)} {r['run']} j={r['j']} p={r['p']:.0e} "
                      f"seed={r['seed']}: dphi {r['sum_dphi'] / r['shots']:.2e}, wrong {r['n_wrong']}, "
                      f"S flipped {r['system_flipped']} [{r['seconds']:.0f}s]", flush=True)
            if a.task is None and (i + 1) % 50 == 0:
                summary()
    if a.task is None:
        summary()


if __name__ == "__main__":
    main()
