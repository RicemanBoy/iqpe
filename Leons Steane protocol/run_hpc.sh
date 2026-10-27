#!/bin/bash
# SLURM array job for the full-iQPE sweep (sweep.py, pure numpy, no Qiskit).
#
#   sbatch run_hpc.sh                     # 32 tasks x 16 cores
#   python sweep.py --summary             # afterwards (needs matplotlib): table + results/sweep_full.png
#
# Each task takes every 32nd job of the plan in sweep.py (PLAN) and skips
# jobs whose result file already exists, so the array can be resubmitted
# after a time-out or to continue a partial run.  Copy results/ back into
# FAUbox (phase5_full_iqpe/results/) to merge with the local results.
#
# Environment: python >= 3.10 with numpy (matplotlib only for --summary).
# Adapt partition / time / module lines to the cluster.
#SBATCH --job-name=iqpe-full
#SBATCH --array=0-31
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%A_%a.out

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
mkdir -p logs results

# module load python   # adapt
# source activate <env>

# one numpy thread per worker process (several BLAS threads per process slow everything down)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

python -u sweep.py --task "$SLURM_ARRAY_TASK_ID" --ntasks "${SLURM_ARRAY_TASK_COUNT:-32}" \
    --procs "${SLURM_CPUS_PER_TASK:-16}"
