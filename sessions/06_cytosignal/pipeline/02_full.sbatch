#!/bin/bash
# CytoSignal full run on the entire breast_atera dataset (~400K cells).
#
# *** Do not submit this until the pilot (02_pilot.sbatch) has completed
# and you've read its actual wall-clock/peak-mem (`sacct -j <jobid> --format=
# Elapsed,MaxRSS`). Scale --time and --mem from those numbers -- the defaults
# below are placeholders, not a benchmark. ***
#
# Submit:   sbatch 02_full.sbatch
# Monitor:  squeue -u $USER   |   tail -f logs/full_*.out
#
# Edit --partition below for your cluster before submitting.
#SBATCH --job-name=cytosignal-breast-full
#SBATCH --partition=<your_cpu_partition>
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=180G
#SBATCH --output=logs/full_%j.out
#SBATCH --error=logs/full_%j.err
set -euo pipefail

mkdir -p logs
echo "host: $(hostname)"
echo "=== start: $(date) ==="

# --- 1. export from the grafiti python venv (full dataset, no --subset) ---
source "${GRAFITI_VENV:-/path/to/your/grafiti/venv}/bin/activate"
python 00_export_breast_atera.py --output-dir ./exported_full
deactivate

# --- 2. run CytoSignal in R (conda env with cytosignal installed -- see README) ---
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CYTOSIGNAL_ENV:-cytosignal-r}"
echo "R env: $(conda info --envs | grep '\*')"; echo "R: $(command -v R)"

export DATA_DIR=./exported_full
export OUTPUT_DIR=./cytosignal_output_full
export N_CORES=16
Rscript 01_cytosignal_breast_atera.R

echo "=== end: $(date) ==="
