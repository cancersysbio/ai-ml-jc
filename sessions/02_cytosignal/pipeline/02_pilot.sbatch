#!/bin/bash
# CytoSignal pilot run (~8000 cells) -- validates the pipeline end-to-end and
# gives a real wall-clock number to extrapolate the full run's --time from.
#
# Submit:   sbatch 02_pilot.sbatch
# Monitor:  squeue -u $USER   |   tail -f logs/pilot_*.out
#
# Edit --partition below for your cluster before submitting.
#SBATCH --job-name=cytosignal-breast-pilot
#SBATCH --partition=<your_cpu_partition>
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=logs/pilot_%j.out
#SBATCH --error=logs/pilot_%j.err
set -euo pipefail

mkdir -p logs
echo "host: $(hostname)"
echo "=== start: $(date) ==="

# --- 1. export from the grafiti python venv ---
source "${GRAFITI_VENV:-/path/to/your/grafiti/venv}/bin/activate"
python 00_export_breast_atera.py --output-dir ./exported_pilot --subset 8000
deactivate

# --- 2. run CytoSignal in R (conda env with cytosignal installed -- see README) ---
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CYTOSIGNAL_ENV:-cytosignal-r}"
echo "R env: $(conda info --envs | grep '\*')"; echo "R: $(command -v R)"

export DATA_DIR=./exported_pilot
export OUTPUT_DIR=./cytosignal_output_pilot
export N_CORES=4
Rscript 01_cytosignal_breast_atera.R

echo "=== end: $(date) ==="
