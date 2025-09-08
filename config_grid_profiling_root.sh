#!/bin/bash
#
# 1) count JSON files
# NUM=$(ls outputs/profiling_json_files/config_DUNE_With_*.json | wc -l)
# echo "$NUM JSON files found"
# 2) submit one array task per file
# sbatch --array=0-$((NUM-1)) config_grid_profiling_root.sh
#
#SBATCH --job-name=GundamFits
#SBATCH --output=slurm_files/fit_log/GundamFits_%A_%a.out
#SBATCH --error=slurm_files/fit_log/GundamFits_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=Flynn.Y.Guo@stonybrook.edu
#SBATCH --mem=50G
#SBATCH -c 16
#SBATCH --gres=gpu

# ------------------------------------------------------------------
JSON_DIR="$PWD/outputs/profiling_json_files/v2"    # where the JSONs live
# ------------------------------------------------------------------

# ──────────────  ENVIRONMENT  ──────────────────────────────────────
source /home/rrazakami/workspace/ROOT/root_binary/bin/thisroot.sh
source /home/fyguo/GUNDAM/Install/gundam/setup.sh

export OA_INPUTS=/storage/shared/DUNE/OA-inputs/lbl/gudi-inputs/v5
export LD_LIBRARY_PATH="/home/fyguo/Work/Projects/GUNDAM/GundamInputDUNE_LBL/Devel/gundamOscAnaTools/resources/TabulateNuOscillator/build-x86_64/_deps/oscprob-src/lib/:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="/home/fyguo/Work/Projects/GUNDAM/GundamInputDUNE_LBL/Devel/gundamOscAnaTools/resources/TabulateNuOscillator/build-x86_64/lib/:$LD_LIBRARY_PATH"
export NUOSCILLATOR_ROOT=/home/fyguo/Work/Projects/GUNDAM/GundamInputDUNE_LBL/Devel/gundamOscAnaTools/resources/TabulateNuOscillator/build-x86_64
source ${NUOSCILLATOR_ROOT}/bin/setup.NuOscillator.sh
# ───────────────────────────────────────────────────────────────────

echo "[$(date '+%F %T')] host=$(hostname) task=$SLURM_ARRAY_TASK_ID"

# ---------- build the list of JSON files ---------------------------
mapfile -t JSON_LIST < <(ls "$JSON_DIR"/config_DUNE_With_*.json | sort)

TOTAL=${#JSON_LIST[@]}
if (( TOTAL == 0 )); then
  echo "ERROR: no JSON files in $JSON_DIR"; exit 1
fi
if (( SLURM_ARRAY_TASK_ID >= TOTAL )); then
  echo "ERROR: task index ${SLURM_ARRAY_TASK_ID} ≥ total ${TOTAL}"; exit 1
fi

JSON_FILE="${JSON_LIST[$SLURM_ARRAY_TASK_ID]}"
echo "Running: gundamFitter -c $JSON_FILE -a --gpu -t 16"
# gundamFitter -c "$JSON_FILE" -a --gpu -t 16
# gundamFitter -d -a -c "$JSON_FILE" --inject-toy-parameter ./injector/parameterInjector1.yaml -t 16 --toy 0
gundamFitter -a -c "$JSON_FILE" --inject-toy-parameter ./injector/parameterInjector1.yaml -t 16 --toy 0
echo "[$(date '+%F %T')] task $SLURM_ARRAY_TASK_ID finished"
