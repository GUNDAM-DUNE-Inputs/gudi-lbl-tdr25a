#!/bin/bash
#SBATCH --job-name=GundamGrid
#SBATCH --output=slurm_files/grid_log/GundamGrid_%A_%a.out
#SBATCH --error=slurm_files/grid_log/GundamGrid_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=Flynn.Y.Guo@stonybrook.edu
#SBATCH --mem=50G
#SBATCH -c 16
#SBATCH --gres=gpu
#SBATCH --array=0-62
#SBATCH --time=2:00:00

TEMPLATE_CFG="./config_DUNE.yaml"
YAML_DIR="./outputs/profiling_yaml/MaCh3_v3"
JSON_DIR="./outputs/profiling_json_files/v2"

source /home/rrazakami/workspace/ROOT/root_binary/bin/thisroot.sh
source /home/fyguo/GUNDAM/Install/gundam/setup.sh

echo "[$(date '+%F %T')] Host: $(hostname)"
echo "Task ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} started"

mapfile -t YAML_LIST < <(ls "${YAML_DIR}"/osc_*.yaml | sort)
OVERRIDE_YAML="${YAML_LIST[${SLURM_ARRAY_TASK_ID}]}"
STEM=$(basename "${OVERRIDE_YAML}" .yaml)

echo "Processing override: ${STEM}"

# ---------- fixed block -------------------------------------------
mkdir -p "${JSON_DIR}"
gundamConfigUnfolder -c "${TEMPLATE_CFG}" \
                     -of "${OVERRIDE_YAML}" \
                     -o  "${JSON_DIR}/config_DUNE_With_${STEM}.json"
# ------------------------------------------------------------------

echo "[$(date '+%F %T')] Task ${SLURM_ARRAY_TASK_ID} finished"
