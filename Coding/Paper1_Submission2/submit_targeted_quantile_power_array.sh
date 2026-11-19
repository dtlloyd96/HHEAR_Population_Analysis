#!/bin/bash
#SBATCH -c 1
#SBATCH -t 0-20:0
#SBATCH -p medium
#SBATCH --mem=8G
#SBATCH -J targeted_qpower

set -euo pipefail

module load gcc/14.2.0 R/4.4.2

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${HHEAR_ROOT:-$HOME/HHEAR/HHEAR_Quantreg}"
CODING_DIR="${PROJECT_ROOT}/Coding"
R_SCRIPT="${CODING_DIR}/Paper1_Submission2/Targeted_Quantile_Power_Array.R"
INPUT_DIR="${PROJECT_ROOT}/Input"
OUT_DIR="${PROJECT_ROOT}/Output/Power_Output"

NSIM=200
ALPHA=0.05
MIN_N=60
SEED=20260504
TAUS="0.05,0.25,0.5,0.75,0.95"
USE_LOG1P=true
RUN_CHILD=true
RUN_MATERNAL=true
RUN_POOLED_STUDY_ANALYSIS=true

cd "${CODING_DIR}"

if [[ "${1:-}" == "combine" ]]; then
  Rscript "${R_SCRIPT}" \
    mode=combine_results \
    project_root="${PROJECT_ROOT}" \
    targeted_rdata="${INPUT_DIR}/Harmonized_Targeted.RData" \
    sdoh_rdata="${INPUT_DIR}/Harmonized_SDOH.RData" \
    out_dir="${OUT_DIR}" \
    nsim="${NSIM}" \
    alpha="${ALPHA}" \
    min_n="${MIN_N}" \
    seed="${SEED}" \
    taus="${TAUS}" \
    use_log1p="${USE_LOG1P}" \
    run_child="${RUN_CHILD}" \
    run_maternal="${RUN_MATERNAL}" \
    run_pooled_study_analysis="${RUN_POOLED_STUDY_ANALYSIS}"
  exit 0
fi

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  Rscript "${R_SCRIPT}" \
    mode=run_task \
    task_id="${SLURM_ARRAY_TASK_ID}" \
    project_root="${PROJECT_ROOT}" \
    targeted_rdata="${INPUT_DIR}/Harmonized_Targeted.RData" \
    sdoh_rdata="${INPUT_DIR}/Harmonized_SDOH.RData" \
    out_dir="${OUT_DIR}" \
    nsim="${NSIM}" \
    alpha="${ALPHA}" \
    min_n="${MIN_N}" \
    seed="${SEED}" \
    taus="${TAUS}" \
    use_log1p="${USE_LOG1P}" \
    run_child="${RUN_CHILD}" \
    run_maternal="${RUN_MATERNAL}" \
    run_pooled_study_analysis="${RUN_POOLED_STUDY_ANALYSIS}"
  exit 0
fi

Rscript "${R_SCRIPT}" \
  mode=build_manifest \
  project_root="${PROJECT_ROOT}" \
  targeted_rdata="${INPUT_DIR}/Harmonized_Targeted.RData" \
  sdoh_rdata="${INPUT_DIR}/Harmonized_SDOH.RData" \
  out_dir="${OUT_DIR}" \
  nsim="${NSIM}" \
  alpha="${ALPHA}" \
  min_n="${MIN_N}" \
  seed="${SEED}" \
  taus="${TAUS}" \
  use_log1p="${USE_LOG1P}" \
  run_child="${RUN_CHILD}" \
  run_maternal="${RUN_MATERNAL}" \
  run_pooled_study_analysis="${RUN_POOLED_STUDY_ANALYSIS}"

TASK_COUNT=$(awk 'NR > 1 {count++} END {print count+0}' "${OUT_DIR}/task_manifest.csv")

ARRAY_JOB_ID=$(sbatch --parsable --array=1-"${TASK_COUNT}" "$0")
COMBINE_JOB_ID=$(sbatch --parsable --dependency=afterok:"${ARRAY_JOB_ID}" "$0" combine)

echo "Manifest built with ${TASK_COUNT} tasks"
echo "Array job id: ${ARRAY_JOB_ID}"
echo "Combine job id: ${COMBINE_JOB_ID}"
