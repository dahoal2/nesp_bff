#!/bin/bash
set -euo pipefail

# -----------------------------
# Fixed directories (your paths)
# -----------------------------
PBS_DIR="/home/565/dh4185/mn51-dh4185/repos_collab/nesp_bff/pbs_scripts"
BASE_DIR="/home/565/dh4185/mn51-dh4185/repos_collab/nesp_bff"

PBS_SCRIPT="${PBS_DIR}/pbs_step1_parallel"
LOGDIR="${PBS_DIR}/job_logs/step1"
mkdir -p "$LOGDIR"

# -----------------------------
# Config
# -----------------------------
DRYRUN="${DRYRUN:-false}"
SLEEP_BETWEEN=1

SCENARIOS=(historical ssp126 ssp370)
#GCMS=(
#  "ACCESS-CM2"
#)
 GCMS=(
   #"ACCESS-CM2"
   "ACCESS-ESM1-5"
   "CESM2"
   "CMCC-ESM2"
   "EC-Earth3"
   "MPI-ESM1-2-HR"
   "NorESM2-MM"
 )

# Time windows
HIST_PERIODS=("1994 2014")
FUTURE_PERIODS=("2020 2040" "2040 2060" "2060 2080")

FREQ="1hr"
NWORKERS=""   # empty -> PBS script uses PBS_NCPUS

# -----------------------------
# Checks
# -----------------------------
command -v qsub >/dev/null 2>&1 || { echo "ERROR: qsub not found in PATH"; exit 1; }
[[ -f "$PBS_SCRIPT" ]] || { echo "ERROR: PBS script not found: $PBS_SCRIPT"; exit 1; }

safe_id() { echo "$1" | sed 's/[^[:alnum:]]/_/g'; }

for SCENARIO in "${SCENARIOS[@]}"; do
  if [[ "$SCENARIO" == "historical" ]]; then
    PERIODS=("${HIST_PERIODS[@]}")
  else
    PERIODS=("${FUTURE_PERIODS[@]}")
  fi

  for GCM in "${GCMS[@]}"; do
    for P in "${PERIODS[@]}"; do
      read -r STARTYEAR ENDYEAR <<< "$P"

      S_SCEN=$(safe_id "${SCENARIO:0:6}")
      S_GCM=$(safe_id "${GCM:0:18}")
      JOBNAME="step1_${S_SCEN}_${S_GCM}_${STARTYEAR}_${ENDYEAR}"

      STDOUT="${LOGDIR}/${JOBNAME}.o"
      STDERR="${LOGDIR}/${JOBNAME}.e"

      ENVVARS="SCENARIO=${SCENARIO},GCM=${GCM},STARTYEAR=${STARTYEAR},ENDYEAR=${ENDYEAR},FREQ=${FREQ}"
      if [[ -n "${NWORKERS}" ]]; then
        ENVVARS="${ENVVARS},NWORKERS=${NWORKERS}"
      fi

      echo "Submitting: $JOBNAME"
      echo "  PBS: $PBS_SCRIPT"
      echo "  Logs: $STDOUT / $STDERR"

      if [[ "$DRYRUN" == true ]]; then
        echo "DRYRUN: qsub -N \"$JOBNAME\" -o \"$STDOUT\" -e \"$STDERR\" -v \"$ENVVARS\" \"$PBS_SCRIPT\""
      else
        qsub -N "$JOBNAME" -o "$STDOUT" -e "$STDERR" -v "$ENVVARS" "$PBS_SCRIPT"
      fi

      sleep "$SLEEP_BETWEEN"
    done
  done
done
