#!/bin/bash
set -euo pipefail

PBS_SCRIPT="/home/565/dh4185/mn51-dh4185/repos_collab/nesp_bff/pbs_scripts/pbs_step3_calc_missing_vars"

LOCS=(
  Melbourne Darwin Cairns Brisbane Longreach Mildura Adelaide Perth Sydney Canberra Hobart
)

for loc in "${LOCS[@]}"; do
  echo "Submitting: $loc"
  qsub -v LOC="$loc" "$PBS_SCRIPT"
done