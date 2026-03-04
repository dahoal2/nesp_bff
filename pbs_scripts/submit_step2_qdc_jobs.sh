#!/bin/bash
set -euo pipefail

PBS="/home/565/dh4185/mn51-dh4185/repos_collab/nesp_bff/pbs_scripts/pbs_qdc_step2"

# Edit if needed — ideally match utils.locations keys
#LOCS=(Darwin)
LOCS=(
  Melbourne Darwin Cairns Brisbane Longreach Mildura Adelaide Perth Sydney Canberra Hobart
)

SCENARIOS=(ssp126 ssp370)

for loc in "${LOCS[@]}"; do
  for ssp in "${SCENARIOS[@]}"; do
    echo "Submitting: $loc $ssp"
    qsub -v LOC="$loc",SCENARIO="$ssp" "$PBS"
  done
done
