#!/bin/bash

# Usage: 
#   ./rechunk_ccam.sh GCM scenario start_year end_year
# e.g.
#   ./rechunk_ccam.sh ACCESS-CM2 historical 1985 2014

set -euo pipefail

# Load NCO module (adjust if needed)
module load nco

# Arguments
GCM=${1:?Error: No GCM specified}
SCENARIO=${2:?Error: No scenario specified}
START_YEAR=${3:?Error: No start year specified}
END_YEAR=${4:?Error: No end year specified}
VAR=${5:-}  # Optional, if empty process all variables

# Variables list
ALL_VARS=(tas hurs huss sfcWind psl uas vas clt rsds rsdsdir)

# If a single variable is specified, overwrite the array
if [ -n "$VAR" ]; then
    VARS=("$VAR")
else
    VARS=("${ALL_VARS[@]}")
fi

# Base input directory (CORDEX structure)
BASE_DIR="/g/data/hq89/CCAM/output/CMIP6/DD/AUS-10i/CSIRO"

# Output directory (ensure it exists)
OUTPUT_DIR="/scratch/eg3/dh4185/rechunked/${GCM}/${SCENARIO}/"
mkdir -p "$OUTPUT_DIR"

echo "Starting rechunk for GCM=$GCM, Scenario=$SCENARIO, Years=$START_YEAR-$END_YEAR"
echo "Output dir: $OUTPUT_DIR"

for VAR in "${VARS[@]}"; do
    echo "Processing variable: $VAR"

    # Construct input file path pattern according to CORDEX directory structure:
    # /g/data/hq89/CCAM/output/CMIP6/DD/AUS-10i/CSIRO/$GCM/$SCENARIO/*/*/*/1hr/$VAR/*/*.nc
    files=( "$BASE_DIR/$GCM/$SCENARIO"/*/*/*/1hr/"$VAR"/*/*.nc )

    # Check if files found
    if [ ${#files[@]} -eq 0 ]; then
        echo "No files found for variable $VAR"
        continue
    fi

    for file in "${files[@]}"; do
        filename=$(basename "$file")

        # Extract start and end dates from filename (expects ..._YYYYMMDD-YYYYMMDD.nc)
        file_start=$(echo "$filename" | sed -n 's/.*_\([0-9]\{8\}\)-\([0-9]\{8\}\).nc/\1/p')
        file_end=$(echo "$filename" | sed -n 's/.*_\([0-9]\{8\}\)-\([0-9]\{8\}\).nc/\2/p')

        # Validate extraction
        if [[ -z "$file_start" || -z "$file_end" ]]; then
            echo "Skipping $filename: Could not parse date range"
            continue
        fi

        # Extract years
        start_year_file=${file_start:0:4}
        end_year_file=${file_end:0:4}

        # Check if file overlaps desired time period
        if (( end_year_file >= START_YEAR && start_year_file <= END_YEAR )); then
            echo "Rechunking $filename ..."

            # Compose output file path
            outfile="$OUTPUT_DIR/$filename"

            if [ -f $outfile ]; then
                echo "File $outfile exists. Skipping."
            else
                echo "Creating it..."
                # Rechunk lat, lon to 15, time chunk 2000 with ncks; overwrite if exists
                ncks -O --cnk_dmn lat,15 --cnk_dmn lon,15 --cnk_dmn time,2000 "$file" "$outfile"
            fi
        else
            echo "Skipping $filename (outside time range)"
        fi
    done
done

echo "Rechunking completed."

