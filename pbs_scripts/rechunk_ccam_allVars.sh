#!/bin/bash

VARS=('tas' 'hurs' 'huss' 'sfcWind' 'psl' 'uas' 'vas' 'clt' 'rsds' 'rsdsdir')
GCM="ACCESS-CM2"
SCENARIO="historical"
START_YEAR=1985
END_YEAR=2014

for VAR in "${VARS[@]}"; do
  qsub -v VAR="$VAR",GCM="$GCM",SCENARIO="$SCENARIO",START_YEAR="$START_YEAR",END_YEAR="$END_YEAR" rechunk_ccam.pbs
done

