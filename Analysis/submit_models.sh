#!/bin/bash

echo "Beginning of script"
date

effects=('corr' 'diff' 'ratio')
exposures=('flood' 'precipitation' 'temperature' 'drought')

for effect in "${effects[@]}"
do
  for exposure in "${exposures[@]}"
  do
    echo "Running $variable $exposure";
    Rscript Analysis/02_prepare_meta_analysis.R -e $effect  -c $exposure;
    Rscript Analysis/03_run_meta_analysis.R -e $effect  -c $exposure;
  done
done

echo "End of script"

date

