#!/bin/bash

echo "Beginning of script"
date

effects=('corr' 'direction_positive')
exposures=('flood' 'precipitation' 'temperature' 'drought')
outclass=('occurrence' 'incidence')

for effect in "${effects[@]}"
do
  for exposure in "${exposures[@]}"
  do
    for out in "${outclass[@]}"
    do
      echo "Running $variable $exposure $out";
      Rscript Analysis/02_prepare_meta_analysis.R -e $effect  -c $exposure -o $out;
      Rscript Analysis/03_run_meta_analysis.R -e $effect  -c $exposure -o $out;
    done
  done
done

echo "End of script"

date

