#!/bin/bash

echo "Beginning of script"
date

effects=('corr' 'direction_positive')
exposures=('flood' 'precipitation' 'temperature' 'drought')
outclass=('occurrence' 'incidence')
levels=('all' 'Asia' 'Africa' 'other' 'Arid' 'Temperate' 'Tropical')
strat=('none' 'continent' 'continent' 'continent' 'new_climatic_zone' 'new_climatic_zone' 'new_climatic_zone')

for effect in "${effects[@]}"
do
  for exposure in "${exposures[@]}"
  do
    for out in "${outclass[@]}"
    do
      for i in {0..6}
      do
      echo "Running $variable $exposure $out";
      Rscript Analysis/02_prepare_meta_analysis.R -e $effect  -c $exposure -o $out -l ${levels[i]} -s ${strat[i]};
      Rscript Analysis/03_run_meta_analysis.R -e $effect  -c $exposure -o $out -l ${levels[i]} -s ${strat[i]};
    done
    done
  done
done

echo "End of script"

date

