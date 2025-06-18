#!/bin/bash

echo "Beginning of script"
date

effects=('direction_positive')
exposures=('flood' 'precipitation' 'temperature' 'drought')
outclass=('occurrence' 'incidence')
levels=('all' 'Asia' 'Africa' 'Other' 'arid' 'temperate' 'tropical' '<50%' '50-95%' '>95%')
strat=('none' 'continent' 'continent' 'continent' 'new_climatic_zone' 'new_climatic_zone' 'new_climatic_zone' 'pct_years_cholera_cat' 'pct_years_cholera_cat' 'pct_years_cholera_cat')

for effect in "${effects[@]}"
do
  for exposure in "${exposures[@]}"
  do
    for out in "${outclass[@]}"
    do
      for i in {0..9}
      do
      echo "------- Running $variable $exposure $out ${strat[i]} ${levels[i]} -------";
      Rscript Analysis/02_prepare_meta_analysis.R -e $effect  -c $exposure -o $out -l ${levels[i]} -s ${strat[i]} -r TRUE;
      Rscript Analysis/03_run_meta_analysis.R -e $effect  -c $exposure -o $out -l ${levels[i]} -s ${strat[i]} -r TRUE;
      done
    done
  done
done

echo "End of script"

date

