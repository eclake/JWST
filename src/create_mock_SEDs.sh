#!/bin/bash

declare -a dropouts=("B" "V" "I" "Z" "Y")
N=10
suffix="_multiple_draws"
beagle_exec="/Users/jchevall/Coding/BEAGLE-merge/build/BEAGLE"

for drop in "${dropouts[@]}"; do
  echo "drop: $drop"
  files=($(ls /Users/jchevall/JWST/Simulations/XDF/${drop}_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian/Summary_MC_*_min_z_?.?${suffix}.param))
  file=${files[0]}
  for i in $(seq 1 $(($N-1))); do 
    echo "i: $i"
    #echo gsed -i "'"s/MC_$(($i-1))/MC_$i/g"'" $file
    gsed -i s/MC_$(($i-1))/MC_$i/g $file
    ${beagle_exec} --mock -p $file
  done
  #echo gsed -i "'"s/MC_$i/MC_0/g"'" $file
  gsed -i s/MC_$i/MC_0/g $file
done
