#!/bin/bash

declare -a dropouts=("B" "V" "I" "Z" "Y")
N=10
suffix="_multiple_draws_mass_SFR_logU"
beagle_exec="/Users/jchevall/Coding/BEAGLE/build/BEAGLE"

for drop in "${dropouts[@]}"; do
  echo "drop: $drop"
  root=/Users/jchevall/JWST/Simulations/XDF/${drop}_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian_weights_mass_SFR_logU
  files=($(ls ${root}/Summary_MC_?_min_z_?_min_M_?e?${suffix}.param))
  file=${files[0]}
  for i in $(seq 0 $(($N-1))); do 
    echo "i: $i"
    #echo gsed -i "'"s/MC_$(($i-1))/MC_$i/g"'" $file
    if [ ${i} -gt 0 ]; then
      gsed -i s/MC_$(($i-1))/MC_$i/g $file
    fi
    beagle_file="${root}/MC_${i}/input-SEDs-IRAC/input_SEDs_MC_${i}.fits"
    echo $beagle_file
    if [ ! -s "${beagle_file}" ]; then
      ${beagle_exec} --mock -p $file
    fi
  done
  #echo gsed -i "'"s/MC_$i/MC_0/g"'" $file
  gsed -i s/MC_$i/MC_0/g $file
done
