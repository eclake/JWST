#!/bin/bash

declare -a dropouts=("B" "V" "I" "Z" "Y")

dest_folder="/Users/jchevall/JWST/Simulations/XDF/ALL_DROPOUTS/MC_0/ETC-simulations/ETC-output"

for drop in "${dropouts[@]}"; do
   folder=/Users/jchevall/JWST/Simulations/XDF/${drop}_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian/MC_0/ETC-simulations/ETC-output
   cd $folder
   for file in *.fits; do 
     #echo "cp ${file} ${dest_folder}/${drop}_${file}"
     cp -p ${file} ${dest_folder}/${drop}_${file}
   done
done
