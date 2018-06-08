#!/bin/bash

declare -a dropouts=("B" "V" "I" "Z" "Y")

suffix="_multiple_draws"
config_json="/Users/jchevall/JWST/Simulations/XDF/config.json"

for drop in "${dropouts[@]}"; do
  files=$(ls /Users/jchevall/JWST/Simulations/XDF/${drop}_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian/Summary_MC_*_min_z_?.?${suffix}.fits)
  for file in $files ; do
    i=$(awk -v a="$file" -v b=".fits" 'BEGIN{print index(a,b)}')
    output=${file:0:$i-1}"_input_for_mock.fits"
    echo "output: $output"
    ./modify_total_mass_in_mock.py -i $file -o $output --config $config_json 
  done
done

# B-band dropouts
#./modify_total_mass_in_mock.py -i /Users/jchevall/JWST/Simulations/XDF/B_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian/Summary_MC_0_min_z_3.0.fits -o /Users/jchevall/JWST/Simulations/XDF/B_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian/Summary_MC_0_min_z_3.0_input_for_mock.fits --config /Users/jchevall/JWST/Simulations/XDF/config.json

# V-band dropouts
#./modify_total_mass_in_mock.py -i /Users/jchevall/JWST/Simulations/XDF/V_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian/Summary_MC_0_min_z_3.0.fits -o /Users/jchevall/JWST/Simulations/XDF/V_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian/Summary_MC_0_min_z_3.0_input_for_mock.fits --config /Users/jchevall/JWST/Simulations/XDF/config.json

# I-band dropouts
#./modify_total_mass_in_mock.py -i /Users/jchevall/JWST/Simulations/XDF/I_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian/Summary_MC_0_min_z_4.0.fits -o /Users/jchevall/JWST/Simulations/XDF/I_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian/Summary_MC_0_min_z_4.0_input_for_mock.fits --config /Users/jchevall/JWST/Simulations/XDF/config.json

# Z-band dropouts
#./modify_total_mass_in_mock.py -i /Users/jchevall/JWST/Simulations/XDF/Z_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian/Summary_MC_0_min_z_5.0.fits -o /Users/jchevall/JWST/Simulations/XDF/Z_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian/Summary_MC_0_min_z_5.0_input_for_mock.fits --config /Users/jchevall/JWST/Simulations/XDF/config.json

# Y-band dropouts
#./modify_total_mass_in_mock.py -i /Users/jchevall/JWST/Simulations/XDF/Y_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian/Summary_MC_0_min_z_5.0.fits -o /Users/jchevall/JWST/Simulations/XDF/Y_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian/Summary_MC_0_min_z_5.0_input_for_mock.fits --config /Users/jchevall/JWST/Simulations/XDF/config.json

