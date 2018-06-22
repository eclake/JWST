#!/bin/bash

declare -a dropouts=("B" "V" "I" "Z" "Y")
MC_run="MC_0"

folder="ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian_weights_mass_SFR_logU"
config_json="/Users/jchevall/JWST/Simulations/XDF/config.json"
script="compute_emission_line_SN.py"

for drop in "${dropouts[@]}"; do
    _folder="/Users/jchevall/JWST/Simulations/XDF/${drop}_DROPOUTS/${folder}/${MC_run}/ETC-simulations"
    ./${script}  --folder ${_folder} --json-file /Users/jchevall/JWST/src/emission_lines_SN_config_PRISM.json
done
