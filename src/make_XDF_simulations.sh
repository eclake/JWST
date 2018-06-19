#!/bin/bash

dropouts=("B" "V" "I" "Z" "Y")

MC_runs=(1 2 3 4 5 6 7 8 9)
MC_runs=(0)

# Effective radii computed by considering (in redshift bins of delta z=1) the
# mean of the log-normal distribution of FWHM in the UDF master catalogue
# (column UVUDF_FWHM_775_noPsf)
# you can use the command line option
# --effective-radius ${effective_radii[$i]}
# to use the values below
effective_radii=(0.10 0.11 0.12 0.13 0.145)

nbexp=108
filters="CLEAR"
gratings="PRISM"
sersic="1.5" # from Shibuya+2015
nproc=22

root_folder="/home/jchevall/JWST/Simulations/XDF"
folder_name="ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian_weights_mass_SFR_logU"

for MC_run in "${MC_runs[@]}"; do
  i=0
  for drop in "${dropouts[@]}"; do
     folder=${root_folder}/${drop}_DROPOUTS/${folder_name}/MC_${MC_run}/ETC-simulations
     if [ -d ${folder}/ETC-output ] ; then
       rm -r ${folder}/ETC-output
     fi

     input_file=${root_folder}/${drop}_DROPOUTS/${folder_name}/MC_${MC_run}/input-SEDs-IRAC/input_SEDs_MC_${MC_run}.fits

    ./make_ETC_simulations.py -i ${input_file} -o ${folder} --exposures ${nbexp} --filters ${filters} --gratings ${gratings} --sersic ${sersic} --effective-radius shibuya+2015 --nproc ${nproc}
    i=$((i + 1))
  done
done

# B-band dropouts
#./make_ETC_simulations.py -i /Users/jchevall/JWST/Simulations/XDF/B_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian/${MC_run}/input-SEDs/input_SEDs_${MC_run}.fits -o /Users/jchevall/JWST/Simulations/XDF/B_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian/${MC_run}/ETC-simulations --exposures 108 --filters CLEAR --gratings PRISM --sersic 1 --effective-radius 0.10

# V-band dropouts: ratio of angular diameter distance at z=4 and z=5: 1.106 (from http://www.astro.ucla.edu/~wright/CosmoCalc.html, H0=69.6, omega_m=0.286, omega_l=0.714)
#./make_ETC_simulations.py -i /Users/jchevall/JWST/Simulations/XDF/V_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian/${MC_run}/input-SEDs/input_SEDs_${MC_run}.fits -o /Users/jchevall/JWST/Simulations/XDF/V_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian/${MC_run}/ETC-simulations --exposures 108 --filters CLEAR --gratings PRISM --sersic 1 --effective-radius 0.11

# I-band dropouts: ratio of angular diameter distance at z=4 and z=6: 1.215
#./make_ETC_simulations.py -i /Users/jchevall/JWST/Simulations/XDF/I_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian/${MC_run}/input-SEDs/input_SEDs_${MC_run}.fits -o /Users/jchevall/JWST/Simulations/XDF/I_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian/${MC_run}/ETC-simulations --exposures 108 --filters CLEAR --gratings PRISM --sersic 1 --effective-radius 0.12

# Z-band dropouts: ratio of angular diameter distance at z=4 and z=7: 1.328
#./make_ETC_simulations.py -i /Users/jchevall/JWST/Simulations/XDF/Z_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian/${MC_run}/input-SEDs/input_SEDs_${MC_run}.fits -o /Users/jchevall/JWST/Simulations/XDF/Z_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian/${MC_run}/ETC-simulations --exposures 108 --filters CLEAR --gratings PRISM --sersic 1 --effective-radius 0.13

# Y-band dropouts: ratio of angular diameter distance at z=4 and z=8: 1.441
#./make_ETC_simulations.py -i /Users/jchevall/JWST/Simulations/XDF/Y_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian/${MC_run}/input-SEDs/input_SEDs_${MC_run}.fits -o /Users/jchevall/JWST/Simulations/XDF/Y_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian/${MC_run}/ETC-simulations --exposures 108 --filters CLEAR --gratings PRISM --sersic 1 --effective-radius 0.145

