#!/bin/bash

# Minimum redshift

#python extract_SEDs.py -r /user_data/jchevall/Emma/BEAGLE_results/WG3_z_lt_4_mock_test_5 -p WG3_z_lt_4_mock_3D_HST_GOODSN_test5.param  -n_samples 1 -n_objects 20 --folder-suffix _new 

python extract_SEDs.py -r /user_data/jchevall/Emma/BEAGLE_results/WG3_z_lt_4_mock_test_5 -p WG3_z_lt_4_mock_3D_HST_GOODSN_test5.param  -n_samples 1 -np 22 --folder-suffix _new 

#python extract_SEDs.py -r /local/jchevall/BEAGLE_results/UVUDF/Emma_IDs/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian -p UVUDF_BEAGLE_Emma.param  -n_samples 4 -n_objects 10 --plot --folder-suffix _new --suffix _min_z_1.0 --params-ranges '{"redshift": {"colName":"redshift", "extName":"galaxy properties", "min":1.0, "step":0.1} }'
#python extract_SEDs.py -r /local/jchevall/BEAGLE_results/UVUDF/Emma_IDs/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian -p UVUDF_BEAGLE_Emma.param  -n_samples 4 -np 22 --folder-suffix _new --suffix _min_z_1.0 --params-ranges '{"redshift": {"colName":"redshift", "extName":"galaxy properties", "min":1.0, "step":0.1} }'

#python extract_SEDs.py -r /local/jchevall/BEAGLE_results/XDF/V_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian -p XDF_BEAGLE_delayed_SFR_V_dropout.param  -n_samples 4 -n_objects 10 --plot --folder-suffix _new --suffix _min_z_3.0 --params-ranges '{"redshift": {"colName":"redshift", "extName":"galaxy properties", "min":3.0, "step":0.1} }'
#python extract_SEDs.py -r /local/jchevall/BEAGLE_results/XDF/V_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian -p XDF_BEAGLE_delayed_SFR_V_dropout.param  -n_samples 4 -np 22 --folder-suffix _new --suffix _min_z_3.0 --params-ranges '{"redshift": {"colName":"redshift", "extName":"galaxy properties", "min":3.0, "step":0.1} }'

#python extract_SEDs.py -r /local/jchevall/BEAGLE_results/UVUDF/Emma_IDs/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian -p UVUDF_BEAGLE_Emma.param  -n_samples 4 -np 22 --suffix _min_z_1.0 --params-ranges '{"redshift": {"colName":"redshift", "extName":"galaxy properties", "min":1.0, "step":0.1} }'

#-#python extract_SEDs.py -r /local/jchevall/BEAGLE_results/XDF/B_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian -n_samples 4 -np 22 --suffix _min_z_2.0 --params-ranges '{"redshift": {"colName":"redshift", "extName":"galaxy properties", "min":2.0, "step":0.1} }'

#-#python extract_SEDs.py -r /local/jchevall/BEAGLE_results/XDF/V_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian -n_samples 4 -np 22 --suffix _min_z_3.0 --params-ranges '{"redshift": {"colName":"redshift", "extName":"galaxy properties", "min":3.0, "step":0.1} }'

#-#python extract_SEDs.py -r /local/jchevall/BEAGLE_results/XDF/I_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian -n_samples 4 -np 22 --suffix _min_z_4.0 --params-ranges '{"redshift": {"colName":"redshift", "extName":"galaxy properties", "min":4.0, "step":0.1} }'

#-#python extract_SEDs.py -r /local/jchevall/BEAGLE_results/XDF/Z_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian -n_samples 4 -np 22 --suffix _min_z_5.0 --params-ranges '{"redshift": {"colName":"redshift", "extName":"galaxy properties", "min":5.0, "step":0.1} }'

#-#python extract_SEDs.py -r /local/jchevall/BEAGLE_results/XDF/Y_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian -n_samples 4 -np 22 --suffix _min_z_5.0 --params-ranges '{"redshift": {"colName":"redshift", "extName":"galaxy properties", "min":5.0, "step":0.1} }'

# Minimum redshift and mass

#python extract_SEDs.py -r /local/jchevall/BEAGLE_results/UVUDF/Emma_IDs/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian -n_samples 4 -np 22 --suffix _min_z_1.0_min_M_4.0e6 --params-ranges '{"redshift": {"colName":"redshift", "extName":"galaxy properties", "min":1.0, "step":0.1}, "mass": {"colName":"M_tot", "extName":"galaxy properties", "min":4.0E+06}}'

#python extract_SEDs.py -r /local/jchevall/BEAGLE_results/XDF/B_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian -n_samples 4 -np 22 --suffix _min_z_2.0_min_M_1.0e7 --params-ranges '{"redshift": {"colName":"redshift", "extName":"galaxy properties", "min":2.0, "step":0.1}, "mass": {"colName":"M_tot", "extName":"galaxy properties", "min":1.0E+07}}'

#python extract_SEDs.py -r /local/jchevall/BEAGLE_results/XDF/V_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian -n_samples 4 -np 22 --suffix _min_z_3.0_min_M_1.5e7  --params-ranges '{"redshift": {"colName":"redshift", "extName":"galaxy properties", "min":3.0, "step":0.1}, "mass": {"colName":"M_tot", "extName":"galaxy properties", "min":1.5E+07}}'

#python extract_SEDs.py -r /local/jchevall/BEAGLE_results/XDF/I_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian -n_samples 4 -np 22 --suffix _min_z_4.0_min_M_2.5e7  --params-ranges '{"redshift": {"colName":"redshift", "extName":"galaxy properties", "min":4.0, "step":0.1}, "mass": {"colName":"M_tot", "extName":"galaxy properties", "min":2.5E+07}}'

#python extract_SEDs.py -r /local/jchevall/BEAGLE_results/XDF/Z_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian -n_samples 4 -np 22 --suffix _min_z_5.0_min_M_2.5e7  --params-ranges '{"redshift": {"colName":"redshift", "extName":"galaxy properties", "min":5.0, "step":0.1}, "mass": {"colName":"M_tot", "extName":"galaxy properties", "min":2.5E+07}}'

#python extract_SEDs.py -r /local/jchevall/BEAGLE_results/XDF/Y_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian -n_samples 4 -np 22 --suffix _min_z_5.0_min_M_2.5e7  --params-ranges '{"redshift": {"colName":"redshift", "extName":"galaxy properties", "min":5.0, "step":0.1}, "mass": {"colName":"M_tot", "extName":"galaxy properties", "min":2.5E+07}}'
