#!/usr/bin/env python

##### PS ########################################################
# Signal to noise of ~10 expected for 10 exposures and an input
# flux of 100 nJy => creating a flat input spectrum
# Forcing the summation parameters to those used for the
# sensitivity computation (2,4)

# Creation of the input spectrum (command-line sequence)
from JWSTpylib import c_spectrum
import os
import random
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import argparse
import re

c_light = 2.99792e+18 # Ang/s
show_plot = False

# Name of the Python script that creates the simulated NIRSpec observations
Pierre_procedure = "/Users/jchevall/JWST/code/JWSTpytools-0.0.3/sensitivity/p_spectrumMOS1x3.py"
datapath = "/Users/jchevall/JWST/code/JWSTpytools-0.0.3/data"
pce = "PCE-NIRS30-IFU30-FPA106"

# How many exposures?
#nbexps = ("108", "36", "36", "36")

def compute_ETC_simulation(input_file, FWA, GWA, nbexp, output_folder, output_prefix):

    sys_command = "python2.7 " + Pierre_procedure + " " + input_file + " " + datapath + " " + pce \
            + " " + FWA + " " + GWA + " PS " + nbexp + " " + output_folder + " " + output_prefix

    os.system(sys_command)

    if show_plot:

        hdulist0 = fits.open(input_file)
        wl = hdulist0[1].data['WAVELENGTH']*1.E+06
        flux = hdulist0[1].data['VALUE']

        file_name = output_prefix + "_snr_PS_" + FWA + "_" + GWA + ".fits"

        hdulist = fits.open(os.path.join(output_folder, file_name))
        wl_simulated = hdulist[1].data['WAVELENGTH']*1.E+06
        flux_simulated = hdulist[1].data['NRSPEC']

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        ax.plot(wl_simulated,
            flux_simulated)

        ax.plot(wl,
            flux)

        ax.set_xlim([0.5, 5.0])
        ax.set_xlabel("$\lambda / \mu\\textnormal{m}$ (observed-frame)")
        ax.set_ylabel("Flux (Jy)")
        plt.show()

        hdulist.close()
        hdulist0.close()

def write_ETC_input_file(wl, flux, redshift, file_name):

    # Redshift the SED and wl
    flux_obs = flux / (1.+redshift)
    wl_obs = wl * (1.+redshift)

    # Convert F_lambda [erg s^-1 cm^-2 A^-1] ----> F_nu [erg s^-1 cm^-2 Hz^-1]
    flux_obs = (wl_obs)**2/c_light*flux_obs

    # Scale to Jy
    flux_obs = flux_obs * 1.e+23

    if show_plot:
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(wl_obs*1e-04,
            flux_obs)

        ax.set_xlim([0.5, 5.0])
        ax.set_xlabel("$\lambda / \mu\\textnormal{m}$ (observed-frame)")
        ax.set_ylabel("Flux (Jy)")
        plt.show()

    # Wl in meters
    wl_obs *= 1.e-10
    
    ETC_spectrum = c_spectrum.c_spectrum()

    # lambda in meters, flux in Jy
    n_wl_obs = wl_obs.size

    # Add the SED, pixel-by-pixel
    for i in range(wl_obs.size):

        if i == n_wl_obs-1:
            wl_obs_start = wl_obs[i] - 0.5 * (wl_obs[i]-wl_obs[i-1])
            wl_obs_end = wl_obs[i] + 0.5 * (wl_obs[i]-wl_obs[i-1])
        else:
            wl_obs_start = wl_obs[i] - 0.5 * (wl_obs[i+1]-wl_obs[i])
            wl_obs_end = wl_obs[i] + 0.5 * (wl_obs[i+1]-wl_obs[i])

        pixel = c_spectrum.c_pixel(wl_obs[i], wl_obs_start, wl_obs_end, value=flux_obs[i])

        ETC_spectrum.m_addPixel(pixel)

    # Write the FITS file containing the input file for Pierre's routines
    ETC_spectrum.m_writeToSimpleFITS(file_name)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-o', '--output-dir',
        help="Directory that will contain the ETC-like simulations.",
        action="store", 
        type=str, 
        dest="output_dir", 
        required=True
    )

    parser.add_argument(
        '-i', '--input-catalogue',
        help="FITS file containing the input catalogue of galaxy SEDs produeced by the Beagle tool.",
        action="store", 
        type=str, 
        dest="input_catalogue", 
        required=True
    )

    parser.add_argument(
        '--filters', 
        help="NIRSpec filters to use to simulate NIRSpec observations. \
                If multiple filters are passed, all of them will be used.",
        action="store", 
        type=str, 
        dest="FWAs", 
        nargs='+',
        choices=["CLEAR", "F100LP", "F170LP", "F290LP"],
        default=("CLEAR",)
    )

    parser.add_argument(
        '--gratings', 
        help="NIRSpec gratings to use to simulate NIRSpec observations. \
                If multiple gratings are passed, all of them will be used.",
        action="store", 
        type=str, 
        dest="GWAs", 
        nargs='+',
        choices=["PRISM", "G140M", "G235M", "G395M"],
        default=("PRISM",)
    )

    parser.add_argument(
        '--exposures', 
        help="Number of exposures.",
        action="store", 
        type=str, 
        dest="nbexps", 
        nargs='+',
        default=("108",)
    )

    parser.add_argument(
        '-N',
        help="The first N objects will be processed.",
        action="store", 
        type=int
    )

    parser.add_argument(
        '--seed', 
        help="Seed of the random number generator.",
        action="store", 
        type=int,
        dest="seed",
        default=123456
    )

    parser.add_argument(
        '--shuffle', 
        help="Shuffle the rows of the input catalogue.",
        action="store_true", 
        dest="shuffle" 
    )


    parser.add_argument(
        '--no-recompute', 
        help="Recompute existing data or not.",
        action="store_true", 
        dest="no_recompute" 
    )

    parser.add_argument(
        '--show-plot', 
        help="Show plots of input / output SED.",
        action="store_true", 
        dest="show_plot" 
    )

    # Get parsed arguments
    args = parser.parse_args()    

    # 
    if not len(args.FWAs) == len(args.GWAs) == len(args.nbexps):
        raise ValueError("The length of the filters, gratings, and number of exposures must be the same!")

    # Set seed for random number generator (used only when shuffling the input
    # catalogue rows)
    random.seed(args.seed)

    # Flag to recompute or not existing data
    recompute = not args.no_recompute

    # Set the global variable "show_plot"
    if args.show_plot:
        show_plot = args.show_plot

    # Create the output folders is necessary
    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)

    # Check whether you need to create the folder that will contain the input
    # FITS file for the ETC simulator
    ETC_input_dir = os.path.join(args.output_dir, 'ETC-input')
    if not os.path.isdir(ETC_input_dir):
        os.makedirs(ETC_input_dir)

    # Check whether you need to create the folder that will contain the output
    # FITS file preoduced by the ETC simulator
    ETC_output_dir = os.path.join(args.output_dir, 'ETC-output')
    if not os.path.isdir(ETC_output_dir):
        os.makedirs(ETC_output_dir)

    # Check if the FITS catalogue contains an "MC_#" string
    tmp = os.path.basename(args.input_catalogue) 
    suffix = re.search('MC_(\d+)', tmp)
    if suffix is not None:
        suffix = '_MC_' + suffix.group(1)
    else:
        suffix = ''

    # Open catalogue of input SEDs
    hdulist = fits.open(args.input_catalogue)

    # Get the wavelength array (units of Ang)
    wl = hdulist['full sed wl'].data['wl'][0,:]

    # Get the redshifts of the different SEDs
    redshifts = hdulist['galaxy properties'].data['redshift']

    # By default you create simulated NIRSpec observations for all the objects
    # in the catalogue, but the user can choose to just run on the first N
    # objects (mainly for testing purposes!)
    n_objects = len(redshifts)
    rows = range(n_objects)    
    if args.shuffle:
        random.shuffle(rows)

    if args.N:
        rows = rows[:args.N]

    # Cycle across all objects in the input FITS catalogue
    for row in rows:

        # Redshfit of the i-th object (i.e., i-th row in the input FITS catalogue)
        redshift = redshifts[row]

        # Name of the FITS file containing the input SED for the ETC simulator
        ETC_input_file = os.path.join(ETC_input_dir, str(row+1) + suffix + '_input_for_ETC.fits')

        # By default you always recompute the input file for the ETC, but in
        # some occasions ypu may just want to create the input file for some
        # missing objects 
        if not os.path.isfile(ETC_input_file) or recompute:

            # SED (units are those putput from Beagle, i.e. erg s^-1 cm^-2 A^-1)
            sed = hdulist['full sed'].data[row,:]

            # Function that creates the FITS file that will later be used as
            # input for the ETC simulator. Note that this function simply
            # convert the flux to observed frame, and from F_lambda into F_nu
            # (in Jansky)
            write_ETC_input_file(wl, sed, redshift, ETC_input_file)


        # Prefix for the output file produced by the ETC simulator
        ETC_output_prefix = str(row+1) + suffix

        # Cycle across each combination of filter, grating, and number of exposures
        for FWA, GWA, nbexp in zip(args.FWAs, args.GWAs, args.nbexps):

            # Name of the file created by the ETC simulator (need the name to
            # check if the file already exists or not)
            ETC_output_file = ETC_output_prefix + "_snr_PS_" + FWA + "_" + GWA + ".fits"
            ETC_output_file = os.path.join(ETC_output_dir, ETC_output_file)
            if not os.path.isfile(ETC_output_file) or recompute:

                # Run the actual scripts that compute the ETC-like simulated NIRSpec observation
                compute_ETC_simulation(ETC_input_file, FWA, GWA, nbexp, ETC_output_dir, ETC_output_prefix)

                # The SED output from the ETC simulator is in units of Jansky
                # (F_nu), while Beagle works in F_lambda. We therefore add two
                # columns to the ETC output containing a Beagle-friendly
                # format.

                # Open the file containing the ETC-like simulation
                hduETC = fits.open(ETC_output_file)

                # Get existing columns
                existing_cols = hduETC[1].columns

                # Add new columns
                new_col = list()

                # Add "FLUX_FLAMBDA" column, expressing the flux in F_lambda, erg s^-1 cm^-2 A^-1
                # The units of the NRSPEC column are Jy
                flux = hduETC[1].data['NRSPEC'] * 1.E-23 * c_light / (hduETC[1].data['WAVELENGTH']*1.E+10)**2
                new_col.append(fits.Column(name='FLUX_FLAMBDA', array=flux, format='E', unit='erg s^-1 cm^-2 A^-1'))

                # Add "NOISE_FLAMBDA" column, expressing the flux in F_lambda, erg s^-1 cm^-2 A^-1
                # The units of the NOISE column are Jy
                noise = hduETC[1].data['NOISE'] * 1.E-23 * c_light / (hduETC[1].data['WAVELENGTH']*1.E+10)**2
                new_col.append(fits.Column(name='NOISE_FLAMBDA', array=noise, format='E', unit='erg s^-1 cm^-2 A^-1'))

                new_col_defs = fits.ColDefs(new_col)

                hduETC[1] = fits.BinTableHDU.from_columns(existing_cols+new_col_defs)

                # Add redshift keyword
                hduETC[1].header['redshift'] = float(redshift)

                # Overwrite the FITS file
                hduETC.writeto(ETC_output_file, clobber=True)

                hduETC.close()


    hdulist.close()
