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
from astropy.cosmology import FlatLambdaCDM
import numpy as np
import matplotlib.pyplot as plt
import argparse
import re

import sys
sys.path.append(os.path.join(os.environ['PYP_BEAGLE'], "PyP-BEAGLE"))
import beagle_multiprocess

from pathos.multiprocessing import ProcessingPool 

c_light = 2.99792e+18 # Ang/s
show_plot = False

# Name of the Python script that creates the simulated NIRSpec observations
jwstpytools = os.environ['JWSTPYTOOLS']
jwstpytools_procedure = os.path.join("p_spectrumMOS1x3_JC.py")
jwstpytools_data = os.path.join(jwstpytools, "data")

# We use the newest Photon Counting Efficiency tables, the same that Pierre provided STScI for their ETC
pce = "PCE-OTE07-NIRS40-IFU31-FPA106-ETC"

ETC_output_dir = ""
ETC_input_dir = ""

# How many exposures?
#nbexps = ("108", "36", "36", "36")

def compute_ETC_simulation(input_file, FWA, GWA, nbexp, output_folder, output_prefix,
        sersic=None, effective_radius=None, seed=None):

    sys_command = "python2.7 " + jwstpytools_procedure + " " + input_file + " " + jwstpytools_data + " " + pce \
            + " " + FWA + " " + GWA + " PS " + nbexp + " " + output_folder + " " + output_prefix 

    if sersic is not None and effective_radius is not None:
            sys_command += " --sersic " + str(sersic) + " --effective-radius " + str(effective_radius)

    if seed is not None:
            sys_command += " --seed " + str(seed)

    print "sys_command: ", sys_command

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


def make_ETC_simulations_single(ETC_simulation_prefix, 
        wl, sed, redshift,
        recompute, 
        FWAs, GWAs, nbexps,
        sersic=None, effective_radius=None, seed=None):

    # SED (units are those putput from Beagle, i.e. erg s^-1 cm^-2 A^-1)
    # Redshfit of the i-th object (i.e., i-th row in the input FITS catalogue)

    # Name of the FITS file containing the input SED for the ETC simulator
    ETC_input_file = os.path.join(ETC_input_dir, ETC_simulation_prefix + '_input_for_ETC.fits')

    # By default you always recompute the input file for the ETC, but in
    # some occasions ypu may just want to create the input file for some
    # missing objects 
    if not os.path.isfile(ETC_input_file) or recompute:

        # Function that creates the FITS file that will later be used as
        # input for the ETC simulator. Note that this function simply
        # convert the flux to observed frame, and from F_lambda into F_nu
        # (in Jansky)
        write_ETC_input_file(wl, sed, redshift, ETC_input_file)

    # Cycle across each combination of filter, grating, and number of exposures
    for FWA, GWA, nbexp in zip(args.FWAs, args.GWAs, args.nbexps):

        # Name of the file created by the ETC simulator (need the name to
        # check if the file already exists or not)
        ETC_output_file = ETC_simulation_prefix + "_snr_PS_" + FWA + "_" + GWA + ".fits"
        ETC_output_file = os.path.join(ETC_output_dir, ETC_output_file)
        if not os.path.isfile(ETC_output_file) or recompute:

            # Run the actual scripts that compute the ETC-like simulated NIRSpec observation
            compute_ETC_simulation(ETC_input_file, FWA, GWA, nbexp, ETC_output_dir, ETC_simulation_prefix, sersic, effective_radius, seed)

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
            hduETC.writeto(ETC_output_file, overwrite=True)

            hduETC.close()

def Shibuya_sizes(redshift, L_UV):

    # See Williams et al 2018, Sec 5.2 (equation 28)
    M_UV_0 = -21.
    L_UV_0 = 10.**(-0.4*(M_UV_0-48.6))

    r_eff_0 = 6.9 * (1.+redshift)**(-1.2)
    r_eff = r_eff_0 * (L_UV/L_UV_0)**0.27

    print "Median radius (kpc): ", np.median(r_eff)

    # We use the same cosmology as in Shibuya to convert the sizes back in arcsec
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)
    r_eff *= cosmo.arcsec_per_kpc_proper(redshift)
    print "Median radius (arcsec): ", np.median(r_eff)

    return r_eff

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
        '--nproc',
        help="Number of processors to use",
        action="store", 
        type=int, 
        dest="nproc",
        default=-1
    )

    parser.add_argument(
        '--sersic',
        help="Sersic index of the source",
        action="store", 
        type=float, 
        dest="sersic"
    )

    parser.add_argument(
        '--effective-radius',
        help="Effective radius (in arcsec) of the source",
        action="store", 
        dest="effective_radius"
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

    # If the args.effective_radius is a number, then use the same radius for all galaxies
    r_eff = np.array((None,)*len(redshifts))
    if args.effective_radius is not None:
        try:
            r_eff = (float(args.effective_radius),)*len(redshifts)
        except:
            if args.effective_radius.lower() == 'shibuya+2015':
                L_UV = 10.**(hdulist['galaxy properties'].data['L_UV'])
                r_eff = Shibuya_sizes(redshift=redshifts, L_UV=L_UV)
            else:
                raise ValueError("Optional argument --effective-radius `" + args.effective_radius + 
                        "` not recognized")
    
    # Get the SEDs
    SEDs_ = list()
    for i in range(hdulist['full sed'].data[:,:].shape[0]):
        SEDs_.append(hdulist['full sed'].data[i,:])

    hdulist.close()

    # By default you create simulated NIRSpec observations for all the objects
    # in the catalogue, but the user can choose to just run on the first N
    # objects (mainly for testing purposes!)
    n_objects = len(redshifts)
    rows = np.array(range(n_objects), dtype=int)

    # If requested, shuffle the roder of the rows
    if args.shuffle:
        random.shuffle(rows)

    # The user can choose to compute the simulated spectra only for a subset of
    # the possible input SEDs
    if args.N:
        rows = rows[:args.N]

    # Re-ordert the SEDs_ list to follow the suffled order of the rows
    SEDs = list()
    for row in rows:
        SEDs.append(SEDs_[row])

    # Re-order the redshifts array to follow the suffled order of the rows
    redshifts = redshifts[rows]
    r_eff = r_eff[rows]

    # Create a list containing the prefix used for the input and output file
    # for the ETC simulator
    ETC_simulation_prefixes = list()
    for row in rows:
        ETC_simulation_prefixes.append(str(row+1) + suffix)

    # If the user does not specify the number of processors to be used, assume that it is a serial job
    if args.nproc <= 0:

        for i in range(len(rows)):
             make_ETC_simulations_single(
                ETC_simulation_prefix=ETC_simulation_prefixes[i],
                wl=wl,
                sed=SEDs[i],
                redshift=redshifts[i],
                recompute=recompute,
                FWAs=args.FWAs,
                GWAs=args.GWAs,
                nbexps=args.nbexps,
                sersic=args.sersic,
                effective_radius=r_eff[i],
                seed=args.seed
                )
    
    # Otherwise you use pathos to run in parallel on multiple CPUs
    else:

        # Set number of parellel processes to use
        pool = ProcessingPool(nodes=args.nproc)

        # Launch the actual calculation on multiple processesors
        pool.map(make_ETC_simulations_single, 
            ETC_simulation_prefixes,
            (wl,)*len(rows),
            SEDs,
            redshifts,
            (recompute,)*len(rows),
            (args.FWAs,)*len(rows),
            (args.GWAs,)*len(rows),
            (args.nbexps,)*len(rows),
            (args.sersic,)*len(rows),
            r_eff,
            (args.seed,)*len(rows)
            )
