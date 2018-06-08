#! /usr/bin/env python

from astropy.io import fits
import os, glob, sys
import matplotlib.pyplot as plt
import numpy as np
import argparse

import autoscale as autoscale

c_light = 2.99792e+18 # Ang/s

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--ID',
        help="ID of the spectrum to plot",
        type=str, 
        dest="ID", 
        required=True
    )

    parser.add_argument(
        '--MC-draw',
        help="Monte Carlo draw",
        dest="MC_draw", 
        type=str,
        default='0'
    )

    parser.add_argument(
        '--range',
        help="Wavelength range to be plotted",
        dest="xrange", 
        type=np.float32,
        nargs=2,
        default=None
    )

    parser.add_argument(
        '--redshift',
        dest="redshift", 
        help="Redshift",
        type=np.float32,
        default=None
        )


    args = parser.parse_args()    

    band, ID = args.ID.split('_')

    folder = "/Users/jchevall/JWST/Simulations/XDF/" + band + "_DROPOUTS/" \
            "ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian/MC_" + args.MC_draw

    # Read magnitude and redshift of the object
    "_WFC3_F160W_APP"

    file_name = os.path.join(folder, "input-SEDs-IRAC", 
            "input_SEDs_MC_" + args.MC_draw + ".fits")

    row = int(ID)-1
    hdulist = fits.open(file_name)
    redshift = hdulist['GALAXY PROPERTIES'].data['redshift'][row]
    m_160 = hdulist['APPARENT MAGNITUDES'].data['_WFC3_F160W_APP'][row]
    SFR = hdulist['STAR FORMATION'].data['SFR_10'][row]
    sSFR = hdulist['STAR FORMATION'].data['sSFR'][row]

    hdulist.close()

    # Load the original (noiseless) SED
    #file_name = os.path.join(folder, "ETC-simulations", "ETC-input", 
    #        ID + "_MC_" + args.MC_draw + "_input_for_ETC.fits")

    #hdulist = fits.open(file_name)
    #wl = hdulist[1].data['WAVELENGTH'] * 1.E+06
    #flux = hdulist[1].data['VALUE']

    #hdulist.close()

    # Load the simulated spectrum 
    file_name = os.path.join(folder, "ETC-simulations", "ETC-output", 
            ID + "_MC_" + args.MC_draw + "_snr_PS_CLEAR_PRISM.fits")

    hdulist = fits.open(file_name)
    wl_simul = hdulist[1].data['WAVELENGTH'] * 1.E+06
    flux_simul = hdulist[1].data['NRSPEC']

    wl = hdulist[1].data['WAVELENGTH'] * 1.E+06
    flux = hdulist[1].data['RSPEC']

    hdulist.close()

    # Plot the spectrum
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    if args.xrange is not None:
        ax.set_xlim(args.xrange)
    else:
        ax.set_xlim((np.amin(wl_simul), np.amax(wl_simul)))

    ax.set_xlabel('$\lambda / \mu\\textnormal{m}$')
    ax.set_ylabel('$f_\\nu / \\textnormal{erg} \, \\textnormal{s}^{-1} \, \\textnormal{cm}^{-2} \, \\textnormal{Hz}^{-1} $')

    alpha = 0.7
    ax.plot(wl,
       flux,
       lw=2.5,
       color="black",
       alpha=alpha
       )

    ax.plot(wl_simul,
       flux_simul,
       lw=1.5,
       color="red",
       alpha=alpha*0.8
       )

    autoscale.autoscale_y(ax)

    y0, y1 = ax.get_ylim()
    ax.set_ylim([0.,y1])

    y = 0.92
    ax.text(0.05, y, 
            "$z=" + "{:.3f}$".format(redshift), 
            transform=ax.transAxes,
            ha='left',
            va='center')

    y -= 0.07
    ax.text(0.05, y, 
            "$H_{160} =" + "{:.2f}$".format(m_160), 
            transform=ax.transAxes,
            ha='left',
            va='center')

    y -= 0.07
    ax.text(0.05, y, 
            "$\\textnormal{SFR} =" + "{:.2f}".format(SFR) + \
            "\; \\textnormal{M}_\odot \, \\textnormal{yr}^{-1}$", 
            transform=ax.transAxes,
            ha='left',
            va='center')

    y -= 0.07
    ax.text(0.05, y, 
            "$\\log(\\textnormal{sSFR}) =" + "{:.2f}".format(sSFR) + \
            "\; \\textnormal{yr}^{-1}$", 
            transform=ax.transAxes,
            ha='left',
            va='center')



    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    plt.title(band + "-dropout, " + ID + "\_MC\_" + args.MC_draw)

    name = os.path.join(folder, 
            ID + "_MC_" + args.MC_draw + "_noiseless_snr_PS_CLEAR_PRISM.pdf")

    fig.savefig(name, dpi=None, facecolor='w', edgecolor='w',
            orientation='portrait', papertype='a4', format="pdf",
            transparent=False, bbox_inches="tight", pad_inches=0.1)

    plt.close(fig)

