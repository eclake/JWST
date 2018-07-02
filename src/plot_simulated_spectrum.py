#! /usr/bin/env python

from astropy.io import fits
import os, glob, sys, fnmatch
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
        '--folder',
        help="Folder containing the simulated NIRSpec spectra",
        type=str,
        dest="folder", 
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

    parser.add_argument(
        '--save',
        dest="savefig", 
        help="Save the figure as a pdf",
        type=bool,
        default=False
        )


    args = parser.parse_args()    

    dropout_band = args.folder.split("_DROPOUTS")[0][-1]

    file_name = os.path.join(args.folder, "MC_" + args.MC_draw, "input-SEDs-IRAC", 
            "input_SEDs_MC_" + args.MC_draw + ".fits")

    row = int(args.ID)-1
    hdulist = fits.open(file_name)
    redshift = hdulist['GALAXY PROPERTIES'].data['redshift'][row]
    for m in hdulist['APPARENT MAGNITUDES'].columns.names:
        if "F160W" in m:
            m_160 = hdulist['APPARENT MAGNITUDES'].data[m][row]
    SFR = hdulist['STAR FORMATION'].data['SFR_10'][row]
    sSFR = hdulist['STAR FORMATION'].data['sSFR'][row]

    hdulist.close()

    _file = None
    for file in os.listdir(args.folder):
        if fnmatch.fnmatch(file, "Summary_MC_" + args.MC_draw + "*.fits") and "_input_for_mock" not in file:
            _file = file
            break
    hdulist = fits.open(os.path.join(args.folder, _file))
    XDF_ID = hdulist['META DATA'].data['ID'][row]

    hdulist.close()


    # Load the original (noiseless) SED
    #file_name = os.path.join(folder, "ETC-simulations", "ETC-input", 
    #        ID + "_MC_" + args.MC_draw + "_input_for_ETC.fits")

    #hdulist = fits.open(file_name)
    #wl = hdulist[1].data['WAVELENGTH'] * 1.E+06
    #flux = hdulist[1].data['VALUE']

    #hdulist.close()

    # Load the simulated spectrum 
    file_name = os.path.join(args.folder, "MC_" + args.MC_draw, "ETC-simulations", "ETC-output", 
            args.ID + "_MC_" + args.MC_draw + "_snr_PS_CLEAR_PRISM.fits")

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

    plt.title(args.ID + "\_MC\_" + args.MC_draw + " " + XDF_ID)

    plt.show()

    if args.savefig:
        name = os.path.join(args.folder, "MC_" + args.MC_draw, "ETC-simulations", "ETC-output",
                args.ID + "_MC_" + args.MC_draw + "_noiseless_snr_PS_CLEAR_PRISM.pdf")

        fig.savefig(name, dpi=None, facecolor='w', edgecolor='w',
                orientation='portrait', papertype='a4', format="pdf",
                transparent=False, bbox_inches="tight", pad_inches=0.1)

    plt.close(fig)

