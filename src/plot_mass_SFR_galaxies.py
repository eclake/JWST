import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import matplotlib.colors as colors
import matplotlib.patches as patches
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from itertools import izip as zip, count
import scipy.stats as stats

folder = "/Users/jchevall/JWST/Simulations/Sep_2016/mass_SFR_weighted_Gaussian"

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

dropout_bands = ('B', 'V', 'I', 'Z', 'Y')
dropout_redshift = (4, 5, 6, 7, 8)
dropout_colors = ''
c_light = 2.99792e+18 # Ang/s
Ha_threshold = np.log10(6.E-19)
m_UV_lim = 28.

def mass_SFR_relation(z, logM):

    t = cosmo.age(z).value
    logSFR = (0.84-0.026*t)*logM - (6.51-0.11*t)

    return logSFR

def plot_single_mass_SFR(fileName, ax, xPar, yPar, xLog, yLog):

    cm = plt.get_cmap('Spectral') 

    cNorm  = colors.Normalize(vmin=0, vmax=len(dropout_bands))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)

    hdulist = fits.open(fileName)

    IDs = hdulist[1].data['ID']
    redshift = hdulist[1].data['z']

    if xPar.lower() == "uv_1500_flambda":
        tmp = hdulist[1].data[xPar] / (1.+redshift)
        wl_central = 1500. * (1+redshift)
        xData = -2.5*np.log10(wl_central**2/c_light*tmp) - 48.6
    elif xPar.lower() == "f160w":
        ff = "/Users/jchevall/Coding/BEAGLE/files/data/XDF/XDF_DROPOUTS_aper_corr.fits"
        hdutmp = fits.open(ff)
        tmpID = hdutmp[1].data['ID']
        tmpF160 = hdutmp[1].data['mag_F160W']
        F160 = np.zeros(len(IDs))
        for i, ID in enumerate(IDs):
            indx = [ii for ii, j in zip(count(), tmpID) if j == ID]
            j = indx[0]
            F160[i] = tmpF160[j]
        hdutmp.close()
        xData = F160
    else:
        xData = hdulist[1].data[xPar]

    if xLog:
        xData = np.log10(xData)

    yData = hdulist[1].data[yPar]

    if yLog:
        yData = np.log10(yData)

    MM = np.array(ax.get_xlim())

    for i, dropout in enumerate(dropout_bands):

        colorVal = scalarMap.to_rgba(len(dropout_bands)-i)
        indices = [j for j, s in enumerate(IDs) if 'XDF'+dropout+'-' in s]

        if yPar == "Halpha_flux":
            n = len(np.where((xData[indices] < m_UV_lim) & (yData[indices] > Ha_threshold))[0])
            n_tot = len(np.where((xData[indices] < m_UV_lim))[0])

            ax.text(0.08, 0.95-i*0.05, "{:.2f}".format(1.*n/n_tot),
                    horizontalalignment='left',
                    verticalalignment='center',
                    fontsize=14, color=colorVal, 
                    weight="heavy",
                    transform=ax.transAxes)


        if xPar == "M_star" and yPar == "SFR":
            lgS = mass_SFR_relation(dropout_redshift[i], MM)
            ax.plot(MM,
                    lgS,
                    color=colorVal,
                    lw=2,
                    zorder=2)

        ax.plot(xData[indices],
                yData[indices],
                marker='o',
                ls="",
                color=colorVal,
                zorder=1
                )


        ax.text(0.05, 0.95-i*0.05, dropout,
                horizontalalignment='center',
                verticalalignment='center',
                fontsize=18, color=colorVal, 
                weight="heavy",
                transform=ax.transAxes)

    if yPar == "Halpha_flux":
        xx = np.array(ax.get_xlim())
        yy = np.array(ax.get_ylim())
        ax.plot(xx,
                xx*0.+Ha_threshold,
                color="black",
                lw=1)

        pp = patches.Rectangle(
                (xx[0], Ha_threshold), 28.-xx[0], yy[1]-yy[0],
                alpha=0.3,
                color="grey",
                fill=True,
                transform=ax.transData)

        ax.add_patch(pp)

    hdulist.close()

def plot_mass_SFR(draw='all', xlim=None, ylim=None, 
        xPar="M_star", yPar="SFR", 
        xLog=True, yLog=True,
        xlabel=None, ylabel=None):

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    if xlabel is not None:
        ax.set_xlabel(xlabel)
    elif xPar == "M_star":
        ax.set_xlabel('$\log(\\textnormal{M} / \\textnormal{M}_\odot)$')
    else:
        ax.set_xlabel(xPar.replace("_", " "))

    if ylabel is not None:
        ax.set_ylabel(ylabel)
    elif yPar == "SFR":
        ax.set_ylabel('$\log(\psi/\\textnormal{M}_\odot \, \\textnormal{yr}^{-1})$')
    else:
        ax.set_ylabel(yPar.replace("_", " "))

    if xlim is not None:
        ax.set_xlim(xlim)

    if ylim is not None:
        ax.set_ylim(ylim)

    file = os.path.join(folder, "Summary_MC_"+str(draw)+".fits")

    plot_single_mass_SFR(file, ax, xPar=xPar, yPar=yPar, xLog=xLog, yLog=yLog)

    fig.suptitle('MC-'+str(draw), fontsize=20)

    fileName = xPar+"_"+yPar+"_MC_"+str(draw)+".pdf"
    fileName = os.path.join(folder, fileName)

    fig.savefig(fileName, dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', papertype='a4', format="pdf",
        transparent=False, bbox_inches="tight", pad_inches=0.1)

    plt.close(fig)


def plot_hist(draw=None, xlim=None, nbins=20, 
        xPar=None, 
        xLog=True,
        xlabel=None):

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    cm = plt.get_cmap('Spectral') 

    cNorm  = colors.Normalize(vmin=0, vmax=len(dropout_bands))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)

    if xlabel is not None:
        ax.set_xlabel(xlabel)
    else:
        ax.set_xlabel(xPar.replace("_", " "))

    ax.set_ylabel("N objects")

    fileName = os.path.join(folder, "Summary_MC_"+str(draw)+".fits")

    hdulist = fits.open(fileName)

    IDs = hdulist[1].data['ID']
    data = hdulist[1].data[xPar]

    if xLog:
        if "EW" in xPar:
            data = -data
        data = np.log10(data)

    if xlim is None:
        xlim = (np.amin(data), np.amax(data))

    for i, dropout in enumerate(dropout_bands):

        colorVal = scalarMap.to_rgba(len(dropout_bands)-i)
        indices = [j for j, s in enumerate(IDs) if 'XDF'+dropout+'-' in s]

        ax.hist(data[indices], nbins, range=xlim, 
                color=colorVal, alpha=0.3, lw=0)

        ax.text(0.05, 0.95-i*0.05, dropout,
                horizontalalignment='center',
                verticalalignment='center',
                fontsize=18, color=colorVal, 
                weight="heavy",
                transform=ax.transAxes)

        if "Halpha" in xPar:
            yy = ax.get_ylim()
            h = sorted(data)
            pdf = stats.norm.pdf(h, 2.57, 0.25) * 0.5*yy[1]
            ax.plot(h, pdf, color="black") 

    hdulist.close()

    fig.suptitle('MC-'+str(draw), fontsize=20)

    fileName = xPar+"_histogram_MC_"+str(draw)+".pdf"
    fileName = os.path.join(folder, fileName)

    fig.savefig(fileName, dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', papertype='a4', format="pdf",
        transparent=False, bbox_inches="tight", pad_inches=0.1)

    plt.close(fig)



if __name__ == "__main__":

    # ***************************************************
    # Plotting the mass-SFR relation

    for draw in range(10):
        plot_mass_SFR(draw=draw, xlim=[6.5, 11], ylim=[-5, 4])

        #plot_mass_SFR(draw=draw, 
        #        xPar="f160w", yPar="Halpha_flux", 
        #        xlim=[24, 32], ylim=[-20,-16],
        #        xlabel="F160W (AB mag)", ylabel="$\\textnormal{H}\\alpha \; \\textnormal{flux}$",
        #        xLog=False)

        #plot_mass_SFR(draw=draw, 
        #        xPar="z", yPar="Halpha_flux", 
        #        xlim=[3,9], ylim=[-20,-16],
        #        xlabel="redshift", ylabel="$\\textnormal{H}\\alpha \; \\textnormal{flux}$",
        #        xLog=False)

        #plot_hist(draw, xPar="Halpha_EW", nbins=15,
        #        xlabel="$\log(\\textnormal{H}\\alpha \; \\textnormal{EW} / \\textnormal{\AA})$",
        #        xLog=True)

        #plot_mass_SFR(draw=draw, 
        #        xPar="uv_1500_flambda", yPar="Halpha_flux", 
        #        xlim=[24, 31.5], ylim=[-20,-16],
        #        xlabel="$m_\\textnormal{AB} \; 1500 \; \\textnormal{\AA}$", 
        #        ylabel="$\log(\\textnormal{H}\\alpha / \\textnormal{erg}\;\\textnormal{cm}^2 \; \\textnormal{s}^{-1})$",
         #       xLog=False)



