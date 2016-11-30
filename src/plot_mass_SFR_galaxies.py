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
from walker_random_sampling import WalkerRandomSampling
import argparse

folder = ""

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

def sample_mass_SFR(logM, SFR, z):

    mean = mass_SFR_relation(z, logM)

    pdf = np.zeros(len(logM))
    for i, value in enumerate(SFR):
        pdf[i] = stats.norm.pdf(SFR[i], loc=mean[i], scale=0.3)

    return pdf


def plot_X_vs_Y(draw, xPar, yPar,
        xlim=None, ylim=None, 
        xLog=True, yLog=True,
        xlabel=None, ylabel=None):

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    if xlabel is not None:
        ax.set_xlabel(xlabel)
    elif "M_star" in xPar["colName"]:
        ax.set_xlabel('$\log(\\textnormal{M}_\\ast / \\textnormal{M}_\odot)$')
    elif "M_tot" in xPar["colName"]:
        ax.set_xlabel('$\log(\\textnormal{M}^\\textnormal{tot}_\\ast / \\textnormal{M}_\odot)$')
    else:
        ax.set_xlabel(xPar["colName"].replace("_", " "))

    if ylabel is not None:
        ax.set_ylabel(ylabel)
    elif "SFR" in yPar["colName"]:
        ax.set_ylabel('$\log(\psi/\\textnormal{M}_\odot \, \\textnormal{yr}^{-1})$')
    else:
        ax.set_ylabel(yPar["colName"].replace("_", " "))

    if xlim is not None:
        ax.set_xlim(xlim)

    if ylim is not None:
        ax.set_ylim(ylim)

    file = os.path.join(folder, "Summary_"+str(draw)+".fits")

    cm = plt.get_cmap('Spectral') 

    cNorm  = colors.Normalize(vmin=0, vmax=10)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)

    hdulist = fits.open(file)

    IDs = hdulist["META DATA"].data['ID']
    redshift = hdulist['GALAXY PROPERTIES'].data['redshift']

    if xPar["colName"].lower() == "uv_1500_flambda":
        tmp = hdulist[xPar["extName"]].data[xPar["colName"]] / (1.+redshift)
        wl_central = 1500. * (1+redshift)
        xData = -2.5*np.log10(wl_central**2/c_light*tmp) - 48.6
    elif xPar["colName"].lower() == "f160w":
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
        xData = hdulist[xPar["extName"]].data[xPar["colName"]]

    if xLog:
        xData = np.log10(xData)

    yData = hdulist[yPar["extName"]].data[yPar["colName"]]

    if yLog:
        yData = np.log10(yData)

    MM = np.array(ax.get_xlim())

    mean_redshift = np.mean(redshift)

    colorVal = scalarMap.to_rgba(10-mean_redshift)
    indices = range(len(redshift))

    if yPar["colName"] == "Halpha_flux":
        n = len(np.where((xData[indices] < m_UV_lim) & (yData[indices] > Ha_threshold))[0])
        n_tot = len(np.where((xData[indices] < m_UV_lim))[0])

        ax.text(0.08, 0.95-i*0.05, "{:.2f}".format(1.*n/n_tot),
                horizontalalignment='left',
                verticalalignment='center',
                fontsize=14, color=colorVal, 
                weight="heavy",
                transform=ax.transAxes)


    if xPar["colName"] == "M_tot" and "SFR" in yPar["colName"]:
        lgS = mass_SFR_relation(mean_redshift, MM)

        X, Y = np.mgrid[6.:11.:50j, -4.:4.:50j]
        positions = np.vstack([X.ravel(), Y.ravel()])
        n_ = positions.shape[1]
        redshift_ = mean_redshift
        indices_ = np.arange(n_)

        pdf_ = sample_mass_SFR(positions[0,:], positions[1,:], np.full(n_, redshift_))
        wrand = WalkerRandomSampling(pdf_, keys=indices_)
        sel = wrand.random(1000)

        ax.plot(positions[0,sel],
                positions[1,sel],
                marker='o',
                ls="",
                ms=3,
                alpha=0.6,
                color=colorVal,
                zorder=1
                )

        ax.plot(MM,
                lgS,
                color=colorVal,
                lw=2,
                zorder=3)

        ax.plot(xData[indices],
                yData[indices],
                marker='o',
                ls="",
                color=colorVal,
                zorder=2
                )


        ax.text(0.05, 0.95, "{:.2f}".format(mean_redshift),
                horizontalalignment='center',
                verticalalignment='center',
                fontsize=18, color=colorVal, 
                weight="heavy",
                transform=ax.transAxes)

    if yPar["colName"] == "Halpha_flux":
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

    fig.suptitle(str(draw).replace("_", " "), fontsize=20)

    fileName = xPar["colName"]+"_"+yPar["colName"]+"_"+str(draw)+".pdf"
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

    fileName = os.path.join(folder, "Summary_"+str(draw)+".fits")

    hdulist = fits.open(fileName)

    IDs = hdulist['META DATA'].data['ID']
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

    fig.suptitle(str(draw).replace("_", " "), fontsize=20)

    fileName = xPar+"_histogram_"+str(draw)+".pdf"
    fileName = os.path.join(folder, fileName)

    fig.savefig(fileName, dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', papertype='a4', format="pdf",
        transparent=False, bbox_inches="tight", pad_inches=0.1)

    plt.close(fig)



if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-r', '--results-dir',
        help="Directory containing BEAGLE results",
        action="store", 
        type=str, 
        dest="results_dir", 
        required=True
    )

    parser.add_argument(
        '--suffix',
        help="Suffix",
        action="store", 
        type=str, 
        dest="suffix", 
        default=""
    )

    args = parser.parse_args()

    folder = args.results_dir

    # ***************************************************
    # Plotting the mass-SFR relation

    #plot_mass_SFR(draw="MAP", yPar="SFR_100", xlim=[6.5, 11], ylim=[-5, 4])

    for draw in range(10):
        draw = "MC_"+str(draw)+args.suffix

        xPar = {"colName":"M_tot","extName":"GALAXY PROPERTIES"}
        yPar = {"colName":"SFR","extName":"STAR FORMATION"}
        plot_X_vs_Y(draw=draw, xPar=xPar, yPar=yPar, xlim=[7.5, 12], ylim=[-2, 4])

        #plot_mass_SFR(draw=draw, 
        #        xPar="M_tot", yPar="SFR_100", 
        #        xlim=[6.5, 11], ylim=[-3,3])

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



