
##### PS ########################################################
# Signal to noise of ~10 expected for 10 exposures and an input
# flux of 100 nJy => creating a flat input spectrum
# Forcing the summation parameters to those used for the
# sensitivity computation (2,4)

# Creation of the input spectrum (command-line sequence)
import os
import sys
sys.path.append(os.path.join(os.environ['PYP_BEAGLE'], "PyP-BEAGLE"))

import json
import ConfigParser
import argparse
import copy
import beagle_multiprocess
from beagle_utils import get_files_list

from pathos.multiprocessing import ProcessingPool 
from multiprocessing import Manager

from JWSTpylib import c_spectrum
import os
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import matplotlib.colors as colors
from walker_random_sampling import WalkerRandomSampling
import scipy.stats as stats
from collections import OrderedDict
import WeightedKDE
from itertools import repeat

from TwoDimPlot import CredibleIntervalContour

from beagle_filters import PhotometricFilters
from beagle_photometry import ObservedCatalogue
from beagle_utils import extract_row, set_plot_ticks

c_light = 2.99792e+18 # Ang/s
passive_fraction = 0.10
sSFR_threshold = -9.0
M_tot_threshold = 7.3
seed = 12345

Jy = np.float32(1.E-23)
microJy = np.float32(1.E-23 * 1.E-06)
nanoJy = np.float32(1.E-23 * 1.E-09)

write_Summary = True
show_plot = False
write_to_disk = True
compute_lines = False

# ########################################################################################
# ########################################################################################
def write_SED(hdulist, indx, fileName, show_plot=False): 

    # Now write in a separate FITS file the pararmeters corresponding to the MAP row
    new_hdulist = fits.HDUList(fits.PrimaryHDU())

    for hdu in hdulist:
         
        if hdu.data is not None:

            if hdu.is_image:
                new_hdu = fits.PrimaryHDU()
                new_hdu.name = hdu.name
                new_hdu.data = hdu.data[indx,:]
            else:
                new_hdu = fits.BinTableHDU.from_columns(hdu.columns, nrows=1)
                new_hdu.name = hdu.name
                if 'SED WL' in hdu.name:
                    new_hdu.data = hdu.data
                else:
                    new_hdu.data[0] = hdu.data[indx]

            new_hdulist.append(new_hdu)

    new_hdulist.writeto(fileName, clobber=True)

    new_hdulist.close()


# ########################################################################################
# ########################################################################################

def mass_sfr(logM, SFR, z, distribution='t-Student'):
    #Speagle+ 14 - log(M*) = (0.84-0.026*t)logM - (6.51-0.11*t)
    #cosmology used in paper to calculate ages (t) is (h,omega_m,omega_lambda) = (0.7,0.3,0.7)
    #Shivaei+ 15 measure scatter in log(SFR(Halpha))-log(M*) to be 0.3 dex (corrected for uncertainties)
    #We employ cauchy scatter - although I think it should only be scatter in one direction...
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    t = cosmo.age(z).value
    mean = (0.84-0.026*t)*logM - (6.51-0.11*t)

    pdf = np.zeros(len(logM))

    for i, value in enumerate(SFR):

        if distribution.lower() == 't-student':
            pdf[i] = stats.t.pdf(SFR[i], df=3, loc=mean[i], scale=0.3)
        elif distribution.lower() == 'gaussian':
            pdf[i] = stats.norm.pdf(SFR[i], loc=mean[i], scale=0.3)
        elif distribution.lower() == 'cauchy':
            pdf[i] = stats.cauchy.pdf(SFR[i], loc=mean[i], scale=0.3)
        else:
            raise ValueError('The input distribution '+distribution+' is not supported!')

    return pdf

def mass_sfr_test(M_tot, SFR, redshift, distribution='t-Student'):
    #Speagle+ 14 - log(M*) = (0.84-0.026*t)logM - (6.51-0.11*t)
    #cosmology used in paper to calculate ages (t) is (h,omega_m,omega_lambda) = (0.7,0.3,0.7)
    #Shivaei+ 15 measure scatter in log(SFR(Halpha))-log(M*) to be 0.3 dex (corrected for uncertainties)
    #We employ cauchy scatter - although I think it should only be scatter in one direction...
    logM = np.log10(M_tot)
    logSFR = np.log10(SFR)

    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    t = cosmo.age(redshift).value
    mean = (0.84-0.026*t)*logM - (6.51-0.11*t)

    pdf = np.zeros(len(logM))

    for i, value in enumerate(logSFR):

        if distribution.lower() == 't-student':
            pdf[i] = stats.t.pdf(logSFR[i], df=3, loc=mean[i], scale=0.3)
        elif distribution.lower() == 'gaussian':
            pdf[i] = stats.norm.pdf(logSFR[i], loc=mean[i], scale=0.3)
        elif distribution.lower() == 'cauchy':
            pdf[i] = stats.cauchy.pdf(logSFR[i], loc=mean[i], scale=0.3)
        else:
            raise ValueError('The input distribution '+distribution+' is not supported!')

    return pdf

# ########################################################################################
# ########################################################################################

def print_ID(IDs, fileNames):

    print "\nID: ", ID


def extract_SED_from_FITS(ID, fileName, count):

    commonPrefix = ID + suffix

    recompute = True

    print "\n fileName: ", fileName

    # The input BEAGLE FITS file contains several (thousands) of SED in the
    # "FULL SED" extension, we will just select the one corresponding to the
    # Maximum-a-Posteriori to create the simulated NIRSpec observation

    # *********************************************
    # Create the "MAP" FITS file
    # *********************************************

    #for i in range(n_MonteCarlo):

        # This FITS file will only contain the row corresponding to the MAP value of the posterior PDF
     #   objPrefix = commonPrefix + '_MC_' + str(i)
      #  newfile_name = os.path.join(output_folder, objPrefix + '.fits.gz')
       # print "newfile_name: ", newfile_name

        #if not os.path.isfile(newfile_name):
         #   recompute = True
          #  break

    if recompute:

        # Open the original BEAGLE FITS file
        hdulist = fits.open(os.path.join(folder, fileName))

        # Get the posterior probability
        post = hdulist['posterior pdf'].data['probability']

        # Now select in sSFR, so as to have a certain fraction of passive and star forming galaxies
        sSFR = hdulist['star formation'].data['sSFR']
        SFR = hdulist['star formation'].data['SFR']
        SFR_100 = hdulist['star formation'].data['SFR_100']
        M_star = hdulist['galaxy properties'].data['M_star']
        M_tot = hdulist['galaxy properties'].data['M_tot']
        template_redshifts = hdulist['posterior pdf'].data['redshift']

        if compute_lines:
            wl = np.ravel(hdulist['FULL SED WL'].data[0][:])

        dropout_redshift = np.float32(obs_catalogue['dropout_redshift'][obs_catalogue['ID']==ID])

        #pdf = mass_sfr(np.log10(M_tot), np.log10(SFR), template_redshifts, distribution='gaussian')
        pdf = mass_sfr(np.log10(M_tot), np.log10(SFR), template_redshifts, distribution='gaussian')

        data = np.zeros((2,len(M_tot)))
        data[0,:] = M_tot
        data[1,:] = SFR
        #kde_pdf = WeightedKDE.gaussian_kde(, weights=pdf)

        #reweighted_pdf = post*pdf
        reweighted_pdf = pdf

        # Redshifts corresponding to the different Beagle solutions
        #template_redshifts = hdulist['posterior pdf'].data['redshift']

        # Pick only the solutions consistent with Bouwens selection, i.e. +/- 0.5 from the redshift in the catalogue

        z_low = dropout_redshift - 0.5
        z_up = dropout_redshift + 0.5
        ok = np.where((template_redshifts >= z_low) & (template_redshifts <= z_up) & (np.log10(M_tot) >= M_tot_threshold))[0]

        # Extend the search at +/- 1 from the redshift in the catalogue
        while (len(ok) < 2*n_MonteCarlo):
            z_low -= redshift_step
            z_up += redshift_step
            ok = np.where((template_redshifts >= z_low) & (template_redshifts <= z_up) & (np.log10(M_tot) >= M_tot_threshold))[0]

        indices = np.arange(len(post))

        wrand = WalkerRandomSampling(reweighted_pdf[ok], keys=indices[ok])
        rows = wrand.random(n_MonteCarlo)

        selIndices = rows

        n = len(selIndices)
        summary = OrderedDict()
        summary['ID'] = np.empty(n, dtype='S20')
        summary['row'] = np.zeros(n, dtype=np.int)
        summary['z'] = np.zeros(n, dtype=np.float32)
        summary['M_star'] = np.zeros(n, dtype=np.float32)
        summary['M_tot'] = np.zeros(n, dtype=np.float32)
        summary['SFR'] = np.zeros(n, dtype=np.float32)
        summary['sSFR'] = np.zeros(n, dtype=np.float32)
        summary['SFR_100'] = np.zeros(n, dtype=np.float32)
        summary['posterior_PDF'] = np.zeros(n, dtype=np.float32)
        summary['reweighted_PDF'] = np.zeros(n, dtype=np.float32)
        summary['UV_1500_FLAMBDA'] = np.zeros(n, dtype=np.float32)

        # Compute Ha and Hb integrated fluxes and EW
        for key, value in data_lines.iteritems():
            summary[key+"_flux"] = np.zeros(n, dtype=np.float32)
            summary[key+"_EW"] = np.zeros(n, dtype=np.float32)

        for i, indx in enumerate(selIndices):

            if write_Summary:
                summary['ID'][i] = ID
                summary['row'][i] = indx
                summary['z'][i] = template_redshifts[indx]
                summary['M_star'][i] = M_star[indx]
                summary['M_tot'][i] = M_tot[indx]
                summary['SFR'][i] = SFR[indx]
                summary['sSFR'][i] = sSFR[indx]
                summary['SFR_100'][i] = SFR_100[indx]
                summary['posterior_PDF'][i] = post[indx]
                summary['reweighted_PDF'][i] = reweighted_pdf[indx]

                if compute_lines:

                    SED = hdulist['FULL SED'].data[indx,:]

                    i0 = np.searchsorted(wl, 1450)
                    i1 = np.searchsorted(wl, 1550)
                    UV_1500 = np.trapz(SED[i0:i1+1], x=wl[i0:i1+1]) / (wl[i1]-wl[i0])
                    summary['UV_1500_FLAMBDA'][i] = UV_1500

                    for key, value in data_lines.iteritems():

                        # Compute the flux integrated around the line center, for all selected rows
                        i0 = np.searchsorted(wl, value["center"]-width)
                        i1 = np.searchsorted(wl, value["center"]+width)

                        flux = np.trapz(SED[i0:i1+1], x=wl[i0:i1+1])
                        summary[key+"_flux"][i] = flux

                        # To compute the EW, you firstly computed the integrated flux in a window on the left of the EL
                        il0 = np.searchsorted(wl, value["left_cont"][0])
                        il1 = np.searchsorted(wl, value["left_cont"][1])
                        if il0 == il1:
                            il0 -= 1
                        flux_left = np.trapz(SED[il0:il1+1], x=wl[il0:il1+1]) / (wl[il1]-wl[il0])

                        # Repeat the same on the right of the line
                        ir0 = np.searchsorted(wl, value["right_cont"][0])
                        ir1 = np.searchsorted(wl, value["right_cont"][1])
                        if ir0 == ir1:
                            ir1 += 1
                        flux_right = np.trapz(SED[ir0:ir1+1], x=wl[ir0:ir1+1]) / (wl[ir1]-wl[ir0])

                        flux /= 0.5*(flux_left+flux_right)
                        flux -= 1.
                        flux = -flux

                        summary[key+"_EW"][i] = flux

            if write_to_disk:
                objPrefix = commonPrefix + '_MC_' + str(i)
                newfile_name = os.path.join(output_folder, objPrefix + '.fits.gz')
                print "newfile_name: ", newfile_name
                write_SED(hdulist, indx, newfile_name)


        hdulist.close()

        if write_Summary:

            return summary
                   
def extract_MAP_from_FITS(ID, fileName):

    commonPrefix = ID + suffix

    recompute = True

    print "\n fileName: ", fileName

    # The input BEAGLE FITS file contains several (thousands) of SED in the
    # "FULL SED" extension, we will just select the one corresponding to the
    # Maximum-a-Posteriori to create the simulated NIRSpec observation

    # *********************************************
    # Create the "MAP" FITS file
    # *********************************************

    #for i in range(n_MonteCarlo):

        # This FITS file will only contain the row corresponding to the MAP value of the posterior PDF
     #   objPrefix = commonPrefix + '_MC_' + str(i)
      #  newfile_name = os.path.join(output_folder, objPrefix + '.fits.gz')
       # print "newfile_name: ", newfile_name

        #if not os.path.isfile(newfile_name):
         #   recompute = True
          #  break

    if recompute:

        # Open the original BEAGLE FITS file
        hdulist = fits.open(os.path.join(folder, fileName))

        # Get the posterior probability
        post = hdulist['posterior pdf'].data['probability']

        # Now select in sSFR, so as to have a certain fraction of passive and star forming galaxies
        sSFR = hdulist['star formation'].data['sSFR']
        SFR = hdulist['star formation'].data['SFR']
        SFR_100 = hdulist['star formation'].data['SFR_100']
        M_star = hdulist['galaxy properties'].data['M_star']
        M_tot = hdulist['galaxy properties'].data['M_tot']
        template_redshifts = hdulist['posterior pdf'].data['redshift']

        wl = np.ravel(hdulist['FULL SED WL'].data[0][:])

        # Redshifts corresponding to the different Beagle solutions
        #template_redshifts = hdulist['posterior pdf'].data['redshift']

        # Pick only the solutions consistent with Bouwens selection, i.e. +/- 0.5 from the redshift in the catalogue
        dropout_redshift = np.float32(obs_catalogue['dropout_redshift'][obs_catalogue['ID']==ID])

        z_low = dropout_redshift - 0.5
        z_up = dropout_redshift + 0.5
        ok = np.where((template_redshifts >= z_low) & (template_redshifts <= z_up))[0]

        # Extend the search at +/- 1 from the redshift in the catalogue
        while (len(ok) < 1):
            z_low -= redshift_step
            z_up += redshift_step
            ok = np.where((template_redshifts >= z_low) & (template_redshifts <= z_up))[0]

        indices = np.arange(len(post))
        indx = np.argmax(post[ok])
        MAP_indx = indices[ok[indx]]

        summary = OrderedDict()

        # Compute Ha and Hb integrated fluxes and EW
        for key, value in data_lines.iteritems():
            summary[key+"_flux"] = np.zeros(n, dtype=np.float32)
            summary[key+"_EW"] = np.zeros(n, dtype=np.float32)

        if write_Summary:
            summary['ID'] = ID
            summary['row'] = MAP_indx
            summary['z'] = template_redshifts[MAP_indx]
            summary['M_star'] = M_star[MAP_indx]
            summary['M_tot'] = M_tot[MAP_indx]
            summary['SFR'] = SFR[MAP_indx]
            summary['sSFR'] = sSFR[MAP_indx]
            summary['SFR_100'] = SFR_100[MAP_indx]
            summary['posterior_PDF'] = post[MAP_indx]

            SED = hdulist['FULL SED'].data[MAP_indx,:]

            i0 = np.searchsorted(wl, 1450)
            i1 = np.searchsorted(wl, 1550)
            UV_1500 = np.trapz(SED[i0:i1+1], x=wl[i0:i1+1]) / (wl[i1]-wl[i0])
            summary['UV_1500_FLAMBDA'] = UV_1500

            for key, value in data_lines.iteritems():

                # Compute the flux integrated around the line center, for all selected rows
                i0 = np.searchsorted(wl, value["center"]-width)
                i1 = np.searchsorted(wl, value["center"]+width)

                flux = np.trapz(SED[i0:i1+1], x=wl[i0:i1+1])
                summary[key+"_flux"] = flux

                # To compute the EW, you firstly computed the integrated flux in a window on the left of the EL
                il0 = np.searchsorted(wl, value["left_cont"][0])
                il1 = np.searchsorted(wl, value["left_cont"][1])
                if il0 == il1:
                    il0 -= 1
                flux_left = np.trapz(SED[il0:il1+1], x=wl[il0:il1+1]) / (wl[il1]-wl[il0])

                # Repeat the same on the right of the line
                ir0 = np.searchsorted(wl, value["right_cont"][0])
                ir1 = np.searchsorted(wl, value["right_cont"][1])
                if ir0 == ir1:
                    ir1 += 1
                flux_right = np.trapz(SED[ir0:ir1+1], x=wl[ir0:ir1+1]) / (wl[ir1]-wl[ir0])

                flux /= 0.5*(flux_left+flux_right)
                flux -= 1.
                flux = -flux

                summary[key+"_EW"] = flux

            if write_to_disk:
                objPrefix = commonPrefix + '_MAP' 
                newfile_name = os.path.join(output_folder, objPrefix + '.fits.gz')
                print "newfile_name: ", newfile_name
                write_SED(hdulist, MAP_indx, newfile_name)


        hdulist.close()

        if write_Summary:

            return summary
                   

def plot_SEDs(ID, hdulist, rows_indices, x_log=False, print_title=True):
        
        # From the (previously loaded) observed catalogue select the row
        # corresponding to the input ID
        observation = extract_row(observed_catalogue.data, ID)

        # Check if you need to apply an aperture correction to the catalogue fluxes
        if 'aper_corr' in observed_catalogue.data.dtype.names:
            aper_corr = 10.**(-0.4*observation[0]['aper_corr'])
        else:
            aper_corr = 1.

        # Put observed photometry and its error in arrays
        obs_flux, obs_flux_err = observed_catalogue.extract_fluxes(filters, ID)
        ok = np.where(obs_flux_err > 0.)[0]
        obs_flux[ok] *= 1.E+09
        obs_flux_err[ok] *= 1.E+09

        if x_log:
            wl_eff = np.log10(filters.data['wl_eff'])
        else:
            wl_eff = np.array(filters.data['wl_eff'])

        # Sort wl_eff array
        sor = np.argsort(wl_eff)
        obs_flux, obs_flux_err, wl_eff = obs_flux[sor], obs_flux_err[sor], wl_eff[sor]
        ok = np.where(obs_flux_err > 0.)[0]

        # Initialize figure
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        # Define axes labels
        if x_log:
            ax.set_xlabel("$\\log (\lambda_\\textnormal{eff} / \\textnormal{\AA}$ (observed-frame))")
        else:
            ax.set_xlabel("$\lambda_\\textnormal{eff} / \\textnormal{\AA}$ (observed-frame)")

        ax.set_ylabel("$f_{\\nu}/\\textnormal{nanoJy}$")

        # Set better location of tick marks
        set_plot_ticks(ax, n_x=5)

        kwargs = {'alpha':0.8}
        
        dwl = wl_eff[-1] - wl_eff[0]
        ax.set_xlim(wl_eff[0]-dwl*0.1, wl_eff[-1]+dwl*0.1)

        # Plot the data with errors bars
        plt.errorbar(wl_eff[ok], 
                obs_flux[ok], 
                yerr = obs_flux_err[ok],
                color = "dodgerblue",
                ls = " ",
                marker = "D",
                markeredgewidth = 0.,
                markersize = 8,
                elinewidth=1.0,
                capsize=3,
                **kwargs)

        cm = plt.get_cmap('Spectral') 
        cNorm  = colors.Normalize(vmin=0, vmax=len(rows_indices))
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)

        model_mag = np.zeros(filters.n_bands)
        # Plot the model
        model_sed = hdulist['marginal photometry']
        for j, row in enumerate(rows_indices):
            color = scalarMap.to_rgba(len(rows_indices)-j)

            for i, band_name in enumerate((filters.data['label'][sor])):
                model_mag[i] = model_sed.data['_'+band_name+'_'][row] / nanoJy

            plt.scatter(wl_eff[ok],
                    model_mag[ok],
                    color = color,
                    marker = "*",
                    s = 34,
                    lw=0.2,
                    edgecolor='black',
                    alpha = 0.7,
                    zorder=3
                    )

            chi_square = np.sum(((model_mag[ok]-obs_flux[ok]) / obs_flux_err[ok])**2)

            plt.scatter(0.03,
                    0.95-j*0.03,
                    color = color,
                    marker = "*",
                    s = 34,
                    lw=0.2,
                    edgecolor='black',
                    transform=ax.transAxes)

            chi_square = np.sum(((model_mag[ok]-obs_flux[ok]) / obs_flux_err[ok])**2)


            ax.text(0.04, 0.95-j*0.03, "{:.2f}".format(chi_square), 
                    horizontalalignment='left',
                    verticalalignment='center',
                    fontsize=10, weight='heavy', 
                    color='black', 
                    transform=ax.transAxes)


        # Title of the plot is the object ID
        if print_title: plt.title(str(ID))

        name = os.path.join(plots_folder, str(ID)+'_SED.pdf')

        fig.savefig(name, dpi=None, facecolor='w', edgecolor='w',
                orientation='portrait', papertype='a4', format="pdf",
                transparent=False, bbox_inches="tight", pad_inches=0.1)

        plt.close(fig)

def draw_rows_from_posterior(ID, fileName, 
        n_samples=1,
        params_ranges=None,
        weight_func=None, 
        weight_func_args=None,
        extensions=None,
        make_plot=False):

    print "Extracting rows from object ID: ", ID

    # Open the original BEAGLE FITS file
    hdulist = fits.open(fileName)

    # Get the posterior probability
    post = hdulist['posterior pdf'].data['probability']

    weights = np.full(1, len(post), dtype=np.float32)

    # Can pass a "weight_func" that will take as input some parameters of the Beagle output and compute a weight, that will be used to multiply the posterior PDF
    if weight_func is not None:
        func_args = OrderedDict()
        for key, value in weight_func_args.iteritems():
            # Load the data from the Beagle FITS file into a dictionary
            if "colName" in value:
                data = hdulist[value["extName"]].data[value["colName"]]
                func_args[key] = data
            # Assume that the other entries of the dictionary correspond to kwargs of the weight_func
            else:
                func_args[key] = value
    
        # Compute the actual weights
        weights = weight_func(**func_args)

    reweighted_pdf = post*weights
    #reweighted_pdf = weights
    #reweighted_pdf = post
    
    # Besides the weight function, the user can select some allowed ranges for some parameters
    mask = np.ones(len(post), dtype=bool)
    n = -1
    if params_ranges is not None:
        i=0
        while n < 2*n_samples:
            mask = np.ones(len(post), dtype=np.bool)
            for key, value in params_ranges.iteritems():
                data = hdulist[value["extName"]].data[value["colName"]]
                # The "step" is used to widen the range at each iteration, in
                # order to have enough valid elements to then randomly draw
                # n_samples from them
                if "step" in value:
                    step = value["step"]
                else:
                    step = 0
                if "min" in value:
                    m = value["min"]-i*step
                    indx = np.where(data < m)[0]
                    if len(indx) > 0:
                        mask[indx] = False
                if "max" in value:
                    m = value["max"]+i*step
                    indx = np.where(data > m)[0]
                    if len(indx) > 0:
                        mask[indx] = False

            n = np.sum(mask)
            i += 1
            
    
    #data = np.zeros((2,len(post)))
    #data[0,:] = hdulist['galaxy properties'].data['M_tot']
    #data[1,:] = hdulist['star formation'].data['SFR']
    #credInterval = CredibleIntervalContour(data, post)
    
    #isoLevels = credInterval.GetProbabilityForCredibleRegion((0.68,0.95))
    # All probabilities larger than than the isocontour probability are ok, hence we select those smaller, to set then mask=False for those indices
    #indx = np.where(post < isoLevels[0])[0]
    #mask[indx] = False
    
    indices = np.arange(len(post))

    # Randomly draw the indices of the rows, using as weights the reweighted posterior PDF 
    wrand = WalkerRandomSampling(reweighted_pdf[mask], keys=indices[mask])
    rows_indices = wrand.random(n_samples)

    #print "---> ", post[rows_indices]

    if make_plot: 
        plot_SEDs(ID, hdulist, rows_indices)

    data = OrderedDict()
    data["rows_indices"] = rows_indices
    if extensions is not None:
        for ext in extensions:
            data[ext] = hdulist[ext].data[rows_indices]

    hdulist.close()

    return data


# ########################################################################################
# ########################################################################################

observed_catalogue = ObservedCatalogue()
filters = PhotometricFilters()
plots_folder = ""

if __name__ == '__main__':

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
            '--plot', 
            dest='plot', 
            action='store_true')

    parser.add_argument(
        '-p', '--parameter-file',
        help="Parameter file used in the BEAGLE run",
        action="store", 
        type=str, 
        dest="param_file"
    )

    # Number of processors to use in the multi-processor parts of the analysis
    parser.add_argument(
        '-n_samples',
        help="Number of Monte Carlo draws",
        action="store", 
        type=int, 
        dest="n_samples",
        default=1
    )

    parser.add_argument(
        '-n_objects',
        help="Number of objects",
        action="store", 
        type=int, 
        dest="n_objects",
        default=-1
    )

    parser.add_argument(
        '-np',
        help="Number of processors to use",
        action="store", 
        type=int, 
        dest="np",
        default=-1
    )

    parser.add_argument(
        '--folder-suffix',
        help="Suffix for output folder",
        action="store", 
        type=str, 
        dest="folder_suffix",
        default=""
    )

    parser.add_argument(
        '--suffix',
        help="Suffix for output filename",
        action="store", 
        type=str, 
        dest="suffix",
        default=""
    )

    parser.add_argument(
        '--params-ranges',
        help="JSON string containing the hard limits on the params",
        action="store", 
        type=str, 
        dest="params_ranges"
    )

    # Get parsed arguments
    args = parser.parse_args()    

    # Initialize seed for random number generator
    ranSeed = 123456
    np.random.seed(ranSeed)

    # Read parameter file
    make_plot=False
    if args.plot:
        make_plot=True
        config = ConfigParser.SafeConfigParser()
        name = os.path.join(args.results_dir, 'BEAGLE-input-files', args.param_file) 
        config.read(name)

        # Load filters file
        filters_file = os.path.expandvars(config.get('main', 'FILTERS FILE'))
        filters.load(filters_file)

        # Load observed catalogue
        file_name = os.path.expandvars(config.get('main', 'PHOTOMETRIC CATALOGUE'))
        observed_catalogue.load(file_name)

    # Get list of Beagle output FITS files in the esults folder, and extract the IDs from the files
    file_list, IDs = get_files_list(results_dir=args.results_dir)
    file_list = [os.path.join(args.results_dir, file) for file in file_list]

    # Finally, write the FITS tables, one per draw
    fold = args.results_dir.split('BEAGLE_results')[1][1:] + args.folder_suffix
    folder = os.path.join('/home/jchevall/JWST/Simulations', fold)
    print "folder:", folder
    if not os.path.exists(folder):
        os.makedirs(folder)

    plots_folder = os.path.join(folder, 'plots')
    if not os.path.exists(plots_folder):
        os.makedirs(plots_folder)

    if args.n_objects < 0:
        args.n_objects = len(IDs)

    weight_func_args = OrderedDict()
    weight_func_args['M_tot'] = {"colName":"M_tot", "extName":"galaxy properties"}
    weight_func_args['SFR'] = {"colName":"SFR", "extName":"star formation"}
    weight_func_args['redshift'] = {"colName":"redshift", "extName":"galaxy properties"}
    weight_func_args['distribution'] = 'gaussian'

    params_ranges=None
    if args.params_ranges is not None:
        params_ranges = json.loads(args.params_ranges)

    #print "params_ranges: ", params_ranges

    #stop

    #params_ranges = OrderedDict()

    #params_ranges['redshift'] = {"colName":"redshift", "extName":"galaxy properties", "min":3.5, "max":4.5, "step":0.1}

    #params_ranges['redshift'] = {"colName":"redshift", "extName":"galaxy properties", "min":5.0, "step":0.1}

    #params_ranges['mass'] = {"colName":"M_tot", "extName":"galaxy properties", "min":10.**(7.3)}


    #
    extensions = ('GALAXY PROPERTIES', 'STAR FORMATION', 'POSTERIOR PDF')

    if args.np <= 0:
        results = list()
        for (ID, file) in zip(IDs[0:args.n_objects], file_list[0:args.n_objects]):
            res = draw_rows_from_posterior(ID, file, 
                    args.n_samples,
                    params_ranges,
                    mass_sfr_test,
                    weight_func_args,
                    extensions,
                    make_plot
                    )
            results.append(res)
    else:

        # Set number of parellel processes to use
        pool = ProcessingPool(nodes=args.np)

        # Launch the actual calculation on multiple processesors
        results = pool.map(draw_rows_from_posterior, IDs, file_list, 
                (args.n_samples,)*len(IDs),
                (params_ranges,)*len(IDs),
                (mass_sfr_test,)*len(IDs),
                (weight_func_args,)*len(IDs),
                (extensions,)*len(IDs)
                )

    # Initialize an empty HDU
    hdulist = fits.HDUList(fits.PrimaryHDU())

    # Create columns to hold the mata data, i.e. the IDs of the objects and index of the rows extracted
    # NB: the row index starts from 0
    c1 = fits.Column(name='ID', format='20A')
    c2 = fits.Column(name='row_index', format='I')
    coldefs = fits.ColDefs([c1, c2])
    t1 = fits.BinTableHDU.from_columns(coldefs, nrows=len(IDs))
    t1.name = 'META DATA'

    hdulist.append(t1)

    # Create extension to hold the posterior pdf data    
    for key, value in results[0].iteritems():
        if 'rows_indices' not in key:
            t2 = fits.BinTableHDU.from_columns(value.columns, nrows=len(IDs))
            t2.name = key
            hdulist.append(t2)

    data = OrderedDict()

    # Parse the results
    # Cycle over all objects
    for j, res in enumerate(results):
        # Cycle over all rows extracted for each object
        for key, value in res.iteritems():
            for i, d in enumerate(value):
                draw_key = 'draw_'+str(i)
                if draw_key not in data:
                    data[draw_key] = copy.deepcopy(hdulist)

                if 'rows_indices' in key:
                    data[draw_key]["META DATA"].data['ID'][j] = IDs[j]
                    data[draw_key]["META DATA"].data['row_index'][j] = d
                else:
                    data[draw_key][key].data[j] = d


    for i, (key, value) in enumerate(data.iteritems()):
        name = os.path.join(folder, 'Summary_MC_'+str(i)+args.suffix+'.fits')
        value.writeto(name, clobber=True)
    stop

# Name of the original XDF catalogue on which the BEAGLE fittting is based
catalogue_name = os.path.expandvars("$BEAGLE_DATA/XDF/XDF_DROPOUTS_aper_corr.fits")
obs_hdu = fits.open(catalogue_name)[1]
obs_catalogue = obs_hdu.data

fit_folder = "ineb_Jan16_logU_xid_delayed_SFR_max_age-Gaussian"
fit_folder = "B_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian"

# Folder containing the input BEAGLE FITS files that will be post-processed to
# obtain simulated NIRSpec observations 
folder = os.path.join("/local/jchevall/BEAGLE_results/XDF", fit_folder)
#folder="/user_data/jchevall/BEAGLE/files/results/XDF/ineb_Jan16_logU_xid_constant_SFR"

suffix = "_BEAGLE"

# Files in this list will be processed into NIRSpec observations
IDs = list()
fileNames = list()
fileName = os.path.join(folder, 'selected_objects.txt')
if os.path.isfile(fileName):
    with open(fileName) as f:
        for line in f:
            objFileName = line.strip()
            fileNames.append(objFileName)
            IDs.append(objFileName.split(suffix)[0])
else:
    fileNames, IDs = get_files_list(results_dir=folder)

n = len(IDs)

data_lines = OrderedDict()

line = {"center":6563., "left_cont":[6535., 6540.], "right_cont":[6590., 6595.]}
data_lines['Halpha'] = line

line = {"center":4861., "left_cont":[4820., 4825.], "right_cont":[4890., 4895.]}
data_lines['Hbeta'] = line

width = 15.

# Output folder, containing the simulated NIRSpec observations
output_folder = os.path.join("/home/jchevall/JWST/Simulations/XDF", fit_folder, "mass_SFR_weighted_Gaussian_no_posterior")
if not os.path.isdir(output_folder):
    os.makedirs(output_folder)

# Number of Monte Carlo draws
n_MonteCarlo = 1
redshift_step = 0.1
ranSeed = 123456
np.random.seed(ranSeed)


if write_Summary:
    summary_data = OrderedDict()
    summary_data['ID'] = np.empty(n, dtype='S20')
    summary_data['row'] = np.zeros(n, dtype=np.int)
    summary_data['z'] = np.zeros(n, dtype=np.float32)
    summary_data['M_star'] = np.zeros(n, dtype=np.float32)
    summary_data['M_tot'] = np.zeros(n, dtype=np.float32)
    summary_data['SFR'] = np.zeros(n, dtype=np.float32)
    summary_data['sSFR'] = np.zeros(n, dtype=np.float32)
    summary_data['SFR_100'] = np.zeros(n, dtype=np.float32)
    summary_data['posterior_PDF'] = np.zeros(n, dtype=np.float32)
    summary_data['reweighted_PDF'] = np.zeros(n, dtype=np.float32)
    summary_data['UV_1500_FLAMBDA'] = np.zeros(n, dtype=np.float32)

    # Compute Ha and Hb integrated fluxes and EW
    for key, value in data_lines.iteritems():
        summary_data[key+"_flux"] = np.zeros(n, dtype=np.float32)
        summary_data[key+"_EW"] = np.zeros(n, dtype=np.float32)

    global_summary = OrderedDict()
    for i in range(n_MonteCarlo):
        key = 'MC_' + str(i)
        print "key: ", key
        global_summary[key] = copy.deepcopy(summary_data)


counters = range(len(IDs))    

pool = ProcessingPool(nodes=23)
#results = extract_MAP_from_FITS(IDs[0], fileNames[0])
if 1 == 0:
    results = pool.map(extract_MAP_from_FITS, IDs, fileNames)
    if write_Summary:
        for j, summary in enumerate(results):
            for k, val in summary.iteritems():
                global_summary[key][k][j] = val

        fileName = 'Summary_MAP.fits'
        fileName = os.path.join(output_folder, fileName)
        tab = Table(global_summary[key])
        tab.write(fileName, format='fits', overwrite=True)

if 0 == 0:
    results = pool.map(extract_SED_from_FITS, IDs, fileNames, counters)

    if write_Summary:
        for i in range(n_MonteCarlo):
            key = 'MC_' + str(i)
            for j, summary in enumerate(results):
                for k, val in summary.iteritems():
                    global_summary[key][k][j] = val[i]

            fileName = 'Summary_' + key + '.fits'
            fileName = os.path.join(output_folder, fileName)
            tab = Table(global_summary[key])
            tab.write(fileName, format='fits', overwrite=True)

