
##### PS ########################################################
# Signal to noise of ~10 expected for 10 exposures and an input
# flux of 100 nJy => creating a flat input spectrum
# Forcing the summation parameters to those used for the
# sensitivity computation (2,4)

# Creation of the input spectrum (command-line sequence)
import sys
sys.path.append("/home/jchevall/Coding/PyP-BEAGLE/PyP-BEAGLE")

import copy
import beagle_multiprocess

from pathos.multiprocessing import ProcessingPool 
from multiprocessing import Manager

from JWSTpylib import c_spectrum
import os
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from walker_random_sampling import WalkerRandomSampling
import scipy.stats as stats
from collections import OrderedDict

c_light = 2.99792e+18 # Ang/s
passive_fraction = 0.10
sSFR_threshold = -9.0
seed = 12345

write_Summary = True
show_plot = False
write_to_disk = True

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
        M_star = hdulist['galaxy properties'].data['M_star']
        template_redshifts = hdulist['posterior pdf'].data['redshift']

        wl = np.ravel(hdulist['FULL SED WL'].data[0][:])

        pdf = mass_sfr(np.log10(M_star), np.log10(SFR), template_redshifts, distribution='gaussian')

        reweighted_pdf = post*pdf

        # Redshifts corresponding to the different Beagle solutions
        #template_redshifts = hdulist['posterior pdf'].data['redshift']

        # Pick only the solutions consistent with Bouwens selection, i.e. +/- 0.5 from the redshift in the catalogue
        dropout_redshift = np.float32(obs_catalogue['dropout_redshift'][obs_catalogue['ID']==ID])

        z_low = dropout_redshift - 0.5
        z_up = dropout_redshift + 0.5
        ok = np.where((template_redshifts >= z_low) & (template_redshifts <= z_up))[0]

        # Extend the search at +/- 1 from the redshift in the catalogue
        while (len(ok) < 2*n_MonteCarlo):
            z_low -= redshift_step
            z_up += redshift_step
            ok = np.where((template_redshifts >= z_low) & (template_redshifts <= z_up))[0]

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
        summary['SFR'] = np.zeros(n, dtype=np.float32)
        summary['sSFR'] = np.zeros(n, dtype=np.float32)
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
                summary['SFR'][i] = SFR[indx]
                summary['sSFR'][i] = sSFR[indx]
                summary['posterior_PDF'][i] = post[indx]
                summary['reweighted_PDF'][i] = reweighted_pdf[indx]

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
                   



# ########################################################################################
# ########################################################################################


# Name of the original XDF catalogue on which the BEAGLE fittting is based
catalogue_name = os.path.expandvars("$BEAGLE_DATA/XDF/XDF_DROPOUTS_aper_corr.fits")
obs_hdu = fits.open(catalogue_name)[1]
obs_catalogue = obs_hdu.data

# Folder containing the input BEAGLE FITS files that will be post-processed to
# obtain simulated NIRSpec observations 
folder = "/local/jchevall/BEAGLE_results/XDF/ineb_Jan16_logU_xid"

suffix = "_BEAGLE"

# Files in this list will be processed into NIRSpec observations
fileName = os.path.join(folder, 'selected_objects.txt')

data_lines = OrderedDict()

line = {"center":6563., "left_cont":[6535., 6540.], "right_cont":[6590., 6595.]}
data_lines['Halpha'] = line

line = {"center":4861., "left_cont":[4820., 4825.], "right_cont":[4890., 4895.]}
data_lines['Hbeta'] = line

width = 15.

# Output folder, containing the simulated NIRSpec observations
output_folder = "/home/jchevall/JWST/Simulations/Sep_2016/mass_SFR_weighted_Gaussian"
if not os.path.isdir(output_folder):
    os.makedirs(output_folder)

# Number of Monte Carlo draws
n_MonteCarlo = 10
redshift_step = 0.1
ranSeed = 123456
np.random.seed(ranSeed)

num_lines = sum(1 for line in open(fileName))
n = num_lines

if write_Summary:
    summary_data = OrderedDict()
    summary_data['ID'] = np.empty(n, dtype='S20')
    summary_data['row'] = np.zeros(n, dtype=np.int)
    summary_data['z'] = np.zeros(n, dtype=np.float32)
    summary_data['M_star'] = np.zeros(n, dtype=np.float32)
    summary_data['SFR'] = np.zeros(n, dtype=np.float32)
    summary_data['sSFR'] = np.zeros(n, dtype=np.float32)
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

IDs = list()
fileNames = list()

f = open(fileName, 'r')
for line in f:

    objFileName = line.strip()
    fileNames.append(objFileName)
    IDs.append(objFileName.split(suffix)[0])

counters = range(len(IDs))    

pool = ProcessingPool(nodes=23)
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
        tab.write(fileName, format='fits')

