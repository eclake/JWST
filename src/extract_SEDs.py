#!/usr/bin/env python

##### PS ########################################################
# Signal to noise of ~10 expected for 10 exposures and an input
# flux of 100 nJy => creating a flat input spectrum
# Forcing the summation parameters to those used for the
# sensitivity computation (2,4)

# Creation of the input spectrum (command-line sequence)
import matplotlib
matplotlib.use('GTKAgg')
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

from credible_intervals import CredibleInterval

from beagle_filters import PhotometricFilters
from beagle_photometry import ObservedCatalogue
from beagle_utils import extract_row, set_plot_ticks

c_light = 2.99792e+18 # Ang/s
nanoJy = np.float32(1.E-23 * 1.E-09)

show_plot = False

# taken from https://github.com/eclake/WG3_NIRSpec/blob/a3ade5286aea8d73bc07a6dbbfd24969a4f6eab2/SF%20Q%20separation.ipynb
# Kevin Hainline employed the Whitaker+ 11 colour-cut criteria
def UVJ_separation(redshift, UV, VJ): # returns 1 for quiescent, 0 for star-forming
    if redshift > 4:
        return -1
    if redshift < 0.5:
        if (UV > 0.88*VJ+0.69 and UV > 1.3 and VJ < 1.6):
            #print 'yes'
            return 1
        else:
            return 0
    if redshift >= 0.5 and redshift < 1.5:
        if (UV > 0.88*VJ+0.59 and UV > 1.3 and VJ < 1.6):
            #print 'yes'
            return 1
        else:
            return 0
    if redshift >= 1.5 and redshift < 2.0:
        if (UV > 0.88*VJ+0.59 and UV > 1.3 and VJ < 1.5):
            #print 'yes'
            return 1
        else:
            return 0
    if redshift >= 2 and redshift < 4:
        if (UV > 0.88*VJ+0.59 and UV > 1.2 and VJ < 1.4):
            #print 'yes'
            return 1
        else:
            return 0



# ****************************************************************************************
def draw_rows_from_interval(data, probability, n_draws=1, level=0.68):

    CredInterv = CredibleInterval(data=data, probability=probability)
    indices = np.arange(len(probability))

    if data.ndim == 1:
        region = CredInterv.Get1DCredibleRegion(levels=(level,))
        ok = np.where((data >= region[0]) & (data <= region[1]))[0] 
    elif data.ndim == 2:
        prob_level = CredInterv.GetProbabilityFor2DCredibleRegion(levels=(level,))
        kde_pdf_grid =  CredInterv.kde_pdf(data)
        ok = np.where((kde_pdf_grid >= prob_level))[0]

    return np.random.choice(indices[ok], size=n_draws)

# ****************************************************************************************
def mass_sfr(mass, SFR, z, distribution='t-Student'):
    #Speagle+ 14 - log(M*) = (0.84-0.026*t)logM - (6.51-0.11*t)
    #cosmology used in paper to calculate ages (t) is (h,omega_m,omega_lambda) = (0.7,0.3,0.7)
    #Shivaei+ 15 measure scatter in log(SFR(Halpha))-log(M*) to be 0.3 dex (corrected for uncertainties)
    #We employ cauchy scatter - although I think it should only be scatter in one direction...
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    t = cosmo.age(z).value
    mean = (0.84-0.026*t)*mass - (6.51-0.11*t)

    pdf = np.zeros(len(mass))

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

# ****************************************************************************************
def mass_sfr_test(mass, SFR, redshift, distribution='t-Student'):
    #Speagle+ 14 - log(M*) = (0.84-0.026*t)logM - (6.51-0.11*t)
    #cosmology used in paper to calculate ages (t) is (h,omega_m,omega_lambda) = (0.7,0.3,0.7)
    #Shivaei+ 15 measure scatter in log(SFR(Halpha))-log(M*) to be 0.3 dex (corrected for uncertainties)
    #We employ cauchy scatter - although I think it should only be scatter in one direction...
    logM = np.log10(mass)
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
        UVJ_data=None,
        extensions=None,
        make_plot=False):

    print "Extracting rows from object ID: ", ID, "(filename: ", fileName, ")"

    # Open the original BEAGLE FITS file
    hdulist = fits.open(fileName)

    # Get the posterior probability
    post = hdulist['posterior pdf'].data['probability']

    # Photometric redshifts computed by Beagle
    redshifts = hdulist['galaxy properties'].data['redshift']

    weights = np.full(1, len(post), dtype=np.float32)

    # Compute the posterior median for the photometric redshift
    CredInterv_redshfit = CredibleInterval(data=redshifts, probability=post)
    median_redshift = CredInterv_redshfit.GetMedian()

    # Check of a UVJ dictionary was passed to the function, in which case you
    # will split among quiescent / SF galaxies based on the UVJ diagram
    galaxy_types = np.zeros(n_samples)
    if UVJ_data is not None:
        UVJ_bands = OrderedDict()
        for key, value in UVJ_data.iteritems():
            UVJ_bands[key] = hdulist[value["extName"]].data[value["colName"]]

        U_V_color = UVJ_bands['U']-UVJ_bands['V']
        V_J_color = UVJ_bands['V']-UVJ_bands['J']

        data = np.zeros((2,len(post)))
        data[0,:] = U_V_color
        data[1,:] = V_J_color

        rows = draw_rows_from_interval(data=data, probability=post, n_draws=n_samples, level=0.68)

        if show_plot:
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)

            plt.scatter(U_V_color,
                    V_J_color,
                    color = 'red',
                    marker = "o",
                    s=24
                    )

            plt.scatter(U_V_color[rows],
                    V_J_color[rows],
                    color = 'blue',
                    marker = "x",
                    s=24
                    )

            plt.show()

        for i, row in enumerate(rows):
            if UVJ_separation(redshifts[row], U_V_color[row], V_J_color[row]) == 1:
                galaxy_types[i] = 1
                

    quiescent_indices = np.where(galaxy_types==1)[0]
    n_quiescent = len(quiescent_indices)

    star_forming_indices = np.where(galaxy_types!=1)[0]
    n_star_forming = len(star_forming_indices)
    
    rows_indices = np.zeros(n_samples, dtype=int)

    if n_quiescent > 0:
        rows_indices[quiescent_indices] = rows[quiescent_indices]

    if n_star_forming > 0:

        # Can pass a "weight_func" that will take as input some parameters of the
        # Beagle output and compute a weight, that will be used to multiply the
        # posterior PDF
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

        #
        if UVJ_data is not None:
            for i, (redshift, UV, VJ) in enumerate(zip(redshifts, U_V_color, V_J_color)):
                if UVJ_separation(redshift, UV, VJ) == 1:
                    mask[i] = False
        
            
        indices = np.arange(len(post))

        # Randomly draw the indices of the rows, using as weights the reweighted posterior PDF 
        wrand = WalkerRandomSampling(reweighted_pdf[mask], keys=indices[mask])
        rows_indices[star_forming_indices] = wrand.random(n_star_forming)

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
        '-nproc',
        help="Number of processors to use",
        action="store", 
        type=int, 
        dest="nproc",
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

    # Arguments used by the "weight function", which is multiplied by the
    # posterior pdf to then select solutions from all the possible ones which
    # are output from Beagle
    weight_func_args = OrderedDict()
    weight_func_args['mass'] = {"colName":"M_star", "extName":"galaxy properties"}
    weight_func_args['SFR'] = {"colName":"SFR", "extName":"star formation"}
    weight_func_args['redshift'] = {"colName":"redshift", "extName":"galaxy properties"}
    weight_func_args['distribution'] = 'gaussian'

    #
    UVJ_data = OrderedDict()
    UVJ_data['U'] = {"colName":"_Bessel_U_ABS", "extName":"absolute magnitudes"}
    UVJ_data['V'] = {"colName":"_Bessel_V_ABS", "extName":"absolute magnitudes"}
    UVJ_data['J'] = {"colName":"_TwoMass_J_ABS", "extName":"absolute magnitudes"}

    # Restrict the allowed solutions to have some parameters within defined ranges
    params_ranges=None
    if args.params_ranges is not None:
        params_ranges = json.loads(args.params_ranges)

    # Which FITS extension (among those available from the Beagle output FITS
    # file) do we print to the output Summary catalogue?
    extensions = ('GALAXY PROPERTIES', 'STAR FORMATION', 'POSTERIOR PDF')

    if args.nproc <= 0:
        results = list()
        for (ID, file) in zip(IDs[0:args.n_objects], file_list[0:args.n_objects]):
            res = draw_rows_from_posterior(ID, file, 
                    n_samples=args.n_samples,
                    params_ranges=params_ranges,
                    weight_func=mass_sfr_test,
                    weight_func_args=weight_func_args,
                    UVJ_data=UVJ_data,
                    extensions=extensions,
                    make_plot=make_plot
                    )
            results.append(res)
    else:

        # Set number of parellel processes to use
        pool = ProcessingPool(nodes=args.nproc)

        # Launch the actual calculation on multiple processesors
        results = pool.map(draw_rows_from_posterior, IDs, file_list, 
                (args.n_samples,)*len(IDs),
                (params_ranges,)*len(IDs),
                (mass_sfr_test,)*len(IDs),
                (weight_func_args,)*len(IDs),
                (UVJ_data,)*len(IDs),
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

