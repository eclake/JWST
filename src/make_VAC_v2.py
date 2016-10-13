from collections import OrderedDict
from scipy.interpolate import interp1d
from astropy.io import fits
import os
import numpy as np
from astropy.table import Table, Column
from make_value_added_catalogue import extract_data, get_mode_rows

def get1DInterval(param_values, probability, levels):

    """ 
    Compute several quantities from a 1D probability density function

    Parameters
    ----------
    param_values : numpy array
        Contains the values of the parameter.

    probability : numpy array 
        Contains the probability associated with each value of the
        parameter, hence must have same dimension as 'param_values'.

    levels : numpy array or list containing float
        Contains the (percentage) levels used to compute the credible
        regions, e.g. levels=[68.,95.] will compute 68 % and 95 % (central)
        credible regions

    Returns
    -------
    mean : float
        Mean of the parameter, computed as
        sum(probability*param_values)/sum(probability)
    median : float
        Median of the parameter, computed from the cumulative integral of
        the PDF
    interval : list of float
        2-dimensional list containing the lower and upper value of the
        parameter corresponding to the different `levels`

    """

    sort_ = np.argsort(param_values)
    cumul_pdf = np.cumsum(probability[sort_])
    cumul_pdf /= cumul_pdf[len(cumul_pdf)-1]

    # Get the interpolant of the cumulative probability
    f_interp = interp1d(cumul_pdf, param_values[sort_])

    #
    mean = np.sum(probability * param_values) / np.sum(probability)

    median = f_interp(0.5)

    interval = list()
    for lev in levels:

        low, high = f_interp([0.5*(1.-lev/100.), 1.-0.5*(1.-lev/100.)])
        interval.append([low,high])

    return mean, median, interval

# Catalogue containing the UVUDF objects in Andy's catalogue
#fileName = "/Users/jchevall/JWST/UDF_value_added_catalogues/matchBouwens_UVUDF.txt"
fileName = "/home/jchevall/JWST/UDF_value_added_catalogues/matchBouwens_UVUDF.txt"
obs_cat = Table.read(fileName, format='ascii')
IDs = obs_cat['UVUDFname']

# Directory containing Beagle results
#results_dir = "/Users/jchevall/Coding/BEAGLE/files/results/UVUDF/Emma_constant_SFR"
results_dir = "/user_data/jchevall/Emma/BEAGLE_results/UVUDF_constant_SFR"

data_lines = OrderedDict()

line = {"center":6563., "left_cont":[6552., 6557.], "right_cont":[6569., 6574.]}
#line = {"center":6563., "left_cont":[6545., 6550.], "right_cont":[6576., 6581.]}
data_lines['Halpha'] = line

line = {"center":4861., "left_cont":[4850., 4855.], "right_cont":[4866., 4871.]}
#line = {"center":4861., "left_cont":[4843., 4848.], "right_cont":[4874., 4879.]}
data_lines['Hbeta'] = line

width = 3.

n = len(IDs)

dictKeys, dictTypes, colFormat = list(), list(), list()

dictKeys.append("UVUDF_ID")
dictTypes.append(np.int)
colFormat.append("4d")

for key, value in data_lines.iteritems():
    dictKeys.append(key+'_flux')
    dictTypes.append(np.float32)
    colFormat.append(".3E")

    dictKeys.append(key+'_flux_err')
    dictTypes.append(np.float32)
    colFormat.append(".3E")

    dictKeys.append(key+'_EW')
    dictTypes.append(np.float32)
    colFormat.append(".3f")

    dictKeys.append(key+'_EW_err')
    dictTypes.append(np.float32)
    colFormat.append(".3f")

galProp = OrderedDict()
galProp["redshift"] = {"extension":'GALAXY PROPERTIES', 'format':".3f"}
galProp["M_star"] = {"extension":'GALAXY PROPERTIES', 'format':".3E"}
galProp["M_tot"] = {"extension":'GALAXY PROPERTIES', 'format':".3E"}
galProp["SFR"] = {"extension":'STAR FORMATION', 'format':".3E"}
galProp["sSFR"] = {"extension":'STAR FORMATION', 'format':".3E"}

for key, value in galProp.iteritems():
    dictKeys.append(key)
    dictTypes.append(np.float32)
    colFormat.append(value['format'])

    dictKeys.append(key+'_err')
    dictTypes.append(np.float32)
    colFormat.append(value['format'])

newCols = OrderedDict()

redshift_index = 2

for key, Type in zip(dictKeys, dictTypes):
    newCols[key] = np.full(n, -99, Type)

for i, ID in enumerate(IDs):

    newCols['UVUDF_ID'][i] = ID
    ID = str(ID)
    print "ID: ", ID

    fits_fileName = ID + "_BEAGLE.fits.gz"
    fits_fileName = os.path.join(results_dir, fits_fileName)

    if not os.path.isfile(fits_fileName) or os.path.getsize(fits_fileName) == 0:
        continue

    data = extract_data(ID, results_dir, n_par=6, redshift_index=redshift_index, redshift_type='high')

    for k in data:
        key = k
        break

    mode_index = key.split('_')[1]
    #rows = get_mode_rows(ID, results_dir, param_names=["mass", "redshift", "tauV_eff", "metallicity", "nebular_logU"], param_indices=[1, 2, 3, 4, 6], mode_index=mode_index)
    rows = get_mode_rows(ID, results_dir, param_names=["mass", "redshift"], param_indices=[1, 2], mode_index=mode_index)

    #print "rows: ", rows[0:-1]

    hdulist = fits.open(fits_fileName) 
    SED = hdulist['FULL SED'].data[rows,:]
    wl = np.ravel(hdulist['FULL SED WL'].data[0][:])
    probability = np.ravel(hdulist['POSTERIOR PDF'].data['probability'][rows])

    # Compute Ha and Hb integrated fluxes and EW
    for key, value in data_lines.iteritems():

        # Compute the flux integrated around the line center, for all selected rows
        i0 = np.searchsorted(wl, value["center"]-width)
        i1 = np.searchsorted(wl, value["center"]+width)
    
        # NB: Note that to get the integrtaed flux in units of erg s^-1 you need to multiply the flux below by delta lambda (wl[i1]-wl[i0])
        print "wl: ", wl[i0], wl[i1]
        flux = np.ravel(np.trapz(SED[:,i0:i1+1], x=wl[i0:i1+1], axis=1)) / (wl[i1]-wl[i0])
        
        # Get the posterior mean, median, and 68% credible interval
        mean, median, interval = get1DInterval(flux, probability, levels=[68.])

        # Put the data into the corresponding columns
        newCols[key+'_flux'][i] = median
        newCols[key+'_flux_err'][i] = 0.5*(interval[0][1]-interval[0][0])
    
        # To compute the EW, you firstly computed the integrated flux in a window on the left of the EL
        il0 = np.searchsorted(wl, value["left_cont"][0])
        il1 = np.searchsorted(wl, value["left_cont"][1])
        if il0 == il1:
            il0 -= 1
        flux_left = np.ravel(np.trapz(SED[:,il0:il1+1], x=wl[il0:il1+1], axis=1)) / (wl[il1]-wl[il0])

        # Repeat the same on the right of the line
        ir0 = np.searchsorted(wl, value["right_cont"][0])
        ir1 = np.searchsorted(wl, value["right_cont"][1])
        if ir0 == ir1:
            ir1 += 1
        flux_right = np.ravel(np.trapz(SED[:,ir0:ir1+1], x=wl[ir0:ir1+1], axis=1)) / (wl[ir1]-wl[ir0])

        # Standard definition of EW = int 1-f_lambda/f0 dlambda
        #f0 = 1.-SED[:,i0:i1+1]/(0.5*(flux_left+flux_right))
        #flux = np.ravel(np.trapz(f0, x=wl[i0:i1+1], axis=1)) / (wl[i1]-wl[i0])

        # Approximate definition
        # NB: Note that to get the EW in the correct units of Ang^-1 you need to multiply the flux below by delta lambda (wl[i1]-wl[i0])
        f0 = SED[:,i0:i1+1]
        flux = np.ravel(np.trapz(f0, x=wl[i0:i1+1], axis=1)) / (wl[i1]-wl[i0])
        flux /= 0.5*(flux_left+flux_right)
        flux -= 1.
        flux = -flux

        # Get the posterior mean, median, and 68% credible interval
        mean, median, interval = get1DInterval(flux, probability, levels=[68.])

        # Put the data into the corresponding columns
        newCols[key+'_EW'][i] = median
        newCols[key+'_EW_err'][i] = 0.5*(interval[0][1]-interval[0][0])

        #print "---> ", newCols[key+'_flux'][i]
        #print "---> ", newCols[key+'_flux_err'][i]

        #print "---> ", newCols[key+'_EW'][i]
        #print "---> ", newCols[key+'_EW_err'][i]


    for key, value in galProp.iteritems():

        x = hdulist[value["extension"]].data[key][rows]

        # Get the posterior mean, median, and 68% credible interval
        mean, median, interval = get1DInterval(x, probability, levels=[68.])
    
        # Put the data into the corresponding columns
        newCols[key][i] = median
        newCols[key+'_err'][i] = 0.5*(interval[0][1]-interval[0][0])
        

    hdulist.close()

    pause

myCols = list()
for i, (key, col) in enumerate(newCols.iteritems()):
    tmpCol = Column(col, name=key, dtype=dictTypes[i], format='%'+colFormat[i])
    myCols.append(tmpCol)

newTable = Table(myCols)

file_name = os.path.basename(fileName).split('.')[0] + "_Beagle_constant_SFR.txt"
newTable.write(file_name, format="ascii.commented_header")
