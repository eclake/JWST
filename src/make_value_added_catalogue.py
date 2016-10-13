from scipy import spatial
import os
import math
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.table import Table, Column
import argparse
import ConfigParser
import numpy as np
from collections import OrderedDict
from collections import defaultdict
import itertools
import matplotlib.pyplot as plt

deltaEvidence_lim = 6.
deltaZ_lim = 0.5
__logBase10of2 = 3.010299956639811952137388947244930267681898814621085413104274611e-1


def RoundToSigFigs( x, sigfigs ):
    """
    Rounds the value(s) in x to the number of significant figures in sigfigs.

    Restrictions:
    sigfigs must be an integer type and store a positive value.
    x must be a real value or an array like object containing only real values.
    """
    if not ( type(sigfigs) is int or np.issubdtype(sigfigs, np.integer)):
        raise TypeError( "RoundToSigFigs: sigfigs must be an integer." )

    if not np.all(np.isreal( x )):
        raise TypeError( "RoundToSigFigs: all x must be real." )

    if sigfigs <= 0:
        raise ValueError( "RoundtoSigFigs: sigfigs must be positive." )

    mantissas, binaryExponents = np.frexp( x )

    decimalExponents = __logBase10of2 * binaryExponents
    intParts = np.floor(decimalExponents)

    mantissas *= 10.0**(decimalExponents - intParts)

    return np.around( mantissas, decimals=sigfigs - 1 ) * 10.0**intParts

def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        dist = np.abs(value-array[idx-1])
        return idx-1, dist
    else:
        dist = np.abs(value-array[idx])
        return idx, dist

def get_mode_rows(ID, results_dir, param_names, param_indices, mode_index):

    file_name = os.path.join(results_dir, ID + "_BEAGLE_MNpost_separate.dat")

    n_empty = 0
    prev_is_empty = False
    values = defaultdict(list)
    with open(file_name , 'r') as f:

        for line in f:

            is_empty = False
            if not line.strip():
                is_empty = True

            if is_empty and not prev_is_empty:
                n_empty += 1

            prev_is_empty = is_empty

            if not is_empty and n_empty == (int(mode_index)):
                split = line.split()
                for name, indx in zip(param_names, param_indices):
                    value = split[1+indx]
                    values[name].append(float(value))

    n = len(values[param_names[0]])
    n_par = len(param_names)
    np_values = np.zeros([n, n_par])
    for i, (key, value) in enumerate(values.iteritems()):
        np_values[:,i] = value

    tree = spatial.cKDTree(np_values)


    file_name = os.path.join(results_dir, ID + "_BEAGLE.fits.gz")

    hdulist = fits.open(file_name)
    n = len(hdulist['POSTERIOR PDF'].data['probability'])
    full_values = np.zeros([n, n_par])
    for i, name in enumerate(param_names):
        full_values[:,i] = hdulist['POSTERIOR PDF'].data[name]

    #sor = np.argsort(full_values)

    full_indx = np.arange(len(full_values))
    #full_values_sor = full_values[sor]
    #full_indx = full_indx[sor]

    rows = list()

    distances, indiced = tree.query(full_values)

    max_distance = 1.E-06
    
    #print "distance: ", distances
    #plt.hist(np.log10(distances), 50, normed=1, facecolor='green', alpha=0.75)
    #plt.show()
    #pause    

    loc = np.where(distances <= max_distance)[0]
    rows = full_indx[loc]

    hdulist.close()

    return rows
    
def extract_data(ID, results_dir, n_par, redshift_index, redshift_type=None):

    # This number include the posterior mean, maximum likelihood and
    # maximum a posteriori for each parameter + the headers
    n_lines_per_mode = 8 + n_par*3

    # Useful information for the first mode start at line 11 (in Python
    # we count from 0)
    first_line = 10

    # Now we read the evidence, post mean, maximum likelihood and map for each mode
    file_name = os.path.join(results_dir, ID + "_BEAGLE_MNstats.dat")
    f = open(file_name , 'r')
    outData = OrderedDict()
    post_sig = list()
    post_mean = list()
    max_likelihood = list()
    max_a_post = list()
    logEvidence = list()
    mode_count = 1

    for i, line in enumerate(f):
        
        # Evidence
        if i == first_line:
            logEv = float(line.split()[2])
        # Posterior mean for each parameter 
        elif((i >= first_line+3) and (i < first_line+3+n_par)):
            post_mean.append(float(line.split()[1])) 
            post_sig.append(float(line.split()[2])) 
        # Maximum likelihood for each parameter 
        elif((i >= first_line+6+n_par) and (i < first_line+6+2*n_par)):
            max_likelihood.append(float(line.split()[1])) 
        # Maximum a posteriori for each parameter 
        elif((i >= first_line+9+2*n_par) and (i < first_line+9+3*n_par)):
            max_a_post.append(float(line.split()[1])) 

        # Once you've read the data for the first mode, put them into
        # the MultiNestObject!
        if i == (first_line + n_lines_per_mode):
            key = "mode_"+str(mode_count)
            outData[key] = {"evidence":logEv, 
                    "posterior_mean":post_mean, 
                    "posterior_sigma":post_sig, 
                    "max_likelihood":max_likelihood,
                    "max_a_post":max_a_post}

            post_mean = list()
            post_sig = list()
            max_likelihood = list()
            max_a_post = list()
            logEvidence.append(logEv)
            first_line += n_lines_per_mode + 5
            mode_count = mode_count + 1

    f.close()
    sor = np.arange(len(logEvidence))
    data = OrderedDict()
    if len(logEvidence) > 1:
        sor = np.argsort(logEvidence)
        key_1 = outData.keys()[sor[-1]]
        key_2 = outData.keys()[sor[-2]]

        deltaEvidence = outData[key_1]["evidence"] - outData[key_2]["evidence"]
        deltaZ = outData[key_1]["posterior_mean"][redshift_index-1] - outData[key_2]["posterior_mean"][redshift_index-1]
        #print 'deltaZ: ', deltaZ, outData[key_1]["posterior_mean"][1], outData[key_2]["posterior_mean"][1]

        if deltaEvidence < deltaEvidence_lim and abs(deltaZ) > deltaZ_lim:
            if redshift_type is not None:
                z_1 = outData[key_1]["posterior_mean"][redshift_index-1]
                z_2 = outData[key_2]["posterior_mean"][redshift_index-1]
                if z_1 > z_2:
                    key_max = key_1
                    key_min = key_2
                else:
                    key_max = key_2
                    key_min = key_1
    
                if redshift_type == "high":
                    data[key_max] = outData[key_max]
                    data[key_min] = outData[key_min]
                    return data
                elif  redshift_type == "low":
                    data[key_min] = outData[key_min]
                    data[key_max] = outData[key_max]
                    return data

            data[key_1] = outData[key_1]
            data[key_2] = outData[key_2]
            return data

        data[key_1] = outData[key_1]
        return data

    return outData

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-i', '--input',
        help="Name of the input catalogue", 
        action="store", 
        type=str, 
        dest="inputCat", 
        required=True
    )

    parser.add_argument(
        '-r', '--results-directory',
        help="Name of the direvotry containing Beagle results",
        action="store", 
        type=str, 
        dest="resultsDir", 
        required=True
    )

    # Get parsed arguments
    args = parser.parse_args()

    # Read the input catalogue
    inputData = Table.read(args.inputCat, format='ascii')

    # Read parameter file
    config = ConfigParser.SafeConfigParser()

    # Search for the parameter file in the results directory
    for root, dirs, files in os.walk(args.resultsDir):
        for file in files:
            if file.endswith('.param'):
                param_file = os.path.join(root, file)
                break

    config.read(param_file)

    file_name = os.path.expandvars(config.get('main', 'PHOTOMETRIC CATALOGUE'))
    BeagleData = fits.open(file_name)[1].data

    inputCoord = SkyCoord(ra=inputData['RA']*u.degree, dec=inputData['DEC']*u.degree)  
    BeagleCoord = SkyCoord(ra=BeagleData['RA']*u.degree, dec=BeagleData['DEC']*u.degree)  
    idx, d2d, d3d = inputCoord.match_to_catalog_sky(BeagleCoord)  

    n = len(inputData.field(0))
    # Columns to be added to the catalogues
    paramKeys = ["z_beagle", "z_beagle_err", "M", "M_err", "sSFR", "sSFR_err"]
    dictKeys = list()
    dictKeys.append("UVUDF_ID")
    dictKeys.append("distance")
    for i in (1,2):
        for key in paramKeys:
            dictKeys.append(key+'_'+str(i))

    dictKeys.append("deltaEvidence")
    dictKeys.append("KF_flag")
    dictKeys.append("P1/P2")

    dictTypes = [
            np.int, 
            np.float32,
            np.float32, np.float32, np.float32, np.float32, np.float32, np.float32,
            np.float32, np.float32, np.float32, np.float32, np.float32, np.float32,
            np.float32, np.int, np.float32
            ]

    colFormat = ["4d",
            ".3f",
            ".3f",".3f",".2f",".2f",".2f",".2f",
            ".3f",".3f",".2f",".2f",".2f",".2f",
            ".2f",
            "1d",
            ".3E"
            ]

    newCols = OrderedDict()

    for key, Type in zip(dictKeys, dictTypes):
        newCols[key] = np.full(n, -99, Type)
    
    zCol, mCol, sCol = 2, 1, 6

    for i, (indx, dist) in enumerate(zip(idx, d2d.arcsecond)):
        #print '\n i: ', i
        if dist > -2.:
            ID = BeagleData['ID'][indx]
            newCols["UVUDF_ID"][i] = int(ID)
            newCols["distance"][i] = dist
            #print "ID: ", ID
            data = extract_data(str(ID), n_par=7)

            for j, (k, mod) in enumerate(data.iteritems()):
                newCols["z_beagle_"+str(j+1)][i] = 1.
                newCols["z_beagle_"+str(j+1)][i] = mod["posterior_mean"][zCol-1]
                newCols["z_beagle_err_"+str(j+1)][i] = mod["posterior_sigma"][zCol-1]

                newCols["M_"+str(j+1)][i] = mod["posterior_mean"][mCol-1]
                newCols["M_err_"+str(j+1)][i] = mod["posterior_sigma"][mCol-1]

                newCols["sSFR_"+str(j+1)][i] = mod["posterior_mean"][sCol-1]
                newCols["sSFR_err_"+str(j+1)][i] = mod["posterior_sigma"][sCol-1]

                if j == 0:
                    deltaEvidence = mod["evidence"]
                else:
                    deltaEvidence -= mod["evidence"]

            if j > 0:
                newCols["deltaEvidence"][i] = deltaEvidence
                newCols["P1/P2"][i] = np.exp(-deltaEvidence)

                if 2.*deltaEvidence < 2.:
                    newCols["KF_flag"][i] = 0
                elif 2.*deltaEvidence < 6.:
                    newCols["KF_flag"][i] = 1
                elif 2.*deltaEvidence < 10.:
                    newCols["KF_flag"][i] = 2
                else:
                    newCols["KF_flag"][i] = 3
    
    myCols = list()
    for i, (key, col) in enumerate(newCols.iteritems()):
        tmpCol = Column(col, name=key, dtype=dictTypes[i], format='%'+colFormat[i])
        myCols.append(tmpCol)

    newTable = Table(myCols)

    file_name = os.path.basename(args.inputCat).split('.')[0] + '_Beagle.txt'
    print "file_name: ", file_name
    newTable.write(file_name, format="ascii.commented_header")
    
    #print "idx: ", idx, d2d.arcsecond
    #print "---> ", inputData.field(0)
    #print "---> ", BeagleData.field(0)[idx]
