#!/usr/bin/env python

from scipy import spatial
import os
import math
from astropy.coordinates import SkyCoord, Angle
from astropy import units as u
from astropy.io import fits
from astropy.table import Table, Column
import argparse
import ConfigParser
import numpy as np
from collections import OrderedDict
from scipy.interpolate import interp1d
from collections import defaultdict
import itertools
import matplotlib.pyplot as plt
import time
import pickle

#From script written by Jacopo Chevallard
#Adapted by Emma Curtis-Lake to provide the output required for target selection
#from photo-z analysis to NIRCam photometry
#Emma is a bit of a nerd about whitespace (but not always consistent) - in this file
#prefer no tabs, just 4 spaces for indents!

import sys
#sys.path.append(os.path.join(os.environ['PYP_BEAGLE'], "PyP-BEAGLE"))
from pyp_beagle import BeagleDirectories, extract_IDs
#import beagle_multiprocess

from pathos.multiprocessing import ProcessingPool 

deltaEvidence_lim = 6.
deltaZ_lim = 0.5
__logBase10of2 = 3.010299956639811952137388947244930267681898814621085413104274611e-1

results_dir = ""

def ReadInfoFromFitsFiles(ID):
  #Read in all the relevant information from the fits file once to pass to other functions
  #Make into the data format read in usually by astropy.io.fits so that the data is accessed the same way
  full_path = os.path.join(args.results_dir, str(ID)+'_'+suffix)
  data = []
  if not os.path.isfile(full_path):
    return None
  else:
    with fits.open(full_path) as f:
      posteriorCols = []
      for col in f['POSTERIOR PDF'].data.dtype.names:
        posteriorCols.append(fits.Column(name=col,array=f['POSTERIOR PDF'].data[col], format='D'))
      posteriorHdu = fits.BinTableHDU.from_columns(posteriorCols)
      posteriorHdu.header['EXTNAME'] = 'POSTERIOR PDF'
      propertyCols = [fits.Column(name='M_star',array=f['GALAXY PROPERTIES'].data['M_star'],format='D')]
      propertyHdu = fits.BinTableHDU.from_columns(propertyCols)
      propertyHdu.header['EXTNAME'] = 'GALAXY PROPERTIES'
      primaryHdu = fits.PrimaryHDU()
      data = fits.HDUList([primaryHdu,propertyHdu,posteriorHdu])
        
  return data
    

def IntegratedProbAboveZ(ID, zLim, rows=None, plot=False, f=None):

  openFile=False
  if f is None: #read in the information from file
    suffix = BeagleDirectories.suffix + '.fits.gz'

    full_path = os.path.join(args.results_dir, str(ID)+'_'+suffix)
    if not os.path.isfile(full_path):
      return None
    else:
      f = fits.open(full_path)
      openFile = True
  if rows is None:
    rows = np.fromiter((x for x in range(len(f['POSTERIOR PDF'].data['probability']))),np.int)
  probability = f['POSTERIOR PDF'].data['probability'][rows]
  z = f['POSTERIOR PDF'].data['redshift'][rows]
   
  if openFile:
    f.close()
      
  sortIdx = np.argsort(z)
  z_sorted = z[sortIdx]

  #This is edited from the PyP-BEAGLE beagle_summary_catalogue.py
  #get1DInterval function
  # ******************************************************************
  # Here you must simply use `cumsum`, and not `cumtrapz` as in
  # beagle_utils.prepare_violin_plot, since the output of MultiNest are a set
  # of weights (which sum up to 1) associated to each set of parameters (the
  # `p_j` of equation 9 of Feroz+2009), and not a probability density (as the
  # MultiNest README would suggest).
  # ******************************************************************
  cumul_pdf = np.cumsum(probability[sortIdx])
  cumul_pdf /= cumul_pdf[len(cumul_pdf)-1]
  if plot:
    plt.figure()
    plt.plot(z_sorted, cumul_pdf)
    plt.show()
   
  #add to the redshift array if the zLimits chosen are outside the limits
  #in redshift sampled
  if z_sorted[0] > np.min(zLim):
    z_sorted = np.concatenate(([np.min(zLim)],z_sorted),axis=0)
    cumul_pdf = np.concatenate(([0.],cumul_pdf),axis=0)
  if z_sorted[-1] < np.max(zLim):
    z_sorted = np.concatenate((z_sorted,[np.max(zLim)]),axis=0)
    cumul_pdf = np.concatenate((cumul_pdf,[1.]),axis=0)

  # Get the interpolant of the cumulative probability
  f_interp = interp1d(z_sorted, cumul_pdf)
   
  #We output the integrated probability that the object is above a given
  #redshift value
  p1 = f_interp(zLim)
  pOut = 1.-p1

  return pOut

def min_chi2(ID, rows=None, f=None):

  openFile=False
  if f is None:
    #Sorry, not the most efficient to have to open the file again
    suffix = BeagleDirectories.suffix + '.fits.gz'

    full_path = os.path.join(args.results_dir, str(ID)+'_'+suffix)
    if not os.path.isfile(full_path):
      return None
    else:
      f = fits.open(full_path)
      openFile = True

  if rows is None:
    rows = np.fromiter((x for x in range(len(f['POSTERIOR PDF'].data['chi_square']))),np.int)
  c = f['POSTERIOR PDF'].data['chi_square'][rows]

  if openFile:
    f.close()
      
  if (len(c) < 1):
  	return -9999
  else:
  	return np.min(c)

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

def get_mode_rows(ID, results_dir, param_names, param_indices, mode_index, hdulist=None):

  file_name = os.path.join(results_dir + np.str(ID) + "_BEAGLE_MNpost_separate.dat")

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

  openFile=False
  if hdulist is None:
    #Sorry, not the most efficient to have to open the file again
    suffix = BeagleDirectories.suffix + '.fits.gz'

    full_path = os.path.join(args.results_dir, str(ID)+'_'+suffix)
    if not os.path.isfile(full_path):
      return None
    else:
      hdulist = fits.open(full_path)
      openFile=True
           
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

  #Emma comment - I believe this is matching the rows in the
  #separate mode file to entries in the posterior output
  #extension.
  loc = np.where(distances <= max_distance)[0]
  rows = full_indx[loc]

  if openFile:
    hdulist.close()

  return rows
    
def extract_data(ID, n_par, redshift_index, redshift_type=None):

  file_name = os.path.join(results_dir, str(ID) + "_BEAGLE_MNstats.dat")

  if not os.path.isfile(file_name):
    return None

  # This number include the posterior mean, maximum likelihood and
  # maximum a posteriori for each parameter + the headers
  n_lines_per_mode = 8 + n_par*3

  # Useful information for the first mode start at line 11 (in Python
  # we count from 0)
  first_line = 10

  # Now we read the evidence, post mean, maximum likelihood and map for each mode
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

def get1DInterval(ID, param_names, levels=[68., 95.], rows=None, f=None):

  #Check if there are actually any rows in the mode and return -99 if not
  if rows is None:
    rows = np.fromiter((x for x in range(len(f['POSTERIOR PDF'].data['probability']))),np.int)
  if len(rows) > 10:
            
    
    openFile=False
    if f is None:
      #Sorry, not the most efficient to have to open the file again
      suffix = BeagleDirectories.suffix + '.fits.gz'

      full_path = os.path.join(args.results_dir, str(ID)+'_'+suffix)
      if not os.path.isfile(full_path):
        return None,1
      else:
        f = fits.open(full_path)
        openFile=True
    
    param_values = OrderedDict()
        
    probability = f['POSTERIOR PDF'].data['probability'][rows]
    for name in param_names:
      if name == 'Mstar':
        param_values[name] = np.log10(f['GALAXY PROPERTIES'].data['M_star'][rows])
      else:
        param_values[name] = f['POSTERIOR PDF'].data[name][rows]

    output = OrderedDict()
    samplerFlag = 0 #This will be set to 1 if any of the limits are
                    #beyond the probability sampled by multinest
    for key, value in param_values.iteritems():
      print key
      sort_ = np.argsort(value)

      cumul_pdf = np.cumsum(probability[sort_])
      cumul_pdf /= cumul_pdf[len(cumul_pdf)-1]

      # Get the interpolant of the cumulative probability
      f_interp = interp1d(cumul_pdf, value[sort_])

      # You shoud integrate rather than summing here
      mean = np.sum(probability * value) / np.sum(probability)

      median = f_interp(0.5)

      interval = OrderedDict()
      for lev in levels:
        lower_lev = 0.5*(1.-lev/100.)
        upper_lev = 1.-0.5*(1.-lev/100.)
        #check if either limit is lower or higher than current samples
        if lower_lev < cumul_pdf[0]:
          lower_lev = cumul_pdf[0]
          samplerFlag = 1
        if upper_lev > cumul_pdf[-1]:
          upper_lev = cumul_pdf[-1]
          samplerFlag = 1
        low, high = f_interp([lower_lev, upper_lev])
        interval[str(lev)] = np.array([low,high])

      output[key] = {'mean':mean, 'median':median, 'regions':interval}
      if openFile:
        f.close()
        
    
  else:
    output = OrderedDict()
    for name in param_names:
      interval = OrderedDict()
      for lev in levels:
        interval[str(lev)]=np.array([-99,-99])
      output[name] = {'mean':-99, 'median':-99, 'regions':interval}
      samplerFlag = -99
        
        
  return output,samplerFlag
    
def do_all_the_stuff(ID):

  outputDict = {}
  outputPickleFile=pickleFileDir+str(ID)+'.p'
  #Check if the pickle file exists for this object
  if not os.path.isfile(outputPickleFile):

    #read BEAGLE.fits.gz file to pass to other functions here
    fdata = ReadInfoFromFitsFiles(ID)
    print '****', ID, fdata

    #This extracts information in the multinest output file with "_BEAGLE_MNstats.dat"
    d = extract_data(ID,
            n_par=n_par,
            redshift_index=paramDict["redshift"]
            )
    outputDict['data'] = d
    print d
        
    if d is not None:
    
      if args.credible_regions is not None:
        print args.credible_regions
        c,flag = get1DInterval(ID,
                  param_names=param_names,
                  levels=args.credible_regions,
                  f=fdata
                  )
        outputDict['data_cred_region'] = c
        outputDict['cred_region_flag'] = flag

      if args.zLim is not None:
        p = IntegratedProbAboveZ(ID,
                  args.zLim,
                  f=fdata
                  )
        outputDict['data_zLim_probs']=p

      if args.chi2 or args.qFlag:
        chi2Modes = {}
        c = min_chi2(ID,f=fdata)
        chi2Modes["mean"] = (c)
      #And split information by mode
      if args.Mstar:
        mStarModes = {}
        m,flag = get1DInterval(ID,
                               param_names=['Mstar'],
                               levels=args.credible_regions,
                               f=fdata)
        mStarModes["mean"] = m
      if args.qFlag:
        qModes = {}
        p_dz_02 = IntegratedProbAboveZ(ID, outputDict['data_cred_region']["redshift"]["mean"]-0.2, f=f) - \
                  IntegratedProbAboveZ(ID, outputDict['data_cred_region']["redshift"]["mean"]+0.2, f=f)
        temp,flags = get1DInterval(ID, param_names=["redshift"], levels=[99.], f=fdata)
        if flags != -99:
          print temp["redshift"]["regions"]
          z99up_low = temp["redshift"]["regions"][str(99.)][1] - temp["redshift"]["regions"][str(99.)][0]
          qModes["mean"] = (chi2Modes["mean"]/args.nFilt) * z99up_low/p_dz_02
        else:
          qModes["mean"] = -99
            
      for k in d.keys():
        mode_index = k.split('_')[1]
        rows = get_mode_rows(ID, args.results_dir, param_names=["mass","redshift"], param_indices=[1,2], mode_index=mode_index, hdulist=fdata)
        if args.Mstar:
          m,flag = get1DInterval(ID,
                       param_names=['Mstar'],
                       levels=args.credible_regions,
                       rows=rows,
                       f=fdata)
          mStarModes[k] = m
        if args.chi2 or args.qFlag:
          c = min_chi2(ID, rows=rows)
          chi2Modes[k] = c
        if args.qFlag:
          p_dz_02 = IntegratedProbAboveZ(ID, d[k]["posterior_mean"][paramDict["redshift"]-1]-0.2,f=fdata) - \
                    IntegratedProbAboveZ(ID, d[k]["posterior_mean"][paramDict["redshift"]-1]+0.2,f=fdata)
          temp,flags = get1DInterval(ID, param_names=["redshift"], levels=[99.], rows=rows,f=fdata)
          if flags != -99:
            z99up_low = temp["redshift"]["regions"][str(99.)][1] - temp["redshift"]["regions"][str(99.)][0]
            qModes[k] = (chi2Modes[k]/args.nFilt) * z99up_low/p_dz_02
          else:
            qModes[k] = -99
                
      if args.Mstar:
        outputDict['mStar']=mStarModes
      if args.chi2:
        outputDict['chi2']=chi2Modes
      if args.qFlag:
        outputDict['qFlag']=qModes
          
  
    else:
      if args.Mstar:
        outputDict['mStar']=None
      if args.chi2:
        outputDict['chi2']=None
      if args.qFlag:
        outputDict['qFlag']=None
            
    #Write the output for this object to a pickle file
    pickle.dump(outputDict,open(outputPickleFile,'w'))
    return

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
        '--input-coord',
        help="",
        action="store", 
        type=str, 
        nargs=2,
        dest="input_coord",
        default=("deg", "deg")
  )

  parser.add_argument(
        '--beagle-coord',
        help="",
        action="store", 
        type=str, 
        nargs=2,
        dest="beagle_coord",
        default=("deg", "deg")
  )

  parser.add_argument(
        '-r', '--results-dir',
        help="Directory containing BEAGLE results",
        action="store", 
        type=str, 
        dest="results_dir", 
        required=True
  )

  parser.add_argument(
        '-p', '--parameter-file',
        help="Parametr file used in the BEAGLE run",
        action="store", 
        type=str, 
        dest="param_file",
        required=True
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
        '--credible-regions',
        help="Credible regions to calculate",
        action="store", 
        type=float, 
        nargs='+',
        dest="credible_regions"
  )

  parser.add_argument(
        '--zLim-probs',
        help="List of redshifts above which to supply the integrated probability",
        action="store",
        type=float,
        nargs='+',
        dest="zLim"
  )

  parser.add_argument(
        '--chi2',
        help="add average of chi2 of last 10 samples?  Not the best metric but useful for comparing to minimum chi2 based codes (default False)",
        action="store_true",
        default=False,
        dest="chi2"
  )

  parser.add_argument(
        '--m-star',
        help="include Mstar in output",
        action="store_true",
        default=False,
        dest="Mstar"
  )

  parser.add_argument(
        '--q-flag',
        help="add EAZY style q-flag in output - requires credible intervals with a level of 99 (default False)",
        action="store_true",
        default=False,
        dest="qFlag"
  )

  parser.add_argument(
        '--n-filt',
        help="number of filters, required for the reduced chi2 in the quality flag (default 9)",
        action="store",
        type=int,
        default=9,
        dest="nFilt"
  )
    
  parser.add_argument(
        '--interim-results',
        help="Specify a folder name where pickle files will be stored for each object containing all the information derived within this script.  The script will search the folder for the [ID].p file and read that in rather than re-calculate things for this file.",
        action="store",
        type=str,
        dest="interimResults",
        required=True
  )

  # ID list
  parser.add_argument(
        '-idlist','--id_number_list',
        help="List of ID numbers for a subsample?",
        action="store",
        type=str,
        dest="id_number_list",
        required=False
  )


#    if args.qFlag: #We need the upper an lower limits of the 99% credible interval
#        if args.credible_regions is not None:
#            print args.credible_regions
#            #check if level 99 is included and add it if not
#            present = False
#            for i in range(len(args.credible_regions)):
#                if args.credible_regions[i] == 99.:
#                    present = True
#            if not present:
#                args.credible_regions.append(99)
#        else:
#            args.credible_regions = [99.]
    

  # Get parsed arguments
  args = parser.parse_args()

  results_dir = args.results_dir

  # Read the input catalogue
  inputData = Table.read(args.inputCat)

  # Create a subsample given an ID number list
  if (args.id_number_list):
    ID_input_file = np.loadtxt(args.id_number_list)
    if (len(ID_input_file.shape) > 1):
      input_IDs = ID_input_file[:,0].astype(int)
    else:
      input_IDs = ID_input_file.astype(int)
    #ID_numbers = ID_input_file
    n_input = len(input_IDs)
  else:
	# ID in the input catalogue
	input_IDs = extract_IDs(inputData, key='ID')
	n_input = len(inputData.field(0))
	
  # ID in the input catalogue
  #input_IDs = extract_IDs(inputData, key='ID')

  # Read parameter file
  config = ConfigParser.SafeConfigParser()

  # Search for the parameter file in the results directory
  if os.path.isabs(args.param_file):
    param_file = args.param_file
  else:
    param_file = os.path.join(args.results_dir, 'BEAGLE-input-files', args.param_file)

  config.read(param_file)

  # Columns to be added to the catalogues
  param_names = ["redshift", "mass"]
  
  pickleFileDir = args.interimResults
  # check if folder exists and make it if not
  if not os.path.exists(args.interimResults):
    os.system("mkdir "+args.interimResults)

  dictKeys = OrderedDict()

#---------------------------------------------------------------------------
#-------- Setting up the column names and types --------------
#---------------------------------------------------------------------------

  dictKeys["ID"] = {"type":"S15", "format":"s"}
  #dictKeys["distance"] = {"type":np.float32, "format":".3f"}

  for name in param_names:
    key = name+"_beagle_mean"
    dictKeys[key] = {"type":np.float32, "format":".3f"}
    for j in range(2):
      suff = str(j+1)
      key = name+"_beagle_"+suff
      dictKeys[key] = {"type":np.float32, "format":".3f"}

      key = name+"_beagle_err_"+suff
      dictKeys[key] = {"type":np.float32, "format":".3f"}

    if args.credible_regions is not None:
      for region in args.credible_regions:
        key = name + "_" + str(np.rint(region)) + "_low"
        dictKeys[key] = {"type":np.float32, "format":".3f"}
        key = name + "_" + str(np.rint(region)) + "_up"
        dictKeys[key] = {"type":np.float32, "format":".3f"}

  if args.zLim is not None:
    for zLim in args.zLim:
      key = "redshift_p_gt_" + str(zLim)
      dictKeys[key] = {"type":np.float32, "format":".3f"}
  
  if args.chi2 is not None:
    key = "chi2"
    dictKeys[key] = {"type":np.float32, "format":".3f"}
            
  if args.Mstar:
    # We add the modal mean Mstar values as well as the total average and credible intervals
    name = "Mstar"
    for j in range(2):
      suff = str(j+1)
      key = name+"_beagle_"+suff
      dictKeys[key] = {"type":np.float32, "format":".3f"}
    if args.credible_regions is not None:
      key = name + "_beagle_mean"
      dictKeys[key] = {"type":np.float32, "format":".3f"}
      for region in args.credible_regions:
        key = name + "_" + str(np.rint(region)) + "_low"
        dictKeys[key] = {"type":np.float32, "format":".3f"}
        key = name + "_" + str(np.rint(region)) + "_up"
        dictKeys[key] = {"type":np.float32, "format":".3f"}
                
  if args.qFlag:
    #We add a quality flag for the main average redshift as well as for each mode
    dictKeys["q_beagle"] = {"type":np.float32, "format":".3f"}
    for j in range(2):
      suff = str(j+1)
      dictKeys["q_beagle_"+suff] = {"type":np.float32, "format":".3f"}
#
#      key = name+"_beagle_err_"+suff
#      dictKeys[key] = {"type":np.float32, "format":".3f"}
            
            
#      key = "chi2"
#      dictKeys[key] = {"type":np.float32, "format":".3f"}
#      dictKeys["q_beagle"] = {"type":np.float32, "format":".3f"}
#      for j in range(2):
#        suff = str(j+1)
#        key = "q_beagle_"+suff
#        dictKeys[key] = {"type":np.float32, "format":".3f"}

  paramDict = OrderedDict()
  # Determine number of free parameters by counting columns in Beagle output file
  suffix = BeagleDirectories.suffix + '.fits.gz'
  for file in sorted(os.listdir(args.results_dir)):
    full_path = os.path.join(args.results_dir, file)
    if file.endswith(suffix) and os.path.getsize(full_path) > 0:
      with fits.open(full_path) as f:
        n_par = len(f["POSTERIOR PDF"].data.dtype.names)-4
        #print n_par, f['POSTERIOR PDF'].data.dtype.names
                
        for name in param_names:
          for i, col_name in enumerate(f["POSTERIOR PDF"].data.dtype.names):
            if name == col_name:
              paramDict[name] = (i+1)-4
              break


      break

  dictKeys["deltaEvidence"] = {"type":np.float32, "format":".2f"}
  dictKeys["KF_flag"] = {"type":np.int, "format":"1d"}
  dictKeys["P1/P2"] = {"type":np.float32, "format":".3e"}
  dictKeys["cred_region_flag"] = {"type":np.int, "format":"1d"}

#---------------------------------------------------------------------------
# setting out the columns that data will be added to once calculated.
#---------------------------------------------------------------------------

  newCols = OrderedDict()

  for key, value in dictKeys.iteritems():
    Type = value["type"]
    if isinstance(Type, str):
      newCols[key] = np.full(n_input, "-99", Type)
    else:
      newCols[key] = np.full(n_input, -99, Type)

#---------------------------------------------------------------------------
# Now performing the calculations
#---------------------------------------------------------------------------

  input_idx = range(n_input)

  # Put oriignal (input catalgoue) IDs in the output catalogue
  newCols["ID"] = np.array(input_IDs)

  #print "-------> ", get1DInterval(Beagle_IDs[match_ok[0]], param_names=param_names, levels=[68., 95., 99.7])

  #pause
  # If the user does not specify the number of processors to be used, assume that it is a serial job
  if args.credible_regions is not None:
    data_cred_region = list()
    cred_region_flag = list()
  if args.zLim is not None:
    data_zLim_probs = list()
  if args.Mstar is not None:
    mStar = list()
  if args.qFlag is not None:
    p_dz_02 = list()
    qFlag = list()

  data = list()
  zProb = list()
  chi2 = list()
  t = list()
  if args.nproc <= 0:

    for indx in range(len(input_IDs)):
        
      t0 = time.clock()
      ID = input_IDs[indx]
      #Perform all the calculations - results get saved to a pickle file
      do_all_the_stuff(ID)
      t.append(time.clock()-t0)

        
  # Otherwise you use pathos to run in parallel on multiple CPUs
  else:
    # Set number of parellel processes to use
    pool = ProcessingPool(nodes=args.nproc)

    # Launch the actual calculation on multiple processesors
    data = pool.map(do_all_the_stuff,
                    input_IDs)

#-----------------------------------------------------------------------------
# If you've got to here, nothing broke so time to make the output file!
#-----------------------------------------------------------------------------
      
#-----------------------------------------------------------------------------
# All the data is put into the output format to be written to a fits file
#-----------------------------------------------------------------------------
  for i, indx in enumerate(input_idx):


    ID = input_IDs[indx]
    outputPickleFile=pickleFileDir+str(ID)+'.p'
    #Check if the pickle file exists for this object
    if os.path.isfile(outputPickleFile):
      outputDict = pickle.load(open(outputPickleFile,'r'))
      d = outputDict['data']
      
      if d is None:
        continue

      for j, (key, value) in enumerate(d.iteritems()):
        suff = str(j+1)
        for name, row_index in paramDict.iteritems():
          newCols[name+"_beagle_"+suff][indx] = value["posterior_mean"][row_index-1]
          newCols[name+"_beagle_err_"+suff][indx] = value["posterior_sigma"][row_index-1]
            
        if args.Mstar:
          newCols["Mstar_beagle_"+suff][indx] = outputDict["mStar"][key]["Mstar"]["mean"]
                
        if args.qFlag:
          print "q_beagle_"+suff, indx, key, outputDict["qFlag"]
          newCols["q_beagle_"+suff][indx] = outputDict["qFlag"][key]

        if j == 0:
          deltaEvidence = value["evidence"]
        else:
          deltaEvidence -= value["evidence"]

      if j > 0:
        newCols["deltaEvidence"][indx] = deltaEvidence
        newCols["P1/P2"][indx] = np.exp(-deltaEvidence)

        if 2.*deltaEvidence < 2.:
          newCols["KF_flag"][indx] = 0
        elif 2.*deltaEvidence < 6.:
          newCols["KF_flag"][indx] = 1
        elif 2.*deltaEvidence < 10.:
          newCols["KF_flag"][indx] = 2
        else:
          newCols["KF_flag"][indx] = 3
        
      if args.credible_regions is not None:
        c = outputDict["data_cred_region"]
        flag = outputDict["cred_region_flag"]
        if c is not None:
          for name in param_names:
            key = name+"_beagle_mean"
            newCols[key][indx] = c[name]["mean"]
            for region in args.credible_regions:
              key = name + "_" + str(np.rint(region)) + "_low"
              print name, c[name]["regions"]
              newCols[key][indx] = c[name]["regions"][str(region)][0]
              key = name + "_" + str(np.rint(region)) + "_up"
              newCols[key][indx] = c[name]["regions"][str(region)][1]
          newCols['cred_region_flag'][indx] = flag
      if args.Mstar:
        m = outputDict["mStar"]
        key = "Mstar_beagle_mean"
        newCols[key][indx] = m["mean"]["Mstar"]["mean"]
        for region in args.credible_regions:
          key = "Mstar_"+str(np.rint(region))+"_low"
          newCols[key][indx] = m["mean"]["Mstar"]["regions"][str(region)][0]
          key = "Mstar_"+str(np.rint(region)) + "_up"
          newCols[key][indx] = m["mean"]["Mstar"]["regions"][str(region)][1]
                
      if args.zLim is not None:
        p = outputDict["data_zLim_probs"]
        if p is not None:
          for j,zLim in enumerate(args.zLim):
            key = "redshift_p_gt_" + str(zLim)
            newCols[key][indx] = p[j]

      if args.chi2:
        c = outputDict["chi2"]
        print indx, c
        newCols["chi2"][indx] = c["mean"]
            
      if args.qFlag:
        newCols["q_beagle"][indx] = outputDict["qFlag"]["mean"]
            
            

  myCols = list()
  for i, (key, col) in enumerate(newCols.iteritems()):
    tmpCol = Column(col, name=key, dtype=dictKeys[key]["type"], format='%'+dictKeys[key]['format'])
    myCols.append(tmpCol)

  newTable = Table(myCols)

  file_name = os.path.splitext(args.inputCat)[0] + '_Beagle_VAC.txt'
  print "Output (ASCII) file_name: ", file_name
  newTable.write(file_name, format="ascii.commented_header")

  file_name = os.path.splitext(args.inputCat)[0] + '_Beagle_VAC.fits'
  print "Output (FITS) file_name: ", file_name
  newTable.write(file_name, format="fits", overwrite=True)


