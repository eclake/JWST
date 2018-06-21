#########################################################################
# IDENT        p_spectrumMOS1x3.py
# LANGUAGE     Python
# AUTHOR       P.FERRUIT          
# PURPOSE      Computation of a simulated spectrum for the observation
#              of a source with JWST/NIRSpec in its MOS mode (MOS1x3).
#
# VERSION
#	1.0.0 PF 26.05.2015 Creation
#	1.0.1 PF 04.06.2015 Update
#		1) Eliminated the  spectral element size parameter that did not
#		make sense because we are working at pixel level. The default
#		number was 1 so the change should be transparent if the
#		user did not set the parameter.
#	1.0.2 PF 12.04.2016 Update
#		1) Added a column tracking the absolute noise (sigma) level
#		in the same units than the spectrum.
#		2) Now using the absolute noise levels to generate the noisy
#		spectrum. This way, we do not have incorrect noisy spectrum
#		values (= zero) in parts of thespectra where the objects
#		intensity is zero.
#
version = '1.0.2'
#########################################################################
import os
import math
import sys
import datetime
import numpy
from astropy.io import fits as pyfits 
import pylab
import argparse
from JWSTpylib import c_straylightOTE
from JWSTpylib import c_zodiacalLight as zod
from JWSTpylib import c_throughputFile
from JWSTpylib import c_transmissionFunction
from JWSTpylib import c_sourceFile
from JWSTpylib import c_sourceFunction
from JWSTpylib import c_spectrum
from JWSTpylib import f_interpolation
from JWSTpylib.sensitivity import c_detectorCosmetics as c_detectorCosmetics
from JWSTpylib.sensitivity import c_summationBoxLosses as c_summationBoxLosses
from JWSTpylib.sensitivity import c_pce as c_pce
from JWSTpylib.sensitivity import c_tableSNR as c_tableSNR

# =======================================================================
# Addition from Jacopo Chevallard - 26/01/2017
from compute_MSA_slit_throughput import MSAThroughput 
# =======================================================================

# Speed of light in m s-1
c = 2.997925e8
# Planck constant in J s
h = 6.62620e-34
# Boltzmann constant in J K-1
k = 1.38062e-23


# =======================================================================
# Declaring and loading the input arguments
# =======================================================================
parser = argparse.ArgumentParser()
parser.add_argument('filename', type=str, help='Input spectrum filename [PS: in Jy; ES in Jy arcsec-2].')
parser.add_argument('data', type=str, help='Path to data folder from the JWSTpytools distribution.')
parser.add_argument('pcefolder', type=str, help='Name of the folder containing the PCE files.')
parser.add_argument('FWA', type=str, help='FWA configuration.')
parser.add_argument('GWA', type=str, help='GWA configuration.')
listOfValidSpatialTypes = ['PS', 'ES']
parser.add_argument('spatype', type=str, help='Spatial type of the source (PS = point source; ES = extended [uniform] source).', choices=listOfValidSpatialTypes)
parser.add_argument('nexp', type=int, help='Number of MULTIACCUM22x4 exposures.')
parser.add_argument('out', type=str, help='Path where the output files and figures will be generated.')
parser.add_argument('prefix', type=str, help='Prefix used to generate the names of the output files and figures.')

parser.add_argument('-tn', '--totalnoise', type=float, help='Detector total noise level for a MULTIACCUM22x4 exposure (1 sigma, in electrons).', default=7.0)
parser.add_argument('-dc', '--darkcurrent', type=float, help='Average dark current rate per pixel (in electrons/s).', default=0.01)
parser.add_argument('-bpf', '--badpixfrac', type=float, help='Fraction of bad pixels (unitless).', default=0.035)
parser.add_argument('-opf', '--openpixfrac', type=float, help='Fraction of open pixels (unitless).', default=0.0000055)
parser.add_argument('-nb', '--nb', type=int, help='Number of background elements.', default=2)
parser.add_argument('-sy', '--sumy', type=int, help='Projected size of the summation box element along the spatial direction (integer, in pixels).', default=4)
parser.add_argument('-rx', '--rejx', type=int, help='Size of the rejection box along the spectral direction for estimating the impact of bad pixels (integer, in pixels).', default=2)
parser.add_argument('-ry', '--rejy', type=int, help='Size of the rejection box along the spatial direction for estimating the impact of bad pixels (integer, in pixels).', default=3)
parser.add_argument('-gd', '--gammadark', type=float, help='Accuracy of the dark current subtraction.', default=0.07)
parser.add_argument('-gff', '--gammaflatfield', type=float, help='Accuracy of the flat-field correction.', default=0.02)

# =======================================================================
# Addition from Jacopo Chevallard - 26/01/2017
# =======================================================================
parser.add_argument('--sersic', type=float, dest='sersic', help='Sersic index of the source, used to calculate slit losses.')
parser.add_argument('--effective-radius', type=float, dest='effective_radius', help='Effective radius (in arcsec) of the source, used to calculate slit losses.')

# =======================================================================
# Addition from Jacopo Chevallard - 21/06/2018
# =======================================================================
parser.add_argument('--seed', type=int, dest='seed', help='Seed for the random number generator.')

args = parser.parse_args()
argv = sys.argv
narg = len(argv)-1
print "====================================="
print argv[0]
print "Version: ", version
print (datetime.datetime.now()).isoformat()
print "====================================="

# =======================================================================
# Addition from Jacopo Chevallard - 26/01/2017
if args.seed is not None:
    numpy.random.seed(seed=args.seed)
# =======================================================================

# =======================================================================
# Initialising the c_tableSNR object
# =======================================================================
table = c_tableSNR.c_tableSNR()

# =======================================================================
# Parsing the input arguments
# =======================================================================
print "# === Configuration and input paths"
fullInputFilename = args.filename
inputPath,inputFilename = os.path.split(fullInputFilename)
table.m_setKeyword('REFSRC', inputFilename)
print "# Path to the input filename: {:s}".format(inputPath)
print "# Input filename: {:s}".format(inputFilename)
dataPath = args.data
table.m_setKeyword('DATAPATH', dataPath)
print "# Path to the JWSTpytools data folder: {:s}".format(dataPath)
pceFolder = args.pcefolder
table.m_setKeyword('PCEDIR', pceFolder)
print "# Name of the folder containing the PCE file: {:s}".format(pceFolder)
pcePrefix = 'PCE-NIRSpec'
print "# Prefix used for the naming of the PCE file (hard-coded): {:s}".format(pcePrefix)
FWA = args.FWA
table.m_setKeyword('FWA', FWA)
print "# FWA configuration: {:s}".format(FWA)
GWA = args.GWA
table.m_setKeyword('GWA', GWA)
print "# GWA configuration: {:s}".format(GWA)
aperture = 'MOS1x3'
table.m_setKeyword('APERTURE', aperture)
print "# Aperture type (hard-coded): {:s}".format(aperture)

print "# === Source categories"
sourceSpatialType = args.spatype
table.m_setKeyword('SPATYPE', sourceSpatialType)
print "# Source spatial category: {:s}".format(sourceSpatialType)
sourceSpectralType = 'CONT'
table.m_setKeyword('SPECTYPE', sourceSpectralType)
print "# Source spectral category (hard-coded): {:s}".format(sourceSpectralType)

# =======================================================================
# Addition from Jacopo Chevallard - 26/01/2017
# You also modified the file /Users/jchevall/JWST/code/JWSTpylib-1.0.4/JWSTpylib/sensitivity/c_tableSNR.py to add the table entries below

if args.sersic is not None and args.effective_radius is not None: 
    table.m_setKeyword('SERSIC', args.sersic)
    table.m_setKeyword('R_EFF', args.effective_radius)
# =======================================================================

print "# === Exposure parameters"

numberOfExposures = args.nexp
if (numberOfExposures < 1):
	print "Non-valid number of MULTIACCUM22x4 exposures on input ({:d}). It should be strictly larger than 0.".format(numberOfExposures)
	raise ValueError
table.m_setKeyword('NEXP', numberOfExposures)
print "# Number of MULTIACCUM22x4 exposures: {:d}".format(numberOfExposures)

sum_spec = 1
table.m_setKeyword('SUMBOXX', sum_spec)
print "# Size of the summation along the spectral direction (hard-coded): {:d} pixels.".format(sum_spec)

sum_spa = int(args.sumy)
if (sum_spa < 1):
	print "Non-valid size for the size of spatial direction summation box ({:d} pixels). It should be strictly larger than 0.".format(sum_spa)
	raise ValueError
if ((sum_spa != 4) and (sum_spa != 3)):
	print "# ==WARNING== summation box loss correction is NOT valid for sum_spa != 3 or 4. It will be disabled."
table.m_setKeyword('SUMBOXY', sum_spa)
print "# Size of the summation box along the spatial direction: {:d} pixels.".format(sum_spa)

numberOfBackgroundElements = args.nb
if (numberOfBackgroundElements < 1):
	print "Non-valid number of background elements on input ({:d}). It should be strictly larger than 0.".format(numberOfBackgroundElements)
	raise ValueError
table.m_setKeyword('NBCKGD', numberOfBackgroundElements)
print "# Number of background elements: {:d}".format(numberOfBackgroundElements)
print "# CAUTION: the size of each background element is identifical to the size of the summation box."

# =======================================================================
# Output file parameters
# =======================================================================
print "# === Path and prefix when generating the output files."
outputPath = args.out
print "# Output path: {:s}".format(outputPath)
prefix = args.prefix
print "# Output file name prefix: ", prefix

# =======================================================================
# Detector noise & cosmetics
# =======================================================================
print "# === Detector noise and cosmetics parameters."
totalNoise = float(args.totalnoise)
if (totalNoise < 0.0):
	print "Non-valid total noise value on input ({:5.2f}). It should be larger or equal to 0.".format(totalNoise)
	raise ValueError
table.m_setKeyword('TOTNOISE', totalNoise)
print "# Total noise for a MULTIACCUM-22x4 exposure: {:5.2f} electrons (1 sigma).".format(totalNoise)

darkCurrent = float(args.darkcurrent)
if (darkCurrent < 0.0):
	print "Non-valid dark current value on input ({:5.2f}). It should be larger or equal to 0.".format(darkCurrent)
	raise ValueError
table.m_setKeyword('DARK', pceFolder)
print "# Average dark current rate: {:5.2f} electrons/s.".format(darkCurrent)

rejx = int(args.rejx)
if (rejx < 1):
	print "Non-valid size along the spectral direction for the rejection box ({:d} pixels). It should be strictly larger than 0.".format(rejx)
	raise ValueError
table.m_setKeyword('REJX', rejx)
print "# Size of the rejection box along the spectral direction: {:d} pixels.".format(rejx)

rejy = int(args.rejy)
if (rejy < 1):
	print "Non-valid size along the spatial direction for the rejection box ({:d} pixels). It should be strictly larger than 0.".format(rejy)
	raise ValueError
table.m_setKeyword('REJY', rejx)
print "# Size of the rejection box along the spatial direction: {:d} pixels.".format(rejy)

badPixFrac = float(args.badpixfrac)
if ((badPixFrac < 0.0) or (badPixFrac > 1.0)):
	print "Non-valid bad pixel fraction on input ({:5.2f}). It should be in the range[0,1].".format(badPixFrac)
	raise ValueError
table.m_setKeyword('BDPIXF', badPixFrac)
print "# Fraction of bad (dead + hot) pixels: {:5.2f} %.".format(1e2*badPixFrac)
openPixFrac = float(args.openpixfrac)
if ((openPixFrac < 0.0) or (openPixFrac > 1.0)):
	print "Non-valid open pixel fraction on input ({:5.2f}). It should be in the range[0,1].".format(openPixFrac)
	raise ValueError
table.m_setKeyword('OPPIXF', openPixFrac)
print "# Fraction of open pixels: {:5.2f} %.".format(1e2*openPixFrac)

cosmetics = c_detectorCosmetics.c_detectorCosmetics()
cosmetics.m_setFractions(badPixFrac, openPixFrac)
badPixelProbabilities = cosmetics.m_computeProbabilities(rejx, rejy)
print "# Probability that a 1x3 nodding pattern is fully usable               : {:5.4f}".format(badPixelProbabilities[0])
print "# Probability that one element of a 1x3 nodding pattern is not usable  : {:5.4f}".format(badPixelProbabilities[1])
print "# Probability that two elements of a 1x3 nodding pattern are no usable : {:5.4f}".format(badPixelProbabilities[2])
print "# Probability that a 1x3 nodding pattern is not uable at all           : {:5.4f}".format(badPixelProbabilities[3])
exposureTimeFactor = 1.0 * badPixelProbabilities[0] + 2.0/3.0 * badPixelProbabilities[1] + 1.0/3.0 * badPixelProbabilities[2]
table.m_setKeyword('TFACTOR', exposureTimeFactor)
print "# Exposure time reduction factor: {:5.3f}".format(exposureTimeFactor)
print "====================================="
print

# =======================================================================
# Summation box losses
# =======================================================================
if (sourceSpatialType == 'PS'):
	print "# PS case - Preparing the summation box losses (only valid for the 2x3 or 2x4 options)"
	boxLosses = c_summationBoxLosses.c_summationBoxLosses()
	boxLosses.m_setSummationBoxSize(sum_spec, sum_spa)
	losses = boxLosses.m_computeLosses(sourceSpectralType)
else:
	print "# ES case - No summation box losses (ignoring the shadow of the bars along the spatial direction)."
	losses = c_transmissionFunction.c_transmissionFunctionFromArrays(1e-6*numpy.array([0.6,5.3]), numpy.array([0.0,0.0]))

fracSummationBox = 1.0 - losses

# =======================================================================
# Loading the effiency of JWST NIRSpec (should include everything:
#	OTE + NIRSpec (including optical train, slit losses and detectors)
# =======================================================================
print "# Loading the JWST/NIRSpec PCE and spectral resolution values from disk."
pceobject = c_pce.c_pce()
pceobject.m_setTelescopeEffectiveArea(25.0)
table.m_setKeyword('ATEL', pceobject.Atel)
print "# Telescope ceffective area: {:5.2f} m2".format(pceobject.Atel)
pceobject.m_setPixelArea(2.35e-13)
table.m_setKeyword('APIX', pceobject.Apix)
print "# Pixel area on the sky: {:8.4e} sr".format(pceobject.Apix)
print "# Pixel area on the sky: {:8.4e} mas2".format(pceobject.Apix * (206265e3)**2)
print "# Pixel size on the sky: {:5.1f} mas".format(math.sqrt(pceobject.Apix * (206265e3)**2))
pceobject.m_setDataPath(dataPath)
pceobject.m_getSpectroscopicPCE(pceFolder, pcePrefix, aperture, sourceSpatialType, FWA, GWA)
table.m_setKeyword('ORDER', pceobject.order)

pcebackground = c_pce.c_pce()
pcebackground.m_setTelescopeEffectiveArea(25.0)
pcebackground.m_setPixelArea(2.35e-13)
pcebackground.m_setDataPath(dataPath)
pcebackground.m_getSpectroscopicPCE(pceFolder, pcePrefix, aperture, 'ES', FWA, GWA)

# =======================================================================
# Wavelength range
# =======================================================================
if (GWA == 'PRISM'):
	lbda_min = 0.60001e-6
	lbda_max = 5.0e-6
elif (GWA == 'G140M' or GWA == 'G140H'):
	if (FWA == 'F070LP'):
		lbda_min = 0.7e-6
		lbda_max = 1.2e-6
	else:
		lbda_min = 1.0e-6
		lbda_max = 1.8e-6
elif (GWA == 'G235M' or GWA == 'G235H'):
	lbda_min = 1.7e-6
	lbda_max = 3.1e-6
elif (GWA == 'G395M'):
	lbda_min = 2.9e-6
	lbda_max = 5.2e-6
elif (GWA == 'G395H'):
	if (FWA == 'F110W'):		
		lbda_min = 1.0e-6
		lbda_max = 1.25e-6
	else:
		lbda_min = 2.9e-6
		lbda_max = 5.2e-6
	
# =======================================================================
# Performance computation
# =======================================================================
# Exposure times in s
ng = 22
nf = 4
tf = 10.74
table.m_setKeyword('FLAGWIN', False)
table.m_setKeyword('WSIZEX', 2048)
table.m_setKeyword('WSIZEY', 2048)
table.m_setKeyword('TF', tf)
table.m_setKeyword('NR1', 1)
table.m_setKeyword('NR2', 2)
table.m_setKeyword('NF', 4)
table.m_setKeyword('NG', 22)
table.m_setKeyword('NINT', 1)

print "# Number of groups in an individual exposure: ng = {:d}".format(ng)
print "# Number of frames per group: nf = {:d}".format(nf)
print "# Time needed to reset or read a full frame: {:8.5f} seconds.".format(tf)
teff = (ng*nf-1) * tf
print "# Effective integration time for each individual exposure: teff = (ng*nf-1)*tf = {:8.5f} seconds.".format(teff)
ttot = (ng*nf+2) * tf
print "# Total exposure time (assuming that 2 resets are performed at the beginning): ttot = (ng*nf+2) * tf = {:8.5f} seconds.".format(ttot)
tc = 1.5e4
# In the technical note of P. Jakobsen this effect is described but finally
# not taken into account
# texp = nexp * tc * (1. - math.exp(-teff / tc))
texp = numberOfExposures * teff * exposureTimeFactor
table.m_setKeyword('TTOT', texp)
print "# Effective total on-source time (after correction by the detector cosmetics factor): {:8.5f} seconds.".format(texp)
# Detector parameters
readout = math.sqrt(totalNoise**2 - darkCurrent * teff)
table.m_setKeyword('RDNOISE', readout)
print "# Effective readout noise (artificial construct): {:5.2f} electrons/pixel.".format(readout)
# Gamma coefficients
gamma_ff = args.gammaflatfield
table.m_setKeyword('GAMMAFF', gamma_ff)
print "# Accuracy of the flat-field correction: {:5.2f} %%".format(1e2 * gamma_ff)
gamma_dark = args.gammadark
table.m_setKeyword('GAMMAD', gamma_dark)
print "# Accuracy of the dark current subtraction: {:5.2f} %%".format(1e2 * gamma_dark)
# pylab.plot(1e6*wave, delta_lbda, wave, delta_lbda_bis)
# pylab.show()
slit_width = 2
table.m_setKeyword('SWIDTH', slit_width)
print "# Slit width: {:d} pixels".format(slit_width)
npix = sum_spec * sum_spa

# =======================================================================
# Loading the input spectrum
# =======================================================================
print "# Loading the input spectrum."
inputSpectrum = c_spectrum.c_spectrum()
inputSpectrum.m_readFromSimpleFITS(fullInputFilename)
inputCentralWavelength, inputRebinGrid, inputRebinGridStepSize, inputValues = inputSpectrum.m_getRebinGrids()

# =======================================================================
# Preparing the output rebin grid
# =======================================================================
print "# Generating the output rebin grids."
temporarySpectrum = c_spectrum.c_spectrum()
temporarySpectrum.m_createSpectrumFromResolutionCurve(lbda_min, lbda_max, pceobject.resolution, elemsize=2.2)
outputCentralWavelength, outputRebinGrid, outputRebinGridStepSize, outputValues = temporarySpectrum.m_getRebinGrids()
wave = numpy.copy(outputCentralWavelength)
delta_lbda = numpy.copy(outputRebinGridStepSize)
boxCorrection = fracSummationBox(wave)


# Configuration parameters
# delta_lbda = wave / (2.2 * pceobject.resolution(wave))

# =======================================================================
# Rebinning the input spectrum
# =======================================================================
rebinnedValues = f_interpolation.f_rebinLinear1D(inputValues * inputRebinGridStepSize, outputRebinGrid, source=inputRebinGrid) / outputRebinGridStepSize


# =======================================================================
# Addition from Jacopo Chevallard - 26/01/2017
# Calculating slit losses for extended soruces, using tabulated data computed by M. Maseda
slit_losses = numpy.ones(len(wave))
if args.sersic is not None and args.effective_radius is not None: 
    m = MSAThroughput(os.path.join(dataPath, "slit_losses"))
    # Get the slit losses wrt to a centered point source, for a given Sersic
    # index and effective_radius radius, but averaged over all positions within
    # the open slit
    slit_losses = m.get_throughput(wl=wave*1.E+06, Sersic=args.sersic, 
            effective_radius=args.effective_radius)

rebinnedValues *= slit_losses
# =======================================================================

# =======================================================================
# Object
# =======================================================================

# =======================================================================
# Modified by Jacopo Chevallard - 26/01/2017 - to add slit losses
if (sourceSpatialType == 'PS'):
	# input spectrum in Jy
	object_elec = 1e-26 / h * rebinnedValues * pceobject.Atel * sum_spec * delta_lbda / wave * texp * pceobject.pce(wave) * boxCorrection
else:
	area = pceobject.Apix * slit_width * sum_spa * 206265**2
	# input spectrum in Jy arcsec-2 - no box correction for the ES case
	object_elec = 1e-26 / h * rebinnedValues *  pceobject.Atel * sum_spec * delta_lbda / wave * texp * area * pceobject.pce(wave)
# =======================================================================

# =======================================================================
# Zodiacal light
# =======================================================================
zodlight = 1.2 * zod.c_zodiacalLight(256.)
straylight = c_straylightOTE.c_straylightOTEFromArrays(numpy.array([0.6,1.0,2.0,3.0,5.0]), numpy.array([0.091e6, 0.091e6, 0.091e6, 0.070e6, 0.070e6]))
background = zodlight + straylight

# =======================================================================
# Noise
# =======================================================================
background_elec = background(1e6*wave) * pcebackground.Atel * texp * slit_width * sum_spa * pcebackground.Apix * 1e6 * delta_lbda * sum_spec * pcebackground.pce(wave)
gamma = gamma_ff**2 / (npix * numberOfExposures)
varback = background_elec + background_elec * background_elec * gamma
vardet = numberOfExposures * npix * (darkCurrent * teff + readout**2 + (darkCurrent * teff * gamma_dark)**2)
varobject = object_elec + object_elec * object_elec * gamma
variance = vardet + varback + varobject
if (sourceSpatialType == 'PS'):
	# Converting the noise (1 sigma) into Jy
	conversionFactor = 1e-26 / h * pceobject.Atel * sum_spec * delta_lbda / wave * texp * pceobject.pce(wave) * boxCorrection
	absoluteNoise = numpy.sqrt(variance) / conversionFactor
	validIndices = numpy.where(conversionFactor == 0.0)
	absoluteNoise[validIndices] = 0.0
else:
	area = pceobject.Apix * slit_width * sum_spa * 206265**2
	# Converting the noise into Jy arcsec-2
	conversionFactor = 1e-26 / h * pceobject.Atel * sum_spec * delta_lbda / wave * texp * area * pceobject.pce(wave)
	absoluteNoise = numpy.sqrt(variance) / conversionFactor
	validIndices = numpy.where(conversionFactor == 0.0)
	absoluteNoise[validIndices] = 0.0

# =======================================================================
# Comparing different contributions to the noise variance
# =======================================================================
line0 = pylab.plot(1e6*wave, varobject)
pylab.setp(line0, 'linestyle','-')
pylab.setp(line0, 'linewidth',2.0)
line1 = pylab.plot(1e6*wave, varback)
pylab.setp(line1, 'linestyle','-')
pylab.setp(line1, 'linewidth',2.0)
line2 = pylab.plot(1e6*wave, 0. * wave + vardet)
pylab.setp(line2, 'linestyle','--')
pylab.setp(line2, 'linewidth',2.0)
pylab.legend((line0[0],line1[0], line2[0]), ('Variance of the object', 'Variance of the background related noise\n(zodiacal light + OTE straylight)', 'Variance of the detector related noise'), loc='upper right', prop={"size":10})
pylab.grid(True)
pylab.xlabel('Wavelength (microns)')
pylab.ylabel('Variance (electrons**2)')
pylab.title('Contribution of background and detector noise to the noise variance\n (MOS mode - {:s}/{:s} - summation over {:d}x{:d} pixels)'.format(FWA, GWA, sum_spec, sum_spa), fontsize=11)
filename = '{:s}_variance_{:s}_{:s}_{:s}.pdf'.format(prefix, sourceSpatialType, FWA, GWA)
pylab.savefig(os.path.join(outputPath, filename))
# pylab.show()
pylab.close()

# =======================================================================
# SGenerating the noisy rebinned spectrum
# =======================================================================
snr = object_elec / numpy.sqrt(variance)
numberOfWavelengths = outputCentralWavelength.size
noise = rebinnedValues * numpy.random.normal(0.0, 1.0, numberOfWavelengths) / snr
absnoise = absoluteNoise * numpy.random.normal(0.0, 1.0, numberOfWavelengths)
noisyRebinnedValues = rebinnedValues + absnoise

# line0 = pylab.plot(1e6*wave, noise)
# pylab.setp(line0, 'linestyle','-')
# pylab.setp(line0, 'linewidth',2.0)
# line1 = pylab.plot(1e6*wave, absnoise)
# pylab.setp(line1, 'linestyle','-')
# pylab.setp(line1, 'linewidth',2.0)
# pylab.legend((line0[0],line1[0]), ('Noise through SNR', 'Noise through absolute noise'), loc='upper right', prop={"size":10})
# pylab.grid(True)
# pylab.xlabel('Wavelength (microns)')
# pylab.ylabel('Sigma (Jy or Jy arcsec-2)')
# pylab.title('COmparison of the noise computed using two different methods\n (MOS mode - {:s}/{:s} - summation over {:d}x{:d} pixels)'.format(FWA, GWA, sum_spec, sum_spa), fontsize=11)
# filename = '{:s}_variance_{:s}_{:s}_{:s}.pdf'.format(prefix, sourceSpatialType, FWA, GWA)
# pylab.savefig(os.path.join(outputPath, filename))
# # pylab.show()
# pylab.close()

# =======================================================================
# Signal to noise ratio
# =======================================================================
table.wavelength = numpy.copy(outputCentralWavelength)
table.resolution = numpy.copy(pceobject.resolution(outputCentralWavelength))
table.deltaWavelength = numpy.copy(outputRebinGridStepSize)
table.minimumWavelength = outputRebinGrid[:-1]
table.maximumWavelength = outputRebinGrid[1:]
table.electronRate = object_elec / texp
table.numberOfAccumulatedElectrons = numpy.copy(object_elec)
table.SNR = numpy.copy(snr)

# =======================================================================
# Modified by Jacopo Chevallard - 26/01/2017
table.noise = absoluteNoise / slit_losses
table.noiselessSpectrum = numpy.copy(rebinnedValues) / slit_losses
table.noisySpectrum = numpy.copy(noisyRebinnedValues) / slit_losses
# =======================================================================

filename = '{:s}_snr_{:s}_{:s}_{:s}.fits'.format(prefix, sourceSpatialType, FWA, GWA)
author = "{:s} - {:s}".format(argv[0], version)
reference = inputFilename
description = "JWST/NIRSpec simulated spectrum"
table.m_writeToFITS(os.path.join(outputPath, filename), author, reference, description)
