#!/usr/bin/env python

from collections import OrderedDict
import argparse
import os
from natsort import natsorted, ns
import json
from astropy.io import fits
from bisect import bisect_left
import numpy as np

def get_line_SN_OLD(file_name, line_wl):

    hdulist = fits.open(file_name)
    redshift = hdulist[1].header['redshift']

    if isinstance(line_wl, (list, tuple, np.ndarray)):

        SN = np.zeros(len(line_wl))
        SN_1pixel = np.zeros(len(line_wl))

        for ii, wl in enumerate(line_wl):

            wl *= (1.+redshift) * 1.E-10

            # Get the nearest wavelength to the redshifted line, to identify the pixel
            # where the emission line falls
            i = bisect_left(hdulist[1].data['wavelength'], wl)
            
            if wl > hdulist[1].data['maxw'][-1] or i >= len(hdulist[1].data['wavelength']):
                SN[ii] = -99.
                SN_1pixel[ii] = -99.
                continue
            
            if wl > hdulist[1].data['minw'][i] and wl <= hdulist[1].data['maxw'][i]:
                indx = i
            elif wl > hdulist[1].data['minw'][i-1] and wl <= hdulist[1].data['maxw'][i-1]:
                indx = i-1

            # Now get the nearest pixel to the line, *excluding* the pixel where the
            # line falls
            n = len(hdulist[1].data['wavelength'])
            if i == (n-1):
                j = i-1
            elif (wl-hdulist[1].data['wavelength'][i-1]) < (wl-hdulist[1].data['wavelength'][i+1]):
                j = i-1
            else:
                j = i+1

            # Get the flux in the pixel containing the line, by reading the "rebinned spectrum" flux
            signal = hdulist[1].data['RSPEC'][indx]

            # The noise estimate will be the sum in quadrature of the noise in the
            # pixel where the line falls and in the nearest pixel. This accounts for
            # the fact that the instrumental resolution of NIRSpec will always sprad an
            # unresolved spectral line (all emission lines from HII regions at R=100
            # and R=1000) between 1.3 and 2 pixel
            noise = np.sqrt((hdulist[1].data['RSPEC'][indx]/hdulist[1].data['SNR'][indx])**2 \
                    + (hdulist[1].data['RSPEC'][j]/hdulist[1].data['SNR'][j])**2)

            SN[ii] = signal/noise
            SN_1pixel[ii] = hdulist[1].data['SNR'][indx]

    else:

        wl = line_wl * (1.+redshift) * 1.E-10

        # Get the nearest wavelength to the redshifted line, to identify the pixel
        # where the emission line falls
        i = bisect_left(hdulist[1].data['wavelength'], wl)
        
        if wl > hdulist[1].data['maxw'][-1] or i >= len(hdulist[1].data['wavelength']):
            hdulist.close()
            return -99., -99.
        
        if wl > hdulist[1].data['minw'][i] and wl <= hdulist[1].data['maxw'][i]:
            indx = i
        elif wl > hdulist[1].data['minw'][i-1] and wl <= hdulist[1].data['maxw'][i-1]:
            indx = i-1

        # Now get the nearest pixel to the line, *excluding* the pixel where the
        # line falls
        n = len(hdulist[1].data['wavelength'])
        if i == (n-1):
            j = i-1
        elif (wl-hdulist[1].data['wavelength'][i-1]) < (wl-hdulist[1].data['wavelength'][i+1]):
            j = i-1
        else:
            j = i+1

        # Get the flux in the pixel containing the line, by reading the "rebinned spectrum" flux
        signal = hdulist[1].data['RSPEC'][indx]

        # The noise estimate will be the sum in quadrature of the noise in the
        # pixel where the line falls and in the nearest pixel. This accounts for
        # the fact that the instrumental resolution of NIRSpec will always sprad an
        # unresolved spectral line (all emission lines from HII regions at R=100
        # and R=1000) between 1.3 and 2 pixel
        noise = np.sqrt((hdulist[1].data['RSPEC'][indx]/hdulist[1].data['SNR'][indx])**2 \
                + (hdulist[1].data['RSPEC'][j]/hdulist[1].data['SNR'][j])**2)

        SN = signal/noise

        SN_1pixel = hdulist[1].data['SNR'][indx]

    hdulist.close()

    return SN, SN_1pixel


def get_lines_SN(file_name, lines):


    #print "file_name: ", file_name

    hdulist = fits.open(file_name)
    redshift = hdulist[1].header['redshift']

    wl = hdulist[1].data['wavelength'] * 1.E+10
    minw = hdulist[1].data['minw'] * 1.E+10
    maxw = hdulist[1].data['maxw'] * 1.E+10
    nwl = len(wl)
    dwl = hdulist[1].data['deltaw'] * 1.E+10
    fluxes = hdulist[1].data['FLUX_FLAMBDA']
    errors = hdulist[1].data['NOISE_FLAMBDA']

    SN = OrderedDict()
    for key, value in lines.iteritems():
        SN[key] = -99.99

    # Cycle across all the lines for which we compute the EWs
    for key, value in lines.iteritems():

        #print "key: ", key

        #### 
        wl_range = np.array(value["wl_range"]) * (1.+redshift)
        if wl_range[0] > wl[-1] or wl_range[1] < wl[0]:
            continue
        i0 = np.searchsorted(wl, wl_range[0]) ; i0 -= 1
        i1 = np.searchsorted(wl, wl_range[1])
        n_wl = i1-i0+1

        # Compute the average left continuum
        il0, flux_left = None, None
        if 'continuum_left' in value:
            continuum_left = np.array(value["continuum_left"]) * (1.+redshift)
            if continuum_left[0] >= wl[0] and continuum_left[1] <= wl[-1]:
                il0 = np.searchsorted(wl, continuum_left[0]) ; il0 -= 1
                il1 = np.searchsorted(wl, continuum_left[1])
                if il0 == il1:
                    il0 -= 1
        else:
            if i0 > 4:
                il0 = max(0, i0-5)
                il1 = il0 + 3

        if il0 is not None:
            while il1 >= i0:
                il1 -= 1
                il0 -= 1

        if il0 is not None:
            dwl_left = maxw[il1] - minw[il0]
            flux_left = np.sum(fluxes[il0:il1+1]*dwl[il0:il1+1]) / dwl_left
            err_flux_left = np.sqrt(np.sum((errors[il0:il1+1]*dwl[il0:il1+1])**2)) / dwl_left
            wl_left = 0.5*(maxw[il1] + minw[il0])
            #print "wl left, dwl: ", wl[il0:il1+1], dwl[il0:il1+1]
            #print "flux_left: ", wl_left, flux_left, flux_left/err_flux_left

        # Compute the average right continuum
        ir0, flux_right = None, None
        if 'continuum_right' in value:
            continuum_right = np.array(value["continuum_right"]) * (1.+redshift)
            if continuum_right[0] >= wl[0] and continuum_right[1] <= wl[-1]:
                ir0 = np.searchsorted(wl, continuum_right[0]) ; ir0 -= 1
                ir1 = np.searchsorted(wl, continuum_right[1])
                if ir0 == ir1:
                    ir0 -= 1
        else:
            if nwl-i1 > 4:
                ir1 = min(nwl-1, i1+5)
                ir0 = ir1 - 3

        if ir0 is not None:
            while ir0 <= i1:
                ir1 += 1
                ir0 += 1

        if ir0 is not None:
            dwl_right = maxw[ir1] - minw[ir0]
            flux_right = np.sum(fluxes[ir0:ir1+1]*dwl[ir0:ir1+1]) / dwl_right
            err_flux_right = np.sqrt(np.sum((errors[ir0:ir1+1]*dwl[ir0:ir1+1])**2)) / dwl_right
            wl_right = 0.5*(maxw[ir1] + minw[ir0])
            #print "wl right, dwl: ", wl[ir0:ir1+1], dwl[ir0:ir1+1]
            #print "flux_right: ", wl_right, flux_right, flux_right/err_flux_right

        # Approximate the continuum with a straght line
        if flux_left is not None and flux_right is not None:
            grad = (flux_right-flux_left)/(wl_right-wl_left)
            grad_error = np.sqrt(err_flux_left**2 + err_flux_right**2) / (wl_right-wl_left)
        else:
            grad = 0.
            grad_error = 0.
        #print "grad SN: ", grad, grad/grad_error

        if flux_left is not None and flux_right is not None:
            intercept = flux_right - grad*wl_right
            intercept_error = np.sqrt(err_flux_right**2 + (grad_error*wl_right)**2) 
        elif flux_left is not None:
            intercept = flux_left
            intercept_error = err_flux_left
        elif flux_right is not None:
            intercept = flux_right
            intercept_error = err_flux_right

        #print "intercept SN: ", intercept, intercept/intercept_error

        # Linear function approximating the continuum
        integrand = fluxes[i0:i1+1] - (grad*wl[i0:i1+1]+intercept)
        #print "fluxes---> ", fluxes[i0:i1+1]
        #print "continuum-----> ", (grad*wl[i0:i1+1]+intercept)
        integrand_errors = errors[i0:i1+1]
        integrated_flux = np.sum(integrand*dwl[i0:i1+1])
        integrated_flux_error = np.sqrt(np.sum((integrand_errors*dwl[i0:i1+1])**2))

        if integrated_flux > 0.:
            SN[key] = integrated_flux/integrated_flux_error

        #integrand_errors = np.sqrt(errors[i0:i1+1]**2 + (grad_error*wl[i0:i1+1])**2 + intercept_error**2)
        #integrated_flux_error = np.sqrt(np.sum((integrand_errors*dwl[i0:i1+1])**2))
        #print "SN: ", SN[key],  integrated_flux/integrated_flux_error

        #print "--------> ", wl[i0], wl[i1], integrated_flux, integrated_flux_error, SN
    return SN

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--folder',
        help="Parent folder containing the NIRSpec simulations",
        dest="folder", 
        type=str,
        required=True
    )

    parser.add_argument(
        '--json-file',
        help="JSON file containing the list of lines (with relative wl) for which the S/N computation is required",
        dest="json_file", 
        type=str,
        required=True
    )

    args = parser.parse_args()    

    # Extract the Monte Carlo draw from the folder name
    for f in args.folder.split('/'):
        if 'MC_' in f:
            MC_draw = f.split('_')[1]

    # List all simulated spectra in <folder>/ETC-output
    folder = os.path.join(args.folder, "ETC-output")
    file_names = list()
    gratings = list()
    filters = list()
    IDs = list()
    for f in natsorted(os.listdir(folder)):
        if f.endswith(".fits"):
            file_names.append(f)
            s = os.path.splitext(f)[0].split('_')
            IDs.append(s[0])
            gratings.append(s[-1])
            filters.append(s[-2])

    
    cols = list()

    cols.append(fits.Column(name='ID', format='20A', array=IDs))
    cols.append(fits.Column(name='filter', format='20A', array=filters))
    cols.append(fits.Column(name='grating', format='20A', array=gratings))

    N = len(file_names)

    with open(args.json_file) as f:
        lines = json.load(f, object_pairs_hook=OrderedDict)

    data = dict()
    line_wl = list()
    for key, value in lines.iteritems():
        data[key] = np.zeros(N)
        #data[key+"_1pixel"] = np.zeros(N)
    
    for i, f in enumerate(file_names):

        file_name = os.path.join(args.folder, "ETC-output", f)
        SN = get_lines_SN(file_name, lines)
        for j, (key, value) in enumerate(lines.iteritems()):
            data[key][i] = SN[key]
            #data[key+"_1pixel"][i] = SN_1pixel[j]

    for key, value in lines.iteritems():
        cols.append(fits.Column(name=str(key), format='E', array=data[key]))
        #cols.append(fits.Column(name=str(key+"_1pixel"), format='E', array=data[key+"_1pixel"]))

    columns = fits.ColDefs(cols)
    new_hdu = fits.BinTableHDU.from_columns(columns)
    new_hdu.name = 'S_to_N'

    new_hdulist = fits.HDUList(fits.PrimaryHDU())
    new_hdulist.append(new_hdu)


    file_name = os.path.join(args.folder, "Emission_lines_observational_SN_MC_" + MC_draw + ".fits")
    new_hdulist.writeto(file_name, overwrite=True)
