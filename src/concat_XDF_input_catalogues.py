import os
import logging
from collections import OrderedDict
import json
import numpy as np
from scipy.interpolate import interp1d
from astropy.io import fits

if __name__ == '__main__':

    output_folder = "/Users/jchevall/JWST/Simulations/XDF/ALL_DROPOUTS"

    dropout_bands = ("B", "V", "I", "Z", "Y")
    catalogue_name = 'input_SEDs_MC_0.fits'

    # Initialize a new (empty) primary HDU for your output FITS file
    new_hdulist = fits.HDUList(fits.PrimaryHDU())

    hdu_names = ('GALAXY PROPERTIES', 'STAR FORMATION', 'STAR FORMATION BINS', 
            'DUST ATTENUATION', 'NEBULAR EMISSION', 'ABSOLUTE MAGNITUDES', 'APPARENT MAGNITUDES')

    # Create a new multi-extension FITS file that will contain the FITS tables
    # of each dropout catalogue

    # First, determine the total number of rows needed
    nrows = 0
    for band in dropout_bands:
        folder = '/Users/jchevall/JWST/Simulations/XDF/' + band + \
            '_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian/MC_0/input-SEDs'

        file_name = os.path.join(folder, catalogue_name)
        hdulist = fits.open(file_name)

        nrows += len(hdulist[1].data.field(0))
        hdulist.close()

    # You consider the first file in the list and use as a "mold" to create
    # the structure (binary tables and their columns) of the output FITS file
    folder = '/Users/jchevall/JWST/Simulations/XDF/' + dropout_bands[0] + \
        '_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian/MC_0/input-SEDs'
    file_name = os.path.join(folder, catalogue_name)
    hdulist = fits.open(file_name)

    # Now you cycle over all extension and columns that you want to put in
    # the summary catalogue
    for hdu_name in hdu_names:
        new_columns = list()

        # The first column of each output extension contains the object ID
        new_columns.append(fits.Column(name='ID', format='20A'))

        # By default you take all columns in that extensions
        col_names = hdulist[hdu_name].columns.names

        for col in hdulist[hdu_name].columns:

            new_columns.append(fits.Column(name=col.name,
                format=col.format, unit=col.unit))

        # Create the "column definition"
        cols_ = fits.ColDefs(new_columns)

        # Create the actual binary table, with the correct number of rows
        # to accomodate all objects
        new_hdu = fits.BinTableHDU.from_columns(cols_, nrows=nrows)
        new_hdu.name = hdu_name

        # And finally append the newly created bunary table to the hdulist
        # that will be printed to the ouput FITS file
        new_hdulist.append(new_hdu)

        hdulist.close()

    # Now you can go through each file, and compute the required quantities
    i0 = 0
    i1 = 0
    for band in dropout_bands:
        folder = '/Users/jchevall/JWST/Simulations/XDF/' + band + \
            '_DROPOUTS/ineb_Jan16_logU_xid_delayed_SFR-Gaussian_max_age-Gaussian/MC_0/input-SEDs'

        file_name = os.path.join(folder, catalogue_name)
        hdulist = fits.open(file_name)
        nrows = len(hdulist[1].data.field(0))
        ID = np.array([band + '_' + str(i+1) for i in range(nrows)])
        i1 += nrows

        print "nrows: ", nrows, len(new_hdulist[hdu_name].data['ID'][i0:i1])

        for hdu_name in hdu_names:

            new_hdulist[hdu_name].data['ID'][i0:i1] = ID

            for col in hdulist[hdu_name].columns.names:
                if col not in new_hdulist[hdu_name].columns.names:
                    continue
                new_hdulist[hdu_name].data[col][i0:i1] = hdulist[hdu_name].data[col]

        hdulist.close()

        i0 += nrows

    name = os.path.join(output_folder, catalogue_name)
    new_hdulist.writeto(name, clobber=True)


