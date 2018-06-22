#!/usr/bin/env python

from collections import OrderedDict
import argparse
import json
from astropy.io import fits
import numpy as np

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-i', '--input',
        help="",
        type=str,
        action="store", 
        required=True,
        dest="input"
    )

    parser.add_argument(
        '-o', '--output',
        help="",
        type=str,
        action="store", 
        required=True,
        dest="output"
    )

    parser.add_argument(
        '--config', 
        help="",
        type=str,
        action="store", 
        required=True,
        dest="config_file"
    )

    args = parser.parse_args()

    hdulist = fits.open(args.input)

    with open(args.config_file) as f:    
        hdu_col = json.load(f, object_pairs_hook=OrderedDict)

    new_columns = list()

    for key, val in hdu_col.iteritems():
        if "extName" in val:
            ext = val["extName"]
        else:
            ext = "POSTERIOR PDF"

        if "colName" in val:
            col = val["colName"]
        else:
            col = key

        if "units" in val:
            units = val["units"]
        else:
            units = ""

        if "log10" in val:
            log10 = val["log10"]
        else:
            log10 = False

        column = hdulist[ext].columns[col]

        data = hdulist[ext].data[col]
        if log10:
            data = np.log10(data)

        print "key: ", key
        new_col = fits.Column(name=str(key), format=column.format, unit=units, 
                array=data)
                    
        new_columns.append(new_col)

    colsDefs = fits.ColDefs(new_columns)

    newhdulist = fits.HDUList(fits.PrimaryHDU())

    newhdu = fits.BinTableHDU.from_columns(colsDefs)

    newhdulist.append(newhdu)

    newhdu.writeto(args.output, overwrite=True)

    hdulist.close()

