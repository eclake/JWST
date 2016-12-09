## Getting ready

In order to use the Python routines included in the ``src`` folder you must:

1. clone the [PyP-BEAGLE](https://github.com/jacopo-chevallard/PyP-BEAGLE) repository

2. set the ``PYP_BEAGLE`` environment variable to point to the location of the cloned ``PyP-BEAGLE`` repository

3. get the ``JWSTpytools-0.0.3`` package

4. get the ``JWSTpylib-1.0.4`` package and install it with the command ``python setup.py install``

5. set the ``JWSTPYTOOLS`` environment variable to point to the location of the ``JWSTpytools-0.0.3`` package

## Usage examples

### Simulating NIRSpec observations

* get the help message by typing
 
 ```
 ./make_ETC_simulations.py --help
 ```
 and you should get all the necessary information to start creating simulated NIRSpec observations starting from an input catalogue of galaxy SEDs created with Beagle. 

### Extracting SEDs from Beagle output files 

* get the help message by typing
 
 ```
 ./extract_SEDs.py --help
 ```

* extract 1 realization from Beagle output, matching the redshift-dependent Speagle+14 mass-SFR relation
 
 ```
 ./extract_SEDs.py -r <path to directory containing Beagle results> -o <path to output directory>
 ```

* extract 1 realization from Beagle output, matching the redshift-dependent Speagle+14 mass-SFR relation **and** using the UVJ color diagram to separet SF and quiescent galaxies
 
 ```
 ./extract_SEDs.py -r <path to directory containing Beagle results> -o <path to output directory> --UVJ
 ```

* extract 1 realization from Beagle output, matching the redshift-dependent Speagle+14 mass-SFR relation but forcing the drawn solutions to have redshift > 4 (the ``step`` guarantees that is no solution with redshift > 4 is found, then the minimum redshift is iteratively lowered by ``step`` until some solutions are found)
 
 ```
 ./extract_SEDs.py -r <path to directory containing Beagle results> -o <path to output directory> --UVJ --params-ranges '{"redshift": {"colName":"redshift", "extName":"galaxy properties", "min":4.0, "step":0.1} }'
 ```

* extract 1 realization from Beagle output, matching the redshift-dependent Speagle+14 mass-SFR relation but forcing the drawn solutions to have redshift > 4 and mass > 1.0E+08
 
 ```
 ./extract_SEDs.py -r <path to directory containing Beagle results> -o <path to output directory> --UVJ --params-ranges '{"redshift": {"colName":"redshift", "extName":"galaxy properties", "min":4.0, "step":0.1}, "mass": {"colName":"M_tot", "extName":"galaxy properties", "min":1.0E+08}}'
 ```

### Building Value Added Catalogues

The script will match the coordinates (read from the RA and DEC columns in the FITS tables) of the sources in the input photometric catalogue with those of the sources of the photometric catalogue fitted by Beagle, and output a "value added catalogue".

* get the help message by typing
 
 ```
 ./make_value_added_catalogue.py --help
 ```

* If the input catalogue and the catalogue fitted with Beagle have ra/dec expressed in degrees, use the command

 ```
 ./make_value_added_catalogue.py -i <path to input photometric catalogue> -r <path to directory containing Beagle results> -p <name of Beagle parameter file>
 ```

* If the input catalogue ra/dec are expressed in degrees, and the catalogue fitted with Beagle have ra expressed in hours and dec in degrees, use the command 

 ```
 ./make_value_added_catalogue.py -i <path to input photometric catalogue> -r <path to directory containing Beagle results> -p <name of Beagle parameter file> --beagle-coord hourangle dec
 ```
