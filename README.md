## Getting ready

In order to use the Python routines included in the ``src`` folder you must:

1. clone the [PyP-BEAGLE](https://github.com/jacopo-chevallard/PyP-BEAGLE) repository

2. set the ``PYP_BEAGLE`` environment variable to point to the location of the cloned ``PyP-BEAGLE`` repository

3. get the ``JWSTpytools-0.0.3`` package

4. get the ``JWSTpylib-1.0.4`` package and install it with the command ``python setup.py install``

5. set the ``JWSTPYTOOLS`` environment variable to point to the location of the ``JWSTpytools-0.0.3`` package

## Simulating NIRSpec observations

* from the ``src`` folder, type 
 ```
 ./make_ETC_simulations.py --help
 ```
 and you should get all the necessary information to start creating simulated NIRSpec observations starting from an input catalogue of galaxy SEDs created with Beagle. 
