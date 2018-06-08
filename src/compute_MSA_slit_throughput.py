import os
import numpy as np
from astropy.io import fits
from scipy import interpolate
import re

def isclose(a, b, rel_tol=1e-07, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

class MSAThroughput(object):

    def __init__(self, folder, aperture_type="OPEN", throughput_type="relative"):

        aperture_types = ("PITCH", "OPEN", "MID", "CONST", "TIGHT")
        if aperture_type not in aperture_types:
             raise Exception('Aperture type`' + aperture_type + '` not supported!')

        file_list = list()

        # Read in the model
        for file in sorted(os.listdir(folder)):
            if file.endswith(".fits") and aperture_type in file:
                file_list.append(file)

        sersic = list()
        for file in file_list:
            f = file.split('_')[0]
            f = map(int, re.findall('\d+', f))
            sersic.append(f[0])

        self.throughput = None

        self.x = np.array(sersic)

        for ix, file in enumerate(file_list):

            #print "file: ", file

            full_path = os.path.join(folder, file)
            hdulist = fits.open(full_path)

            wl = hdulist[1].data["wavelength"]
            unique_wl = reduce(lambda l, x: l if x in l else l+[x], wl, [])
            self.z = np.array(unique_wl)*1.E+06
            r_eff = hdulist[1].data["r_eff"]
            if throughput_type == "relative":
                throughput = hdulist[1].data["correctedthroughput"]
            elif throughput_type == "total":
                throughput = hdulist[1].data["throughput"]
    
            unique_r_eff = reduce(lambda l, x: l if x in l else l+[x], r_eff, [])

            self.y = np.array(unique_r_eff)

            if self.throughput is None:
                self.throughput = np.zeros((self.x.size, self.y.size, self.z.size))
            
            for iy, r in enumerate(unique_r_eff):
                ok = np.where(r_eff == r)[0]
                t = throughput[ok]
                    
                self.throughput[ix,iy,:] = t    
                #print "ix, iy: ", ix, iy
                #print "throughput: ", t

            hdulist.close()

        self.abscissa = (self.x, self.y, self.z)

    def get_throughput(self, wl, Sersic=4, effective_radius=0.1):

        if isclose(Sersic, self.x[0]) or isclose(Sersic, self.x[-1]):
            pass
        elif Sersic < self.x[0]:
            raise ValueError('Input Sersic index cannot be < ' + str(self.x[0]) + '!')
        elif Sersic > self.x[-1]:
            raise ValueError('Input Sersic index cannot be > ' + str(self.x[-1]) + '!')

        if isclose(effective_radius, self.y[0]) or isclose(effective_radius, self.y[-1]):
            pass
        elif effective_radius < self.y[0]:
            raise ValueError('Input effective radius cannot be < ' + str(self.y[0]) + '!')
        elif effective_radius > self.y[-1]:
            raise ValueError('Input effective radius cannot be > ' + str(self.y[-1]) + '!')

        # Linear interpolate at the right Sersic index
        iy = np.searchsorted(self.y, effective_radius)
        ix = np.searchsorted(self.x, Sersic)

        if isclose(Sersic, self.x[ix]):
            if isclose(effective_radius, self.y[iy]):
                fy = self.throughput[ix,iy,:]
            else:
                fy = (self.y[iy+1]-effective_radius)/(self.y[iy+1]-self.y[iy])*self.throughput[ix,iy,:] \
                        + (effective_radius-self.y[iy])/(self.y[iy+1]-self.y[iy])*self.throughput[ix,iy+1,:]
        else:
            if isclose(effective_radius, self.y[iy]):
                fy = (self.x[ix+1]-Sersic)/(self.x[ix+1]-self.x[ix])*self.throughput[ix,iy,:] \
                        + (Sersic-self.x[ix])/(self.x[ix+1]-self.x[ix])*self.throughput[ix+1,iy,:]
            else:
                fy1 = (self.x[ix+1]-Sersic)/(self.x[ix+1]-self.x[ix])*self.throughput[ix,iy,:] \
                        + (Sersic-self.x[ix])/(self.x[ix+1]-self.x[ix])*self.throughput[ix+1,iy,:]

                fy2 = (self.x[ix+1]-Sersic)/(self.x[ix+1]-self.x[ix])*self.throughput[ix,iy+1,:] \
                        + (Sersic-self.x[ix])/(self.x[ix+1]-self.x[ix])*self.throughput[ix+1,iy+1,:]

                fy = (self.y[iy+1]-effective_radius)/(self.y[iy+1]-self.y[iy])*fy1 \
                        + (effective_radius-self.y[iy])/(self.y[iy+1]-self.y[iy])*fy2


        f = interpolate.interp1d(self.z, fy, kind='cubic')
    
        # Avoid extrapolation errors
        i0=0
        i1=len(wl)+1
        f_wl = np.zeros(len(wl))
        if wl[0] < self.z[0]:
            i0 = np.searchsorted(wl, self.z[0])

        if wl[-1] > self.z[-1]:
            i1 = np.searchsorted(wl, self.z[-1])

        f_wl[i0:i1] = f(wl[i0:i1])

        if wl[0] < self.z[0]:
            f_wl[0:i0] = f(self.z[0])

        if wl[-1] > self.z[-1]:
            f_wl[i1:] = f(self.z[-1])

        return f_wl

        #grid_z2 = griddata(points, values, (sersic, effective_radius, wl), method='cubic')

if __name__ == '__main__':

    m = MSAThroughput("/Users/jchevall/JWST/code/JWSTpytools-0.0.3/data/slit_losses", throughput_type="total")

    # Create table of slit losses at 4.5 micron for Emma
    dropout_bands = ('B', 'V', 'I', 'Z', 'Y')
    redshifts = (8, 4, 5, 6, 7, 8)
    sizes = (0.06, 0.10, 0.11, 0.12, 0.13, 0.145)
    sersic = 1.0
    wl = np.array([0.6, 2.5,  5.3])

    print "redshift wl/micron size/arcsec throughput"
    for size, redshift in zip(sizes, redshifts):
        throughput = m.get_throughput(wl=wl, Sersic=sersic, effective_radius=size)
        print redshift, wl, size, throughput

    # Create table of slit losses at redshifted Halpha wl
    wl = np.array((0.6563, ))
    print "redshift wl/micron size/arcsec throughput"
    for size, redshift in zip(sizes, redshifts):
        Ha_wl = wl*(1.0+redshift)
        throughput = m.get_throughput(wl=Ha_wl, Sersic=sersic, effective_radius=size)
        print redshift, Ha_wl[0], size, throughput[0]


    #wl = np.arange(0.7, 4.5, 0.25)
    #through = m.get_throughput(wl=wl, Sersic=1, effective_radius=0.1)
    #print wl, through
    #for w, t in zip(wl, through): 
    #    print "wl, t: ", w, t

