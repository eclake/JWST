# External modules.
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import brentq
from scipy.integrate import simps

# Internal modules.
import WeightedKDE

# This value is the minimum relative probability used when setting the axes
# limits and building the grid of points where to compute the PDF

class CredibleInterval:

    def __init__(self, data, probability):

        # Copy the posterior probability, likelihood and parameter values
        self.data = np.array(data)

        self.probability = np.array(probability)

        if self.data.ndim == 1:
            self.ComputeKDE1D()
        elif self.data.ndim == 2:
            self.ComputeKDE2D()

    def ComputeKDE1D(self,
                    minRelProb=1.E-03,
                    nXgrid=200):

        # Minimum and maximum value of the x and y parameters
        loc = np.where(self.probability >= np.amax(self.probability)*minRelProb)[0]

        self.min_x = np.amin(self.data[loc])
        self.max_x = np.amax(self.data[loc])

        # Compute a continuous marginal probability by means of (weighted)
        # kernel density estimation
        self.kde_pdf = WeightedKDE.gaussian_kde( self.data, weights=self.probability)

        # Now consider a regular grid of the parameter x
        self.x_grid = np.linspace(self.min_x, self.max_x, nXgrid)

        self.kde_pdf_grid = self.kde_pdf(self.x_grid)

        # Reshape the pdf on the grid
        self.kde_pdf_grid = np.array(self.kde_pdf_grid)

        # Compute the PDF normalization so the integral over the whole parameter range is 1
        # Use Simpson integration
        self.kde_pdf_norm = simps(self.kde_pdf_grid, self.x_grid)

        # Normalize the PDF
        self.kde_pdf_grid /= self.kde_pdf_norm

    def GetCumulativeIntegral(self):

        self.cumul_pdf = np.cumsum(self.kde_pdf_grid)

        # Be sure the it sums up to 1
        self.cumul_pdf /= self.cumul_pdf[-1]

        # Compute interpolant of cumulative PDF (linear interpolation)
        self.interp_cumul_pdf = interp1d(self.cumul_pdf, self.x_grid)

    def GetMedian(self):
        """
        Arguments:

        """

        # Check if the interpolant of the cumulative integral of the marginal PDF already exists
        try:
            interp_cumulpdf = self.interp_cumul_pdf
        except AttributeError:
            self.GetCumulativeIntegral()
            interp_cumulpdf = self.interp_cumul_pdf

        # Posterior median is the value of param that corresponds to a cumulative
        # probability = 0.5
        median = self.interp_cumul_pdf(0.5)

        return median

    def Get1DCredibleRegion(self, levels):
        """ 

        Arguments:

        """

        # Check if the interpolant of the cumulative integral of the marginal PDF already exists
        try:
            interp_cumul_pdf = self.interp_cumul_pdf
        except AttributeError:
            self.GetCumulativeIntegral()
            interp_cumul_pdf = self.interp_cumul_pdf

        # Check if the posterior median has been computed or not
        try:
            posterior_median = self.posterior_median
        except AttributeError:
            posterior_median = self.interp_cumul_pdf(0.5)

        # Plot credible regions intervals above data.
        limits = list()
        for level in levels:
            half_level = 0.5 * (1.0-level)

            low = interp_cumul_pdf(half_level)
            up = interp_cumul_pdf(1.-half_level)

            limits.append((low, up))

        return limits

    def ComputeKDE2D(self,
                    minRelProb=1.E-03,
                    nXgrid=40,
                    nYgrid=40):

        # Minimum and maximum value of the x and y parameters
        loc = np.where(self.probability >= np.amax(self.probability)*minRelProb)[0]

        self.min_x = np.amin(self.data[0,loc])
        self.max_x = np.amax(self.data[0,loc])

        self.min_y = np.amin(self.data[1,loc])
        self.max_y = np.amax(self.data[1,loc])

        # Compute a continuous marginal probability by means of (weighted)
        # kernel density estimation
        self.kde_pdf = WeightedKDE.gaussian_kde( self.data, weights=self.probability)

        # Now consider a regular grid of the parameters x and y and compute the joint
        # probability on the grid by using the above kernel pdf
        self.x_grid = np.linspace(self.min_x, self.max_x, nXgrid)
        self.y_grid = np.linspace(self.min_y, self.max_y, nYgrid)

        self.xx_grid, self.yy_grid = np.meshgrid(self.x_grid, self.y_grid)

        self.kde_pdf_grid =  self.kde_pdf((np.ravel(self.xx_grid), np.ravel(self.yy_grid)))

        # Reshape the pdf on the grid
        self.kde_pdf_grid = np.reshape(self.kde_pdf_grid, self.xx_grid.shape)

        # Compute the PDF normalization so the integral over the whole parameter range is 1
        # Use Simpson integration
        self.kde_pdf_norm = simps(simps(self.kde_pdf_grid,self.y_grid), self.x_grid)

        #print "Integral of the 2-D joint PDF: ", self.kde_pdf_norm

        # Normalize the PDF
        self.kde_pdf_grid /= self.kde_pdf_norm


    def GetProbabilityFor2DCredibleRegion(self,
            levels):

        """ Plot the posterior mean

        Arguments:

        CredibleRegions     -- Class containing the definitions for plotting credible regions

        """

        # Plot credible regions intervals above data.

        isocontourLevel = list()
        colors = list()
        labels= list()
        linewidths = list()
        max_pdf = max(np.ravel(self.kde_pdf_grid)) 

        for level in levels:
            
            # Parametric integral: you compute the integral for all values of
            # probability >= t
            def f_integral(t):
                temp_grid = np.copy(self.kde_pdf_grid)
                temp_grid[temp_grid < t] = 0.
                i = simps(simps(temp_grid,self.y_grid), self.x_grid)
                return i-level

            # Root finding algorithm to search for the zero of the function,
            # i.e. to seach for the value of the parameter "t" for which the
            # integral of the PDF is equal to "levels[j]" (and hence the
            # difference (integral-levels[j]) = 0)
            isocontourLevel.append(brentq(f_integral, 0., max_pdf))


        return isocontourLevel

