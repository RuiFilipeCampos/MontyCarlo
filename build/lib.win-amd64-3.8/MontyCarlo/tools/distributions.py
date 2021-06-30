from scipy.interpolate import CubicSpline, RectBivariateSpline, BarycentricInterpolator, interp1d
from numpy import *
from numpy.random import rand
from scipy.integrate import *



class UnivariateDistribution:
    def __init__(self, xAxis = [], yAxis = []):
        xAxis, yAxis = map(array, (xAxis, yAxis))
        self.xAxis = xAxis
        self.yAxis = yAxis
        
        normalization = trapz(yAxis, x = xAxis)
        distribution = yAxis/normalization
        
        cum = cumtrapz(distribution, x = xAxis, initial = 0)
        self.cum    = interp1d(xAxis, cum)
        self.invCum = interp1d(cum, xAxis)
    


class Distribution: #RITA method(kind of)
    def __init__(self, xAxis = [], yAxisCol = [], EAxis = [], cumMoments = False):
        """
        Normalizes the distributions. Interpolates the probability distribution.
        Calculates and interpolates the inverse cumulutative function on the provided grid 'xAxis'.
        
        Usage:
        > x = arange(0, 1, .1)
        > energy = arange(0, 50)
        > probList = [prob1, prob2, prob3, ..., prob50]
        
        > dist = Distribution(xAxis = x, EAxis = energy, yAxisCol = probList)
        
        This will perform all operations and store them in memory.
        To get the distribution, one must update the state by giving an energy value:
        
        > dist.update_dist(1e6)
        
        This will perform aliasing according to a log interpolation in E and update
        the state, i.e. dist.prob and dist.invCum are available(and callable).
        """
        #note: yCol = Collection of y's
        
        assert len(EAxis) == len(yAxisCol)
        assert len(xAxis) == len(yAxisCol[0])
        
        #note: yAxis is an array of arrays with len(E)
        EAxis, xAxis, yAxisCol = map(array, [EAxis, xAxis, yAxisCol])

        self.EAxis    = EAxis
        self.xAxis    = xAxis
        self.yAxisCol = yAxisCol

        adhocCS        = [trapz(yAxis, x = xAxis) for yAxis in yAxisCol]
        self.probCol   = [yAxis/cs                for cs, yAxis in zip(adhocCS, yAxisCol)]
        self.probCol   = array(self.probCol)
        
        self.cumCol    = [cumtrapz(prob, x = xAxis, initial = 0) for prob in self.probCol]
        self.invCumCol = [interp1d(cum, xAxis)                   for cum  in self.cumCol]
        
        #can be improved, xAxis can be made denser if interpolation of prob is used
        if cumMoments is True: 
            self.probCol_m1 = []
            self.probCol_m2 = []
            
            for prob in self.probCol:
                self.probCol_m1.append(cumtrapz(prob*xAxis,    x = xAxis, initial = 0))
                self.probCol_m2.append(cumtrapz(prob*xAxis**2, x = xAxis, initial = 0))

            self.T1Col = [interp1d(xAxis, m1) for m1 in self.probCol_m1]
            self.T2Col = [interp1d(xAxis, m2) for m2 in self.probCol_m2]

        self.cumMoments = cumMoments
 
    def update_dist(self, E):
        """Update distributions according to energy."""
        
        k       = searchsorted(self.EAxis, E, side='left')
        Ek, Ek_ = self.EAxis[k-1], self.EAxis[k]
        pr      = (log(Ek_) - log(E))/(log(Ek_) - log(Ek))
        
        if rand() < pr: k -= 1

        self.E      = E
        self.prob   = self.probCol[k]
        self.invCum = self.invCumCol[k]
        
        if self.cumMoments is True:
            self.T1 = self.T1Col[k]
            self.T2 = self.T2Col[k]
