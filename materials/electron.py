from scipy.interpolate import CubicSpline, RectBivariateSpline, BarycentricInterpolator, interp1d
from numpy import *
from numpy.random import rand
from scipy.integrate import *
#from mpmath import exp

def get_array(dic, key):
    """ find closest key and interpolates corresponding data"""
    try:
        return array(dic[key])
    except:
        keys = list(dic.keys())
        i = searchsorted(keys, key, side = "left")
        key0, keyf = keys[i-1], keys[i]
        return .5*(array(dic[key0]) + array(dic[keyf]))

def to_mu(theta):
    return (1-cos(theta))/2

def to_theta(mu):
    return arccos(1-2*mu)




class Distribution: #RITA method
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





    

class Electron:
    def __init__(self,
                 X, ang,
                 sp_args,
                 elastic_args, elasticDCS,
                 N,
                 density):

        self.brem    = Brem(X, ang, N)
        self.elastic = Elastic(elastic_args, elasticDCS, N, density)
        
        E = sp_args[0]
        sp = (sp_args[1] + sp_args[2])*density*1e3 ########## ALTERED DATA #$###############################
        self.sp = CubicSpline(E, sp)
        







class Elastic(Distribution):
    def __init__(self, elastic_args, elasticDCS, N, density):
        #E, sigma, sigma1, sigma2 = elastic_args #eV, cm^2, cm^2, cm^2
        #the grid of energies is approximetly logaritmic, with 15 points per decade

        #self.E = E

        
        self.N = N
        self.density = density  #g/cm^3



        self.E   = array(list(elasticDCS.keys()))
        E = self.E
        
        self._sigma, self._sigma1, self._sigma2 = [], [], []
        self._DCS = []
        for tupl in elasticDCS.values():

            self._sigma, self._sigma1, self._sigma2 = map(append,[ self._sigma,
                                                                   self._sigma1,
                                                                   self._sigma2], tupl[0:-1])
            self._DCS += [tupl[-1]]

        #self._sigma, self._sigma1, self._sigma2 = map(array, [self._sigma, self._sigma1, self._sigma2])
        self._DCS = array(self._DCS)

        
        
        #self.sigma  = CubicSpline(E, self._sigma )
        #self.sigma1 = CubicSpline(E, self._sigma1)
        #self.sigma2 = CubicSpline(E, self._sigma2)

        self.lambda_  = CubicSpline(E, 1/(self._sigma *self.N))
        self.lambda_1 = CubicSpline(E, 1/(self._sigma1*self.N))
        #self.lambda_2 = CubicSpline(E, 1/(self._sigma2*self.N))
            

        
        #self._DCS_ = array(list(elasticDCS.values())) #cm^2 - contains integrated cross sections (?)

        dN = 1/605
        self.mu = arange(0, 1+dN, dN)

        super().__init__(xAxis      = self.mu,
                         yAxisCol   = self._DCS,
                         EAxis      = self.E,
                         cumMoments = True)




    def __call__(self, E):
        """Return inverse of mean free path of hard elastic colision."""
        return 1/self.lambda_hard(E)

    def lambda_hard(self, E):
        """Returns the mean free path between two hard elastic collisions."""
        return max(self.lambda_(E), 0.2*self.lambda_1(E))
    
    def mu_c(self): 
        """Find cut off angle. (state dependent)"""
        
        return self.invCum(1 - self.lambda_(self.E)/self.lambda_hard(self.E))


    
    def F(self, mu_c, s): 
        """Simplified probability distribution for the polar displacement of a random hinge. Returns a and b."""      

        T1, T2 = self.T1(mu_c), self.T2(mu_c)
        
        #soft elastic first and second transport mean free paths
        lambda_1s = 1/(2*T1)
        lambda_2s = 1/(6*(T1-T2))

        #first and second moments of exact probability distribution
        mu  = .5*(1-exp(-s/lambda_1s))
        mu2 = mu - (1/6)*(1-exp(-s/lambda_2s))

        assert 1-2*mu > 0
        
        #paramaters of the simplified distribution
        #this distribution replicates the first and second moments of the exact distribution  
        b = (2*mu-3*mu2)/(1-2*mu)
        a = 1-2*mu+b
        
        return (a, b)

















class Brem:
    #k = load("K.npy")
    def __init__(self, _X, ang, N): #*10**-27

        k   = arange(0, 32)/32
        E   = array(list(_X.keys()))
        _X_ = [array(tup[1]) for tup in _X.values()]
        
        self.X = Distribution(xAxis    = k,
                              yAxisCol = _X_,
                              EAxis    = E)

        lambda_      = [1/(tup[0]*N*1e-27) for tup in _X.values()]
        self.lambda_ = CubicSpline(E, lambda_)
        self.ang     = ang


    def __call__(self, E):
        return 1/self.lambda_(E)
        


##class X(Distribution):
##    def __init__(self, xAxis = [], yAxisCol = [], EAxis = []):
##        super().__init__(xAxis    = xAxis,
##                         yAxisCol = yAxisCol,
##                         EAxis    = EAxis)







    
        
