"""
Objective:
-> create a callable that receives x**2 and outputs F(x**2)**2.

should be low level, and it should be possible easily save
and read it.
"""

#external imports
from numba import njit
from numpy import *

#internal imports
from ....tools.data import getAxis


#dealing with paths >.<
from ....settings import __montecarlo__
__coherent__ = __montecarlo__/'materials'/'photon'/'coherent'
directory = str(__coherent__)



#This is the function that will be fitted.
@njit
def func(x, a1, a2, a3, a4, a5):
    """x = 10**-10 * sin(theta/2)/lambda"""
    A = 1 + a1*x**2 + a2 * x**3 + a3 * x**4
    B = 1 + a4*x**2 + a5 * x**4
    return A/B**2


#This is the error to be minimized (weighted squared sum of residuals)
@njit
def error(a1, a2, a3, a4, a5, w = 1, xAxis = [], yAxis = []):
    #a1, a2, a3, a4, a5 = param
    #yAxisTrue = func(xAxis, *param)
    x = xAxis
    A = 1 + a1*x**2 + a2 * x**3 + a3 * x**4
    B = 1 + a4*x**2 + a5 * x**4
    yAxisTrue = A/B**2
    return sum(w*(yAxisTrue - yAxis)**2)


class FormFactorWriter:
    #I've used several optimizers from scipy, some required
    #slice objects, others just the bounds.
    """
	DATA BASE DOC: 
        https://www-nds.iaea.org/epics/DOCUMENTS/ENDL2002.pdf
    
    REFERENCES:
        -----------------------------------------------------------------------------------------------
        ??? :: Atomic Form Factors, Incoherent Scattering Functions,and Photon Scattering Cross Sections
        https://drive.google.com/file/d/1hbhDTCn1NGYZIB31K-OW21RsbFEPO9ta/
        ->  by: Hubble et al.
        -----------------------------------------------------------------------------------------------
        photonCS1994 :: Analytical cross sections for Monte Carlo simulation of photon transport.
        https://drive.google.com/file/d/1Rt2DqkhwJINQC1S469adqn5Whz9110ep/view?usp=sharing
        ->  savat et. al
        -----------------------------------------------------------------------------------------------
    
	"""

    ranges = (slice(0, 10, 10), #a1
              slice(0, 10, 10), #a2
              slice(0, 10, 50), #a3
              slice(0, 15, 10), #a4
              slice(0, 10, 10)) #a5


    def __init__(self, DATA, Z):
        self.bounds = [(s.start, s.stop) for s in FormFactorWriter.ranges]

        #_, CS = CS #ignoring interpolation flag
        
        self.xAxis, self.yAxis = DATA.X, DATA.Y
        self.Z = Z
        
        self.xAxis = self.xAxis*10**6
        self.xAxis = self.xAxis[self.yAxis < 0.99*Z]
        self.yAxis = self.yAxis[self.yAxis < 0.99*Z]

        #the is not flexible enough for all data
        #the arrays have been cut in accordance to the
        #paper

        if self.Z < 11:
            self.xAxis = self.xAxis[1e-5 < self.yAxis]
            self.yAxis = self.yAxis[1e-5 < self.yAxis]

        else:
            self.xAxis = self.xAxis[2 < self.yAxis]
            self.yAxis = self.yAxis[2 < self.yAxis]


        assert len(self.xAxis) == len(self.yAxis)
        
        self.yAxis = self.yAxis/Z

        w = lambda x, y: 1/y**2 * 1/(1-y)
        self.w = w(self.xAxis, self.yAxis)

        #self.basinhopping() #basinhopping seems to work pretty well
        #self._measureGoodness()
        #if self.relative_error >= 1:
            #self.bounds = [(a-1, a+1) for a in self.param]
            #self.brute(Ns = 30)
        #self.brute(Ns = 30)
        #self.simplex()
        
        
        self.brute()
        self.__save__()
        
        


  
    def error(self, param):
        return error(*param, w = self.w, xAxis = self.xAxis, yAxis = self.yAxis)
    
    
    def simplex(self):
        """
        NUMERICAL RECIPES: https://www.google.com/url?q=https://drive.google.com/file/d/1K6tg_um1G-5X-5hdgwKb8vKyBD5oz4vX/view?usp%3Dsharing&sa=D&ust=1606856195190000&usg=AFQjCNH2Sc8A5HVZrZb72T-DrPM2-1kRtQ
        """
        
        from scipy.optimize import fmin as simplex
        
        self.param = simplex(self.error, 5*[5], xtol=1e-50, ftol=1e-50, maxfun=1e6)
        print(self.param)
        
        self.fit = lambda x: func(x, *self.param)
        self._measureGoodness()

    def basinhopping(self):
        from scipy.optimize import basinhopping

        res = basinhopping(self.error, 
                           tuple(5*[5]),
                           niter = 3*1000,
                           minimizer_kwargs = dict(method = 'Nelder-Mead'))
        
        self.res   = res        #saving report
        self.param = self.res.x #(a1, a2, a3, a4, a5)

        print(self.res)

        self.fit = lambda x: func(x, *self.param)
        self._measureGoodness()

        
    def _measureGoodness(self):
        self.R2 = sum((self.fit(self.xAxis) - self.yAxis)**2)
        
        self.relative_error = abs(self.fit(self.xAxis)- self.yAxis)/self.yAxis 
        self.relative_error = 100 * sum( self.relative_error ) / len(self.xAxis)
        #print("Relative Error", self.error)
        #print("R2:", self.sq)
        
        

        
    
    def __save__(self):
        
        Z = self.Z
        
        
        
        ### - creating jitted function
        


        ### - saving jitted function
        import dill as pickle
        import os
        
        path = directory + r"\pickles" + "\\" + str(self.Z)

        with open(path, "wb") as file:
            #to_save = (F2, self.xAxis**2) #not so sure about the sq
            pickle.dump(tuple(self.param), file)



        
        print("")
        print(f"> Coherent form factor of element {self.Z} has been pickled.")
        # print(f"> Report has been saved to {path}.")
        # print("")

        ### - writing report so I know what went down
        # from time import asctime
        # report = "###################################################\n"
        # report = report + asctime() + "\n"
        # report = report + self.res.__str__() + "\n"
        # report = report + "RELATIVE ERROR: " + str(self.relative_error) + "\n"
        # report = report + "R2: " + str(self.R2) + "\n"
        # report = report + "###################################################\n"
        # report = report + "\n\n\n"
        
        # path = path + ".report.txt"

        # with open(path, "a") as file:
        #     file.write(report)
            

        

        a1, a2, a3, a4, a5 = self.param

        from numba import njit
        if self.Z > 10:
            @njit
            def F2(x2):
                x = sqrt(x2) #changing variables

                #calculating f(Z, x)
                A = 1 + a1*x**2 + a2 * x**3 + a3 * x**4
                B = 1 + a4*x**2 + a5 * x**4
                f = Z*A/B**2


                if f < 2.:
                    alpha = 1/137.03599908421  #fine structure constant
                    a = (Z - 5/6)*alpha

                    Q = x/20.6074/2/a          #changing variables to suit formula in paper

                    gamma = (1 - a**2)**.5

                    #calculating Fk
                    num = sin(2*gamma*arctan(Q))
                    den = gamma*Q*(1 + Q**2)**gamma
                    Fk  = num/den
                    return max(f, Fk)**2
                else: return f**2
        else:
            #FF of elements below 11 are not defined by a piecewise func (see paper)
            @njit
            def F2(x2):
                x = sqrt(x2)
                A = 1 + a1*x**2 + a2 * x**3 + a3 * x**4
                B = 1 + a4*x**2 + a5 * x**4
                return (Z*A/B**2)**2
        
        self.F2 = F2
        
        path = directory + r"\pickles" + "\\" + str(self.Z)
        
        self.plot(path)




    def plot(self, path):

        import matplotlib.pyplot as plt
        plt.figure(figsize=(11, 11))
        
        ax = plt.gca()

        x0, xf = self.xAxis[0], self.xAxis[-1]
        xAxis = arange(x0, xf, (xf-x0)/(4*len(self.xAxis)))
 
        
        
        self.Y = [self.Z**-1 * (self.F2(x**2))**.5 for x in xAxis]
        
        ax.plot(xAxis, self.Y, label="fit")
        ax.scatter(self.xAxis, self.yAxis, s=3, color="r", label="data from EPDL")
        

        ax.set_yscale('log')
        ax.set_xscale('log')
        
        plt.title("Fit of Form Factor for element Z = " + str(self.Z))
        plt.xlabel("x**2 (units??)")
        plt.ylabel("|F(x)|**2 (units??)")
        
        ax.legend()
        ax.text(0.01, 0.01, 
                f"""
                Relative Error: {self.relative_error}
                R2: {self.R2}
                """, transform=ax.transAxes
                )
        
        
        plt.savefig(path)




    # #OTHERS
    # def plot(self):

    #     import matplotlib.pyplot as plt
    #     plt.figure()
        
    #     ax = plt.gca()
    #     ax.plot(self.xAxis, self.fit(self.xAxis))
    #     ax.scatter(self.xAxis, self.yAxis)


    #     ax.set_yscale('log')
    #     ax.set_xscale('log')

    def brute(self, Ns = 35):
        from scipy.optimize import brute, fmin
        res = brute(self.error, self.bounds, Ns=Ns, finish = fmin)
        self.param = res
        print(self.param)
        self.fit = lambda x: func(x, *self.param)
        self._measureGoodness()


    def shgo(self):
        from scipy.optimize import shgo
        
        res = shgo(self.error, 
                   self.bounds, 
                   constraints=None, 
                   n=60, 
                   sampling_method='sobol',
                   iters=30)
        
        self.param = res.x
        self.fit = lambda x: func(x, *self.param)
        self._measureGoodness()