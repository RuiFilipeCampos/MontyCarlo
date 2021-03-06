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

    ranges = (slice(0, 10, 10), #a1
              slice(0, 10, 10), #a2
              slice(0, 10, 50), #a3
              slice(0, 15, 10), #a4
              slice(0, 10, 10)) #a5


    def __init__(self, CS, Z):
        self.bounds = [(s.start, s.stop) for s in FormFactorWriter.ranges]

        _, CS = CS #ignoring interpolation flag
        
        self.xAxis, self.yAxis = getAxis(CS)
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

        self.basinhopping() #basinhopping seems to work pretty well
        self.__save__()
        
        


  
    def error(self, param):
        return error(*param, w = self.w, xAxis = self.xAxis, yAxis = self.yAxis)


    def basinhopping(self):
        from scipy.optimize import basinhopping

        res = basinhopping(self.error, 
                           tuple(5*[0.5]),
                           niter = 500,
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
        print(f"> Report has been saved to {path}.")
        print("")

        ### - writing report so I know what went down
        from time import asctime
        report = "###################################################\n"
        report = report + asctime() + "\n"
        report = report + self.res.__str__() + "\n"
        report = report + "RELATIVE ERROR: " + str(self.relative_error) + "\n"
        report = report + "R2: " + str(self.R2) + "\n"
        report = report + "###################################################\n"
        report = report + "\n\n\n"
        
        path = path + ".report.txt"

        with open(path, "a") as file:
            file.write(report)



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
                    Fk  = num/denom
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
        
        










    #OTHERS
    def plot(self):

        import matplotlib.pyplot as plt
        plt.figure()
        
        ax = plt.gca()
        ax.plot(self.xAxis, self.fit(self.xAxis))
        ax.scatter(self.xAxis, self.yAxis)


        ax.set_yscale('log')
        ax.set_xscale('log')

    def brute(self, Ns = 10):
        from scipy.optimize import brute
        res = brute(self.error, self.bounds, Ns=Ns)
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