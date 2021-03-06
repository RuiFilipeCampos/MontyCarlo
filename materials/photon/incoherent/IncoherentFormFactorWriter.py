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


@njit
def func(x, a1, a2, a3, a4, a5):
    """x = 10**-10 * sin(theta/2)/lambda"""
    A = 1 + a1*x**2 + a2 * x**3 + a3 * x**4
    B = 1 + a4*x**2 + a5 * x**4
    return 1 - A/B**2

@njit
def error(a1, a2, a3, a4, a5, w = 1, xAxis = [], yAxis = []):
    #a1, a2, a3, a4, a5 = param
    #yAxisTrue = func(xAxis, *param)
    x = xAxis
    A = 1 + a1*x**2 + a2 * x**3 + a3 * x**4
    B = 1 + a4*x**2 + a5 * x**4
    yAxisTrue = 1 - A/B**2
    return sum(w*(yAxisTrue - yAxis)**2)


class IncoherentFormFactorWriter:
    ranges = (slice(0, 10, 10), 
              slice(0, 10, 10), 
              slice(0, 10, 50),  
              slice(0, 15, 10),  
              slice(0, 10, 10))
    def __init__(self, CS, Z):
        self.bounds = [(s.start, s.stop) for s in IncoherentFormFactorWriter.ranges]
        
        _, CS = CS
        
        self.xAxis, self.yAxis = getAxis(CS)
        self.Z = Z
        
        self.xAxis = self.xAxis*10**6
        
        self.yAxis = self.yAxis[5.00E-03 < self.xAxis]
        self.xAxis = self.xAxis[5.00E-03 < self.xAxis]

        self.xAxis = self.xAxis[9.9999E-1 > self.yAxis]
        self.yAxis = self.yAxis[9.9999E-1 > self.yAxis]
        
        
        

        assert len(self.xAxis) == len(self.yAxis)
        
        self.yAxis = self.yAxis/Z

        #w = sy.lambdify([x, y], w)
        w = lambda x, y: 1/y * 1/(1-y)
        self.w = w(self.xAxis, self.yAxis)
        self.basinhopping()
        self.__save__()

        self.plot()
        
        
    def __save__(self):
        a1, a2, a3, a4, a5 = self.param
        import dill as pickle

        
        
        import os
        directory = os.path.dirname(__file__)
        path = directory + r"\pickles" + "\\" + str(self.Z)
        with open(path, "wb") as file:
            pickle.dump(tuple(self.param), file)
            
        from time import asctime
        report = "###################################################\n"
        report = report + asctime() + "\n"
        report = report + self.res.__str__() + "\n"
        report = report + "RELATIVE ERROR: " + str(self.relative_error) + "\n"
        report = report + "R2: " + str(self.R2) + "\n"
        report = report + "\n\n\n"
        
        path = path + ".report.txt"

        
        with open(path, "a") as file:
            file.write(report)
        
        print("")
        print(f"> Incoherent form factor of element {self.Z} has been pickled.")
        print(f"> Report has been saved to {path}.")
        print("")
            
            

        
    def error(self, param):
        return error(*param, w = self.w, xAxis = self.xAxis, yAxis = self.yAxis)
    
    
    def basinhopping(self):
        from scipy.optimize import basinhopping
        res = basinhopping(self.error, 
                           tuple(5*[0.5]),
                           niter = 500,
                           minimizer_kwargs = dict(method = 'Nelder-Mead'))
        
        self.res = res
        print(self.res)
        self.param = self.res.x
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
    
    def brute(self, Ns = 10):
        from scipy.optimize import brute
        res = brute(self.error, self.bounds, Ns=Ns)
        self.param = res
        print(self.param)
        self.fit = lambda x: func(x, *self.param)
        self._measureGoodness()
    
    
        
        
    def _measureGoodness(self):
        self.R2 = sum((self.fit(self.xAxis) - self.yAxis)**2)
        
        self.relative_error = abs(self.fit(self.xAxis)- self.yAxis)/self.yAxis 
        self.relative_error = 100 * sum( self.relative_error ) / len(self.xAxis)
        #print("Relative Error", self.error)
        #print("R2:", self.sq)
        
        
    def plot(self):

        import matplotlib.pyplot as plt
        plt.figure()
        
        ax = plt.gca()
        ax.plot(self.xAxis, self.fit(self.xAxis))
        ax.scatter(self.xAxis, self.yAxis)


        ax.set_yscale('log')
        ax.set_xscale('log')
        
    
