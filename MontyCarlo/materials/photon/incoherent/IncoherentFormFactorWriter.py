"""
Objective:
-> create a callable that receives x**2 and outputs F(x**2)**2. REESCREVER

should be low level, and it should be possible easily save
and read it.


REFERENCES:
        -----------------------------------------------------------------------------------------------
        photonCS1994 :: Analytical cross sections for Monte Carlo simulation of photon transport.
        https://drive.google.com/file/d/1Rt2DqkhwJINQC1S469adqn5Whz9110ep/view?usp=sharing
        ->  savat et. al
        -----------------------------------------------------------------------------------------------
"""

#external imports
from numba import njit
from numpy import *

#internal imports
from ....tools.data import getAxis


























































@njit
def func(x, a1, a2, a3, a4, a5):
    """
    Fit function.
    
    x = 10**-10 * sin(theta/2)/lambda
    
    Equation 29 @ photonCS1994
    https://drive.google.com/file/d/1Rt2DqkhwJINQC1S469adqn5Whz9110ep
    
    """
    A = 1 + a1*x**2 + a2 * x**3 + a3 * x**4
    B = 1 + a4*x**2 + a5 * x**4
    return 1 - A/B**2

@njit
def error(a1, a2, a3, a4, a5, w = 1, xAxis = [], yAxis = []):
    """
    Error to be minimized.
    

    Equation 30 @ photonCS1994
    https://drive.google.com/file/d/1Rt2DqkhwJINQC1S469adqn5Whz9110ep
    
    """
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
        
        """
        ??? @ photonCS1994
        https://drive.google.com/file/d/1Rt2DqkhwJINQC1S469adqn5Whz9110ep

        """
        
        
        
        self.bounds = [(s.start, s.stop) for s in IncoherentFormFactorWriter.ranges]
        
        _, CS = CS
        
        self.xAxis, self.yAxis = getAxis(CS)
        self.Z = Z
        
        self.xAxis = self.xAxis*10**6
        
        
        # 
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
        
        path = directory + r"\pickles" + "\\" + str(self.Z) + "report\\" #+ str(self.Z) + ".txt"
        if not os.path.exists(path):
            os.makedirs(path)
        
        path = directory + r"\pickles" + "\\" + str(self.Z) + "report\\" + str(self.Z) + ".txt"
        with open(path, "a") as file:
            file.write(report)
        
        path = directory + r"\pickles" + "\\" + str(self.Z) + "graph"
        self.plot(path)
        
        print("")
        print(f"> Incoherent form factor of element {self.Z} has been pickled.")
        print(f"> Report has been saved to {path}.")
        print("")
            
            

        
    def error(self, param):
        return error(*param, w = self.w, xAxis = self.xAxis, yAxis = self.yAxis)
    
    
    def basinhopping(self):
        """
        REFERENCES:
            WIKI :: https://en.wikipedia.org/wiki/Basin-hopping
            SCIPY DOCS :: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.basinhopping.html
            -----------------------------------------------------------------------------------------------
            Olson-2012
            Basin Hopping as a General and Versatile Optimization Framework for the Characterization 
            of Biological Macromolecules, Advances in Artificial Intelligence
            https://drive.google.com/file/d/1Yqop2eV52wUH-L8B7AwlGMNbY5WGawmp/
            -----------------------------------------------------------------------------------------------
        
        DESCRIPTION OF BASINHOPPING:
            Basin-hopping is a two-phase method that combines a global stepping algorithm with local 
            minimization at each step. Designed to mimic the natural process of energy minimization 
            of clusters of atoms, it works well for similar problems with “funnel-like, but rugged” 
            energy landscapes
            
        """
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
        """
        Equation 32 @ photonCS1994
        https://drive.google.com/file/d/1Rt2DqkhwJINQC1S469adqn5Whz9110ep
        
        """
        
        
        
        self.R2 = sum((self.fit(self.xAxis) - self.yAxis)**2)
        
        self.relative_error = abs(self.fit(self.xAxis)- self.yAxis)/self.yAxis 
        self.relative_error = 100 * sum( self.relative_error ) / len(self.xAxis)
        
        print("___________________________")
        print("Relative Error", self.relative_error)
        print("R2:", self.R2)
        print("¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯")
        
        
    def plot(self, path):

        import matplotlib.pyplot as plt
        plt.figure(figsize=(8, 8))
        
        ax = plt.gca()

        x0, xf = self.xAxis[0], self.xAxis[-1]
        xAxis = arange(x0, xf, (xf-x0)/(4*len(self.xAxis)))
        ax.plot(xAxis, self.fit(xAxis), label="fit")
        ax.scatter(self.xAxis, self.yAxis, s=3, color="r", label="data from EPDL")
        

        ax.set_yscale('log')
        ax.set_xscale('log')
        
        plt.title("Fit of Incoherent Form Factor for element Z = " + str(self.Z))
        plt.xlabel("x (units??)")
        plt.ylabel("S(x) (units??)")
        
        ax.legend()
        ax.text(0.3, 0.01, 
                f"""
                Relative Error: {self.relative_error}
                R2: {self.R2}
                """, transform=ax.transAxes
                )
        
        
        plt.savefig(path)
        
    
