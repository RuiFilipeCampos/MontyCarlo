# -*- coding: utf-8 -*-
"""
Created on Sat Nov 28 16:58:34 2020

@author: Rui Campos
"""

#External Imports
from numpy import *
from numba import *


#dealing with paths >.<
from ....settings import __montecarlo__
__coherent__ = __montecarlo__/'materials'/'photon'/'coherent'
directory = str(__coherent__)





def getFF(Z):
    import os
        
    path = directory + "\pickles\\" + str(Z)
        
    if not os.path.isfile(path):
            
        print("> *** photon/coherent: Writing form factor...")
            
        from ...database import EPDL
        from .FormFactorWriter import FormFactorWriter
        
        FFdata = EPDL[Z-1][(7, 93, 0, 0, 0, 941)]
        
        #xAxis, yAxis = getAxis(FFdata)
        
        FormFactorWriter(FFdata, Z)
        
            
    
    
    
    import dill as pickle

    with open(path, 'rb') as file:
        param = pickle.load(file)

    a1, a2, a3, a4, a5 = param


    #depending on Z, the fit is calculated differently
    if Z < 11:
        @njit
        def _eval(x2):
            x = sqrt(x2)
            A = 1 + a1*x**2 + a2 * x**3 + a3 * x**4
            B = 1 + a4*x**2 + a5 * x**4
            f = Z*A/B**2
            return f**2

    else:

        @njit
        def _eval(x2):
            x = sqrt(x2)
            A = 1 + a1*x**2 + a2 * x**3 + a3 * x**4
            B = 1 + a4*x**2 + a5 * x**4
            f = Z*A/B**2

            if f > 2: return f**2
			
            alpha = 1/137.03599908421  #fine structure constant
            a = (Z - 5/6)*alpha

            Q = x/20.6074/2/a #changing variables to suit formula in paper

            gamma = (1 - a**2)**.5

            #calculating Fk
            num = sin(2*gamma*arctan(Q))
            den = gamma*Q*(1 + Q**2)**gamma
            Fk  = num/den
            return max(f, Fk)**2
        
    return _eval

def composeFF(formula):
    
    fits = [getFF(Z) for Z in formula]
    coefs = list(formula.keys())
    N = len(fits)

    #@njit
    def FF(x):
        return sum(fits[i](x)*coefs[i] for i in range(N))
    return FF
        

