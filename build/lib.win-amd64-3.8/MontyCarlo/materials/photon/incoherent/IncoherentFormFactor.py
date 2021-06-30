# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 10:29:02 2020

@author: Rui Campos
"""

#External Imports
from numpy import *
from numba import *


#dealing with paths >.<
from ....settings import __montecarlo__
__incoherent__ = __montecarlo__/'materials'/'photon'/'incoherent'
directory = str(__incoherent__)

def getIFF(Z):
    import os
    path = directory + "\pickles\\" + str(Z)
    if not os.path.isfile(path):
        print("> Writing Incoherent Form Factor...")
        
        from ...database import EPDL
        IFFdata = EPDL[Z-1][(7, 93, 0, 0, 0, 942)]
        
        from . import IncoherentFormFactorWriter
        IncoherentFormFactorWriter.IncoherentFormFactorWriter(IFFdata, Z)
        
    import dill as pickle
    
    with open(path, 'rb') as file:
        param = pickle.load(file)

    a1, a2, a3, a4, a5 = param


    @njit
    def _eval(x):
        A = 1 + a1*x**2 + a2 * x**3 + a3 * x**4
        B = 1 + a4*x**2 + a5 * x**4
        return (1 - A/B**2)/Z

    return _eval


def composeIFF(formula):
    coefs = list(formula.values())
    fits  = [getIFF(Z) for Z in formula]
    
    
    Z = coefs[0]
    fit = fits[0]
    
    @njit
    def new_result(x):
        return Z*fit(x)
    
    for i in range(1, len(coefs)):
        Z = coefs[i]
        fit = fits[i]
        
        previous_result = new_result
        
        
        @njit
        def new_result(x):
            return previous_result(x) + Z*fit(x)

    return new_result

