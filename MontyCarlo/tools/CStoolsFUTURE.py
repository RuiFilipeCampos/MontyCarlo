#External Imports
from numpy import *
#from scipy.interpolate import interpl1d

#Internal Imports
from .data import getAxis
from .others import searchsorted
#from .interpol import LinearInterpolation
from .interpolation import LinearInterpolation
from .others import searchsorted

from ..materials import database


def getData(Z):
    x = database.EPDL[Z][(7, 71, 0, 0, 0, 0)]
    Iflag, data = x[0], x[1]
    xAxis, yAxis = getAxis(data) #verificar unidades
    return Iflag, xAxis, yAxis
    
def getInterpol(Z):
    iflag, xAxis, yAxis = getData(Z) #verificar unidades
    cythonInterpolator = LinearInterpolation(xAxis, yAxis)
    return cythonInterpolator._eval

def getCompoundCS(formula, density):
    w = formula.keys()
    for Z in w:
        formula[Z] = getInterpol(Z)

    def CS(E):
        return sum(Z*formula[Z](E) for Z in formula)

    return CS



