# cython: profile=True

"""

MAIN SOURCE OF INFO:
    PENELOPE 11  - https://drive.google.com/file/d/1aIfsVbZUbIwjDzlvOwWBA_lN2XTu2VGd/
    PENELOP 2018 - https://drive.google.com/file/d/1rb_wKkICOyL5UMuG4chuRxBQHuqR_8q1/
"""


#from numpy import * #array, geomspace, flip, load, searchsorted

#from . import electron as el
from .photon.photon import Photon
from .pyRelax import Atom

#msg = "Can't modify initialized material!"

class Molecule:
    """
    Container class for relaxation models and other information.
    """
    def __init__(self, formula):
        self.formula = {Atom(Z):formula[Z] for Z in formula}

class Material:
    """
    Objects of this class are given to Vol.
    """
    def __init__(self, formula, density):
        
        
        print("")
        print("_____________________________________")
        print("INTITIALIZING MATERIAL")
        print("FORMULA: ", formula)
        print("DENSITY:", density)


        self.formula = formula
        self.density = density
        
        print("COMPILING PHOTON DATA...")
        print("")

        self.photon   = Photon(formula, density)
        
        print("")
        print("CONSTRUCTING RELAXATION MODELS...")
        
        self.molecule = Molecule(formula)
        self.photon.molecule = self.molecule
        
        print("")
        print("WARNING: Electron reference is set to None")
        
        self.electron = False
        
        print("")
        print("DONE")
        print("----------------------------------")
        