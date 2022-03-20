__doc__ = """Just a sketch for now.
"""
__author__ = "Rui Campos"



class Element(dict):
    
    def __init__(self, Z):
        self.__dict__[Z] = 1
        self.Z = Z
        
        
    def __mul__(self, other):
        #usual checks
        if not (isinstance(other, float) or isinstance(other, int)):
            return NotImplemented
        
        self.__dict__[Z] *= other
       
    def __add__(self, other):
        if isinstance(other, float) or isinstance(other, int):
            return NotImplemented
        
        return Molecule()
    
