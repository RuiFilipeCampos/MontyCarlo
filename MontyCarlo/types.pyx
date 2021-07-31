__doc__ = """
types.pyx
    A place to keep the frequently used data structures. (I will eventually move everything here).
"""

__author__ = "Rui Campos"

print("Importing .types")


import numpy as np

nan = np.nan

cdef struct double3:
    double x, y, z
    
    
class MAP(dict):
    """A dictionary with the following syntax sugar:
    
    > map = MAP()
    > map.x = 1
    
    is equivelent to:
    
    > map = MAP()
    > map['x'] = 1
    
    """
    
    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            raise AttributeError

    def __setattr__(self, key, value):
        self[key] = value

    def __delattr__(self, key):
        try:
            del self[key]
        except KeyError:
            raise AttributeError
