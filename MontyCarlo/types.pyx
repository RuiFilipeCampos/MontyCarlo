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
