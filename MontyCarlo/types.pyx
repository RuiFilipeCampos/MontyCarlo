"""
types.pyx
    A place to keep the frequently used data structures. Will eventually move everything here.
"""

import numpy as np

nan = np.nan

cdef struct double3:
    double x, y, z
