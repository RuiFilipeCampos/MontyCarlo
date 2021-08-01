

__doc__ = """A place for tools to live before being moved to a proper module.
"""

__author__ = "Rui Campos"


import numpy as np
from collections import deque

cdef object remove_duplicates(ndarray X, ndarray Y):
    """Removes duplicates from the (x, y) tuple of arrays.
    """

    if not (len(X) == len(Y)):
        raise RuntimeError(f"Arrays have different shapes.")

    new_x, new_y = deque(), deque()

    duplicates = deque([None])
    
    cdef int i
    cdef int N = len(X)

    for i in range(N - 1):
        if X[i+1] - X[i] == 0: # is duplicate
            if X[i] == duplicates[-1]: # and has already been cached
                continue

            duplicates.append(X[i])
            new_x.append(X[i])
            new_y.append(Y[i])
            continue

        if X[i] == duplicates[-1]:
            continue

        new_x.append(X[i])
        new_y.append(Y[i])

    if not (X[-1] == duplicates[-1]):
        new_x.append(X[i])
        new_y.append(Y[i])

    return (np.array(new_x), np.array(new_y))



class python_hooks:
    """Namespace for collecting python wrappers.
    """

    @staticmethod
    def remove_duplicates(x, y):
        return remove_duplicates(x, y)