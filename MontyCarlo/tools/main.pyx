

__doc__ = """A place for tools to live before being moved to a proper module.
"""

__author__ = "Rui Campos"





from numpy cimport ndarray
cdef object remove_duplicates(ndarray x, ndarray Y):
    cdef ndarray u, c, dup
    u, c = np.unique(x, return_counts=True)
    dup = u[c > 1]
    
    cdef bint keep = True
    cdef list new_y = []
    cdef int i
    cdef double y
    
    for i, y in enumerate(Y):
        
        if x[i] in dup:
            if keep:
                new_y.append(y)
                keep = False
                continue
            else: continue
        
        if keep is False:
            keep = True
        new_y.append(y)
    
    return u, np.array(new_y)