

cimport numpy as cnp



def searchsorted(arr, element, x0, xf):
    return _searchsorted(arr, element, x0, xf)

cdef int _searchsorted(cnp.ndarray[cnp.float_t, ndim=1] arr, 
                       double element, 
                       int x0, 
                       int xf):
    
    cdef int k = (x0 + xf)/2

    if k == 0: return k
    
    if arr[k - 1] < element <= arr[k]:
        return k

    if arr[k] < element:
        #discard left side
        return _searchsorted(arr, element, k, xf)

    if arr[k] > element:
        #discard right side
        return _searchsorted(arr, element, x0, k)