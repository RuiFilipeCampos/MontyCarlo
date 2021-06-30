# cython: profile=False



cdef int _sortedListDOUBLE(list L, double l, int start, int end):
    if start > end:
        return -1
    
    cdef int mid_point = int (start + (end - start)/2)
    
    cdef int i
    
    if l < L[mid_point]:
        i =  _sortedListDOUBLE(L, l, start, mid_point-1)
        if i == -1: return mid_point - 1
        else:       return i
    
    if L[mid_point] < l:
        i = _sortedListDOUBLE(L, l, mid_point+1, end)
        if i == -1: return mid_point
        else:       return i
        
    else:
        return mid_point
    


def sortedListDOUBLE(L, a, start, end):
    return _sortedListDOUBLE(L, a, start, end)





cimport numpy as cnp





cimport cython

@cython.boundscheck(False)  # SET TO TRUE
@cython.wraparound(False) 
@cython.initializedcheck(False)
@cython.cdivision(True)
cdef int _sortedArrayDOUBLE(double[:] arr, double element, int start, int end):
    cdef int mid
    
    while start <= end:
        #print(f"Subarray: {list(arr[start:end+1])}")
        mid = start + (end - start)//2
        
        if element == arr[mid]:
            return mid
        
        if element < arr[mid]:
            end = mid - 1
        else:
            start = mid + 1
    return end









# #cimport cython
# #@cython.boundscheck(False)  # Deactivate bounds checking
# #@cython.wraparound(False)   # Deactivate negative indexing.
# cdef int _sortedArrayDOUBLE(double [:] L, double l, int start, int end):
#     if start > end:
#         return -1
    
#     cdef int mid_point = int(start + (end - start)/2)
    
#     cdef int i
    
#     if l < L[mid_point]:
#         i =  _sortedArrayDOUBLE(L, l, start, mid_point-1)
#         if i == -1: return mid_point - 1
#         else:       return i
    
#     if L[mid_point] < l:
#         i = _sortedArrayDOUBLE(L, l, mid_point+1, end)
#         if i == -1: return mid_point
#         else:       return i
        
#     else:
#         return mid_point
    


def sortedArrayDOUBLE(double [:] L, double l):
    return _sortedArrayDOUBLE(L, l, 0, len(L))