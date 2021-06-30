# import numpy as np
# import pyximport
# pyximport.install(reload_support= True,
#                   inplace       = True,
#                   setup_args = {'include_dirs':np.get_include()})


print("____INIT_____")


#### PREPARING ENERGY AXIS AND GRID COMMON TO ALL MODULES
from .settings import __montecarlo__
import numpy as np 


__path__ = __montecarlo__/'materials'/'electron'


PATH = __path__/'elastic'
PATH = str(PATH)


#mu = np.load(path + "/muGRID.npy")

LEeax = np.load(PATH + "/LEeax.npy")
HEeax = np.load(PATH + "/HEeax.npy")
eax =  np.append(LEeax, HEeax[1:])


eax = np.ascontiguousarray(eax) # for python data processing
EAX = eax  ### array  defined in pxd, cdef double EAX[2695]

 



cdef int get_exp(double x):
    cdef int exp;
    frexp(x, &exp);
    return exp;




hashed =  np.array([get_exp(E) for E in eax], dtype = int)
indexes = np.arange(0, len(hashed))
Imax = int(max(hashed))



lims = [np.array([0, 0, 0], dtype = int)]


cdef int i 
for i in range(Imax + 1): #every possible value of the hash, index = hash
    selected = indexes[hashed == i]
    n = len(selected)
    if n == 0: # either out of bounds or no values in this range
        #if no values in this range, interpolate using last interval
        n_last = lims[-1][2]
        if n_last == 0: #out of bounds
            lims.append(np.array([0, 0, 0], dtype = int))
            continue
        i_last = lims[-1][1]
        lims.append(np.array([i_last, i_last, 1], dtype = int))
        continue

    lims.append(np.array([selected[0], selected[-1] , n], dtype = int))

LIMS = np.array(lims[1:], dtype = int) ### memory view defined in pxd, cdef double[::1] EAX



