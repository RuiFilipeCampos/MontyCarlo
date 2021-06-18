# distutils: include_dirs = /mixmax
# distutils: language = c++

from .mixmax.interface cimport main
from .mixmax.interface cimport mixmax_engine






cdef mixmax_engine gen = mixmax_engine(0,0,0,123);

genPTR = &gen;


cdef double res
def py_main():
    return genPTR.get_next_float()


