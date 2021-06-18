#distutils: language = c++

from .mixmax.interface cimport mixmax_engine

ctypedef double (*func)


cdef class rng:
    cdef func gen

cdef mixmax_engine* genPTR