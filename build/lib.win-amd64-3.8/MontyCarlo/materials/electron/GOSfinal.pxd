from ..cppRelaxAPI cimport Atom as rAtom
from ..materials cimport Atom as _Atom


from ..cppRelaxAPI cimport PARTICLES
from ..materials cimport Molecule
from ..materials cimport Atom


from ..._random.interface cimport mixmax_engine


cdef class gosMolecule:
    cdef double number_density
    cdef gosAtom[::1] gosATOMS
    cdef public double[:, ::1] totalCS;
    cdef void sample(gosMolecule self, mixmax_engine* genPTR, int index, double E, PARTICLES *particles)
    cdef int Nat
    
    
    
cdef class gosAtom:
    cdef gosShell[::1] gosSHELLS
    cdef double[:, ::1] totalCS
    cdef double Z
    cdef int Nat
    cdef public int Nsh
    cdef rAtom rATOM
    cdef _Atom ATOMptr

    cdef void sample(gosAtom self, mixmax_engine* genPTR, int index, double E, PARTICLES *particles)

cimport numpy as cnp

cdef class gosShell:
    cdef public cnp.ndarray p1, p2, p3
    cdef double fk

    cdef int Nsh
    cdef rAtom rATOM
    cdef double Wk;
    cdef _Atom ATOMptr
    cdef:
        double[::1, :] totalCS;
        double[::1, :] pCLOSE 
        double[::1, :] pLFAR  
        double[::1, :] pTFAR  
    cdef void sample(self, mixmax_engine* genPTR, int index, double E, PARTICLES *particles)
