# distutils: language = c++

#cdef double[::1] EAX


cdef extern from "<math.h>" nogil:
    double frexp(double x, int* exponent)


#cdef int get_exp(double x)





from .electron.main cimport Electron
from .positron.main cimport Positron

from .photon.photon cimport Photon


from ..external.mixmax_interface cimport mixmax_engine
from .cppRelaxAPI cimport Atom as rAtom
from .cppRelaxAPI cimport Shell as rShell


from libcpp.vector cimport vector
from .cppRelaxAPI cimport PARTICLES
from .pyRelax cimport Atom as crAtom #cython relaxation atom - Atom instance from cython class that wrapps the relaxation model atom instance
from ..tools.interpol1 cimport CSa


cdef class Shell:
    cdef double Z
    cdef double binding_energy
    cdef int index
    cdef rAtom rATOM
    cdef double[:] PHELa, PHELb
    #cdef void ionize(self, mixmax_engine *genPTR, PARTICLES* particles)
    cdef void ionize(self, mixmax_engine *genPTR, PARTICLES* particles, double E)
    cdef double sample_compton_profile(self, mixmax_engine *genPTR, double pz_max)
    cdef double[::1] cumul
    cdef double[:, ::1] cCUMUL, cINVCUMUL
    
cdef class Atom(crAtom):
    cdef double[:, ::1] ALIAS
    cdef double CUT_OFF
   # cdef int Z
    cdef double[::1] PHELa, PHELb
    cdef Shell[::1] arrSHELLS
    cdef CSa S
    
    cdef void PHELchoose(Atom self,int index, double E, mixmax_engine* genPTR, PARTICLES* particles)
    cdef void ionize(self, int shell_index, mixmax_engine* genPTR, PARTICLES* particles)


cdef class Molecule:
    cdef Atom choose_atom(self, mixmax_engine *genPTR)
    cdef double[::1, :] atomALIAS
    
    
    
    
    
    cdef void ionize_particular(self, mixmax_engine *genPTR, PARTICLES* particles, int atom_index, int shell_index)
    cdef int Nat
    #cdef Atom* arrATOMS
    cdef Atom[::1] arrATOMS
    cdef double[::1] arrNi
    cdef double[::1] PHELa, PHELb
    cdef void PHELionize(Molecule self, int index, double E, mixmax_engine *genPTR, PARTICLES* particles)
    cdef Atom get(self, double Z)
    
    cdef double[:, ::1] ALIAS
    cdef double Nsh
    cdef void ionize(self, mixmax_engine *genPTR, PARTICLES* particles)
    

cdef class Material:
    cdef Electron electron
    cdef Photon photon
    cdef Positron positron
    cdef public object formula
    cdef double C1, C2, Wcc, Wcr
    cdef public double density
    cdef Molecule molecule
    cdef public object name
    cdef public dict python
    cdef public double Am, N

        