#distutils: language = c++
#distutils: sources = cppRelax.cpp

from .._random.random cimport mixmax_engine
from libcpp.deque cimport deque
from libcpp.vector cimport vector
from libcpp.string cimport string



cdef extern from "cppRelax.cpp":
	pass

cdef extern from "cppRelax.h" namespace "model":
	cdef struct PARTICLES:
		vector[double] PHOTONS; 
		vector[double] ELECTRONS;

ctypedef void (*f)


# Declare the class with cdef
cdef extern from "cppRelax.h" namespace "model":

		
	cdef struct Transition:
		double p;
		Transition* t;
		double E;
		Shell* j;
		Shell* k;
		void (*perform)(void *_this, PARTICLES *particles);
        
        
	void setADRESS_RAD(Transition * tr);
	void setADRESS_NONRAD(Transition * tr);


		
		
		
	cdef cppclass Shell:
		
		Shell() except +


		Shell(double* frac, bint dontSIMULATE) except +
		string list_transitions() except +
		Transition* sample_transition() except +
		double *frac
		double *original_frac;
		bint dontSIMULATE;
		bint LAST;
		int Nt;

		Transition *trSTART;

		void  setTRANSITIONS(Transition* arr, int N)  except +;

	cdef cppclass Atom:
	
		Atom() except +

		
		
		Atom(Shell* arr, int Nsh, int Nt) except +;
		void run(int shell_index, PARTICLES* particles,  mixmax_engine *genPTR)	 except +

		Shell* fetchFI(int index)  except +
		
		string list_array() except +;



