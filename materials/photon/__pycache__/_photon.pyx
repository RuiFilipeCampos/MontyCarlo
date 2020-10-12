# cython: profile=True

#External Imports

#Internal Imports
from .coherent.coherent                   import Coherent
from .incoherent.incoherent               import Incoherent
from .photoelectric.photoelectric         import Photoelectric
from .pairproduction.pairproduction       import Pairproduction
from .tripletproduction.tripletproduction import Tripletproduction



cdef class Photon:
	def __init__(self, Material upper_dir, dict bookmarked_text):


		self.upper_dir = upper_dir

		cdef tuple CS, F, I, R, E, S
		#COHERENT SCATTERING
		CS = bookmarked_text[(7, 71, 0, 0, 0, 0)]   #integrated cross section
		F  = bookmarked_text[(7, 93, 0, 0, 0, 941)] #form factor
		I  = bookmarked_text[(7, 93, 0, 0, 0, 943)] #imaginary anomalous scattering
		R  = bookmarked_text[(7, 93, 0, 0, 0, 944)] #real anomalous scattering
		E  = bookmarked_text[(7, 71, 0, 0, 7, 10)]  #average energy of the scattered photon
		
		self.coherent = Coherent(upper_dir, CS, F, R, I, E)

		#INCOHERENT SCATTERING
		CS = bookmarked_text[(7, 72, 0, 0, 0, 0)]   #integrated cross section
		S  = bookmarked_text[(7, 93, 0, 0, 0, 942)] #incoherent form factor
		
		self.incoherent = Incoherent(upper_dir, CS, S)



		#PHOTOELECTRIC
		CS = bookmarked_text[(7, 73, 0, 0, 0, 0)]
		self.photoelectric = Photoelectric(upper_dir, CS)



		#PAIR PRODUCTION
		CS = bookmarked_text[(7, 74, 0, 0, 0, 0)]
		self.pairproduction = Pairproduction(upper_dir, CS)



		#TRIPLET PRODUCTION
		CS = bookmarked_text[(7, 75, 0, 0, 0, 0)]
		self.tripletproduction = Tripletproduction(upper_dir, CS)
							 
		
	def __add__(self, Photon other):
		coherent   = self.coherent   + other.coherent
		incoherent = self.incoherent + other.incoherent

		return CompoundPhoton(self.upper_dir, 
							  coherent, 
							  incoherent, 
							  self.photoelectric, 
							  self.pairproduction,
			           		  self.tripletproduction)
		
	def __mul__(self, float other):
		coherent   = self.coherent*other
		incoherent = self.incoherent*other

		return CompoundPhoton(self.upper_dir, 
							  coherent, 
							  incoherent, 
							  self.photoelectric, 
							  self.pairproduction,
			           		  self.tripletproduction)


	cpdef final_init(self, double density):
		self.coherent.final_init(density)
		self.incoherent.final_init(density)

	def __repr__(self):
		return "<" + self.upper_dir._name + ": /photon>"









cdef class CompoundPhoton(Photon):
	def __init__(self, Material upper_dir, 
			           object coherent,
			           object incoherent,
			           object photoelectric,
			           object pairproduction,
			           object tripletproduction):

		self.upper_dir         = upper_dir
		self.coherent          = coherent
		self.incoherent        = incoherent
		self.photoelectric     = photoelectric
		self.pairproduction    = pairproduction
		self.tripletproduction = tripletproduction




