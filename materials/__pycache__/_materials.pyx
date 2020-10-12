# cython: profile=True

from numpy import * #array, geomspace, flip, load, searchsorted

from . import electron as el



from .photon.photon import Photon
from ..settings import __montecarlo__

__materials__ = __montecarlo__/'materials'
directory = str(__materials__)



msg = "Can't modify initialized material!"

cdef class Material:
	def __init__(self, double Z, double A, object photon, object electron):
		self.Z = Z
		self.A = A
		self.photon   = photon
		self.electron = False
		self._name = "Compound"

	cpdef create(self, double density):
		"""Create a material by providing this function with a density value."""
		self.created = True
		self.photon.final_init(density)
		#self.electron.final_init(density)
		return self 

	def __mul__(first, second):

		if isinstance(first, float) and isinstance(second, Material):
			if second.created: raise RuntimeError(msg)
			photon = second.photon*first
			A      = second.A     *first
			return Material(second.Z, A, photon, second.electron)

		if isinstance(second, float) and isinstance(first, Material):
			if first.created: raise RuntimeError(msg)
			photon = first.photon*second
			A      = first.A     *second
			return Material(first.Z, A, photon, first.electron)

		if isinstance(first, int) and isinstance(second, Material):
			if second.created: raise RuntimeError(msg)
			photon = second.photon*first
			A      = second.A     *first
			return Material(second.Z, A, photon, second.electron)
			
		if isinstance(second, int) and isinstance(first, Material):
			if first.created: raise RuntimeError(msg)
			photon = first.photon*second
			A      = first.A     *second
			return Material(first.Z, A, photon, first.electron)

		print("> NotImplemented in materials.pyx")
		return NotImplemented

		
		


	def __add__(first, second):
		if isinstance(first, Material) and isinstance(second, Material):
			if first.created or second.created:
				raise RuntimeError(msg)

			A = first.A + second.A
		
			photon   = first.photon + second.photon
			electron = first.electron# + other.electron
			return Material(first.Z, A, photon, electron)

		return NotImplemented

	def __repr__(self):
		return "<" + self._name + ">"
























		


class Element(Material):
	def __init__(self, Z, full_init = False):
		self.Z = Z
		self._name = f"Element {self.Z}"
		
		#Evaluated Atomic Database Library
		EADL_path      = directory + "\\EADL\\" + str(Z) + ".txt"
		with open(EADL_path, "r") as file:
			text = file.readlines()
			text = [line.strip('\n') for line in text]
			line = text[0]
			self.A = float(line[13:24])
			
		
		#Evaluated Photon Database Library
		EPDL_path      = directory + "\\EPDL\\" + str(Z) + ".txt"
		self.EPDL_dict = self.get_bookmarked_text(EPDL_path)
		self.photon    = Photon(self, self.EPDL_dict)
		

		#Evaluated Electron Database Library
		EEDL_path      = directory + r"\\EEDL\\" + str(Z) + ".txt"
		self.EEDL_dict = self.get_bookmarked_text(EEDL_path)
		self.electron  = False
		#self.electron  = ph.Electron(self.EEDL_dict)
		
		
		

	def get_bookmarked_text(self, path):
		"""
		Reads EPICS file format and returns a dict with flags as keys
		and data as a list of strings.
		"""
		with open(path, "r") as file:
			text = file.readlines()
			text = [line.strip('\n') for line in text]

			bookmarks = [0]
			
			for n, line in enumerate(text):
				if line == "                                                                       1":
					bookmarks += [n + 1]
					
			#gather all bookmarked text into a dict
			bookmark_ = bookmarks[0:-1]
			_ookmarks = bookmarks[1:]
		
			bookmarked_text = {}
			for i, j in zip(bookmark_, _ookmarks):
				line1, line2 = text[i], text[i+1]

				#on line 1
				Yi = float(line1[7:9])    #particle identifier
				Yo = float(line1[10:12])  #secondary particle designator
				Iflag = float(line1[31])  #interpolation flag
				
				#on line 2
				C  = float(line2[0:2])    #reaction descriptor
				I  = float(line2[2:5])    #reaction property
				S  = float(line2[5:8])    #reaction modifier
				X1 = float(line2[22:32])  #subshell designator
				
				flags = (Yi, C, S, X1, Yo, I)

				flags = tuple(map(int, flags))
				bookmarked_text[flags] = (Iflag, text[i+2:j-1])
				
		return bookmarked_text


print("> Imported materials!")