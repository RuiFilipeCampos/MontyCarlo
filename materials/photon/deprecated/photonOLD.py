# cython: profile=True

#External Imports

#Internal Imports
from .coherent.coherent                   import Coherent
from .incoherent.incoherent               import Incoherent
from .photoelectric.photoelectric         import Photoelectric
from .pairproduction.pairproduction       import Pairproduction
from .tripletproduction.tripletproduction import Tripletproduction



class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'

print(color.BOLD + 'Hello World !' + color.END)



class Photon:
	"""
	DATA BASE DOC: https://www-nds.iaea.org/epics/DOCUMENTS/ENDL2002.pdf
	"""
	def __init__(self, upper_dir, bookmarked_text):

		self.upper_dir = upper_dir


		print("> ** Photon: Reading data for coherent scattering.")
		#COHERENT SCATTERING
		CS = bookmarked_text[(7, 71, 0, 0, 0, 0)]   #integrated cross section
		F  = bookmarked_text[(7, 93, 0, 0, 0, 941)] #form factor
		I  = bookmarked_text[(7, 93, 0, 0, 0, 943)] #imaginary anomalous scattering
		R  = bookmarked_text[(7, 93, 0, 0, 0, 944)] #real anomalous scattering
		E  = bookmarked_text[(7, 71, 0, 0, 7, 10)]  #average energy of the scattered photon


		


      		self.coherent = Coherent(upper_dir, CS, F, R, I, E)
      		
		print(" ")

		print("> ** /photon: Reading data for incoherent scattering.")
		#INCOHERENT SCATTERING
		CS = bookmarked_text[(7, 72, 0, 0, 0, 0)]   #integrated cross section
		S  = bookmarked_text[(7, 93, 0, 0, 0, 942)] #incoherent form factor
		
		self.incoherent = Incoherent(upper_dir, CS, S)
		print(" ")

		print("> ** /photon: Reading data for photoelectric effect.")
		#PHOTOELECTRIC
		CS = bookmarked_text[(7, 73, 0, 0, 0, 0)]
		self.photoelectric = Photoelectric(upper_dir, CS, bookmarked_text)
		     #bookmarked_text[(7, 73, 91, ?, 0, 0)]
		print(" ")

		print("> ** /photon: Reading data for pair production.")
		#PAIR PRODUCTION
		CS = bookmarked_text[(7, 74, 0, 0, 0, 0)]
		self.pairproduction = Pairproduction(upper_dir, CS)
		print(" ")


		print("> ** /photon: Reading data for triplet production.")
		#TRIPLET PRODUCTION
		CS = bookmarked_text[(7, 75, 0, 0, 0, 0)]
		self.tripletproduction = Tripletproduction(upper_dir, CS)
		print(" ")

		self._plotTotalCrossSections()

	def _plotTotalCrossSections(self):

		import matplotlib.pyplot as plt
		fig = plt.figure()
		ax  = plt.gca()

		ax.set_title(f"Integrated Cross Sections for Element {self.upper_dir.Z}.")
		ax.set_xlabel("Energy(MeV)")
		ax.set_ylabel("Cross Section (barn)")
		ax.set_yscale("log")
		ax.set_xscale("log")

		xAxis, yAxis = self.coherent.CS.xAxis, self.coherent.CS.yAxis
		ax.plot(xAxis, yAxis, label="Coherent")

		xAxis, yAxis = self.incoherent.CS.xAxis, self.incoherent.CS.yAxis
		ax.plot(xAxis, yAxis, label="incoherent")

		xAxis, yAxis = self.photoelectric.CS.xAxis, self.photoelectric.CS.yAxis
		ax.plot(xAxis, yAxis, label="photoelectric")

		xAxis, yAxis = self.pairproduction.CS.xAxis, self.pairproduction.CS.yAxis
		ax.plot(xAxis, yAxis, label="pairproduction")

		xAxis, yAxis = self.tripletproduction.CS.xAxis, self.tripletproduction.CS.yAxis
		ax.plot(xAxis, yAxis, label="tripletproduction")
		ax.legend()
		plt.show()




















	def __add__(self, other):

		return CompoundPhoton(self.upper_dir, 
							  self.coherent          + other.coherent, 
							  self.incoherent        + other.incoherent, 
							  self.photoelectric     + other.photoelectric, 
							  self.pairproduction    + other.pairproduction,
			           		  self.tripletproduction + other.tripletproduction)
		
	def __mul__(self, other):
		
		return CompoundPhoton(self.upper_dir, 
							  self.coherent*other, 
							  self.incoherent*other, 
							  self.photoelectric*other, 
							  self.pairproduction*other,
			           		  self.tripletproduction*other)


	def final_init(self, density):
		self.coherent.final_init(density)
		self.incoherent.final_init(density)
		self.photoelectric.final_init(density)
		self.pairproduction.final_init(density)
		self.tripletproduction.final_init(density)

	def __repr__(self):
		return "<" + self.upper_dir._name + ": /photon>"









class CompoundPhoton(Photon):
	def __init__(self,  upper_dir, 
			            coherent,
			            incoherent,
			            photoelectric,
			            pairproduction,
			            tripletproduction):

		self.upper_dir         = upper_dir
		self.coherent          = coherent
		self.incoherent        = incoherent
		self.photoelectric     = photoelectric
		self.pairproduction    = pairproduction
		self.tripletproduction = tripletproduction




