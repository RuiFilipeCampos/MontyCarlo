# cython: profile=False

print("Importing .materials.database")

__doc__ = """
DATABASE DOCUMENTATION: 
	https://www-nds.iaea.org/epics/DOCUMENTS/ENDL2002.pdf
	https://www-nds.iaea.org/epics/
"""

__author__ = "Rui Campos"



#External Imports
from numpy import * #array, geomspace, flip, load, searchsorted

#Internal Imports
from ..tools.CubicInverseTransform import makeAlias
from ..settings import __montecarlo__
from ..tools.data import getAxis
from ..tools.interpol1 import LinLinInterpolation


__materials__ = __montecarlo__/'materials'
directory = str(__materials__)

# META
__directory__ = __montecarlo__/'materials'



##################### CONSTANTS
Na = 6.02214129e23 #1/mol
c = 2.99792458e8
alpha = 1/137.035999074
hbar = 6.58211928e16 #eVs
e = 1.602176565e-19 #C
m_electron = 9.10938291e-31 #kg
E0_electron = 510.998928 #keV





class EADL:
	def __init__(self):
		self.container = [EADLelement(Z) for Z in range(1, 101)]

	def __getitem__(self, Z):
		return self.container[Z-1]
			



class EADLelement:
	def __init__(self, Z):
		self.Z = Z
	
		file = str(Z) + ".txt"
		self.path = str(__materials__/'EADL'/file)
		del file
	
		#self.path = directory + "\\EADL\\" + str(Z) + ".txt"
		self.Aw, self.EADL_dict = self.getBookmarkedText()
		self.container = {}
		
		for Id in self.EADL_dict:
			if Id[0:3] == (0, 92, 91) and Id[4:] == (7, 931):
				j_, fr_, Er_ = [], [], []
				for line in self.EADL_dict[Id]:
					j, fr, Er = [float(x) for x in line.split()]
					j_ += [j]
					fr_ += [fr]
					Er_ += [Er]

				self.container[Id] = tuple(map(array, (j_, fr_, Er_) ))
			   
			if Id[0:3] == (0, 92, 91) and Id[4:] == (9, 932):
				j_, k_, fnr_, Enr_ = [], [], [], []
			   
				for line in self.EADL_dict[Id]:
					j, k, fnr, Enr = [float(x) for x in line.split()]
				   
					j_ += [j]
					k_ += [k]
					fnr_ += [fnr]
					Enr_ += [Enr]
			
				
				self.container[Id] = tuple(map(array, (j_, k_, fnr_, Enr_) ))
		
		
		

		
		
	def getBookmarkedText(self):
		path = self.path

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
			
			line0 = text[0]
			Z  = int(line0[0:3])
			Aw = float(line0[13:24])
			
	
			bookmarked_text = {}
			
			for i, j in zip(bookmark_, _ookmarks):
					line1, line2 = text[i], text[i+1]
					
				   #on line 1
					Yi = float(line1[7:9])    #particle identifier
					Yo = float(line1[10:12])  #secondary particle designator
					
					#on line 2
					C  = float(line2[0:2])    #reaction descriptor
					I  = float(line2[2:5])    #reaction property
					S  = float(line2[5:8])    #reaction modifier
					X1 = float(line2[22:32])  #subshell designator
					
					flags = (Yi, C, S, X1, Yo, I)

					flags = tuple(map(int, flags))
					bookmarked_text[flags] = text[i+2:j-1]
						
		return Aw, bookmarked_text



def get_bookmarked_text_EADL(path):

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

		line0 = text[0]
		Z  = int(line0[0:3])
		Aw = float(line0[13:24])

		bookmarked_text = {}

		for i, j in zip(bookmark_, _ookmarks):
				line1, line2 = text[i], text[i+1]
				
			 	#on line 1
				Yi = float(line1[7:9])    #particle identifier
				Yo = float(line1[10:12])  #secondary particle designator
				
				#on line 2
				C  = float(line2[0:2])    #reaction descriptor
				I  = float(line2[2:5])    #reaction property
				S  = float(line2[5:8])    #reaction modifier
				X1 = float(line2[22:32])  #subshell designator
				
				flags = (Yi, C, S, X1, Yo, I)

				flags = tuple(map(int, flags))
				bookmarked_text[flags] = text[i+2:j-1]
					
	return Aw, bookmarked_text
















def getEADL(Z):
	"""
	EADL DOC: 
		https://drive.google.com/file/d/1i5ndh-G6eD1ginpxNLzBtskJh3XTU5p9/
	"""
	
	file = str(Z) + ".txt"
	EADL_path = str(__materials__/'EADL'/file)
	del file
	
	Aw, EADL_dict = get_bookmarked_text_EADL(EADL_path)

	EADL_dict['Aw'] = Aw
	EADL_dict['Z']  = Z 
	
	return EADL_dict
	


		

	


 






























def getEPDL(Z):
	"""
	DATABASE DOCS:
		https://drive.google.com/file/d/1_Dtsfd4A18m1BsZPCyWfm6PjcWxjyf1n/
		
	
	cullen1997
	https://drive.google.com/file/d/1WqppBrR-C3yiRuhp7P16c0yCPDFr2T9t/view?usp=sharing
	"""
	file = str(Z) + ".txt"
	EPDL_path = str(__materials__/'EPDL'/file)
	del file

	#EPDL_path      = directory + "\\EPDL\\" + str(Z) + ".txt"

	EPDL_dict = get_bookmarked_text(EPDL_path)
	
	

	for Id in EPDL_dict:
		EPDL_dict[Id] = Table(EPDL_dict[Id], Id, Z)
		
	
	return EPDL_dict

def getEEDL(Z):
	"""
	DATABASE DOCS:
		https://drive.google.com/file/d/1ef8Ww_0PPWuLEuwpfv9KOF4wyd75_eOp/

	"""
	Z = int(Z)

	file = str(Z) + ".txt"
	EEDL_path = str(__materials__/'EEDL'/file)
	del file

	#EEDL_path = directory + r"\\EEDL\\" + str(Z) + ".txt"

	EEDL_dict = get_bookmarked_text(EEDL_path)
	
	
	for Id in EEDL_dict:
		EEDL_dict[Id] = EEDLtable(EEDL_dict[Id], Id, Z)
		
	return EEDL_dict

















def get_bookmarked_text(path):
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



# def get_bookmarked_text2(path):
#         """
#         Reads EPICS file format and returns a dict with flags as keys
#         and data as a list of strings.
#         """
#         with open(path, "r") as file:
#                 text = file.readlines()
#                 text = [line.strip('\n') for line in text]

#                 bookmarks = [0]
				
#                 for n, line in enumerate(text):
#                         if line == "                                                                       1":
#                                 bookmarks += [n + 1]
								
#                 #gather all bookmarked text into a dict
#                 bookmark_ = bookmarks[0:-1]
#                 _ookmarks = bookmarks[1:]
		
#                 bookmarked_text = {}
#                 for i, j in zip(bookmark_, _ookmarks):
#                         line1, line2 = text[i], text[i+1]

#                         #on line 1
#                         Z = float(line1[0:3])     #
#                         A = float(line1[3:6])  #secondary particle designator
#                         Yi = float(line1[7:9])  #interpolation flag
#                         Yo = float(line1[10:12])
#                         AW = float(line1[13:24])
						
#                         #on line 2
#                         #C  = float(line2[0:2])    #reaction descriptor
#                         #I  = float(line2[2:5])    #reaction property
#                         #S  = float(line2[5:8])    #reaction modifier
#                         #X1 = float(line2[22:32])  #subshell designator
						
#                         flags = int(Z)

#                         #flags = tuple(map(int, flags))
#                         bookmarked_text[flags] =  text[i+2:j-1]
						
#         return bookmarked_text











class Table:
	
	IDtranslation = {(7, 71, 0, 0, 0, 0): "Coherent",
					 (7, 71, 0, 0, 7, 10): "FormFactor",
					 (7, 72, 0, 0, 0, 0): "Incoherent",
					 (7, 72, 0, 0, 7, 10): "IncoherentFormFactor",
					 (7, 72, 0, 0, 9, 10): "???",
					 (7, 73, 0, 0, 0, 0): "???",
					 (7, 73, 91, 1, 0, 0): "???",
					 (7, 73, 91, 1, 0, 11): "???",
					 (7, 73, 91, 1, 9, 10): "???",
					 (7, 74, 0, 0, 0, 0): "???",
					 (7, 74, 0, 0, 8, 10): "???",
					 (7, 74, 0, 0, 9, 10): "???",
					 (7, 75, 0, 0, 0, 0): "???",
					 (7, 75, 0, 0, 8, 10): "???",
					 (7, 75, 0, 0, 9, 10): "???",
					 (7, 93, 0, 0, 0, 941): "???",
					 (7, 93, 0, 0, 0, 942): "???",
					 (7, 93, 0, 0, 0, 943): "???",
					 (7, 93, 0, 0, 0, 944): "???"}
	
	def __init__(self, EPDLtable, Id, Z):
		
   
		self.Z = Z
		self.Id = Id
		#print(EPDLtable)
		self.Iflag, rawData = EPDLtable[0], EPDLtable[1]
		
		try:    self.X, self.Y = getAxis(rawData)
		except: self.rawData = rawData
	
	def getLinLinInterpol(self):
		return LinLinInterpolation(self.X, self.Y)

	def getLogLogInterpol(self):
		
		f = self.getLinLinInterpol()
		
		x_min, x_max = min(self.X), max(self.X)
		X_logspaced = logspace(x_min, x_max, num=100)
		
		Y_logspaced = [f(x) for x in X_logspaced]
		Y_logspaced = array(X_logspaced)
		
		return LogLogInterpol(X_logspaced, Y_logspaced)

	def scatter(self):
		import matplotlib.pyplot as plt
		fig = plt.figure(figsize=(10, 10))
		plt.scatter(self.X, self.Y, s=3)
		plt.xscale("log")
		plt.yscale("log")
		plt.title("EADL table id = " + str(self.Id) + f" of element {self.Z}")
		plt.grid(which = 'both')
		plt.show()




class EEDLtable:
	
	
	def __init__(self, table, Id, Z):
		
		self.Z = Z
		self.Id = Id
		#print(EPDLtable)
		self.Iflag, rawData = table[0], table[1]
		Ncolumns = len(rawData[0].split())
		self.Ncolumns = Ncolumns
		if Ncolumns == 2:
			self.X, self.Y = getAxis(rawData)
			return

	
		axis1 = []
		axis2 = []
		axis3 = []
		
		for line in rawData:
			numbers = [float(x) for x in line.split()]
			axis1 += [numbers[0]]
			axis2 += [numbers[1]]
			axis3 += [numbers[2]]
		
		self.E = array(axis1)
		
		
		self.Y1 = {E:[] for E in self.E}
		self.Y2 = {E:[] for E in self.E}
		for i, E in enumerate(self.E):
			self.Y1[E].append(axis2[i])
			self.Y2[E].append(axis3[i])
			

	
			


	def scatter(self, E):
		import matplotlib.pyplot as plt
		fig = plt.figure(figsize=(10, 10))
		plt.scatter(self.Y1[E], self.Y2[E], s=3)
		plt.xscale("log")
		plt.yscale("log")
		plt.title("EADL table id = " + str(self.Id) + f" of element {self.Z}")
		plt.grid(which = 'both')
		plt.show()

			






class _EPDL:
	__cache__ = {}
	def __getitem__(self, Z):
		try:
			return self.__cache__[Z+1]
		except KeyError:
			self.__cache__[Z+1] = getEPDL(Z+1)
			return self.__cache__[Z+1]
	def __call__(self, Z):
		try:
			return self.__cache__[Z-1]
		except KeyError:
			self.__cache__[Z-1] = getEPDL(Z)
			return self.__cache__[Z-1]
	
		
	
class _EEDL:
	__cache__ = {}
	def __getitem__(self, Z):
		try:
			return self.__cache__[Z+1]
		except KeyError:
			self.__cache__[Z+1] = getEEDL(Z+1)
			return self.__cache__[Z+1]
		
	def __call__(self, Z):
		try:
			return self.__cache__[Z-1]
		except KeyError:
			self.__cache__[Z-1] = getEEDL(Z)
			return self.__cache__[Z-1]
	













from scipy.interpolate import CubicSpline


import pickle
from os import path

if not path.exists(directory + "\\data.pkl"):
	print("""
	________________
	> Reading EADL. """)
	#EADL = EADL()
	EADL = [getEADL(Z) for Z in range(1, 101)]
	print("""> Done! EPDL available @ database.EPDL
	¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯""")


	
	EPDL = _EPDL()
	
	print("""
	________________
	> Reading EPDL. """)
	#EPDL = [getEPDL(Z) for Z in range(1, 101)] 
	print("""> Done! EPDL available @ database.EPDL
	¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯""")
		
	EEDL = _EEDL()
	
	print("""
	________________
	> Reading EEDL. """)
	
	#EEDL = EEDL()
	#EEDL = [getEEDL(Z) for Z in range(1, 101)]
	print("""> Done! EPDL available @ database.EEDL
	¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯
	""")
	
	
	print("Imported databases")


	data = [EADL, EPDL, EEDL]


	with open(directory + "\\data.pkl", 'wb') as output:
		pickle.dump(data, output)
		
	print("Saved database.")
else:
	with open(directory + "\\data.pkl", 'rb') as pkl_file:
		data = pickle.load(pkl_file)
	
	EADL, EPDL, EEDL = data
	
	
	
### designators

epdl = EPDL[99]



designator_to_index = dict()

i = 0
for key in epdl:
	if key[0] == 7 and key[1] == 73 and key[2] == 91 and key[4] == 0 and key[5] == 0: #(7, 73, 91, designator, 0, 0)
		
		designator_to_index[key[3]] = i
		i += 1
	
	


del epdl
del key


import numpy as np 
class MoleculeDATA:
	def __init__(self, formula):
		self.ATOMS = []
		
		for Z, x in formula.items():
			self.ATOMS.append(AtomDATA(Z, x))
			
			
			
		# Making aliases for choosing shell based on occupation numbers
		
		self.Nsh = 0
		i = 0
		I = []
		info = []
		probs = []
		
		atom_index = 0
		for atom in self.ATOMS:
			self.Nsh += len(atom.number_of_electrons)
			
			
			for shell_index, Nel in enumerate(atom.number_of_electrons):
				info.append([atom_index, shell_index])    
				I.append(i)
				probs.append(atom.x*Nel)
				i += 1
			atom_index += 1
			
			
		probs = np.array(probs)
		probs = probs/sum(probs)
		index = array(I)
		
		ALIAS = makeAlias(index, probs)
		
		self.ALIAS = []
		for k, line in enumerate(ALIAS):
			
			
			a = int(line[0])
			p = line[1]
			b = int(line[2])
			
			new_line =  [p] + info[a] + info[b]
			
			self.ALIAS.append(new_line)
			
		self.ALIAS = np.array(self.ALIAS)


	def test_alias(self, N = 10_000):
		from numpy.random import rand
		samples = []
		#choose randomly in array
		for _ in range(N):
			r = rand()*self.Nsh
			N = int(r)
			
			if r - N < self.ALIAS[N, 0]: #accept
			   # sample = ( self.ALIAS[N, 1], self.ALIAS[N, 2])    
				samples.append(N)
			
			else:
				#sample = ( self.ALIAS[N, 3], self.ALIAS[N, 4])
				for n, line in enumerate(self.ALIAS):
					if line[1] == self.ALIAS[N, 3] and line[2] == self.ALIAS[N, 4]:
						samples.append(n)
				
		return samples
			
				
		
			
			
			
			
		#     #construct prob array with indexes like: (Z, shell_index), maybe in form of arr, array([Z, shell_index ])        
		# probs = self.number_of_electrons / sum(self.number_of_electrons)
		# index = np.arange(0, len(probs))
		# self.ALIAS = makeAlias(index, probs)
			
			
			
	def __repr__(self):
		rep = "<Molecule DATA :: "
		
		for atom in self.ATOMS:
			rep += f"{atom.x}x(Z = {atom.Z}) "
		
		rep += ">"
		return rep
		
	def __getitem__(self, *args, **kwargs):
		s = args[0]
		if isinstance(s, int) or isinstance(s, float):
			return self.ATOMS[int(s)]
		
		if isinstance(s, slice):
			if s.start == None and s.stop == None:
				Z = s.step
				for atom in self.ATOMS:
					if atom.Z == Z:
						return atom
			
		else: return NotImplemented
	
	def __str__(self):
		to_print = self.__repr__() + "\n"
		to_print += "\n"
		
		for atom in self.ATOMS:
			to_print += atom.__str__()
			to_print += "\n \n"
		
		return to_print
	
class AtomDATA:
	def __init__(self, Z, x):
		self.x = x
		self.Z = Z

		file_name = str(Z) + '.txt'
		self.path = str(__directory__/'EADL'/file_name)
		del file_name

		self.Aw, self.EADL_dict = self.getBookmarkedText()

		for key, item in self.EADL_dict.copy().items():
			content = self.EADL_dict[key]

			replace = []
			for line in content:
				numerical_line = [float(x) for x in line.split()]
				replace.append(numerical_line)

			self.EADL_dict[key] = np.array(replace)


		self.number_of_electrons = self.EADL_dict[(0, 91, 0, 0, 0, 912)][:, 1]
		self.binding_energy = self.EADL_dict[(0, 91, 0, 0, 0, 913)][:, 1]*1e6
		self.kinetic_energy = self.EADL_dict[(0, 91, 0, 0, 0, 914)][:, 1]

		J0 = []
		
		
		J0path = str(__directory__/'comptonJ0.txt')
		with open(J0path) as file:
			text = file.readlines()[2:]
			text = [line.strip('\n') for line in text]
			for line in text:
				numbers = [float(x) for x in line.split()]
				if numbers[0] == self.Z:
					J0.append(numbers[3])
		self.J0 = np.array(J0)
		
		


		
				
		
		
		
		
	def __repr__(self):
		return f"<Atom DATA :: {self.Z} >"

	def __str__(self):
		to_print = self.__repr__() + "\n"
		for i, Nel in enumerate(self.number_of_electrons):
			BE = self.binding_energy[i]
			J0 = self.J0[i]
			to_print += f"    <Shell #{i} ::  Nel = {Nel} | BE = {BE} eV | J0 = {J0}> \n"
		return to_print
		
		
		
		
		
		
		
		
		
		
   
		
		
		

		
		
	def getBookmarkedText(self):
		path = self.path

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
				
				line0 = text[0]
				Z  = int(line0[0:3])
				Aw = float(line0[13:24])
				
		
				bookmarked_text = {}
				
				for i, j in zip(bookmark_, _ookmarks):
						line1, line2 = text[i], text[i+1]
						
					   #on line 1
						Yi = float(line1[7:9])    #particle identifier
						Yo = float(line1[10:12])  #secondary particle designator
						
						#on line 2
						C  = float(line2[0:2])    #reaction descriptor
						I  = float(line2[2:5])    #reaction property
						S  = float(line2[5:8])    #reaction modifier
						X1 = float(line2[22:32])  #subshell designator
						
						flags = (Yi, C, S, X1, Yo, I)

						flags = tuple(map(int, flags))
						bookmarked_text[flags] = text[i+2:j-1]
						
		return Aw, bookmarked_text
		
		
		
	
class ShellDATA:
	def __init__(self, index, designator, number_el, binding_energy, J0):
		self.designator = designator
		self.index = index
		self.J0 = J0
		
		
		pass












