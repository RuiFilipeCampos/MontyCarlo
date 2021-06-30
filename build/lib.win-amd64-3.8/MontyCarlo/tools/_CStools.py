#External Imports
from numpy import *
#from scipy.interpolate import interpl1d

#Internal Imports
from .data import getAxis
from .others import searchsorted
#from .interpol import LinearInterpolation






class TotalCrossSection:
	def __init__(self, upper_dir, CS):

		#creating new element##########################################
		self.upper_dir = upper_dir
		self.CS, self.Iflag, self.data = CS, CS[0], CS[1]
		self.xAxis, self.yAxis = getAxis(self.data)

		assert self.Iflag == 2 or self.Iflag == 0

		Z = upper_dir.Z


		self.interpolations = {LinearInterpolation(Z, self.xAxis, self.yAxis)}
		print(f"> Created a TotalCrossSection with element {upper_dir.Z}.")
		##################################################################


	def __mul__(self, other):
		for x in self.interpolations:
			x.multiply(other)
		return self

	def __add__(self, other):
		new_self = TotalCrossSection.__new__(TotalCrossSection)

		intersection = self.interpolations & other.interpolations

		for x in intersection:
			x.multiply(2)


		new_self.interpolations = self.interpolations | other.interpolations
		new_self.upper_dir = self.upper_dir
		print(new_self.interpolations)

		return new_self

	def final_init(self, density):
		print(f"> Finalized a TotalCrossSection with {len(self.interpolations)} elements.")
		for interpol in self.interpolations:
			interpol.multiply(6.02214076E23 * density/self.upper_dir.A * 1E-24)
			interpol.final_init()

	def __call__(self, x):
		result = 0
		for interpol in self.interpolations:
			result += interpol._eval(x)
		return result


	def _debug(self):
		import numpy as np
		import matplotlib.pyplot as plt

		x = np.arange(0, 20, .1)
		y = [self(x) for x in x]
		plt.plot(x, y)
		plt.show()







class LinearInterpolation:
	"""
	Composition phase: holds yAxis for multiplication.
	"""
	def __init__(self, Z, xAxis, yAxis):
		self.xAxis = xAxis
		self.yAxis = yAxis
		self.Z = Z

	def multiply(self, other):
		self.yAxis = self.yAxis*other

	def __hash__(self):
		return self.Z




	def final_init(self):
		"""
		Creates intervals.
		"""

		self.N = len(self.xAxis) - 1
		self.x_max = self.xAxis[-1]
		self.x_min = self.xAxis[0]

		self.intervals = []

		for i in range(self.N):
			x0, xf = self.xAxis[i], self.xAxis[i+1]
			y0, yf = self.yAxis[i], self.yAxis[i+1]
			self.intervals += [Interval(x0, xf, y0, yf)]



	def _eval(self, x):
		if x <  self.x_min: return 0
		if x == self.x_min: return self.x_min
		if x == self.x_max: return self.x_max

		k = searchsorted(self.xAxis, x, 0, self.N)
		return self.intervals[k-1]._eval(x)





class Interval:
	def __init__(self, x0, xf, y0, yf):
		self.m = (yf-y0)/(xf-x0)
		self.x0, self.xf = x0, xf
		self.y0, self.yf = y0, yf

	def _eval(self, x):
		return self.y0 + self.m * (x - self.x0)

	def __repr__(self):
		return str((self.x0, self.yf))
	def __str__(self):
		return self.__repr__()







#External Imports
from numpy import *
#from scipy.interpolate import interpl1d

#Internal Imports
from .data import getAxis
from .others import searchsorted
#from .interpol import LinearInterpolation






class TotalCrossSection:
	def __init__(self, upper_dir, CS):

		#creating new element##########################################
		self.upper_dir = upper_dir
		self.CS, self.Iflag, self.data = CS, CS[0], CS[1]
		self.xAxis, self.yAxis = getAxis(self.data)

		assert self.Iflag == 2 or self.Iflag == 0

		Z = upper_dir.Z


		self.interpolations = {LinearInterpolation(Z, self.xAxis, self.yAxis)}
		print(f"> Created a TotalCrossSection with element {upper_dir.Z}.")
		##################################################################


	def __mul__(self, other):
		for x in self.interpolations:
			x.multiply(other)
		return self

	def __add__(self, other):
		new_self = TotalCrossSection.__new__(TotalCrossSection)

		intersection = self.interpolations & other.interpolations

		for x in intersection:
			x.multiply(2)


		new_self.interpolations = self.interpolations | other.interpolations
		new_self.upper_dir = self.upper_dir
		print(new_self.interpolations)

		return new_self

	def final_init(self, density):
		print(f"> Finalized a TotalCrossSection with {len(self.interpolations)} elements.")
		for interpol in self.interpolations:
			interpol.multiply(6.02214076E23 * density/self.upper_dir.A * 1E-24)
			interpol.final_init()

	def __call__(self, x):
		result = 0
		for interpol in self.interpolations:
			result += interpol._eval(x)
		return result


	def _debug(self):
		import numpy as np
		import matplotlib.pyplot as plt

		x = np.arange(0, 20, .1)
		y = [self(x) for x in x]
		plt.plot(x, y)
		plt.show()







class LinearInterpolation:
	"""
	Composition phase: holds yAxis for multiplication.
	"""
	def __init__(self, Z, xAxis, yAxis):
		self.xAxis = xAxis
		self.yAxis = yAxis
		self.Z = Z

	def multiply(self, other):
		self.yAxis = self.yAxis*other

	def __hash__(self):
		return self.Z




	def final_init(self):
		"""
		Creates intervals.
		"""

		self.N = len(self.xAxis) - 1
		self.x_max = self.xAxis[-1]
		self.x_min = self.xAxis[0]

		self.intervals = []

		for i in range(self.N):
			x0, xf = self.xAxis[i], self.xAxis[i+1]
			y0, yf = self.yAxis[i], self.yAxis[i+1]
			self.intervals += [Interval(x0, xf, y0, yf)]



	def _eval(self, x):
		if x <  self.x_min: return 0
		if x == self.x_min: return self.x_min
		if x == self.x_max: return self.x_max

		k = searchsorted(self.xAxis, x, 0, self.N)
		return self.intervals[k-1]._eval(x)





class Interval:
	def __init__(self, x0, xf, y0, yf):
		self.m = (yf-y0)/(xf-x0)
		self.x0, self.xf = x0, xf
		self.y0, self.yf = y0, yf

	def _eval(self, x):
		return self.y0 + self.m * (x - self.x0)

	def __repr__(self):
		return str((self.x0, self.yf))
	def __str__(self):
		return self.__repr__()





















































