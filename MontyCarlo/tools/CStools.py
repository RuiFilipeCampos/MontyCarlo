__doc__ = """DEPRECATED
"""

__author__ = "Rui Campos"



#External Imports
from numpy import *
#from scipy.interpolate import interpl1d

#Internal Imports
from .data import getAxis
from .others import searchsorted
#from .interpol import LinearInterpolation
from .interpolation import LinearInterpolation
from .others import searchsorted




class TotalCrossSection:
	def __init__(self, upper_dir, CS):

		#creating new element##########################################
		self.upper_dir = upper_dir
		self.CS, self.Iflag, self.data = CS, CS[0], CS[1]
		self.xAxis, self.yAxis = getAxis(self.data)

		assert self.Iflag == 2 or self.Iflag == 0

		Z = upper_dir.Z


		self.interpolations = {LinearInterpolation(Z, self.xAxis, self.yAxis)}
		print(f"> *** \t Created a TotalCrossSection with element {upper_dir.Z}.")
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
		return new_self

	def final_init(self, density):
		
		for interpol in self.interpolations:
			interpol.multiply(6.02214076E23 * density/self.upper_dir.A * 1E-24)
			interpol.final_init()
		print(f"> ***  \t Finalized a TotalCrossSection with {len(self.interpolations)} elements.")

	def __call__(self, x):
		result = 0
		for interpol in self.interpolations:
			result += interpol.eval(x)
		return result

	def printTable(self):
		import pandas as pd
		dic = {'Energy(MeV)':self.xAxis, 'CS(barns)':self.yAxis}
		df  = pd.DataFrame(dic)
		print(df)

	def _debug(self):
		import numpy as np
		import matplotlib.pyplot as plt

		
		#y = [self(x) for x in self.xAxis]
		plt.plot(self.xAxis, self.yAxis)
		plt.show()

