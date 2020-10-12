
from numpy import *


class LinearInterpolation:
	def __init__(self, xAxis, yAxis):

		self.xAxis = xAxis

		M = []

		for i in range(len(xAxis)):
			m = (yAxis[i] - yAxis[i-1])/(xAxis[i] - yAxis[i-1])
			M += [m]

	def __call__(self, x):
		k = searchsorted(xAxis, x)
		return y[k] + m[k]*(x[k+1] - x)

