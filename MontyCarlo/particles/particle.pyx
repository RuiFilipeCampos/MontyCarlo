# distutils: language = c++

print("Importing `.particles.particle` ...")


# Conditional compilation 
DEF FULL_RECORD = True



# External Imports
cimport cython
import numpy as np
cimport numpy as cnp
from libc.math cimport sin, cos, log, pow, sqrt
import numpy.random
from collections import deque
from cython.operator import dereference as deref
from cython.operator import preincrement as inc
#from numpy.math cimport INFINITY

rand = numpy.random.rand

# Internal Imports
from ..settings import DEBUG
from ..tools.performance import timer
from ..geometry.main    cimport Volume
from ..types cimport double3



cdef struct STATE:
		mixmax_engine* genPTR
		void *current_region
		double3 pos
		double3 dire
		double3 axis
		double E
		double L 
		double last_displacement




@cython.boundscheck(False)
@cython.initializedcheck(False)
@cython.cdivision(True)
@cython.freelist(100_000_000)
cdef class Particle:
		"""Base class for all particle types.

		Attributes:
			STATE state :: struct that holds the particle state
			int nSECONDARY :: number of secondary particles produced
			cdef public object FILE :: file object used for debugging purposes
			cdef long double imfp_T :: total inverse mean free path
			cdef object secondary   :: python container for the secondary particles 

		"""


		cdef inline double ENERGY(self):
				raise RuntimeError("ENERGY METHOD CALLED FROM Particles")


		cpdef get_record_pos(self):
			"""Return the particles trajectory, making it available on python.

			Note:
					If `self.state_record` does not exist (see: conditional compilation on this module), it
					will cause segmentation error. 

			Returns:
					Python deque holding the particles trajectory. (??): deque([(x,y,z), (x,y,z)])
			"""


			cdef vector[STATE].iterator it = self.state_record.begin()
			
			traj = deque()

			cdef int i
			for i in range(self.state_record.size()):
				traj.append(deref(it).pos)
				inc(it)

			return traj

		cdef void record(self):
				"""
				Add current position and energy to the track record.
				"""
				
				IF FULL_RECORD:
						self.state_record.push_back(self.state)
				#ELSE:
				#    self.pos_record.push_back(self.state.pos)
				#    self.E_record.push_back(self.state.E)

				
		def getEnergy(self):
				return self.state.E
				
		cdef void move(self, double L):
				"""Move particle a displacement `L` in the current direction.

				Note:
						This function is not inlined. Using it frequently will cause performance issues.

				Args:
						L (double): The particle displacement.

				Returns:
						(void)
				"""

				# self.state.pos.x += self.state.EZ.x*L
				# self.state.pos.y += self.state.EZ.y*L
				# self.state.pos.z += self.state.EZ.z*L

				self.state.pos.x += self.state.dire.x*L
				self.state.pos.y += self.state.dire.y*L
				self.state.pos.z += self.state.dire.z*L

		cdef inline void deposit(self):
				"""Deprecated: Energy deposition is performed by `Volume` instances.
				"""
				pass
		
	
		cdef double rsqrt(self, double x):
				"""Fourth order approximation on the inverse square root.

				Note:
						Used to normalize the direction and polarization vectors.

				Args:
						x (double)

				Returns:
						(double) approximation to 1/sqrt(x)
				"""

				x = 1 - x
				cdef double y = 1 + x*.5  # 0th and 1st order
				x *= x              
				y += x * 0.375            # 2nd order
				x *= x
				y += x * 0.3125
				x *= x
				y += x * 0.2734375
				return y
		


		cdef void normalize(self):
				"""Normalize the direction and polarization vectors.

				Note:
						Uses a first order approximation on 1/sqrt(x). As long as this is done
						frequently, the numeric error is kept under control with minimal 
						performance hits.
				"""

				cdef double x
				cdef double norm
				
				x = 1 - (self.state.dire.x*self.state.dire.x + self.state.dire.y*self.state.dire.y + self.state.dire.z*self.state.dire.z )
				norm = 1 + x*.5  # 0th and 1st order
				
				self.state.dire.x *= norm
				self.state.dire.y *= norm
				self.state.dire.z *= norm


		cdef void invert_axis(self):
				self.state.axis.x *= -1
				self.state.axis.y *= -1
				self.state.axis.z *= -1


		cdef void invert_dire(self):
				self.state.dire.x *= -1
				self.state.dire.y *= -1
				self.state.dire.z *= -1


		cdef void throwAZIMUTH(self):
				cdef double x, y, a, w1, w2, w3

				while True:
						while True:
								x = 2*self.state.genPTR.get_next_float() - 1
								y = 2*self.state.genPTR.get_next_float() - 1
								a = x**2 + y**2
								if a < 1:
										break


						w1 = 1 - 2*a 
						a = 2 * sqrt(1 - a)
						w2 = x*a
						w3 = y*a
						a =  w1*self.state.dire.x +  w2*self.state.dire.y + w3*self.state.dire.z
						#w and ez are unit vectors, this should be enough to guarentee that they are not in the same direc
						if a < 1: 
								
								x = 1/sqrt(1 - a**2)
								self.state.axis.x = (w1 - a*self.state.dire.x)*x
								self.state.axis.y = (w2 - a*self.state.dire.y)*x
								self.state.axis.z = (w3 - a*self.state.dire.z)*x
								return
								
								# x = 1/(self.eyx**2 + self.eyy**2 + self.eyz**2)**.5
								
								# self.eyx *= x
								# self.eyy *= x
								# self.eyz *= x
								

								 
				# else: 
				#     import time
						
				#     print("22222111ASdsasadjlkdsja lksadalkdsfjlkdasf lkjdsahf kjsad hfkjda hsflkj hdsaflkj hadlkjf halskjd")
						
				#     time.sleep(1000000)
				#     raise RuntimeError("Exceeded max")
		
		
		
		cdef void rotateTHETA(self, double cos):
				# spent an entire afternoon looking for a serious bug
				# as it happens, a 1.0000000000000002 value of cos will result in NaN values on the rotated vector
				
				
				if cos > 1 :
						cos =  1
				elif cos < -1: 
						cos = -1
				
				

				
				cdef double cos2sq = .5*(1 + cos)
				cdef double cos2 = sqrt(cos2sq)
				
				cdef double sin2sq = 1 - cos2sq
				cdef double sin2 = sqrt(sin2sq)
				cdef double sin2cos2 = 2*sin2*cos2
				
				cdef double eyx2 = self.state.axis.x**2
				cdef double eyy2 = self.state.axis.y**2
				cdef double eyz2 = self.state.axis.z**2
				
				cdef double ezx = self.state.dire.x
				cdef double ezy = self.state.dire.y
				cdef double ezz = self.state.dire.z
				
				#### calculated formulas using sympy's quaternions
				#### copied it here and organized best I could for performance
				self.state.dire.x = \
				cos2sq* ezx + \
				sin2sq* (eyx2 * ezx \
								+ 2*self.state.axis.x  * (self.state.axis.y*ezy + self.state.axis.z*ezz) \
								- ezx * (eyy2  + eyz2))         \
				+ sin2cos2*(self.state.axis.y*ezz - self.state.axis.z*ezy)
				
				

				self.state.dire.y = \
				cos2sq*ezy + \
				sin2sq*( ezy*( -eyx2 \
								+ eyy2 \
								- eyz2)
								+ 2*self.state.axis.y*(self.state.axis.x*ezx \
														+ self.state.axis.z*ezz)) \
				+ sin2cos2*(-self.state.axis.x*ezz + self.state.axis.z*ezx)

				self.state.dire.z = cos2sq*ezz + \
				sin2sq*(- ezz*   (eyx2 + eyy2 - eyz2) \
								+ 2*self.state.axis.z*(self.state.axis.x*ezx  + self.state.axis.y*ezy)) \
				+ sin2cos2*(self.state.axis.x*ezy - self.state.axis.y*ezx)


			 # print("")
			#  print("ey", self.eyx, self.eyy, self.eyz)
			 # print("norm:", self.eyx**2 + self.eyy**2 + self.eyz**2 )
			#  print("ez", self.ezx, self.ezy, self.ezz)
			#  print("norm:", self.ezx**2 + self.ezy**2 + self.ezz**2 )




# ('ey', 0.304192241115606, 0.7845300950417162, 0.5403513768085548)
# ('norm:', 1.0)
# ('ez', 0.09865239061376316, 0.5382383859203108, -0.8369988923219004)
# ('norm:', 1.0)


				#cdef double x
 
				#x = (self.ezx**2 + self.ezy**2 + self.ezz**2 )
				
				
				#nevermind the names, just using already allocated variables
				cos2sq = 1 - (self.state.dire.x**2 + self.state.dire.y**2 + self.state.dire.z**2 )
				ezz = 1 + cos2sq*.5  # 0th and 1st order
				cos2sq *= cos2sq              
				ezz += cos2sq * 0.375            # 2nd order
				cos2sq *= cos2sq
				ezz += cos2sq * 0.3125
				cos2sq *= cos2sq
				ezz += cos2sq * 0.2734375
				
				self.state.dire.x *= ezz
				self.state.dire.y *= ezz
				self.state.dire.z *= ezz
				
				
				# print("")
				# print("ey", self.eyx, self.eyy, self.eyz)
				# print("norm:", self.eyx**2 + self.eyy**2 + self.eyz**2 )
				# print("ez", self.ezx, self.ezy, self.ezz)
				# print("norm:", self.ezx**2 + self.ezy**2 + self.ezz**2 )
				
				
				# print("rotateTHETA out")




		cdef void rotateTHETAvers2(self, double cos):
				# spent an entire afternoon looking for a serious bug
				# as it happens, a 1.0000000000000002 value of cos will result in NaN values on the rotated vector
				if cos > 1 :  
						cos =  1
				elif cos < -1: 
						cos = -1
				
				
				
				cdef double cos2sq = .5*(1 + cos)
				cdef double cos2 = sqrt(cos2sq)
				
				cdef double sin2sq = 1 - cos2sq
				cdef double sin2 = sqrt(sin2sq)
				cdef double sin2cos2 = 2*sin2*cos2
				
				cdef double eyx2 = self.state.axis.x**2 # eyx**2
				cdef double eyy2 = self.state.axis.y**2 # eyy**2
				cdef double eyz2 = self.state.axis.z**2 # eyz**2
				
				cdef double ezx = self.state.dire.x # ezx
				cdef double ezy = self.state.dire.y  #ezy
				cdef double ezz = self.state.dire.z # ezz
				
				#### calculated formulas using sympy's quaternions
				#### copied it here and organized best I could for performance
				self.ezx = \
				cos2sq* ezx + \
				sin2sq* (eyx2 * ezx \
								+ 2*self.eyx  * (self.eyy*ezy + self.eyz*ezz) \
								- ezx * (eyy2  + eyz2))         \
				+ sin2cos2*(self.eyy*ezz - self.eyz*ezy)
				
				

				self.ezy = \
				cos2sq*ezy + \
				sin2sq*( ezy*( -eyx2 \
								+ eyy2 \
								- eyz2)
								+ 2*self.eyy*(self.eyx*ezx \
														+ self.eyz*ezz)) \
				+ sin2cos2*(-self.eyx*ezz + self.eyz*ezx)

				self.ezz = cos2sq*ezz + \
				sin2sq*(- ezz*   (eyx2 + eyy2 - eyz2) \
								+ 2*self.eyz*(self.eyx*ezx  + self.eyy*ezy)) \
				+ sin2cos2*(self.eyx*ezy - self.eyy*ezx)


			 # print("")
			#  print("ey", self.eyx, self.eyy, self.eyz)
			 # print("norm:", self.eyx**2 + self.eyy**2 + self.eyz**2 )
			#  print("ez", self.ezx, self.ezy, self.ezz)
			#  print("norm:", self.ezx**2 + self.ezy**2 + self.ezz**2 )




# ('ey', 0.304192241115606, 0.7845300950417162, 0.5403513768085548)
# ('norm:', 1.0)
# ('ez', 0.09865239061376316, 0.5382383859203108, -0.8369988923219004)
# ('norm:', 1.0)


				#cdef double x
 
				#x = (self.ezx**2 + self.ezy**2 + self.ezz**2 )
				
				
				#nevermind the names, just using already allocated variables
				cos2sq = 1 - (self.ezx**2 + self.ezy**2 + self.ezz**2 )
				ezz = 1 + cos2sq*.5  # 0th and 1st order
				cos2sq *= cos2sq              
				ezz += cos2sq * 0.375            # 2nd order
				cos2sq *= cos2sq
				ezz += cos2sq * 0.3125
				cos2sq *= cos2sq
				ezz += cos2sq * 0.2734375
				
				self.ezx *= ezz
				self.ezy *= ezz
				self.ezz *= ezz
				
				
				# print("")
				# print("ey", self.eyx, self.eyy, self.eyz)
				# print("norm:", self.eyx**2 + self.eyy**2 + self.eyz**2 )
				# print("ez", self.ezx, self.ezy, self.ezz)
				# print("norm:", self.ezx**2 + self.ezy**2 + self.ezz**2 )
				
				
				# print("rotateTHETA out")
				










		cdef void rotateAZIMUTH(self, double cos):
				# spent an entire afternoon looking for a serious bug
				# as it happens, a 1.0000000000000002 value of cos will result in NaN values on the rotated vector
				if cos > 1 :  
						
						cos =  1
				elif cos < -1: 
						
						cos = -1
				
				
				#cos = (<double> 1e-15 * (<int> (1e15* cos)))
				cdef double tempEYx = self.eyx
				cdef double tempEYy = self.eyy
				cdef double tempEYz = self.eyz
				
				cdef double tempEZx = self.ezx
				cdef double tempEZy = self.ezy
				cdef double tempEZz = self.ezz

				
				self.eyx = tempEZx
				self.eyy = tempEZy
				self.eyz = tempEZz
				
				self.ezx = tempEYx
				self.ezy = tempEYy
				self.ezz = tempEYz





				#print(f"rotateTHETA in -- cos = {cos}")
				
				#cdef double eyx = self.eyx
				#cdef double eyy = self.eyy
				#cdef double eyz = self.eyz
				#print("ey", self.eyx, self.eyy, self.eyz)
			 # print("norm:", self.eyx**2 + self.eyy**2 + self.eyz**2 )
			 # print("ez", self.ezx, self.ezy, self.ezz)
			#  print("norm:", self.ezx**2 + self.ezy**2 + self.ezz**2 )
				
				cdef double cos2sq = .5*(1 + cos)
				cdef double cos2 = sqrt(cos2sq)
				
				cdef double sin2sq = 1 - cos2sq
				cdef double sin2 = sqrt(sin2sq)
				cdef double sin2cos2 = 2*sin2*cos2
				
				cdef double eyx2 = self.eyx**2
				cdef double eyy2 = self.eyy**2
				cdef double eyz2 = self.eyz**2
				
				cdef double ezx = self.ezx
				cdef double ezy = self.ezy
				cdef double ezz = self.ezz
				
				#### calculated formulas using sympy's quaternions
				#### copied it here and organized best I could for performance
				self.ezx = \
				cos2sq* ezx + \
				sin2sq* (eyx2 * ezx \
								+ 2*self.eyx  * (self.eyy*ezy + self.eyz*ezz) \
								- ezx * (eyy2  + eyz2))         \
				+ sin2cos2*(self.eyy*ezz - self.eyz*ezy)
				
				

				self.ezy = \
				cos2sq*ezy + \
				sin2sq*( ezy*( -eyx2 \
								+ eyy2 \
								- eyz2)
								+ 2*self.eyy*(self.eyx*ezx \
														+ self.eyz*ezz)) \
				+ sin2cos2*(-self.eyx*ezz + self.eyz*ezx)

				self.ezz = cos2sq*ezz + \
				sin2sq*(- ezz*   (eyx2 + eyy2 - eyz2) \
								+ 2*self.eyz*(self.eyx*ezx  + self.eyy*ezy)) \
				+ sin2cos2*(self.eyx*ezy - self.eyy*ezx)


			 # print("")
			#  print("ey", self.eyx, self.eyy, self.eyz)
			 # print("norm:", self.eyx**2 + self.eyy**2 + self.eyz**2 )
			#  print("ez", self.ezx, self.ezy, self.ezz)
			#  print("norm:", self.ezx**2 + self.ezy**2 + self.ezz**2 )




# ('ey', 0.304192241115606, 0.7845300950417162, 0.5403513768085548)
# ('norm:', 1.0)
# ('ez', 0.09865239061376316, 0.5382383859203108, -0.8369988923219004)
# ('norm:', 1.0)


				#cdef double x
 
				#x = (self.ezx**2 + self.ezy**2 + self.ezz**2 )
				
				
				#nevermind the names, just using already allocated variables
				cos2sq = 1 - (self.ezx**2 + self.ezy**2 + self.ezz**2 )
				ezz = 1 + cos2sq*.5  # 0th and 1st order
				cos2sq *= cos2sq              
				ezz += cos2sq * 0.375            # 2nd order
				cos2sq *= cos2sq
				ezz += cos2sq * 0.3125
				cos2sq *= cos2sq
				ezz += cos2sq * 0.2734375
				
				self.ezx *= ezz
				self.ezy *= ezz
				self.ezz *= ezz



				self.eyx = self.ezx
				self.eyy = self.ezy
				self.eyz = self.ezz
				
				self.ezx = tempEZx
				self.ezy = tempEZy
				self.ezz = tempEZz






		cpdef add_to_cell(self, object old_cell, object old_points,  int Npoints):
				
				
				cdef int N = self.X.size()
				cdef int i
				
				old_cell.append(N)
				
				cdef int x
				for x in range(0, N):
						old_cell.append(x + Npoints)
			 # old_cell.extend(x + Npoints for x in range(0, N))
				
				cdef cnp.ndarray points = np.zeros((N, 3))
				
				
				
				for i in range(N):
						points[i, 0] = self.X[i]
						points[i, 1] = self.Y[i]
						points[i, 2] = self.Z[i]
				
				old_points.extend(points)
				#old_cell.extend(the_cell)
						
				return N
				
		def get_track(self):
				
				cdef int N = self.X.size()
				cdef int i
				
				cdef cnp.ndarray points = np.zeros((N, 3))
				for i in range(N):
						points[i, 0] = self.X[i]
						
						points[i, 1] = self.Y[i]
						
						points[i, 2] = self.Z[i]
				return points
						


		def track(self):
				return list(self.X), list(self.Y), list(self.Z)

		cdef void update_references(self):
				raise RuntimeError("update_references from particle.pyx was "          \
												 + "not overriden by update_references from photon.pyx")


		cdef void _run(self, mixmax_engine *genPTR):
				raise RuntimeError(">>>> _run called from particle")