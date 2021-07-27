#distutils: language = c++

print("Importing .geometry.CSG")


DEF VERBOSE = False
DEF VERBOSE_TALLY = False
DEF DEBUG_MODE = False

from libcpp.vector cimport vector
from libcpp.list cimport list as cpplist;

cdef extern from "IntervalArithmetics.cpp":
	pass


cdef extern from "IntervalArithmetics.h":
	cdef cppclass Interval:
		double t1, t2; 
		Interval();
		Interval(double t1, double t2);

	ctypedef cpplist[Interval] intLIST;
	intLIST intIntersect(intLIST& left, intLIST& right);
	intLIST intPlus(intLIST& left, intLIST& right);
	intLIST intMinus(intLIST& left, intLIST& right);

	cdef cppclass intIterator:
		intIterator();
		intIterator(intLIST crosses);
		void inc();
		double deref();
		double current();



cimport numpy as cnp
import numpy as np 


from libc.string cimport memcpy 
from cython.operator import dereference as deref
from cython.operator import preincrement as inc

from .main cimport Volume

from ..particles.particle cimport STATE

from numpy.math cimport INFINITY as INF

from libcpp.deque cimport deque

from libc.math cimport fmin, fmax, sqrt, cos, sin

from libc.stdlib cimport malloc, free

from . import sdf as plt_geo

cdef double nan = np.nan;
ctypedef BVH Vol
ctypedef BVH V

cdef double eps = .1


from ..types cimport double3


#cdef class Proxy(CSGvol):
#	cdef intLIST crosses
#	cdef cpplist[Interval] it
#	cdef void* vol
#	cdef double particle_position
#
#
#	def __init__(self, vol):
#		self.vol = <void*> vol
#		self.is_proxy = True
#
#	cdef double localSDF(self, STATE& state):
#		self.particle_position += state.last_displacement
#
#
#
#
#	cdef double SDF(self, double3 pos):
#		pass


def lock(msg):
	def _lock(method):
		def new_method(self, *args, **kwargs):
			if self.lock:
				raise RuntimeWarning(msg)
			return method(self, *args, **kwargs)
		return new_method
	return _lock












cdef class BVH(Volume):
	# Workspace
	cdef int Nws
	cdef list tmp_ws
	cdef void** ws
	cdef void** original_ws
	cdef bint has_name


	# Boundary Crossing
	cdef int i                       #position in outers work space
	cdef bint keep                   #which intersected volume will keep its intersections cached for the next iteration

	# user related
	cdef bint lock 					#prevent user from modifying volume after exit code
	cdef str name
	cdef bint render

	# ray marching
	cdef double sdf                  # nearest distance to this volumes surface

	# Ray Tracing
	cdef bint cache                  # is this volume storing cached intersections?
	cdef intIterator cross           # custom c++ iterator for aiding in simulation with cached intersections
	cdef double particle_position;   # must keep track of particles position along the ray 




	def __init__(self):
		self.cache = False
		self.lock = False



	#@lock("Context has been opened more than once.")
	def __enter__(self):
		# Workspace
		self.Nws = 1
		self.tmp_ws = [self]
		return self

	def configure(self, name, render = True):
		self.name = name
		self.render = render

	# Constructing BVH
	#@lock("Modifiying volume after being closed")
	def __contains__(self, other):
		if isinstance(other, BVH):
			self.tmp_ws.append(other)
			self.Nws += 1
			other.setOuter(self, len(self.tmp_ws) - 1)
			return True
		
		cdef double3 pos
		pos.x = other[0]
		pos.y = other[1]
		pos.z = other[2]
		return self.is_inside(pos)

	cpdef setOuter(self, BVH other, int index):
		"""
		other -> outer volume
		index -> self's position in outers workspace
		"""
		self.outer = other
		self.i = index

	def __iter__(self):
		yield from self.tmp_ws[1:]

	def __len__(self):
		return self.Nws

	def set_name(self, str name):
		self.name = name
		self.has_name = True


	# Exit Code
	def __exit__(self, *args, **kwargs):
		self.lock = True

		cdef BVH region
		cdef int i
		
		self.ws = <void**> malloc(self.Nws * sizeof(void*))

		print(self.tmp_ws)
		print(self.Nws)

		for i, region in enumerate(self.tmp_ws):
			self.ws[i]          = <void*> region

		if self.render:
			@plt_geo.sdf3
			def this():
				def SDF(double[:,:] P):
					cdef double3 p
					cdef int N = len(P)
					cdef cnp.ndarray sd = np.zeros(N)
					cdef int i
					for i in range(N):
						p.x = P[i, 0]
						p.y = P[i, 1]
						p.z = P[i, 2]
						sd[i] = self.SDF(p)
					return sd
				return SDF

			generator = this()
			import os
			generator.save(f"geo/{self.name}.stl")

	def get_mesh(self):
		import pyvista as pv
		return pv.read(f"geo/{self.name}.stl")

	def plot(self):
		import pyvista as pv
		mesh = pv.read(f"geo/{self.name}.stl")
		mesh.plot()



			


	cdef bint move(self, STATE& state, double SP):
		raise RuntimeError("'move' called from its virtual in 'Volume.BVH' ")

	cdef void depositUNIFORM(self, STATE& state, double SP):
		raise RuntimeError("depositUNIFORM called from BVH (virtual)")
		print("depositUNIFORM called from BVH (virtual)")
		import time
		time.sleep(10_000)


	cdef void depositDISCRETE(self, STATE& state):
		raise RuntimeError("'depositDISCRETE' called from its virtual in 'Volume.BVH' ")
	
	cdef void depositLOCAL(self, double3& pos, double E):
		raise RuntimeError("'depositLOCAL' called from its virtual in 'Volume.BVH' ")

	cdef void depositRANDOM(self, STATE& state, double E, double tau):
		raise RuntimeError("'depositRANDOM' called from its virtual in 'Volume.BVH' ")

	cdef double main_intersect(self, STATE& state):
		raise RuntimeError("'main_intersect' called from its virtual in 'Volume.BVH' ")


	cdef void localSDF(self, STATE& state):
		raise RuntimeError("'localSDF' called from its virtual in 'Volume.BVH' ")


	cdef void* searchO(self, STATE& state):
		cdef int i

		for i in range(1, self.Nws):
			if self.ws[i] == state.current_region: continue

			if (<BVH> self.ws[i]).is_inside(state.pos):
				print(i)
				return self.ws[i]
		print(0)
		return <void*> self

		#for i in range(1, (<BVH> state.current_region).i):
		#	if (<BVH> self.ws[i]).is_inside(state.pos):
		#		return self.ws[i]
#
		#for i in range((<BVH> state.current_region).i+1, self.Nws):
		#	if (<BVH> self.ws[i]).is_inside(state.pos):
		#		return self.ws[i]


	cdef double SDF(self, double3 pos):
		raise RuntimeError("'SDF' called from its virtual in 'Volume.BVH' ")

	cdef bint is_inside(self, double3 pos):
		raise RuntimeError("'is_inside' called from its virtual in 'Volume.BVH' ")

	cdef void exit(self):
		cdef int i
		for i in range(self.Nws):
			#if (<BVH> self.ws[i]).keep: continue
			(<BVH> self.ws[i]).cache = False























#   ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄ 
#  ▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌
#  ▐░█▀▀▀▀▀▀▀▀▀ ▐░█▀▀▀▀▀▀▀▀▀ ▐░█▀▀▀▀▀▀▀▀▀ 
#  ▐░▌          ▐░▌          ▐░▌          
#  ▐░▌          ▐░█▄▄▄▄▄▄▄▄▄ ▐░▌ ▄▄▄▄▄▄▄▄ 
#  ▐░▌          ▐░░░░░░░░░░░▌▐░▌▐░░░░░░░░▌
#  ▐░▌           ▀▀▀▀▀▀▀▀▀█░▌▐░▌ ▀▀▀▀▀▀█░▌
#  ▐░▌                    ▐░▌▐░▌       ▐░▌
#  ▐░█▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄█░▌▐░█▄▄▄▄▄▄▄█░▌
#  ▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌
#   ▀▀▀▀▀▀▀▀▀▀▀  ▀▀▀▀▀▀▀▀▀▀▀  ▀▀▀▀▀▀▀▀▀▀▀ 
#                                         


cdef class CSGvol(BVH):
	# Ray Marching
	cdef double global_sdf           # nearest distance to nearest surface contained by and including this volumes surface 
	cdef int i0                      # position of closest volume in selfs workspace


	cdef double E

	# CONSTRUCTING A VOLUME
	cdef bint is_inside(self, double3& pos):
		raise RuntimeError("asdsadsad")

	cdef void depositUNIFORM(self, STATE& state, double SP):
		pass

	cdef void depositLOCAL(self, double3& pos, double E):
		pass


	def __init__(self):
		# Opening lock, volume can be modified
		super(CSGvol, self).__init__()
		self.E = 0


	cdef void depositRANDOM(self, STATE& state, double E, double tau):
		pass


	cdef void depositLocaly(self, double3& pos, double E):
		pass

	#@lock("Modifiying volume after being closed")
	def rotate(self, axis, angle):
		return NotImplemented

	#@lock("Modifiying volume after being closed")
	def translate(self, direction, displacement):
		return NotImplemented

	#@lock("Modifiying volume after being closed")
	def __or__(self, other):
		return Union(self, other)

	#@lock("Modifiying volume after being closed")
	def __add__(self, other):
		return Union(self, other)

	#@lock("Modifiying volume after being closed")
	def __and__(self, CSGvol other):
		return Intersection(self, other)

	#@lock("Modifiying volume after being closed")
	def __sub__(self, CSGvol other):
		return Subtraction(self, other)


	cdef intLIST intersect(self, double3& pos, double3& dire):
		raise RuntimeError(".intersect called from virtual")

	


	cdef double SDF(self, double3 pos):
		raise RuntimeError("SDF FROM VOL WAS CALLED")

	cdef void localSDF(self, STATE& state):

		if self.cache:
			IF DEBUG_MODE:
				print(f"Volume[{self.i}] has cached intersections...")
				input(f"Volume[{self.i}].particle_position: {self.particle_position}")
				input(f"state.last_displacement: {state.last_displacement}")

			self.particle_position += state.last_displacement

			IF DEBUG_MODE:
				print(f"Updated particle position in Volume[{self.i}]...")
				input(f"Volume[{self.i}].particle_position: {self.particle_position}")

			self.sdf = self.cross.current() - self.particle_position

			IF DEBUG_MODE:
				input(f"Safest distance to Volume[{self.i}]: {self.sdf}")

			return


		IF DEBUG_MODE: print(f"Volume[{self.i}] does not have cached intersections...")

		self.sdf = self.SDF(state.pos)

		IF DEBUG_MODE: input(f"Safest distance to Volume[{self.i}]: {self.sdf}")
		

	cdef double main_intersect(self, STATE& state):

		IF VERBOSE: print("MAIN_INTERSECT")
		IF VERBOSE: print(self.cache)

		if self.cache:
			return self.cross.current() - self.particle_position

		self.particle_position = 0

		cdef intLIST temp = self.intersect(state.pos, state.dire)
		cdef intIterator temp2 = intIterator(temp)

		IF VERBOSE: print(temp.size())

		cdef int i
		IF VERBOSE:
			print("PRINTING INTERSECTIONS")
			for i in range(temp.size()):
				print(temp2.current())
				temp2.inc()
				print(temp2.current())
				temp2.inc()

		self.cross = intIterator(temp)
		self.cache = True
		return self.cross.current()

	cdef void globalSDF(self, STATE& state):
		IF DEBUG_MODE: input("-- CALCULATING SAFEST DISTANCE")

		IF DEBUG_MODE: input(f"cache[0] = {self.cache}")

		if self.cache:
			IF DEBUG_MODE:
				print(f"Volume[0] has cached intersections...")
				input(f"Volume[0].particle_position: {self.particle_position}")
				input(f"state.last_displacement: {state.last_displacement}")

			self.particle_position += state.last_displacement

			IF DEBUG_MODE:
				print("Updated particle position in Volume[0]...")
				input(f"Volume[0].particle_position: {self.particle_position}")

			self.sdf = self.cross.current() - self.particle_position

			IF DEBUG_MODE:
				print("Updated particle position in Volume[0]...")
				input(f"Safest distance to Volume[0]: {self.sdf}")

		else:
			IF DEBUG_MODE:
				print(f"Volume[0] does not have cached intersections...")

			self.sdf = -self.SDF(state.pos)

			IF DEBUG_MODE:
				input(f"Safest distance to Volume[0]: {self.sdf}")

		self.i0 = 0
		self.global_sdf = self.sdf

		cdef int i
		for i in range(1, self.Nws):
			(<V> self.ws[i]).localSDF(state)
			if (<V> self.ws[i]).sdf < self.global_sdf:
				self.global_sdf = (<V> self.ws[i]).sdf
				self.i0 = i

		IF DEBUG_MODE: input(f"The closest volume is: Volume[{self.i0}] @ {(<V> self.ws[self.i0]).sdf}cm")




	cdef inline void final(self, STATE& state):
		IF DEBUG_MODE: input(f"FINAL DISPLACEMENT: L = {state.L}")

		state.last_displacement = state.L

		state.pos.x += state.dire.x*state.L
		state.pos.y += state.dire.y*state.L
		state.pos.z += state.dire.z*state.L

		state.L = 0

		IF VERBOSE:
			for i in range(self.Nws):
				print(f"is_inside[{i}] = {(<V> self.ws[i]).is_inside(state.pos)}")



	cdef inline void virtual_event(self, STATE& state, double dr):
		IF DEBUG_MODE: input(f"VIRTUAL: dr = {dr}, L = {state.L}")
		state.last_displacement = dr


		IF VERBOSE: print("norm of dire:", state.dire.x**2 + state.dire.y**2 + state.dire.z**2)
		state.pos.x += state.dire.x*dr
		state.pos.y += state.dire.y*dr
		state.pos.z += state.dire.z*dr

		state.L -= dr

		IF VERBOSE:
			for i in range(self.Nws):
				print(f"is_inside[{i}] = {(<V> self.ws[i]).is_inside(state.pos)}")



	cdef int intEVENT(self, STATE& state):
		"""
		Instruct nearest surface to intersect with the ray defined
		in the particles state:

			ray = state.pos + t*ray.dire  for t >= 0

		Act according to the number of intersections found:

			0 intersections
			---------------
			Particle hits L before hitting second nearest surface?

				YES:
					- Perform final displacement event.


				NO: 
					- Perform virtual displacement event;
					- Substitute nearest surface by a point at
					  infinity.


			N intersections
			---------------
			- insert point at infinity in the intersecions map

			We have three quantities:
				L: distance to end of step
				second_nearest: distance to second nearest surface
				d: distance to intersection

				if min() == L:
					go to L; end step;

				elif min() == second_nearest:
					virtual step
					make proxy
					continue

				elif min() == d:
					go to d, finish loop



		"""
		IF DEBUG_MODE: input("INTERSECTION EVENT")
		# get the closest volume, easier to write this way...
		cdef void* vol = self.ws[self.i0]
		cdef double cross = (<V> vol).main_intersect(state)

		

		# look for distance to the second nearest surface
		cdef double second_nearest = INF
		
		cdef int i
		print(state.pos)

		for i in range(0, self.i0):
			print(i, (<V> self.ws[i]).sdf, (<V> self.ws[i]))
			if (<V> self.ws[i]).sdf < second_nearest:
				second_nearest = (<V> self.ws[i]).sdf

		for i in range(self.i0+1, self.Nws):
			print(i, (<V> self.ws[i]).sdf, (<V> self.ws[i]))
			if (<V> self.ws[i]).sdf < second_nearest:
				second_nearest = (<V> self.ws[i]).sdf



		IF VERBOSE: print(f"cross = {cross} | second_nearest = {second_nearest}")

		# no intersections found
		if cross == INF:
			IF VERBOSE: print("no intersections found")
			# min() = L
			if second_nearest > state.L:
				self.final(state)
				self.exit()
				return 0
			# min() = second_nearest
			self.virtual_event(state, second_nearest)
			#self.ws[self.i0] = proxy
			#self.reset = True
			return 1


		if cross < second_nearest:

			# min() == L
			if cross > state.L:
				IF VERBOSE: print("min() == L")
				self.final(state)
				self.exit()
				return 0

			# min() == it.first 
			IF VERBOSE: print("particle will definetly hit the surface")
			# particle will definetly hit the surface
			self.virtual_event(state, cross)

			IF VERBOSE: print(f"before incrementing: current = {(<V> vol).cross.current()}")
			(<V> vol).cross.inc()
			IF VERBOSE: print("icremented successfully")
			IF VERBOSE: print(f"after incrementing: current = {(<V> vol).cross.current()}")
			self.boundary_crossing(state)
			return 2

		# min() == L
		if second_nearest > state.L:
			IF VERBOSE: print("min() == L 222")
			self.final(state)
			self.exit()
			return 0



	cdef void exitINNER_TO_OUTER(self):
		cdef int i
		for i in range(1, self.Nws):
			(<V> self.ws[i]).cache = False
			(<V> self.ws[i]).keep = False


	cdef void exitINNER_TO_INNER(self):
		cdef int i
		for i in range(0, self.Nws):
			(<V> self.ws[i]).cache = False
			(<V> self.ws[i]).keep = False

	cdef void exitOUTER_TO_INNER(self):
		cdef int i
		for i in range(0, self.Nws):
			(<V> self.ws[i]).cache = False
			(<V> self.ws[i]).keep = False

	cdef void boundary_crossing(self, STATE& state):

		# from inner to somewhere in outer
		if self.i0 == 0:
			state.current_region = (<V> self.outer).searchO(state)

			# staying in outer, must keep cached intersections
			if state.current_region == <void*> self.outer:
				self.keep = True
				self.exitINNER_TO_OUTER()
				return

			# entering some adjacent volume, must intersect it then
			(<V> state.current_region).main_intersect(state)
			(<V> state.current_region).keep = True
			self.exitINNER_TO_INNER()
			return 


		# from outer to inner
		state.current_region = self.ws[self.i0]
		self.exitOUTER_TO_INNER()
		(<V> state.current_region).keep = True
		(<V> state.current_region).cache = True
		return


	cdef bint move(self, STATE& state, double SP):
		cdef int case
		cdef int i
		IF DEBUG_MODE: input("\n ----MOVING PARTICLE----")
		IF DEBUG_MODE: input(f"CURRENT_POSITION: {state.pos}")
		IF DEBUG_MODE: input(f"Is inside current region? {(<V> state.current_region).is_inside(state.pos)}")
		IF DEBUG_MODE:
			for i in range(self.Nws):
				print(f"is_inside[{i}] = {(<V> self.ws[i]).is_inside(state.pos)}")

		IF DEBUG_MODE: input("STARTING EVENT LOOP")
		while True:
			self.globalSDF(state)


			IF DEBUG_MODE: print(f"The safest distance is {self.global_sdf}cm | Physics proposed {state.L}cm ")

			if self.global_sdf > state.L:
				self.final(state)
				
				for i in range(self.Nws):
					(<BVH> self.ws[i]).keep = False
					(<BVH> self.ws[i]).cache = False
				self.exit()


				return False

			if self.global_sdf < eps:
				case = self.intEVENT(state)


				if case == 2: # boundary crossing
					

					return True

				if case == 0: # final displacement

					return False

				if case == 1: # virtual displacement, proxys have been set
					continue
			self.virtual_event(state, self.global_sdf)





























































																															

#   ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄        ▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄       ▄▄ 
#  ▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░▌      ▐░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░▌     ▐░░▌
#   ▀▀▀▀█░█▀▀▀▀ ▐░█▀▀▀▀▀▀▀█░▌▐░█▀▀▀▀▀▀▀█░▌▐░▌░▌     ▐░▌▐░█▀▀▀▀▀▀▀▀▀ ▐░█▀▀▀▀▀▀▀▀▀ ▐░█▀▀▀▀▀▀▀█░▌▐░█▀▀▀▀▀▀▀█░▌▐░▌░▌   ▐░▐░▌
#       ▐░▌     ▐░▌       ▐░▌▐░▌       ▐░▌▐░▌▐░▌    ▐░▌▐░▌          ▐░▌          ▐░▌       ▐░▌▐░▌       ▐░▌▐░▌▐░▌ ▐░▌▐░▌
#       ▐░▌     ▐░█▄▄▄▄▄▄▄█░▌▐░█▄▄▄▄▄▄▄█░▌▐░▌ ▐░▌   ▐░▌▐░█▄▄▄▄▄▄▄▄▄ ▐░█▄▄▄▄▄▄▄▄▄ ▐░▌       ▐░▌▐░█▄▄▄▄▄▄▄█░▌▐░▌ ▐░▐░▌ ▐░▌
#       ▐░▌     ▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░▌  ▐░▌  ▐░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░▌       ▐░▌▐░░░░░░░░░░░▌▐░▌  ▐░▌  ▐░▌
#       ▐░▌     ▐░█▀▀▀▀█░█▀▀ ▐░█▀▀▀▀▀▀▀█░▌▐░▌   ▐░▌ ▐░▌ ▀▀▀▀▀▀▀▀▀█░▌▐░█▀▀▀▀▀▀▀▀▀ ▐░▌       ▐░▌▐░█▀▀▀▀█░█▀▀ ▐░▌   ▀   ▐░▌
#       ▐░▌     ▐░▌     ▐░▌  ▐░▌       ▐░▌▐░▌    ▐░▌▐░▌          ▐░▌▐░▌          ▐░▌       ▐░▌▐░▌     ▐░▌  ▐░▌       ▐░▌
#       ▐░▌     ▐░▌      ▐░▌ ▐░▌       ▐░▌▐░▌     ▐░▐░▌ ▄▄▄▄▄▄▄▄▄█░▌▐░▌          ▐░█▄▄▄▄▄▄▄█░▌▐░▌      ▐░▌ ▐░▌       ▐░▌
#       ▐░▌     ▐░▌       ▐░▌▐░▌       ▐░▌▐░▌      ▐░░▌▐░░░░░░░░░░░▌▐░▌          ▐░░░░░░░░░░░▌▐░▌       ▐░▌▐░▌       ▐░▌
#        ▀       ▀         ▀  ▀         ▀  ▀        ▀▀  ▀▀▀▀▀▀▀▀▀▀▀  ▀            ▀▀▀▀▀▀▀▀▀▀▀  ▀         ▀  ▀         ▀ 



cdef cnp.ndarray new_rotationT(_axis, angle):

		cdef cnp.ndarray axis = np.array(_axis)

		axis = axis/np.sqrt(np.sum(axis**2))

		cdef double ux = axis[0]
		cdef double uy = axis[1]
		cdef double uz = axis[2]


		cdef cnp.ndarray T = np.zeros((4,4))

		cdef double _cos = cos(angle)
		cdef double oneMcos = 1 - _cos
		cdef double _sin = sqrt(1 - _cos*_cos)

		T[0, 0] = _cos + ux*ux*oneMcos
		T[0, 1] = ux*uy*oneMcos - uz*_sin
		T[0, 2] = ux*uz*oneMcos + uy*_sin

		T[1, 0] = uy*ux*oneMcos + uz*_sin
		T[1, 1] = _cos + uy*uy*oneMcos
		T[1, 2] = uy*uz*oneMcos - ux*_sin

		T[2, 0] = uz*ux*oneMcos-uy*_sin
		T[2, 1] = uz*uy*oneMcos + ux*_sin
		T[2, 2] = _cos + uz*uz*oneMcos

		T[3, 3] = 1

		return T

cdef cnp.ndarray Carr_to_NParr(double* arr):
	numbers = np.zeros(16)
	cdef int i
	for i in range(16):
		numbers[i] = arr[i]

	numbers.shape = (4, 4)
	return numbers

cdef class Transform(CSGvol):
	cdef Primitive primitive
	cdef double[16] T, iT

	def __init__(self, Primitive primitive, cnp.ndarray T, cnp.ndarray iT):
		self.primitive = primitive

		cdef int i
		cdef double t, it
		for i, (t, it) in enumerate(zip(T.flat, iT.flat)):
			self.T[i] = t
			self.iT[i] = it

	@property
	def matrix(self):
		return Carr_to_NParr(self.T)

	@property
	def inv_matrix(self):
		return Carr_to_NParr(self.iT)

	def translate(self, double dx, double dy, double dz):
		self.T[3]  += dx
		self.T[7]  += dy
		self.T[11] += dz

		self.iT[3]  -= dx
		self.iT[7]  -= dy
		self.iT[11] -= dz
		return self




	def rotate(self, axis, angle):
		cdef cnp.ndarray nT = new_rotationT(axis, angle)
		cdef cnp.ndarray T = Carr_to_NParr(self.T)

		nT = nT@T

		cdef cnp.ndarray inT = np.linalg.inv(nT)

		cdef int i
		cdef double t, it

		for i, (t, it) in enumerate(zip(nT.flat, inT.flat)):
			self.T[i] = t
			self.iT[i] = it

		return self



	cdef void inv_pos(self, double3& rpos):
		cdef double3 pos = rpos
		rpos.x = self.iT[0]*pos.x + self.iT[1]*pos.y + self.iT[2] *pos.z  + self.iT[3]
		rpos.y = self.iT[4]*pos.x + self.iT[5]*pos.y + self.iT[6] *pos.z  + self.iT[7]
		rpos.z = self.iT[8]*pos.x + self.iT[9]*pos.y + self.iT[10]*pos.z  + self.iT[11]

	cdef void inv_dire(self, double3& rdire):
		cdef double3 dire = rdire
		rdire.x = self.iT[0]*dire.x + self.iT[1]*dire.y + self.iT[2] *dire.z  + self.iT[3]
		rdire.y = self.iT[4]*dire.x + self.iT[5]*dire.y + self.iT[6] *dire.z  + self.iT[7]
		rdire.z = self.iT[8]*dire.x + self.iT[9]*dire.y + self.iT[10]*dire.z  + self.iT[11]




	cdef intLIST intersect(self, double3& pos, double3& dire):
		cdef double3 rpos
		rpos.x = self.iT[0]*pos.x + self.iT[1]*pos.y + self.iT[2] *pos.z  + self.iT[3]
		rpos.y = self.iT[4]*pos.x + self.iT[5]*pos.y + self.iT[6] *pos.z  + self.iT[7]
		rpos.z = self.iT[8]*pos.x + self.iT[9]*pos.y + self.iT[10]*pos.z  + self.iT[11]

		cdef double3 rdire
		rdire.x = self.iT[0]*dire.x + self.iT[1]*dire.y + self.iT[2] *dire.z  + self.iT[3]
		rdire.y = self.iT[4]*dire.x + self.iT[5]*dire.y + self.iT[6] *dire.z  + self.iT[7]
		rdire.z = self.iT[8]*dire.x + self.iT[9]*dire.y + self.iT[10]*dire.z  + self.iT[11]

		return self.primitive.intersect(rpos, rdire)

	cdef bint is_inside(self, double3& pos):
		cdef double3 rpos
		rpos.x = self.iT[0]*pos.x + self.iT[1]*pos.y + self.iT[2] *pos.z  + self.iT[3]
		rpos.y = self.iT[4]*pos.x + self.iT[5]*pos.y + self.iT[6] *pos.z  + self.iT[7]
		rpos.z = self.iT[8]*pos.x + self.iT[9]*pos.y + self.iT[10]*pos.z  + self.iT[11]

		return self.primitive.is_inside(rpos)


cdef class Isometry(Transform):


	def __init__(self, Primitive primitive, cnp.ndarray T, cnp.ndarray iT):
		self.primitive = primitive

		cdef int i
		cdef double t, it
		for i, (t, it) in enumerate(zip(T.flat, iT.flat)):
			self.T[i] = t
			self.iT[i] = it


	cdef void inv_pos(self, double3& rpos):
		IF VERBOSE: print("Isometry.inv_pos")
		cdef double3 pos = rpos
		rpos.x = self.iT[0]*pos.x + self.iT[1]*pos.y + self.iT[2]*pos.z + self.iT[3]
		rpos.y = self.iT[4]*pos.x + self.iT[5]*pos.y + self.iT[6]*pos.z + self.iT[7]
		rpos.z = self.iT[8]*pos.x + self.iT[9]*pos.y + self.iT[10]*pos.z + self.iT[11]





	cdef void inv_dire(self, double3& pos):
		IF VERBOSE: print("Isometry.inv_dire")

		cdef double3 tmp_pos = pos
		pos.x = self.iT[0]*tmp_pos.x + self.iT[1]*tmp_pos.y + self.iT[2]*tmp_pos.z 
		pos.y = self.iT[4]*tmp_pos.x + self.iT[5]*tmp_pos.y + self.iT[6]*tmp_pos.z 
		pos.z = self.iT[8]*tmp_pos.x + self.iT[9]*tmp_pos.y + self.iT[10]*tmp_pos.z



	cdef double SDF(self, double3& pos):
		cdef double3 rpos
		rpos.x = self.iT[0]*pos.x + self.iT[1]*pos.y + self.iT[2] *pos.z  + self.iT[3]
		rpos.y = self.iT[4]*pos.x + self.iT[5]*pos.y + self.iT[6] *pos.z  + self.iT[7]
		rpos.z = self.iT[8]*pos.x + self.iT[9]*pos.y + self.iT[10]*pos.z  + self.iT[11]
		return self.primitive.SDF(rpos)






cdef class Identity(Isometry):
	
	def __init__(self, Primitive primitive):
		self.primitive = primitive
		cdef int i
		for i in range(16):
			self.T[i] = 0
			self.iT[i] = 0

		self.T[0] = 1
		self.T[5] = 1
		self.T[10] = 1
		self.T[15] = 1

		self.iT[0] = 1
		self.iT[5] = 1
		self.iT[10] = 1
		self.iT[15] = 1




	cdef void inv_pos(self, double3& pos):
		pass

	cdef void inv_dire(self, double3& dire):
		pass


	def translate(self, dx, dy, dz):
		return Translation(self.primitive, dx, dy, dz)

	def rotate(self, axis, angle):
		return Rotation(self.primitive, axis, angle)

	def scale(self, s):
		self.primitive.scale(s)

	cdef intLIST intersect(self, double3& pos, double3& dire):
		return self.primitive.intersect(pos, dire)

	cdef double SDF(self, double3& pos):
		return self.primitive.SDF(pos)

	#def matrix(self):
	#	cdef cnp.ndarray m = np.zeros((4, 4))
	#	cdef int i
	#	for i in range(4):
	#		m[i, i] = 1
	#	return m


cdef class NonIsometry(Transform):
	pass





cdef class Translation(Isometry):

	def __init__(self, Primitive primitive, dx, dy, dz):

		self.primitive = primitive

		self.T[3]  = dx
		self.T[7]  = dy
		self.T[11] = dz

		self.iT[3]  = -dx
		self.iT[7]  = -dy
		self.iT[11] = -dz

		self.T[0] = 1
		self.T[5] = 1
		self.T[10] = 1
		self.T[15] = 1

		self.iT[0] = 1
		self.iT[5] = 1
		self.iT[10] = 1
		self.iT[15] = 1


	cdef void inv_pos(self, double3& pos):
		pos.x += self.iT[3]
		pos.y += self.iT[7]
		pos.z += self.iT[11]

	cdef void inv_dire(self, double3& pos):
		pass

	cdef double SDF(self, double3 pos):
		pos.x += self.iT[3]
		pos.y += self.iT[7]
		pos.z += self.iT[11]

		return self.primitive.SDF(pos)

	cdef intLIST intersect(self, double3 pos, double3& dire):
		pos.x += self.iT[3]
		pos.y += self.iT[7]
		pos.z += self.iT[11]

		return self.primitive.intersect(pos, dire)


	def rotate(self, axis, angle):
		"""
		R * T yields a general isometry.
		"""

		cdef cnp.ndarray rot = new_rotationT(axis, angle)
		cdef cnp.ndarray T = rot @ self.matrix
		cdef cnp.ndarray iT = np.linalg.inv(T) 

		return Isometry(self.primitive, T, iT)







cdef class Rotation(Isometry):

	def __init__(self,Primitive primitive, axis, angle):
		self.primitive = primitive



		cdef cnp.ndarray T = new_rotationT(axis, angle)
		cdef cnp.ndarray iT = np.linalg.inv(T)

		cdef int i
		cdef double t, it
		for i, (t, it) in enumerate(zip(T.flat, iT.flat)):
			self.T[i] = t
			self.iT[i] = it


	cdef void inv_pos(self, double3& rpos):
		cdef double3 pos = rpos
		rpos.x = self.iT[0]*pos.x + self.iT[1]*pos.y + self.iT[2]*pos.z
		rpos.y = self.iT[4]*pos.x + self.iT[5]*pos.y + self.iT[6]*pos.z
		rpos.z = self.iT[8]*pos.x + self.iT[9]*pos.y + self.iT[10]*pos.z





	cdef void inv_dire(self, double3& pos):
		cdef double3 tmp_pos = pos
		pos.x = self.iT[0]*tmp_pos.x + self.iT[1]*tmp_pos.y + self.iT[2]*tmp_pos.z
		pos.y = self.iT[4]*tmp_pos.x + self.iT[5]*tmp_pos.y + self.iT[6]*tmp_pos.z
		pos.z = self.iT[8]*tmp_pos.x + self.iT[9]*tmp_pos.y + self.iT[10]*tmp_pos.z




	def translate(self, dx, dy, dz):
		cdef cnp.ndarray T = np.zeros((4, 4))

		T[0, 0] = self.T[0]
		T[0, 1] = self.T[1]
		T[0, 2] = self.T[2]
		T[0, 3] = dx

		T[1, 0] = self.T[3]
		T[1, 1] = self.T[4]
		T[1, 2] = self.T[5]
		T[1, 3] = dy


		T[2, 0] = self.T[6]
		T[2, 1] = self.T[7]
		T[2, 2] = self.T[8]
		T[2, 3] = dz

		T[3, :] = np.array([0, 0, 0, 1])

		cdef cnp.ndarray iT = np.linalg.inv(T)

		return Isometry(self, T, iT)

	def rotate(self, axis, angle):
		cdef cnp.ndarray rot = new_rotationT(axis, angle)
		cdef cnp.ndarray T = Carr_to_NParr(self.T)

		T = rot @ T
		iT = np.linalg.inv(T)

		cdef int i
		cdef double t, it
		for i, (t, it) in enumerate(zip(T.flat, iT.flat)):
			self.T[i] = t
			self.iT[i] = it

		return self












































#   ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄        ▄  ▄▄▄▄▄▄▄▄▄▄▄ 
#  ▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░▌      ▐░▌▐░░░░░░░░░░░▌
#  ▐░█▀▀▀▀▀▀▀█░▌▐░█▀▀▀▀▀▀▀█░▌▐░█▀▀▀▀▀▀▀▀▀ ▐░█▀▀▀▀▀▀▀█░▌▐░█▀▀▀▀▀▀▀█░▌ ▀▀▀▀█░█▀▀▀▀  ▀▀▀▀█░█▀▀▀▀ ▐░█▀▀▀▀▀▀▀█░▌▐░▌░▌     ▐░▌▐░█▀▀▀▀▀▀▀▀▀ 
#  ▐░▌       ▐░▌▐░▌       ▐░▌▐░▌          ▐░▌       ▐░▌▐░▌       ▐░▌     ▐░▌          ▐░▌     ▐░▌       ▐░▌▐░▌▐░▌    ▐░▌▐░▌          
#  ▐░▌       ▐░▌▐░█▄▄▄▄▄▄▄█░▌▐░█▄▄▄▄▄▄▄▄▄ ▐░█▄▄▄▄▄▄▄█░▌▐░█▄▄▄▄▄▄▄█░▌     ▐░▌          ▐░▌     ▐░▌       ▐░▌▐░▌ ▐░▌   ▐░▌▐░█▄▄▄▄▄▄▄▄▄ 
#  ▐░▌       ▐░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌     ▐░▌          ▐░▌     ▐░▌       ▐░▌▐░▌  ▐░▌  ▐░▌▐░░░░░░░░░░░▌
#  ▐░▌       ▐░▌▐░█▀▀▀▀▀▀▀▀▀ ▐░█▀▀▀▀▀▀▀▀▀ ▐░█▀▀▀▀█░█▀▀ ▐░█▀▀▀▀▀▀▀█░▌     ▐░▌          ▐░▌     ▐░▌       ▐░▌▐░▌   ▐░▌ ▐░▌ ▀▀▀▀▀▀▀▀▀█░▌
#  ▐░▌       ▐░▌▐░▌          ▐░▌          ▐░▌     ▐░▌  ▐░▌       ▐░▌     ▐░▌          ▐░▌     ▐░▌       ▐░▌▐░▌    ▐░▌▐░▌          ▐░▌
#  ▐░█▄▄▄▄▄▄▄█░▌▐░▌          ▐░█▄▄▄▄▄▄▄▄▄ ▐░▌      ▐░▌ ▐░▌       ▐░▌     ▐░▌      ▄▄▄▄█░█▄▄▄▄ ▐░█▄▄▄▄▄▄▄█░▌▐░▌     ▐░▐░▌ ▄▄▄▄▄▄▄▄▄█░▌
#  ▐░░░░░░░░░░░▌▐░▌          ▐░░░░░░░░░░░▌▐░▌       ▐░▌▐░▌       ▐░▌     ▐░▌     ▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░▌      ▐░░▌▐░░░░░░░░░░░▌
#   ▀▀▀▀▀▀▀▀▀▀▀  ▀            ▀▀▀▀▀▀▀▀▀▀▀  ▀         ▀  ▀         ▀       ▀       ▀▀▀▀▀▀▀▀▀▀▀  ▀▀▀▀▀▀▀▀▀▀▀  ▀        ▀▀  ▀▀▀▀▀▀▀▀▀▀▀ 


cdef double delta = 1e-10


cdef class CSGop(CSGvol):
	cdef CSGvol R, L
	cdef double (*rule)(double, double)

	def __init__(self, CSGvol L, CSGvol R):
		self.L = L
		self.R = R

	cdef bint is_inside(self, double3 pos):
		raise RuntimeError("'is_inside' called from virtual Volume.BVH.CSGvol.CSGop")

	def translate(self, dx, dy, dz):
		if isinstance(self.L, Transform):
			self.L.primitive.translate(dx, dy, dz)
		else: self.L.translate(dx, dy, dz)

		if isinstance(self.R, Transform):
			self.R.primitive.translate(dx, dy, dz)
		else: self.R.translate(dx, dy, dz)


	def rotate(self, axis, angle):
		self.L.rotate(axis, angle)
		self.R.rotate(axis, angle)

	cdef intLIST intersect(self, double3& pos, double3& dire):
		raise RuntimeError("Called from virtual;")
		
		#return self.rule(self.L.intersect(pos, dire), self.R.intersect(pos, dire))



cdef class Subtraction(CSGop):

	"""
	SDF:  max(-SDF_L, SDF_R)
		
	RAY TRACE:
		
		   L/R     IN    OUT     BORDER  
		 -------- ---- -------- -------- 
		  IN       OUT      IN       BORDER      
		  OUT      OUT     OUT      OUT  
		  BORDER   OUT  BORDER      
		  
		  1 = out
		  -1 = in
		  0 = border
			
	"""
	def __init__(self, CSGvol L, CSGvol R):
		# child nodes
		super(Subtraction, self).__init__(L, R)
		#self.mesh = L.mesh - R.mesh




	cdef double SDF(self, double3& pos):
		return fmax(self.L.SDF(pos), -self.R.SDF(pos))

	cdef bint is_inside(self, double3& pos):
		return self.L.is_inside(pos) and (not self.R.is_inside(pos))

	def __repr__(self):
		return "<Subtraction>"

	cdef intLIST intersect(self, double3& pos, double3& dire):
		IF VERBOSE: print("SUBTRACTING: \n --left-- \n")

		cdef intLIST L = self.L.intersect(pos, dire)

		if L.size() == 0:
			return L

		IF VERBOSE: print("\n --right-- \n")

		cdef intLIST R = self.R.intersect(pos, dire)

		if R.size() == 0:
			return L

		return intMinus(L, R)








			   
cdef class Union(CSGop):

	"""
	SDF:
		min(SDF_a, SDF_b)
		
	RAY TRACE:
		
		   L/R     IN    OUT     BORDER  
		 -------- ---- -------- -------- 
		  IN       IN   IN       IN      
		  OUT      IN   OUT      BORDER  
		  BORDER   IN   BORDER      
		  
		  1 = out
		  -1 = in
		  0 = border
			
	"""
	def __init__(self, CSGvol L, CSGvol R):
		# child nodes
		super(Union, self).__init__(L, R)
		#self.mesh = L.mesh + R.mesh

	def __repr__(self):
		return "<Union>"

	cdef double SDF(self, double3& pos):
		return fmin(self.L.SDF(pos), self.R.SDF(pos))

	cdef intLIST intersect(self, double3& pos, double3& dire):
		IF VERBOSE: print("UNION: \n --left-- \n")

		cdef intLIST L = self.L.intersect(pos, dire)

		if L.size() == 0:
			return self.R.intersect(pos, dire)

		IF VERBOSE: print("\n --right-- \n")

		cdef intLIST R = self.R.intersect(pos, dire)

		if R.size() == 0:
			return L

		return intPlus(L, R)

	cdef bint is_inside(self, double3& pos):
		return self.L.is_inside(pos) or self.R.is_inside(pos)


cdef class Intersection(CSGop):

	"""
	SDF:
		max(SDF_a, SDF_b)
		
	RAY TRACE:
		
		   L/R     IN    OUT     BORDER  
		 -------- ---- -------- -------- 
		  IN       IN      OUT       BORDER      
		  OUT      OUT     OUT      OUT  
		  BORDER   BORDER  OUT      
		  
		  1 = out
		  -1 = in
		  0 = border
			
	"""
	def __init__(self, CSGvol L, CSGvol R):
		# child nodes
		super(Intersection, self).__init__(L, R)
		#self.mesh = L.mesh.boolean_cut(R.mesh)


	def __repr__(self):
		return "<Intersection>"
 
	cdef double SDF(self, double3& pos):
		return fmax(self.L.SDF(pos), self.R.SDF(pos))

	cdef intLIST intersect(self, double3& pos, double3& dire):
		IF VERBOSE: print("INTERSECTING: ")
		cdef intLIST L = self.L.intersect(pos, dire)

		if L.size() == 0:
			return L

		cdef intLIST R = self.R.intersect(pos, dire)
		if R.size() == 0:
			return R

		return intIntersect(L, R)

	cdef bint is_inside(self, double3& pos):
		return self.L.is_inside(pos) and self.R.is_inside(pos)



















































#   ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄       ▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄               ▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄ 
#  ▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░▌     ▐░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░▌             ▐░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌
#  ▐░█▀▀▀▀▀▀▀█░▌▐░█▀▀▀▀▀▀▀█░▌ ▀▀▀▀█░█▀▀▀▀ ▐░▌░▌   ▐░▐░▌ ▀▀▀▀█░█▀▀▀▀  ▀▀▀▀█░█▀▀▀▀  ▀▀▀▀█░█▀▀▀▀  ▐░▌           ▐░▌ ▐░█▀▀▀▀▀▀▀▀▀ ▐░█▀▀▀▀▀▀▀▀▀ 
#  ▐░▌       ▐░▌▐░▌       ▐░▌     ▐░▌     ▐░▌▐░▌ ▐░▌▐░▌     ▐░▌          ▐░▌          ▐░▌       ▐░▌         ▐░▌  ▐░▌          ▐░▌          
#  ▐░█▄▄▄▄▄▄▄█░▌▐░█▄▄▄▄▄▄▄█░▌     ▐░▌     ▐░▌ ▐░▐░▌ ▐░▌     ▐░▌          ▐░▌          ▐░▌        ▐░▌       ▐░▌   ▐░█▄▄▄▄▄▄▄▄▄ ▐░█▄▄▄▄▄▄▄▄▄ 
#  ▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌     ▐░▌     ▐░▌  ▐░▌  ▐░▌     ▐░▌          ▐░▌          ▐░▌         ▐░▌     ▐░▌    ▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌
#  ▐░█▀▀▀▀▀▀▀▀▀ ▐░█▀▀▀▀█░█▀▀      ▐░▌     ▐░▌   ▀   ▐░▌     ▐░▌          ▐░▌          ▐░▌          ▐░▌   ▐░▌     ▐░█▀▀▀▀▀▀▀▀▀  ▀▀▀▀▀▀▀▀▀█░▌
#  ▐░▌          ▐░▌     ▐░▌       ▐░▌     ▐░▌       ▐░▌     ▐░▌          ▐░▌          ▐░▌           ▐░▌ ▐░▌      ▐░▌                    ▐░▌
#  ▐░▌          ▐░▌      ▐░▌  ▄▄▄▄█░█▄▄▄▄ ▐░▌       ▐░▌ ▄▄▄▄█░█▄▄▄▄      ▐░▌      ▄▄▄▄█░█▄▄▄▄        ▐░▐░▌       ▐░█▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄█░▌
#  ▐░▌          ▐░▌       ▐░▌▐░░░░░░░░░░░▌▐░▌       ▐░▌▐░░░░░░░░░░░▌     ▐░▌     ▐░░░░░░░░░░░▌        ▐░▌        ▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌
#   ▀            ▀         ▀  ▀▀▀▀▀▀▀▀▀▀▀  ▀         ▀  ▀▀▀▀▀▀▀▀▀▀▀       ▀       ▀▀▀▀▀▀▀▀▀▀▀          ▀          ▀▀▀▀▀▀▀▀▀▀▀  ▀▀▀▀▀▀▀▀▀▀▀ 

cdef class InfiniteVolume(CSGvol):
	def __init__(self, vaccum = False):
		super(InfiniteVolume, self).__init__()
		if vaccum:
			self.opaque = True


	cdef bint is_inside(self, double3& pos):
		return True

	cdef double SDF(self, double3& pos):
		return -INF

	cdef intLIST intersect(self, double3& pos, double3& dire):
		cdef intLIST result
		return result
		#return INF

cdef class Primitive(CSGvol):
	cdef Transform tr


	def __init__(self):
		super(Primitive, self).__init__()
		self.tr = Identity(self)

	def translate(self, dx, dy, dz):
		self.tr = self.tr.translate(dx, dy, dz)
		#self.mesh.translate([dx, dy, dz])

	def rotate(self, axis, angle):
		self.tr = self.tr.rotate(axis, angle)

	@property
	def matrix(self):
		return self.tr.matrix

	@property
	def inv_matrix(self):
		return self.tr.inv_matrix


cdef class Sphere(Primitive):
	cdef double r


	def __init__(self, double r):
		super(Sphere, self).__init__()
		self.r = r
		#self.mesh_gen = plt_geo.sphere(r)
		#import pyvista as pv
		#self.mesh = pv.Sphere(self.r)

	cdef bint is_inside(self, double3& _pos):
		cdef double3 pos = _pos
		self.tr.inv_pos(pos)

		return pos.x*pos.x + pos.y*pos.y + pos.z*pos.z <= self.r*self.r

	def __repr__(self):
		return f"<Sphere: r = {self.r}>"

	cdef double SDF(self, double3 _pos):
		cdef double3 pos = _pos
		self.tr.inv_pos(pos)
		return sqrt(pos.x*pos.x + pos.y*pos.y + pos.z*pos.z) - self.r

	def scale(self, s):
		self.r *= s
		return self




	cdef intLIST intersect(self, double3& _pos, double3& _dire):
		IF VERBOSE: print("INTERSECTING SPHERE")
		# direction * origin

		cdef double3 pos = _pos
		cdef double3 dire = _dire

		self.tr.inv_pos(pos)
		self.tr.inv_dire(dire)


		cdef double b = pos.x*dire.x + pos.y*dire.y + pos.z*dire.z
		# b*b - (|o|**2 - r**2)
		cdef double DELTA = b*b - pos.x*pos.x - pos.y*pos.y - pos.z*pos.z + self.r*self.r

		cdef intLIST result
		cdef Interval I

		IF VERBOSE: print(f"b = {b}, DELTA = {DELTA})")

		if DELTA <= 0:
			IF VERBOSE: print("RETURNING EMPTY")
			return result


		DELTA = sqrt(DELTA)
		b *= -1

		IF VERBOSE: print(f"proposed t2 = {b + DELTA}")
		if b + DELTA >= -1e-12:
			I.t2 = b + DELTA
		else:
			IF VERBOSE: print("RETURNING EMPTY")
			return result

		IF VERBOSE: print(f"proposed t1 = {b - DELTA}")

		if b - DELTA >= -1e-12:
			I.t1 = b - DELTA
		else:
			I.t1 = -10

		result.push_back(I)
		return result

		
#cdef class Box(Primitive):
#	cdef double x, y, z
#	def __init__(self, x, y, z):
#		self.x = x/2
#		self.y = y/2
#		self.z = z/2
#
#	cdef double SDF(self, double3 _pos):
#		cdef double3 pos = _pos
#		self.tr.inv_pos(pos)
#
#		pos.x = abs(pos.x) - self.x
#		pos.y = abs(pos.y) - self.y
#		pos.z = abs(pos.z) - self.z
#
#		return sqrt(fmax(pos.x, 0.)**2 + fmax(pos.y, 0.)**2 + fmax(pos.z, 0.)**2) + fmin(0., fmax(pos.x, fmax(pos.y, pos.z)))





































cdef class Tally(BVH):

	cdef void depositUNIFORM(self,STATE& state, double SP):
		raise RuntimeError("depositUNIFORM called from Tally (virtual)")
		import time
		print("depositUNIFORM called from Tally (virtual)")
		time.sleep(10_000)

	cdef void depositLOCAL(self, double3& pos, double E):
		raise RuntimeError("'depositLOCAL' called from its virtual in 'Volume.BVH' ")


	def __init__(self):
		pass



cdef class Z_TALLY(Tally):
	cdef vector[double] bins, counts, tmp
	cdef double3 last_pos
	cdef double L
	cdef double DZ
	cdef int Nbins


	cdef void depositRANDOM(self, STATE& state, double E, double tau):
		#print("depositing random ")
		self.bins[<int> (  (state.pos.z - state.dire.z*tau*state.genPTR.get_next_float())/self.DZ   )] += E
		state.E -= E


	def get_bins(self):
		import numpy as np
		z_axis = np.arange(0, self.Nbins)*self.DZ
		bins = np.zeros(self.Nbins)
		cdef int i
		for i in range(self.Nbins):
			#if self.counts[i] == 0: 

				#continue
			bins[i] = self.bins[i]#/self.counts[i]

		return z_axis, bins


	def __init__(self, DZ = 0.5, zmax = 1000):
		self.DZ = DZ
		self.Nbins = <int> (zmax/self.DZ) # 5 meters
		self.bins.reserve(self.Nbins)
		self.counts.reserve(self.Nbins)
		cdef int i
		for i in range(self.Nbins):
			self.bins[i] = 0
			self.counts[i] = 0


	cpdef reset(self):
		for i in range(self.Nbins):
			self.bins[i] = 0


	cdef void depositLOCAL(self, double3& pos, double E):
		if pos.z < 0:
			print("OUT OF BOUNDS, LOCAL DEPOSIT: z =", pos.z, E)
			return
		#if pos.y**2 + pos.x**2 > 1: return 
		#if self.id == state.id:
		#	self.bins[<int> (state.pos.z/self.DZ)] += E
		#	return


		#self.id = state.id
		#self.counts[<int> (pos.z/self.DZ)] += 1
		self.bins[<int> (pos.z/self.DZ)] += E


		
		#self.counts[<int> (pos.z/self.DZ)] += 1

	cdef bint move(self, STATE& state, double SP):

		IF DEBUG_MODE: input("\n ----MOVING PARTICLE----")
		IF DEBUG_MODE: input(f"CURRENT_POSITION: {state.pos}")
		IF DEBUG_MODE: input(f"Is inside current region? {(<V> state.current_region).is_inside(state.pos)}")
		IF DEBUG_MODE:
			for i in range(self.Nws):
				print(f"is_inside[{i}] = {(<V> self.ws[i]).is_inside(state.pos)}")

		IF DEBUG_MODE: input("STARTING EVENT LOOP")





		IF VERBOSE_TALLY: print(state.pos, state.L)
		#if state.pos.z < -0.01:
			#import time
			#time.sleep(10000)

		if SP != 0: self.last_pos = state.pos

		self.L = state.L

		while True:
			IF DEBUG_MODE: input(f"The safest distance is {state.pos.z}cm | Physics proposed {state.L}cm ")

			if state.pos.z > state.L:
				self.final(state)
				if SP != 0:
					self.deposit(state, SP)
				self.cache = False
				return False

			if state.pos.z < eps:
				IF DEBUG_MODE: input(f"Intersection Event")

				case = self.intEVENT(state)

				if SP != 0:
					self.deposit(state, SP)

				if case == 2: # boundary crossing
					return True

				if case == 0: # final displacement
					self.cache = False
					return False

			IF DEBUG_MODE: input(f"Virtual Event")
			self.virtual_event(state, state.pos.z)


	cdef inline int intEVENT(self, STATE& state):
		IF VERBOSE_TALLY: print("intEVENT", f"sdf = {state.pos.z}")
		IF VERBOSE_TALLY: print(f"cache = {self.cache}")
		IF VERBOSE_TALLY: print("direction:", state.dire)
		if self.cache:
			self.final(state)
			return 0

		cdef double cos = state.dire.z
		cdef double t
		if abs(cos) > 0.0001:
			t = -state.pos.z / cos
			if state.L < t:
				self.final(state)
				return 0

			self.virtual_event(state, t)
			self.boundary_crossing(state)
			self.cache = True
			state.current_region = <void*> self.outer
			return 2
		self.final(state)
		return 0


	cdef void boundary_crossing(self, STATE& state):
		state.current_region = <void*> self.outer


	cdef void deposit(self, STATE& state, double SP):


		cdef int i0 = <int> (self.last_pos.z/self.DZ)
		cdef int i1 = <int> (state.pos.z/self.DZ)

		if i0 < 0 or i1 < 0:
			print("OUT OF BOUNDS:", self.last_pos, state.pos)
			return

		IF VERBOSE_TALLY: print(i0, i1)

		if i0 == i1:
			self.bins[i0] += (self.L-state.L)*SP
			self.counts[i0] += 1
			return

		cdef double dz
		cdef double cos

		dz = self.last_pos.z - i0*self.DZ
		cos = abs(state.dire.z)

		self.bins[i0] += SP*(self.DZ - dz)/cos
		self.counts[i0] += 1

		self.bins[i1] += SP*(state.pos.z - i1*self.DZ)/cos
		self.counts[i1] += 1

		#else: # change indexes around
		#	dz = self.last_pos.z - i1
		#	cos = abs(state.dire.z)
#
		#	self.bins[i1] += (self.DZ - dz)/cos*SP
#
		#	self.bins[i0] += (state.pos.z - i0)/cos*SP


		if abs(i1 - i0) == 1:
			return

		cdef double dE = self.DZ*SP/cos

		if i0 < i1:
			for i in range(i0+1, i1):
				self.bins[i] += dE
				self.counts[i] += 1
		else:
			for i in range(i1+1, i0):
				self.bins[i] += dE
				self.counts[i] += 1



	cdef double SDF(self, double3 pos):
		return -pos.z


	cdef void localSDF(self, STATE& state):
		"""
		JOB OF 'localSDF':
			(1) Set attribute .sdf equal to a *safe distance*

		note: 'safe distance' does not imply the result of a signed
			   distance function (despite the name I gave it).
			   The rest of the code will work as long as this distance can
			   be travelled without hitting some other surface.
		"""
		if self.cache:
			self.sdf = INF
			return

		self.sdf = -state.pos.z


	cdef double main_intersect(self, STATE& state):
		"""
		JOB OF 'main_intersect':
			(1) if previously marked as cache:
					PROVIDE DISTANCE TO THE NEXT CACHED INTERSECTION

			(2) if not marked as cache:
				(1.1) Intersect ray with self
				(1.2) collect intersections into an 'intLIST'
				(1.3) create and store an 'intIerator'
				(1.4) return distance to the next intersection
		"""
		IF VERBOSE_TALLY: print("main_intersect", f"cache = {self.cache}")

		if self.cache:
			# cache implies previous boundary crossing
			# after a BC, only way forward is to infinity
			return INF

		self.cache = True

		cdef double cos = state.dire.z
		cdef intLIST crosses
		cdef Interval I

		if abs(cos) > 0:
			I.t1 = -10
			I.t2 = -state.pos.z / cos
			IF VERBOSE_TALLY: print(I.t2)
			crosses.push_back(I)

			self.cross = intIterator(crosses)
			IF VERBOSE_TALLY: print("next_crosS:", self.cross.current())
			return self.cross.current()

		return INF

	cdef inline void final(self, STATE& state):



		state.last_displacement = state.L

		state.pos.x += state.dire.x*state.L
		state.pos.y += state.dire.y*state.L
		state.pos.z += state.dire.z*state.L

		state.L = 0

	cdef bint is_inside(self, double3& pos):
		return pos.z >= 0


	cdef inline void virtual_event(self, STATE& state, double dr):

		state.pos.x += state.dire.x*dr
		state.pos.y += state.dire.y*dr
		state.pos.z += state.dire.z*dr

		state.last_displacement = dr

		state.L -= dr
