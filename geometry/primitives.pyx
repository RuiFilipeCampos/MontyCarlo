# cython: profile=True

#external imports
import numba as nb

#internal imports
from ..tools.vectors cimport *


cdef bint OR(bint a, bint b):
	return a or b

cdef bint AND(bint a, bint b):
	return a and b

cdef bint NOT(bint a, bint b):
	return not a


cdef class Volume:
	"""
	Logic that all volume types inherit.
	primitives (Sphere Cube, etc) inherit this logic
	composition volumes (Volume) inherit this logic
	"""

	cpdef Volume fill(self, object material):
		self.material = material
		self.N = 0

	def __lshift__(self, Volume other):
		if not self.material:
			raise RuntimeError("You are trying to insert a volume inside a volume that does not have a material.")
		if not other.material:
			raise RuntimeError("That volume you are trying to insert is not filled!")
		self.inner += [other]
		other.outer = self

		return self

	def __str__(self):
		name = self.__repr__()
		outer_name = self.outer.__repr__()
		return name + "\n Containing: " + str(self.inner) + "\n Contained by: " + outer_name

	cdef public bint _eval(self, double x, double y, double z):
		"""To be overwriten"""
		return False







		

	
	
	#def __invert__(self): #******  THERE IS TROUBLE HERE
	#	"""Negation of a region. Example: inverted_region = ~region"""
	#	@nb.njit(nb.b1(nb.float64, nb.float64))
	#	def new_condition(double x, double y, double z):
	#		return not self.condition(x, y, z)
	#
	#	vaccum = UnboundedVolume(%NOT, self)
	#	self.outer = vaccum
	#	return vaccum
		
	#def __sub__(self, other):
	#	"""
	#	Subtraction of two regions:
	#		new_region = region2 - region1
	#	region1 has been striped of all the points that it had in common with region1
	#	"""
	#	@nb.njit(nb.b1(nb.float64, nb.float64))
	#	def new_condition(double x, double y, double z):
	#		return (not other.condition(x, y, z)) and self.condition(x, y, z)
	#	#eturn (~other) & self
	#	return VolumeNode(self, other, new_condition)



	def __contains__(self, Vector pos):
		"""Checks if position is in this region. Example:
		> Vector(x, y, z) in region
		returns True or False.
		"""

		return self._eval(pos.x, pos.y, pos.z)

	def __iter__(self):
		yield from self.inner




	cdef public double getOuterCrossing(self, Vector pos, Vector ez, double L):
		cdef list crossings = self.intercept(pos, ez, L)
		if crossings:
			return min(crossings)
		return 0.



	cdef public tuple getCrossings(self, Vector pos, Vector ez, double L):

		#first, test intersection with bounding volume
		cdef double dL = self.getOuterCrossing(pos, ez, L)
		cdef Volume proposed_vol

		if dL == 0.: proposed_vol = self
		else:        proposed_vol = self.outer

		#then test intersection with inner volumes
		cdef double new_dL
		cdef Volume region
		for region in self:
			new_dL = region.getOuterCrossing(pos, ez, L)
			if 0. < new_dL < dL:
				proposed_vol = region
				dL = new_dL

		#no intersections, in or out
		if dL == 0.:
			return pos + ez*L, self

		#if inner, just return it
		if proposed_vol is not self.outer:
			return pos + ez*dL, proposed_vol

		#if outer volume, search inside for a region change (send out a probe)
		cdef Vector probe  = pos + ez*(dL + 0.0001)

		for region in self.outer:
			if probe in region:
				return pos + ez*dL, region
		else:   return pos + ez*dL, self.outer


	#LOGIC	
	def __and__(self, Volume other):
		"""Intersection of two regions. Example: new_region = region1 & region2"""
		cdef Volume vol = VolumeNode(self, other)
		vol.condition = AND
		return vol

	
	def __or__(self, Volume other):
		"""Union of two regions. Example: new_region = region1 | region2"""
		cdef Volume vol = VolumeNode(self, other)
		vol.condition = OR
		return vol

	def __invert__(self):
		self.outer = UnboundedVolume(self)
		return self.outer




cdef class VolumeNode(Volume):
	"""Generic volume. Consists of the composition of two volumes/primitives."""
	def __init__(self, Volume volume1, Volume volume2):
		self.imp = 1
		self.vol1 = volume1
		self.vol2 = volume2
		self.inner = []

	cpdef list intercept(self, Vector pos, Vector ez, double L):
		return self.vol1.intercept(pos, ez, L) + self.vol2.intercept(pos, ez, L)

	cdef bint _eval(self, double x, double y, double z):
		return self.condition(self.vol1._eval(x, y, z), 
				              self.vol2._eval(x, y, z))

	def __repr__(self):
		return "<Compound Volume>"


































cdef class UnboundedVolume(Volume):
	def __init__(self, Volume ROI):
		self.imp = 0
		self.condition = NOT
		self.inner = [ROI]

	cdef bint _eval(self, double x, double y, double z):
		return not self.inner[0]._eval(x, y, z)

	cpdef list intercept(self, Vector pos, Vector ez, double L):
		return []

	def __repr__(self):
		return "<Unbounded Volume>"






































cdef class Sphere(Volume):
	def __init__(self, double radius, double x0, double y0, double z0):
		self.inner = []
		self.imp = 1
		self.radius = radius
		self.x0, self.y0, self.z0 = x0, y0, z0
		self.center = Vector(x0, y0, z0)

	



	#def move(Vector displacement):
	#	self.center = self.center + displacement
	#	self.x0, self.y0, self.z0 = self.center.x, self.center.y, self.center.z

	cdef bint _eval(self, double x, double y, double z):
		return (self.x0 - x)**2 + (self.y0 - y)**2 + (self.z0 - z)**2 <= self.radius**2

	#cpdef list intercept(self, Vector P, Vector Q):
	#	return self._intercept(P, Q)

	cpdef list intercept(self, Vector pos, Vector ez, double L):
		"""
		Parametric equation for ray:
		position_in_ray = ray.P + ray.Q * t

		Equation of sphere:
		(P - self.center) | (P - self.center) = self.radius**2

		Substituting position_in_ray

		(ray.P + ray.Q * t - self.center)**2 = self.radius**2

		
		ray.P**2 + self.center**2 - 2* self.center | ray.P - self.radius**2
		+ ray.Q **2 * t**2 
		+ 2 * (ray.Q | ray.P -  self.center | ray.Q)  t 
		= 0 

		okaaaaay....

		We get, A*t**2 + B*t + C = 0 (scalars only)

		t = [-b +/- sqrt(discriminant)] / (2*a)
		discriminant = b**2 - 4*a*c
		discriminant < 0 --> forget about iiitttt
		discriminent = 0 --> tangent (ignore??)
		discriminent > 0 --> intercection, get both points, both are needed
		"""
		#cdef Vector U =

		cdef double b = 2 * ( ez | (pos - self.center) )
		cdef double c = (pos >> self.center)**2 - self.radius**2

		cdef double discriminant = b**2 - 4*c

		if discriminant <= 0: return []

		cdef double t1 =  (-b  + discriminant**.5 ) * .5
		cdef double t2 =  (- b - discriminant**.5 ) * .5

		cdef bint cond1 = 0. < t1 <= L
		cdef bint cond2 = 0. < t2 <= L

		if cond1 and cond2:
			return [min(t1, t2)]
		if cond1 and not cond2:
			return [t1]
		if not cond1 and cond2:
			return [t2]
		return []




	def __repr__(self):
		return f"<Primitive Sphere: r = {self.radius}, center = ({self.x0},{self.y0},{self.z0})>"

