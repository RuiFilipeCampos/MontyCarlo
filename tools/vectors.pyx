# cython: profile=True

from pyquaternion import Quaternion


cdef class Vector:
	"""
	This class (Vectors) was created to substitute Numpy arrays.
	It is specialized to store 3 numbers of type double that behave
	like a vector.

	>> Available operations:
	Vector addition: P + Q (commutative)
	Vector subtraction: P - Q and Q - P
	Scalar multiplication: P * s
	>>>>>> Note that the reverse operation s * P is not implemented. It would
	>>>>>> since it would require one extra comparison (and calling native Python)
	Dot Product: P|Q and Q|P
	Distance from vector P to vector Q: P>>Q and Q<<P (they return the same result.)

	Two extra methods are implemented to synergize with the rest of the code:

	self.rotateCos(axis, cos)
	self.rotateAngle(axis, theta)

	These actually call a Quartenion class outside of Cython. This behaviour will
	change in the future.
	"""
	#cdef public double x, y, z
	def __init__(self, double x, double y, double z):
		self.x, self.y, self.z, self.real = x, y, z, 0.





	def __add__(self, Vector other):
		"""Elementwise addition: P + Q"""
		return Vector(self.x + other.x, self.y + other.y, self.z + other.z)

	def __radd__(self, Vector other):
		"""Elementwise addition: Q + P """
		return self.__add__(other)

	def __sub__(self, Vector other):
		"""Elementwise subtraction: P - Q """
		return Vector(self.x - other.x, self.y - other.y, self.z - other.z)

	def __rsub__(self, Vector other):
		"""Elementwise subtraction: Q - P """
		return Vector(other.x - self.x , other.y - self.y ,  other.z - self.z)

	def __mul__(self, double other):
		"""
		Scalar multiplication: P * s

		Note: s * P is not defined. Implementing it is possible but requires adding
		extra comparisons. 
		"""
		return Vector(self.x*other, self.y*other, self.z*other)

	def __or__(self, Vector other):
		"""Dot product: P | Q """
		return self.x*other.x + self.y*other.y + self.z*other.z

	def __ror__(self, Vector other):
		"""Dot product: Q | P """
		return self.__or__(other)

	def __rshift__(self, Vector other):
		"""Distance between two vectors: dist = P >> Q."""
		return  ( (self.x - other.x)**2 + (self.y - other.y)**2 + (self.z - other.z)**2)**.5 


	def __lshift__(self, Vector other):
		"""Distance between two vectors: dist = P << Q."""
		return  ( (self.x - other.x)**2 + (self.y - other.y)**2 + (self.z - other.z)**2)**.5 

	def __repr__(self):
		return f"Vector({self.x}, {self.y}, {self.z})"

	def __truediv__(self, double other):
		return Vector(self.x/other, self.y/other, self.z/other)



	cdef public double norm(self):
		return (self.x**2 + self.y**2 + self.z**2)**.5

	cdef public Vector normalize(self):
		return self/self.norm()




	def __xor__(self, Vector q):
		cdef Vector result = Vector(self.real*q.x  + self.x*q.real + self.y*q.z    - self.z*q.y,
									self.real*q.y  - self.x*q.z    + self.y*q.real + self.z*q.x,
									self.real*q.z  + self.x*q.y    - self.y*q.x    + self.z*q.real)
		
		result.real = self.real*q.real - self.x*q.x - self.y*q.y - self.z*q.z
		return result


	cdef Vector _mult(self, Vector q):

		cdef Vector result = Vector(self.real*q.x  + self.x*q.real + self.y*q.z    - self.z*q.y,
									self.real*q.y  - self.x*q.z    + self.y*q.real + self.z*q.x,
									self.real*q.z  + self.x*q.y    - self.y*q.x    + self.z*q.real)
		
		result.real = (self.real*q.real) - (self | q)
		return result

	cdef Vector _conj(self):
		cdef Vector result = Vector(-self.x, -self.y, -self.z)
		result.real = self.real
		return result



	def rotateAngle(self, Vector axis, double theta):
		""" Rotates vec by an angle theta along the provided axis """
		cdef tuple _axis, _vec
		_axis = (axis.x, axis.y, axis.z)
		_vec  = (self.x, self.y, self.z)

		cdef object qq
		qq = Quaternion(axis = _axis, angle = theta).rotate(_vec)

		return Vector(qq[0], qq[1], qq[2]) 


	#rotate = lambda vec, axis, theta: Quaternion(axis=axis, angle=theta).rotate(vec)


	cdef public Vector rotateCos(self, Vector axis, double cos):
		"""Rotate vector along axis with angle defined by c = cos(theta)"""
		self.real = 0.
		
		cdef double cos2, sin2
		cos2 = (.5*(1 + cos))**.5
		sin2 = (.5*(1 - cos))**.5

		cdef Vector axis_angle = axis * sin2
		axis_angle.real = cos2

		return (axis_angle ^ self ^ axis_angle._conj())

		#q = Quaternion(cos2, axis_angle.x, axis_angle.y, axis_angle.z)
		#q = (q * vector * q.conjugate).imaginary

		#return Vector(q[0], q[1], q[2]



