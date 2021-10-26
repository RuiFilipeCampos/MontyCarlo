










cdef class Sphere(Primitive):
	cdef double r


	def __init__(self, *args, **kwargs):
		super(Sphere, self).__init__(*args, **kwargs)
		self.r = kwargs['radius']


	cdef bint is_inside(self, double3& _pos):
		cdef double3 pos = _pos
		self.tr.inv_pos(pos)
		return pos.x*pos.x + pos.y*pos.y + pos.z*pos.z <= self.r*self.r

	def __repr__(self):
		return f"<Sphere: radius={self.r}cm>"

	cdef double SDF(self, double3& _pos):

		cdef double3 pos = _pos

		self.tr.inv_pos(pos)
		return sqrt(
			pos.x*pos.x +
		    pos.y*pos.y +
			pos.z*pos.z
			) - self.r


	def public__SDF(self, double x, double y, double z):

		print(f"Arguments: x={x}, y={y}, z={z}")
		cdef double3 pos;
		pos.x = x;
		pos.y = y;
		pos.z = z;
		print(pos)
		return self.SDF(pos)

	def scale(self, s):
		self.r *= s
		return self

	cdef intLIST intersect(self, double3& _pos, double3& _dire):
		IF VERBOSE: print("SPHERE::INTERSECTING")


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
