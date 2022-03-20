










cdef class Sphere(Primitive):
	cdef double r
	cdef double r2


	def __init__(self, *args, **kwargs):

		self.r = kwargs['radius']
		self.r2 = self.r*self.r

		super(Sphere, self).__init__(*args, **kwargs)
		


	def __repr__(self):
		return f"<Sphere: radius={self.r}cm>"

	cdef bint is_inside(self, double3& _pos):
		cdef double3 pos = _pos
		self.apply_map(pos, self.inverse_transform)
		return pos.x*pos.x + pos.y*pos.y + pos.z*pos.z <= self.r2


	cdef double SDF(self, double3& pos):

		cdef double3 _pos = pos


		# applying the inverse transformations 
		self.apply_map(
			_pos, self.inverse_transform
		)

		# .... not sure what this was? :
        # self.to_primitive_space(pos)

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



	cdef intLIST intersect(self, double3& _pos, double3& _dire):
		IF VERBOSE: print("SPHERE::INTERSECTING")


		cdef double3 pos = _pos
		cdef double3 dire = _dire

		self.tr.inv_pos(pos)
		self.tr.inv_dire(dire)


		cdef double b = pos.x*dire.x + pos.y*dire.y + pos.z*dire.z
		# b*b - (|o|**2 - r**2)
		cdef double DELTA = b*b - pos.x*pos.x - pos.y*pos.y - pos.z*pos.z + self.r2

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
