

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



	cdef double SDF(self, double3 pos):
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

	cdef double SDF(self, double3 pos):
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





















