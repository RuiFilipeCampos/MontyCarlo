









def Translate(*args, delta_x=0, delta_y=0, delta_z=0):
    for volume in args:
        volume.translate(delta_x, delta_y, delta_Z)
        Translate(*volume.ws, 
            delta_x=delta_x, 
            delta_y=delta_y, 
            delta_z=delta_z,
        )
    return args



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