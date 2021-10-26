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

import numpy as np

ctypedef (*map_t)(double3& pos, double[:] transformation)


cdef class Primitive(CSGvol):

    cdef map_t apply_transform
	cdef double[16] direct_transform
    cdef double[16] inverse_transform


	def __init__(self, *args, **kwargs):

        tmp_direct_transform = np.array(
            [
                [1, 0, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1],
            ], dtype = np.float64,
        )


        translated = False
        rotated = False

        if kwargs['transforms']:
            for transformation in kwargs['transforms']:

                if transformation[0] == "translate":
                    translated = True
                    tmp_direct_transform[3]  += transformation[1]
                    tmp_direct_transform[7]  += transformation[2]
                    tmp_direct_transform[11] += transformation[3]

                elif transformation[0] == "rotate":
                    rotated = True
                    axis = np.array([
                        transformation[1][0], 
                        transformation[1][1], 
                        transformation[1][2], 
                    ], 
                    dtype = np.float64 )

                    axis = axis/np.sqrt(np.sum(axis**2))

                    ux = axis[0]
                    uy = axis[1]
                    uz = axis[2]

                    angle = transformation[2]

                    rotation_matrix = np.zeros((4,4))

                    _cos = cos(angle)
                    oneMcos = 1 - _cos
                    _sin = sqrt(1 - _cos*_cos)

                    rotation_matrix[0, 0] = _cos + ux*ux*oneMcos
                    rotation_matrix[0, 1] = ux*uy*oneMcos - uz*_sin
                    rotation_matrix[0, 2] = ux*uz*oneMcos + uy*_sin

                    rotation_matrix[1, 0] = uy*ux*oneMcos + uz*_sin
                    rotation_matrix[1, 1] = _cos + uy*uy*oneMcos
                    rotation_matrix[1, 2] = uy*uz*oneMcos - ux*_sin

                    rotation_matrix[2, 0] = uz*ux*oneMcos-uy*_sin
                    rotation_matrix[2, 1] = uz*uy*oneMcos + ux*_sin
                    rotation_matrix[2, 2] = _cos + uz*uz*oneMcos

                    rotation_matrix[3, 3] = 1

                    tmp_direct_transform = rotation_matrix @ tmp_direct_transform

        if translated and rotated:
            self.apply_transform = &apply_general_transform
        
        elif translated:
            self.apply_transform = &apply_translation
        
        elif rotated:
            self.apply_transform = &apply_rotation

		super(Primitive, self).__init__(*args, **kwargs)
