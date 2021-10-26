

import numpy as np


cdef class Primitive(CSGvol):

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

        tmp_inverse_transform = np.linalg.inv(tmp_direct_transform)
        
        cdef int i, y, x

        i = 0
        for y in range(4):
            for x in range(4):
                self.direct_transform[i]  = tmp_direct_transform[x, y]
                self.inverse_transform[i] = tmp_inverse_transform[x, y]
                i += 1

        if translated and rotated:
            self.apply_transform = &apply_general_transform
        
        elif translated:
            self.apply_transform = &apply_translation
        
        elif rotated:
            self.apply_transform = &apply_rotation

		super(Primitive, self).__init__(*args, **kwargs)



    def translate(self, double x, double y, double z):  
