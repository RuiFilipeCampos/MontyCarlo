

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
                    tmp_direct_transform[3] += transformation[1]
                    tmp_direct_transform[7] += transformation[2]
                    tmp_direct_transform[11] += transformation[3]

                if transformation[0] == "rotate":
                    rotated = True
                    pass

        if translated and rotated:
            self.apply_transform = &apply_general_transform
        
        elif translated:
            self.apply_transform = &apply_translation
        
        elif rotated:
            self.apply_transform = &apply_rotation

		super(Primitive, self).__init__(*args, **kwargs)
