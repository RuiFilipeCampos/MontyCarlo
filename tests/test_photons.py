__doc__ = """Unit-testing the `.tools` package.
"""

__author__ = "Rui Campos"



import _cmd
import sys
del sys.argv[1]


import numpy as np
import unittest as ut

# Importing 
from MontyCarlo.types import PySTATE
from MontyCarlo.particles.photons import python_hooks
Photon = python_hooks.Photon


class input_val:
    """A namespace indicating input values.
    """
    pass

class ground_truth:
    """A namespace indicating groundtruth.
    """
    pass
    

class output_val:
    """A namespace indicating calculated values.
    """
    pass


class test_Photon(ut.TestCase):
    """Unit testing photons.
    """
    from MontyCarlo.geometry.CSG import Sphere
    from MontyCarlo.geometry.CSG import InfiniteVolume
    from MontyCarlo.materials.materials import Mat

    py_state = PySTATE()
    photon = Photon(py_state)
    water = Mat({1:2, 8:1}, 1)

    with InfiniteVolume() as OUTER:
        OUTER.fill(water)
        OUTER.configure("no_name", render = False)
        with Sphere(1) as sphere:
            sphere in OUTER
            sphere.fill(water)
            sphere.configure("no_name", render = False)

    photon.current_region = sphere


    def test_updates(self):
        """Checks for segmentation errors when calling update methods.
        """
        cls = test_Photon
        cls.photon(method = "update_references")
        cls.photon(method = "update_imfp")

    def test_compton(self):
        pass

    def test_coherent(self):
        pass


if __name__ == '__main__':
    ut.main()