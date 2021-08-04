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
    from MontyCarlo._init import eax

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

    def test_find_index(self):
        import numpy.random as npr

        cls = test_Photon
        photon = cls.photon
        eax = cls.eax
        N = len(eax) - 1
        points = [1e3, 1.1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1.9e8, 1e9, ]

        for E0, Ef in zip(points[:-1], points[1:]):

            for i in range(10_000):
                E = E0 + npr.rand()*(Ef - E0)

                photon.E = E
                i = photon.find_index()

                error_msg = f"""
                INVALID INDEX
                --------------
                photon.find_index() failed for E = {E} eV. 

                It found index i = {i}. Which is out of range for the array `eax`.
                """


                self.assertTrue(0 <= i <= N, msg = error_msg)

                error_msg = f"""
                FOUND WRONG INDEX
                -----------------
                photon.find_index() failed for E = {E}eV. 

                It found index i = {i}. Corresponding to the following interval:

                {eax[i]} <= {E} < {eax[i+1]}

                Note: eax[i] <= E < eax[i+1]
                """

                self.assertTrue( eax[i] <= E < eax[i+1], msg = error_msg)
                #print(f"SUCCESS: {eax[i]} <= {E} < {eax[i+1]}")





    def test_compton(self):
        pass

    def test_coherent(self):
        pass


if __name__ == '__main__':
    ut.main()