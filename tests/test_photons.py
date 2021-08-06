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


    # A basic set-up for holding one particle -----------------------
    from MontyCarlo.geometry.CSG import Sphere
    from MontyCarlo.geometry.CSG import InfiniteVolume
    from MontyCarlo.materials.materials import Mat
    
    from MontyCarlo._init import eax

    print("Creating photon...")
    photon = Photon()

    print("Creating water...")
    water = Mat({1:2, 8:1}, 1)

    print("Creating geometry...")
    with InfiniteVolume() as OUTER:
        OUTER.fill(water)
        OUTER.configure("no_name", render = False)
        with Sphere(1) as sphere:
            sphere in OUTER
            sphere.fill(water)
            sphere.configure("no_name", render = False)

    print("Setting current region...")
    photon.current_region = sphere
 
    print("UPDATING .........................")
    photon.update_references()
    photon.update_imfp()

    print("WORKED")
    # ----------------------------------------------------------------


    def test_updates(self):
        """Checks for segmentation errors when calling update methods.
        """
        print("\n\nTESTING UPDATES")
        cls = test_Photon
        cls.photon.update_references()
        cls.photon.update_imfp()

    def test_find_index(self):
        print("\n\n TESTING `find_index`")
        cls = test_Photon
        import numpy.random as npr

        photon = cls.photon
        eax    = cls.eax


        N = len(eax) - 1
        points = [1e3, 1.1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1.9e8, 1e9, ]

        for E0, Ef in zip(points[:-1], points[1:]):
            print(f"Testing `find_index` in range [{E0}, {Ef}]")

            for i in range(50_000):
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
        print("\n\nTESTING INCOHERENT")
        from collections import deque
        cls = test_Photon
        photon = cls.photon
        print("Seeding photon:")
        photon.set_seed(1234)
        points = [1e3, 1.1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1.9e8, 1e9 ]
        
        for E in points:
            print(f"Running incoherent for energy {E}eV")
            photon.k = E/0.5110e6
            photon.E = E
            photon.secondary = deque()
            photon._incoherent()
            

    def test_coherent(self):
        print("\n\nTESTING COHERENT")
        cls = test_Photon
        photon = cls.photon
        print("Seeding photon:")
        photon.set_seed(1234)
        points = [1e3, 1.1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1.9e8, 1e9 ]
        
        for E in points:
            print(f"Running coherent for energy {E}eV")
            k = E/0.5110e6
            photon.k = E/0.5110e6
            photon._coherent()
            self.assertEqual(photon.k, k, msg = "Coherent is not conserving energy...")
            
            


if __name__ == '__main__':
    ut.main()
