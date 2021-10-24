__doc__ = """Unit-testing the `.tools` package.
"""

__author__ = "Rui Campos"



import _cmd
import sys
del sys.argv[1]


import numpy as np
import unittest as ut


from MontyCarlo.geometry.CSG import Sphere

class test_Sphere(ut.TestCase):
    
    def test_updates(self):
        sphere = Sphere(1)
        
        for x, y, z in zip(
            range(10), range(10), range(10)
        ):
            print(
                sphere.SDF(.3*x, .3*y, .3*z)
            )