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
    """Unit testing the `.tools.main.remove_duplicates` function.
    """

    def test_all(self):
        """Test all cases.
        """
        py_state = PySTATE()
        photon = Photon(py_state)


if __name__ == '__main__':
    ut.main()