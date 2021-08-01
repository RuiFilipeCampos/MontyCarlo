
__doc__ = """
"""

__author__ = "Rui Campos"

import _cmd


import numpy as np
import unittest as ut

from MontyCarlo.tools.main import remove_duplicates


class input_val:
	"""A namespace indicating input values.
	"""

class ground_truth:
	"""A namespace indicating groundtruth.
	"""
	pass

class output_val:
	"""A namespace indicating calculated value.
	"""


class AreDBPresent(ut.TestCase):
    """Checks if databases are where they should be.
    """

    def EADL(self):
        pass

    def EEDL(self):
        pass

    def EPDL(self):
        pass

    def electron_elastic(self):
        pass

    def positron_elastic(self):
        pass



if __name__ == '__main__':
	ut.main()

