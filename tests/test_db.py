
__doc__ = """
"""

__author__ = "Rui Campos"

import _cmd


import numpy as np
import unittest as ut


from pathlib import Path
__directory__ = Path(repr(__file__)[1:])
__directory__ = __directory__.parent


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
        work_dir = __directory__.parent
        work_dir = _dir/'MontyCarlo'/'materials'/'electron'/'elastic'

        self.assertTrue(work_dir.exists(), msg = "`electron/elastic` directory does not exist.")
        self.assertTrue(work_dir.is_dir(), msg = "`electron/elastic` is not a directory.")


        for i in range(100):
            element_dir = work_dir/str(i)

            self.assertTrue(element_dir.exists(), msg = f"`electron/elastic/{i}` directory does not exist.")
            self.assertTrue(element_dir.is_dir(), msg = f"`electron/elastic/{i}` is not a directory.")



    def positron_elastic(self):
        pass



if __name__ == '__main__':
	ut.main()

