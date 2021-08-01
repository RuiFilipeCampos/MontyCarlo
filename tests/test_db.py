
__doc__ = """
"""

__author__ = "Rui Campos"

import _cmd
import sys
del sys.argv[1]

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
    pass





top_level_files = ["HEeax.npy", "LEeax.npy", "muGRID.npy"]
element_level_files = ["DCS.npy", "HEtransportTCS.npy", "LEtransportTCS.npy"]



class AreDBPresent(ut.TestCase):
    """Checks if databases are where they should be.
    """



    def test_electron_elastic(self):
        """Checks if the database for elastic scattering of electrons is present.
        """

        work_dir = __directory__.parent
        work_dir = work_dir/'MontyCarlo'/'materials'/'electron'/'elastic'

        self.assertTrue(work_dir.exists(), msg = f"`{work_dir}` directory does not exist.")
        self.assertTrue(work_dir.is_dir(), msg = f"`{work_dir}` is not a directory.")

        for file_name in top_level_files:
            file = work_dir/'file_name'
            self.assertTrue(file.exists(), msg = "Missing file: " + str(file))

        for i in range(100):
            element_dir = work_dir/str(i)

            self.assertTrue(element_dir.exists(), msg = f"`electron/elastic/{i}` directory does not exist.")
            self.assertTrue(element_dir.is_dir(), msg = f"`electron/elastic/{i}` is not a directory.")

            for file_name in element_level_files:
                file = element_dir/file_name
                self.assertTrue(file.exists(), msg = "Missing file: " + str(file))


    def positron_elastic(self):
        pass

    def EADL(self):
        pass

    def EEDL(self):
        pass

    def EPDL(self):
        pass

if __name__ == '__main__':
	ut.main()

