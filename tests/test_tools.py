__doc__ = """Unit-testing the `.tools` package.
"""

__author__ = "Rui Campos"



import _cmd
import sys
del sys.argv[1]


import numpy as np
import unittest as ut

# Importing 
from MontyCarlo.tools.main import python_hooks
remove_duplicates = python_hooks.remove_duplicates


class input_val:
    """A namespace indicating input values.
    """
    test_cases = [[ 1,  1,  2, 3, 4, 5, 6],
                  [-1, -1, -4, 2, 5, 7, 7],
                  [1., 1., 1., 2., 2., 2., 3., 4., 5., 5., 5., 5., 6., 6., 6, 6, 6, 7, 1000, 1000, 5000]]


class ground_truth:
    """A namespace indicating groundtruth.
    """
    test_cases = [[ 1,  2, 3, 4, 5, 6],
                  [-1, -4, 2, 5, 7],
                  [1., 2., 3., 4.,  5., 6, 7, 1000, 5000]]
    

class output_val:
    """A namespace indicating calculated values.
    """
    pass


functions = [lambda x: x**2,
             lambda x: x + 1,
             lambda x: x/10,
             np.cos]

class test_remove_duplicates(ut.TestCase):
    """Unit testing the `.tools.main.remove_duplicates` function.
    """

    def test_all(self):
        """Test all cases.
        """

        for i, x_axis in enumerate(input_val.test_cases):
            input_val.X    = np.array(x_axis, dtype = float)
            ground_truth.X = np.array(ground_truth.test_cases[i], dtype = float)

            for func in functions:
                input_val.Y    = func(input_val.X)
                ground_truth.Y = func(ground_truth.X)

                output_val.X, output_val.Y = remove_duplicates(input_val.X, input_val.Y)

                self.assertEqual(
                                 list(output_val.X), 
                                 list(ground_truth.X), 
                                 f"Should be {ground_truth.X}"
                                )

                self.assertEqual(
                                 list(output_val.Y), 
                                 list(ground_truth.Y),
                                 f"Should be {ground_truth.Y}"
                                )


if __name__ == '__main__':
    ut.main()