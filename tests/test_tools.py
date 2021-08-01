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

class ground_truth:
    """A namespace indicating groundtruth.
    """
    pass

class output_val:
    """A namespace indicating calculated value.
    """


def check_for_duplicates(arr):
    """Checks if `arr` contains duplicate values.
    """

    print(f"""
    RESULT: {arr}
    """)


    found = []
    for x in arr:
        if x not in found:
            found.append(x)
            continue

        raise RuntimeError("Found duplicates.")


print("""











""")

input_val.test_cases = [[ 1,  1,  2, 3, 4, 5, 6],
                        [-1, -1, -4, 2, 5, 7, 7],
                        [1., 1., 1., 2., 2., 2., 3., 4., 5., 5., 5., 5., 6., 6., 6, 6, 6, 7, 1000, 1000, 5000]]

ground_truth.test_cases = [[ 1,  2, 3, 4, 5, 6],
                           [-1, -4, 2, 5, 7],
                           [1., 2., 3., 4.,  5., 6, 7, 1000, 5000]]

functions = [lambda x: x**2,
             lambda x: x + 1,
             lambda x: x/10]

class main(ut.TestCase):

    def test_all(self):

        for x_axis in input_val.test_cases:
            input_val.X = np.array(x_axis, dtype = float)

            for func in functions:
                input_val.Y = f(func)

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




    def test1(self):

        input_val.X = np.array([ 1,  1,  2, 3, 4, 5, 6], dtype = float)
        input_val.Y = np.array([-1, -1, -4, 2, 5, 7, 7], dtype = float)

        ground_truth.X = np.array([  1,  2, 3, 4, 5, 6], dtype = float)
        ground_truth.Y = np.array([ -1, -4, 2, 5, 7, 7], dtype = float)

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

    def test2(self):

        input_val.X = np.array([ 1., 1., 1., 2., 2., 2., 3., 4., 5., 5., 5., 5., 6., 6., 6, 6, 6, 7, 1000, 1000, 5000], dtype = float)
        input_val.Y = input_val.X**2

        ground_truth.X = np.array([1, 2, 3, 4, 5, 6, 7, 1000, 5000], dtype = float)
        ground_truth.Y = ground_truth.X**2

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
    print("""







STARTING UNIT TEST
""")
    ut.main()
    print("""









""")