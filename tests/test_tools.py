__doc__ = """
"""

__author__ = "Rui Campos"

import _cmd


import numpy as np 
from MontyCarlo.tools.main import remove_duplicates



def check_for_duplicates(arr):
	"""Checks if `arr` contains duplicate values.
	"""

	print(f"""
	RESULT:	{arr}
	""")


	found = []
	for x in arr:
		if x not in found:
			found.append(x)
			continue

		raise RuntimeError("Found duplicates.")



X = np.array([1, 1, 2, 3, 4, 5, 6])
Y = np.array([-1, -1, -4, 2, 5, 7])

print(f"""
	Array being tested:
	{X}
	""")

x, y = remove_duplicates(X, Y)

check_for_duplicates(x)