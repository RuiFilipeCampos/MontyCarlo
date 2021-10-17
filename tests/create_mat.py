__doc__ = """
"""

import _cmd
import sys
del sys.argv[1]

import MontyCarlo as myco
myco.Mat({1:2, 8:1}, 1)
myco.Mat({4:2, 9:1}, 2)