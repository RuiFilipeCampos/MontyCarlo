__doc__ = """Handling command line arguments for all unit testing modules. 

This is necessary since there are several possible scenarios when unit testing Monty Carlo:

    (1) MyCo is installed through pip
    (2) MyCo has been built inplace at `/MontyCarlo`
    (3) ...
"""

__author__ = "Rui Campos"


import sys
from pathlib import Path

if len(sys.argv) == 1:
	raise RuntimeError("Missing arguments.")

arg1 = sys.argv[1]

if arg1 == "--pip_installed" or arg1 == "-p":
	pass

elif arg1 == "--built_inplace" or arg1 == "-b":
	__path__ = Path(__file__)
	__folder__ = __path__.parent
	sys.path.append(str(__folder__.parent))
	print(f"The path '{sys.path[-1]}' has been appended to `sys.path`. If the path looks weird, you probably ran the script as a python file. Please run it as a module:")
	print("`python -m unit_test --built_inplace`")

else:
	raise RuntimeError("No arguments at all.")



