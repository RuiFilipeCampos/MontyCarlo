from ._init import eax 




print("IMPORTING GEOMETRY")
from .geometry.main import *
from .geometry.CSG import Sphere

print("IMPORTING MATERIAL")
from .materials.materials import *
from .sources import *

from matplotlib.pyplot import *
from numpy import *

from .materials import database

# # import numpy as np
# # import pyximport
# # pyximport.install(reload_support= True,
# #                   inplace       = True,
# #                   setup_args = {'include_dirs':np.get_include()})


# print("____INIT_____")

# print("from .materials.materials import eax")
# from .materials.materials import eax


