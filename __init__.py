#from ._init import eax 

from . import _init



#print("IMPORTING GEOMETRY")
from .geometry.main import *
from .geometry.CSG import Sphere



#print("IMPORTING MATERIAL")
from .materials.materials import *


from .sources import *


from matplotlib.pyplot import *
from numpy import *

from .materials import database


