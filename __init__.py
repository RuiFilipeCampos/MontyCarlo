import os


from .settings import __montecarlo__


with open(str(__montecarlo__/'DATA_STATUS.txt'), "r") as file:
    for status in file:
        continue

status = int(status)
if status == 0:
    from . import install
    with open(str(__montecarlo__/'DATA_STATUS.txt'), "w") as file:
        file.write("1")
    status = None




from ._init import eax
from .plotter import *
#from MontyCarlo.init import eax 




#print("IMPORTING GEOMETRY")
from MontyCarlo.geometry.main import *
from .geometry.CSG import Sphere



#print("IMPORTING MATERIAL")
from .materials.materials import *


from .sources import *


from matplotlib.pyplot import *
from numpy import *

from .materials import database


