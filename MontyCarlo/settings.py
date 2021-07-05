__doc__ = """
settings.py
  This was originally meant to be an editable file that would change simulation paramaters. 
  It will be deprecated in the future. VERY important parameters defined here:
  
  - photon cut off            :: below this energy value the photon simulation is stopped
  - electron/positron cut off :: below this energy value the simulation of charged particles are stopped
  - DEBUG                     :: parameter that defines run time debug status
  - SURFACE_THICKNESS         :: completely ignore this, it was for a very old version of the geometry module (exclusive sphere tracer)
  - Wcc                       :: completely ignore this, this is very old, the Wcc and Wcr parameters are directly defined when creating a material
  
How this file will change:
  Similar to Wcc and Wcr parameters, the particle energy cut offs will be defined in the Material instance or Volume instance. If they are not defined, the cut offs should
  then be automatically set equal to their condensed history counterparts(Wcr and Wcc).
  
  Printing won't be as messy!
  
  __montecarlo__ is also a very old name, should be changed to __path__ or __PATH__. 
"""


print("SHOWING setting.py")

#__PATH__ = repr(__file__)[:-14][1:]

from pathlib import Path
__montecarlo__ = Path(repr(__file__)[1:])
__montecarlo__ = __montecarlo__.parent


print("_________________________________________________________________")
print("INSTALL PATH: ", __montecarlo__)
print("----")

DEBUG = False # runtime debug mode

__photonCUTOFF__   = 10e3  
__electronCUTOFF__ = 100e3

SURFACE_THICKNESS = .001e-4  #500 um --- this has no meaning, keeping it in case some older module is calling it
Wcc = 10 #eV < --- no meaning

print(f"DEBUG: {DEBUG}")
print("---------------")

print(f"PHOTON CUT OFF: {__photonCUTOFF__*0.511} MeV")
print("----")

print(f"ELECTRON CUT OFF: {__electronCUTOFF__} MeV")
print("---------------------------------------------")

print("");print("");



Wcc = Wcc*1e-6 



print("BRANCH: default1 - pre-alpha")
