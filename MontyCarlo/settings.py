print("SHOWING setting.py")

__PATH__ = repr(__file__)[:-14][1:]

print("_________________________________________________________________")
print("INSTALL PATH: ", __PATH__)
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

from pathlib import Path

__montecarlo__ = Path(__PATH__)

print("BRANCH: default1 - pre-alpha")
