print(">>>> SETTINGS.Py")


#################################################################
#################################################################
__PATH__ = r"C:\Users\Rui Campos\Dropbox\MontyCarlo\MontyCarlo" #<------- ONLY CHANGE THIS 
#################################################################
#################################################################


print("_________________________________________________________________")
print("INSTALL PATH: ", __PATH__)
print("----")

DEBUG = False

__photonCUTOFF__   = 10e3
__electronCUTOFF__ = 100e3

SURFACE_THICKNESS = .001e-4  #500 um


#__photonCUTOFF__ = __photonCUTOFF__*1e-6 / 0.511


Wcc = 10 #eV




print(f"SURFACE_THICKNESS: {SURFACE_THICKNESS}")
print("---------------")


print(f"DEBUG: {DEBUG}")
print("---------------")

print(f"PHOTON CUT OFF: {__photonCUTOFF__*0.511} MeV")
print("----")

print(f"ELECTRON CUT OFF: {__electronCUTOFF__} MeV")
print("---------------------------------------------")

print(f"Wcc: {Wcc} eV")
print("---------------")

print("");print("");






Wcc = Wcc*1e-6 

from pathlib import Path
__montecarlo__ = Path(__PATH__)
