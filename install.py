
print("Downloading EADL...")
from .materials.EADL import download_database
print("Downloading EPDL...")
from .materials.EPDL import download_database
print("Downloading EEDL...")
from .materials.EEDL import download_database

print("Downloading electron elastic data...")
from .materials.electron import download_database
print("Downloading electron positron data...")
from .materials.positron import download_database

from .settings import __montecarlo__

#print(f"""
#IN THIS VERSION OF MONTY CARLO YOU HAVE TO MANUALLY FINISH THE INSTALL.
#
#PLEASE GO TO THE FOLLOWING FOLDERS, FIND THE .tar FILES AND EXTRACT THEM:
#
#{str(__montecarlo__/'materials'/'electron')} 
#
#{str(__montecarlo__/'materials'/'positron')} 
#
#
#""")
#
#a = "n"
#while a != "y":
#	a = input(">>>> continue? [y]/[N]")
#
#del a#