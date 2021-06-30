

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
