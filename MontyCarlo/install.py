__doc__ = """
install.py
    Second phase of installation. Triggered on first import. Downloads the following databases:
      - EADL :: Evaluated Atomic Data Library :: downloaded with a webscraper <---- this will eventually change and the data base will be downloaded from my drive
      - EPDL :: Evaluated Photon Data Library :: downloaded with a webscraper <---- this will eventually change and the data base will be downloaded from my drive
      - EEDL :: Evaluated Electron Data Library :: downloaded with a webscraper <---- this will eventually change and the data base will be downloaded from my drive
      - Data regarding the elastic scattering of electrons, compiled using ELSPA :: downloaded from my personal drive
      - Data regarding the elastic scattering of positrons, compiled using ELSPA :: downloaded from my personal drive
"""

print("Downloading EADL...")
from .materials.EADL import download_database
print("""Done!
""")

print("Downloading EPDL...")
from .materials.EPDL import download_database
print("""Done!
""")

print("Downloading EEDL...")
from .materials.EEDL import download_database
print("""Done!
""")

print("Downloading electron elastic data...")
from .materials.electron import download_database
print("""Done!
""")

print("Downloading electron positron data...")
from .materials.positron import download_database
print("""Done!
""")
