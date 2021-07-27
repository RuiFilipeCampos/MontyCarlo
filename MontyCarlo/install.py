__doc__ = """
install.py
    Second phase of installation. Triggered on first import. Downloads the following databases:
      - EADL :: Evaluated Atomic Data Library   (*.txt)
      - EPDL :: Evaluated Photon Data Library   (*.txt)
      - EEDL :: Evaluated Electron Data Library (*.txt)
      - Data regarding the elastic scattering of electrons, compiled using ELSPA (*.npy) 
      - Data regarding the elastic scattering of positrons, compiled using ELSPA (*.npy)
"""


print("Importing .install")


print("""
The data being download can be seen in `https://github.com/RuiFilipeCampos/MontyCarlo/tree/pre-alpha/0.0.41/docs` 
in the `EADL`, `EPDL`, `EEDL`, and `elastic` folders. 
""")


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
