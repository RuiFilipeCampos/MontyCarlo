__doc__ = """Download elastic scattering database from `https://ruifilipecampos.github.io/MontyCarlo/`.
"""
__author__ = "Rui Campos"


import requests
import os

"""
from ...settings import __montecarlo__

PATH = __montecarlo__/'materials'/'electron'/'elastic'
elastic_folder = str(PATH)
"""

os.mkdir(elastic_folder)

url = r"https://ruifilipecampos.github.io/MontyCarlo/elastic/electron"
top_level_files = ["HEeax.npy", "LEeax.npy", "muGRID.npy"]

for filename in top_level_files:
  file = requests.get(url + filename)
  print(file)
  
