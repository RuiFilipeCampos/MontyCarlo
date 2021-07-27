__doc__ = """Download elastic scattering database from `https://ruifilipecampos.github.io/MontyCarlo/`.
"""
__author__ = "Rui Campos"


import requests
import os
from pathlib import Path

"""
from ...settings import __montecarlo__

PATH = __montecarlo__/'materials'/'electron'/'elastic'
elastic_folder = str(PATH)
"""

__path__ = Path(__file__)
__folder__ = __path__.parent
elastic_folder = str(__folder__/'elastic')


os.mkdir(elastic_folder)

url = r"https://ruifilipecampos.github.io/MontyCarlo/elastic/electron/"

top_level_files = ["HEeax.npy", "LEeax.npy", "muGRID.npy"]

for filename in top_level_files:
  with requests.get(url + filename) as file:
	  with open(str(__folder__/'elastic'/filename), 'wb') as local_file:
	  	local_file.write(file.content)


element_level_files = ["DCS.npy", "HEtransportTCS.npy", "LEtransportTCS.npy"]

for i in range(1, 100):
	element_folder = str(__folder__/'elastic'/str(i))

	os.mkdir(element_folder)

	for filename in top_level_files:
	  with requests.get(url + filename) as file:
		  with open(str(__folder__/'elastic'/str(i)/filename), 'wb') as local_file:
		  	local_file.write(file.content)



  
