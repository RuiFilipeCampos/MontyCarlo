import requests
import os
from ...settings import __montecarlo__

PATH = __montecarlo__/'materials'/'electron'/'elastic'
elastic_folder = str(PATH)


os.mkdir(elastic_folder)

url = r"https://ruifilipecampos.github.io/MontyCarlo/elastic/electron"
top_level_files = ["HEeax.npy", "LEeax.npy", "muGRID.npy"]

for filename in top_level_files:
  file = requests.get(url + filename)
  


"""
import gdown

url = 'https://drive.google.com/uc?id=1wod967aO8K90AtW9nvgtL8irRr0Go9iE'
output = str(PATH) + ".tar"
gdown.download(url, output, quiet=False) 

from pyunpack import Archive

Archive(output).extractall(str(__montecarlo__/'materials'/'electron'))
"""

