import requests

from ...settings import __montecarlo__

PATH = __montecarlo__/'materials'/'electron'/'elastic'

import gdown

url = 'https://drive.google.com/uc?id=1wod967aO8K90AtW9nvgtL8irRr0Go9iE'
output = str(PATH) + ".tar"
gdown.download(url, output, quiet=False) 


