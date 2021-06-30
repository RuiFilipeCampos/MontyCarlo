import requests

from ...settings import __montecarlo__

PATH = __montecarlo__/'materials'/'positron'/'elastic'

import gdown

url = 'https://drive.google.com/uc?id=1V2cB4s4YenPFBaHoOnnn3rYFpnlWulDl'
output = str(PATH) + ".tar"
gdown.download(url, output, quiet=False)



from pyunpack import Archive

Archive(output).extractall(str(__montecarlo__/'materials'/'positron'))
