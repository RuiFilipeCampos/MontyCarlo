import requests
from ...settings import __montecarlo__

PATH = __montecarlo__/'materials'/'EPDL'

for N in range(1, 101):
    #url = "https://www-nds.iaea.org/epics/ENDL2017/EADL.ELEMENTS/ZA00" + str(N) + "000"
    url = "https://ruifilipecampos.github.io/MontyCarlo/EPDL/" + str(N) + ".txt"
    file = requests.get(url)
    filename = str(N) + ".txt"
    open(str(PATH/filename), 'wb').write(file.content)


"""
for N in range(10, 100):
    url = "https://www-nds.iaea.org/epics/ENDL2017/EADL.ELEMENTS/ZA0" + str(N) + "000"
    file = requests.get(url)
    filename = str(N) + ".txt"
    open(str(PATH/filename), 'wb').write(file.content)

url = "https://www-nds.iaea.org/epics/ENDL2017/EADL.ELEMENTS/ZA100000"
file = requests.get(url)
filename = "100.txt"
open(str(PATH/filename), 'wb').write(file.content)
"""
