import requests
from ...settings import __montecarlo__

PATH = __montecarlo__/'materials'/'EPDL'
for N in range(1, 10):
    url = "https://www-nds.iaea.org/epics/ENDL2017/EPDL.ELEMENTS/ZA00" + str(N) + "000"
    file = requests.get(url)
    filename = str(N) + ".txt"
    open(str(PATH/filename), 'wb').write(file.content)


for N in range(10, 100):
    url = "https://www-nds.iaea.org/epics/ENDL2017/EPDL.ELEMENTS/ZA0" + str(N) + "000"
    file = requests.get(url)
    filename = str(N) + ".txt"
    open(str(PATH/filename), 'wb').write(file.content)

url = "https://www-nds.iaea.org/epics/ENDL2017/EPDL.ELEMENTS/ZA100000"
file = requests.get(url)
filename = "100.txt"
open(str(PATH/filename), 'wb').write(file.content)
