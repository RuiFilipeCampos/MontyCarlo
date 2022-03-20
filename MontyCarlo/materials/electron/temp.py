

from ...settings import __montecarlo__


__path__ = __montecarlo__/'materials'/'electron'


__path__ = str(__path__)

from numpy import pi


def getData(Z):
    from numpy import array, zeros
    from os import listdir
    path = __path__ + f"/{Z}/"
    files = listdir(path)
    print(path)
    #files = listdir(path)
    for f in files:
        if ".pkl" in f:
            files.remove(f)

    
    eax = []
    
    for file in files:
        e = file[4:12].replace("p", ".")
        e = float(e)
        eax.append(e)
    
    n = len(eax)
    print(n)
    
    files = [f for _, f in sorted(zip(eax, files))]
    
    allDCS = zeros((n, 606))
    SIGMA0 = zeros(n)
    SIGMA1 = zeros(n)
    SIGMA2 = zeros(n)
    
    for j, file in enumerate(files):

        with open(path + file, "r") as f:
            lines = f.readlines()
        
        DCS = zeros(606)

        for line in lines:
            if line[1] == "#":
                if "Total elastic cross section" in line:
                    SIGMA0[j] =  float(line.split()[6]) 
                if "1st transport cross section" in line:
                    SIGMA1[j] = float(line.split()[6])
                if "2nd transport cross section" in line:
                    SIGMA2[j] = float(line.split()[6]) 
                    i = 0
                    continue
 
                continue
            numb = line.split()
            dcs = numb[2]
            DCS[i] = float(dcs)
            i += 1

        allDCS[j] = DCS*4*pi
    eax.sort()
    return allDCS, SIGMA0, SIGMA1, SIGMA2, array(eax)
    
from numpy import *
    
allDCS, SIGMA0, SIGMA1, SIGMA2, eax = getData(1)
TCS = array([SIGMA0, SIGMA1, SIGMA2])

import os
def getGrid():
    from numpy import array
    from os import listdir
    path = __path__ + f"/{1}/"

    files = listdir(path)
    for f in files:
        if ".pkl" in f:
            files.remove(f)

    
    
    allMU  = []

    
    with open(path + files[0], "r") as f:
        lines = f.readlines()

    MU = []
  
    for line in lines:
        if line[1] == "#":
            continue
                
        numb = line.split()
        mu = numb[1]
        MU.append(float(mu))
    return array(MU)
#os.makedirs(__path__ + "/elastic")

mu = getGrid()


path = __path__ + f"/elastic/"
with open(path + "muGRID.npy", 'wb') as f:
    save(f, mu)
        
        
        
        
#for Z in range(1, 100):

#    path = __path__ + f"/elastic/{Z}/"
    #os.makedirs(path)
#    with open(path + "HEtransportTCS.npy", 'wb') as f:
 #       save(f, TCS)
#    with open(path + "DCS.npy", 'wb') as f:
 #       save(f, allDCS)
    
    
    

    
