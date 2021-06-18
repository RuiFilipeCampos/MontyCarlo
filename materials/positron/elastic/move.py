

#from ....settings import __montecarlo__
import subprocess
import os
from multiprocessing import Pool

#__path__ = __montecarlo__/'materials'/'positron'/'elastic'


#__path__ = str(__path__)


from os import listdir

import shutil


for Z in range(2, 100):
    path =  f"./{Z}/"
    final_path =  f"./positron_elastic/{Z}/"
    files = listdir(path)
    os.mkdir(final_path)
    for file in files[0:-1]:
    	if "dcs_" in file:
    		#print(file)
        	shutil.move(path + file, final_path + file)
    
    

