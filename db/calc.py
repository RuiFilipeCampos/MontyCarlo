import subprocess
import os
from multiprocessing import Pool







def replace_first_line( src_filename, target_filename, replacement_line):
    f = open(src_filename)
    first_line, remainder = f.readline(), f.read()
    t = open(target_filename,"w")
    t.write(replacement_line + "\n")
    t.write(remainder)
    t.close()
    



def make(Z):
    subprocess.call(f"md {Z}", shell=True)
    rep = f"IZ      {Z}         atomic number                               [none]"
    replace_first_line("minimal.in", "minimal.in", rep)
    subprocess.call(r"elscata < minimal.in", shell=True)
    subprocess.call(f"move dcs*.DAT {Z}/", shell=True)
    

for Z in range(1, 25):
    make(Z)

#if __name__ == '__main__':
#    with Pool(5) as p:
#        p.map(make, range(1, 10))




#a = shlex.split(r"elscata < minimal.in")

#print(a)
#subprocess.run(a)

#return_code = subprocess.call(r"elscata < minimal.in", shell=True)  
#print(return_code)

