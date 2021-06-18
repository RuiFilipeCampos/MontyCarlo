import os

path = directory + "\pickles\\" + str(upper_dir.Z)

if not os.path.isfile(path):
        print("> *** photon/coherent: Writing form factor...")
        from .FormFactorWriter import FormFactorWriter
        FormFactorWriter(CS, upper_dir.Z)
        
import dill as pickle

with open(path, 'rb') as file:
        a1, a2, a3, a4, a5 = pickle.load(file)
        print(f"> *** photon/coherent: Unpickled paramaters of form factor of element {self.Z}")

#a1, a2, a3, a4, a5 = self.param
#Z = self.Z
