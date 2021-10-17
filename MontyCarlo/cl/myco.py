BOILER_PLATE = """


# Imports will be cleaner eventually...
from MontyCarlo import *
from MontyCarlo.sources import *
from MontyCarlo.geometry.CSG import *


# Define your materials. They are compiled once per project, so don't worry if it takes a while
# the first time you run the script, the second time `Mat` will just load from cache.

water = Mat({1:2, 8:1}, 1, name = "Water")

air = Mat({6:1.50187000E-04,
           7:7.84430000E-01,
           8:2.10748000E-01,
           18:4.67111000E-03}, 
           1.20479e-03,
           C1 = .2, 
           C2 = .2,
           name = "AirDryNearSeaLevel")

gold = Mat({79:1}, 1.93200000E+01, name = "Gold")


# Define geometry
# The indentation tells you exactly how the BVH is being constructed.
# The `.configure` method will be replaced with a smarter system.
with InfiniteVolume() as outer:
    outer.configure("OUTER", render = False)
    outer.fill(gold) 
    
    with Sphere(100) as outer_sphere:
        outer_sphere in outer
        outer_sphere.configure("outer_sphere", render = True)
        outer_sphere.fill(air)

        with Sphere(50) as inner_sphere:
            inner_sphere in outer_sphere
            inner_sphere.configure("inner_sphere", render = True)
            inner_sphere.fill(water)



photon_beam = Beam(
                   "photon",       # kind of particle 
                   inner_sphere,   # initial volume
                   E = 10e6,       # initial eneryg in eV
                   N = 1_000,     # number of particles in the source, careful with this number, might break your run and fill your ram
                   pos = (0, 0, 0) # initial position 
                  ) 



# let Plotter handle the run
plotter = Plotter(photon_beam)


# then ask it for a fig
fig = plotter.new_plot() b 

# use this method to draw the geometry onto the figure (this will be better)
plotter.add_geometry(fig, outer)


fig.show()
"""




import argparse
import sys
from pathlib import Path

import cmd, sys
import colorama
from colorama import Fore, Back, Style
colorama.init()


def is_project(path):
    x = path.glob(".myco")
    x = list(x)

    if len(x):
        return x[0].is_file()
    else:
        return False
class NavigationCMD(cmd.Cmd):
    intro = "Welcome to MontyCarlo! Type help or ? to list the available commands ^.^ \n"
    file = None
    project = None
    path = Path(sys.argv[0]).parent
    prompt = f"{Fore.RED}myco{Style.RESET_ALL}@{Fore.CYAN}{path}{Style.RESET_ALL}> "

    def do_cd(self, args):
        

        if args in ["..", ""]:
            NavigationCMD.path = NavigationCMD.path.parent
            NavigationCMD.prompt = f"{Fore.RED}myco{Style.RESET_ALL}@{Fore.CYAN}{NavigationCMD.path}{Style.RESET_ALL}> "
            return 

        new_path = NavigationCMD.path/args
        if new_path.is_dir():
            NavigationCMD.path = NavigationCMD.path/args
            NavigationCMD.prompt = f"{Fore.RED}myco{Style.RESET_ALL}@{Fore.CYAN}{NavigationCMD.path}{Style.RESET_ALL}> "
        else:
            print("Error: NOT A PATH")

    def do_ld(self, args):
        x = NavigationCMD.path.glob('*')
        all_paths = list(x)
        directories = []
        files = []


        for path in all_paths:
            if path.is_dir():
                directories.append(path.name)
            else:
                files.append(path.name)
        

        directories.sort()
        files.sort()
        
        print("# Sorted alphabetically with capital letters coming first.")
        print("\nDIRECTORIES:")
        for directory in directories:
            print("    ", directory)
        print("")
        print("FILES:")
        for file in files:
            print("    ", file)
        print("")


    def do_ls(self, args):
        self.do_ld(args)




class MontyCarloShell(NavigationCMD):

    
    def do_install(self, args):
        print("# Downloading databases !")
        exec("import MontyCarlo")

        pass

    def do_uninstall(self, args):
        print("# Uninstall is not yet available !")
        pass

    def do_open(self, args):
        project_path = NavigationCMD.path/args
        if project_path.is_dir():
            if is_project(project_path):
                MontyCarloProjectShell().cmdloop(project_path)
            else:
                print("Not a project !")
        else:
            print("Not a directory !")

    ##### NAVIGATION

    def do_create(self, args):
        """Create a new MontyCarlo project.
        """

        

        if args == "":
            print("----Error------------------------------")
            print("Please provide a name for your project:")
            print("    (MontyCarlo) create name_of_my_project")
            print("--------------------------------------- \n\n")
            return
        
        root_folder = NavigationCMD.path/args
        root_folder.mkdir(parents=True, exist_ok=True)


        # DIRECTORIES
        (root_folder/'mat').mkdir(parents=True, exist_ok=True)
        (root_folder/'geo').mkdir(parents=True, exist_ok=True)
        (root_folder/'out').mkdir(parents=True, exist_ok=True)

        # FILES
        (root_folder/'.myco').touch()
        (root_folder/'main.py').touch()
        (root_folder/'main.py').write_text(BOILER_PLATE)



        print("The name of your project is:", args)
        print("\n")



    def do_exit(self, arg):
        return True



class MontyCarloProjectShell(cmd.Cmd):
    intro = "Your project is now open! \n"
    path = Path(sys.argv[0]).parent


    def do_exit(self, args):
        MontyCarloShell().cmdloop()


    def cmdloop(self, path):
        MontyCarloProjectShell.path = path
        MontyCarloProjectShell.prompt = f"{Fore.RED}myco{Style.RESET_ALL}@{Fore.GREEN}YourProject{Style.RESET_ALL}> "


        super(MontyCarloProjectShell, self).cmdloop()

    def do_run(self, args):
        code = (MontyCarloProjectShell.path/'main.py').read_text()

        class NAMESPACE:
            exec(code)
        
        del NAMESPACE
        
def main():
    MontyCarloShell().cmdloop()
