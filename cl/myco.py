import argparse
import sys
from pathlib import Path

print(Path(sys.argv[0]))
import cmd, sys
from turtle import *
import colorama
from colorama import Fore, Back, Style
colorama.init()

# Set the color semi-permanently
print(Fore.CYAN)
print("Text will continue to be cyan")
print("until it is reset or changed")
print(Style.RESET_ALL)

# Colorize a single line and then reset
print(Fore.RED + 'You can colorize a single line.' + Style.RESET_ALL)

# Colorize a single word in the output
print('Or a single ' + Back.GREEN + 'words' + Style.RESET_ALL + ' can be highlighted')

# Combine foreground and background color
print(Fore.BLUE + Back.WHITE)
print('Foreground, background, and styles can be combined')
print("==========            ")

print(Style.RESET_ALL)
print('If unsure, reset everything back to normal.')
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

    def do_open(self, args):
        project_path = NavigationCMD.path/args
        if project_path.is_dir():
            if is_project(project_path):
                MontyCarloProjectShell().cmdloop(0, 0)
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
        (root_folder/'.myco').touch()



        print("The name of your project is:", args)
        print("\n")


    def do_run():
        pass

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

MontyCarloShell().cmdloop()
