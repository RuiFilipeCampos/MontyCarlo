import argparse

myco = argparse.ArgumentParser()


subparser = myco.add_subparsers(dest='command')


# Creating a project
new = subparser.add_parser("new", help = "Create a new MontyCarlo project.")
new.add_argument("project_name", type=str)







def main():
  args = myco.parse_args()


  if args.command is None:
    pass

  elif args.command == "new":
    pass

