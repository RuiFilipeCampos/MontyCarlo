__doc__ = """Setting up a new version.
"""

__author__ = "Rui Campos"

version = "0.0.41" # to be imported by setup*.py scripts


with open("setup.cfg", "w") as setup_cfg:
  text = setup_cfg.readlines()
  print(text)

### to be added: open files and write version
