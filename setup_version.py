__doc__ = """Setting up a new version.
"""

__author__ = "Rui Campos"



version = "0.0.5dev1" # to be imported by setup*.py scripts


with open("setup.cfg", "r") as setup_cfg:
  text = setup_cfg.readlines()

text[3] = f"version = {version}\n"


print(text)
with open("setup.cfg", "w") as setup_cfg:
  new_text = ""
  for line in text:
  	new_text += line
  setup_cfg.write(new_text)


### to be added: open files and write version
