'''
DOCSTRING TO BE WRITTEN
'''

import os

fragment_defs   = {}
fragment_metadata = {}

### Runs through all definition files and adds their definitions to the global definition dictionary
### Finds list of all definition files in directory
files = os.listdir(os.path.dirname(__file__))
defs_files = []
for file in files:
    if file.endswith(".py") and not file == "__init__.py":
        ### Adds base path and removes ".py" extension
        defs_files.append("COBY.fragment_definitions." + file[:-3])

for defs_file in defs_files:
    defs = __import__(defs_file, fromlist=[''])

    if hasattr(defs, "fragment_defs"):
        fragment_defs.update(defs.fragment_defs)

    if hasattr(defs, "fragment_metadata"):
        fragment_metadata.update(defs.fragment_metadata)
