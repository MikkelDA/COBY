'''
DOCSTRING TO BE WRITTEN
'''

import os

### 'lipid_scaffolds' somewhat follows the old insane method of creating generalized structures from which many variations of which can be made
### 'lipid_scaffolds' is limited to single-residue molecules. If you have multi-residue lipids then use 'lipid_defs' instead
lipid_scaffolds = {}
### 'lipid_defs', 'solvent_defs' and 'ion_defs' use a new method whereby multiple residues can be designated
lipid_defs      = {}
solvent_defs    = {}
ion_defs        = {}
prot_defs       = {}

### Runs through all definition files and adds their definitions to the global definition dictionaries
### Finds list of all definition files in directory
files = os.listdir(os.path.dirname(__file__))
defs_files = []
for file in files:
    if file.endswith(".py") and not file == "__init__.py":
        ### Adds base path and removes ".py" extension
        defs_files.append("COBY.molecule_definitions." + file[:-3])

for defs_file in defs_files:
    defs = __import__(defs_file, fromlist=[''])

    if hasattr(defs, "lipid_scaffolds"):
        lipid_scaffolds.update(defs.lipid_scaffolds)
    if hasattr(defs, "lipid_defs"):
        lipid_defs.update(defs.lipid_defs)
    if hasattr(defs, "solvent_defs"):
        solvent_defs.update(defs.solvent_defs)
    if hasattr(defs, "ion_defs"):
        ion_defs.update(defs.ion_defs)
    if hasattr(defs, "prot_defs"):
        prot_defs.update(defs.prot_defs)
