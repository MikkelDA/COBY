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
pos_ion_defs    = {}
neg_ion_defs    = {}
prot_defs       = {}
### Defs for molecule fragment builder
fragment_defs   = {}

lipid_metadata    = {}
solvent_metadata  = {}
pos_ion_metadata  = {}
neg_ion_metadata  = {}
prot_metadata     = {}
fragment_metadata = {}

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
    if hasattr(defs, "pos_ion_defs"):
        pos_ion_defs.update(defs.pos_ion_defs)
    if hasattr(defs, "neg_ion_defs"):
        neg_ion_defs.update(defs.neg_ion_defs)
    if hasattr(defs, "prot_defs"):
        prot_defs.update(defs.prot_defs)
    
    if hasattr(defs, "lipid_metadata"):
        lipid_metadata.update(defs.lipid_metadata)
    if hasattr(defs, "solvent_metadata"):
        solvent_metadata.update(defs.solvent_metadata)
    if hasattr(defs, "pos_ion_metadata"):
        pos_ion_metadata.update(defs.pos_ion_metadata)
    if hasattr(defs, "neg_ion_metadata"):
        neg_ion_metadata.update(defs.neg_ion_metadata)
    if hasattr(defs, "prot_metadata"):
        prot_metadata.update(defs.prot_metadata)
