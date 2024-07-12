ion_defs     = {}

############
### IONS ###
############
params = "default"
ion_defs[params] = {}
ion_defs[params]["positive"] = {
    ### Monovalent
    "NA":  {"beads": [{"name": "NA",  "x": 0, "y":  0, "z": 0, "charge": 1}]},
    "K":   {"beads": [{"name": "K",   "x": 0, "y":  0, "z": 0, "charge": 1}]},
    "TMA": {"beads": [{"name": "TMA", "x": 0, "y":  0, "z": 0, "charge": 1}]},
    ### Divalent
    "CA":  {"beads": [{"name": "CA",  "x": 0, "y":  0, "z": 0, "charge": 2}]},
}

ion_defs[params]["negative"] = {
    ### Monovalent
    "CL":   {"beads": [{"name": "CL",  "x": 0, "y":  0, "z": 0, "charge": -1}]},
    "BR":   {"beads": [{"name": "BR",  "x": 0, "y":  0, "z": 0, "charge": -1}]},
    "IOD":  {"beads": [{"name": "ID",  "x": 0, "y":  0, "z": 0, "charge": -1}]},
    "ACE":  {"beads": [{"name": "CL",  "x": 0, "y":  0, "z": 0, "charge": -1}]},
    "BF4":  {"beads": [{"name": "BF4", "x": 0, "y":  0, "z": 0, "charge": -1}]},
    "PF6":  {"beads": [{"name": "PF6", "x": 0, "y":  0, "z": 0, "charge": -1}]},
    "SCN":  {"beads": [{"name": "SCN", "x": 0, "y":  0, "z": 0, "charge": -1}]},
    "CLO4": {"beads": [{"name": "CLO", "x": 0, "y":  0, "z": 0, "charge": -1}]},
    "NO3":  {"beads": [{"name": "NO3", "x": 0, "y":  0, "z": 0, "charge": -1}]},
}
