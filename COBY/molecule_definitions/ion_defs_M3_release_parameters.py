ion_defs     = {}

############
### IONS ###
############
params = "default"
ion_defs[params] = {}
ion_defs[params]["positive"] = {
    ### Monovalent
    "NA":  {"beads": "NA",  "charge": 1, "x": (0,), "y": (0,), "z": (0,)},
    "K":   {"beads": "K",   "charge": 1, "x": (0,), "y": (0,), "z": (0,)},
    "TMA": {"beads": "TMA", "charge": 1, "x": (0,), "y": (0,), "z": (0,)},
    ### Divalent
    "CA":  {"beads": "CA",  "charge": 2, "x": (0,), "y": (0,), "z": (0,)},
}

ion_defs[params]["negative"] = {
    ### Monovalent
    "CL":   {"beads": "CL",  "charge": -1, "x": (0,), "y": (0,), "z": (0,)},
    "BR":   {"beads": "BR",  "charge": -1, "x": (0,), "y": (0,), "z": (0,)},
    "IOD":  {"beads": "ID",  "charge": -1, "x": (0,), "y": (0,), "z": (0,)},
    "ACE":  {"beads": "CL",  "charge": -1, "x": (0,), "y": (0,), "z": (0,)},
    "BF4":  {"beads": "BF4", "charge": -1, "x": (0,), "y": (0,), "z": (0,)},
    "PF6":  {"beads": "PF6", "charge": -1, "x": (0,), "y": (0,), "z": (0,)},
    "SCN":  {"beads": "SCN", "charge": -1, "x": (0,), "y": (0,), "z": (0,)},
    "CLO4": {"beads": "CLO", "charge": -1, "x": (0,), "y": (0,), "z": (0,)},
    "NO3":  {"beads": "NO3", "charge": -1, "x": (0,), "y": (0,), "z": (0,)},
}
