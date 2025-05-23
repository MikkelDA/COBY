pos_ion_defs = {}
neg_ion_defs = {}
pos_ion_metadata = {}
neg_ion_metadata = {}

params = "default"
pos_ion_metadata[params] = {
    "Description": "Ions in this parameter library are based on the release parameters for Martini 3 as described in https://doi.org/10.1038/s41592-021-01098-3.",
}
params = "default"
neg_ion_metadata[params] = {
    "Description": "Ions in this parameter library are based on the release parameters for Martini 3 as described in https://doi.org/10.1038/s41592-021-01098-3.",
}

############
### IONS ###
############
params = "default"
pos_ion_defs[params] = {
    ### Monovalent
    "NA":  {
        "beads": [{"name": "NA",  "x": 0, "y":  0, "z": 0, "charge": 1}],
        "tags": ("martini3", "ion", "positive", "monovalent"),
    },
    "K":   {
        "beads": [{"name": "K",   "x": 0, "y":  0, "z": 0, "charge": 1}],
        "tags": ("martini3", "ion", "positive", "monovalent"),
    },
    "TMA": {
        "beads": [{"name": "TMA", "x": 0, "y":  0, "z": 0, "charge": 1}],
        "tags": ("martini3", "ion", "positive", "monovalent"),
    },
    ### Divalent
    "CA": {
        "beads": [{"name": "CA",  "x": 0, "y":  0, "z": 0, "charge": 2}],
        "tags": ("martini3", "ion", "positive", "divalent"),
    },
}

neg_ion_defs[params] = {
    ### Monovalent
    "CL": {
        "beads": [{"name": "CL",  "x": 0, "y":  0, "z": 0, "charge": -1}],
        "tags": ("martini3", "ion", "negative", "monovalent"),
    },
    "BR": {
        "beads": [{"name": "BR",  "x": 0, "y":  0, "z": 0, "charge": -1}],
        "tags": ("martini3", "ion", "negative", "monovalent"),
    },
    "IOD": {
        "beads": [{"name": "ID",  "x": 0, "y":  0, "z": 0, "charge": -1}],
        "tags": ("martini3", "ion", "negative", "monovalent"),
    },
    "ACE": {
        "beads": [{"name": "CL",  "x": 0, "y":  0, "z": 0, "charge": -1}],
        "tags": ("martini3", "ion", "negative", "monovalent"),
    },
    "BF4": {
        "beads": [{"name": "BF4", "x": 0, "y":  0, "z": 0, "charge": -1}],
        "tags": ("martini3", "ion", "negative", "monovalent"),
    },
    "PF6": {
        "beads": [{"name": "PF6", "x": 0, "y":  0, "z": 0, "charge": -1}],
        "tags": ("martini3", "ion", "negative", "monovalent"),
    },
    "SCN": {
        "beads": [{"name": "SCN", "x": 0, "y":  0, "z": 0, "charge": -1}],
        "tags": ("martini3", "ion", "negative", "monovalent"),
    },
    "CLO4": {
        "beads": [{"name": "CLO", "x": 0, "y":  0, "z": 0, "charge": -1}],
        "tags": ("martini3", "ion", "negative", "monovalent"),
    },
    "NO3": {
        "beads": [{"name": "NO3", "x": 0, "y":  0, "z": 0, "charge": -1}],
        "tags": ("martini3", "ion", "negative", "monovalent"),
    },
}
