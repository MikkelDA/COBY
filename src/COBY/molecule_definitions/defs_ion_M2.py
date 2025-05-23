pos_ion_defs = {}
neg_ion_defs = {}
pos_ion_metadata = {}
neg_ion_metadata = {}

params = "M2"
pos_ion_metadata[params] = {
    "Description": "Ions in this parameter library are based on the release parameters for Martini 3 as described in https://doi.org/10.1038/s41592-021-01098-3.",
}
params = "M2"
neg_ion_metadata[params] = {
    "Description": "Ions in this parameter library are based on the release parameters for Martini 3 as described in https://doi.org/10.1038/s41592-021-01098-3.",
}

############
### IONS ###
############
params = "M2"
pos_ion_defs[params] = {
    ### Monovalent
    "NA": {
        "beads": [{"name": "NA+",  "x": 0, "y":  0, "z": 0, "charge": 1}],
        "tags": ("martini3", "ion", "positive", "monovalent"),
    },
    "NC3+": {
        "beads": [{"name": "NC3",  "x": 0, "y":  0, "z": 0, "charge": 1}],
        "tags": ("martini3", "ion", "positive", "monovalent"),
    },
    "CA+": {
        "beads": [{"name": "CA+",  "x": 0, "y":  0, "z": 0, "charge": 1}],
        "tags": ("martini3", "ion", "positive", "monovalent"),
    },
}

neg_ion_defs[params] = {
    ### Monovalent
    "CL": {
        "beads": [{"name": "CL-",  "x": 0, "y":  0, "z": 0, "charge": -1}],
        "tags": ("martini3", "ion", "negative", "monovalent"),
    },
}
