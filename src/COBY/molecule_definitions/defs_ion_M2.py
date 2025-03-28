pos_ion_defs = {}
neg_ion_defs = {}

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
