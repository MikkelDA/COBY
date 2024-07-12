solvent_defs = {}
ion_defs = {}
### molar_mass: [g/mol]
### density: [g/cm^3]

################################
### ATOMISTIC WATER AND IONS ###
################################

### I suggest adding the following subarguments to solvent arguments when using atomistic solvent
### "kick:0.0425 bead_radius:0.17"
### E.g. "solv:SOL:params:AA_tip4p pos:SOD:params:AA neg:CLA:params:AA kick:0.0425 bead_radius:0.17"
### The bead radius is a (very) rough estimation from table 3 in https://pubs.acs.org/doi/10.1021/jp953141%2B
### The kick is simply bead_radius/4

### Atom positions
#      O
#  H1     H2
params = "AA_tip3p"
solvent_defs[params] = {
    "SOL" : {
        "mapping_ratio": 1, "density": 0.99669, "molar_mass": 18.01528,
        "beads": [
            {"name": "OW",  "x":  0,     "y": 0,     "z": 0, "charge": -0.834},
            {"name": "HW1", "x":  0.074, "y": 0.064, "z": 0, "charge": 0.417},
            {"name": "HW2", "x": -0.074, "y": 0.064, "z": 0, "charge": 0.417},
        ],
    },
}

### Atom positions
#      O
#      M
#  H1     H2
params = "AA_tip4p"
solvent_defs[params] = {
    "SOL" : {
        "mapping_ratio": 1, "density": 0.99669, "molar_mass": 18.01528,
        "beads": [
            {"name": "OW",  "x":  0,     "y": 0,     "z": 0, "charge":  0},
            {"name": "HW1", "x":  0.074, "y": 0.064, "z": 0, "charge":  0.52},
            {"name": "HW2", "x": -0.074, "y": 0.064, "z": 0, "charge":  0.52},
            {"name": "MW",  "x":  0,     "y": 0.032, "z": 0, "charge": -1.04},
        ],
    },
}

### Atom positions
#   LP1 LP2
#      O
#  H1     H2
params = "AA_tip5p"
solvent_defs[params] = {
    "SOL" : {
        "mapping_ratio": 1, "density": 0.99669, "molar_mass": 18.01528,
        "beads": [
            {"name": "OW",  "x":  0,     "y":  0,     "z": 0, "charge":  0},
            {"name": "HW1", "x":  0.074, "y":  0.064, "z": 0, "charge":  0.241},
            {"name": "HW2", "x": -0.074, "y":  0.064, "z": 0, "charge":  0.241},
            {"name": "LP1", "x":  0.02,  "y":  0.02,  "z": 0, "charge": -0.241},
            {"name": "LP1", "x": -0.02,  "y": -0.02,  "z": 0, "charge": -0.241},
        ],
    },
}

params = "AA"
ion_defs[params] = {}
ion_defs[params]["positive"] = {
    ### Monovalent
    "SOD" : {"beads": [{"name": "SOD", "x": 0, "y":  0, "z": 0, "charge": 1}]},
}

ion_defs[params]["negative"] = {
    ### Monovalent
    "CLA" : {"beads": [{"name": "CLA", "x": 0, "y":  0, "z": 0, "charge": -1}]},
}
