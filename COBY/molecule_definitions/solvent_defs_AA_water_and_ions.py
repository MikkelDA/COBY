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
        "beads": ("OW", "HW1",   "HW2", ), 
        "x":     ( 0,    0.074,  -0.074,),
        "y":     ( 0,    0.064,   0.064,),
        "z":     ( 0,    0,       0,    ),
        "charges": ((0, -0.834), (1, 0.417), (2, 0.417)),
        "mapping_ratio": 1, "density": 0.99669, "molar_mass": 18.01528,
    },
}

### Atom positions
#      O
#      M
#  H1     H2
params = "AA_tip4p"
solvent_defs[params] = {
    "SOL" : {
        "beads": ("OW", "HW1",   "HW2",  "MW",  ), 
        "x":     ( 0,    0.074,  -0.074,  0,    ),
        "y":     ( 0,    0.064,   0.064,  0.032,),
        "z":     ( 0,    0,       0,      0,    ),
        "charges": ((1, 0.52), (2, 0.52), (3, -1.04)),
        "mapping_ratio": 1, "density": 0.99669, "molar_mass": 18.01528,
    },
}

### Atom positions
#   LP1 LP2
#      O
#  H1     H2
params = "AA_tip5p"
solvent_defs[params] = {
    "SOL" : {
        "beads": ("OW", "HW1",   "HW2",  "LP1", "LP2",), 
        "x":     ( 0,    0.074,  -0.074,  0.02, -0.02,),
        "y":     ( 0,    0.064,   0.064, -0.02, -0.02,),
        "z":     ( 0,    0,       0,      0,     0,   ),
        "charges": ((1, 0.241), (2, 0.241), (3, -0.241), (4, -0.241)),
        "mapping_ratio": 1, "density": 0.99669, "molar_mass": 18.01528,
    },
}

params = "AA"
ion_defs[params] = {}
ion_defs[params]["positive"] = {
    ### Monovalent
    "SOD":  {"beads": "SOD",  "charge": 1, "x": (0,), "y": (0,), "z": (0,)},
}

ion_defs[params]["negative"] = {
    ### Monovalent
    "CLA":   {"beads": "CLA",  "charge": -1, "x": (0,), "y": (0,), "z": (0,)},
}
