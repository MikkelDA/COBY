solvent_defs = {}

################
### SOLVENTS ###
################
### molar_mass: [g/mol]
### density: [g/cm^3]
params = "default"
solvent_defs[params] = {
    "W" : {"beads": "W",  "x": (0,), "y": (0,), "z": (0,), "mapping_ratio": 4, "density": 0.99669, "molar_mass": 18.01528},
    "SW": {"beads": "SW", "x": (0,), "y": (0,), "z": (0,), "mapping_ratio": 3, "density": 0.99669, "molar_mass": 18.01528},
    "TW": {"beads": "TW", "x": (0,), "y": (0,), "z": (0,), "mapping_ratio": 2, "density": 0.99669, "molar_mass": 18.01528},
}
### Amino acids
solvent_defs[params].update({
    ### ### 1 atom residues
    "GLY":  {"beads": ("BB",), "mapping_ratio": 1, "x": (0,), "y": (0,), "z": (0,)},
    
    ### ### 2 atom residues
    ### Residues without variants
    "ALA":  {"beads": ("BB", "SC1"), "mapping_ratio": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
    "CYS":  {"beads": ("BB", "SC1"), "mapping_ratio": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
    "VAL":  {"beads": ("BB", "SC1"), "mapping_ratio": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
    "LEU":  {"beads": ("BB", "SC1"), "mapping_ratio": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
    "ILE":  {"beads": ("BB", "SC1"), "mapping_ratio": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
    "MET":  {"beads": ("BB", "SC1"), "mapping_ratio": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
    "PRO":  {"beads": ("BB", "SC1"), "mapping_ratio": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
    "HYP":  {"beads": ("BB", "SC1"), "mapping_ratio": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
    "ASN":  {"beads": ("BB", "SC1"), "mapping_ratio": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
    "GLN":  {"beads": ("BB", "SC1"), "mapping_ratio": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
    "THR":  {"beads": ("BB", "SC1"), "mapping_ratio": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
    "SER":  {"beads": ("BB", "SC1"), "mapping_ratio": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
    ### ASP variants
    "ASP":  {"beads": ("BB", "SC1"), "mapping_ratio": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0), "charges": (1, -1)},
    "ASPP": {"beads": ("BB", "SC1"), "mapping_ratio": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)}, # Neutral ASP
    "ASH":  {"beads": ("BB", "SC1"), "mapping_ratio": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)}, # Neutral ASP
    ### Glutamate variants
    "GLU":  {"beads": ("BB", "SC1"), "mapping_ratio": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0), "charges": (1, -1)},
    "GLUP": {"beads": ("BB", "SC1"), "mapping_ratio": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)}, # Neutral GLU
    "GLH":  {"beads": ("BB", "SC1"), "mapping_ratio": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)}, # Neutral GLU
    
    ### ### 3 atom residues
    "ARG":  {"beads": ("BB", "SC1", "SC2"), "mapping_ratio": 1, "x": (0.25, 0, -0.25), "y": (0, 0, 0.125), "z": (0, 0, 0), "charges": (2, 1)},
    ### Lysine variants
    "LYS":  {"beads": ("BB", "SC1", "SC2"), "mapping_ratio": 1, "x": (0.25, 0, -0.25), "y": (0, 0, 0.125), "z": (0, 0, 0), "charges": (2, 1)},
    "LSN":  {"beads": ("BB", "SC1", "SC2"), "mapping_ratio": 1, "x": (0.25, 0, -0.25), "y": (0, 0, 0.125), "z": (0, 0, 0)}, # Neutral LYS
    "LYN":  {"beads": ("BB", "SC1", "SC2"), "mapping_ratio": 1, "x": (0.25, 0, -0.25), "y": (0, 0, 0.125), "z": (0, 0, 0)}, # Neutral LYS
    
    ### ### 4 atom residues
    "PHE":  {"beads": ("BB", "SC1", "SC2", "SC3"), "mapping_ratio": 1, "x": (0.25, 0, -0.25, -0.25), "y": (0, 0, 0.125, -0.125), "z": (0, 0, 0, 0)},
    ### Histidines
    "HIS":  {"beads": ("BB", "SC1", "SC2", "SC3"), "mapping_ratio": 1, "x": (0.25, 0, -0.25, -0.25), "y": (0, 0, 0.125, -0.125), "z": (0, 0, 0, 0)},
    "HIE":  {"beads": ("BB", "SC1", "SC2", "SC3"), "mapping_ratio": 1, "x": (0.25, 0, -0.25, -0.25), "y": (0, 0, 0.125, -0.125), "z": (0, 0, 0, 0)},
    "HSE":  {"beads": ("BB", "SC1", "SC2", "SC3"), "mapping_ratio": 1, "x": (0.25, 0, -0.25, -0.25), "y": (0, 0, 0.125, -0.125), "z": (0, 0, 0, 0)},
    "HSD":  {"beads": ("BB", "SC1", "SC2", "SC3"), "mapping_ratio": 1, "x": (0.25, 0, -0.25, -0.25), "y": (0, 0, 0.125, -0.125), "z": (0, 0, 0, 0)},
    "HID":  {"beads": ("BB", "SC1", "SC2", "SC3"), "mapping_ratio": 1, "x": (0.25, 0, -0.25, -0.25), "y": (0, 0, 0.125, -0.125), "z": (0, 0, 0, 0)},
    "HSP":  {"beads": ("BB", "SC1", "SC2", "SC3"), "mapping_ratio": 1, "x": (0.25, 0, -0.25, -0.25), "y": (0, 0, 0.125, -0.125), "z": (0, 0, 0, 0), "charges": ((2, 0.5), (3, 0.5))},
    "HIP":  {"beads": ("BB", "SC1", "SC2", "SC3"), "mapping_ratio": 1, "x": (0.25, 0, -0.25, -0.25), "y": (0, 0, 0.125, -0.125), "z": (0, 0, 0, 0), "charges": ((2, 0.5), (3, 0.5))},
    
    ### ### 5 atom residues
    "TYR":  {"beads": ("BB", "SC1", "SC2", "SC3", "SC4"), "mapping_ratio": 1, "x": (0.25, 0.25, 0, 0, -0.25), "y": (0.125, 0, -0.125, 0.125, 0), "z": (0, 0, 0, 0, 0)},
    
    ### ### 6 atom residues
    "TRP":  {"beads": ("BB", "SC1", "SC2", "SC3", "SC4", "SC5"), "mapping_ratio": 1, "x": (0.25, 0.25, 0, 0, -0.25, -0.25), "y": (0.125, 0, -0.125, 0.125, 0, 0.125), "z": (0, 0, 0, 0, 0, 0)},
})

### Example of a multi-residue solvent molecule
### Two different ways to designate charges are shown
### ### Can be done within each residue dictionary or for the whole molecule
### ### Note that python-numbering starts from 0, meaning that 0 is the first bead (and 1 is the second)
#solvent_defs[params].update({
#     "Example_molecule": { # Name of molecule (used when selecting molecule in solvation/flooding commands)
#         "residues": [ # The list of residue dictionaries
#             { # Residue 0 (python indexing starts from zero)
#                 "resname": "RES0",
#                 "beads":   ("R0A0", "R0A1", "R0A2", "R0A3", "R0A4"),
#                 "x":       (     0,  0.125,   0.25,  0.375,    0.5),
#                 "y":       (     0,      0,  0.125,      0, -0.125),
#                 "z":       (     0,  0.125,      0,  0.125,      0),
#                 "charges": (1, 1), # Bead nr 1 has a charge of 1 (E.g. Bead "R0A1")
#             },
#             { # Residue 1
#                 "resname": "RES1",
#                 "beads":   ("R1A0", "R1A1", "R1A2", "R1A3", "R1A4"),
#                 "x":       (     0,  0.125,   0.25,  0.375,    0.5),
#                 "y":       (-0.125,      0,      0, -0.125,      0),
#                 "z":       ( 0.125,   0.25,  0.125,   0.25,  0.125),
#                 "charges": (4, -1), # Bead nr 4 has a charge of -1 (E.g. Bead "R1A4")
#             },
#         ],
#         "mapping_ratio": 1, # Not strictly needed here as the mapping ratio is assumed to be 1
#         ### Alternative way to set charges (residue nr, bead nr, charge)
#         # "charges": ((0, 1, 1), (1, 4, -1)), # Beads "R0A1" and "R1A4"
#     },
#})

