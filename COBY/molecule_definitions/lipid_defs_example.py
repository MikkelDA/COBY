lipid_defs = {}
params = "default"
lipid_defs[params] = {}
### Example of a multi-residue lipid molecule
### Two different ways to designate charges are shown
### ### Can be done within each residue dictionary or for the whole molecule
### ### Note that python-numbering starts from 0, meaning that 0 is the first bead (and 1 is the second)
#lipid_defs[params].update({
#     "Example_molecule": { # Name of molecule (used when selecting molecule in solvation/flooding commands)
#         "residues": [ # The list of residue dictionaries
#             { # Residue 0 (python indexing starts from zero)
#                 "resname": "RES0",
#                 "beads":   ("R0A0", "R0A1", "R0A2", "R0A3", "R0A4"),
#                 "x":       (     0,  0.125,   0.25,  0.125,      0),
#                 "y":       (     0,      0,      0,      0,      0),
#                 "z":       (     0,   0.25,    0.5,   0.75,      1),
#                 "charges": (1, 1), # Bead nr 1 has a charge of 1 (E.g. Bead "R0A1")
#             },
#             { # Residue 1
#                 "resname": "RES1",
#                 "beads":   ("R1A0", "R1A1", "R1A2", "R1A3", "R1A4"),
#                 "x":       (     0,  0.125,   0.25,  0.125,      0),
#                 "y":       (  0.25,   0.25,   0.25,   0.25,   0.25),
#                 "z":       (     0,   0.25,    0.5,   0.75,      1),
#                 "charges": (4, -1), # Bead nr 4 has a charge of -1 (E.g. Bead "R1A4")
#             },
#         ],
#         "mapping_ratio": 1, # Not strictly needed here as the mapping ratio is assumed to be 1
#         ### Alternative way to set charges (residue nr, bead nr, charge)
#         # "charges": ((0, 1, 1), (1, 4, -1)), # Beads "R0A1" and "R1A4"
#     },
#})

