
### Below is shown an example of a multi-resdiue lipid
### ### Note that it is not necessary to include "charge" in all bead dictionaries as it is assumed to be zero if not defined
# lipid_defs = {}
# params = "default"
# lipid_defs[params] = {}
# lipid_defs[params].update({
#     "Example_lipid": { # Name of molecule (used when selecting molecule in solvation/flooding commands)
#         "residues": [ # The list of residue dictionaries
#             { # Residue 0 (python indexing starts from zero)
#                 "resname": "RES0",
#                 "beads": [
#                     {"name": "R0A0", "x": 0,     "y": 0, "z": 0,               },
#                     {"name": "R0A1", "x": 0.125, "y": 0, "z": 0.25, "charge": 1},
#                     {"name": "R0A2", "x": 0.25,  "y": 0, "z": 0.5,             },
#                     {"name": "R0A3", "x": 0.125, "y": 0, "z": 0.75,            },
#                     {"name": "R0A4", "x": 0,     "y": 0, "z": 1,               },
#                 ],
#             },
#             { # Residue 1
#                 "resname": "RES1",
#                 "beads": [
#                     {"name": "R1A0", "x": 0,     "y": 0.25, "z": 0,                },
#                     {"name": "R1A1", "x": 0.125, "y": 0.25, "z": 0.25,             },
#                     {"name": "R1A2", "x": 0.25,  "y": 0.25, "z": 0.5,              },
#                     {"name": "R1A3", "x": 0.125, "y": 0.25, "z": 0.75,             },
#                     {"name": "R1A4", "x": 0,     "y": 0.25, "z": 1,    "charge": -1},
#                 ],
#             },
#         ],
#         "mapping_ratio": 1, # Not strictly needed here as the mapping ratio is assumed to be 1, but shown for clarity
#         "tags": ("martini3", "custom_lipid"), # Tags used with COBY.Library() method for examining the molecule parameter libraries
#     },
# })
