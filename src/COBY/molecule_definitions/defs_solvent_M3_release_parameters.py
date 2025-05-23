solvent_defs = {}
solvent_metadata = {}

params = "default"
solvent_metadata[params] = {
    "Description": "Solvent and solutes in this parameter library are based on the release parameters for Martini 3 as described in https://doi.org/10.1038/s41592-021-01098-3.",
}

################
### SOLVENTS ###
################
### molar_mass: [g/mol]
### density: [g/cm^3]
params = "default"
solvent_defs[params] = {}

### Water
solvent_defs[params].update({
    "W" : {
        "mapping_ratio": 4, "density": 0.99669, "molar_mass": 18.01528,
        "beads": [{"name": "W", "x": 0, "y": 0, "z": 0, "charge": 0}],
        "tags": ("martini3", "water"),
    },
    "SW" : {
        "mapping_ratio": 3, "density": 0.99669, "molar_mass": 18.01528,
        "beads": [{"name": "SW", "x": 0, "y": 0, "z": 0, "charge": 0}],
        "tags": ("martini3", "water"),
    },
    "TW" : {
        "mapping_ratio": 2, "density": 0.99669, "molar_mass": 18.01528,
        "beads": [{"name": "TW", "x": 0, "y": 0, "z": 0, "charge": 0}],
        "tags": ("martini3", "water"),
    },
})

### Amino acids
solvent_defs[params].update({
    ### ### 1 atom residues
    "GLY" : {
        "mapping_ratio": 1,
        "beads": [
            {"name": "BB", "x": 0, "y": 0, "z": 0, "charge": 0},
        ],
        "tags": ("martini3", "amino_acid"),
    },
    
    ### ### 2 atom residues
    ### Residues without variants
    **{
        resname: {
            "mapping_ratio": 1,
            "beads": [
                {"name": "BB",  "x":  0.25, "y": 0, "z": 0, "charge": 0},
                {"name": "SC1", "x": -0.25, "y": 0, "z": 0, "charge": 0},
            ],
            "tags": ("martini3", "amino_acid"),
        } for resname in [
            ### Non-charged 2-bead residues without variants
            "ALA", "CYS", "VAL", "LEU", "ILE", "MET", "PRO", "HYP", "ASN", "GLN", "THR", "SER",
            ### Neutral ASP variants
            "ASPP", "ASH",
            ### Neutral GLU variants
            "GLUP", "GLH",
        ]
    },

    ### Charged ASP and GLU
    **{
        resname: {
            "mapping_ratio": 1,
            "beads": [
                {"name": "BB",  "x":  0.25, "y": 0, "z": 0, "charge": 0},
                {"name": "SC1", "x": -0.25, "y": 0, "z": 0, "charge": -1},
            ],
            "tags": ("martini3", "amino_acid"),
        } for resname in ["ASP", "GLU"]
    },

    ### ### 3 atom residues
    ### Charged ARG and LYS
    **{
        resname: {
            "mapping_ratio": 1,
            "beads": [
                {"name": "BB",  "x":  0.25, "y": 0,     "z": 0, "charge": 0},
                {"name": "SC1", "x":  0,    "y": 0,     "z": 0, "charge": 0},
                {"name": "SC2", "x": -0.25, "y": 0.125, "z": 0, "charge": 1},
            ],
            "tags": ("martini3", "amino_acid"),
        } for resname in ["ARG", "LYS"]
    },
    ### Non-charged variants of lysine
    **{
        resname: {
            "mapping_ratio": 1,
            "beads": [
                {"name": "BB",  "x":  0.25, "y": 0,     "z": 0, "charge": 0},
                {"name": "SC1", "x":  0,    "y": 0,     "z": 0, "charge": 0},
                {"name": "SC2", "x": -0.25, "y": 0.125, "z": 0, "charge": 0},
            ],
            "tags": ("martini3", "amino_acid"),
        } for resname in ["LSN", "LYN"]
    },

    ### ### 4 atom residues
    ### Non-charged PHE and HIS variants
    **{
        resname: {
            "mapping_ratio": 1,
            "beads": [
                {"name": "BB",  "x":  0.25, "y": 0,      "z": 0, "charge": 0},
                {"name": "SC1", "x":  0,    "y": 0,      "z": 0, "charge": 0},
                {"name": "SC2", "x": -0.25, "y": 0.125,  "z": 0, "charge": 0},
                {"name": "SC3", "x": -0.25, "y": -0.125, "z": 0, "charge": 0},
            ],
            "tags": ("martini3", "amino_acid"),
        } for resname in ["PHE", "HIS", "HIE", "HSE", "HSD", "HID"]
    },
    ### Charged variants of histidine
    **{
        resname: {
            "mapping_ratio": 1,
            "beads": [
                {"name": "BB",  "x":  0.25, "y": 0,      "z": 0, "charge": 0},
                {"name": "SC1", "x":  0,    "y": 0,      "z": 0, "charge": 0},
                {"name": "SC2", "x": -0.25, "y": 0.125,  "z": 0, "charge": 0.5},
                {"name": "SC3", "x": -0.25, "y": -0.125, "z": 0, "charge": 0.5},
            ],
            "tags": ("martini3", "amino_acid"),
        } for resname in ["HSP", "HIP"]
    },
    
    ### ### 5 atom residues
    "TYR" : {
        "mapping_ratio": 1,
        "beads": [
            {"name": "BB",  "x":  0.25, "y":  0.125, "z": 0, "charge": 0},
            {"name": "SC1", "x":  0.25, "y":  0,     "z": 0, "charge": 0},
            {"name": "SC2", "x":  0,    "y": -0.125, "z": 0, "charge": 0},
            {"name": "SC3", "x":  0,    "y":  0.125, "z": 0, "charge": 0},
            {"name": "SC4", "x": -0.25, "y":  0,     "z": 0, "charge": 0},
        ],
        "tags": ("martini3", "amino_acid"),
    },
    
    ### ### 6 atom residues
    "TRP" : {
        "mapping_ratio": 1,
        "beads": [
            {"name": "BB",  "x":  0.25, "y":  0.125, "z": 0, "charge": 0},
            {"name": "SC1", "x":  0.25, "y":  0,     "z": 0, "charge": 0},
            {"name": "SC2", "x":  0,    "y": -0.125, "z": 0, "charge": 0},
            {"name": "SC3", "x":  0,    "y":  0.125, "z": 0, "charge": 0},
            {"name": "SC4", "x": -0.25, "y":  0,     "z": 0, "charge": 0},
            {"name": "SC5", "x": -0.25, "y":  0.125, "z": 0, "charge": 0},
        ],
        "tags": ("martini3", "amino_acid"),
    },
})


### Below is shown an example of a multi-resdiue solvent
### ### Note that it is not necessary to include "charge" in all bead dictionaries as it is assumed to be zero if not defined
# solvent_defs[params].update({
#     "Example_solvent": { # Name of molecule (used when selecting molecule in solvation/flooding commands)
#         "residues": [ # The list of residue dictionaries
#             { # Residue 0 (python indexing starts from zero)
#                 "resname": "RES0",
#                 "beads": [
#                     {"name": "R0A0", "x": 0,     "y":  0,     "z": 0,                },
#                     {"name": "R0A1", "x": 0.125, "y":  0,     "z": 0.125, "charge": 1},
#                     {"name": "R0A2", "x": 0.25,  "y":  0.125, "z": 0,                },
#                     {"name": "R0A3", "x": 0.375, "y":  0,     "z": 0.125,            },
#                     {"name": "R0A4", "x": 0.5,   "y": -0.125, "z": 0,                },
#                 ]
#             },
#             { # Residue 1
#                 "resname": "RES1",
#                 "beads": [
#                     {"name": "R1A0", "x": 0,     "y":  0,     "z": 0,                 },
#                     {"name": "R1A1", "x": 0.125, "y":  0,     "z": 0.125,             },
#                     {"name": "R1A2", "x": 0.25,  "y":  0.125, "z": 0,                 },
#                     {"name": "R1A3", "x": 0.375, "y":  0,     "z": 0.125,             },
#                     {"name": "R1A4", "x": 0.5,   "y": -0.125, "z": 0,     "charge": -1},
#                 ]
#             },
#         ],
#         "mapping_ratio": 1, # Not strictly needed here as the mapping ratio is assumed to be 1, but shown for clarity
#         "tags": ("martini3", "custom_solvent"), # Tags used with COBY.Library() method for examining the molecule parameter libraries
#     },
# })
