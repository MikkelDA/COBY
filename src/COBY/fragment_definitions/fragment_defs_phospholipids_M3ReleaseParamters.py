from COBY.fragment_definitions.function_tailbuilder import tail_builder

fragment_defs = {}
fragment_metadata = {}

moltype = "phospholipid"
fragment_metadata[moltype] = {
    "Description": "Fragments for this moleculetype (moltype) are based on the phospholipids from the release parameters for Martini 3 phospholipids as described in https://doi.org/10.1038/s41592-021-01098-3.",
}

### General lipid type
### Includes glycerophoshpolipids (GL), etherphospholipids (ET) and plasmalogens (PL), depending on used linker
moltype = "phospholipid"
fragment_defs = {}
fragment_defs[moltype] = {
    ### The parts that the "moltype" can accept
    "accepted_parts":("head", "linker", "tail1", "tail2"),
    ### The order in which parts should be written (to fit with the order in .itp files)
    "order": ("head", "linker", "tail1", "tail2"),
    ### Empty dictionary for parts
    "parts": {},
    ### The default parts used in case the part is not specified
    "default_parts": {
        "linker": "GL", # Assume glycerophospholipids if linker not specified
    },
    ### Ascii art for reference
    ### Prefixing 'r' before a string forces raw formatting which prevents backslashes from causing problems in ascii art
    "ascii": [
        r"  head      ",
        r"  |         ",
        r"  linker    ",
        r"  |     \   ",
        r"tail1  tail2",
    ],
}

molpart = "head"
fragment_defs[moltype]["parts"] = {}
fragment_defs[moltype]["parts"][molpart] = {}
fragment_defs[moltype]["parts"][molpart].update({
    "PC": { # Phosphatidylcholine
        "beads": [
            {"name": "NC3", "charge": 1,  "x": 0, "y": 0, "z": 0.3, "resname": False, "resnr": 0},
            {"name": "PO4", "charge": -1, "x": 0, "y": 0, "z": 0,   "resname": False, "resnr": 0},
        ],
        "join_to": ("linker", (0, 0, -0.3)), # GL1 bead
    },
    
    "PE": { # Phosphatidylethanolamine
        "beads": [
            {"name": "NH3", "charge": 1,  "x": 0, "y": 0, "z": 0.3, "resname": False, "resnr": 0},
            {"name": "PO4", "charge": -1, "x": 0, "y": 0, "z": 0,   "resname": False, "resnr": 0},
        ],
        "join_to": ("linker", (0, 0, -0.3)), # GL1 bead
    },
    
    "PG": { # Phosphatidylglycerol
        "beads": [
            {"name": "GL0", "charge": 0,  "x": 0, "y": 0, "z": 0.3, "resname": False, "resnr": 0},
            {"name": "PO4", "charge": -1, "x": 0, "y": 0, "z": 0,   "resname": False, "resnr": 0},
        ],
        "join_to": ("linker", (0, 0, -0.3)), # GL1 bead
    },
    
    "PA": { # Phosphatidic acid
        "beads": [
            {"name": "PO4", "charge": -1, "x": 0, "y": 0, "z": 0,   "resname": False, "resnr": 0},
        ],
        "join_to": ("linker", (0, 0, -0.3)), # GL1 bead
    },
    
    "PS": { # Phosphatidylserine
        "beads": [
            {"name": "CNO", "charge": 0,  "x": 0, "y": 0, "z": 0.3, "resname": False, "resnr": 0},
            {"name": "PO4", "charge": -1, "x": 0, "y": 0, "z": 0,   "resname": False, "resnr": 0},
        ],
        "join_to": ("linker", (0, 0, -0.3)), # GL1 bead
    },
})    

molpart ="linker"
fragment_defs[moltype]["parts"][molpart] = {}
fragment_defs[moltype]["parts"][molpart].update({
    "GL": { # Glycerol
        "beads": [
            {"name": "GL1", "charge": 0, "x": 0,     "y": 0, "z": 0, "resname": False, "resnr": 0},
            {"name": "GL2", "charge": 0, "x": 0.125, "y": 0, "z": 0, "resname": False, "resnr": 0},
        ],
        "join_from": {
            "head":  (0,    0, 0), # GL1 bead
            "tail1": (0,    0, 0), # GL1 bead
            "tail2": (0.25, 0, 0), # GL2 bead
        },
    },
})

molpart ="tail1"
fragment_defs[moltype]["parts"][molpart] = {}
fragment_defs[moltype]["parts"][molpart].update({
    "function": tail_builder,
    "kwargs": {
        "moltype": moltype, # Only used for error messages
        "molpart": molpart, # Only used for error messages
        "suffix": "A",
        "code_translator": {
            "C": "C", # Regular bead
            "c": "C", # Small bead
            "D": "D", # Regular bead with double bond
        },
    },
    "join_to": ("linker", (0, 0, 0.3)), # GL1 bead
})

molpart ="tail2"
fragment_defs[moltype]["parts"][molpart] = {}
fragment_defs[moltype]["parts"][molpart].update({
    "function": tail_builder,
    "kwargs": {
        "moltype": moltype, # Only used for error messages
        "molpart": molpart, # Only used for error messages
        "suffix": "B",
        "code_translator": {
            "C": "C", # Regular bead
            "c": "C", # Small bead
            "D": "D", # Regular bead with double bond
        },
    },
    "join_to": ("linker", (0, 0, 0.3)), # GL2 bead
})
