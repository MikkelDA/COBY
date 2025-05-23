
from COBY.fragment_definitions.function_tailbuilder import tail_builder

import copy

fragment_defs = {}
fragment_metadata = {}

moltype = "monoglyceride_LTF"
fragment_metadata[moltype] = {
    "Description": "Fragments for this moleculetype (moltype) are based on the monoglycerides from the lipid task force parameters for Martini 3 as described in https://doi.org/10.26434/chemrxiv-2024-8bjrr.",
}

### General lipid type
### Includes glycerophoshpolipids (GL), etherphospholipids (ET) and plasmalogens (PL), depending on used linker
moltype = "sphingolipid_LTF"
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
        "linker": "SM", # Currently the only option, so no point in requiring it to be specified
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

### All heads simply copied from "fragment_defs_phospholipids.py"
molpart = "head"
fragment_defs[moltype]["parts"] = {}
fragment_defs[moltype]["parts"][molpart] = {}
fragment_defs[moltype]["parts"][molpart].update({
    "CER": { # Ceramide COH head
        "beads": [
            {"name": "COH", "charge": 0, "x": 0, "y": 0, "z": 0, "resname": False, "resnr": 0},
        ],
        "join_to": ("linker", (0, 0, -0.3)), # GL1 bead
    },
    
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
            {"name": "PO4", "charge": -1, "x": 0, "y": 0, "z": 0, "resname": False, "resnr": 0},
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

### Generating phosphoinositol (PI) head groups using scaffold-like method
### Beads (and coordinates further down) are ordered in accordance with the order that they appear in topology files.
inositol_head_beads = {
    "PI": ["C1", "C2", "C3", "C4", "PO4", "",   "",   "",  ],
    "P1": ["C1", "C2", "C3", "C4", "PO4", "P3", "",   "",  ],
    "P2": ["C1", "C2", "C3", "C4", "PO4", "P3", "P4", "",  ],
    "P3": ["C1", "C2", "C3", "C4", "PO4", "P3", "P4", "P5",],
    "P4": ["C1", "C2", "C3", "C4", "PO4", "",   "P4", "",  ],
    "P5": ["C1", "C2", "C3", "C4", "PO4", "",   "",   "P5",],
    "P6": ["C1", "C2", "C3", "C4", "PO4", "",   "P4", "P5",],
    "P7": ["C1", "C2", "C3", "C4", "PO4", "P3", "",   "P5",],
}

inositol_head_charges = {
    "PI": {"PO4": -1},
    "P1": {"PO4": -1, "P3": -2.0},
    "P2": {"PO4": -1, "P3": -1.5, "P4": -1.5},
    "P3": {"PO4": -1, "P3": -1.3, "P4": -1.4, "P5": -1.3},
    "P4": {"PO4": -1, "P4": -2.0},
    "P5": {"PO4": -1, "P5": -2.0},
    "P6": {"PO4": -1, "P4": -1.5, "P5": -1.5},
    "P7": {"PO4": -1, "P3": -1.5, "P5": -1.5},
}

### Phosphoinositol structure
#   P5     P4
#   |       \    P3
#  C2 — — — C3 /
#    \  C4  /
#     \    /
#       C1
#       |
#      PO4
#
# Axes: PO4 is origo to fit with linker join position
#      /|\+z
#       |
#  -x   |   x+
# <----PO4---->

#              C1      C2     C3      C4  PO4     P3   P4    P5
inostiol_x = (  0,   -0.2,   0.2,      0,   0, 0.275,   0.175,   -0.2)
inostiol_y = (  0,      0,     0,      0,   0,     0,       0,      0)
inostiol_z = (1/3,    2/3, 1.75/3, 1.5/3,   0, 2.25/3,  2.5/3, 2.75/3)

for headname, beadslist in inositol_head_beads.items():
    xs = [coord for i, coord in enumerate(inostiol_x) if beadslist[i]]
    ys = [coord for i, coord in enumerate(inostiol_y) if beadslist[i]]
    zs = [coord for i, coord in enumerate(inostiol_z) if beadslist[i]]
    beadslist_reduced = [bead for i, bead in enumerate(beadslist) if beadslist[i]]
    
    charges = [inositol_head_charges[headname][bead] if bead in inositol_head_charges[headname] else 0 for bead in beadslist_reduced]
    
    beads = [
        {"name": beadname, "charge": charge, "x": x, "y": y, "z": z, "resname": False, "resnr": 0}
        for beadname, charge, x, y, z in zip(beadslist_reduced, charges, xs, ys, zs)
    ]

    fragment_defs[moltype]["parts"][molpart].update({
        headname: {
            "beads": beads,
            "join_to": ("linker", (0, 0, -0.3)), # OH1 bead
        },
    })

molpart ="linker"
fragment_defs[moltype]["parts"][molpart] = {}
fragment_defs[moltype]["parts"][molpart].update({
    "SM": { # Sphingomyelin
        "beads": [
            {"name": "OH1", "charge": 0, "x": 0,    "y": 0, "z": 0, "resname": False, "resnr": 0},
            {"name": "AM2", "charge": 0, "x": 0.25, "y": 0, "z": 0, "resname": False, "resnr": 0},
        ],
        "join_from": {
            "head":  (0,    0, 0), # OH1 bead
            "tail1": (0,    0, 0), # OH1 bead
            "tail2": (0.25, 0, 0), # AM2 bead
        },
    },
})

molpart ="tail1" # Tail 1 is long chain base (LCB)
fragment_defs[moltype]["parts"][molpart] = {}
fragment_defs[moltype]["parts"][molpart].update({
    "function": tail_builder,
    "kwargs": {
        "moltype": moltype, # Only used for error messages
        "molpart": molpart, # Only used for error messages
        "suffix": "A",
        "code_translator": {
            "C": "C", # Regular bead - Chain with no double bonds
            "c": "C", # Small bead   - Chain with no double bonds

            "D": "D", # Regular bead - Chain with 1 cis double bond

            "T": "T", # Regular bead - Chain with 1 trans double bond
            "t": "t", # Small bead   - Chain with 1 trans double bond

            "F": "F", # Regular bead - Chain with more than 1 double bond
        },
    },
    "join_to": ("linker", (0, 0, 0.3)), # OH1 bead
})

molpart ="tail2" # Tail 2 is fatty acid (FA)
fragment_defs[moltype]["parts"][molpart] = {}
fragment_defs[moltype]["parts"][molpart].update({
    "function": tail_builder,
    "kwargs": {
        "moltype": moltype, # Only used for error messages
        "molpart": molpart, # Only used for error messages
        "suffix": "B",
        "code_translator": {
            "C": "C", # Regular bead - Chain with no double bonds
            "c": "C", # Small bead   - Chain with no double bonds

            "D": "D", # Regular bead - Chain with 1 cis double bond

            "T": "T", # Regular bead - Chain with 1 trans double bond
            "t": "t", # Small bead   - Chain with 1 trans double bond

            "F": "F", # Regular bead - Chain with more than 1 double bond
        },
    },
    "join_to": ("linker", (0, 0, 0.3)), # AM2 bead
})


