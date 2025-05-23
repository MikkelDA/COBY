from COBY.fragment_definitions.function_tailbuilder import tail_builder

fragment_defs = {}
fragment_metadata = {}

moltype = "BMP_LTF"
fragment_metadata[moltype] = {
    "Description": "Fragments for this moleculetype (moltype) are based on the BMPs (bis(monoacylglycero)phosphates) from the lipid task force parameters for Martini 3 as described in https://doi.org/10.26434/chemrxiv-2024-8bjrr.",
}

### General lipid type
moltype = "BMP_LTF" # bis(monoacylglycero)phosphates
fragment_defs = {}
fragment_defs[moltype] = {
    ### The parts that the "moltype" can accept
    "accepted_parts":("head", "linker1", "tail1", "linker2", "tail2"),
    ### The order in which parts should be written (to fit with the order in .itp files)
    "order": ("head", "linker1", "tail1", "linker2", "tail2"),
    ### Empty dictionary for parts
    "parts": {},
    ### The default parts used in case the part is not specified
    "default_parts": {
        "head": "PO4", # Currently the only option, so no point in requiring it to be specified
    },
    ### Ascii art for reference
    ### Prefixing 'r' before a string forces raw formatting which prevents backslashes from causing problems in ascii art
    "ascii": [
        r"       head       ",
        r"      /    \      ",
        r"linker1    linker2",
        r"   |          |   ",
        r" tail1      tail2 ",
    ],
}

molpart = "head"
fragment_defs[moltype]["parts"] = {}
fragment_defs[moltype]["parts"][molpart] = {}
fragment_defs[moltype]["parts"][molpart].update({
    "PO4": { # Phosphatidylcholine
        "beads": [
            {"name": "PO4", "charge": -1, "x": 0, "y": 0, "z": 0, "resname": False, "resnr": 0},
        ],
        "join_from": {
            "linker1": (0, 0, 0), # PO4 bead
            "linker2": (0, 0, 0), # PO4 bead
        },
    },
})

molpart ="linker1"
fragment_defs[moltype]["parts"][molpart] = {}
fragment_defs[moltype]["parts"][molpart].update({
    "2": {
        "beads": [
            {"name": "OH1", "charge": -1, "x": 0, "y": 0, "z": 0,    "resname": False, "resnr": 0},
            {"name": "GL1", "charge": -1, "x": 0, "y": 0, "z": -0.3, "resname": False, "resnr": 0},
        ],
        "join_to": ("head", (-0.2, 0, 1/4)), # PO4 bead
        "join_from": {
            "tail1": (0, 0, -0.3), # GL1 bead
        },
    },
    "3": {
        "beads": [
            {"name": "GL1", "charge": -1, "x": 0,    "y": 0, "z": 0, "resname": False, "resnr": 0},
            {"name": "OH1", "charge": -1, "x": 0.25, "y": 0, "z": 0, "resname": False, "resnr": 0},
        ],
        "join_to": ("head", (-0.2, 0, 1/4)), # PO4 bead
        "join_from": {
            "tail1": (0, 0, 0), # GL1 bead
        },
    },
})

molpart ="linker2"
fragment_defs[moltype]["parts"][molpart] = {}
fragment_defs[moltype]["parts"][molpart].update({
    "2": {
        "beads": [
            {"name": "OH2", "charge": -1, "x": 0, "y": 0, "z": 0,    "resname": False, "resnr": 0},
            {"name": "GL2", "charge": -1, "x": 0, "y": 0, "z": -0.3, "resname": False, "resnr": 0},
        ],
        "join_to": ("head", (0.2, 0, 1/4)), # PO4 bead
        "join_from": {
            "tail2": (0, 0, -0.3), # GL2 bead
        },
    },
    "3": {
        "beads": [
            {"name": "GL2", "charge": -1, "x": 0,    "y": 0, "z": 0, "resname": False, "resnr": 0},
            {"name": "OH2", "charge": -1, "x": -0.25, "y": 0, "z": 0, "resname": False, "resnr": 0},
        ],
        "join_to": ("head", (0.2, 0, 1/4)), # PO4 bead
        "join_from": {
            "tail2": (0, 0, 0), # GL2 bead
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
            "C": "C", # Regular bead - Chain with no double bonds
            "c": "C", # Small bead   - Chain with no double bonds

            "D": "D", # Regular bead - Chain with 1 cis double bond

            "T": "T", # Regular bead - Chain with 1 trans double bond
            "t": "t", # Small bead   - Chain with 1 trans double bond

            "F": "F", # Regular bead - Chain with more than 1 double bond
        },
    },
    "join_to": ("linker1", (0, 0, 0.3)), # GL1 bead
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
            "C": "C", # Regular bead - Chain with no double bonds
            "c": "C", # Small bead   - Chain with no double bonds

            "D": "D", # Regular bead - Chain with 1 cis double bond

            "T": "T", # Regular bead - Chain with 1 trans double bond
            "t": "t", # Small bead   - Chain with 1 trans double bond

            "F": "F", # Regular bead - Chain with more than 1 double bond
        },
    },
    "join_to": ("linker2", (0, 0, 0.3)), # GL1 bead
})


