from COBY.fragment_definitions.function_tailbuilder import tail_builder

fragment_defs = {}
fragment_metadata = {}

moltype = "cardiolipin_LTF"
fragment_metadata[moltype] = {
    "Description": "Fragments for this moleculetype (moltype) are based on the cardiolipins from the lipid task force parameters for Martini 3 as described in https://doi.org/10.26434/chemrxiv-2024-8bjrr.",
}

### General lipid type
moltype = "cardiolipin_LTF"
fragment_defs = {}
fragment_defs[moltype] = {
    ### The parts that the "moltype" can accept
    "accepted_parts":("core", "linker1", "tail11", "tail21", "linker2", "tail12", "tail22"),
    ### The order in which parts should be written (to fit with the order in .itp files)
    "order": ("core", "linker1", "tail11", "tail21", "linker2", "tail12", "tail22"),
    ### Empty dictionary for parts
    "parts": {},
    ### The default parts used in case the part is not specified
    "default_parts": {
        "core": "GLC",     # Currently the only option, so no point in requiring it to be specified
        "linker1": "PhGL", # Currently the only option, so no point in requiring it to be specified
        "linker2": "PhGL", # Currently the only option, so no point in requiring it to be specified
    },
    ### Ascii art for reference
    ### Prefixing 'r' before a string forces raw formatting which prevents backslashes from causing problems in ascii art
    "ascii": [
        r"            core            ",
        r"           /    \           ",
        r"     linker1    linker2     ",
        r"    /   |        |     \    ",
        r"tail21 tail11  tail12 tail22",
    ],
}

# """
# core = GLC
# linker1 = PO41, GL11, GL21
# linker2 = PO42, GL12, GL22
#             GLC
#            /   \
#         PO41   PO42
#           |     |
#    GL11-GL12   GL12-GL21
#      |    |     |    |
#    C1B1 C1A1   C1A2 C1B2
#      |    |     |    |
#      *    *     *    *
#      |    |     |    |
#      *    *     *    *
# """

molpart ="core"
fragment_defs[moltype]["parts"][molpart] = {}
fragment_defs[moltype]["parts"][molpart].update({
    "GLC": { # Glycerol Core
        "beads": [
            {"name": "GLC", "charge": 0, "x": 0, "y": 0, "z": 0, "resname": False, "resnr": 0},
        ],
        "join_from": {
            "linker1": (0, 0, 0), # PO41 bead
            "linker2": (0, 0, 0), # PO42 bead
        },
    },
})

molpart ="linker1"
fragment_defs[moltype]["parts"][molpart] = {}
fragment_defs[moltype]["parts"][molpart].update({
    "PhGL": { # Phosphoglycerol
        "beads": [
            {"name": "PO41", "charge": 0, "x": 0,     "y": 0, "z": 0,    "resname": False, "resnr": 0},
            {"name": "GL11", "charge": 0, "x": 0,     "y": 0, "z": -0.3, "resname": False, "resnr": 0},
            {"name": "GL21", "charge": 0, "x": -0.25, "y": 0, "z": -0.3, "resname": False, "resnr": 0},
        ],
        "join_to": ("core", (0.15, 0, 0.3)), # GLC bead
        "join_from": {
            "tail11": (0,     0, -0.3), # GL11 bead
            "tail21": (-0.25, 0, -0.3), # GL21 bead
        },
    },
})

molpart ="linker2"
fragment_defs[moltype]["parts"][molpart] = {}
fragment_defs[moltype]["parts"][molpart].update({
    "PhGL": { # Phosphoglycerol
        "beads": [
            {"name": "PO42", "charge": 0, "x": 0,    "y": 0, "z": 0,    "resname": False, "resnr": 0},
            {"name": "GL12", "charge": 0, "x": 0,    "y": 0, "z": -0.3, "resname": False, "resnr": 0},
            {"name": "GL22", "charge": 0, "x": 0.25, "y": 0, "z": -0.3, "resname": False, "resnr": 0},
        ],
        "join_to": ("core", (-0.15, 0, 0.3)), # GLC bead
        "join_from": {
            "tail12": (0,    0, -0.3), # GL12 bead
            "tail22": (0.25, 0, -0.3), # GL22 bead
        },
    },
})

molpart ="tail11"
fragment_defs[moltype]["parts"][molpart] = {}
fragment_defs[moltype]["parts"][molpart].update({
    "function": tail_builder,
    "kwargs": {
        "moltype": moltype, # Only used for error messages
        "molpart": molpart, # Only used for error messages
        "suffix": "A1",
        "code_translator": {
            "C": "C", # Regular bead - Chain with no double bonds
            "c": "C", # Small bead   - Chain with no double bonds

            "D": "D", # Regular bead - Chain with 1 cis double bond

            "T": "T", # Regular bead - Chain with 1 trans double bond
            "t": "t", # Small bead   - Chain with 1 trans double bond

            "F": "F", # Regular bead - Chain with more than 1 double bond
        },
    },
    "join_to": ("linker1", (0, 0, 0.3)), # GL11 bead
})
molpart ="tail21"
fragment_defs[moltype]["parts"][molpart] = {}
fragment_defs[moltype]["parts"][molpart].update({
    "function": tail_builder,
    "kwargs": {
        "moltype": moltype, # Only used for error messages
        "molpart": molpart, # Only used for error messages
        "suffix": "B1",
        "code_translator": {
            "C": "C", # Regular bead - Chain with no double bonds
            "c": "C", # Small bead   - Chain with no double bonds

            "D": "D", # Regular bead - Chain with 1 cis double bond

            "T": "T", # Regular bead - Chain with 1 trans double bond
            "t": "t", # Small bead   - Chain with 1 trans double bond

            "F": "F", # Regular bead - Chain with more than 1 double bond
        },
    },
    "join_to": ("linker1", (0, 0, 0.3)), # GL21 bead
})

molpart ="tail12"
fragment_defs[moltype]["parts"][molpart] = {}
fragment_defs[moltype]["parts"][molpart].update({
    "function": tail_builder,
    "kwargs": {
        "moltype": moltype, # Only used for error messages
        "molpart": molpart, # Only used for error messages
        "suffix": "A2",
        "code_translator": {
            "C": "C", # Regular bead - Chain with no double bonds
            "c": "C", # Small bead   - Chain with no double bonds

            "D": "D", # Regular bead - Chain with 1 cis double bond

            "T": "T", # Regular bead - Chain with 1 trans double bond
            "t": "t", # Small bead   - Chain with 1 trans double bond

            "F": "F", # Regular bead - Chain with more than 1 double bond
        },
    },
    "join_to": ("linker2", (0, 0, 0.3)), # GL12 bead
})
molpart ="tail22"
fragment_defs[moltype]["parts"][molpart] = {}
fragment_defs[moltype]["parts"][molpart].update({
    "function": tail_builder,
    "kwargs": {
        "moltype": moltype, # Only used for error messages
        "molpart": molpart, # Only used for error messages
        "suffix": "B2",
        "code_translator": {
            "C": "C", # Regular bead - Chain with no double bonds
            "c": "C", # Small bead   - Chain with no double bonds

            "D": "D", # Regular bead - Chain with 1 cis double bond

            "T": "T", # Regular bead - Chain with 1 trans double bond
            "t": "t", # Small bead   - Chain with 1 trans double bond

            "F": "F", # Regular bead - Chain with more than 1 double bond
        },
    },
    "join_to": ("linker2", (0, 0, 0.3)), # GL22 bead
})


