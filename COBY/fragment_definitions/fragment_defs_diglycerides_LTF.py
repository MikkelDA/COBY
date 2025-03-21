
from COBY.fragment_definitions.function_tailbuilder import tail_builder

fragment_defs = {}

### General lipid type
moltype = "diglyceride_LTF"
fragment_defs = {}
fragment_defs[moltype] = {
    ### The parts that the "moltype" can accept
    "accepted_parts":("linker", "tail1", "tail2"),
    ### The order in which parts should be written (to fit with the order in .itp files)
    "order": ("linker", "tail1", "tail2"),
    ### Empty dictionary for parts
    "parts": {},
    ### The default parts used in case the part is not specified
    "default_parts": {
        "linker": "GL", # Currently the only option, so no point in requiring it to be specified
    },
    ### Ascii art for reference
    ### Prefixing 'r' before a string forces raw formatting which prevents backslashes from causing problems in ascii art
    "ascii": [
        r"  linker    ",
        r"  |     \   ", # Double \ because of string formatting
        r"tail1  tail2",
    ],
}

molpart ="linker"
fragment_defs[moltype]["parts"][molpart] = {}
fragment_defs[moltype]["parts"][molpart].update({
    "GL": {
        "beads": [
            {"name": "COH", "charge": 0, "x": 0,    "y": 0, "z": 0.3, "resname": False, "resnr": 0},
            {"name": "GL1", "charge": 0, "x": 0,    "y": 0, "z": 0,   "resname": False, "resnr": 0},
            {"name": "GL2", "charge": 0, "x": 0.25, "y": 0, "z": 0,   "resname": False, "resnr": 0},
        ],
        "join_from": {
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
            "C": "C", # Regular bead - Chain with no double bonds
            "c": "C", # Small bead   - Chain with no double bonds

            "D": "D", # Regular bead - Chain with 1 cis double bond

            "T": "T", # Regular bead - Chain with 1 trans double bond
            "t": "t", # Small bead   - Chain with 1 trans double bond

            "F": "F", # Regular bead - Chain with more than 1 double bond
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
            "C": "C", # Regular bead - Chain with no double bonds
            "c": "C", # Small bead   - Chain with no double bonds

            "D": "D", # Regular bead - Chain with 1 cis double bond

            "T": "T", # Regular bead - Chain with 1 trans double bond
            "t": "t", # Small bead   - Chain with 1 trans double bond

            "F": "F", # Regular bead - Chain with more than 1 double bond
        },
    },
    "join_to": ("linker", (-0.125, 0, 0.3)), # GL1 bead
})


