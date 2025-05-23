from COBY.fragment_definitions.function_tailbuilder import tail_builder

fragment_defs = {}
fragment_metadata = {}

moltype = "triglyceride_LTF"
fragment_metadata[moltype] = {
    "Description": "Fragments for this moleculetype (moltype) are based on the monoglycerides from the lipid task force parameters for Martini 3 as described in https://doi.org/10.26434/chemrxiv-2024-8bjrr.",
}

### General lipid type
moltype = "triglyceride_LTF"
fragment_defs = {}
fragment_defs[moltype] = {
    ### The parts that the "moltype" can accept
    "accepted_parts":("linker1", "tail1", "linker2", "tail2", "linker3", "tail3"),
    ### The order in which parts should be written (to fit with the order in .itp files)
    "order": ("linker1", "tail1", "linker2", "tail2", "linker3", "tail3"),
    ### Empty dictionary for parts
    "parts": {},
    ### The default parts used in case the part is not specified
    "default_parts": {
        "linker1": "GL", # Currently the only option, so no point in requiring it to be specified
        "linker2": "GL", # Currently the only option, so no point in requiring it to be specified
        "linker3": "GL", # Currently the only option, so no point in requiring it to be specified
    },
    ### Ascii art for reference
    ### Prefixing 'r' before a string forces raw formatting which prevents backslashes from causing problems in ascii art
    "ascii": [
        r"linker1-linker2-linker3",
        r"  |       |       |    ",
        r"tail1   tail2   tail3  ",
    ],
}

molpart ="linker1"
fragment_defs[moltype]["parts"][molpart] = {}
fragment_defs[moltype]["parts"][molpart].update({
    "GL": {
        "beads": [
            {"name": "GL1", "charge": 0, "x": 0, "y": 0, "z": 0, "resname": False, "resnr": 0},
        ],
        "join_from": {
            "tail1":   (0, 0, 0), # GL1 bead
            "linker2": (0, 0, 0), # GL1 bead
            "linker3": (0, 0, 0), # GL1 bead
        },
    },
})

molpart ="linker2"
fragment_defs[moltype]["parts"][molpart] = {}
fragment_defs[moltype]["parts"][molpart].update({
    "GL": {
        "beads": [
            {"name": "GL2", "charge": 0, "x": 0, "y": 0, "z": 0, "resname": False, "resnr": 0},
        ],
        "join_to": ("linker1", (-0.25, 0, 0)), # GL1 bead
        "join_from": {
            "tail2":   (0, 0, 0), # GL2 bead
            "linker3": (0, 0, 0), # GL2 bead
        },
    },
})

molpart ="linker3"
fragment_defs[moltype]["parts"][molpart] = {}
fragment_defs[moltype]["parts"][molpart].update({
    "GL": {
        "beads": [
            {"name": "GL3", "charge": 0, "x": 0, "y": 0, "z": 0, "resname": False, "resnr": 0},
        ],
        "join_to": ("linker2", (-0.25, 0, 0)), # GL1 bead
        "join_from": {
            "tail3": (0, 0, 0), # GL2 bead
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
    "join_to": ("linker2", (0, 0, 0.3)), # GL2 bead
})
molpart ="tail3"
fragment_defs[moltype]["parts"][molpart] = {}
fragment_defs[moltype]["parts"][molpart].update({
    "function": tail_builder,
    "kwargs": {
        "moltype": moltype, # Only used for error messages
        "molpart": molpart, # Only used for error messages
        "suffix": "C",
        "code_translator": {
            "C": "C", # Regular bead - Chain with no double bonds
            "c": "C", # Small bead   - Chain with no double bonds

            "D": "D", # Regular bead - Chain with 1 cis double bond

            "T": "T", # Regular bead - Chain with 1 trans double bond
            "t": "t", # Small bead   - Chain with 1 trans double bond

            "F": "F", # Regular bead - Chain with more than 1 double bond
        },
    },
    "join_to": ("linker3", (0, 0, 0.3)), # GL3 bead
})


