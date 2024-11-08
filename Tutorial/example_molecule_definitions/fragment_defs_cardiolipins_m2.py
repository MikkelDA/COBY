def tail_builder(tailcode, **kwargs):
    '''
    tailcode: [str]
        A string of letters indicating the tail bead combintation.

    kwargs["suffix"]: [str]
        Suffix added onto the end of tail bead names.

    kwargs["code_translator"]: [dict]
        Dictionary of accepted tailcode to bead convertions.

    kwargs["moltype"]: [str]
        The current lipid type.

    kwargs["molpart"]: [str]
        The current lipid part.

    return: [list of tuples]
        Returns a list of tuples formatted as (beadname, resname, resnr, xcoord, ycoord, zcoord, charge).
    '''
    ### Checks errormessage-specific input
    if "moltype" in kwargs:
        moltype = kwargs["moltype"]
    else:
        moltype = "False (not specified)"
    if "molpart" in kwargs:
        molpart = kwargs["molpart"]
    else:
        molpart = "False (not specified)"
    
    ### Checks if need kwargs are defined
    assert "suffix" in kwargs.keys(), "The key 'suffix' must be designated in the kwargs for lipid part '{molpart}' for lipid type '{moltype}'. It is allowed to be an empty string (e.g. '') in case no suffix is wanted.".format(molpart=molpart, moltype=moltype)
    assert "code_translator" in kwargs.keys(), "The key 'code_translator' must be designated in the kwargs for lipid part '{molpart}' for lipid type '{moltype}'.".format(molpart=molpart, moltype=moltype)
    assert len(kwargs["code_translator"]) > 0, "The key 'code_translator', for lipid part '{molpart}' for lipid type '{moltype}', must contain at least one 'code:bead' pair.".format(molpart=molpart, moltype=moltype)

    suffix = kwargs["suffix"]
    translator = kwargs["code_translator"]

    ### Checks if all letters in tailcode are in the translator
    assert all([i in translator.keys() for i in tailcode]), "\n".join([
        "At least one letter of the provided tailcode is not present in the code translator for lipid part '{molpart}' for lipid type '{moltype}'.".format(molpart=molpart, moltype=moltype),
        "Tailcode:",
        "    "+tailcode,
        "Code translator items (code:bead):",
        *["    " + key + ":" + val for key, val in translator.items()],
    ])

    ### List of dictionaries
    tail_beads = [
        {
            "name": translator[letter] + str(li+1) + str(suffix),
            "x": 0,
            "y": 0,
            "z": -li*0.3,
            "charge": 0,
            "resname": False,
            "resnr": 0,
        }
        for li, letter in enumerate(tailcode)
    ]

    return tail_beads

fragment_defs = {}

### General lipid type
moltype = "m2_cardiolipin"
fragment_defs = {}
fragment_defs[moltype] = {
    ### The parts that the "moltype" can accept
    "accepted_parts": ("core", "linker1", "tail11", "tail21", "linker2", "tail12", "tail22"),
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
            {"name": "GL0", "charge": 0,  "x": 0, "y": 0, "z": 0, "resname": False, "resnr": 0},
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
            {"name": "PO41", "charge": 0,  "x":  0,    "y": 0, "z":  0,   "resname": False, "resnr": 0},
            {"name": "GL11", "charge": 0,  "x":  0,    "y": 0, "z": -0.3, "resname": False, "resnr": 0},
            {"name": "GL21", "charge": 0,  "x": -0.25, "y": 0, "z": -0.3, "resname": False, "resnr": 0},
        ],
        "join_to": ("core", (0.15, 0, 1/3)), # GL1 bead
        "join_from": {
            "tail11": ( 0,    0, -1/3), # GL11 bead
            "tail21": (-0.25, 0, -1/3), # GL21 bead
        },
    },
})

molpart ="linker2"
fragment_defs[moltype]["parts"][molpart] = {}
fragment_defs[moltype]["parts"][molpart].update({
    "PhGL": { # Phosphoglycerol
        "beads": [
            {"name": "PO42", "charge": 0,  "x": 0,    "y": 0, "z":  0,   "resname": False, "resnr": 0},
            {"name": "GL12", "charge": 0,  "x": 0,    "y": 0, "z": -0.3, "resname": False, "resnr": 0},
            {"name": "GL22", "charge": 0,  "x": 0.25, "y": 0, "z": -0.3, "resname": False, "resnr": 0},
        ],
        "join_to": ("core", (-0.15, 0, 1/3)), # GL1 bead
        "join_from": {
            "tail12": (0,    0, -1/3), # GL12 bead
            "tail22": (0.25, 0, -1/3), # GL22 bead
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
            "C": "C", # Regular bead
            "c": "C", # Small bead
            "D": "D", # Regular bead with double bond
        },
    },
    "join_to": ("linker1", (0, 0, 1/3)), # GL11 bead
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
            "C": "C", # Regular bead
            "c": "C", # Small bead
            "D": "D", # Regular bead with double bond
        },
    },
    "join_to": ("linker1", (0, 0, 1/3)), # GL21 bead
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
            "C": "C", # Regular bead
            "c": "C", # Small bead
            "D": "D", # Regular bead with double bond
        },
    },
    "join_to": ("linker2", (0, 0, 1/3)), # GL12 bead
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
            "C": "C", # Regular bead
            "c": "C", # Small bead
            "D": "D", # Regular bead with double bond
        },
    },
    "join_to": ("linker2", (0, 0, 1/3)), # GL22 bead
})


