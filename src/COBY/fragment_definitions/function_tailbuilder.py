def tail_builder(tailcode, **kwargs):
    '''
    tailcode: [str]
        A string of letters indicating the tail bead combintation.

    kwargs["suffix"]: [str]
        Suffix added onto the end of tail bead names.

    kwargs["code_translator"]: [dict]
        Dictionary of accepted tailcode to bead convertions.

    kwargs["lipidtype"]: [str]
        The current lipid type.

    kwargs["lipidpart"]: [str]
        The current lipid part.

    return: [list of tuples]
        Returns a list of tuples formatted as (beadname, resname, resnr, xcoord, ycoord, zcoord, charge).
    '''
    ### Checks errormessage-specific input
    if "lipidtype" in kwargs:
        lipidtype = kwargs["lipidtype"]
    else:
        lipidtype = "False (not specified)"
    if "lipidpart" in kwargs:
        lipidpart = kwargs["lipidpart"]
    else:
        lipidpart = "False (not specified)"
    
    ### Checks if need kwargs are defined
    assert "suffix" in kwargs.keys(), "The key 'suffix' must be designated in the kwargs for lipid part '{lipidpart}' for lipid type '{lipidtype}'. It is allowed to be an empty string (e.g. '') in case no suffix is wanted.".format(lipidpart=lipidpart, lipidtype=lipidtype)
    assert "code_translator" in kwargs.keys(), "The key 'code_translator' must be designated in the kwargs for lipid part '{lipidpart}' for lipid type '{lipidtype}'.".format(lipidpart=lipidpart, lipidtype=lipidtype)
    assert len(kwargs["code_translator"]) > 0, "The key 'code_translator', for lipid part '{lipidpart}' for lipid type '{lipidtype}', must contain at least one 'code:bead' pair.".format(lipidpart=lipidpart, lipidtype=lipidtype)

    suffix = kwargs["suffix"]
    translator = kwargs["code_translator"]

    ### Checks if all letters in tailcode are in the translator
    assert all([i in translator.keys() for i in tailcode]), "\n".join([
        "At least one letter of the provided tailcode is not present in the code translator for lipid part '{lipidpart}' for lipid type '{lipidtype}'.".format(lipidpart=lipidpart, lipidtype=lipidtype),
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
