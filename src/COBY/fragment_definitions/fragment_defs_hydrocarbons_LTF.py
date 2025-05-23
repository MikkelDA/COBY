from COBY.fragment_definitions.function_tailbuilder import tail_builder

fragment_defs = {}
fragment_metadata = {}

moltype = "hydrocarbon_LTF"
fragment_metadata[moltype] = {
    "Description": "Fragments for this moleculetype (moltype) are based on the hydrocarbons from the lipid task force parameters for Martini 3 as described in https://doi.org/10.26434/chemrxiv-2024-8bjrr.",
}

### General lipid type
moltype = "hydrocarbon_LTF"
fragment_defs = {}
fragment_defs[moltype] = {
    ### The parts that the "moltype" can accept
    "accepted_parts":("tail",),
    ### The order in which parts should be written (to fit with the order in .itp files)
    "order": ("tail",),
    ### Empty dictionary for parts
    "parts": {},
    ### The default parts used in case the part is not specified
    "default_parts": {},
    ### Ascii art for reference
    ### Prefixing 'r' before a string forces raw formatting which prevents backslashes from causing problems in ascii art
    "ascii": [
        r"tail",
    ],
}

molpart ="tail"
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
})
