__version__="1.0.11"

major_changes = [
    "Residue numbers for structures imported using the 'protein' argument can now be conserved in the outputted structure.",
    [
        "System components are separated with 'TER' lines in outputted .pdb files.",
        "Outputted .gro files are unaffected as they do not have a similar functinality to 'TER'. This means that outputted .gro and .pdb files may not be identical in terms of atom/residue numbering.",
        "This new behaviour can be turned off using the 'keep_residue_numbering'/'krn' argument as shown below:",
        [
            "COBY(..., keep_residue_numbering = \"False\", ...)",
            "python -m COBY ... -keep_residue_numbering False ...)",
        ],
    ],
]

minor_changes = [
    "Structure file importers now check if all x/y/z coordinate values are numbers and not 'nan' and crash the program if a 'nan' is found",
    "Added warnings when using 'charge:lib' or 'charge:[val]' alongside 'moleculetypes'",
]

bug_fixes = [
    "Fixed crash in the 'Library' program when looking through tags.",
    "Fixed large proteins not having their residue numbers correctly assigned due to the use of 'is not' instead of '!='.",
    [
        "This problem lead to residue referencing when centering proteins not working properly, but did not impact the residue numbering of written files.",
    ],
]

documentation_changes = [
]

tutorial_changes = [
]

other_changes = [
]

def version_change_writer(iterable, recursion_depth = 0):
    list_of_strings = []
    for i in iterable:
        if type(i) == str:
            ### Headers
            if recursion_depth == 0:
                list_of_strings.append(i)
            ### Changes. -1 to have no spaces for first recursion. Two spaces "  " to fit with GitHub list formatting.
            else:
                list_of_strings.append("  " * (recursion_depth - 1) + "-" + " " + i)

        elif type(i) in [list, tuple]:
            list_of_strings.extend(version_change_writer(i, recursion_depth + 1))
    return list_of_strings

### Extra empty "" is to add a blank line between sections
all_changes = []
if len(major_changes) > 0:
    all_changes += ["Major changes:", major_changes, ""]

if len(minor_changes) > 0:
    all_changes += ["Minor changes:", minor_changes, ""]

if len(bug_fixes) > 0:
    all_changes += ["Bug fixing:", bug_fixes, ""]

if len(documentation_changes) > 0:
    all_changes += ["Documentation changes:", documentation_changes, ""]

if len(tutorial_changes) > 0:
    all_changes += ["Tutorial changes:", tutorial_changes, ""]

if len(other_changes) > 0:
    all_changes += ["Other changes:", other_changes, ""]

if len(all_changes) > 0:
    all_changes = all_changes[:-1] # Removes the last ""

version_changes_list = version_change_writer(all_changes)
__version_changes__ = "\n".join(version_changes_list)
__changes__         = __version_changes__
__changelog__       = __version_changes__
