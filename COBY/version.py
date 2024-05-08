__version__="0.1.2"

new_features = [
    "New:",
    [
        "It is now possible to pip install COBY"
    ],
]

feature_changes = [
    "Changes:",
    [
        "CGSB has been renamed to COBY (COarse grained system B[Y]uilder)",
        "Changed 'CGSB Logo' tutorial in 'Tutorial_advanced' to 'COBY Logo'.",
    ],
]

bug_fixes = [
    "Bug fixing:",
    [
        "Changed assert statements for 'itp_input' so that it tells you if a file does not exist instead of showing 'incorrect argument'.",
        "Fixed error with backup functionaly when output files are in the same folder as the script from which COBY is run.",
        "Added momentary fix to never ending 'CREATING LIPID GRID'. Not ideal fix but fine for now.",
        "Fixed membranes being able to span outside PBC when manually designating center and x/y-length.",
    ],
]

def version_change_writer(iterable, recursion_depth = 0):
    list_of_strings = []
    for i in iterable:
        if type(i) == str:
            list_of_strings.append("    " * recursion_depth + i)
        elif type(i) in [list, tuple]:
            list_of_strings.extend(version_change_writer(i, recursion_depth + 1))
    return list_of_strings

all_changes = []
if len(new_features) > 1:
    all_changes += new_features + [""]
if len(feature_changes) > 1:
    all_changes += feature_changes + [""]
if len(bug_fixes) > 1:
    all_changes += bug_fixes

version_changes_list = version_change_writer(all_changes)
version_changes_str = "\n".join(version_changes_list)

def version_changes():
    print(version_changes_str)

