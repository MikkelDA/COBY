__version__="0.1.1"

new_features = [
    "New features:",
    [
        "Added 'CGSB.version' value to check current CGSB version.",
        "Added 'CGSB.version_changes()' function to see changes from previous version to current version.",
    ],
]

feature_changes = [
    "Feature changes:",
    [
        "Solvations can now be done without ions.",
        "Output filenames now have their extension cases checked (e.g., '.GRO' is read as '.gro', such that the file is not renamed to '.GRO.gro').",
    ],
]

bug_fixes = [
    "Bug fixing:",
    [
        "Fixed error when manually designating lipid insertion algorithm.",
        "No longer adds '.gro' or '.pdb' if filenames already contains it for 'out_gro' and 'out_pdb' arguments.",
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

