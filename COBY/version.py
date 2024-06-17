__version__="0.1.6"

new_features = [
    "Topology reading:",
    [
        "Re-implemented ifdef/ifndef processeing.",
        "Added 'define' subargument to 'itp_input' argument where defines used in the topology files can be specified.",
        "An '#ifdef' statement will only be processed if the variable has been defined.",
        "An '#ifndef' statement will only be processed if the variable has not been defined.",
        "Variables can also be defined in the topology file with '#define VARNAME'",
    ]
]

feature_changes = [
]

bug_fixes = [
    "Fixed bug with 'include:file.itp' for 'itp_input' argument where it would only include the last itp file given with 'include'",
]

minor_changes = [
]

tutorial_changes = [
]

def version_change_writer(iterable, recursion_depth = 0):
    list_of_strings = []
    for i in iterable:
        if type(i) == str:
            list_of_strings.append("    " * recursion_depth + i)
        elif type(i) in [list, tuple]:
            list_of_strings.extend(version_change_writer(i, recursion_depth + 1))
    return list_of_strings

### Extra empty "" is to add a blank line between sections
all_changes = []
if len(new_features) > 0:
    all_changes += ["New:", new_features, ""]

if len(feature_changes) > 0:
    all_changes += ["Changes:", feature_changes, ""]

if len(bug_fixes) > 0:
    all_changes += ["Bug fixing:", bug_fixes, ""]

if len(minor_changes) > 0:
    all_changes += ["Minor changes:", minor_changes, ""]

if len(tutorial_changes) > 0:
    all_changes += ["Tutorial changes:", tutorial_changes]

version_changes_list = version_change_writer(all_changes)
version_changes_str = "\n".join(version_changes_list)

def version_changes():
    print(version_changes_str)

