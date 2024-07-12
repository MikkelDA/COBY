__version__="0.2.1"

new_features = [
]

feature_changes = [
]

minor_changes = [
    "Changed the definitions import argument 'import_definitions' / 'import_defs' to 'import_library'. No functional changes beyond the name change."
]

bug_fixes = [
    "Fixed an unintended error warning when loading COBY without lipid definitions (there are none in the inbuilt library)."
]

tutorial_changes = [
    "Corrected 'import_definitions' / 'import_defs' to 'import_library' per the argument name change."
]

def version_change_writer(iterable, recursion_depth = 0):
    list_of_strings = []
    for i in iterable:
        if type(i) == str:
            if recursion_depth == 0:
                list_of_strings.append("    " * recursion_depth + i)
            else:
                list_of_strings.append("    " * recursion_depth + "-" + " " + i)

        elif type(i) in [list, tuple]:
            list_of_strings.extend(version_change_writer(i, recursion_depth + 1))
    return list_of_strings

### Extra empty "" is to add a blank line between sections
all_changes = []
if len(new_features) > 0:
    all_changes += ["New features:", new_features, ""]

if len(feature_changes) > 0:
    all_changes += ["Feature changes:", feature_changes, ""]

if len(minor_changes) > 0:
    all_changes += ["Minor changes:", minor_changes, ""]

if len(bug_fixes) > 0:
    all_changes += ["Bug fixing:", bug_fixes, ""]

if len(tutorial_changes) > 0:
    all_changes += ["Tutorial changes:", tutorial_changes]

version_changes_list = version_change_writer(all_changes)
version_changes_str = "\n".join(version_changes_list)

def version_changes():
    print(version_changes_str)

### Abbreviation
changes = version_changes

