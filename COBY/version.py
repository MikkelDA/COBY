__version__="0.2.6"

major_changes = [
]

minor_changes = [
    "Changed parameter library name for lipid task force lipids to 'LTF' from 'beta3'",
    "Molecule fragment builder arguments written in topology files can now be read even if there are spaces between ';' and '@COBY' such as for '; @COBY'.",
    [
        "Also cleaned up the 'molecule_fragment_builder.py' file.",
    ],
]

bug_fixes = [
    "Molecule fragment builder arguments written in topology files no longer causes COBY to crash.",
    "'moleculetype' assignment is now handled correctly for molecule fragment builder arguments.",
]

documentation_changes = [
]

tutorial_changes = [
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

if len(all_changes) > 0:
    all_changes = all_changes[:-1] # Removes the last ""

version_changes_list = version_change_writer(all_changes)
version_changes_str = "\n".join(version_changes_list)

def version_changes():
    print(version_changes_str)

### Abbreviations
changes   = version_changes
changelog = version_changes

