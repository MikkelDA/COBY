__version__="1.0.16"

major_changes = [
    "Redid the parameter library names for for Martini 2/3 lipid, solvent/ion and protein parameter libraries.",
    [
        "Changed from arbitrary names to standardised {force field abbreviation}_{year of publication} format.",
        "Changed default lipid library from 'M2_2021' (previously called 'default') to 'M3_2025' (previously called 'M3_LTF')",
    ],
]

minor_changes = [
    "Added possibility to do random protein rotations when using the 'rotate' and 'r[x/y/z]' subarguments. See the documentation section 6.18.2 for examples.",
]

bug_fixes = [
    "Added the ability to set 'keep_residue_numbering' when running COBY from the terminal.",
]

documentation_changes = [
    "Added list of all included parameter libraries to section 6.12 'Setting parameter libraries: params'",
]

tutorial_changes = [
    "Changed a few tutorials to match the new parameter library names",
]

other_changes = [
    "Added matplotlib and pyrecorder as dependencies in the pyproject.toml.",
    "Removed checks to the membrane grid plotter that tested for a pyrecorder installation.",
    "Added urls for the COBY documentation and the COBY paper to the [project.urls] section in the pyproject.toml.",
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
