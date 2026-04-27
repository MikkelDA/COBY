__version__="1.0.12"

major_changes = [
    "Added two new subarguments for the 'membrane' argument.",
    [
        "'grid_maker_maximum_radius'/'gm_maxr': Sets the maximum radius for separate radius grouping. The default value is 4.0 [nm]. Radius groups that are above this maximum value are merged together.",
        "'grid_maker_minimum_radius'/'gm_minr': Sets the minimum radius for separate radius grouping. The default value is 2.5 [nm]. Radius groups that are below this minimum value are merged together.",
        "The effective change is that (with default values), lipids with radii below 2.5 nm (or above 4.0 nm) are always placed in the same lipid insertion step. This should fix various bugs related to lipid insertion."
    ],
]

minor_changes = [
    "COBY now prints the current version when being run. This is also written to the log file.",
    "The log file is now written while COBY is running.",
    "Added 'there are no tags in this parameter library' print for COBY.Library when no tags are present.",
]

bug_fixes = [
    "Various bugs (more like unintended edge-case behaviour) related lipid radius grouping and lipid insertion.",
]

documentation_changes = [
    "Added documentation for 'grid_maker_maximum_radius/gm_maxr' and 'grid_maker_minimum_radius/gm_minr'.",
]

tutorial_changes = [
    "Various fixes for 3_Manuscript_tutorial.ipynb:",
    [
        "Added DOI to the top of the Notebook.",
        "Fixed incorrect subargument in the 'c) Multiple solvent spaces' tutorial.",
        "Removed duplicate 'b)' section located in the 'c)' section.",
        "Corrected system names so that they correspond to the correct subfigure in the publication.",
    ],
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
