__version__="1.0.10"

major_changes = [
    "Changed the way molcule alignment works with the 'molecule_import' argument.",
    [
        "Added subargument 'alignment' which can be either False (default), 'manual' or 'principal'",
        [
            "False: The alignment algorithms are turned off.",
            "'manual': Aligns the molecule along the line that runs through the center of two groups of beads. This is how alignment used to work before if 'upbeads' and 'downbeads' were supplied.",
            [
                "Requires both 'upbeads' and 'downbeads' to be given",
            ],
            "'principal': Aligns the molecule along the first principal axis based on a principal component analysis.",
            [
                "It is required to specify either 'upbeads' or 'downbeads' (but not both). The algorithm will use that bead group to ensure that the molecule points in the correct direction.",
            ],
        ],
    ],
]

minor_changes = [
    "Added unit to final 'Time spent running COBY' print.",
    "Adde 'dir' subargument to 'plot_grid' argument, which allows one to change the name of the directory where the plotting files are placed. Default is 'grid_plots'.",
]

bug_fixes = [
    "Fixed potential crash caused by incorrect checking and reconfiguring of the 'grid resolution' during solvations.",
    "Solvent charge is now rounded to 5 decimals to prevent errors from float values.",
]

documentation_changes = [
    "Fixed an error in an example for 'hole:rectangle' using 'xradius' and 'yradius' instead of 'xlength' and 'ylength'.",
    "Changed existing documentation and added extra documentation for alignment with 'molecule_import'.",
    "Added documentation for new 'dir' subargument for the 'plot_grid' argument.",
]

tutorial_changes = [
    "Changed the tutorial '2_Advanced_tutorial - Import a lipid' to fit with the new syntax.",
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
