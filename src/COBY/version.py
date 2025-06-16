__version__="1.0.2"

major_changes = [
    "Made it possible to apply all 'grid_maker' arguments to specific leaflets instead of always being applied to both leaflets.",
]

minor_changes = [
    "Slightly modified the algorithm for finding initial lipid placements.",
    "Added 'gm_rlo'/'grid_maker_lipid_distribution' subargument to the Membrane argument allowing for the reversal of the lipid placement order.",
]

bug_fixes = [
]

documentation_changes = [
    "Added documentation for gm_rlo'/'grid_maker_lipid_distribution' subargument."
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
