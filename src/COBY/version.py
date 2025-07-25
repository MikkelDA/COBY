__version__="1.0.7"

major_changes = [
    "Added 'scale', 'scale_x', 'scale_y' and 'scale_z' subsubarguments for the 'lipid', 'solvent', 'solute', 'pos_ion' and 'neg_ion' subarguments",
    [
        "'scale_x:val', 'scale_y:val' and 'scale_z:val' can be used to scale the coordinates along any given axis of any given molecule.",
        "'scale:xval:yval:zval' can be used to scale all the x/y/z-coordinates of any given molecule in a single subsubargument. All three values must be given.",
    ],
    "Added 'rotate', 'rotate_x', 'rotate_y' and 'rotate_z' subsubarguments for the 'lipid', 'solvent', 'solute', 'pos_ion' and 'neg_ion' subarguments",
    [
        "'rotate_x:val', 'rotate_y:val' and 'rotate_z:val' can be used to scale the coordinates along any given axis of any given molecule.",
        "'rotate:xval:yval:zval' can be used to rotate all the x/y/z-coordinates of any given molecule in a single subsubargument. All three values must be given.",
    ],
    "Added 'interleaflet_buffer' subargument for membrane argument. Allows one to set the buffer distance between the middle of a membrane and the lipids. Can be given to the whole membrane or each leaflet individually.",
    "Added 'rotation_angles' subargument to solvation and flooding arguments. Allows one to set the allowed range of rotations for molecules along each axis.",
    [
        "Accepts 'xval1:xval2:yval1:yval2:zval1:zval2' or any combination of 'x:val1:val2', 'y:val1:val2' and 'z:val1:val2' (example: rotation_angles:z:-40:40:y:-50:80)",
        "The values must be within 360 of each other (i.e. it is required that 'maxval - minval <= 360'). All val1 and val2 values are -180 and 180, respectively, by default."
    ]
]

minor_changes = [
]

bug_fixes = [
    "Fixed a bug where the lipid grouping algorithm would separate lipids into different groups when they should be contained within a single group.",
    "Fixed a bug where no radius groups were created causing no lipids to be placed but did not result in a crash.",
]

documentation_changes = [
    "Added documentation for the new 'scale', 'rotate', 'interleaflet_buffer' and 'rotation_angles' functionalities."
]

tutorial_changes = [
]

other_changes = [
    "Added a 'known issues' section to the front page README file.",
    [
        "Added 'Recursion depth crash on Mac' to the list of known issues."
    ]
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
