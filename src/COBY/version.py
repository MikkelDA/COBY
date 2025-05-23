__version__="1.0.0"

major_changes = [
    "The paper for COBY has now been published. It can be found at https://doi.org/10.1021/acs.jcim.5c00069."
]

minor_changes = [
    "Added new 'lipid_extra' subargument for membranes. It allows the user to add extra lipids to a membrane/leaflet after the initial lipid calculations have been done. See the documentation for further explanation of how it works.",
    "COBY.Library() can now also be used to investigate fragment molecule types.",
    "Added metadata dictionaries for lipid, solvent, ion, protein parameter libraries and fragment molecule types.",
    [
        "This metadata is printed when accessing the parameter library / fragment library with COBY.Library().",
    ],
    "Added -doi/-DOI / COBY.doi/COBY.DOI and -citation/COBY.citation arguments/classes. The first prints the DOI for the COBY paper and the second prints a citation for the paper in the bibtex style.",
        [
            "Also cleaned up a bit how the printing of version and version changes is handled. The strings for version",
        ],
]

bug_fixes = [
    "Lipid types that end up with zero (or a negative) number of lipids are now properly removed (as in, they won't be written to the .top file).",
    "Fixed some instances where COBY.Library() would crash when using 'r'/'return'",
]

documentation_changes = [
    "Added documentation for new 'lipid_extra' membrane subargument.",
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
