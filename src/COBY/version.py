__version__="0.2.7"

major_changes = [
]

minor_changes = [
    "Cleaned up the code/math a bit for the lipid placement optimizer. Also added handling for when two points are placed directly on top of each other during the placement optimization.",
    "Slightly increased performance of solvent placement algorithm by reducing the number of calls to debug prints.",
    "Added the ability to give lipid-specific APL values. Lipid APL values default to the leaflet-wide APL value.",
    [
        "Example: membrane = 'apl:0.6 lipid:POPC:5:apl:0.5 lipid:CHOL:2:apl:0.4 lipid:DOPE:2'",
        [
            "POPC and CHOL would have APLs of 0.5 and 0.4, respectively, while DOPE would have an APL of 0.6, inherited from the leaflet (membrane) APL.",
            "The leaflet APL (used in later calculations) is then adjusted to become sum([apl*ratio for apl, ratio in zip(apls, ratios)]) / sum(ratios) = 0.5",
        ],
    ],
]

bug_fixes = [
    "COBY.Library can now be called from the terminal without any other arguments being given.",
    "Fixed crash caused by file backup method that occured on Windows systems due to improper OS-agnostic path handling.",
    "Fixed crash caused by creating a membranes with lipids having a ratio of zero, e.g. 'lipid:POPC:1 lipid:POPE:0'. COBY now properly ignores the lipid types that end up with zero lipids.",
]

documentation_changes = [
    "Fixed wrong documentation for 'upbead'/'downbead' for molecule_import method.",
    "Added documentation for new apl subsubargument to be used in lipid subarguments",
]

tutorial_changes = [
    "Renamed tutorial notebook files.",
]

other_changes = [
    "Updated both the general GitHub README.md and the GitHub Tutorial README.md.",
    "Added uncompressed program folder to GitHub",
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
version_changes_str = "\n".join(version_changes_list)

def version_changes():
    print(version_changes_str)

### Abbreviations
changes   = version_changes
changelog = version_changes

