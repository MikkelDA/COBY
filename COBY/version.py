__version__="0.2.0"

new_features = [
    "Implemented new molecule fragment builder argument (molecule_builder) (see documentation for details).",
    [
        "Allows for the dynamic construction of molecules based on different fragments.",
        "The fragments can either be specified from a prebuilt library (e.g. heads and linkers for lipids) or dynamically built (e.g. tails for lipids).",
        "Currently only contains glycerophospholipids with PC, PE, PA, PG and PS heads, but more lipid types will be added in the future.",
        "No non-lipid molecules are currently implemented but the system allows for them as molecules built with molecule_builder are added to all libraries (lipid, solvent and ions).",
        "More lipid builder fragments can be loaded with the argument 'import_defintions' / 'import_defs'.",
    ],
    "Implemented new molecule builder program COBY.Crafter separate from COBY.COBY (see documentation for details).",
    [
        "Accessed from python with 'COBY.Crafter' (as opposed to 'COBY.COBY').",
        "Accessed from terminals by using '--program Crafter' (--program defaults to COBY if it is not specified)",
        "Instead of creating simulation-ready systems, 'COBY.Crafter' instead creates systems containing a single molecule.",
        "Molecules can be obtained from any inbuilt library (lipid, solvent/solute, neg/pos ions), the lipid fragment builder or read from topology files",
    ],
]

feature_changes = [
    "Topology reader now reads COBY lipid fragment builder arguments written in topology files (see documentation for details).",
    [
        "Arguments must start with ';@COBY' and be followed by the lipid fragment builder argument.",
        "Lipids from topology files are by default placed in the 'TOP' parameter library though this can be changed with the 'params' subargument in the lipid fragment builder argument.",
    ],
    "Changed how molecule definition (solvent_defs, ion_defs, lipid_defs) dictionaries are formatted to be more similar to the new fragment definition dictionaries.",
]

minor_changes = [
    "Changed how some terminal printouts appear.",
    "Added lipid scaffold for all tail combinations for M2/M3 release phospholipids and sphingolipids with simple heads.",
    "Added documentation for nanodisc functionality.",
]

bug_fixes = [
    "Fixed definitions importer (import_defintions / import_defs) so it properly updates molecule libraries.",
    [
        "It no longer replaces the old library but instead performs deep merging of the old and new nested dictionaries ensuring that previous data is maintained.",
        "If duplicate data (with non-identical) values is detected, then it simply overwrites the prior data with the new but throws a warning to let you know.",
    ],
    "Various minor bug fixes here and there.",
]

tutorial_changes = [
    "Added the 'Tutorial_paper.ipynb' notebook which includes the COBY commands used to create the systems from our upcoming paper.",
    [
        "Also added files needed for the creation of the paper figures.",
    ],
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

