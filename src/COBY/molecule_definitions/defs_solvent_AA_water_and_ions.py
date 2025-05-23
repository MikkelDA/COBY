solvent_defs = {}
pos_ion_defs = {}
neg_ion_defs = {}
solvent_metadata = {}
pos_ion_metadata = {}
neg_ion_metadata = {}

note_text = [
    "Please note the following when solvating with atomistic molecules:",
    "    " + "I suggest adding the following subarguments to solvent arguments when using atomistic solvent",
    "    " + "    " + "kick:0.03 bead_radius:0.12",
    "    " + "E.g. 'solv:TIP3:params:AA pos:SOD:params:AA neg:CLA:params:AA kick:0.04 bead_radius:0.12'",
    "    " + "    " + "The bead radius is a rough estimation obtained by comparing the output of COBY and CharmmGUI's 'Solution Builder' for systems containing AA proteins in solution.",
    "    " + "    " + "The kick is simply bead_radius/4",
]
params = "AA"
solvent_metadata[params] = {
    "Description": "This parameter library contains atomistic tip3p, tip4p and tip5p water where the molecule names are 'TIP3', 'TIP4' and 'TIP5', respectively.",
    "Note": note_text,
}
params = "AA_tip3p"
solvent_metadata[params] = {
    "Description": "This parameter library contains atomistic tip3p water where the molecule name is 'SOL'.",
    "Note": note_text,
}
params = "AA_tip4p"
solvent_metadata[params] = {
    "Description": "This parameter library contains atomistic tip4p water where the molecule name is 'SOL'.",
    "Note": note_text,
}
params = "AA_tip5p"
solvent_metadata[params] = {
    "Description": "This parameter library contains atomistic tip5p water where the molecule name is 'SOL'.",
    "Note": note_text,
}
params = "AA_Charmm"
solvent_metadata[params] = {
    "Description": "This parameter library contains atomistic tip3p water where the molecule name is 'SOL' and the atom names use the ones from the Charmm force field.",
    "Note": note_text,
}
params = "AA"
pos_ion_metadata[params] = {
    "Description": "This parameter library contains atomistic positively charge ions.",
    "Note": note_text,
}
params = "AA"
neg_ion_metadata[params] = {
    "Description": "This parameter library contains atomistic positively charge ions.",
    "Note": note_text,
}


### molar_mass: [g/mol]
### density: [g/cm^3]

################################
### ATOMISTIC WATER AND IONS ###
################################

### I suggest adding the following subarguments to solvent arguments when using atomistic solvent
### "kick:0.03 bead_radius:0.12"
### E.g. "solv:TIP3:params:AA pos:SOD:params:AA neg:CLA:params:AA kick:0.04 bead_radius:0.12"
### The bead radius is a rough estimation obtained by comparing the output of COBY and CharmmGUI's "Solution Builder" for systems containing AA proteins in solution.
### The kick is simply bead_radius/4

### All AA water models under seperate parameter libraries with the same residue name of "SOL"
### Atom positions
#      O
#  H1     H2
params = "AA_tip3p"
solvent_defs[params] = {
    "SOL": {
        "mapping_ratio": 1, "density": 0.99669, "molar_mass": 18.01528,
        "beads": [
            {"name": "OW",  "x":  0,     "y": 0,     "z": 0, "charge": -0.834},
            {"name": "HW1", "x":  0.074, "y": 0.064, "z": 0, "charge":  0.417},
            {"name": "HW2", "x": -0.074, "y": 0.064, "z": 0, "charge":  0.417},
        ],
        "tags": ("atomistic", "water"),
    },
}

### Atom positions
#      O
#      M
#  H1     H2
params = "AA_tip4p"
solvent_defs[params] = {
    "SOL": {
        "mapping_ratio": 1, "density": 0.99669, "molar_mass": 18.01528,
        "beads": [
            {"name": "OW",  "x":  0,     "y": 0,     "z": 0, "charge":  0},
            {"name": "HW1", "x":  0.074, "y": 0.064, "z": 0, "charge":  0.52},
            {"name": "HW2", "x": -0.074, "y": 0.064, "z": 0, "charge":  0.52},
            {"name": "MW",  "x":  0,     "y": 0.032, "z": 0, "charge": -1.04},
        ],
        "tags": ("atomistic", "water"),
    },
}

### Atom positions
#   LP1 LP2
#      O
#  H1     H2
params = "AA_tip5p"
solvent_defs[params] = {
    "SOL": {
        "mapping_ratio": 1, "density": 0.99669, "molar_mass": 18.01528,
        "beads": [
            {"name": "OW",  "x":  0,     "y":  0,     "z": 0, "charge":  0},
            {"name": "HW1", "x":  0.074, "y":  0.064, "z": 0, "charge":  0.241},
            {"name": "HW2", "x": -0.074, "y":  0.064, "z": 0, "charge":  0.241},
            {"name": "LP1", "x":  0.02,  "y":  0.02,  "z": 0, "charge": -0.241},
            {"name": "LP1", "x": -0.02,  "y": -0.02,  "z": 0, "charge": -0.241},
        ],
        "tags": ("atomistic", "water"),
    },
}

### All AA water models under a single parameter library with residue names being "TIP3", "TIP4" and "TIP5"
params = "AA"
solvent_defs[params] = {
    "TIP3": {
        "mapping_ratio": 1, "density": 0.99669, "molar_mass": 18.01528,
        "beads": [
            {"name": "OW",  "x":  0,     "y": 0,     "z": 0, "charge": -0.834},
            {"name": "HW1", "x":  0.074, "y": 0.064, "z": 0, "charge":  0.417},
            {"name": "HW2", "x": -0.074, "y": 0.064, "z": 0, "charge":  0.417},
        ],
        "tags": ("atomistic", "water"),
    },
    "TIP4": {
        "mapping_ratio": 1, "density": 0.99669, "molar_mass": 18.01528,
        "beads": [
            {"name": "OW",  "x":  0,     "y": 0,     "z": 0, "charge":  0},
            {"name": "HW1", "x":  0.074, "y": 0.064, "z": 0, "charge":  0.52},
            {"name": "HW2", "x": -0.074, "y": 0.064, "z": 0, "charge":  0.52},
            {"name": "MW",  "x":  0,     "y": 0.032, "z": 0, "charge": -1.04},
        ],
        "tags": ("atomistic", "water"),
    },
    "TIP5": {
        "mapping_ratio": 1, "density": 0.99669, "molar_mass": 18.01528,
        "beads": [
            {"name": "OW",  "x":  0,     "y":  0,     "z": 0, "charge":  0},
            {"name": "HW1", "x":  0.074, "y":  0.064, "z": 0, "charge":  0.241},
            {"name": "HW2", "x": -0.074, "y":  0.064, "z": 0, "charge":  0.241},
            {"name": "LP1", "x":  0.02,  "y":  0.02,  "z": 0, "charge": -0.241},
            {"name": "LP1", "x": -0.02,  "y": -0.02,  "z": 0, "charge": -0.241},
        ],
        "tags": ("atomistic", "water"),
    },
}


### Uses the Charmm36 atom names
params = "AA_Charmm"
solvent_defs[params] = {
    "TIP3": {
        "mapping_ratio": 1, "density": 0.99669, "molar_mass": 18.01528,
        "beads": [
            {"name": "OH2", "x":  0,     "y": 0,     "z": 0, "charge": -0.834},
            {"name": "H1",  "x":  0.074, "y": 0.064, "z": 0, "charge":  0.417},
            {"name": "H2",  "x": -0.074, "y": 0.064, "z": 0, "charge":  0.417},
        ],
        "tags": ("atomistic", "water"),
    },
}


params = "AA"
pos_ion_defs[params] = {
    ### Monovalent
    "SOD": {
        "beads": [{"name": "SOD", "x": 0, "y":  0, "z": 0, "charge": 1}],
        "tags": ("atomistic", "ions"),
    },
}

params = "AA"
neg_ion_defs[params] = {
    ### Monovalent
    "CLA": {
        "beads": [{"name": "CLA", "x": 0, "y":  0, "z": 0, "charge": -1}],
        "tags": ("atomistic", "ions"),
    },
}
