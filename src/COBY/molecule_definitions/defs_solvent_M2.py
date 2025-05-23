solvent_defs = {}
solvent_metadata = {}

params = "M2"
solvent_metadata[params] = {
    "Description": "Solvent in this parameter library is based on the Martini 2 parameters as described in https://doi.org/10.1021/jp071097f.",
}

################
### SOLVENTS ###
################
### molar_mass: [g/mol]
### density: [g/cm^3]
params = "M2"
solvent_defs[params] = {}

### Water
solvent_defs[params].update({
    "W" : {
        "mapping_ratio": 4, "density": 0.99669, "molar_mass": 18.01528,
        "beads": [{"name": "W", "x": 0, "y": 0, "z": 0, "charge": 0}],
        "tags": ("martini2", "water"),
    },
    "WF" : {
        "mapping_ratio": 4, "density": 0.99669, "molar_mass": 18.01528,
        "beads": [{"name": "WF", "x": 0, "y": 0, "z": 0, "charge": 0}],
        "tags": ("martini2", "water"),
    },
})
