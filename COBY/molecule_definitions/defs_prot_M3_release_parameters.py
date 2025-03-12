prot_defs    = {}

###########################
### PROTEIN CHARGE DATA ###
###########################
params = "default"
prot_defs[params] = {}
prot_defs[params]["charges"] = {
    ### ### 1 atom residues
    "GLY":  {"BB": 0},
    
    ### ### 2 atom residues
    ### Residues without variants
    "ALA":  {"BB": 0, "SC1": 0},
    "CYS":  {"BB": 0, "SC1": 0},
    "VAL":  {"BB": 0, "SC1": 0},
    "LEU":  {"BB": 0, "SC1": 0},
    "ILE":  {"BB": 0, "SC1": 0},
    "MET":  {"BB": 0, "SC1": 0},
    "PRO":  {"BB": 0, "SC1": 0},
    "HYP":  {"BB": 0, "SC1": 0},
    "ASN":  {"BB": 0, "SC1": 0},
    "GLN":  {"BB": 0, "SC1": 0},
    "THR":  {"BB": 0, "SC1": 0},
    "SER":  {"BB": 0, "SC1": 0},
    ### ASP variants
    "ASP":  {"BB": 0, "SC1": -1},
    "ASPP": {"BB": 0, "SC1": 0}, # Neutral ASP
    "ASH":  {"BB": 0, "SC1": 0}, # Neutral ASP
    ### Glutamate variants
    "GLU":  {"BB": 0, "SC1": -1},
    "GLUP": {"BB": 0, "SC1": 0}, # Neutral GLU
    "GLH":  {"BB": 0, "SC1": 0}, # Neutral GLU
    
    ### ### 3 atom residues
    "ARG":  {"BB": 0, "SC1": 0, "SC2": 1},
    ### Lysine variants
    "LYS":  {"BB": 0, "SC1": 0, "SC2": 1},
    "LSN":  {"BB": 0, "SC1": 0, "SC2": 0}, # Neutral LYS
    "LYN":  {"BB": 0, "SC1": 0, "SC2": 0}, # Neutral LYS
    
    ### ### 4 atom residues
    "PHE":  {"BB": 0, "SC1": 0, "SC2": 0, "SC3": 0},
    ### Histidines
    "HIS":  {"BB": 0, "SC1": 0, "SC2": 0, "SC3": 0},
    "HIE":  {"BB": 0, "SC1": 0, "SC2": 0, "SC3": 0},
    "HSE":  {"BB": 0, "SC1": 0, "SC2": 0, "SC3": 0},
    "HSD":  {"BB": 0, "SC1": 0, "SC2": 0, "SC3": 0},
    "HID":  {"BB": 0, "SC1": 0, "SC2": 0, "SC3": 0},
    "HSP":  {"BB": 0, "SC1": 0, "SC2": 0.5, "SC3": 0.5},
    "HIP":  {"BB": 0, "SC1": 0, "SC2": 0.5, "SC3": 0.5},
    
    ### ### 5 atom residues
    "TYR":  {"BB": 0, "SC1": 0, "SC2": 0, "SC3": 0, "SC4": 0},
    
    ### ### 6 atom residues
    "TRP":  {"BB": 0, "SC1": 0, "SC2": 0, "SC3": 0, "SC4": 0, "SC5": 0},
}
