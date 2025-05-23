
lipid_scaffolds = {}
lipid_metadata = {}

params = "LTF"
lipid_metadata[params] = {
    "Description": "Lipids in this parameter library are based on the parameters from the lipid task force as described in https://doi.org/10.26434/chemrxiv-2024-8bjrr.",
}

### ### ### ### ### ################################################################## ### ### ### ### ###
### ### ### ### ### ### Phospholipids (glycerol and ether linker) and diglycerides ### ### ### ### ### ###
### ### ### ### ### ################################################################## ### ### ### ### ###

### Some tails have the same beads. The difference is in the topology (regular vs small beads).

phospholipid_tails = {
    ### No double bonds
    "DT": "C1A   C2A     -     -     -     -   C1B   C2B    -     -     -     - ",
    "DJ": "C1A   C2A     -     -     -     -   C1B   C2B    -     -     -     - ",

    "DU": "C1A   C2A   C3A     -     -     -   C1B   C2B   C3B    -     -     - ",
    "DM": "C1A   C2A   C3A     -     -     -   C1B   C2B   C3B    -     -     - ",
    
    "MP": "C1A   C2A   C3A   C4A     -     -   C1B   C2B   C3B    -     -     - ",
    "MS": "C1A   C2A   C3A   C4A     -     -   C1B   C2B   C3B    -     -     - ",
    
    "PM": "C1A   C2A   C3A     -     -     -   C1B   C2B   C3B   C4B     -    - ",
    "SM": "C1A   C2A   C3A     -     -     -   C1B   C2B   C3B   C4B     -    - ",

    "DP": "C1A   C2A   C3A   C4A     -     -   C1B   C2B   C3B   C4B     -    - ",
    "DS": "C1A   C2A   C3A   C4A     -     -   C1B   C2B   C3B   C4B     -    - ",
    "PS": "C1A   C2A   C3A   C4A     -     -   C1B   C2B   C3B   C4B     -    - ",

    "DK": "C1A   C2A   C3A   C4A   C5A     -   C1B   C2B   C3B   C4B   C5B    - ",
    "DB": "C1A   C2A   C3A   C4A   C5A     -   C1B   C2B   C3B   C4B   C5B    - ",

    "DX": "C1A   C2A   C3A   C4A   C5A     -   C1B   C2B   C3B   C4B   C5B   C6B",
    "DC": "C1A   C2A   C3A   C4A   C5A   C6A   C1B   C2B   C3B   C4B   C5B   C6B",

    ### One double bond on each tail
    "DR": "C1A   D2A   C3A     -     -     -   C1B   D2B   C3B    -     -     - ",

    "DY": "C1A   C2A   D3A   C4A     -     -   C1B   C2B   D3B   C4B    -     - ",
    "DO": "C1A   D2A   C3A   C4A     -     -   C1B   D2B   C3B   C4B    -     - ",
    "DV": "C1A   C2A   D3A   C4A     -     -   C1B   C2B   D3B   C4B    -     - ",
    "YO": "C1A   D2A   C3A   C4A     -     -   C1B   C2B   D3B   C4B    -     - ",
    
    "OE": "C1A   C2A   D3A   C4A   C5A     -   C1B   D2B   C3B   C4B    -     - ",

    "DG": "C1A   C2A   D3A   C4A   C5A     -   C1B   C2B   D3B   C4B   C5B    - ",
    "DE": "C1A   C2A   D3A   C4A   C5A     -   C1B   C2B   D3B   C4B   C5B    - ",

    "DN": "C1A   C2A   C3A   D4A   C5A   C6A   C1B   C2B   C3B   D4B   C5B   C6B",

    ### Two double bonds on each tail
    "DL": "C1A   D2A   D3A   C4A     -     -   C1B   D2B   D3B   C4B    -     - ",

    ### Three double bonds on each tail
    "DF": "C1A   D2A   D3A   D4A     -     -   C1B   D2B   D3B   D4B    -     - ",
    "DA": "C1A   D2A   D3A   D4A   C5A     -   C1B   D2B   D3B   D4B   C5B    - ",

    ### Five double bonds on each tail
    "DD": "D1A   D2A   D3A   D4A   D5A     -   D1B   D2B   D3B   D4B   D5B    - ",

    ### One double bond on tail A (SN2), no double bonds on tail B (SN1)
    "PY": "C1A   C2A   D3A   C4A     -     -   C1B   C2B   C3B   C4B    -     - ",
    "PO": "C1A   D2A   C3A   C4A     -     -   C1B   C2B   C3B   C4B    -     - ",
    "SO": "C1A   D2A   C3A   C4A     -     -   C1B   C2B   C3B   C4B    -     - ",
    "PE": "C1A   C2A   D3A   C4A   C5A     -   C1B   C2B   C3B   C4B    -     - ",

    ### Two double bonds on tail A (SN2), no double bonds on tail B (SN1)
    "PL": "C1A   D2A   D3A   C4A     -     -   C1B   C2B   C3B   C4B    -     - ",
    "SL": "C1A   D2A   D3A   C4A     -     -   C1B   C2B   C3B   C4B    -     - ",
    "PI": "C1A   C2A   D3A   D4A   C5A     -   C1B   C2B   C3B   C4B    -     - ",

    ### Three double bonds on tail A (SN2), no double bonds on tail B (SN1)
    "PF": "C1A   D2A   D3A   D4A     -     -   C1B   C2B   C3B   C4B    -     - ",
    "PQ": "C1A   D2A   D3A   D4A   C5A     -   C1B   C2B   C3B   C4B    -     - ",
    "PA": "C1A   D2A   D3A   D4A   C5A     -   C1B   C2B   C3B   C4B    -     - ",
    "SA": "C1A   D2A   D3A   D4A   C5A     -   C1B   C2B   C3B   C4B    -     - ",

    ### Five double bonds on tail A (SN2), no double bonds on tail B (SN1)
    "PD": "D1A   D2A   D3A   D4A   D5A     -   C1B   C2B   C3B   C4B    -     - ",
    "SD": "D1A   D2A   D3A   D4A   D5A     -   C1B   C2B   C3B   C4B    -     - ",

    ### Two double bonds on tail A (SN2), one double bond on tail B (SN1)
    "OL": "C1A   D2A   D3A   C4A     -     -   C1B   D2B   C3B   C4B    -     - ",

    ### Five double bonds on tail A (SN2), one double bond on tail B (SN1)
    "OD": "D1A   D2A   D3A   D4A   D5A     -   C1B   D2B   C3B   C4B    -     - ",

    ### Three double bonds on tail A (SN2), two double bond on tail B (SN1)
    "LF": "C1A   D2A   D3A   D4A     -     -   C1B   D2B   D3B   C4B    -     - ",
}



### Phospholipids with simple head structures
#
#                     B15-B16-B17-B18-B19-B20
#                    /
#    H02   H05    L08
#   /  |   |  \    |
# H01-H03-H04-H06-L07-A09-A10-A11-A12-A13-A14

phospholipid_heads_simple = {
          #   H01   H02   H03   H04   H05   H06
    "PA": "    -     -     -     -     -    PO4",
    "PC": "    -     -     -    NC3    -    PO4",
    "PE": "    -     -     -    NH3    -    PO4",
    "PG": "    -     -     -    GL0    -    PO4",
    "PS": "    -     -     -    CNO    -    PO4",
    "DG": "    -     -     -     -     -     OH", # Diglyceride
}

lipid_type, params = "diacyl_glycerols_LTF_release", "LTF"
lipid_scaffolds[(lipid_type, params)] = {}   #  H01   H02   H03   H04   H05   H06   L07   L08   A09   A10   A11   A12   A13   A14   B15   B16   B17   B18   B19   B20
lipid_scaffolds[(lipid_type, params)]["x"] = (    0,   .5,    0,    0,   .5,    0,    0,   .5,    0,    0,    0,    0,    0,    0,    1,    1,    1,    1,    1,    1)
lipid_scaffolds[(lipid_type, params)]["y"] = (    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0)
lipid_scaffolds[(lipid_type, params)]["z"] = (   10,    9,    9,    8,    8,    7,    6,    6,    5,    4,    3,    2,    1,    0,    5,    4,    3,    2,    1,    0)
lipid_scaffolds[(lipid_type, params)]["charges"] = (("NC3", 1), ("NH3", 1), ("PO4", -1))
lipid_scaffolds[(lipid_type, params)]["bd"] = (0.25, 0.25, 0.3)
lipid_scaffolds[(lipid_type, params)]["lipids"] = {}

### Linker does not change lipid name (e.g. POPC is still called POPC if it has an ether linker)
### Adds glycerol linker as standard lipid when requesting phospholipids (e.g. Accessible using POPC)
phospholipid_linkers = {
    "GL": "GL1   GL2", # Glycerol
}
diacyl_glycerols_lipids = {
    (tail_key + head_key, "beads"): "   ".join([head_beads, link_beads, tail_beads])
    for head_key, head_beads in phospholipid_heads_simple.items()
    for link_key, link_beads in phospholipid_linkers.items()
    for tail_key, tail_beads in phospholipid_tails.items()
}
lipid_scaffolds[(lipid_type, params)]["lipids"].update(diacyl_glycerols_lipids)

### Adds both glycerol and ether linkers as ".GL" / ".ET" extensions to lipid names (e.g. Accessible using POPC.GL and POPC.ET, respectively)
phospholipid_linkers = {
    "GL": "GL1   GL2", # Glycerol
    "ET": "ET1   ET2", # Ether
}
diacyl_glycerols_lipids = {
    (tail_key + head_key + "." + link_key, "beads"): "   ".join([head_beads, link_beads, tail_beads])
    for head_key, head_beads in phospholipid_heads_simple.items()
    for link_key, link_beads in phospholipid_linkers.items()
    for tail_key, tail_beads in phospholipid_tails.items()
    if head_key != "DG" # No diglyceride for ether linkers so no point in having ".GL" versions of them
}
lipid_scaffolds[(lipid_type, params)]["lipids"].update(diacyl_glycerols_lipids)


### Inositol phospholipids
# Very rough ascii-art due to small differences in positions
#
# H06\                 B17-B18-B19-B20-B21-B22
#    /H02\           /
# H07 |   H01     L10
#     |V04/  \     |
#     |  /    H05-L09-A11-A12-A13-A14-A15-A16
# H08-H03        

inositol_head_beads = {
    "PI": "    C1    C2    C3    C4   PO4     -     -     -",
    "P1": "    C1    C2    C3    C4   PO4    P3     -     -",
    "P2": "    C1    C2    C3    C4   PO4    P3    P4     -",
    "P3": "    C1    C2    C3    C4   PO4    P3    P4    P5",
    "P4": "    C1    C2    C3    C4   PO4     -    P4     -",
    "P5": "    C1    C2    C3    C4   PO4     -     -    P5",
    "P6": "    C1    C2    C3    C4   PO4     -    P4    P5",
    "P7": "    C1    C2    C3    C4   PO4    P3     -    P5",
}

inositol_head_charges = {
    "P1": [("P3", -2)],
    "P2": [("P3", -1.5), ("P4", -1.5)],
    "P3": [("P3", -1.3), ("P4", -1.4), ("P5", -1.3)],
    "P4": [("P4", -2)],
    "P5": [("P5", -2)],
    "P6": [("P4", -1.5), ("P5", -1.5)],
    "P7": [("P3", -1.5), ("P5", -1.5)],
}

lipid_type, params = "phosphoinositollipids_LTF_release", "LTF"
lipid_scaffolds[(lipid_type, params)] = {}   #  H01   H02   H03   V04   H05   H06   H07   H08   L09   L10   A11   A12   A13   A14   A15   A16   B17   B18   B19   B20   B21   B22
lipid_scaffolds[(lipid_type, params)]["x"] = (  0.5,  0.8, -0.3, 0.25,    0,    1,   .5, -0.3,    0,   .5,    0,    0,    0,    0,    0,    0,    1,    1,    1,    1,    1,    1)
lipid_scaffolds[(lipid_type, params)]["y"] = (    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0)
lipid_scaffolds[(lipid_type, params)]["z"] = (    8,  9.1,  8.8,  8.5,    7,   10,   10,   10,    6,    6,    5,    4,    3,    2,    1,    0,    5,    4,    3,    2,    1,    0)
lipid_scaffolds[(lipid_type, params)]["charges"] = ("PO4", -1)
lipid_scaffolds[(lipid_type, params)]["bd"] = (0.25, 0.25, 0.3)
lipid_scaffolds[(lipid_type, params)]["lipids"] = {}

### Adds glycerol linker as standard lipid when requesting phospholipids (e.g. Accessible using POPC)
phospholipid_linkers = {
    "GL": "GL1   GL2", # Glycerol
}
phosphoinositollipids_lipids = {
    (tail_key + head_key, "beads"): "   ".join([head_beads, link_beads, tail_beads])
    for head_key, head_beads in inositol_head_beads.items()
    for link_key, link_beads in phospholipid_linkers.items()
    for tail_key, tail_beads in phospholipid_tails.items()
}
lipid_scaffolds[(lipid_type, params)]["lipids"].update(phosphoinositollipids_lipids)

### Adds both glycerol and ether linkers as ".GL" / ".ET" extensions to lipid names (e.g. Accessible using POPC.GL and POPC.ET, respectively)
phospholipid_linkers = {
    "GL": "GL1   GL2", # Glycerol
    "ET": "ET1   ET2", # Ether
}
phosphoinositollipids_lipids = {
    (tail_key + head_key + "." + link_key, "beads"): "   ".join([head_beads, link_beads, tail_beads])
    for head_key, head_beads in inositol_head_beads.items()
    for link_key, link_beads in phospholipid_linkers.items()
    for tail_key, tail_beads in phospholipid_tails.items()
}
lipid_scaffolds[(lipid_type, params)]["lipids"].update(phosphoinositollipids_lipids)

### Creates lipid-specific charges
phosphoinositollipids_LTF_release_charges = {
    (tail_key + head_key + link_suffix, "charges"): head_charges
    for head_key, head_charges in inositol_head_charges.items()
    for link_suffix in ["", ".GL", "ET"]
    for tail_key in phospholipid_tails.keys()
}
lipid_scaffolds[(lipid_type, params)]["lipids"].update(phosphoinositollipids_LTF_release_charges)


### ### ### ### ### ########################################## ### ### ### ### ### 
### ### ### ### ### ### Phospholipids (plasmalogen linker) ### ### ### ### ### ###
### ### ### ### ### ########################################## ### ### ### ### ### 

### Some tails have the same beads. The difference is in the topology (regular vs small beads).
phospholipid_tails = {
    ### One double bond on tail A (SN2), no double bonds on tail B (SN1)
    "O": "C1A   D2A   C3A   C4A     -     -   C1B   C2B   C3B   C4B     -     -",

    ### Three double bonds on tail A (SN2), no double bonds on tail B (SN1)
    "A": "C1A   D2A   D3A   D4A   C5A     -   C1B   C2B   C3B   C4B     -     -",

    ### Five double bonds on tail A (SN2), no double bonds on tail B (SN1)
    "D": "D1A   D2A   D3A   D4A   D5A     -   C1B   C2B   C3B   C4B     -     -",
}

### Phospholipids with simple head structures
#
#                     B15-B16-B17-B18-B19-B20
#                    /
#    H02   H05    L08
#   /  |   |  \    |
# H01-H03-H04-H06-L07-A09-A10-A11-A12-A13-A14

phospholipid_heads_simple = {
    "C": "    -     -     -   NC3     -   PO4",
    "E": "    -     -     -   NH3     -   PO4",
}

phospholipid_linkers = {
    "PL": "GL1   PL2", # Plasmalogens
}

plasmalogen_lipids = {
    (tail_key + link_key + head_key, "beads"): "   ".join([head_beads, link_beads, tail_beads])
    for head_key, head_beads in phospholipid_heads_simple.items()
    for link_key, link_beads in phospholipid_linkers.items()
    for tail_key, tail_beads in phospholipid_tails.items()
}

lipid_scaffolds[(lipid_type, params)]["lipids"].update(plasmalogen_lipids)



### ### ### ### ### #################### ### ### ### ### ###
### ### ### ### ### ### Hydrocarbons ### ### ### ### ### ###
### ### ### ### ### #################### ### ### ### ### ###

### Hydrocarbons
"""
A1-A2-A3-A4-A5-A6
"""

lipid_type, params = "hydrocarbons_LTF_release", "LTF"
lipid_scaffolds[(lipid_type, params)] = {}   #  A1    A2    A3    A4    A5    A6
lipid_scaffolds[(lipid_type, params)]["x"] = (   0,    0,    0,    0,    0,    0)
lipid_scaffolds[(lipid_type, params)]["y"] = (   0,    0,    0,    0,    0,    0)
lipid_scaffolds[(lipid_type, params)]["z"] = (   5,    4,    3,    2,    1,    0)
lipid_scaffolds[(lipid_type, params)]["bd"] = (0.25, 0.25, 0.3)
lipid_scaffolds[(lipid_type, params)]["lipids"] = {
    ("C16", "beads"): "C1A   C2A   C3A   C4A    -     - ",
    ("C20", "beads"): "C1A   C2A   C3A   C4A   C5A    - ",
    ("C24", "beads"): "C1A   C2A   C3A   C4A   C5A   C6A",

    ("IC1", "beads"): "C1A   C2A   D3A   C4A   C5A    - ",
    ("IC2", "beads"): "C1A   C2A   D3A   D4A   C5A    - ",
    ("IC3", "beads"): "C1A   C2A   D3A   D4A   D5A    - ",

    ### IC4 and IC5 have identical structures but different topologies
    ("IC4", "beads"): "C1A   D2A   D3A   D4A   D5A    - ",
    ("IC5", "beads"): "C1A   D2A   D3A   D4A   D5A    - ",

    ("IC6", "beads"): "C1A   D2A   C3A   C4A   C5A    - ",
    ("IC7", "beads"): "C1A   D2A   C3A   D4A   D5A    - ",
    ("IC8", "beads"): "C1A   D2A   C3A   D4A   C5A    - ",
    ("IC9", "beads"): "C1A   D2A   C3A   C4A   D5A    - ",
}


### ### ### ### ### ################### ### ### ### ### ###
### ### ### ### ### ### Fatty Acids ### ### ### ### ### ###
### ### ### ### ### ################### ### ### ### ### ###

### Hydrocarbons
#
# H1-A2-A3-A4-A5-A6-A7

lipid_type, params = "fattyacids_LTF_release", "LTF"
lipid_scaffolds[(lipid_type, params)] = {}   #  H1     A2    A3    A4    A5    A6    A7
lipid_scaffolds[(lipid_type, params)]["x"] = (   0,     0,    0,    0,    0,    0,    0)
lipid_scaffolds[(lipid_type, params)]["y"] = (   0,     0,    0,    0,    0,    0,    0)
lipid_scaffolds[(lipid_type, params)]["z"] = (   6,     5,    4,    3,    2,    1,    0)
lipid_scaffolds[(lipid_type, params)]["charges"] = ("COO", -1)
lipid_scaffolds[(lipid_type, params)]["bd"] = (0.25, 0.25, 0.3)
lipid_scaffolds[(lipid_type, params)]["lipids"] = {
    ("OA", "beads"):  "COO   C1A   D2A   C3A   C4A    -     - ",
    ("LA", "beads"):  "COO   C1A   D2A   D3A   C4A    -     - ",
    ("LNA", "beads"): "COO   C1A   D2A   D3A   D4A    -     - ",
    
    ("MA", "beads"):  "COO   C1A   C2A   C3A    -     -     - ",

    ### PA and SA have identical structures but different topologies
    ("PA", "beads"):  "COO   C1A   C2A   C3A   C4A    -     - ",
    ("SA", "beads"):  "COO   C1A   C2A   C3A   C4A    -     - ",
}


### ### ### ### ### ###################### ### ### ### ### ###
### ### ### ### ### ### Monoglycerides ### ### ### ### ### ###
### ### ### ### ### ###################### ### ### ### ### ###

### Monoglycerides
#
# H1-L2-A3-A4-A5-A6-A7-A8

lipid_type, params = "monoglycerides_LTF_release", "LTF"
lipid_scaffolds[(lipid_type, params)] = {}   #  H1    L2    A3    A4    A5    A6    A7    A8
lipid_scaffolds[(lipid_type, params)]["x"] = (   0,    0,    0,    0,    0,    0,    0,    0)
lipid_scaffolds[(lipid_type, params)]["y"] = (   0,    0,    0,    0,    0,    0,    0,    0)
lipid_scaffolds[(lipid_type, params)]["z"] = (   7,    6,    5,    4,    3,    2,    1,    0)
lipid_scaffolds[(lipid_type, params)]["bd"] = (0.25, 0.25, 0.3)
lipid_scaffolds[(lipid_type, params)]["lipids"] = {
    ("MO",  "beads"): "DOH   GL1   C1A   D2A   C3A   C4A    -     - ",
    ("ML",  "beads"): "DOH   GL1   C1A   D2A   D3A   C4A    -     - ",
    ("MLN", "beads"): "DOH   GL1   C1A   D2A   D3A   D4A    -     - ",

    ### MS and MP have identical structures but different topologies
    ("MS",  "beads"): "DOH   GL1   C1A   C2A   C3A   C4A    -     - ", 
    ("MP",  "beads"): "DOH   GL1   C1A   C2A   C3A   C4A    -     - ",

    ### MM and MLA have identical structures but different topologies
    ("MM",  "beads"): "DOH   GL1   C1A   C2A   C3A    -     -     - ",
    ("MLA", "beads"): "DOH   GL1   C1A   C2A   C3A    -     -     - ",
}

### ### ### ### ### #################### ### ### ### ### ###
### ### ### ### ### ### Diglycerides ### ### ### ### ### ###
### ### ### ### ### #################### ### ### ### ### ###

### Diglycerides are included as a headgroup with the phospholipids

### ### ### ### ### ##################### ### ### ### ### ###
### ### ### ### ### ### Triglycerides ### ### ### ### ### ###
### ### ### ### ### ##################### ### ### ### ### ###

### Triglycerides
#
# L01-A02-A03-A04-A05-A06-A07
#  |
# L08-B09-B10-B11-B12-B13-B14
#  |                    
# L15-C16-C17-C18-C19-C20-C21

lipid_type, params = "triglycerides_LTF_release", "LTF"
lipid_scaffolds[(lipid_type, params)] = {}   #  L01   A02   A03   A04   A05   A06   A07   L08   B09   B10   B11   B12   B13   B14   L15   C16   C17   C18   C19   C20   C21
lipid_scaffolds[(lipid_type, params)]["x"] = (    0,    0,    0,    0,    0,    0,    0,    1,    1,    1,    1,    1,    1,    1,    2,    2,    2,    2,    2,    2,    2)
lipid_scaffolds[(lipid_type, params)]["y"] = (    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0)
lipid_scaffolds[(lipid_type, params)]["z"] = (    6,    5,    4,    3,    2,    1,    0,    6,    5,    4,    3,    2,    1,    0,    6,    5,    4,    3,    2,    1,    0)
lipid_scaffolds[(lipid_type, params)]["bd"] = (0.25, 0.25, 0.3)
lipid_scaffolds[(lipid_type, params)]["lipids"] = {
    
        # L01   A02   A03   A04   A05   A06   A07
        # L08   B09   B10   B11   B12   B13   B14
        # L15   C16   C17   C18   C19   C20   C21

    ("TO", "beads"):  " ".join([
        " GL1   C1A   D2A   C3A   C4A    -     -  ",
        " GL2   C1B   D2B   C3B   C4B    -     -  ",
        " GL3   C1C   D2C   C3C   C4C    -     -  ",
    ]),

    ("TL", "beads"):  " ".join([
        " GL1   C1A   D2A   D3A   C4A    -     -  ",
        " GL2   C1B   D2B   D3B   C4B    -     -  ",
        " GL3   C1C   D2C   D3C   C4C    -     -  ",
    ]),

    ("TLN", "beads"):  " ".join([
        " GL1   C1A   D2A   D3A   D4A    -     -  ",
        " GL2   C1B   D2B   D3B   D4B    -     -  ",
        " GL3   C1C   D2C   D3C   D4C    -     -  ",
    ]),

    ### TS and TP have identical structures but different topologies
    ("TS", "beads"):  " ".join([
        " GL1   C1A   C2A   C3A   C4A    -     -  ",
        " GL2   C1B   C2B   C3B   C4B    -     -  ",
        " GL3   C1C   C2C   C3C   C4C    -     -  ",
    ]),
    ("TP", "beads"):  " ".join([
        " GL1   C1A   C2A   C3A   C4A    -     -  ",
        " GL2   C1B   C2B   C3B   C4B    -     -  ",
        " GL3   C1C   C2C   C3C   C4C    -     -  ",
    ]),

    ### TM and TLA have identical structures but different topologies
    ("TM", "beads"):  " ".join([
        " GL1   C1A   C2A   C3A    -     -     -  ",
        " GL2   C1B   C2B   C3B    -     -     -  ",
        " GL3   C1C   C2C   C3C    -     -     -  ",
    ]),
    ("TLA", "beads"):  " ".join([
        " GL1   C1A   C2A   C3A    -     -     -  ",
        " GL2   C1B   C2B   C3B    -     -     -  ",
        " GL3   C1C   C2C   C3C    -     -     -  ",
    ]),

}


### ### ### ### ### ########################################### ### ### ### ### ### 
### ### ### ### ### ### Sphingolipids (including ceramides) ### ### ### ### ### ###
### ### ### ### ### ########################################### ### ### ### ### ### 

### Some tails have the same beads. The difference is in the topology (regular vs small beads).
sphingolipid_tails = {
    ### No double bonds
    "U":  "T1A   C2A   C3A   C4A    -     -   C1B   C2B   C3B    -     -     - ",
    "M":  "T1A   C2A   C3A   C4A    -     -   C1B   C2B   C3B    -     -     - ",

    "P":  "T1A   C2A   C3A   C4A    -     -   C1B   C2B   C3B   C4B    -     - ",
    "S":  "T1A   C2A   C3A   C4A    -     -   C1B   C2B   C3B   C4B    -     - ",

    "K":  "T1A   C2A   C3A   C4A    -     -   C1B   C2B   C3B   C4B   C5B    - ",
    "B":  "T1A   C2A   C3A   C4A    -     -   C1B   C2B   C3B   C4B   C5B    - ",

    "X":  "T1A   C2A   C3A   C4A    -     -   C1B   C2B   C3B   C4B   C5B   C6B",
    "C":  "T1A   C2A   C3A   C4A    -     -   C1B   C2B   C3B   C4B   C5B    - ",

    ### One double bond on fatty acid (FA), no double bonds on long chain base (LCB)
    "O":  "T1A   C2A   C3A   C4A    -     -   C1B   D2B   C3B   C4B    -     - ",

    "N":  "T1A   C2A   C3A   C4A    -     -   C1B   C2B   C3B   D4B   C5B   C6B",
}

### Sphingolipids with simple head structures
#
#                     B15-B16-B17-B18-B19-B20
#                    /
#    H02   H05    L08
#   /  |   |  \    |
# H01-H03-H04-H06-L07-A09-A10-A11-A12-A13-A14

sphingolipid_heads_simple = {
           #   H01   H02   H03   H04   H05   H06
    "SM":  "    -     -     -    NC3    -    PO4",
    "CER": "    -     -     -     -     -    COH",
}

lipid_type, params = "sphingolipids_LTF_release", "LTF"
lipid_scaffolds[(lipid_type, params)] = {}   #  H01   H02   H03   H04   H05   H06   L07   L08   A09   A10   A11   A12   A13   A14   B15   B16   B17   B18   B19   B20
lipid_scaffolds[(lipid_type, params)]["x"] = (    0,   .5,    0,    0,   .5,    0,    0,   .5,    0,    0,    0,    0,    0,    0,    1,    1,    1,    1,    1,    1)
lipid_scaffolds[(lipid_type, params)]["y"] = (    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0)
lipid_scaffolds[(lipid_type, params)]["z"] = (   10,    9,    9,    8,    8,    7,    6,    6,    5,    4,    3,    2,    1,    0,    5,    4,    3,    2,    1,    0)
lipid_scaffolds[(lipid_type, params)]["charges"] = (("NC3", 1), ("PO4", -1))
lipid_scaffolds[(lipid_type, params)]["bd"] = (0.25, 0.25, 0.3)
lipid_scaffolds[(lipid_type, params)]["lipids"] = {}

sphingolipid_linkers = {
    "SM": "OH1   AM2",
}
sphingolipid_LTF_release_lipids = {
    (tail_key + head_key, "beads"): "   ".join([head_beads, link_beads, tail_beads])
    for head_key, head_beads in sphingolipid_heads_simple.items()
    for link_key, link_beads in sphingolipid_linkers.items()
    for tail_key, tail_beads in sphingolipid_tails.items()
}

lipid_scaffolds[(lipid_type, params)]["lipids"].update(sphingolipid_LTF_release_lipids)

### ### ### ### ### #################### ### ### ### ### ### 
### ### ### ### ### ### Cardiolipins ### ### ### ### ### ###
### ### ### ### ### #################### ### ### ### ### ### 

### Cardiolipins
#
# IDXX: ID indicates part of lipid. XX indicates bead number
# LC = Linker "core"
# L1/L2 = Linker 1/2
# A1/B1 = Tail A/B on linker 1
# A2/B2 = Tail A/B on linker 2
#
#              LC01 # Linker "core"
#             /    \
#          L102     L217 # Linker PO4 beads
#            |       | 
#    L104--L103     L218--L219 # Linker GL beads
#      |     |       |     | 
#    B111  A105     A220  B226 # Tails
#      |     |        |    | 
#    B112  A106     A221  B227
#      |     |        |    | 
#    B113  A107     A222  B228
#      |     |        |    | 
#    B114  A108     A223  B229
#      |     |        |    | 
#    B115  A109     A224  B230
#      |     |        |    | 
#    B116  A110     A225  B231

lipid_type, params = "cardiolipins_LTF_release", "LTF"
lipid_scaffolds[(lipid_type, params)] = {}   # LC01  L102  L103  L104  A105  A106  A107  A108  A109  A110  B111  B112  B113  B114  B115  B116  L217  L218  L219  A220  A221  A222  A223  A224  A225  B226  B227  B228  B229  B230  B231
lipid_scaffolds[(lipid_type, params)]["x"] = (    0,   -1,   -1,   -2,   -1,   -1,   -1,   -1,   -1,   -1,   -2,   -2,   -2,   -2,   -2,   -2,    1,    1,    2,    1,    1,    1,    1,    1,    1,    2,    2,    2,    2,    2,    2)
lipid_scaffolds[(lipid_type, params)]["y"] = (    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0)
lipid_scaffolds[(lipid_type, params)]["z"] = (    8,    7,    6,    6,    5,    4,    3,    2,    1,    0,    5,    4,    3,    2,    1,    0,    7,    6,    6,    5,    4,    3,    2,    1,    0,    5,    4,    3,    2,    1,    0,)
lipid_scaffolds[(lipid_type, params)]["charges"] = (("PO41", -1), ("PO42", -1))
lipid_scaffolds[(lipid_type, params)]["bd"] = (0.25, 0.25, 0.3)
lipid_scaffolds[(lipid_type, params)]["lipids"] = {
                    #    LC01  L102  L103  L104  A105  A106  A107  A108  A109  A110  B111  B112  B113  B114  B115  B116  L217  L218  L219  A220  A221  A222  A223  A224  A225  B226  B227  B228  B229  B230  B231
    ("TMCL",  "beads"): " GLC  PO41  GL11  GL21  C1A1  C2A1  C3A1    -     -     -   C1B1  C2B1  C3B1    -     -     -   PO42  GL12  GL22  C1A2  C2A2  C3A2    -     -     -   C1B2  C2B2  C3C2    -     -     - ",
    ("TOCL",  "beads"): " GLC  PO41  GL11  GL21  C1A1  C2A1  C3A1  C4A1    -     -   C1B1  C2B1  C3B1  C4B1    -     -   PO42  GL12  GL22  C1A2  C2A2  C3A2  C4A2    -     -   C1B2  C2B2  C3C2  C4B2    -     - ",
}


### ### ### ### ### ########################################### ### ### ### ### ### 
### ### ### ### ### ### Bis(monoacylglycero)phosphate (BMP) (Both 2,2' and 3,3') ### ### ### ### ### ###
### ### ### ### ### ########################################### ### ### ### ### ### 

### BMP
#
#           H01
#          /   \
# L103-L102     L211-L212 # L102/L211 are OH1/OH2 for 3,3' and GL1/GL2 for 2,2'. L103/L212 are not used for 3,3' and OH1/OH2 for 2,2'.
#        |        | 
#      E104     E213 # E104/E213 are GL1/GL2 for 3,3' and C1A/C1B for 2,2'.
#        |        | 
#       A05      B13
#        |        | 
#       A06      B14
#        |        | 
#       A07      B15
#        |        | 
#       A08      B16
#        |        | 
#       A09      B17
#        |        | 
#       A10      B18

lipid_type, params = "BMP_LTF_release", "LTF"
lipid_scaffolds[(lipid_type, params)] = {}   #  H01  L102  L103  E104   A05   A06   A07   A08   A09   A10  L211  L212  E213   B13   B14   B15   B16   B17   B18
lipid_scaffolds[(lipid_type, params)]["x"] = (    0,   -1,   -2,   -1,   -1,   -1,   -1,   -1,   -1,   -1,    1,    2,    1,    1,    1,    1,    1,    1,    1)
lipid_scaffolds[(lipid_type, params)]["y"] = (    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0)
lipid_scaffolds[(lipid_type, params)]["z"] = (    8,    7,    7,    6,    5,    4,    3,    2,    1,    0,    7,    7,    6,    5,    4,    3,    2,    1,    0)
lipid_scaffolds[(lipid_type, params)]["charges"] = ("PO4", -1)
lipid_scaffolds[(lipid_type, params)]["bd"] = (0.25, 0.25, 0.3)
lipid_scaffolds[(lipid_type, params)]["lipids"] = {
                    #     H01  L102  L103  E104   A05   A06   A07   A08   A09   A10  L211  L212  E213   B13   B14   B15   B16   B17   B18
    ("DO2B",  "beads"): " PO4   OH1    -    GL1   C1A   D2A   C3A   C4A    -     -    OH2    -    GL2   C1B   D2B   C3B   C4B    -     - ",
    ("DO3B",  "beads"): " PO4   GL1   OH1   C1A   D2A   C3A   C4A    -     -     -    GL2   OH2   C1B   D2B   C3B   C4B    -     -     - ",
}


