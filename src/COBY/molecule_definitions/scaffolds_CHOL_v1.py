lipid_scaffolds = {}

### Sterols # These are from the cholesterol paper (https://doi.org/10.1021/acs.jctc.3c00547)
### The cholesterol entry is duplicated so that it is present in both the Martini 3 release and LTF parameter libraries, since it is compatible with both.

lipid_type, params = "sterol", "M3_2021"
lipid_scaffolds[(lipid_type, params)] = {}
lipid_scaffolds[(lipid_type, params)]["x"] = (       0,   1,   0,   0,  1,   0, 0.5, 0.5,   0,  0)
lipid_scaffolds[(lipid_type, params)]["y"] = (       0,   0,   0,   0,  0,   0,   0,   0,   0,  0)
lipid_scaffolds[(lipid_type, params)]["z"] = (     5.3, 4.5, 3.9, 3.3,  3, 2.6, 4.5, 2.6, 1.4,  0)
lipid_scaffolds[(lipid_type, params)]["bd"] = (0.25, 0.25, 0.3)
lipid_scaffolds[(lipid_type, params)]["tags"] = ("martini3", "sterol", "head_ROH")
lipid_scaffolds[(lipid_type, params)]["lipids"] = {
    ("CHOL", "beads"): " ROH  R1  R2  R3  R4  -  R5  R6  C1  C2 ",
    ("CHOL", "tags"): ("core_CHOL",),
}

lipid_type, params = "sterol", "M3_2025"
lipid_scaffolds[(lipid_type, params)] = {}
lipid_scaffolds[(lipid_type, params)]["x"] = (       0,   1,   0,   0,  1,   0, 0.5, 0.5,   0,  0)
lipid_scaffolds[(lipid_type, params)]["y"] = (       0,   0,   0,   0,  0,   0,   0,   0,   0,  0)
lipid_scaffolds[(lipid_type, params)]["z"] = (     5.3, 4.5, 3.9, 3.3,  3, 2.6, 4.5, 2.6, 1.4,  0)
lipid_scaffolds[(lipid_type, params)]["bd"] = (0.25, 0.25, 0.3)
lipid_scaffolds[(lipid_type, params)]["tags"] = ("martini3", "sterol", "head_ROH")
lipid_scaffolds[(lipid_type, params)]["lipids"] = {
    ("CHOL", "beads"): " ROH  R1  R2  R3  R4  -  R5  R6  C1  C2 ",
    ("CHOL", "tags"): ("core_CHOL",),
}
