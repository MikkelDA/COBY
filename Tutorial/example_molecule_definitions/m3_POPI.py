### Martini 2 lipids

lipid_scaffolds = {}

lipid_type, params = "INOSITOLLIPIDS", "IMPORTED"
lipid_scaffolds[(lipid_type, params)] = {}   #  H01   H02   H03   V04   H05   H06   H07   H08   L09   L10   A11   A12   A13   A14   A15   A16   B17   B18   B19   B20   B21   B22
lipid_scaffolds[(lipid_type, params)]["x"] = (  0.5,  0.8, -0.3, 0.25,    0,    1,   .5, -0.3,    0,   .5,    0,    0,    0,    0,    0,    0,    1,    1,    1,    1,    1,    1)
lipid_scaffolds[(lipid_type, params)]["y"] = (    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0)
lipid_scaffolds[(lipid_type, params)]["z"] = (    8,  9.1,  8.8,  8.5,    7,   10,   10,   10,    6,    6,    5,    4,    3,    2,    1,    0,    5,    4,    3,    2,    1,    0)
lipid_scaffolds[(lipid_type, params)]["charges"] = ("PO4", -1)
lipid_scaffolds[(lipid_type, params)]["bd"] = (0.25, 0.25, 0.3)
lipid_scaffolds[(lipid_type, params)]["lipids"] = {
    ("POPI", "beads"): "    C1    C2    C3    C4   PO4     -     -     -   GL1   GL2   C1A   D2A   C3A   C4A     -     -   C1B   C2B   C3B   C4B    -     - ",
}


