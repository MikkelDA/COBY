### Martini 2 lipids

lipid_scaffolds = {}

############################################### PHOSPHOLIPIDS and LYSEPHOSPHOLIPIDS
phospholipid_tails = {
    "DP": "C1A   C2A   C3A   C4A     -     -   C1B   C2B   C3B   C4B     -     - ",
    "PO": "C1A   D2A   C3A   C4A     -     -   C1B   C2B   C3B   C4B     -     - ",
    "PA": "D1A   D2A   D3A   D4A   C5A     -   C1B   C2B   C3B   C4B     -     - ",
    "DO": "C1A   D2A   C3A   C4A     -     -   C1B   D2B   C3B   C4B     -     - ",
    "PU": "D1A   D2A   D3A   D4A   D5A     -   C1B   C2B   C3B   C4B     -     - ",
    "PF": "C1A   D2A   D3A   D4A     -     -   C1B   C2B   C3B   C4B     -     - ",
    "OI": "C1A   D2A   D3A   C4A     -     -   C1B   D2B   C3B   C4B     -     - ",
    "OU": "D1A   D2A   D3A   D4A   D5A     -   C1B   D2B   C3B   C4B     -     - ",
    "OA": "D1A   D2A   D3A   D4A   C5A     -   C1B   D2B   C3B   C4B     -     - ",
}
lysephospholipid_tails = {
    "P": "C1A   C2A   C3A   C4A     -     -     -     -     -     -     -     - ",
    "I": "C1A   D2A   D3A   C4A     -     -     -     -     -     -     -     - ",
}

### Phospholipids with simple head structures
#
#                     B15-B16-B17-B18-B19-B20
#                    /
#    H02   H05    L08
#   /  |   |  \    |
# H01-H03-H04-H06-L07-A09-A10-A11-A12-A13-A14

phospholipid_heads_simple = {
    "PA": "    -     -     -     -     -   PO4   GL1   GL2",
    "PC": "    -     -     -   NC3     -   PO4   GL1   GL2",
    "PE": "    -     -     -   NH3     -   PO4   GL1   GL2",
    "PG": "    -     -     -   GL0     -   PO4   GL1   GL2",
    "PS": "    -     -     -   CNO     -   PO4   GL1   GL2",
    "DG": "    -     -     -    -      -    -    GL1   GL2",
}

lipid_type, params = "diacyl_glycerols_generated", "IngolfssonMembranes"
lipid_scaffolds[(lipid_type, params)] = {}   #  H01   H02   H03   H04   H05   H06   L07   L08   A09   A10   A11   A12   A13   A14   B15   B16   B17   B18   B19   B20
lipid_scaffolds[(lipid_type, params)]["x"] = (    0,   .5,    0,    0,   .5,    0,    0,   .5,    0,    0,    0,    0,    0,    0,    1,    1,    1,    1,    1,    1)
lipid_scaffolds[(lipid_type, params)]["y"] = (    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0)
lipid_scaffolds[(lipid_type, params)]["z"] = (   10,    9,    9,    8,    8,    7,    6,    6,    5,    4,    3,    2,    1,    0,    5,    4,    3,    2,    1,    0)
lipid_scaffolds[(lipid_type, params)]["charges"] = (("NC3", 1), ("NH3", 1), ("PO4", -1))
lipid_scaffolds[(lipid_type, params)]["bd"] = (0.25, 0.25, 0.3)
lipid_scaffolds[(lipid_type, params)]["lipids"] = {}

diacyl_glycerols_generated = {
    (tail_key + head_key, "beads"): "   ".join([head_beads, tail_beads])
    for head_key, head_beads in phospholipid_heads_simple.items()
    for tail_key, tail_beads in {**phospholipid_tails, **lysephospholipid_tails}.items()
}

lipid_scaffolds[(lipid_type, params)]["lipids"].update(diacyl_glycerols_generated)

############################################### INOSITOLLIPIDS
# 1,2,3 - is the inositol and 4 is the phosphate that links to the tail part.
#
#  5
#   \
#  6-2-1-4-8--10-11-12-13-14-15
#    |/    |
#  7-3     9--16-17-18-19-20-21 
lipid_type, params = "inositollipids", "IngolfssonMembranes"
lipid_scaffolds[(lipid_type, params)] = {}   #  H01   H02   H03   H04   H05   H06   L07   L08   A09   A10   A11   A12   A13   A14   B15   B16   B17   B18   B19   B20
lipid_scaffolds[(lipid_type, params)]["x"] = (   .5,  .5,   0,   0,   1, .5,  0,  0,   .5,   0,   0,   0,   0,   0,   0,   1,   1,   1,   1,   1,   1)
lipid_scaffolds[(lipid_type, params)]["y"] = (    0,   0,   0,   0,   0,  0,  0,  0,    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0)
lipid_scaffolds[(lipid_type, params)]["z"] = (    8,   9,   9,   7,  10, 10, 10,  6,    6,   5,   4,   3,   2,   1,   0,   5,   4,   3,   2,   1,   0)
# lipid_scaffolds[(lipid_type, params)]["charges"] = (("NC3", 1), ("NH3", 1), ("PO4", -1))
lipid_scaffolds[(lipid_type, params)]["bd"] = (0.25, 0.25, 0.3)
lipid_scaffolds[(lipid_type, params)]["lipids"] = {      # 1     2    3    4    5   6   7   8    9    10    11    12    13    14   15    16    17    18    19   20 
    # ("DPPI", "beads"): (" C1   C2   C3   PO4   -   -   -  GL1  GL2  C1A  C2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
    # ("DOPI", "beads"): (" C1   C2   C3   PO4   -   -   -  GL1  GL2  C1A  C2A  D3A  C4A  C5A   -   C1B  C2B  D3B  C4B  C5B   - "),
    # ("PI"  , "beads"): (" C1   C2   C3    CP   -   -   -  GL1  GL2  C1A  C2A  C3A  C4A   -    -   CU1  CU2  CU3  CU4  CU5   - "),
    # ("PI34", "beads"): (" C1   C2   C3    CP PO1 PO2   -  GL1  GL2  C1A  C2A  C3A  C4A   -    -   CU1  CU2  CU3  CU4  CU5   - "),
    ("POPI", "beads"): (" C1   C2   C3   PO4   -   -   -  GL1  GL2  C1A  D2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
    ("PIPI", "beads"): (" C1   C2   C3   PO4   -   -   -  GL1  GL2  C1A  D2A  D3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
    ("PAPI", "beads"): (" C1   C2   C3   PO4   -   -   -  GL1  GL2  D1A  D2A  D3A  D4A  C5A   -   C1B  C2B  C3B  C4B   -    - "),
    ("PUPI", "beads"): (" C1   C2   C3   PO4   -   -   -  GL1  GL2  D1A  D2A  D3A  D4A  D5A   -   C1B  C2B  C3B  C4B   -    - "),
    ("POP1", "beads"): (" C1   C2   C3   PO4  P1   -   -  GL1  GL2  C1A  D2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
    ("PAP1", "beads"): (" C1   C2   C3   PO4  P1   -   -  GL1  GL2  D1A  D2A  D3A  D4A  C5A   -   C1B  C2B  C3B  C4B   -    - "),
    ("POP2", "beads"): (" C1   C2   C3   PO4  P1  P2   -  GL1  GL2  C1A  D2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
    ("PAP2", "beads"): (" C1   C2   C3   PO4  P1  P2   -  GL1  GL2  D1A  D2A  D3A  D4A  C5A   -   C1B  C2B  C3B  C4B   -    - "),
    ("POP3", "beads"): (" C1   C2   C3   PO4  P1  P2  P3  GL1  GL2  C1A  D2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
    ("PAP3", "beads"): (" C1   C2   C3   PO4  P1  P2  P3  GL1  GL2  D1A  D2A  D3A  D4A  C5A   -   C1B  C2B  C3B  C4B   -    - "),
}


############################################### CERAMIDES and SPHINGOMYELINS
sphingomyelin_tails = {
    "DP": "T1A   C2A   C3A    -      -     -   C1B   C2B   C3B   C4B    -     - ",
    "PN": "T1A   C2A   C3A    -      -     -   C1B   C2B   C3B   D4B   C5B   C5B",
    "PO": "T1A   C2A   C3A    -      -     -   C1B   D2B   C3B   C4B    -     - ",
    "PB": "T1A   C2A   C3A    -      -     -   C1B   C2B   C3B   C4B   C5B    - ",
    "DB": "T1A   C2A   C3A   C4A     -     -   C1B   C2B   C3B   C4B   C5B    - ",
}

### Phospholipids with simple head structures
#
#                     B15-B16-B17-B18-B19-B20
#                    /
#    H02   H05    L08
#   /  |   |  \    |
# H01-H03-H04-H06-L07-A09-A10-A11-A12-A13-A14

sphingomyelin_heads = {
    "SM": "    -     -     -   NC3    -   PO4   AM1   AM2",
    "CE": "    -     -     -    -     -    -    AM1   AM2",
}

lipid_type, params = "sphingomyelins_generated", "IngolfssonMembranes"
lipid_scaffolds[(lipid_type, params)] = {}   #  H01   H02   H03   H04   H05   H06   L07   L08   A09   A10   A11   A12   A13   A14   B15   B16   B17   B18   B19   B20
lipid_scaffolds[(lipid_type, params)]["x"] = (    0,   .5,    0,    0,   .5,    0,    0,   .5,    0,    0,    0,    0,    0,    0,    1,    1,    1,    1,    1,    1)
lipid_scaffolds[(lipid_type, params)]["y"] = (    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0)
lipid_scaffolds[(lipid_type, params)]["z"] = (   10,    9,    9,    8,    8,    7,    6,    6,    5,    4,    3,    2,    1,    0,    5,    4,    3,    2,    1,    0)
lipid_scaffolds[(lipid_type, params)]["charges"] = (("NC3", 1), ("NH3", 1), ("PO4", -1))
lipid_scaffolds[(lipid_type, params)]["bd"] = (0.25, 0.25, 0.3)
lipid_scaffolds[(lipid_type, params)]["lipids"] = {}

sphingomyelins_generated = {
    (tail_key + head_key, "beads"): "   ".join([head_beads, tail_beads])
    for head_key, head_beads in sphingomyelin_heads.items()
    for tail_key, tail_beads in sphingomyelin_tails.items()
}

lipid_scaffolds[(lipid_type, params)]["lipids"].update(sphingomyelins_generated)

#Prototopology for longer and branched glycosil and ceramide based glycolipids
#
#     17-15-14-16
#         |/
#        13
#         |
# 12-10-9-7-6-4-3-1--18--26-27-28-29-30-31
#  |/   |/  |/  |/    |
#  11   8   5   2    19--20-21-22-23-24-25 

lipid_type, params = "glycolipids", "IngolfssonMembranes"
lipid_scaffolds[(lipid_type, params)] = {}   #  H01   H02   H03   H04   H05   H06   L07   L08   A09   A10   A11   A12   A13   A14   B15   B16   B17   B18   B19   B20
lipid_scaffolds[(lipid_type, params)]["x"] = (    0,  .5,   0,   0,  .5,  0,  0, .5,  0,    0,   .5,    0,    0,    0,   0,    0,    0,    0,   .5,   0,   0,   0,   0,   0,   0,   1,   1,   1,   1,   1,   1)
lipid_scaffolds[(lipid_type, params)]["y"] = (    0,   0,   0,   0,   0,  0,  0,  0,  0,    0,    0,    0,   .5,    1,   1,    1,    1,    0,    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0)
lipid_scaffolds[(lipid_type, params)]["z"] = (    7,   8,   8,   9,  10, 10, 11, 12, 12,   13,   14,   14,   11,   10,  11,    9,   12,    6,    6,   5,   4,   3,   2,   1,   0,   5,   4,   3,   2,   1,   0)
# lipid_scaffolds[(lipid_type, params)]["charges"] = (("NC3", 1), ("NH3", 1), ("PO4", -1))
lipid_scaffolds[(lipid_type, params)]["bd"] = (0.25, 0.25, 0.3)
lipid_scaffolds[(lipid_type, params)]["lipids"] = {      # 1     2    3    4    5   6   7   8   9    10    11    12    13    14   15    16    17    18    19   20   21   22   23   24   25   26   27   28   29   30   31
    # ("GM1" , "beads"): ("GM1  GM2  GM3  GM4  GM5 GM6 GM7 GM8 GM9  GM10  GM11  GM12  GM13  GM14 GM15  GM16  GM17  GM18  GM19 GM20 GM21 GM22 GM23 GM24   -  GM25 GM26 GM27 GM28   -    - "),
    # ("DGDG", "beads"): ("GB2  GB3  GB1  GA1  GA2 GA3   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
    # ("MGDG", "beads"): (" C1   C2   C3    -    -   -   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
    # ("SQDG", "beads"): (" S1   C1   C2   C3    -   -   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
    # ("CER" , "beads"): ("  -    -    -    -    -   -   -   -   -     -     -     -     -     -    -     -     -   AM1   AM2  C1A  C2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
    # ("GCER", "beads"): (" C1   C2   C3    -    -   -   -   -   -     -     -     -     -     -    -     -     -   AM1   AM2  C1A  C2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
    ("DPGS", "beads"): (" C1   C2   C3    -    -   -   -   -   -     -     -     -     -     -    -     -     -   AM1   AM2  T1A  C2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
    ("DBGS", "beads"): (" C1   C2   C3    -    -   -   -   -   -     -     -     -     -     -    -     -     -   AM1   AM2  T1A  C2A  C3A  C4A  C5A   -   C1B  C2B  C3B  C4B  C5B   - "),
    ("POGS", "beads"): (" C1   C2   C3    -    -   -   -   -   -     -     -     -     -     -    -     -     -   AM1   AM2  T1A  C2A  C3A  C4A   -    -   C1B  D2B  C3B  C4B   -    - "),
    ("PNGS", "beads"): (" C1   C2   C3    -    -   -   -   -   -     -     -     -     -     -    -     -     -   AM1   AM2  T1A  C2A  C3A  C4A   -    -   C1B  C2B  C3B  D4B  C5B  C6B"),
    # ("PXSU", "beads"): (" S1   C1   C2   C3    -   -   -   -   -     -     -     -     -     -    -     -     -   AM1   AM2  T1A  C2A  C3A   -    -    -   C1B  C2B  C3B  C4B  C5B  C6B"),
    # ("PNSU", "beads"): (" S1   C1   C2   C3    -   -   -   -   -     -     -     -     -     -    -     -     -   AM1   AM2  T1A  C2A  C3A   -    -    -   C1B  C2B  C3B  D4B  C5B  C6B"),
# HII add GM1(s) with dififrent tails  -  2013.11.10 + extended matrix by x2 columns
    ("DPG1", "beads") : ("GM1  GM2  GM3  GM4  GM5 GM6 GM7 GM8 GM9  GM10  GM11  GM12  GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A   -    -    -   C1B  C2B  C3B  C4B   -    - "),
    ("DBG1", "beads") : ("GM1  GM2  GM3  GM4  GM5 GM6 GM7 GM8 GM9  GM10  GM11  GM12  GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B  C5B   - "),
    ("POG1", "beads") : ("GM1  GM2  GM3  GM4  GM5 GM6 GM7 GM8 GM9  GM10  GM11  GM12  GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A   -    -    -   C1B  D2B  C3B  C4B   -    - "),
    ("PNG1", "beads") : ("GM1  GM2  GM3  GM4  GM5 GM6 GM7 GM8 GM9  GM10  GM11  GM12  GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A   -    -    -   C1B  C2B  C3B  D4B  C5B  C6B"),
    # ("DXG1", "beads") : ("GM1  GM2  GM3  GM4  GM5 GM6 GM7 GM8 GM9  GM10  GM11  GM12  GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A  C4A  C5A   -   C1B  C2B  C3B  C4B  C5B  C6B"),
    # ("XNG1", "beads") : ("GM1  GM2  GM3  GM4  GM5 GM6 GM7 GM8 GM9  GM10  GM11  GM12  GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A  C4A  C5A   -   C1B  C2B  C3B  D4B  C5B  C6B"),
    ("DPG3", "beads") : ("GM1  GM2  GM3  GM4  GM5 GM6  -   -   -    -     -     -    GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A   -    -    -   C1B  C2B  C3B  C4B   -    - "),
    ("DBG3", "beads") : ("GM1  GM2  GM3  GM4  GM5 GM6  -   -   -    -     -     -    GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B  C5B   - "),
    ("POG3", "beads") : ("GM1  GM2  GM3  GM4  GM5 GM6  -   -   -    -     -     -    GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A   -    -    -   C1B  D2B  C3B  C4B   -    - "),
    ("PNG3", "beads") : ("GM1  GM2  GM3  GM4  GM5 GM6  -   -   -    -     -     -    GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A   -    -    -   C1B  C2B  C3B  D4B  C5B  C6B"),
    # ("DXG3", "beads") : ("GM1  GM2  GM3  GM4  GM5 GM6  -   -   -    -     -     -    GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A  C4A  C5A   -   C1B  C2B  C3B  C4B  C5B  C6B"),
    # ("XNG3", "beads") : ("GM1  GM2  GM3  GM4  GM5 GM6  -   -   -    -     -     -    GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A  C4A  C5A   -   C1B  C2B  C3B  D4B  C5B  C6B"),
}

############################################### CHOLESTEROL
lipid_type, params = "sterols", "IngolfssonMembranes"
lipid_scaffolds[(lipid_type, params)] = {}   #  H01   H02   H03   H04   H05   H06   L07   L08   A09   A10   A11   A12   A13   A14   B15   B16   B17   B18   B19   B20
lipid_scaffolds[(lipid_type, params)]["x"] = (     0,  0,  0,  0,  0, 0,   0,  0,  1,  0,  1,  0,  0,  0,  0,  0,  0,  0)
lipid_scaffolds[(lipid_type, params)]["y"] = (     0,  0,  0,  0,  0, 0,   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0)
lipid_scaffolds[(lipid_type, params)]["z"] = (     0,  0,  0,  0,  0, 0, 5.3,4.5,3.9,3.3, 3 ,2.6,1.4,  0,  0,  0,  0,  0)
lipid_scaffolds[(lipid_type, params)]["bd"] = (0.25, 0.25, 0.3)
lipid_scaffolds[(lipid_type, params)]["lipids"] = {
    ("CHOL", "beads"): (" -   -   -   -   -   -  ROH  R1  R2  R3  R4  R5  C1  C2  -   -   -   - "),
}
