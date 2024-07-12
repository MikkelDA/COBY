lipid_scaffolds = {}

### DISCLAIMER
### The following scaffolds are converted to the COBY format from insane (https://github.com/Tsjerk/Insane).

##############
### LIPIDS ###
##############
### Diacyl glycerols
lipid_type, params = "lipid", "default"
lipid_scaffolds[(lipid_type, params)] = {}
lipid_scaffolds[(lipid_type, params)]["x"] = (    0, .5,  0,  0, .5,  0,  0, .5,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1)
lipid_scaffolds[(lipid_type, params)]["y"] = (    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0)
lipid_scaffolds[(lipid_type, params)]["z"] = (   10,  9,  9,  8,  8,  7,  6,  6,  5,  4,  3,  2,  1,  0,  5,  4,  3,  2,  1,  0)
lipid_scaffolds[(lipid_type, params)]["bd"] = (0.25, 0.25, 0.3)
lipid_scaffolds[(lipid_type, params)]["charges"] = (("NC3", 1), ("NH3", 1), ("PO4", -1))
lipid_scaffolds[(lipid_type, params)]["lipids"] = {      # 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20
## Monoacylglycerol
    ("GMO", "beads"):  " -   -   -   -   -   -  GL1 GL2 C1A C2A D3A C4A C5A  -   -   -   -   -   -   - ",
}

linkers = {
    "GL": "GL1 GL2",
}

SM_linkers = {
    "SM": "AM1 AM2",
}

heads = {
          #   H01   H02   H03   H04   H05   H06
    "PA": "    -     -     -     -     -    PO4",
    "PC": "    -     -     -    NC3    -    PO4",
    "PE": "    -     -     -    NH3    -    PO4",
    "PG": "    -     -     -    GL0    -    PO4",
    "PS": "    -     -     -    CNO    -    PO4",
}

SM_heads = {
          #   H01   H02   H03   H04   H05   H06
    "SM": "    -     -     -    NC3    -    PO4",
}

tailcodes = {
    "C": "C",
    "T": "CC",
    "L": "CCC",
    "M": "CCC", # M used sometimes instead of L
    "P": "CCCC",
    "B": "CCCCC",
    "X": "CCCCCC",
    "Y": "CDC",
    "O": "CDCC",
    "V": "CCDC",
    "G": "CCDCC",
    "N": "CCCDCC",
    "I": "CDDC",
    "F": "CDDD",
    "E": "CCDDC",
    "Q": "CDDDC",
    "A": "DDDDC",
    "U": "DDDDD",
    "R": "DDDDDD",
    "J": "TCCC",
}

SM_tailcodes = {
    "P": "TCC",
    "B": "TCCC",
    "X": "TCCCC",
}

lipids = {}
for head_name, head_string in heads.items():
    for linker_name, linker_string in linkers.items():
        for sn1_name, sn1_code in tailcodes.items():
            for sn2_name, sn2_code in tailcodes.items():
                stringlist = []
                stringlist.extend(head_string.split())
                stringlist.extend(linker_string.split())
                for i, code in enumerate(sn1_code):
                    stringlist.append(code + str(i+1) + "A")
                for i in range(6-len(sn1_code)):
                    stringlist.append(" - ")
                for i, code in enumerate(sn2_code):
                    stringlist.append(code + str(i+1) + "B")
                for i in range(6-len(sn2_code)):
                    stringlist.append(" - ")
                    
                if sn1_name == sn2_name:
                    tailname = "D" + sn1_name
                else:
                    tailname = sn2_name + sn1_name
                name = tailname + head_name
                lipids[(name, "beads")] = " ".join(stringlist)
                
for head_name, head_string in SM_heads.items():
    for linker_name, linker_string in SM_linkers.items():
        for sn1_name, sn1_code in SM_tailcodes.items():
            for sn2_name, sn2_code in tailcodes.items():
                stringlist = []
                stringlist.extend(head_string.split())
                stringlist.extend(linker_string.split())
                for i, code in enumerate(sn1_code):
                    stringlist.append(code + str(i+1) + "A")
                for i in range(6-len(sn1_code)):
                    stringlist.append(" - ")
                for i, code in enumerate(sn2_code):
                    stringlist.append(code + str(i+1) + "B")
                for i in range(6-len(sn2_code)):
                    stringlist.append(" - ")
                    
                if sn1_name == sn2_name:
                    tailname = "D" + sn1_name
                else:
                    tailname = sn2_name + sn1_name
                name = tailname + head_name
                lipids[(name, "beads")] = " ".join(stringlist)

lipid_scaffolds[(lipid_type, params)]["lipids"].update(lipids)



### Prototopology for phosphatidylinositol type lipids 5,6,7 are potentail phosphates (PIP1,PIP2 and PIP3)
# 1,2,3 - is the inositol and 4 is the phosphate that links to the tail part.
#  5
#   \
#  6-2-1-4-8--10-11-12-13-14-15
#    |/    |
#  7-3     9--16-17-18-19-20-21 

lipid_type, params = "INOSITOLLIPIDS", "default"
lipid_scaffolds[(lipid_type, params)] = {}
lipid_scaffolds[(lipid_type, params)]["x"] = (   .5,  .5,   0,   0,   1, .5,  0,  0,   .5,   0,   0,   0,   0,   0,   0,   1,   1,   1,   1,   1,   1)
lipid_scaffolds[(lipid_type, params)]["y"] = (    0,   0,   0,   0,   0,  0,  0,  0,    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0)
lipid_scaffolds[(lipid_type, params)]["z"] = (    8,   9,   9,   7,  10, 10, 10,  6,    6,   5,   4,   3,   2,   1,   0,   5,   4,   3,   2,   1,   0)
lipid_scaffolds[(lipid_type, params)]["bd"] = (0.25, 0.25, 0.3)
lipid_scaffolds[(lipid_type, params)]["charges"] = (("NC3", 1), ("PO4", -1))
lipid_scaffolds[(lipid_type, params)]["lipids"] = {      # 1     2    3    4    5   6   7   8    9    10    11    12    13    14   15    16    17    18    19   20 
    ("DPPI", "beads"): " C1   C2   C3    CP   -   -   -  GL1  GL2  C1A  C2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - ",
    ("POPI", "beads"): " C1   C2   C3    CP   -   -   -  GL1  GL2  C1A  D2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - ",
    ("PIPI", "beads"): " C1   C2   C3    CP   -   -   -  GL1  GL2  C1A  D2A  D3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - ",
    ("PAPI", "beads"): " C1   C2   C3    CP   -   -   -  GL1  GL2  D1A  D2A  D3A  D4A  C5A   -   C1B  C2B  C3B  C4B   -    - ",
    ("PUPI", "beads"): " C1   C2   C3    CP   -   -   -  GL1  GL2  D1A  D2A  D3A  D4A  D5A   -   C1B  C2B  C3B  C4B   -    - ",
    ("POP1", "beads"): " C1   C2   C3    CP  P1   -   -  GL1  GL2  C1A  C2A  D3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - ",
    ("POP2", "beads"): " C1   C2   C3    CP  P1  P2   -  GL1  GL2  C1A  C2A  D3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - ",
    ("POP3", "beads"): " C1   C2   C3    CP  P1  P2  P3  GL1  GL2  C1A  C2A  D3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - ",
}

### Prototopology for longer and branched glycosil and ceramide based glycolipids
#
#     17-15-14-16
#         |/
#        13
#         |
# 12-10-9-7-6-4-3-1--18--20-21-22-23-24
#  |/   |/  |/  |/    |
#  11   8   5   2    19--25-26-27-28-29 
lipid_type, params = "GLYCOLIPIDS", "default"
lipid_scaffolds[(lipid_type, params)] = {}
lipid_scaffolds[(lipid_type, params)]["x"] = (    0,  .5,   0,   0,  .5,  0,  0, .5,  0,    0,   .5,    0,    0,    0,   0,    0,    0,    0,   .5,   0,   0,   0,   0,   0,   1,   1,   1,   1,   1)
lipid_scaffolds[(lipid_type, params)]["y"] = (    0,   0,   0,   0,   0,  0,  0,  0,  0,    0,    0,    0,   .5,    1,   1,    1,    1,    0,    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0)
lipid_scaffolds[(lipid_type, params)]["z"] = (    6,   7,   7,   8,   9,  9, 10, 11, 11,   12,   13,   13,   10,    9,  10,    8,   11,    5,    5,   4,   3,   2,   1,   0,   4,   3,   2,   1,   0)
lipid_scaffolds[(lipid_type, params)]["bd"] = (0.25, 0.25, 0.3)
lipid_scaffolds[(lipid_type, params)]["lipids"] = {
                      # 1     2    3    4    5   6   7   8   9    10    11    12    13    14   15    16    17    18    19   20   21   22   23   24   25   26   27   28   29
    ("DPG1", "beads"): " GM1   GM2  GM3  GM4  GM5 GM6 GM7 GM8 GM9  GM10  GM11  GM12  GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A   -    -    -   C1B  C2B  C3B  C4B   -    - ",
    ("DXG1", "beads"): " GM1   GM2  GM3  GM4  GM5 GM6 GM7 GM8 GM9  GM10  GM11  GM12  GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A  C4A  C5A   -   C1B  C2B  C3B  C4B  C5B  C6B",
    ("PNG1", "beads"): " GM1   GM2  GM3  GM4  GM5 GM6 GM7 GM8 GM9  GM10  GM11  GM12  GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A   -    -    -   C1B  C2B  C3B  D4B  C5B  C6B",
    ("XNG1", "beads"): " GM1   GM2  GM3  GM4  GM5 GM6 GM7 GM8 GM9  GM10  GM11  GM12  GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A  C4A  C5A   -   C1B  C2B  C3B  D4B  C5B  C6B",
    ("DPG3", "beads"): " GM1   GM2  GM3  GM4  GM5 GM6  -   -   -    -     -     -    GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A   -    -    -   C1B  C2B  C3B  C4B   -    - ",
    ("DXG3", "beads"): " GM1   GM2  GM3  GM4  GM5 GM6  -   -   -    -     -     -    GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A  C4A  C5A   -   C1B  C2B  C3B  C4B  C5B  C6B",
    ("PNG3", "beads"): " GM1   GM2  GM3  GM4  GM5 GM6  -   -   -    -     -     -    GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A   -    -    -   C1B  C2B  C3B  D4B  C5B  C6B",
    ("XNG3", "beads"): " GM1   GM2  GM3  GM4  GM5 GM6  -   -   -    -     -     -    GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A  C4A  C5A   -   C1B  C2B  C3B  D4B  C5B  C6B",
    ("DPCE", "beads"): "  -     -    -    -    -   -   -   -   -    -     -     -     -     -    -     -     -     AM1   AM2  T1A  C2A  C3A   -    -   C1B  C2B  C3B  C4B   -    -    - ",
    ("DPGS", "beads"): " C1    C2   C3    -    -   -   -   -   -    -     -     -     -     -    -     -     -     AM1   AM2  T1A  C2A  C3A   -    -   C1B  C2B  C3B  C4B   -    -    - ",
    ("DPMG", "beads"): " C1    C2   C3    -    -   -   -   -   -    -     -     -     -     -    -     -     -     GL1   GL2  C1A  C2A  C3A  C4A   -   C1B  C2B  C3B  C4B   -    -    - ",
    ("DPSG", "beads"): " S1    C1   C2   C3    -   -   -   -   -    -     -     -     -     -    -     -     -     GL1   GL2  C1A  C2A  C3A  C4A   -   C1B  C2B  C3B  C4B   -    -    - ",
    ("DPGG", "beads"): "GB2   GB3  GB1  GA1  GA2 GA3   -   -   -    -     -     -     -     -    -     -     -     GL1   GL2  C1A  C2A  C3A  C4A   -   C1B  C2B  C3B  C4B   -    -    - ",
    ("OPMG", "beads"): " C1    C2   C3    -    -   -   -   -   -    -     -     -     -     -    -     -     -     GL1   GL2  C1A  C2A  C3A  C4A   -   C1B  D2B  C3B  C4B   -    -    - ",
    ("OPSG", "beads"): " S1    C1   C2   C3    -   -   -   -   -    -     -     -     -     -    -     -     -     GL1   GL2  C1A  C2A  C3A  C4A   -   C1B  D2B  C3B  C4B   -    -    - ",
    ("OPGG", "beads"): "GB2   GB3  GB1  GA1  GA2 GA3   -   -   -    -     -     -     -     -    -     -     -     GL1   GL2  C1A  C2A  C3A  C4A   -   C1B  D2B  C3B  C4B   -    -    - ",
    ("FPMG", "beads"): " C1    C2   C3    -    -   -   -   -   -    -     -     -     -     -    -     -     -     GL1   GL2  C1A  C2A  C3A  C4A   -   C1B  D2B  D3B  D4B   -    -    - ",
    ("DFMG", "beads"): " C1    C2   C3    -    -   -   -   -   -    -     -     -     -     -    -     -     -     GL1   GL2  C1A  D2A  D3A  D4A   -   C1B  D2B  D3B  D4B   -    -    - ",
    ("FPSG", "beads"): " S1    C1   C2   C3    -   -   -   -   -    -     -     -     -     -    -     -     -     GL1   GL2  C1A  C2A  C3A  C4A   -   C1B  D2B  D3B  D4B   -    -    - ",
    ("FPGG", "beads"): "GB2   GB3  GB1  GA1  GA2 GA3   -   -   -    -     -     -     -     -    -     -     -     GL1   GL2  C1A  C2A  C3A  C4A   -   C1B  D2B  D3B  D4B   -    -    - ",
    ("DFGG", "beads"): "GB2   GB3  GB1  GA1  GA2 GA3   -   -   -    -     -     -     -     -    -     -     -     GL1   GL2  C1A  D2A  D3A  D4A   -   C1B  D2B  D3B  D4B   -    -    - ",
}

lipid_type, params = "QUINONES", "default"
lipid_scaffolds[(lipid_type, params)] = {}
lipid_scaffolds[(lipid_type, params)]["x"] = (    0,  .5,   0,    0,   0,   0,   0,   0,   0,    0,    0,    0)
lipid_scaffolds[(lipid_type, params)]["y"] = (    0,   0,   0,    0,   0,   0,   0,   0,   0,    0,    0,    0)
lipid_scaffolds[(lipid_type, params)]["z"] = (    6,   7,   7,   5.5,  5,  4.5,  4,  3.5, 2.5,   2,  1.5,    1)
lipid_scaffolds[(lipid_type, params)]["bd"] = (0.25, 0.25, 0.3)
lipid_scaffolds[(lipid_type, params)]["lipids"] = {      # 1     2    3    4    5    6    7    8    9    10    11    12
    ("PLQ", "beads"): (" PLQ3 PLQ2 PLQ1 PLQ4 PLQ5 PLQ6 PLQ7 PLQ8 PLQ9 PLQ10 PLQ11 PLQ12"),
}

# Prototopology for cardiolipins
#  
#       4-11-12-13-14-15-16
#       |
#   2---3--5--6--7--8--9-10
#  / 
# 1
#  \
#   17-18-20-21-22-23-24-25
#       |
#      19-26-27-28-29-30-31
#
lipid_type, params = "CARDIOLIPINS", "default"
lipid_scaffolds[(lipid_type, params)] = {}
lipid_scaffolds[(lipid_type, params)]["x"] = (   0.5,   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1)
lipid_scaffolds[(lipid_type, params)]["y"] = (     1,   0,  0,  1,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,   0,  0,  1,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1)
lipid_scaffolds[(lipid_type, params)]["z"] = (     8,   7,  6,  6,  5,  4,  3,  2,  1,  0,  5,  4,  3,  2,  1,  0,   7,  6,  6,  5,  4,  3,  2,  1,  0,  5,  4,  3,  2,  1,  0)
lipid_scaffolds[(lipid_type, params)]["bd"] = (0.25, 0.25, 0.3)
lipid_scaffolds[(lipid_type, params)]["lipids"] = {      #  1    2   3   4   5   6   7   8   9  10  11  12  13  14  15  16   17  18  19  20  21  22  23  24  25  26  27  28  29  30  31
    ("CDL0", "beads"): ("GL5 PO41 GL1 GL2 C1A C2A D3A C4A C5A   - C1B C2B D3B C4B C5B   - PO42 GL3 GL4 C1C C2C D3C C4C C5C   - C1D C2D D3D C4D C5D   -"),
    ("CDL1", "beads"): ("GL5 PO41 GL1 GL2 C1A C2A D3A C4A C5A   - C1B C2B D3B C4B C5B   - PO42 GL3 GL4 C1C C2C D3C C4C C5C   - C1D C2D D3D C4D C5D   -"),
    ("CDL2", "beads"): ("GL5 PO41 GL1 GL2 C1A C2A D3A C4A C5A   - C1B C2B D3B C4B C5B   - PO42 GL3 GL4 C1C C2C D3C C4C C5C   - C1D C2D D3D C4D C5D   -"),
    ("CL4P", "beads"): ("GL5 PO41 GL1 GL2 C1A C2A C3A C4A C5A   - C1B C2B C3B C4B C5B   - PO42 GL3 GL4 C1C C2C C3C C4C C5C   - C1D C2D C3D C4D C5D   -"), 
    ("CL4M", "beads"): ("GL5 PO41 GL1 GL2 C1A C2A C3A   -   -   - C1B C2B C3B   -   -   - PO42 GL3 GL4 C1C C2C C3C   -   -   - C1D C2D C3D   -   -   -"), 
}

# Prototopology for mycolic acid(s)
#
#  1--2--3--4--5--6--7--8
#                       |
# 16-15-14-13-12-11-10--9
# |
# 17-18-19-20-21-22-23-24
#                     /
# 32-31-30-29-28-27-25-26
#

lipid_type, params = "MYCOLIC ACIDS", "default"
lipid_scaffolds[(lipid_type, params)] = {}
lipid_scaffolds[(lipid_type, params)]["x"] = (      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,    0,    1,    1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,  1,   1,   1)
lipid_scaffolds[(lipid_type, params)]["y"] = (      0,   0,   0,   0,   0,   0,   0,   0,   1,   1,   1,   1,   1,   1,   1,    1,    1,    1,    1,   1,   1,   1,   1,   1,   0,   0,   0,   0,   0,  0,   0,   0)
lipid_scaffolds[(lipid_type, params)]["z"] = (      7,   6,   5,   4,   3,   2,   1,   0,   0,   1,   2,   3,   4,   5,   6,    7,    7,    6,    5,   4,   3,   2,   1,   0,   1,   0,   2,   3,   4,  5,   6,   7)
lipid_scaffolds[(lipid_type, params)]["bd"] = (0.25, 0.25, 0.3)
lipid_scaffolds[(lipid_type, params)]["lipids"] = {        # 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15    16    17    18    19   20   21   22   23   24   25   26   27   28   29   30   31   32
    ("AMA", "beads"):   ("  -    -    -  C1A  C2A  C3A  C4A  C5A  M1A  C1B  C2B  C3B  C4B    -    -     -     -     -   M1B  C1C  C2C  C3C    -    -  COH  OOH  C1D  C2D  C3D  C4D  C5D  C6D"),
    ("KMA", "beads"):   ("  -    -    -  C1A  C2A  C3A  C4A  C5A  M1A  C1B  C2B  C3B  C4B    -    -     -     -     -   M1B  C1C  C2C  C3C    -    -  COH  OOH  C1D  C2D  C3D  C4D  C5D  C6D"),
    ("MMA", "beads"):   ("  -    -    -  C1A  C2A  C3A  C4A  C5A  M1A  C1B  C2B  C3B  C4B    -    -     -     -     -   M1B  C1C  C2C  C3C    -    -  COH  OOH  C1D  C2D  C3D  C4D  C5D  C6D"),
}

### Hopanoids
lipid_type, params = "Hopanoids", "default"
lipid_scaffolds[(lipid_type, params)] = {}
lipid_scaffolds[(lipid_type, params)]["x"] = (     0,  0,  0,  0, 0.5,-0.5,   0,   0, 0.5, 0.5,   0,   0,   0,   0,  0,  0,  0,  0)
lipid_scaffolds[(lipid_type, params)]["y"] = (     0,  0,  0,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  0)
lipid_scaffolds[(lipid_type, params)]["z"] = (     0,  0,  0,  0, 0.5, 1.4, 2.6,   3, 3.3, 3.9, 4.5, 5.0, 5.5, 6.0,  0,  0,  0,  0) 
lipid_scaffolds[(lipid_type, params)]["bd"] = (0.25, 0.25, 0.3)
lipid_scaffolds[(lipid_type, params)]["lipids"] = {
    ("HOPR", "beads"): (" -   -   -   R1   R2   R3   R4   R5   R6   R7   R8   -    -    -    -   -   -   - "),
    ("HHOP", "beads"): (" -   -   -   R1   R2   R3   R4   R5   R6   R7   R8   C1   -    -    -   -   -   - "),
    ("HDPT", "beads"): (" -   -   -   R1   R2   R3   R4   R5   R6   R7   R8   C1   -    -    -   -   -   - "),
    ("HBHT", "beads"): (" -   -   -   R1   R2   R3   R4   R5   R6   R7   R8   C1   C2   C3   -   -   -   - "),
}

