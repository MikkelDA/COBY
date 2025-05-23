lipid_scaffolds = {}
lipid_metadata = {}

params = "M2"
lipid_metadata[params] = {
    "Description": "Lipids in this parameter library are based on the Martini 2 parameters as described in https://doi.org/10.1021/jp071097f.",
}

##############
### LIPIDS ###
##############
### Diacyl glycerols
lipid_type, params = "lipid", "M2"
lipid_scaffolds[(lipid_type, params)] = {}
lipid_scaffolds[(lipid_type, params)]["x"] = (    0, .5,  0,  0, .5,  0,  0, .5,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1)
lipid_scaffolds[(lipid_type, params)]["y"] = (    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0)
lipid_scaffolds[(lipid_type, params)]["z"] = (   10,  9,  9,  8,  8,  7,  6,  6,  5,  4,  3,  2,  1,  0,  5,  4,  3,  2,  1,  0)
lipid_scaffolds[(lipid_type, params)]["bd"] = (0.25, 0.25, 0.3)
lipid_scaffolds[(lipid_type, params)]["charges"] = (("NC3", 1), ("NH3", 1), ("PO4", -1))
lipid_scaffolds[(lipid_type, params)]["tags"] = ("martini2",)
lipid_scaffolds[(lipid_type, params)]["lipids"] = {}

### Phospholipids with simple head structures
#
#                     B15-B16-B17-B18-B19-B20
#                    /
#    H02   H05    L08
#   /  |   |  \    |
# H01-H03-H04-H06-L07-A09-A10-A11-A12-A13-A14

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
    "DG": "    -     -     -     -     -     - ", # Diglyceride # Different bead type in GL1 for DG in M2
}

SM_heads = {
          #   H01   H02   H03   H04   H05   H06
    "SM": "    -     -     -    NC3    -    PO4",
}

tailcodes = {
    ### Based on http://cgmartini.nl/index.php/force-field-parameters/lipids2/350-lipid-details
    "C": "C",
    "T": "CC",
    "L": "CCC",
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
                lipids[(name, "tags")] = (
                    "phospholipid",
                    "head_"+head_name,
                    "tail_"+tailname,
                )
                
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
                lipids[(name, "tags")] = (
                    "sphingomyelin",
                    "head_"+head_name,
                    "tail_"+tailname,
                )

lipid_scaffolds[(lipid_type, params)]["lipids"].update(lipids)
