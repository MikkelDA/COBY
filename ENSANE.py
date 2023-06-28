# %time
# %reset -f
import ast
import itertools
import re
import math
import re
import numpy as np
import random
import copy
import scipy.spatial
import alphashape
import shapely
import shapely.plotting
import inspect
import sys
import os
import argparse

import matplotlib.pyplot as plt
from matplotlib import patches

import pickle


from scipy.spatial.distance import cdist

################### ALL INTERNAL CALCULATIONS ARE DONE USING ÅNGSTROM

##############
### LIPIDS ###
##############
### Lipid definitions explanation:
### lipid_type: Designates which of the following values should be used together
### ff: "force-field". Can be used to create multiple libraries of lipids. Can be specified when creation membranes
### x,y,z: 0,0,0 is used as center for lipid placement
### center: z-value used for centering lipids. Also defines the cut between hydrophobic and hydrophilic
### lipids: Lipid bead definitions
### ### (lipid name, beads): Specifies that these are the beads for a specific lipid
### bd: A bead distance multiplier. Ajusts the values of x,y,z coordinates. Remnant from insane.py lipid definitions
### charges: (bead name, charge value)

### DON*T TOUCH
lipid_defs = {}
solvent_defs = {}
ion_defs = {}
prot_defs = {}


# PROTOLIPID (diacylglycerol), 18 beads
#
# 1-3-4-6-7--9-10-11-12-13-14
#  \| |/  |
#   2 5   8-15-16-17-18-19-20
#
#
## Diacyl glycerols
lipid_type, ff = "lipid", "default"
lipid_defs[(lipid_type, ff)] = {}
lipid_defs[(lipid_type, ff)]["x"] = (    0, .5,  0,  0, .5,  0,  0, .5,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1)
lipid_defs[(lipid_type, ff)]["y"] = (    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0)
lipid_defs[(lipid_type, ff)]["z"] = (   10,  9,  9,  8,  8,  7,  6,  6,  5,  4,  3,  2,  1,  0,  5,  4,  3,  2,  1,  0)
lipid_defs[(lipid_type, ff)]["center"] = 7
lipid_defs[(lipid_type, ff)]["bd"] = (0.25, 0.25, 0.3)
lipid_defs[(lipid_type, ff)]["charges"] = (("NC3", 1), ("NH3", 1), ("PO4", -1))
lipid_defs[(lipid_type, ff)]["lipids"] = {      # 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20
## Phospholipids
    ("DTPC", "beads"): (" -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A  -   -   -   -  C1B C2B  -   -   -   - "),
    ("DLPC", "beads"): (" -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A C3A  -   -   -  C1B C2B C3B  -   -   - "),
    ("DPPC", "beads"): (" -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A C3A C4A  -   -  C1B C2B C3B C4B  -   - "),
    ("DBPC", "beads"): (" -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A C3A C4A C5A  -  C1B C2B C3B C4B C5B  - "),
    ("POPC", "beads"): (" -   -   -  NC3  -  PO4 GL1 GL2 C1A D2A C3A C4A  -   -  C1B C2B C3B C4B  -   - "),
    ("DOPC", "beads"): (" -   -   -  NC3  -  PO4 GL1 GL2 C1A D2A C3A C4A  -   -  C1B D2B C3B C4B  -   - "),
    ("DAPC", "beads"): (" -   -   -  NC3  -  PO4 GL1 GL2 D1A D2A D3A D4A C5A  -  D1B D2B D3B D4B C5B  - "),
    ("DIPC", "beads"): (" -   -   -  NC3  -  PO4 GL1 GL2 C1A D2A D3A C4A  -   -  C1B D2B D3B C4B  -   - "),
    ("DGPC", "beads"): (" -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A D3A C4A C5A  -  C1B C2B D3B C4B C5B  - "),
    ("DNPC", "beads"): (" -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A C3A D4A C5A C6A C1B C2B C3B D4B C5B C6B"),
    ("DTPE", "beads"): (" -   -   -  NH3  -  PO4 GL1 GL2 C1A C2A  -   -   -   -  C1B C2B  -   -   -   - "),
    ("DLPE", "beads"): (" -   -   -  NH3  -  PO4 GL1 GL2 C1A C2A C3A  -   -   -  C1B C2B C3B  -   -   - "),
    ("DPPE", "beads"): (" -   -   -  NH3  -  PO4 GL1 GL2 C1A C2A C3A C4A  -   -  C1B C2B C3B C4B  -   - "),
    ("DBPE", "beads"): (" -   -   -  NH3  -  PO4 GL1 GL2 C1A C2A C3A C4A C5A  -  C1B C2B C3B C4B C5B  - "),
    ("POPE", "beads"): (" -   -   -  NH3  -  PO4 GL1 GL2 C1A D2A C3A C4A  -   -  C1B C2B C3B C4B  -   - "),
    ("DOPE", "beads"): (" -   -   -  NH3  -  PO4 GL1 GL2 C1A D2A C3A C4A  -   -  C1B D2B C3B C4B  -   - "),
    ("POPG", "beads"): (" -   -   -  GL0  -  PO4 GL1 GL2 C1A D2A C3A C4A  -   -  C1B C2B C3B C4B  -   - "),
    ("DOPG", "beads"): (" -   -   -  GL0  -  PO4 GL1 GL2 C1A D2A C3A C4A  -   -  C1B D2B C3B C4B  -   - "),
    ("POPS", "beads"): (" -   -   -  CN0  -  PO4 GL1 GL2 C1A D2A C3A C4A  -   -  C1B C2B C3B C4B  -   - "),
    ("DOPS", "beads"): (" -   -   -  CN0  -  PO4 GL1 GL2 C1A D2A C3A C4A  -   -  C1B D2B C3B C4B  -   - "),
    ("DPSM", "beads"): (" -   -   -  NC3  -  PO4 AM1 AM2 T1A C2A C3A  -   -   -  C1B C2B C3B C4B  -   - "),
    ("DBSM", "beads"): (" -   -   -  NC3  -  PO4 AM1 AM2 T1A C2A C3A C4A  -   -  C1B C2B C3B C4B C5B  - "),
    ("BNSM", "beads"): (" -   -   -  NC3  -  PO4 AM1 AM2 T1A C2A C3A C4A  -   -  C1B C2B C3B C4B C5B C6B"),
# PG for thylakoid membrane   
    ("OPPG", "beads"): (" -   -   -  GL0  -  PO4 GL1 GL2 C1A C2A C3A C4A  -   -  C1B D2B C3B C4B  -   - "),
# PG for thylakoid membrane of spinach (PPT with a trans-unsaturated bond at sn1 and a triple-unsaturated bond at sn2, 
# and PPG  with a transunsaturated bond at sn1 and a palmitoyl tail at sn2)
    ("JPPG", "beads"): (" -   -   -  GL0  -  PO4 GL1 GL2 C1A C2A C3A C4A  -   -  D1B C2B C3B C4B  -   - "),
    ("JFPG", "beads"): (" -   -   -  GL0  -  PO4 GL1 GL2 C1A D2A D3A D4A  -   -  D1B C2B C3B C4B  -   - "),
## Monoacylglycerol
    ("GMO", "beads"):  (" -   -   -   -   -   -  GL1 GL2 C1A C2A D3A C4A C5A  -   -   -   -   -   -   - "),
## Templates using the old lipid names and definitions
  ("DHPC.o", "beads"): (" -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A  -   -   -   -  C1B C2B  -   -   -   - "),
  ("DMPC.o", "beads"): (" -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A C3A  -   -   -  C1B C2B C3B  -   -   - "),
  ("DSPC.o", "beads"): (" -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A C3A C4A C5A  -  C1B C2B C3B C4B C5B  - "),
  ("POPC.o", "beads"): (" -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A C3A C4A  -   -  C1B C2B D3B C4B C5B  - "),
  ("DOPC.o", "beads"): (" -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A D3A C4A C5A  -  C1B C2B D3B C4B C5B  - "),
  ("DUPC.o", "beads"): (" -   -   -  NC3  -  PO4 GL1 GL2 C1A D2A D3A C4A  -   -  C1B D2B D3B C4B  -   - "),
  ("DEPC.o", "beads"): (" -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A C3A D4A C5A C6A C1B C2B C3B D4B C5B C6B"),
  ("DHPE.o", "beads"): (" -   -   -  NH3  -  PO4 GL1 GL2 C1A C2A  -   -   -   -  C1B C2B  -   -   -   - "),
  ("DLPE.o", "beads"): (" -   -   -  NH3  -  PO4 GL1 GL2 C1A C2A C3A  -   -   -  C1B C2B C3B  -   -   - "),
  ("DMPE.o", "beads"): (" -   -   -  NH3  -  PO4 GL1 GL2 C1A C2A C3A  -   -   -  C1B C2B C3B  -   -   - "),
  ("DSPE.o", "beads"): (" -   -   -  NH3  -  PO4 GL1 GL2 C1A C2A C3A C4A C5A  -  C1B C2B C3B C4B C5B  - "),
  ("POPE.o", "beads"): (" -   -   -  NH3  -  PO4 GL1 GL2 C1A C2A C3A C4A  -   -  C1B C2B D3B C4B C5B  - "),
  ("DOPE.o", "beads"): (" -   -   -  NH3  -  PO4 GL1 GL2 C1A C2A D3A C4A C5A  -  C1B C2B D3B C4B C5B  - "),
  ("PPCS.o", "beads"): (" -   -   -  NC3  -  PO4 AM1 AM2 C1A C2A C3A C4A  -   -  D1B C2B C3B C4B  -   - "),
  ("DOPG.o", "beads"): (" -   -   -  GL0  -  PO4 GL1 GL2 C1A C2A D3A C4A C5A  -  C1B C2B D3B C4B C5B  - "),
  ("POPG.o", "beads"): (" -   -   -  GL0  -  PO4 GL1 GL2 C1A C2A C3A C4A  -   -  C1B C2B D3B C4B C5B  - "),
  ("DOPS.o", "beads"): (" -   -   -  CN0  -  PO4 GL1 GL2 C1A C2A D3A C4A C5A  -  C1B C2B D3B C4B C5B  - "),
  ("POPS.o", "beads"): (" -   -   -  CN0  -  PO4 GL1 GL2 C1A C2A C3A C4A  -   -  C1B C2B D3B C4B C5B  - "),
   ("CPG.o", "beads"): (" -   -   -  GL0  -  PO4 GL1 GL2 C1A C2A C3A C4A  -   -  C1B C2B D3B C4B  -   - "),
   ("PPG.o", "beads"): (" -   -   -  GL0  -  PO4 GL1 GL2 C1A C2A C3A C4A  -   -  D1B C2B C3B C4B  -   - "),
   ("PPT.o", "beads"): (" -   -   -  GL0  -  PO4 GL1 GL2 C1A D2A D3A D4A  -   -  D1B C2B C3B C4B  -   - "),
  ("DSMG.o", "beads"): (" -   -   -  C6   C4 C1  GL1 GL2 C1A C2A C3A C4A C5A  -  C1B C2B C3B C4B C5B  - "),
  ("DSDG.o", "beads"): ("C61 C41 C11 C62 C42 C12 GL1 GL2 C1A C2A C3A C4A C5A  -  C1B C2B C3B C4B C5B  - "),
  ("DSSQ.o", "beads"): (" -   -   S6 C6   C4 C1  GL1 GL2 C1A C2A C3A C4A C5A  -  C1B C2B C3B C4B C5B  - "),
}

# HII fix for PI templates and new templates PI(s) with diffrent tails, PO-PIP1(3) and POPIP2(4,5)  
#Prototopology for phosphatidylinositol type lipids 5,6,7 are potentail phosphates (PIP1,PIP2 and PIP3)
# 1,2,3 - is the inositol and 4 is the phosphate that links to the tail part.
#  5
#   \
#  6-2-1-4-8--10-11-12-13-14-15
#    |/    |
#  7-3     9--16-17-18-19-20-21 
lipid_type, ff = "INOSITOLLIPIDS", "default"
lipid_defs[(lipid_type, ff)] = {}
lipid_defs[(lipid_type, ff)]["x"] = (   .5,  .5,   0,   0,   1, .5,  0,  0,   .5,   0,   0,   0,   0,   0,   0,   1,   1,   1,   1,   1,   1)
lipid_defs[(lipid_type, ff)]["y"] = (    0,   0,   0,   0,   0,  0,  0,  0,    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0)
lipid_defs[(lipid_type, ff)]["z"] = (    8,   9,   9,   7,  10, 10, 10,  6,    6,   5,   4,   3,   2,   1,   0,   5,   4,   3,   2,   1,   0)
lipid_defs[(lipid_type, ff)]["center"] = 7
lipid_defs[(lipid_type, ff)]["bd"] = (0.25, 0.25, 0.3)
lipid_defs[(lipid_type, ff)]["charges"] = (("NC3", 1), ("PO4", -1))
lipid_defs[(lipid_type, ff)]["lipids"] = {      # 1     2    3    4    5   6   7   8    9    10    11    12    13    14   15    16    17    18    19   20 
    ("DPPI", "beads"): (" C1   C2   C3    CP   -   -   -  GL1  GL2  C1A  C2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
    ("POPI", "beads"): (" C1   C2   C3    CP   -   -   -  GL1  GL2  C1A  D2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
    ("PIPI", "beads"): (" C1   C2   C3    CP   -   -   -  GL1  GL2  C1A  D2A  D3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
    ("PAPI", "beads"): (" C1   C2   C3    CP   -   -   -  GL1  GL2  D1A  D2A  D3A  D4A  C5A   -   C1B  C2B  C3B  C4B   -    - "),
    ("PUPI", "beads"): (" C1   C2   C3    CP   -   -   -  GL1  GL2  D1A  D2A  D3A  D4A  D5A   -   C1B  C2B  C3B  C4B   -    - "),
    ("POP1", "beads"): (" C1   C2   C3    CP  P1   -   -  GL1  GL2  C1A  C2A  D3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
    ("POP2", "beads"): (" C1   C2   C3    CP  P1  P2   -  GL1  GL2  C1A  C2A  D3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
    ("POP3", "beads"): (" C1   C2   C3    CP  P1  P2  P3  GL1  GL2  C1A  C2A  D3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
## Templates using the old lipid names and definitions
  ("PI.o", "beads")  : (" C1   C2   C3    CP   -   -   -  GL1  GL2  C1A  C2A  C3A  C4A   -    -   CU1  CU2  CU3  CU4  CU5   - "),
  ("PI34.o", "beads"): (" C1   C2   C3    CP PO1 PO2   -  GL1  GL2  C1A  C2A  C3A  C4A   -    -   CU1  CU2  CU3  CU4  CU5   - "),
}

#Prototopology for longer and branched glycosil and ceramide based glycolipids
#
#     17-15-14-16
#         |/
#        13
#         |
# 12-10-9-7-6-4-3-1--18--20-21-22-23-24
#  |/   |/  |/  |/    |
#  11   8   5   2    19--25-26-27-28-29 
lipid_type, ff = "GLYCOLIPIDS", "default"
lipid_defs[(lipid_type, ff)] = {}
lipid_defs[(lipid_type, ff)]["x"] = (    0,  .5,   0,   0,  .5,  0,  0, .5,  0,    0,   .5,    0,    0,    0,   0,    0,    0,    0,   .5,   0,   0,   0,   0,   0,   1,   1,   1,   1,   1)
lipid_defs[(lipid_type, ff)]["y"] = (    0,   0,   0,   0,   0,  0,  0,  0,  0,    0,    0,    0,   .5,    1,   1,    1,    1,    0,    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0)
lipid_defs[(lipid_type, ff)]["z"] = (    6,   7,   7,   8,   9,  9, 10, 11, 11,   12,   13,   13,   10,    9,  10,    8,   11,    5,    5,   4,   3,   2,   1,   0,   4,   3,   2,   1,   0)
lipid_defs[(lipid_type, ff)]["center"] = 6
lipid_defs[(lipid_type, ff)]["bd"] = (0.25, 0.25, 0.3)
lipid_defs[(lipid_type, ff)]["lipids"] = {      # 1     2    3    4    5   6   7   8   9    10    11    12    13    14   15    16    17    18    19   20   21   22   23   24   25   26   27   28   29
    ("DPG1", "beads"): ("GM1  GM2  GM3  GM4  GM5 GM6 GM7 GM8 GM9  GM10  GM11  GM12  GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A   -    -    -   C1B  C2B  C3B  C4B   -    - "),
    ("DXG1", "beads"): ("GM1  GM2  GM3  GM4  GM5 GM6 GM7 GM8 GM9  GM10  GM11  GM12  GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A  C4A  C5A   -   C1B  C2B  C3B  C4B  C5B  C6B"),
    ("PNG1", "beads"): ("GM1  GM2  GM3  GM4  GM5 GM6 GM7 GM8 GM9  GM10  GM11  GM12  GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A   -    -    -   C1B  C2B  C3B  D4B  C5B  C6B"),
    ("XNG1", "beads"): ("GM1  GM2  GM3  GM4  GM5 GM6 GM7 GM8 GM9  GM10  GM11  GM12  GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A  C4A  C5A   -   C1B  C2B  C3B  D4B  C5B  C6B"),
    ("DPG3", "beads"): ("GM1  GM2  GM3  GM4  GM5 GM6  -   -   -    -     -     -    GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A   -    -    -   C1B  C2B  C3B  C4B   -    - "),
    ("DXG3", "beads"): ("GM1  GM2  GM3  GM4  GM5 GM6  -   -   -    -     -     -    GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A  C4A  C5A   -   C1B  C2B  C3B  C4B  C5B  C6B"),
    ("PNG3", "beads"): ("GM1  GM2  GM3  GM4  GM5 GM6  -   -   -    -     -     -    GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A   -    -    -   C1B  C2B  C3B  D4B  C5B  C6B"),
    ("XNG3", "beads"): ("GM1  GM2  GM3  GM4  GM5 GM6  -   -   -    -     -     -    GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A  C4A  C5A   -   C1B  C2B  C3B  D4B  C5B  C6B"),
    ("DPCE", "beads"): ("  -    -    -    -    -   -   -   -   -     -     -     -     -     -    -     -     -   AM1   AM2  T1A  C2A  C3A   -    -   C1B  C2B  C3B  C4B   - "),
    ("DPGS", "beads"): (" C1   C2   C3    -    -   -   -   -   -     -     -     -     -     -    -     -     -   AM1   AM2  T1A  C2A  C3A   -    -   C1B  C2B  C3B  C4B   - "),
    ("DPMG", "beads"): (" C1   C2   C3    -    -   -   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -   C1B  C2B  C3B  C4B   - "),
    ("DPSG", "beads"): (" S1   C1   C2   C3    -   -   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -   C1B  C2B  C3B  C4B   - "),
    ("DPGG", "beads"): ("GB2  GB3  GB1  GA1  GA2 GA3   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -   C1B  C2B  C3B  C4B   - "),
#lipids for thylakoid membrane of cyanobacteria: oleoyl tail at sn1 and palmiotyl chain at sn2. SQDG no double bonds
    ("OPMG", "beads"): (" C1   C2   C3    -    -   -   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -   C1B  D2B  C3B  C4B   - "),
    ("OPSG", "beads"): (" S1   C1   C2   C3    -   -   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -   C1B  D2B  C3B  C4B   - "),
    ("OPGG", "beads"): ("GB2  GB3  GB1  GA1  GA2 GA3   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -   C1B  D2B  C3B  C4B   - "),
#lipids for thylakoid membrane of spinach: for the *T both chains are triple unsaturated and the *G have a triple unsaturated chain at sn1 and a palmitoyl chain at sn2. 
    ("FPMG", "beads"): (" C1   C2   C3    -    -   -   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -   C1B  D2B  D3B  D4B   - "),
    ("DFMG", "beads"): (" C1   C2   C3    -    -   -   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  D2A  D3A  D4A   -   C1B  D2B  D3B  D4B   - "),
    ("FPSG", "beads"): (" S1   C1   C2   C3    -   -   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -   C1B  D2B  D3B  D4B   - "),
    ("FPGG", "beads"): ("GB2  GB3  GB1  GA1  GA2 GA3   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -   C1B  D2B  D3B  D4B   - "),
    ("DFGG", "beads"): ("GB2  GB3  GB1  GA1  GA2 GA3   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  D2A  D3A  D4A   -   C1B  D2B  D3B  D4B   - "),
## Templates using the old lipid names and definitions
  ("GM1.o", "beads") : ("GM1  GM2  GM3  GM4  GM5 GM6 GM7 GM8 GM9  GM10  GM11  GM12  GM13  GM14 GM15  GM16  GM17   AM1   AM2  C1A  C2A  C3A  C4A  C5A  C1B  C2B  C3B  C4B   - "), 
  ("DGDG.o", "beads"): ("GB2  GB3  GB1  GA1  GA2 GA3   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -   C1B  C2B  C3B  C4B   - "),
  ("MGDG.o", "beads"): (" C1   C2   C3    -    -   -   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -   C1B  C2B  C3B  C4B   - "),
  ("SQDG.o", "beads"): (" S1   C1   C2   C3    -   -   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -   C1B  C2B  C3B  C4B   - "),
  ("CER.o", "beads") : ("  -    -    -    -    -   -   -   -   -     -     -     -     -     -    -     -     -   AM1   AM2  C1A  C2A  C3A  C4A   -   C1B  C2B  C3B  C4B   - "),
  ("GCER.o", "beads"): (" C1   C2   C3    -    -   -   -   -   -     -     -     -     -     -    -     -     -   AM1   AM2  C1A  C2A  C3A  C4A   -   C1B  C2B  C3B  C4B   - "),
  ("DPPI.o", "beads"): (" C1   C2   C3    -   CP   -   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -   C1B  C2B  C3B  C4B   - "),
}

lipid_type, ff = "QUINONES", "default"
lipid_defs[(lipid_type, ff)] = {}
lipid_defs[(lipid_type, ff)]["x"] = (    0,  .5,   0,    0,   0,   0,   0,   0,   0,    0,    0,    0)
lipid_defs[(lipid_type, ff)]["y"] = (    0,   0,   0,    0,   0,   0,   0,   0,   0,    0,    0,    0)
lipid_defs[(lipid_type, ff)]["z"] = (    6,   7,   7,   5.5,  5,  4.5,  4,  3.5, 2.5,   2,  1.5,    1)
lipid_defs[(lipid_type, ff)]["center"] = 6
lipid_defs[(lipid_type, ff)]["bd"] = (0.25, 0.25, 0.3)
lipid_defs[(lipid_type, ff)]["lipids"] = {      # 1     2    3    4    5    6    7    8    9    10    11    12
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
lipid_type, ff = "CARDIOLIPINS", "default"
lipid_defs[(lipid_type, ff)] = {}
lipid_defs[(lipid_type, ff)]["x"] = (   0.5,   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1)
lipid_defs[(lipid_type, ff)]["y"] = (     1,   0,  0,  1,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,   0,  0,  1,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1)
lipid_defs[(lipid_type, ff)]["z"] = (     8,   7,  6,  6,  5,  4,  3,  2,  1,  0,  5,  4,  3,  2,  1,  0,   7,  6,  6,  5,  4,  3,  2,  1,  0,  5,  4,  3,  2,  1,  0)
lipid_defs[(lipid_type, ff)]["center"] = 7
lipid_defs[(lipid_type, ff)]["bd"] = (0.25, 0.25, 0.3)
lipid_defs[(lipid_type, ff)]["lipids"] = {      #  1    2   3   4   5   6   7   8   9  10  11  12  13  14  15  16   17  18  19  20  21  22  23  24  25  26  27  28  29  30  31
    ("CDL0", "beads"): ("GL5 PO41 GL1 GL2 C1A C2A D3A C4A C5A   - C1B C2B D3B C4B C5B   - PO42 GL3 GL4 C1C C2C D3C C4C C5C   - C1D C2D D3D C4D C5D   -"), # Warning not the same names is in .itp
    ("CDL1", "beads"): ("GL5 PO41 GL1 GL2 C1A C2A D3A C4A C5A   - C1B C2B D3B C4B C5B   - PO42 GL3 GL4 C1C C2C D3C C4C C5C   - C1D C2D D3D C4D C5D   -"), # Warning not the same names is in .itp
    ("CDL2", "beads"): ("GL5 PO41 GL1 GL2 C1A C2A D3A C4A C5A   - C1B C2B D3B C4B C5B   - PO42 GL3 GL4 C1C C2C D3C C4C C5C   - C1D C2D D3D C4D C5D   -"), # Warning not the same names is in .itp 
    ("CL4P", "beads"): ("GL5 PO41 GL1 GL2 C1A C2A C3A C4A C5A   - C1B C2B C3B C4B C5B   - PO42 GL3 GL4 C1C C2C C3C C4C C5C   - C1D C2D C3D C4D C5D   -"), 
    ("CL4M", "beads"): ("GL5 PO41 GL1 GL2 C1A C2A C3A   -   -   - C1B C2B C3B   -   -   - PO42 GL3 GL4 C1C C2C C3C   -   -   - C1D C2D C3D   -   -   -"), 
## Templates using the old lipid names and definitions
  ("CL4.o", "beads") : ("GL5 PO41 GL1 GL2 C1A C2A D3A C4A C5A   - C1B C2B D3B C4B C5B   - PO42 GL3 GL4 C1C C2C D3C C4C C5C   - C1D C2D D3D C4D C5D   -"), 
  ("CL4O.o", "beads"): ("GL5 PO41 GL1 GL2 C1A C2A D3A C4A C5A   - C1B C2B D3B C4B C5B   - PO42 GL3 GL4 C1C C2C D3C C4C C5C   - C1D C2D D3D C4D C5D   -"),
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

lipid_type, ff = "MYCOLIC ACIDS", "default"
lipid_defs[(lipid_type, ff)] = {}
lipid_defs[(lipid_type, ff)]["x"] = (      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,    0,    1,    1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,  1,   1,   1)
lipid_defs[(lipid_type, ff)]["y"] = (      0,   0,   0,   0,   0,   0,   0,   0,   1,   1,   1,   1,   1,   1,   1,    1,    1,    1,    1,   1,   1,   1,   1,   1,   0,   0,   0,   0,   0,  0,   0,   0)
lipid_defs[(lipid_type, ff)]["z"] = (      7,   6,   5,   4,   3,   2,   1,   0,   0,   1,   2,   3,   4,   5,   6,    7,    7,    6,    5,   4,   3,   2,   1,   0,   1,   0,   2,   3,   4,  5,   6,   7)
lipid_defs[(lipid_type, ff)]["center"] = 7
lipid_defs[(lipid_type, ff)]["bd"] = (0.25, 0.25, 0.3)
lipid_defs[(lipid_type, ff)]["lipids"] = {        # 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15    16    17    18    19   20   21   22   23   24   25   26   27   28   29   30   31   32
    ("AMA", "beads"):   ("  -    -    -  C1A  C2A  C3A  C4A  C5A  M1A  C1B  C2B  C3B  C4B    -    -     -     -     -   M1B  C1C  C2C  C3C    -    -  COH  OOH  C1D  C2D  C3D  C4D  C5D  C6D"),
    ("AMA.w", "beads"): ("  -    -    -  C1A  C2A  C3A  C4A  C5A  M1A  C1B  C2B  C3B  C4B    -    -     -     -     -   M1B  C1C  C2C  C3C    -    -  COH  OOH  C1D  C2D  C3D  C4D  C5D  C6D"),
    ("KMA", "beads"):   ("  -    -    -  C1A  C2A  C3A  C4A  C5A  M1A  C1B  C2B  C3B  C4B    -    -     -     -     -   M1B  C1C  C2C  C3C    -    -  COH  OOH  C1D  C2D  C3D  C4D  C5D  C6D"),
    ("MMA", "beads"):   ("  -    -    -  C1A  C2A  C3A  C4A  C5A  M1A  C1B  C2B  C3B  C4B    -    -     -     -     -   M1B  C1C  C2C  C3C    -    -  COH  OOH  C1D  C2D  C3D  C4D  C5D  C6D"),
}

# Sterols
lipid_type, ff = "sterol", "default"
lipid_defs[(lipid_type, ff)] = {}
lipid_defs[(lipid_type, ff)]["x"] = (     0,  0,  0,  0,  0, 0,   0,  0,  1,  0,  1,  0,  0,  0,  0,  0,  0,  0)
lipid_defs[(lipid_type, ff)]["y"] = (     0,  0,  0,  0,  0, 0,   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0)
lipid_defs[(lipid_type, ff)]["z"] = (     0,  0,  0,  0,  0, 0, 5.3,4.5,3.9,3.3, 3 ,2.6,1.4,  0,  0,  0,  0,  0)
lipid_defs[(lipid_type, ff)]["center"] = 4.9
lipid_defs[(lipid_type, ff)]["bd"] = (0.25, 0.25, 0.3)
lipid_defs[(lipid_type, ff)]["lipids"] = {
    ("CHOL", "beads"): (" -   -   -   -   -   -  ROH  R1  R2  R3  R4  R5  C1  C2  -   -   -   - "),
    ("ERGO", "beads"): (" -   -   -   -   -   -  ROH  R1  R2  R3  R4  R5  C1  C2  -   -   -   - "),
}

# Hopanoids
lipid_type, ff = "Hopanoids", "default"
lipid_defs[(lipid_type, ff)] = {}
lipid_defs[(lipid_type, ff)]["x"] = (     0,  0,  0,  0, 0.5,-0.5,   0,   0, 0.5, 0.5,   0,   0,   0,   0,  0,  0,  0,  0)
lipid_defs[(lipid_type, ff)]["y"] = (     0,  0,  0,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  0)
lipid_defs[(lipid_type, ff)]["z"] = (     0,  0,  0,  0, 0.5, 1.4, 2.6,   3, 3.3, 3.9, 4.5, 5.0, 5.5, 6.0,  0,  0,  0,  0) 
lipid_defs[(lipid_type, ff)]["center"] = 4.9
lipid_defs[(lipid_type, ff)]["bd"] = (0.25, 0.25, 0.3)
lipid_defs[(lipid_type, ff)]["lipids"] = {
    ("HOPR", "beads"): (" -   -   -   R1   R2   R3   R4   R5   R6   R7   R8   -    -    -    -   -   -   - "),
    ("HHOP", "beads"): (" -   -   -   R1   R2   R3   R4   R5   R6   R7   R8   C1   -    -    -   -   -   - "),
    ("HDPT", "beads"): (" -   -   -   R1   R2   R3   R4   R5   R6   R7   R8   C1   -    -    -   -   -   - "),
    ("HBHT", "beads"): (" -   -   -   R1   R2   R3   R4   R5   R6   R7   R8   C1   C2   C3   -   -   -   - "),
}

################
### SOLVENTS ###
################
### molar_mass: [g/mol]
### density: [g/cm^3]
ff = "default"
solvent_defs[ff] = {
     "W" : {"beads": "W",  "x": (0,), "y": (0,), "z": (0,), "solvcount": 4, "density": 0.99669, "molar_mass": 18.01528},
     "SW": {"beads": "SW", "x": (0,), "y": (0,), "z": (0,), "solvcount": 3, "density": 0.99669, "molar_mass": 18.01528},
     "TW": {"beads": "TW", "x": (0,), "y": (0,), "z": (0,), "solvcount": 2, "density": 0.99669, "molar_mass": 18.01528},
}
### Amino acids
solvent_defs[ff].update({
     "GLY": {"beads": ("BB",), "solvcount": 1, "x": (0,), "y": (0,), "z": (0,)},
     "ALA": {"beads": ("BB",), "solvcount": 1, "x": (0,), "y": (0,), "z": (0,)},

     "ASN": {"beads": ("BB", "SC1"), "solvcount": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
     "ASP": {"beads": ("BB", "SC1"), "solvcount": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
     "GLU": {"beads": ("BB", "SC1"), "solvcount": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
     "GLN": {"beads": ("BB", "SC1"), "solvcount": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
     "LEU": {"beads": ("BB", "SC1"), "solvcount": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
     "ILE": {"beads": ("BB", "SC1"), "solvcount": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
     "VAL": {"beads": ("BB", "SC1"), "solvcount": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
     "SER": {"beads": ("BB", "SC1"), "solvcount": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
     "THR": {"beads": ("BB", "SC1"), "solvcount": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
     "CYS": {"beads": ("BB", "SC1"), "solvcount": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
     "LYS": {"beads": ("BB", "SC1"), "solvcount": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
     "PRO": {"beads": ("BB", "SC1"), "solvcount": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
     "HYP": {"beads": ("BB", "SC1"), "solvcount": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},

     "ARG": {"beads": ("BB", "SC1", "SC2"), "solvcount": 1, "x": (0.25, 0, -0.25), "y": (0, 0, 0.125), "z": (0, 0, 0)},
     "PHE": {"beads": ("BB", "SC1", "SC2", "SC3"), "solvcount": 1, "x": (0.25, 0, -0.25, -0.25), "y": (0, 0, 0.125, -0.125), "z": (0, 0, 0, 0)},
     "TYR": {"beads": ("BB", "SC1", "SC2", "SC3"), "solvcount": 1, "x": (0.25, 0, -0.25, -0.25), "y": (0, 0, 0.125, -0.125), "z": (0, 0, 0, 0)},
     "TRP": {"beads": ("BB", "SC1", "SC2", "SC3", "SC4"), "solvcount": 1, "x": (0.25, 0.25, 0, 0, -0.25), "y": (0.125, 0, -0.125, 0.125, 0), "z": (0, 0, 0, 0, 0)},
})

############
### IONS ###
############
ff = "default"
ion_defs[ff] = {}
ion_defs[ff]["positive"] = {
     "NA": {"beads": "NA", "charge": 1, "x": (0,), "y": (0,), "z": (0,)},
     "TMA": {"beads": "TMA", "charge": 1, "x": (0,), "y": (0,), "z": (0,)},
     "CA": {"beads": "CA", "charge": 2, "x": (0,), "y": (0,), "z": (0,)},
}

ion_defs[ff]["negative"] = {
     "CL": {"beads": "CL", "charge": -1, "x": (0,), "y": (0,), "z": (0,)},
     "BR": {"beads": "BR", "charge": -1, "x": (0,), "y": (0,), "z": (0,)},
     "IOD": {"beads": "ID", "charge": -1, "x": (0,), "y": (0,), "z": (0,)},
     "ACE": {"beads": "CL", "charge": -1, "x": (0,), "y": (0,), "z": (0,)},
     "BF4": {"beads": "BF4", "charge": -1, "x": (0,), "y": (0,), "z": (0,)},
     "PF6": {"beads": "PF6", "charge": -1, "x": (0,), "y": (0,), "z": (0,)},
     "SCN": {"beads": "SCN", "charge": -1, "x": (0,), "y": (0,), "z": (0,)},
     "CLO4": {"beads": "CLO", "charge": -1, "x": (0,), "y": (0,), "z": (0,)},
     "NO3": {"beads": "NO3", "charge": -1, "x": (0,), "y": (0,), "z": (0,)},
}

###########################
### PROTEIN CHARGE DATA ###
###########################
ff = "default"
prot_defs[ff] = {}
prot_defs[ff]["charges"] = {
    "ARG":1, "LYS":1, "ASP":-1, "GLU":-1, "DLPG":-1, "DMPG":-1, "DPPG":-1, "DSPG":-1, "POPG":-1,
    "SOPG":-1, "DOPG":-1, "POPS":-1, "DOPS":-1, "SAPS":-1, "POPA":-1, "DSSQ":-1,"SAP2_45":-4
}



########################################################################################################
########################################### THE ACTUAL CLASS ###########################################
########################################################################################################

class ENSANE:
    
    def __init__(self, **kwargs):
        try:
            self.lipid_defs = lipid_defs.copy()
        except:
            self.lipid_defs = {}
        try:
            self.solvent_defs = solvent_defs.copy()
        except:
            self.solvent_defs = {}
        try:
            self.ion_defs = ion_defs.copy()
        except:
            self.ion_defs = {}
        try:
            self.prot_defs = prot_defs.copy()
        except:
            self.prot_defs = {}
        
        self.PROTEINS = {}
        self.PROTEINS_cmds = []
        
        self.LEAFLETS = {}
        self.LEAFLETS_cmds = []
        
        self.SOLVATIONS = {}
        self.SOLVATIONS_cmds = []
        
        self.itp_moltypes = {}
        self.ITP_INPUT_cmds = []
        
        self.SOLUTE_INPUT_cmds = []
        
        self.PLOT_cmd = []
        self.plot_data = {}
        self.plots_requested = False
        
        self.PICKLE_cmd = False
        
        self.sys_params = "default"
        self.system_charge = 0
        self.system_name = "PLACEHOLDER_TITLE"
        
        self.output_system_file_name = "output.pdb"
        self.output_topol_file_name = "topol.top"
        
        self.LOG_FILE = []
        self.output_log_file_name = False
        
        self.pbc_set = []
        self.backup = True
        self.pickle = False
        
        self.debug_prints = False
        self.extra_info = True
        self.warnings = True
        self.quiet = False
        
        self.RUN = True
        
        self.commands_handler(kwargs)
    
    ##############################
    ### GIVE COMMANDS TO CLASS ###
    ##############################
    def commands_handler(self, kwargs):
        
        for key, cmd in kwargs.items():
            ### General system inputs
            if any(key.startswith(i) for i in ["protein", "prot"]):
                if type(cmd) != list:
                    cmd = [cmd]
                for subcmd in cmd:
                    self.PROTEINS_cmds.extend([subcmd])
                
            if any(key.startswith(i) for i in ["membrane", "memb"]):
                if type(cmd) != list:
                    cmd = [cmd]
                for subcmd in cmd:
                    self.LEAFLETS_cmds.extend([subcmd])
                
            if any(key.startswith(i) for i in ["solvation", "solv"]):
                if type(cmd) != list:
                    cmd = [cmd]
                for subcmd in cmd:
                    self.SOLVATIONS_cmds.extend([subcmd])
            
            ### Box size
            if key in ["pbc", "box", "pbc_box"]:
                if len(cmd) == 2:
                    self.pbcx, self.pbcy = cmd[0]*10
                    self.pbcz = cmd[1]*10
                    self.pbc_box = [self.pbcx, self.pbcy, self.pbcz]
                elif len(cmd) == 3:
                    self.pbcx, self.pbcy, self.pbcz = [i*10 for i in cmd]
                    self.pbc_box = [self.pbcx, self.pbcy, self.pbcz]
                else:
                    assert False, "Incorrect pbc box dimensions: " + str(key, cmd)
            
            ### Imports
            if key in ["itp_input", "itp_in"]:
                if type(cmd) != list:
                    cmd = [cmd]
                for subcmd in cmd:
                    self.ITP_INPUT_cmds.extend([subcmd])
            
            if key in ["solute_input", "solute_in"]:
                if type(cmd) != list:
                    cmd = [cmd]
                for subcmd in cmd:
                    self.SOLUTE_INPUT_cmds.extend([subcmd])
            
            ### Outputs
            if key in ["output_system", "out_sys", "out"]:
                if any([cmd.endswith(i) for i in ["pdb", "gro"]]):
                    self.output_system_file_name = cmd
                else:
                    if "." in cmd:
                        cmd = cmd.split(".")[0]
                    self.output_system_file_name = cmd + ".pdb"
            
            if key in ["output_topol", "out_top", "top_out", "top_output"]:
                self.output_topol_file_name = cmd
                
            if key in ["output_log", "out_log", "log_out", "log"]:
                self.output_log_file_name = cmd
                
#             if key in ["imp_o", "output_imported"]:
#                 self.output_imported = cmd
            
            ### Extra functionalities
            if key in ["plot"]:
                self.PLOT_cmd = cmd
                if self.PLOT_cmd:
                    self.plots_requested = True
                    
            if key in ["pickle"]:
                self.PICKLE_cmd = cmd
                
            if key in ["backup"]:
                self.backup = cmd

            if key in ["params", "sys_params"]:
                self.sys_ff = cmd
            
            ### Printer settings
            if key in ["quiet"]:
                self.quiet = cmd
                
            if key in ["debug"]:
                self.debug_prints = cmd
                
            if key in ["extra"]:
                self.extra_info = cmd
                
            if key in ["warn"]:
                self.warnings = cmd
                
            
            ### Run the program
            if key in ["run"]:
                self.RUN = cmd
        
        if self.RUN:
            self.run()
        
    def run(self):
        '''
        Runs the entire system creation process
        '''
        ### Initial checks
        assert len(self.pbc_box) > 0, "Box dimensions not set. Please do so using 'box=[x,y,z]'"
        assert len(self.pbc_box) == 3, "Box dimensions improperly defined. 3 dimensions must be given"
        assert all([self.is_number(ax)[0] for ax in self.pbc_box]), "Not all box values are numbers"
        
        assert any([len(cmd) > 0 for cmd in [self.PROTEINS_cmds, self.LEAFLETS_cmds, self.SOLVATIONS_cmds]]), "Running requires at least one command to be supplied of either type 'protein', 'membrane' or 'solvation'"
        
        ### Topology
        self.itp_read_initiater()
        
        self.print_term("------------------------------ PREPROCESSING")
        ### Defs
        self.import_structures_handler()
        self.lipid_defs_preprocessor()
        
        ### Commands
        self.prot_preprocessor()
        self.leaf_preprocessor()
        self.solv_preprocessor()
        self.print_term("------------------------------ PREPROCESSING COMPLETE", "\n")
        
        ### Run the program
        self.prot_placer()
        self.leaflet_rectangle()
        self.overlap_checker()
        self.leaflet_adjuster()
        self.lipid_inserter()
        self.solvater()
        
        self.print_term("--------------------")
        self.print_term("Final system charge:", self.system_charge)
        self.print_term("--------------------", "\n")
        
        self.plotter()
        self.pickler()
        
        ### Write the files
        self.system_file_writer()
        self.topol_file_writer()
        self.log_file_writer()
        
        self.print_term("My task is complete. Did i do a good job?")
    
    #####################################
    ### Specialized printing function ###
    #####################################
    def print_term(self, *string, debug = False, extra = False, warn = False):
        '''
        Specialized printing function.
        Allows for easy customization of requirements for when a print statement should be printed
        '''
        print_true = False
        if type(string) in [tuple, list]:
            string = " ".join([str(i) for i in string])

        ### Check if string should be printed based on settings
        if debug:
            if self.debug_prints:
                print_true = True
        elif extra:
            if self.extra_info:
                print_true = True
        elif warn:
            if self.warnings:
                print_true = True
        else:
            ### If no special case, then assume it should be printed
            print_true = True

        if print_true and not self.quiet:
            ### Print to terminal if not quiet
            print(string)

        ### Appends string to log file, for later writing
        if print_true and self.output_log_file_name:
            self.LOG_FILE.append(string)
    
    ###############
    ### Readers ###
    ###############
    def pdb_reader(self, pdb_file_dest):
        '''
        Reads a pdb file and returns it as a dictionary
        '''
        pdb_column_lengths = [6, 5, 5, 5, 1, 4, 4, 8, 8, 8, 6, 6, 4, 2, 2]
        processed_file = {}
        atom_nr = 0
        with open(pdb_file_dest, "r") as input_file:
            for line_nr, line in enumerate(input_file):
                if any([line.startswith(string) for string in ["ATOM", "HETATM"]]):
                    atom_dict = {}
                    key_indexes = [
                        ("atom_nr"  , 6, 11, int),
                        ("atom_name", 11, 16, str),
                        ("res_name" , 16, 21, str),
                        ("res_nr"   , 22, 26, int),
                        ("x"        , 30, 38, float), # [Å]
                        ("y"        , 38, 46, float), # [Å]
                        ("z"        , 46, 54, float), # [Å]
                    ]
                    for key, i1, i2, func in key_indexes:
                        if len(line) >= i2:
                            atom_dict[key] = func(line[i1:i2].replace(" ",""))
                    processed_file[(atom_nr, atom_dict["atom_nr"], atom_dict["res_nr"])] = atom_dict
                    atom_nr += 1
        return processed_file
    
    def gro_reader(self, gro_file_dest):
        '''
        Reads a gro file and returns it as a dictionary
        '''
        gro_column_lengths = [5, 5, 5, 5, 8, 8, 8, 8, 8, 8]
        processed_file = {}
        atom_nr = 0
        with open(gro_file_dest, "r") as input_file:
            file_length = len(input_file.readlines())
        with open(gro_file_dest, "r") as input_file:
            for line_nr, line in enumerate(input_file):
#                 if file_length - 2 > line_nr > 1:
                if line_nr > 1:
                    EOF_length = len(line) == 31
                    if EOF_length:
                        EOF_x = line[:10]
                        EOF_y = line[10:20]
                        EOF_z = line[20:30]
                        try:
                            EOF_x_lit = ast.literal_eval(EOF_x)
                            EOF_y_lit = ast.literal_eval(EOF_y)
                            EOF_z_lit = ast.literal_eval(EOF_z)
                            if all([type(i) in [int, float] for i in [EOF_x_lit, EOF_y_lit, EOF_z_lit]]):
                                break
                        except:
                            continue
                    if len(line):
                        atom_dict = {}
                        key_indexes = [
                            ("res_nr"   , 0, 5, int),
                            ("res_name" , 5, 10, str),
                            ("atom_name", 10, 15, str),
                            ("atom_nr"  , 15, 20, int),
                            ("x"        , 20, 28, float), # [nm]
                            ("y"        , 28, 36, float), # [nm]
                            ("z"        , 36, 44, float), # [nm]
                        ]
                        for key, i1, i2, func in key_indexes:
                            if len(line) >= i2:
                                atom_dict[key] = func(line[i1:i2].replace(" ",""))
                                if func == float: # Convert [nm] to [Å]
                                    atom_dict[key] *= 10
                        processed_file[(atom_nr, atom_dict["atom_nr"], atom_dict["res_nr"])] = atom_dict
                        atom_nr += 1
        return processed_file
    
    def itp_reader(self, itp_file_dest):
        '''
        Reads itp files
        Understands itp file definitions using "#ifdef" and "#endif"
        Understands references to other files using "#include" (by calling itself recursively on the files)
        Currently only charges are used for anything
        Ignores all types of virtual site definitions as of right now
        '''
        def line_filter(string):
            return list(filter(None, string.split()))

        dict_names = {
            "atoms": ["i", "type", "resnr", "resid", "atom", "cgnr", "charge", "mass"],
            "bonds": ["i", "j", "funct", "length", "fc"],
            "angles": ["i", "j", "k", "funct", "angle", "fc"],
            "dihedrals": ["i", "j", "k", "l", "funct", "angle", "fc", "multiplicity"],
    #         "virtual_sites3": ["site", "i", "j", "k", "funct", "a", "b", "c"], # 3-body virtual site (fd)
    #         "virtual_sitesn": ["site", "funct", "i", "j", "k", "l", "o", "p", "n", "m"], # N-body virtual site (Center Of Mass (COM))
        }

        defs_dict_names = {
            "bondtypes": ["funct", "length", "fc"],
            "angletypes": ["funct", "angle", "fc"],
            "dihedraltypes": ["funct", "angle", "fc"],
        }

        with open(itp_file_dest, "r") as input_file:
            cur_cmd = False
            ifdef = False
            for line_nr, line in enumerate(input_file):
                ### Comments, space-only lines or empty lines
                line = line.split(";")[0]
                if line.isspace() or line == "":
                    continue

                ### If lines are within an "#ifdef" statement then don't include
                elif line.startswith("#ifdef"):
                    ifdef = True
                elif line.startswith("#endif"):
                    ifdef = False
                elif ifdef:
                    continue

                ### Recursively calls the function for '#include' statements
                elif line.startswith("#include"):
                    string = line_filter(line)
                    inc_path = os.path.join("/".join(itp_file_dest.split("/")[:-1]), string[1].replace('"', ''))
                    self.itp_reader(inc_path)
#                     molecule_dict, defs = self.itp_reader(inc_path, defs, molecule_dict)
#                     self.top_mols, self.top_defs = self.itp_reader(top_file, top_defs, top_mols)

                ### Determine current command
                elif "[" in line and "]" in line:
                    cur_cmd = line[line.find("[")+len("["):line.rfind("]")].replace(" ", "")  #"moleculetype"
                    ### Break loop if system name found
                    if cur_cmd == "system":
                        break

                ### Continue if unnecessary data types are currently being run through
                elif cur_cmd in ["defaults", "atomtypes", "nonbond_params", "constraints", "exclusions", "virtual_sites3", "virtual_sitesn"]:
                    continue

                ### '#defines'
                elif line.startswith("#define"):
                    string = line_filter(line)
                    self.itp_defs[cur_cmd][string[1]] = {}
                    for i, st in enumerate(string[:2]):
                        self.itp_defs[cur_cmd][string[1]][defs_dict_names[cur_cmd][i]] = st

                ### Determine new moleculetype
                elif cur_cmd == "moleculetype":
                    string = line_filter(line)
                    self.itp_moltypes[string[0]] = {
                        "atoms": {},
                        "bonds": {},
                        "angles": {},
                        "dihedrals": {},
    #                     "virtual_sitesn": {},
    #                     "virtual_sites3": {},
                    }
                    moltype = string[0]

                ### Atoms are specially treated. Might not be necessary
                elif cur_cmd == "atoms":
                    string = line_filter(line)
                    self.itp_moltypes[moltype][cur_cmd][string[0]] = {}
                    for i, st in enumerate(string):
                        self.itp_moltypes[moltype][cur_cmd][string[0]][dict_names[cur_cmd][i]] = string[i]

                ### Bonds, angles and dihedrals
                elif [cur_cmd == cmd for cmd in ["bonds", "angles", "dihedrals"]]:
                    string = line_filter(line)
                    if cur_cmd == "bonds":
                        ids = 2
                    if cur_cmd == "angles":
                        ids = 3
                    if cur_cmd == "dihedrals":
                        ids = 4
    #                 if cur_cmd.startswith("virtual_sites"):
    #                     ids = 1
                    self.itp_moltypes[moltype][cur_cmd][tuple(string[0:ids])] = {}
                    cur_cmd_type = cur_cmd[:-1] + "types"
                    for i, st in enumerate(string):
                        if st in self.itp_defs[cur_cmd_type]:
                            for keyi, (key, val) in enumerate(self.itp_defs[cur_cmd_type][st].items()):
                                self.itp_moltypes[moltype][cur_cmd][tuple(string[0:ids])][defs_dict_names[cur_cmd_type][keyi]] = val
                            break
                        else:
                            self.itp_moltypes[moltype][cur_cmd][tuple(string[0:ids])][dict_names[cur_cmd][i]] = string[i]
    
    ###############
    ### Writers ###
    ###############
    def pdb_atom_writer(self, ATOM, a_nr, a_name, aLoc, r_name, chain, r_nr, aChar, x, y, z, oc, T, JUMP, E, C):
        '''
        PCL: pdb ATOM column lengths. "+" means that multiple column groups are combined
        '''
        PCL = [6, 5, 1 + 4, 0, 1 + 3 + 1, 1, 4, 1 + 3, 8, 8, 8, 6, 6, 10, 2, 2] 
        if len(r_name) > 4:
            r_name = r_name[:4]
        if len(a_name) > 5:
            a_name = a_name[:5]
        string = '{ATOM:<{L0}}{a_nr:>{L1}}{a_name:>{L2}}{aLoc:>{L3}}{r_name:>{L4}}{chain:^{L5}}{r_nr:>{L6}}{aChar:>{L7}}{x:>{L8}.3f}{y:>{L9}.3f}{z:>{L10}.3f}{oc:>{L11}.2f}{T:>{L12}.2f}{JUMP:>{L13}}{E:>{L14}}{C:>{L15}}'.format(
            ATOM = ATOM, L0 = PCL[0],
            a_nr = a_nr, L1 = PCL[1], # int
            a_name = a_name, L2 = PCL[2],
    #         aLoc = aLoc, L3 = PCL[3],
            aLoc = "", L3 = 0,
            r_name = r_name, L4 = PCL[4],
            chain = chain, L5 = PCL[5],
            r_nr = r_nr, L6 = PCL[6], # int
            aChar = aChar, L7 = PCL[7],
            x = float(x), L8 = PCL[8], # float
            y = float(y), L9 = PCL[9], # float
            z = float(z), L10 = PCL[10], # float
            oc = float(oc), L11 = PCL[11], # float
            T = float(T), L12 = PCL[12], # float
            JUMP = JUMP, L13 = PCL[13],
            E = E, L14 = PCL[14],
            C = C, L15 = PCL[15],
        )
        return string

    def gro_atom_writer(self, r_nr, r_name, a_name, a_nr, x, y, z, vx, vy, vz):
        GCL = [5, 5, 5, 5, 8, 8, 8, 8, 8, 8] # gro atom column lengths
    #     string = '{r_nr:>{L0}}{r_name:<{L1}}{a_name:>{L2}}{a_nr:<{L3}}{x:>{L4}.3f}{y:>{L5}.3f}{z:>{L6}.3f}{vx:>{L7}.4f}{vy:>{L8}.4f}{vz:>{L9}.4f}'.format(
        if len(r_name) > 5:
            r_name = r_name[:5]
        if len(a_name) > 4:
            a_name = a_name[:5]
        string = '{r_nr:>{L0}}{r_name:<{L1}}{a_name:>{L2}}{a_nr:>{L3}}{x:>{L4}.3f}{y:>{L5}.3f}{z:>{L6}.3f}{vx:>{L7}}{vy:>{L8}}{vz:>{L9}}'.format(
            r_nr = r_nr, L0 = GCL[0], # int
            r_name = r_name, L1 = GCL[1], # str
            a_name = a_name, L2 = GCL[2], # str
            a_nr = a_nr, L3 = GCL[3], # int
            x = float(x), L4 = GCL[4], # float
            y = float(y), L5 = GCL[5], # float
            z = float(z), L6 = GCL[6], # float
            vx = vx, L7 = GCL[7], # float # str placeholder
            vy = vy, L8 = GCL[8], # float # str placeholder
            vz = vz, L9 = GCL[9], # float # str placeholder
        )
        return string
    
    def system_file_writer(self):
        self.print_term("Writing structure file:", self.output_system_file_name)
        ### ### Creating beginning of file lines
        if self.output_system_file_name.endswith("pdb"):
            output_system_file_lines = [
                "TITLE     " + self.system_name,
                "REMARK    " + "PLACEHOLDER_REMARK",
                '{Rname:<{RnameL}}{a:>{aL}.3f}{b:>{bL}.3f}{c:>{cL}.3f}{alpha:>{alphaL}.2f}{beta:>{betaL}.2f}{gamma:>{gammaL}.2f} {sGroup:<{sGroupL}}{z:>{zL}}'.format(
                    Rname = "CRYST1", RnameL = 6,
                    a = float(self.pbc_box[0]), aL = 9, b = float(self.pbc_box[1]), bL = 9, c = float(self.pbc_box[2]), cL = 9,
                    alpha = 90, alphaL = 7, beta = 90, betaL = 7, gamma = 90, gammaL = 7,
                    sGroup = "P 1", sGroupL = 11, z = 1, zL = 4,
                ),
                "MODEL        1",
        #         "".join([str(i) for i in range(1, 10)]) + "".join([str(i) for _ in range(8) for i in range(0, 10)]),
            ]
        elif self.output_system_file_name.endswith("gro"):
            output_system_file_lines = [
                self.system_name,
                "PLACEHOLDER_ATOM_COUNT",
        #         "".join([str(i) for i in range(1, 10)]) + "".join([str(i) for _ in range(8) for i in range(0, 10)]),
            ]

        atom_count = 0
        self.molecules_for_top = []

        old_res_nr = 0
        atom_nr = 0
        res_nr = 0
        ### ### Writing proteins lines
        if len(self.PROTEINS_cmds) != 0:
            for protein_nr, protein in self.PROTEINS.items():
                for mol_name in protein["mol_names"]:
                    self.molecules_for_top.append((mol_name, str(1)))
                current_prot_res = 0
                for (i, atom, res), bead_vals in protein["beads_centered"].items():
                    if current_prot_res != res:
                        current_prot_res = res
                        res_nr += 1
                        if res_nr >= 10000:
                            res_nr -= 10000 * (res_nr // 10000)
                    atom_nr += 1
                    atom_count += 1
                    if atom_nr >= 100000:
                        atom_nr -= 100000 * (atom_nr // 100000)
                    x = bead_vals["x"] + (self.pbc_box[0] / 2)
                    y = bead_vals["y"] + (self.pbc_box[1] / 2)
                    z = bead_vals["z"] + (self.pbc_box[2] / 2)
                    a_name = bead_vals["atom_name"]
                    r_name = bead_vals["res_name"]
                    if self.output_system_file_name.endswith("pdb"):
                        output_system_file_lines.append(self.pdb_atom_writer("ATOM", atom_nr, a_name, " ", r_name, "A", res_nr, " ", float(x), float(y), float(z), float(1), float(0), " ", " ", " "))
                    elif self.output_system_file_name.endswith("gro"): ### gro coordinates are in [nm] not [Å]
                        output_systemt_file_lines.append(self.gro_atom_writer(res_nr, r_name, a_name, atom_nr, x / 10, y / 10, z / 10, " ", " ", " "))

        ### ### Writing leaflets lines
        if len(self.LEAFLETS_cmds) != 0:
            for leaflet_nr, leaflet in self.LEAFLETS.items():
                self.molecules_for_top.extend(sorted(leaflet["leaf_lipid_count"], key=lambda m: m[0]))

                ### Ordering lipids per leaflet for concise topology creation
                ordered_lipids = sorted(
                    [(grid_point, grid_vals) for (grid_point, grid_vals) in leaflet["grid_adjusted"].items()],
                    key=lambda gp: gp[1]["lipid"]["name"]
                )

                for i, (grid_point, grid_vals) in enumerate(ordered_lipids):
                    res_nr += 1
                    if res_nr >= 10000:
                        res_nr -= 10000 * (res_nr // 10000)
                    r_name = grid_vals["lipid"]["name"]

                    for x, y, z, bead_name in zip(grid_vals["lipid"]["x"], grid_vals["lipid"]["y"], grid_vals["lipid"]["z"], grid_vals["lipid"]["beads"]):
                        atom_nr += 1
                        atom_count += 1
                        if atom_nr >= 100000:
                            atom_nr -= 100000 * (atom_nr // 100000)
                        x = x + (self.pbc_box[0] / 2)
                        y = y + (self.pbc_box[1] / 2)
                        z = z + (self.pbc_box[2] / 2)
                        a_name = bead_name
                        string = ""
                        if self.output_system_file_name.endswith("pdb"):
                            output_system_file_lines.append(self.pdb_atom_writer("ATOM", atom_nr, a_name, " ", r_name, "A", res_nr, " ", float(x), float(y), float(z), float(1), float(0), " ", " ", " "))
                        elif self.output_system_file_name.endswith("gro"): ### gro coordinates are in [nm] not [Å]
                            output_system_file_lines.append(self.gro_atom_writer(res_nr, r_name, a_name, atom_nr, x / 10, y / 10, z / 10, " ", " ", " "))

        ### ### Writing solvent/ion lines
        if len(self.SOLVATIONS_cmds) != 0:
            for solvation_nr, solvation in self.SOLVATIONS.items():
                solvent_count = [(key_name, count) for (key_name, key_type, key_charge), count in solvation["solv_count"].items()]
                self.molecules_for_top.extend(solvent_count)

                for i, grid_point_3D in enumerate(solvation["grid"]):
                    res_nr += 1
                    if res_nr >= 10000:
                        res_nr -= 10000 * (res_nr // 10000)

                    r_name = grid_point_3D["name"]

                    for (x, y, z), bead_name in zip(grid_point_3D["coords"], grid_point_3D["beads"]):
                        atom_nr += 1
                        atom_count += 1
                        if atom_nr >= 100000:
                            atom_nr -= 100000 * (atom_nr // 100000)
                        x = x + (self.pbc_box[0] / 2)
                        y = y + (self.pbc_box[1] / 2)
                        z = z + (self.pbc_box[2] / 2)
                        a_name = bead_name
                        string = ""
                        if self.output_system_file_name.endswith("pdb"):
                            output_system_file_lines.append(self.pdb_atom_writer("ATOM", atom_nr, a_name, " ", r_name, "A", res_nr, " ", float(x), float(y), float(z), float(1), float(0), " ", " ", " "))
                        elif self.output_system_file_name.endswith("gro"): ### gro coordinates are in [nm] not [Å]
                            output_system_file_lines.append(self.gro_atom_writer(res_nr, r_name, a_name, atom_nr, x / 10, y / 10, z / 10, " ", " ", " "))

        ### ### Creating end of file lines
        if self.output_system_file_name.endswith("pdb"):
            output_system_file_lines.append("TER")
            output_system_file_lines.append("END")
        
        elif self.output_system_file_name.endswith("gro"):
            output_system_file_lines[1] = " " + str(atom_count)
            output_system_file_lines.append( ### gro vectors are in [nm] not [Å]
                '{v1x:>{v1xL}.5f}{v2y:>{v2yL}.5f}{v3z:>{v3zL}.5f}{v1y:>{v1yL}.5f}{v1z:>{v1zL}.5f}{v2x:>{v2xL}.5f}{v2z:>{v2zL}.5f}{v3x:>{v3xL}.5f}{v3y:>{v3yL}.5f}'.format(
                    v1x = self.pbc_box[0] / 10, v1xL = 10,
                    v2y = self.pbc_box[1] / 10, v2yL = 10,
                    v3z = self.pbc_box[2] / 10, v3zL = 10,
                    v1y = 0 / 10, v1yL = 10,
                    v1z = 0 / 10, v1zL = 10,
                    v2x = 0 / 10, v2xL = 10,
                    v2z = 0 / 10, v2zL = 10,
                    v3x = 0 / 10, v3xL = 10,
                    v3y = 0 / 10, v3yL = 10,
                )
            )
            output_system_file_lines.append("")
        
        if self.backup:
            self.backupper(self.output_system_file_name)

        new_file = open(self.output_system_file_name, "w")
        for line in output_system_file_lines:
            new_file.write(line + "\n")
        new_file.close()

        self.print_term("Structure file written", "\n")

    def topol_file_writer(self):
        if self.output_topol_file_name:
            self.print_term("Writing topology file:", self.output_topol_file_name, "\n")
            output_topol_file_lines = [
                "",
                "[ system ]",
                "; name",
                self.system_name or "PLACEHOLDER_TITLE",
                "",
                "[ molecules ]",
                "; name number",
            ]
            molecules_for_top_lengths = [len(str(max(data, key = lambda d: len(str(d))))) for data in zip(*self.molecules_for_top)]
            for n, c in self.molecules_for_top:
                output_topol_file_lines.append(
                    '{NAME:<{Ln}} {COUNT:>{Lc}}'.format(
                        NAME = n, Ln = molecules_for_top_lengths[0],
                        COUNT = c, Lc = molecules_for_top_lengths[1],
                    )
                )
            if self.backup:
                self.backupper(self.output_topol_file_name)

            new_file = open(self.output_topol_file_name, "w")
            for line in output_topol_file_lines:
                new_file.write(line + "\n")
            new_file.close()
            self.print_term("Topology file written", "\n")
    
    def log_file_writer(self):
        if self.output_log_file_name:
            self.print_term("Writing log file file:", self.output_log_file_name, "\n")
            if self.backup:
                self.backupper(self.output_log_file_name)
            new_file = open(self.output_log_file_name, "w")
            for line in self.LOG_FILE:
                new_file.write(line + "\n")
            new_file.close()
            self.print_term("Log file written", "\n")
    
    #########################################
    ### Lipid/solv/ion defs preprocessors ###
    #########################################
    def lipid_defs_preprocessor(self):
        '''
        Preprocesses lipid defintions, by rearranging the data into a more useful system
        '''
        self.print_term("Preprocessing lipid definitions")
        lipid_dict = {}
        ### Loop over generel lipid types (e.g. phospholipids, sterols, etc.)
        for (lipid_type, ff), lipid_type_dict in self.lipid_defs.items():
            if ff not in lipid_dict.keys():
                lipid_dict[ff] = {}
            x, y, z = lipid_type_dict["x"], lipid_type_dict["y"], lipid_type_dict["z"]

            ### Centering x/y/z-values for current general lipid type
            if "center" in lipid_type_dict.keys():
                lcz = lipid_type_dict["center"]
            else:
                lcz = 0

            ### bead distance adjustment values for current general lipid type
            if "bd" in lipid_type_dict.keys():
                bdx, bdy, bdz = lipid_type_dict["bd"]
            else:
                bdx, bdy, bdz = 0.25, 0.25, 0.3

            ### Loop over indexes
            for (lipid_name, data_specifier), lipid_details in lipid_type_dict["lipids"].items():
                if lipid_name not in lipid_dict[ff].keys():
                    lipid_dict[ff][lipid_name] = {
    #                     "x": [(xi - lcx) * bdx * 10 for xi in x],
    #                     "y": [(yi - lcy) * bdy * 10 for yi in y],
    #                     "z": [(zi - lcz) * bdz * 10 for zi in z],
    #                     "charges": lipid_type_charges
    #                     "conformations": [],
                    }

    #                 if len(imp_cmds) != []:
    #                     if lipid_name in imp_structs:
    #                         lipid_dict[ff][lipid_name]["conformations"] = imp_structs[lipid_name]["conformations"]
    #                         for confi, vals in enumerate(lipid_dict[ff][lipid_name]["conformations"]):
    #                             ### Aligning conformations somewhat
    #                             xs, ys, zs = vals["x"], vals["y"], vals["z"]
    #                             headi = 0
    #                             hx, hy, hz = xs[headi], ys[headi], zs[headi]
    #                             mxs = [i - hx for i in xs]
    #                             mys = [i - hy for i in ys]
    #                             mzs = [i - hz for i in zs]
    #                             cx, cy, cz = np.mean(xs), np.mean(ys), np.mean(zs)
    #                             vector = [hx - cx, hy - cy, hz - cz]
    #                             vertical = [0, 0, 1]
    #                             rot_matrix = rotation_matrix_from_vectors(vector, vertical)
    #                             nxs, nys, nzs = list(zip(*[rot_matrix.dot([x, y, z]) for x, y, z in zip(mxs, mys, mzs)]))
    #                             ncx, ncy, ncz = np.mean(nxs), np.mean(nys), np.mean(nzs)
    #                             lipid_dict[ff][lipid_name]["conformations"][confi]["x"] = [i - ncx for i in nxs]
    #                             lipid_dict[ff][lipid_name]["conformations"][confi]["y"] = [i - ncy for i in nys]
    #                             lipid_dict[ff][lipid_name]["conformations"][confi]["z"] = [i - (lcz * bdz * 10) for i in nzs]

                    if "charges" in lipid_type_dict.keys():
                        lipid_dict[ff][lipid_name]["charges"] = lipid_type_dict["charges"]

                if data_specifier == "beads":
                    lipid_dict[ff][lipid_name]["beads"] = list(filter(None, lipid_details.split(" ")))
                    ### Remove coordinate and bead indexes if no bead is assigned
                    lipid_dict[ff][lipid_name]["x"] = [xi for xi, bead in zip(x, lipid_dict[ff][lipid_name]["beads"]) if bead != "-"]
                    lipid_dict[ff][lipid_name]["y"] = [yi for yi, bead in zip(y, lipid_dict[ff][lipid_name]["beads"]) if bead != "-"]
                    lipid_dict[ff][lipid_name]["z"] = [zi for zi, bead in zip(z, lipid_dict[ff][lipid_name]["beads"]) if bead != "-"]
                    lipid_dict[ff][lipid_name]["beads"] = [bead for bead in lipid_dict[ff][lipid_name]["beads"] if bead != "-"]

                    ### Calculate lipid-spicific x/y-center
                    lcx = (max(lipid_dict[ff][lipid_name]["x"]) + min(lipid_dict[ff][lipid_name]["x"])) / 2
                    lcy = (max(lipid_dict[ff][lipid_name]["y"]) + min(lipid_dict[ff][lipid_name]["y"])) / 2

                    ### Center x/y/z and scale using bd
                    lipid_dict[ff][lipid_name]["x"] = [(xi - lcx) * bdx * 10 for xi in lipid_dict[ff][lipid_name]["x"]]
                    lipid_dict[ff][lipid_name]["y"] = [(yi - lcy) * bdy * 10 for yi in lipid_dict[ff][lipid_name]["y"]]
                    lipid_dict[ff][lipid_name]["z"] = [(zi - lcz) * bdz * 10 for zi in lipid_dict[ff][lipid_name]["z"]]

                if data_specifier == "charges":
                    lipid_dict[ff][lipid_name]["charges"] = lipid_details
        
        tot_lipids = sum([len(vals) for vals in lipid_dict.values()])
        
        self.lipid_dict = lipid_dict
        self.print_term("    ", "Number of lipids preprocessed:", tot_lipids)
    
    def center_coords(self, coords):
        ### Calculate lipid-spicific x/y-center
        coord_diff = (max(coords) + min(coords)) / 2
        centered_coords = [coord - coord_diff for coord in coords]
        return centered_coords
    
    #############################################
    ### Protein/leaflet/solvent preprocessors ###
    #############################################
    def prot_preprocessor(self):
        '''
        Preprocesses protein commands for later ease of use
        '''
        self.PROTEINS = {}
        if len(self.PROTEINS_cmds) != 0:
            self.print_term("\nPreprocessing protein requests")
            for cmd_nr, prot_cmd in enumerate(self.PROTEINS_cmds, 1):
                self.print_term("    ", "Starting protein:", cmd_nr)

                ### Defaults
                prot_dict = {
                    "tx": 0,
                    "ty": 0,
                    "tz": 0,
                    "rx": 0,
                    "ry": 0,
                    "rz": 0,
                    "cen_method": ("cog",), # "cog" (center of geometry), "axis", "bead:INT", "res:INT" or "point:x:y:z"
                    "lipids_inside": False, # [bool]
                    "pbc_check": True,
                    "buffer": 1.32, # [Å] default = (vdw of regular beads) / 2
                    "alpha_mult": 1.0,
                    "mol_names": [],
                    "charge": "top", # int/float, "top" or "auto"
                    "tot_charge": 0,
                }

                ### ### Check protein command
                for cmd in prot_cmd.split():
                    sub_cmd = cmd.split(":")

                    ### Read pdb/gro file
                    if sub_cmd[0].lower().endswith(".pdb"):
                        prot_dict["beads"] = self.pdb_reader(sub_cmd[0])

                    elif sub_cmd[0].lower().endswith(".gro"):
                        prot_dict["beads"] = self.gro_reader(sub_cmd[0])
                    
                    ### Center method "mean", "bead:INT" or "res:INT"
                    elif sub_cmd[0].lower() == "cen_method":
                        if sub_cmd[1].lower() in ["cog", "axis"]:
                            prot_dict["cen_method"] = (sub_cmd[1].lower(),)
                        elif sub_cmd[1].lower() in ["bead", "res"]:
                            if len(sub_cmd) == 3:
                                prot_dict["cen_method"] = (sub_cmd[1].lower(), ast.literal_eval(sub_cmd[2]))
                            elif len(sub_cmd) == 4:
                                prot_dict["cen_method"] = (sub_cmd[1].lower(), ast.literal_eval(sub_cmd[2]), ast.literal_eval(sub_cmd[3]))
                        elif sub_cmd[1].lower() in ["point"]:
                            prot_dict["cen_method"] = (sub_cmd[1].lower(), ast.literal_eval(sub_cmd[2]), ast.literal_eval(sub_cmd[3]), ast.literal_eval(sub_cmd[4]))

                    ### True/False whether or not lipds may be placed within protein area
                    elif sub_cmd[0].lower() == "lipids_inside":
                        prot_dict["lipids_inside"] = ast.literal_eval(sub_cmd[1])

                    ### True/False whether to check for pbc crossing
                    elif sub_cmd[0].lower() == "pbc_check":
                        prot_dict["pbc_check"] = ast.literal_eval(sub_cmd[1])

                    ### True/False whether to check for pbc crossing
                    elif sub_cmd[0].lower() == "buffer":
                        prot_dict["buffer"] = ast.literal_eval(sub_cmd[1])

                    ### Integer multiplier for radius used in alphashape function [multiplier]
                    elif sub_cmd[0].lower() == "alpha_mult":
                        prot_dict["alpha_mult"] = ast.literal_eval(sub_cmd[1])

                    ### Integer multiplier for radius used in alphashape function [multiplier]
                    elif sub_cmd[0].lower() in ["mol_names", "mol_name"]:
                        for mol_name in sub_cmd[1:]:
                            prot_dict["mol_names"].append(mol_name)

                    ### Determine charge of protein int/float or "top"
                    elif sub_cmd[0].lower() == "charge":
                        if sub_cmd[1] == "top":
                            prot_dict["charge"] = "top"
                        elif sub_cmd[1] == "auto":
                            prot_dict["charge"] = "auto"
                        else:
                            prot_dict["charge"] = ast.literal_eval(sub_cmd[1])

                    ### xyz axis translations [nm]
                    elif sub_cmd[0].lower() in ["tx", "ty", "tz"]:
                        prot_dict[sub_cmd[0].lower()] = ast.literal_eval(sub_cmd[1]) * 10

                    ### xyz axis rotations [degrees]
                    elif sub_cmd[0].lower() in ["rx", "ry", "rz"]:
                        prot_dict[sub_cmd[0].lower()] = ast.literal_eval(sub_cmd[1])

                    ### Random kick to beads x/y/z positions [nm] #[Å]
                    elif sub_cmd[0].lower() in ["kick", "kickx", "kicky", "kickz"]:
                        if sub_cmd[0].lower() == "kick":
                            prot_dict["kickx"] = ast.literal_eval(sub_cmd[1])
                            prot_dict["kicky"] = ast.literal_eval(sub_cmd[1])
                            prot_dict["kickz"] = ast.literal_eval(sub_cmd[1])
                        else:
                            prot_dict[sub_cmd[0].lower()] = ast.literal_eval(sub_cmd[1])

                    ### Errors out if unknown subcommand used, and prints the subcommand to console
                    else:
                        assert False, "Unknown subcommand given to '-prot'. The subcommand is: '" + str(cmd) + "'"

                        
                def prot_charge_finder(beads):
                    '''
                    Finds the charge of a protein based on residue names and the 'prot_defs' dictionary
                    '''
                    charge = 0
                    counted_residues = []
                    for (bead_i, bead_nr, res_nr), values in beads.items():
                        if res_nr not in counted_residues:
                            if values["res_name"] in self.prot_defs[self.sys_ff]["charges"]:
                                charge += self.prot_defs[self.sys_ff]["charges"][values["res_name"]]
                            else:
                                print(values["res_name"], "residue name not in 'prot_defs' charges for system force field:", self.sys_ff)
                            counted_residues.append(res_nr)
                    return charge
                
                ### Post-preprocessing (topology and charge determination)
                if isinstance(prot_dict["charge"], (float, int)):
                    prot_dict["tot_charge"] += prot_dict["charge"]
                    
                elif prot_dict["charge"] == "auto":
                    prot_dict["tot_charge"] += prot_charge_finder(prot_dict["beads"])
                    
                elif prot_dict["charge"] == "top":
                    if prot_dict["mol_names"] == []:
                        self.print_term("    ", "    ", "No protein names given.  Will estimate charges from residue names. Use 'prot_name' to assign protein names (name used in topology files)")
                        prot_dict["tot_charge"] = prot_charge_finder(prot_dict["beads"])
                    else:
                        for mol_name in prot_dict["mol_names"]:
                            if mol_name not in self.itp_moltypes.keys():
                                print("    ", "    ", "A protein name could not be found in your topology file(s): " +  prot_name)
                            else:
                                prot_dict["tot_charge"] += self.itp_moltypes[mol_name]["charge_sum"]

                prot_dict["mol_names"] = prot_dict["mol_names"] or "_".join(["PROT", str(cmd_nr)])

                self.PROTEINS[cmd_nr] = prot_dict.copy()
            self.print_term("    ", "Number of protein insertions preprocessed:", len(self.PROTEINS))

    def leaf_preprocessor(self):
        '''
        Preprocesses leaflet commands for later ease of use
        '''
        self.LEAFLETS = {}
        ### ### Preprocessing leaflets
        if len(self.LEAFLETS_cmds) != 0:
            self.print_term("\nPreprocessing leaflet requests")
            leaf_nr = 1
            for cmd_nr, leaf_cmd in enumerate(self.LEAFLETS_cmds, 1):
                self.print_term("    ", "Starting leaflet command:", cmd_nr)

                ### Defaults
                leaf_dict = {
                    "lipids": {}, # empty
                    "lipids_preprocessing": [], # empty
                    "shape": "rectangle", # "rectangle"
                    "x": self.pbcx, # [å]
                    "y": self.pbcy, # [å]
                    "center": [0, 0, 0], # [nm]
                    "kickx": 0.25, # [Å]
                    "kicky": 0.25, # [Å]
                    "kickz": 0.25, # [Å]
                    "apl": 0.6 * 100, # [nm^2] converted to [Å^2]
                    "planeWR": 1.32, # [Å] default = (vdw of regular beads) / 2
                    "heightWR": 1.32, # [Å] default = (vdw of regular beads) / 2
                    "lip_round_func": ("int", int), # int, round, math.floor or math.ceil
            #         "bdx": 0.25, # [multiplier]
            #         "bdy": 0.25, # [multiplier]
            #         "bdz": 0.3, # [multiplier]
#                     "rx": 0, # [degrees]
#                     "ry": 0, # [degrees]
#                     "rz": 0, # [degrees]
                    "pbc_check": True, # [bool]
                    "lipid_optim": 'avg_optimal', # str: 'abs_val', 'force_fill', 'fill', 'avg_optimal', 'no'
                    "params": False, # False or str
                    "charge": "top", # "lib" or "top"
                }
                
                # Mono as explicitly upper added by request
                monolayer_upper_designation = ["u", "up", "upper", "m", "mo", "mono", "monolayer"]
                monolayer_lower_designation = ["d", "do", "down", "l", "lo", "lower"]
                bilayer_designation = ["b", "bi", "bilayer"]
                
                layer_defition = "bilayer"

                ### ### Check leaf command
                ### Liberal use of ast.literal_eval() below to convert strings to int/float
                ### ast.literal_eval checks if code is a valid python datatype before interpreting it
            #     for cmd in leaf_cmd.split(","):
                for cmd in leaf_cmd.split():
                    sub_cmd = cmd.split(":")

                    ### Bilayer/monolayer defition
                    if sub_cmd[0].lower() == "type":
                        if sub_cmd[1].lower() in monolayer_upper_designation:
                            layer_defition = "upper"
                        if sub_cmd[1].lower() in monolayer_lower_designation:
                            layer_defition = "lower"
                        if sub_cmd[1].lower() in bilayer_designation:
                            layer_defition = "bilayer"
                    
                    ### Replaced by the above
#                     if sub_cmd[0].lower() in ["u", "up", "upper", "d", "do", "down", "l", "lo", "lower", "b", "bi", "bilayer"]:
#                         layer_defition = sub_cmd[0].lower()

                    ### Leaflet shape
                    elif sub_cmd[0].lower() == "shape":
                        leaf_dict["shape"] = sub_cmd[1]

                    ### Center [nm]
                    elif any(sub_cmd[0].lower() == cen for cen in ["cen", "center"]):
                        leaf_dict["center"] = [ast.literal_eval(i) * 10 for i in sub_cmd[1:]]

                    ### Area per lipid [nm^2] converted to [Å^2]
                    elif sub_cmd[0].lower() == "apl":
                        leaf_dict["apl"] = ast.literal_eval(sub_cmd[1]) * 100

                    ### Wiggle room. Minimum distance from a bead to the edge of a bin [Å]
                    elif sub_cmd[0].lower() == "plane_wr":
                        leaf_dict["planeWR"] = ast.literal_eval(sub_cmd[1])

                    elif sub_cmd[0].lower() == "height_wr":
                        leaf_dict["heightWR"] = ast.literal_eval(sub_cmd[1])

                    ### Rounding method for rounding number of lipids in leaflet [function]
                    elif sub_cmd[0].lower() == "round": # "int", "round", "min" or "max"
                        rounding_types = {"int": ("int", int), "round": ("round", round), "floor": ("math.floor", math.floor), "ceil": ("math.ceil", math.ceil)}
                        leaf_dict["lip_round_func"] = rounding_types[sub_cmd[1]] 

                    elif sub_cmd[0].lower() in ["kick", "kickx", "kicky", "kickz"]:
                        ### Random kick to beads x/y/z positions [nm] #[Å]
                        if sub_cmd[0].lower() == "kick":
                            leaf_dict["kickx"] = ast.literal_eval(sub_cmd[1])
                            leaf_dict["kicky"] = ast.literal_eval(sub_cmd[1])
                            leaf_dict["kickz"] = ast.literal_eval(sub_cmd[1])
                        else:
                            leaf_dict[sub_cmd[0].lower()] = ast.literal_eval(sub_cmd[1])

                    ### Bead distance scaling for z is 0.3 by default and 0.25 by default for x and y [multiplier]
                    elif sub_cmd[0].lower() in ["bdx", "bdy", "bdz"]:
                        leaf_dict[sub_cmd[0].lower()] = ast.literal_eval(sub_cmd[1])

                    ### ### Rectangle specific subcommands:
                    ### xy dimensions of leaflet [nm]. Converted to [Å]
                    elif sub_cmd[0].lower() in ["x", "y"]:
                        leaf_dict[sub_cmd[0].lower()] = ast.literal_eval(sub_cmd[1]) * 10

                    ### xyz axis rotations [degrees]
                    elif sub_cmd[0].lower() in ["rx", "ry", "rz"]:
                        leaf_dict[sub_cmd[0].lower()] = ast.literal_eval(sub_cmd[1])

                    ### True/False whether to check for pbc crossing
                    elif sub_cmd[0].lower() == "pbc_check":
                        leaf_dict["pbc_check"] = ast.literal_eval(sub_cmd[1])

                    ### Lipid optimization method
                    elif sub_cmd[0].lower() == "lipid_optim":
                        valid_lipid_optim = sub_cmd[1] in ['avg_optimal', 'abs_val', 'force_fill', 'fill', 'no']
                        assert valid_lipid_optim, "Invalid lipid selection method: '" + str(sub_cmd[1]) + "'"
                        leaf_dict["lipid_optim"] = sub_cmd[1]

                    ### Designates which force field to collect lipids from
                    elif sub_cmd[0].lower() == "params":
                        leaf_dict["params"] = ast.literal_eval(sub_cmd[1])

                    ### True/False whether to check for pbc crossing
                    elif sub_cmd[0].lower() == "charge":
                        leaf_dict["charge"] = sub_cmd[1]

                    ### Processes only lipids in lipid dictionary
                    elif any([sub_cmd[0] in self.lipid_dict[ff].keys() for ff in self.lipid_dict.keys()]):
                        leaf_dict["lipids_preprocessing"].append(sub_cmd)

                    ### Errors out if unknown subcommand used, and prints the subcommand to terminal
                    else:
                        assert False, "Unknown subcommand given to '-leaf'. The subcommand is: '" + str(cmd) + "'"

                ################################
                ### Lipid data incorporation ###
                ################################
                ### Reconfigures lipid-specific data for leaflet according to subcommands
                assert leaf_dict["lipids_preprocessing"] != {}, "No lipids given to '-leaf' flag. Please specify at least one lipid if you use it."

                for l_name, *rest in leaf_dict["lipids_preprocessing"]:
                    """
                    Checks if parameters specified for lipid.
                    Order: lipid specific params, leaflet specific params, general system params (defaults to M3)
                    """
                    if len(rest) != 1:
                        l_ratio, l_params = rest
                    else:
                        l_ratio = rest[0]
                        l_params = leaf_dict["params"] or self.sys_params # (sys_ff defaults to M3)

                    cur_lipid = self.lipid_dict[l_params][l_name]
                    leaf_dict["lipids"][l_name] = {
                        "ratio": ast.literal_eval(l_ratio),
                        "beads": cur_lipid["beads"],
                        "x"    : cur_lipid["x"],
                        "y"    : cur_lipid["y"],
                        "z"    : cur_lipid["z"],
                        "zmax" : max(cur_lipid["z"]),
                        "zmin" : min(cur_lipid["z"]),
                    }

                    ### Charges
                    if leaf_dict["charge"] == "top":
                        if l_name in self.itp_moltypes.keys():
                            leaf_dict["lipids"][l_name]["charge_sum"] = self.itp_moltypes[l_name]["charge_sum"]
                        elif "charges" in leaf_dict["lipids"][l_name].keys():
                            leaf_dict["lipids"][l_name]["charge_sum"] = sum([charge for bead, charge in cur_lipid["charges"]])
                        else:
                            leaf_dict["lipids"][l_name]["charge_sum"] = 0
                    elif leaf_dict["charge"] == "lib" and "charges" in leaf_dict["lipids"][l_name].keys():
                        leaf_dict["lipids"][l_name]["charge_sum"] = sum([charge for bead, charge in cur_lipid["charges"]])
                    else:
                        leaf_dict["lipids"][l_name]["charge_sum"] = 0

                ### Upper leaflet monolayer
                if layer_defition in monolayer_upper_designation:
                    leaf_dict["HG_direction"] = "up"
                    self.LEAFLETS[leaf_nr] = leaf_dict
                    leaf_nr += 1
                ### Lower leaflet monolayer
                elif layer_defition in monolayer_lower_designation:
                    leaf_dict["HG_direction"] = "down"
                    self.LEAFLETS[leaf_nr] = leaf_dict
                    leaf_nr += 1
                ### Bilayer
                elif layer_defition in bilayer_designation:
                    leaf_dict["HG_direction"] = "up"
                    self.LEAFLETS[leaf_nr] = leaf_dict.copy()
                    leaf_nr += 1
                    leaf_dict["HG_direction"] = "down"
                    self.LEAFLETS[leaf_nr] = leaf_dict.copy()
                    leaf_nr += 1
            
            self.print_term("    ", "Number of leaflets preprocessed:", len(self.LEAFLETS), "(from", cmd_nr, "commands)")

    def solv_preprocessor(self):
        '''
        Preprocesses solvation commands for later ease of use
        '''
        self.SOLVATIONS = {}
        
        if len(self.SOLVATIONS_cmds) != 0:
            self.print_term("\nPreprocessing solvent requests")
            for cmd_nr, solvate_cmd in enumerate(self.SOLVATIONS_cmds, 1):
                self.print_term("    ", "Starting Solvent command:", cmd_nr)
                ### Defaults
                solv_dict = {
                    ### ### Solvent
                    "solvent": {}, # empty
                    "solv_preprocessing": [], # empty
            #         "solvcount": False, # [bool] whether to consider solvcount or not
                    "solv_molarity": 55.56, # [int] or [float] # mol/L
                    "solvfreevol": True,
                    ### molarity of water = 55.56, thus molarity of regular bead water should be a fourth of that (13.89)

                    ### ### Ions
                    "neg_ions": {}, # empty
                    "pos_ions": {}, # empty
                    "neg_ions_preprocessing": [], # empty
                    "pos_ions_preprocessing": [], # empty
            #         "ionscount": False, # [bool] whether to consider solvcount or not
                    "charge": 0, # int or float # Target charge for system
                    "sys_charge": True, # [bool]
                    "salt_molarity": 0.15, # int or float # mol/L
                    "ionsvol": "solv", # "box", "free", "solv"
                    ### starting mol/liter concentration of both negative and positive ions (each will have the concentration) 

                    ### ### General
                    "count": False, # Uses specific molarity value as absolute number of molecules instead of molarity ratio. Will be rounded using int(val+0.5)
                    "kick": 0.5, # [Å]
                    "bdx": 1.0, # [multiplier]
                    "bdy": 1.0, # [multiplier]
                    "bdz": 1.0, # [multiplier]
                    "params": False, # False or str
                    "bead_radius": 2.64, # [Å]
                    "gridres": 1.2, # [Å]
                    "WR": 2.64 * 1, # [Å]
                    "buffer": 2.0, # [Å]
                    "min_cellsize": 1, # Integer value
                }
                ### ### Check leaf command
                ### Liberal use of ast.literal_eval() below to convert strings to int/float
                ### ast.literal_eval checks if code is a valid python datatype before interpreting it
            #     for cmd in solvate_cmd.split(","):
                if solvate_cmd == "":
                    solvate_cmd = "solv:W pos:NA neg:CL"
                for cmd in solvate_cmd.split():
                    sub_cmd = cmd.split(":")
                    ########################
                    ### SOLVENT SPECIFIC ###
                    ########################
                    ### Solvents
                    if sub_cmd[0].lower() == "solv":
                        solv_dict["solv_preprocessing"].append(sub_cmd[1:])

                    ### Solvent concentration
                    elif sub_cmd[0].lower() == "solv_molarity":
                        solv_dict["solv_molarity"] = ast.literal_eval(sub_cmd[1])

                    ### Solvent count
                    elif sub_cmd[0].lower() == "solvcount":
                        solv_dict["solvcount"] = ast.literal_eval(sub_cmd[1])

                    ### Whether to consider free volume (excluding lipid and protein) or purely the box volume
                    elif sub_cmd[0].lower() == "solvfreevol":
                        solv_dict["solvfreevol"] = ast.literal_eval(sub_cmd[1])

                    ####################
                    ### ION SPECIFIC ###
                    ####################
                    ### Negative ions
                    elif sub_cmd[0].lower() == "neg":
                        solv_dict["neg_ions_preprocessing"].append(sub_cmd[1:])

                    ### Positive ions
                    elif sub_cmd[0].lower() == "pos":
                        solv_dict["pos_ions_preprocessing"].append(sub_cmd[1:])

                    ### Solvent count
                    elif sub_cmd[0].lower() == "ionscount":
                        solv_dict["ionscount"] = ast.literal_eval(sub_cmd[1])

                    ### Target charge [int/float] or False for no "neutralization"
                    elif sub_cmd[0].lower() == "charge":
                        solv_dict["charge"] = ast.literal_eval(sub_cmd[1])

                    ### Whether to consider system charge for charge calculations
                    elif sub_cmd[0].lower() == "sys_charge":
                        solv_dict["sys_charge"] = ast.literal_eval(sub_cmd[1])

                    ### Salt concentration
                    elif sub_cmd[0].lower() == "salt_molarity":
                        solv_dict["salt_molarity"] = ast.literal_eval(sub_cmd[1])

                    ### Whether to consider free volume (excluding lipid and protein) or purely the box volume
                    elif sub_cmd[0].lower() == "ionsvol":
                        solv_dict["ionsvol"] = sub_cmd[1]

                    ###############
                    ### GENERAL ###
                    ###############
                    ### Whether to use ratios as absolute number of molecules
                    elif sub_cmd[0].lower() == "count":
                        solv_dict["count"] = ast.literal_eval(sub_cmd[1])

                    ### Radius of non-solvent beads for volume calculations
                    elif sub_cmd[0].lower() == "bead_radius":
                        solv_dict["bead_radius"] = ast.literal_eval(sub_cmd[1])

                    ### Grid resolution
                    elif sub_cmd[0].lower() == "gridres":
                        solv_dict["gridres"] = ast.literal_eval(sub_cmd[1])

                    ### Minimum wiggle room for beads
                    elif sub_cmd[0].lower() == "wr":
                        solv_dict["WR"] = ast.literal_eval(sub_cmd[1])

                    elif sub_cmd[0].lower() in ["kick"]:
                        ### Random kick to beads x/y/z positions [Å]
                        solv_dict["kick"] = ast.literal_eval(sub_cmd[1])
                        ### Bead distance scaling for z is 0.3 by default and 0.25 by default for x and y [multiplier]

                    ### Bead scaling distances
                    elif sub_cmd[0].lower() in ["bdx", "bdy", "bdz"]:
                        leaf_dict[sub_cmd[0].lower()] = ast.literal_eval(sub_cmd[1])

                    ### Designates force field
                    elif sub_cmd[0].lower() == "params":
                        solv_dict["params"] = sub_cmd[1]

                    ### Buffer for solvent placement from leaflet hydrophobic area
                    elif sub_cmd[0].lower() == "buffer":
                        solv_dict["buffer"] = ast.literal_eval(sub_cmd[1])

                    ### Buffer for solvent placement from leaflet hydrophobic area
                    elif sub_cmd[0].lower() == "min_cellsize":
                        solv_dict["min_cellsize"] = int(ast.literal_eval(sub_cmd[1]) + 0.5)

                    ### Errors out if unknown subcommand used, and prints the subcommand to terminal
                    else:
                        assert False, "Unknown subcommand given to '-solvate'. The subcommand is: '" + str(cmd) + "'"
                
                ### Changed such that the following is done before the command is preprocessed
#                 if solv_dict["solv_preprocessing"] == {} and solv_dict["pos_ions_preprocessing"] == {} and solv_dict["neg_ions_preprocessing"] == {}:
#                     solv_dict["solv_preprocessing"].append("solv:W")
#                     solv_dict["pos_ions_preprocessing"].append("solv:NA")
#                     solv_dict["neg_ions_preprocessing"].append("solv:CL")
                    
                ######################################
                ### SOLVENT/ION DATA INCORPORATION ###
                ######################################
                ### Reconfigures lipid-specific data for leaflet according to subcommands
                assert solv_dict["solv_preprocessing"] != {}, "No solvent given to '-solvate' flag. Please specify at least one non-ionic solvent if you use it."

                params = solv_dict["params"] or self.sys_params # (sys_ff defaults to M3)
                ### ### SOLVENT DATA
                for name, *rest in solv_dict["solv_preprocessing"]:
                    ### Default values for solvent
                    solv_dict_default = {
                        "charge": 0,
                        "ratio": 1,
                        "solvcount": 1,
                        "molarity": solv_dict["solv_molarity"],
                    }

                    ### New concentration
                    if len(rest) != 0:
                        solv_dict_default["molarity"] = ast.literal_eval(rest[0])

                    ### If count has been set, then convert to integer value and treat as absolute number of molecules
                    if solv_dict["count"]:
                        solv_dict_default["molarity"] = int(solv_dict_default["molarity"] + 0.5)

                    ### Update with solven library defaults
#                     print(self.solvent_defs)
                    solv_dict_default.update(self.solvent_defs[params][name])

                    for ax in ["x", "y", "z"]:
                        ### Axis geometric center
                        lc = (max(solv_dict_default[ax]) + min(solv_dict_default[ax])) / 2
                        ### Adjust internal coordinate in accordance with axis geometric center and axis bead-distance (bd) scaling
                        solv_dict_default[ax] = [(i - lc) * solv_dict["bd" + ax] * 10 for i in solv_dict_default[ax]]

                    ### Ensuring that beads are stored as a tuple for consistency
                    if type(solv_dict_default["beads"]) == str:
                        solv_dict_default["beads"] = (solv_dict_default["beads"],)

                    ### Adding name:dict combo to solvent dict
                    solv_dict["solvent"][name] = solv_dict_default

                ### ### Positive ions
                for name, *rest in solv_dict["pos_ions_preprocessing"]:
                    ### Default values for solvent
                    solv_dict_default = {
                        "charge": 1,
                        "ratio": 1,
                        "solvcount": 1,
                        "molarity": solv_dict["salt_molarity"],
                    }

                    ### New ratio
                    if len(rest) != 0:
                        solv_dict_default["molarity"] = ast.literal_eval(rest[0])

                    ### If count has been set, then convert to integer value and treat as absolute number of molecules
                    if solv_dict["count"]:
                        solv_dict_default["molarity"] = int(solv_dict_default["molarity"] + 0.5)

                    ### Update with solven library defaults
                    solv_dict_default.update(self.ion_defs[params]["positive"][name])

                    for ax in ["x", "y", "z"]:
                        ### Axis geometric center
                        lc = (max(solv_dict_default[ax]) + min(solv_dict_default[ax])) / 2
                        ### Adjust internal coordinate in accordance with axis geometric center and axis bead-distance (bd) scaling
                        solv_dict_default[ax] = [(i - lc) * solv_dict["bd" + ax] * 10 for i in solv_dict_default[ax]]

                    ### Ensuring that beads is a tuple for consistency
                    if type(solv_dict_default["beads"]) == str:
                        solv_dict_default["beads"] = (solv_dict_default["beads"],)

                    ### Adding name:dict combo to solvent dict
                    solv_dict["pos_ions"][name] = solv_dict_default

                ### ### Negative ions
                for name, *rest in solv_dict["neg_ions_preprocessing"]:
                    ### Default values for solvent
                    solv_dict_default = {
                        "charge": -1,
                        "ratio": 1,
                        "solvcount": 1,
                        "molarity": solv_dict["salt_molarity"],
                    }

                    ### New ratio
                    if len(rest) != 0:
                        solv_dict_default["molarity"] = ast.literal_eval(rest[0])

                    ### If count has been set, then convert to integer value and treat as absolute number of molecules
                    if solv_dict["count"]:
                        solv_dict_default["molarity"] = int(solv_dict_default["molarity"] + 0.5)

                    ### Update with solven library defaults
                    solv_dict_default.update(self.ion_defs[params]["negative"][name])

                    for ax in ["x", "y", "z"]:
                        ### Axis geometric center
                        lc = (max(solv_dict_default[ax]) + min(solv_dict_default[ax])) / 2
                        ### Adjust internal coordinate in accordance with axis geometric center and axis bead-distance (bd) scaling
                        solv_dict_default[ax] = [(i - lc) * solv_dict["bd" + ax] * 10 for i in solv_dict_default[ax]]

                    ### Ensuring that beads is a tuple for consistency
                    if type(solv_dict_default["beads"]) == str:
                        solv_dict_default["beads"] = (solv_dict_default["beads"],)

                    ### Adding name:dict combo to solvent dict
                    solv_dict["neg_ions"][name] = solv_dict_default
                    
                self.SOLVATIONS[cmd_nr] = solv_dict
                
            self.print_term("    ", "Number of solvent commands preprocessed:", len(self.SOLVATIONS))
    
    def itp_read_initiater(self):
        if len(self.ITP_INPUT_cmds) != 0:
            self.print_term("\nLoading topology file(s)")
            self.itp_defs = {
                "bondtypes": {},
                "angletypes": {},
                "dihedraltypes": {},
            }
            self.itp_moltypes = {}
            for itp_i, itp_cmd in enumerate(self.ITP_INPUT_cmds, 0):
                itp_file = itp_cmd.split()[0]
                self.itp_reader(itp_file)
            for moltype, vals in self.itp_moltypes.items():
                charge_sum = 0
                for atom, atom_vals in vals["atoms"].items():
                    if "charge" in atom_vals:
                        charge_sum += ast.literal_eval(atom_vals["charge"])
                self.itp_moltypes[moltype]["charge_sum"] = charge_sum

            self.print_term("    ", "Finished loading topologies. Number of moleculetypes found:", len(self.itp_moltypes))
    
    def import_structures_handler(self):
        if len(self.SOLUTE_INPUT_cmds) != 0:
            for imp_struc in self.SOLUTE_INPUT_cmds:
                structures = []
                
                imp_dict = {
                    "structures": [], # empty
                    "params": "M3", # str
                    "charge": [], # "top" or int or float
                    "structs": "multiple", # multiple/mult/m or single/s
                    ### "single/s" currently not implemented
#                     "top_name": False, # Only relevant if "structs" is single
                }
                
                for cmd in imp_struc.split():
                    sub_cmd = cmd.split(":")
                    if sub_cmd[0].endswith("pdb"):
                        imp_dict["structures"].append(self.pdb_reader(sub_cmd[0]))
                    if sub_cmd[0].endswith("gro"):
                        imp_dict["structures"].append(self.gro_reader(sub_cmd[0]))
                    
                    if sub_cmd[0] == "params":
                        imp_dict["params"] = sub_cmd[1]
                    
                    if sub_cmd[0] == "charge":
                        for charge in sub_cmd[1:]:
                            isnum, isint = self.isnumnber(charge)
                            if sub_cmd[0] == "top":
                                imp_dict["charge"].append(charge)
                            elif isnum:
                                imp_dict["charge"].append(ast.literal_eval(charge))
                    
                    if sub_cmd[0] == "structs":
                        if sub_cmd[1] in ["multiple", "mult", "m"]:
                            imp_dict["structs"] = "multiple"
                        if sub_cmd[1] in ["single", "s"]:
                            imp_dict["structs"] = "single"
                    
#                     if sub_cmd[0] == "top_name":
#                             imp_dict["structs"] = sub_cmd[1]
                
#                 ### Post-preprocess structures
#                 assert len(imp_dict["structures"]) > 0, "No structure files given (.pdb and .gro are accepted files)"
#                 if len(imp_dict["charge"]) == 0: ### If no charge data given
#                     imp_dict["charge"] = ["top"] * len(imp_dict["structures"])
#                 elif len(imp_dict["charge"]) == 1 and imp_dict["charge"][0] == "top": ### All should be topology based
#                     imp_dict["charge"] = ["top"] * len(imp_dict["structures"])
#                 elif len(imp_dict["charge"]) > 0 and not len(imp_dict["charge"]) == len(imp_dict["structures"]): ### Incorrect number of charges given
#                     assert False, "Incorrect number of charges [" + str(len(imp_dict["charge"])) +"] supplied compared to number of structures [" + str(len(imp_dict["charge"])) +"]"
                ### Post-preprocess structures
                assert len(imp_dict["structures"]) > 0, "No structure files given (.pdb and .gro are accepted files)"
                if len(imp_dict["charge"]) == 0: ### If no charge data given
                    imp_dict["charge"] = "top"
                elif len(imp_dict["charge"]) == 1 and imp_dict["charge"][0] == "top": ### All should be topology based
                    imp_dict["charge"] = "top"
                
                def mol_dict_maker():
                    if resname and resnumber: ### Prevents first blank set of data from being written as solvent
                        if type(imp_dict["charge"]) == list:
                            charge_method = imp_dict["charge"][rescounter]
                        else:
                            charge_method = imp_dict["charge"]
                        
                        if charge_method == "top":
                            if resname in self.itp_moltypes.keys():
                                charge = self.itp_moltypes[resname]["charge_sum"]
                            else:
                                self.print_term("    ", "    ", "WARNING: resname not found in topology. Setting charge to 0 for resname:", resname)
                                charge = 0
                        else:
                            charge = imp_dict["charge"][rescounter]
                        
                        mol_dict = {
                            "beads": tuple(beads),
                            "solvcount": solvcount,
                            "x": tuple(self.center_coords(xs)),
                            "y": tuple(self.center_coords(ys)),
                            "z": tuple(self.center_coords(zs)),
                            "charge": charge,
                        }
                        
                        ### Adds force field to solvent and ion defs if they are not present
                        if imp_dict["params"] not in self.solvent_defs.keys():
                            self.solvent_defs[imp_dict["params"]] = {}
                        if imp_dict["params"] not in self.ion_defs.keys():
                            self.ion_defs[imp_dict["params"]] = {}

                        self.solvent_defs[imp_dict["params"]][resname] = mol_dict
                        self.ion_defs[imp_dict["params"]][resname] = mol_dict
                
                rescounter = 0
                for struc_i, structure in enumerate(imp_dict["structures"]):
                    resname = False
                    resnumber = False
                    solvcount = 1
                    ### If each residue should be treated as a separate molecule
                    if imp_dict["structs"] == "multiple":
                        for a_nr, (key, atom) in enumerate(structure.items()):
                            if resname != atom["res_name"] and resnumber != atom["res_nr"]:
                                
                                mol_dict_maker()
    
                                resname = atom["res_name"]
                                resnumber = atom["res_nr"]
                                solvcount = 1
                                beads = []
                                xs = []
                                ys = []
                                zs = []
                                rescounter += 1

                            beads.append(atom["atom_name"])
                            ### Divide by 10 because values are multiplied by 10 in solv_preprocessor()
                            ### Should maybe change that
                            xs.append(atom["x"] / 10)
                            ys.append(atom["y"] / 10)
                            zs.append(atom["z"] / 10)
                        mol_dict_maker()

    #########################
    ### GENERAL FUNCTIONS ###
    #########################
    def rotate_point(self, posx, posy, posz, x_ang_deg, y_ang_deg, z_ang_deg):
        '''
        https://chat.openai.com
        '''
        # Convert angles to radians
        x_ang = math.radians(x_ang_deg)
        y_ang = math.radians(y_ang_deg)
        z_ang = math.radians(z_ang_deg)

        # Calculate rotation matrices
        rx = [[1, 0, 0]       , [0, math.cos(x_ang) , -math.sin(x_ang)], [0, math.sin(x_ang), math.cos(x_ang)]]
        ry = [[math.cos(y_ang), 0, math.sin(y_ang)] , [0, 1, 0]        , [-math.sin(y_ang)  , 0, math.cos(y_ang)]]
        rz = [[math.cos(z_ang), -math.sin(z_ang), 0], [math.sin(z_ang) , math.cos(z_ang), 0], [0, 0, 1]]

        # Combine rotation matrices
        rxy  = [[sum(a*b for a,b in zip(row,col)) for col in zip(*ry)] for row in rx]
        rxyz = [[sum(a*b for a,b in zip(row,col)) for col in zip(*rz)] for row in rxy]

        # Calculate new position of point
        new_posx = rxyz[0][0] * posx + rxyz[0][1] * posy + rxyz[0][2] * posz
        new_posy = rxyz[1][0] * posx + rxyz[1][1] * posy + rxyz[1][2] * posz
        new_posz = rxyz[2][0] * posx + rxyz[2][1] * posy + rxyz[2][2] * posz

        return new_posx, new_posy, new_posz

    def rotation_matrix_from_vectors(self, vec1, vec2):
        """
        https://stackoverflow.com/questions/45142959/calculate-rotation-matrix-to-align-two-vectors-in-3d-space
        Find the rotation matrix that aligns vec1 to vec2
        :param vec1: A 3d "source" vector
        :param vec2: A 3d "destination" vector
        :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
        """
        a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
        v = np.cross(a, b)
        c = np.dot(a, b)
        s = np.linalg.norm(v)
        kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
        return rotation_matrix

    def n_list_mixer(self, *Lists):
        ### Mixes n lists to achieve even spacing between values from initial lists
        mix_determiner = []
        for l in Lists:
            mix_determiner.append([(i / len(l), v) for i, v in enumerate(l)])
        mixed = [v for i, v in sorted(list(itertools.chain(*mix_determiner)))]
        return mixed

    def coord_checker(self, positions, axes, error_count = True):
            '''
            Checks if coordinates are within the box and moves them inside if they are not
            '''
            count = 0
            if type(positions) != list:
                positions = [positions]
            if type(axes) != list:
                axes = [axes]

            if len(positions) == len(axes):
                for pi, (pos, ax) in enumerate(zip(positions, axes)):
                    if pos > ax / 2:
                        positions[pi] -= ax
                        if error_count:
                            count += 1
                        else:
                            self.print_term("    ", "    ", "WARNING: Points are outside pbc. Moved to other side. Expect potential overlap of lipids", warn = True)
                        
                    elif pos < -ax / 2:
                        positions[pi] += ax
                        if error_count:
                            count += 1
                        else:
                            self.print_term("    ", "    ", "WARNING: Points are outside pbc. Moved to other side. Expect potential overlap of lipids", warn = True)
            else:
                self.print_term("    ", "    ", "Incorrect dimensions for 'positions' and 'axes' in in 'coord_checker' function")
                self.print_term("    ", "    ", "'positions' length:", len(positions), "'axes' length:", len(axes))
            
            if error_count:
                return positions, count
            else:
                return positions

    def fix_points(self, grid, leaf_cmd):
        #################
        ### CENTERING ###
        #################
        for grid_point_nr in grid.keys():
            grid[grid_point_nr]["x"] -= 0.5 * leaf_cmd["x"]
            grid[grid_point_nr]["y"] -= 0.5 * leaf_cmd["y"]
            grid[grid_point_nr]["z"] = 0
            grid[grid_point_nr]["internal_z"] = 0

        #################
        ### ROTATIONS ###
        #################
        ### No rotations around x/y axis as of right now
#         x_ang_deg, y_ang_deg, z_ang_deg = leaf_cmd["rx"], leaf_cmd["ry"], leaf_cmd["rz"]
#         x_ang_deg, y_ang_deg, z_ang_deg = 0, 0, leaf_cmd["rz"]
#         if any([ang != 0 for ang in [x_ang_deg, y_ang_deg, z_ang_deg]]):
#             for grid_point_nr in grid.keys():
#                 grid[grid_point_nr]["x"], grid[grid_point_nr]["y"], grid[grid_point_nr]["z"] = rotate_point(
#                     grid[grid_point_nr]["x"],
#                     grid[grid_point_nr]["y"],
#                     grid[grid_point_nr]["z"],
#                     x_ang_deg,
#                     y_ang_deg,
#                     z_ang_deg,
#                 )
#                 grid[grid_point_nr]["x_ang_deg"], grid[grid_point_nr]["y_ang_deg"], grid[grid_point_nr]["z_ang_deg"] = x_ang_deg, y_ang_deg, z_ang_deg

        ####################
        ### TRANSLATIONS ###
        ####################
        cx, cy, cz = leaf_cmd["center"]
        if any([ang != 0 for ang in [cx, cy, cz]]):
            for grid_point_nr in grid.keys():
                grid[grid_point_nr]["x"] += cx
                grid[grid_point_nr]["y"] += cy
                grid[grid_point_nr]["z"] += cz

        #################################################
        ### CHECKS IF COORDINATES ARE OUTSIDE THE BOX ###
        #################################################
        if leaf_cmd["pbc_check"]:
            errors_count = 0
            for grid_point_nr in grid.keys():
                bead_coords = [grid[grid_point_nr]["x"], grid[grid_point_nr]["y"], grid[grid_point_nr]["z"]]
                checked_beads, error = self.coord_checker(bead_coords, self.pbc_box, error_count = True)
                if error > 0:
                    errors_count += 1
                    grid[grid_point_nr]["x"] = bead_coords[0]
                    grid[grid_point_nr]["y"] = bead_coords[1]
                    grid[grid_point_nr]["z"] = bead_coords[2]

        return grid

    def backupper(self, output_file_name):
        """
        Checks if output file already exists and backs it up
        """
        output_file_split = output_file_name.split("/")
        output_path = ""
        output_name = output_file_split[-1]
        if len(output_file_split) > 1:
            for i in range(len(output_file_split) - 1):
                output_path += output_file_split[i] + "/"
        if os.path.exists(output_file_name):
            self.print_term("File " + output_file_name + " already exists. Backing it up")
            number = 1
            while True:
                if os.path.exists(output_path + "#" + output_name + "." + str(number) + "#"):
                    number += 1
                else:
                    os.rename(output_file_name, output_path + "#" + output_name + "." + str(number) + "#")
                    break

    def is_number(self, s):
        sign = "+"
        number = False
        integer = False
        if type(s) == str:
            if s.startswith("-"):
                sign, s = "-", s[1:]
        try:
            s = float(s)
            number = True
            if s.is_integer():
                integer = True
            return number, integer
        except ValueError:
            return False, False
                    
    ######################
    ### PROTEIN PLACER ###
    ######################
    def prot_placer(self):
        '''
        Places all proteins into the systems internal coordinate system
        Checks if all atoms/beads are within the pbc and moves them if they are outside
        '''
        if len(self.PROTEINS) != 0:
            self.print_term("------------------------------ PROTEIN PLACEMENT")
            for protein_nr, protein in self.PROTEINS.items():
                self.print_term("Starting protein nr", protein_nr)
                self.print_term("    ", "Finished placing protein nr", protein_nr)
                
                xs, ys, zs = [[vals[ax] for key, vals in protein["beads"].items()] for ax in ["x", "y", "z"]]
                xmax, ymax, zmax, xmin, ymin, zmin = [func(ax) for func in [max, min] for ax in [xs, ys, zs]]

                beads_centered = copy.deepcopy(protein["beads"])

                #################
                ### CENTERING ###
                #################
                ### Centered on the mean coordinate of all beads (center of geometry)
                if protein["cen_method"][0] == "cog":
                    xcen, ycen, zcen = [np.mean(ax) for ax in [xs, ys, zs]]

                ### Centered on the mean of largest/smallest x/y/z coordinate
                if protein["cen_method"][0] == "axis":
                    xcen, ycen, zcen = np.mean([xmax, xmin]), np.mean([ymax, ymin]), np.mean([zmax, zmin])

                ### Centered on a single bead
                if protein["cen_method"][0] == "bead":
                    if len(protein["cen_method"]) == 2:
                        beads = [protein["cen_method"][1]]
                        
                    elif len(protein["cen_method"]) == 3:
                        if protein["cen_method"][1] < protein["cen_method"][2]:
                            beads = list(range(protein["cen_method"][1],
                                               protein["cen_method"][2]+1))
                        elif protein["cen_method"][1] > protein["cen_method"][2]:
                            beads = list(range(protein["cen_method"][2],
                                               protein["cen_method"][1],
                                               -1))
                        elif protein["cen_method"][1] == protein["cen_method"][2]:
                            beads = [protein["cen_method"][1]]
                            
                    xcen, ycen, zcen = [
                        vals[ax]
                        for key, vals in protein["beads"].items()
                        for ax in ["x", "y", "z"]
                        ### If atom_nr
                        if key[1] in beads
                    ]

                ### Centered on the mean position of all beads in a single residue
                if protein["cen_method"][0] == "res":
                    if len(protein["cen_method"]) == 2:
                        residues = [protein["cen_method"][1]]
                        
                    elif len(protein["cen_method"]) == 3:
                        if protein["cen_method"][1] < protein["cen_method"][2]:
                            residues = list(range(protein["cen_method"][1],
                                               protein["cen_method"][2]+1))
                        elif protein["cen_method"][1] > protein["cen_method"][2]:
                            residues = list(range(protein["cen_method"][2],
                                               protein["cen_method"][1],
                                               -1))
                        elif protein["cen_method"][1] == protein["cen_method"][2]:
                            residues = [protein["cen_method"][1]]

                    zipped = zip(*[
                        (vals["x"], vals["y"], vals["z"])
                        for key, vals in protein["beads"].items()
                        ### If res_nr
                        if key[2] in residues
                    ])
                    xcen, ycen, zcen = [np.mean(ax) for ax in zipped]

                ### Centered on the specific x/y/z coordinates
                if protein["cen_method"][0] == "point":
                    xcen, ycen, zcen = protein["cen_method"][1:]

                ### Doing the actual centering
                for bead in protein["beads"].keys():
                    for i, (ax, mean) in enumerate(zip(["x", "y", "z"], [xcen, ycen, zcen])):
                        beads_centered[bead][ax] -= mean

                self.print_term("    ", "Centering protein using", "'" + " ".join([str(i) for i in protein["cen_method"]])+"'", "at x/y/z:", round(xcen, 3), round(ycen, 3), round(zcen, 3), "(Input file coordinate system [Å])")

                #################
                ### ROTATIONS ###
                #################
                x_ang_deg, y_ang_deg, z_ang_deg = protein["rx"], protein["ry"], protein["rz"]
                if any([ang != 0 for ang in [x_ang_deg, y_ang_deg, z_ang_deg]]):
                    for bead in beads_centered.keys():
                        beads_centered[bead]["x"], beads_centered[bead]["y"], beads_centered[bead]["z"] = self.rotate_point(
                            beads_centered[bead]["x"],
                            beads_centered[bead]["y"],
                            beads_centered[bead]["z"],
                            x_ang_deg,
                            y_ang_deg,
                            z_ang_deg
                        )

                ####################
                ### TRANSLATIONS ###
                ####################
                tx, ty, tz = protein["tx"], protein["ty"], protein["tz"]
                if any([ax != 0 for ax in [tx, ty, tz]]):
                    for bead in beads_centered.keys():
                        beads_centered[bead]["x"] += tx
                        beads_centered[bead]["y"] += ty
                        beads_centered[bead]["z"] += tz

                #################################################
                ### CHECKS IF COORDINATES ARE OUTSIDE THE BOX ###
                #################################################
                if protein["pbc_check"]:
                    errors_count = 0
                    for bead in beads_centered.keys():
                        bead_coords = [beads_centered[bead]["x"], beads_centered[bead]["y"], beads_centered[bead]["z"]]
                        checked_beads, error = self.coord_checker(bead_coords, self.pbc_box, error_count = True)
                        if error > 0:
                            errors_count += 1
                            beads_centered[bead]["x"] = checked_beads[0]
                            beads_centered[bead]["y"] = checked_beads[1]
                            beads_centered[bead]["z"] = checked_beads[2]
                    if errors_count > 0:
                        self.print_term("    ", "    ", "WARNING:", str(errors_count), "Beads are outside pbc. Moved to other side. Expect potential problems from this.", warn = True)
                        self.print_term("    ", "    ", "Please move the protein such that it fits within the pbc.", warn = True)

                xcen_new, ycen_new, zcen_new = 0, 0, 0
                xcen_new, ycen_new, zcen_new = xcen_new + tx, ycen_new + ty, zcen_new + tz
                self.print_term("    ", "New protein center at x/y/z:", round(xcen_new, 3), round(ycen_new, 3), round(zcen_new, 3), "(Internal coordinate system [Å])")

                self.system_charge += protein["tot_charge"]
                
                self.PROTEINS[protein_nr]["beads_centered"] = beads_centered
                
            self.print_term("------------------------------ PLACEMENT COMPLETE", "\n")

    ###################################
    ### SHAPE FUNCTIONS - RECTANGLE ###
    ###################################
    def leaflet_rectangle(self):
        if len(self.LEAFLETS) != 0:
            self.print_term("------------------------------ CREATING LEAFLETS")
            for leaflet_nr, leaflet in self.LEAFLETS.items():
                self.print_term("Starting leaflet nr", leaflet_nr, "with shape:", "'" + leaflet["shape"] + "'")
                self.print_term("    ", "Base parameters:", "x=" + str(leaflet["x"] / 10) + "nm", "y=" + str(leaflet["y"] / 10) + "nm", "APL=" + str(leaflet["apl"] / 100) + "nm^2")

                ### Initial area per lipid calculations and max potential lipids in area
                apl_sqrt = math.sqrt(leaflet["apl"])
                max_lipids_possible_decimal = leaflet["x"] * leaflet["y"] / leaflet["apl"]
                leaflet_area = leaflet["x"] * leaflet["y"]

                ### Rounding max number of possible lipids according to requested rounding method
                if leaflet["lip_round_func"][1] == int: # int rounding
                    max_lipids_possible = int(max_lipids_possible_decimal + 0.5)
                    func_for_printer = str(round(max_lipids_possible_decimal, 3)) + "+0.5"
                else: # round(), math.floor() or math.ceil() rounding
                    max_lipids_possible = leaflet["lip_round_func"][1](max_lipids_possible_decimal)
                    func_for_printer = str(round(max_lipids_possible_decimal, 3))

                self.print_term("    ", "Maximum possible lipids:", round(max_lipids_possible_decimal, 3), "==>", max_lipids_possible)
                self.print_term("    ", "Rounding method:", "'" + str(leaflet["lip_round_func"][0]) + "(" + func_for_printer + ")'")    

                #############################
                ### INITIAL GRID CREATION ###
                #############################
                grid_points = []
                ### Ensures that longest dimension is adjusted (more wiggle room for change along that axis)
                if leaflet["x"] >= leaflet["y"]:
                    ax1 = "x"
                    ax2 = "y"
                elif leaflet["x"] < leaflet["y"]:
                    ax1 = "y"
                    ax2 = "x"

                breaker = False
                for d1 in np.arange(0, leaflet[ax1], apl_sqrt):
                    grid_points.append(0)
                    for d2 in np.arange(0, leaflet[ax2], apl_sqrt):
                        ### Checks if further points should be added
                        if sum(grid_points) < max_lipids_possible:
                            ### Checks if points are exceeding leaflet x/y values
                            if d1 + 0.5 * apl_sqrt <= leaflet[ax1] or leaflet[ax1] > leaflet[ax2]:
                                ### Below was once needed but now causes problems. Keeping in case it is needed later
                                ### if d2 + 0.5 * apl_sqrt <= leaflet[ax2] or leaflet[ax2] > leaflet[ax1]:
                                grid_points[-1] += 1
                        else:
                            breaker = True
                            break
                    if breaker == True:
                        break

                ###################################
                ### GRID OPTIMIZATION/SQUEEZING ###
                ###################################
                if grid_points[0] > grid_points[-1]:
                    last_grid = grid_points[-1]
                    grid_points = grid_points[:-1]
                    ### "Transposes"(ish) grid points
                    transposed_grid_points = [len(grid_points) for _ in range(grid_points[0])]
                    for i in range(last_grid):
                        step = i
                        if step >= len(transposed_grid_points):
                            step -= (len(transposed_grid_points) * (step // len(transposed_grid_points)))
                        transposed_grid_points[step] += 1

                    ### Ensuring that differently sized gridpoints/bins will be spread out evenly
                    if min(transposed_grid_points) != max(transposed_grid_points):
                        first_min = transposed_grid_points.index(min(transposed_grid_points))
                        transposed_grid_points = self.n_list_mixer(transposed_grid_points[:first_min], transposed_grid_points[first_min:])
                else:
                    transposed_grid_points = [len(grid_points) for _ in range(grid_points[0])]

                ######################
                ### POINT CREATION ###
                ######################
                grid = []
                dim2points = len(transposed_grid_points)
                dim2 = leaflet[ax2] / dim2points
                for d2_nr, dim1points in enumerate(transposed_grid_points):
                    dim1 = leaflet[ax1] / dim1points
                    row = d2_nr
                    for d1_nr in range(dim1points):
                        d1pos = dim1 * d1_nr + dim1 * 0.5
                        d2pos = dim2 * d2_nr + dim2 * 0.5
                        col = d1_nr
                        if ax1 == "x":
                            grid.append(((d1pos, d2pos), (dim1, dim2), (row, col)))
                        if ax1 == "y":
                            grid.append(((d2pos, d1pos), (dim2, dim1), (row, col)))

                ######################################
                ### CONVERT TO NESTED DICTIONARIES ###
                ######################################
                grid = {
                    (grid_point_nr, row, col): {
                        "internal_x": xpos,
                        "internal_y": ypos,
                        "x": xpos,
                        "y": ypos,
                        "xdim": xdim,
                        "ydim": ydim,
                        "Overlap": False,
                        "Inside_Poly": False,
                        "InvalidPoint": False,
                        "Removed": False,
                        ### For plotting:
                        "Inside_ConvexHull": False,
                    } for grid_point_nr, ((xpos, ypos), (xdim, ydim), (row, col)) in enumerate(grid)
                }

                self.print_term("    ", "Number of grid-points created:", len(grid))

                grid = self.fix_points(grid, leaflet)

                maxx, maxy, maxz, minx, miny, minz = [
                    func([
                        val + (leaflet[wr + "WR"] * sign)
                        for lipid in leaflet["lipids"].keys()
                        for val, bead in zip(leaflet["lipids"][lipid][ax], leaflet["lipids"][lipid]["beads"])
                        if bead != ""
                    ])
                    for func, sign in [(max, +1), (min, -1)] for ax, wr in [("x", "plane"), ("y", "plane"), ("z", "height")]
                ]
                
                lipid_dimensions = {
                    "maxx": maxx,
                    "maxy": maxy,
                    "maxz": maxz,
                    "minx": minx,
                    "miny": miny,
                    "minz": minz,
                    "zOverLipCen": abs(maxz),
                    "zUnderLipCen": abs(minz),
                    "lipid_height": abs(maxz) + abs(minz),
                    "lipid_radius": max(abs(maxx), abs(maxy), abs(minx), abs(miny)),
                }

                ### Update leaflet dictionary. Some are duplicated here, but potentially overwritten later
                self.LEAFLETS[leaflet_nr].update({
                    "grid": grid,
                    "grid_adjusted": grid,
                    
                    "lipid_dimensions": lipid_dimensions,
                    
                    "leaflet_area": leaflet_area,
                    "refined_area": leaflet_area,
                    
                    "max_lipids_possible": max_lipids_possible,
                    "new_max_lipids_possible": max_lipids_possible,
                    
                    "prot_area": 0,

                    "grid_overlapping": {},
                    "grid_removed": {},
                })
            self.print_term("------------------------------ LEAFLET CREATION COMPLETE", "\n")

    ####################################
    ### GRID-PROTEIN OVERLAP CHECKER ###
    ####################################
    def overlap_checker(self):
        if self.plots_requested:
            self.plot_data = {
                "Grid_ConvexHull_Overlap": [],
                "Grid_ProtBead_Overlap": [],
                "Prot_ConvexHull": [],
                "Prot_ConcaveHull_Polygon": [],
                "Prot_ConcaveHull_Polygon_Buffered": [],
                "Prot_ConcaveHulls_Points": [],
                "Delaunay_triangulation": [],
                "Prot_Flattened_points": [],
                "final_clusters": [],
            }
            self.print_term("Plot data will be obtained as requested", "\n")
        
        if len(self.PROTEINS) != 0 and len(self.LEAFLETS) != 0:
            self.print_term("------------------------------ CHECKING FOR PROTEIN-LEAFLET OVERLAP")
            for leaflet_i, (leaflet_nr, leaflet) in enumerate(self.LEAFLETS.items()):
                if leaflet_i != 0:
                    self.print_term()
                self.print_term("Starting leaflet nr", leaflet_nr)

                if self.plots_requested:
                    self.plot_data["Grid_ConvexHull_Overlap"].append([])
                    self.plot_data["Grid_ProtBead_Overlap"].append([])
                    self.plot_data["Prot_ConvexHull"].append([])
                    self.plot_data["Prot_ConcaveHull_Polygon"].append([])
                    self.plot_data["Prot_ConcaveHull_Polygon_Buffered"].append([])
                    self.plot_data["Prot_ConcaveHulls_Points"].append([])
                    self.plot_data["Delaunay_triangulation"].append([])
                    self.plot_data["Prot_Flattened_points"].append([])
                    self.plot_data["final_clusters"].append([])

                ######################################################
                ### Creating dictionary for removal of grid points ###
                ######################################################

                for protein_i, (protein_nr, protein) in enumerate(self.PROTEINS.items()):
                    if self.plots_requested:
                        self.plot_data["Grid_ConvexHull_Overlap"][leaflet_i].append([])
                        self.plot_data["Grid_ProtBead_Overlap"][leaflet_i].append([])
                        self.plot_data["Prot_ConvexHull"][leaflet_i].append([])
                        self.plot_data["Prot_ConcaveHull_Polygon"][leaflet_i].append([])
                        self.plot_data["Prot_ConcaveHull_Polygon_Buffered"][leaflet_i].append([])
                        self.plot_data["Prot_ConcaveHulls_Points"][leaflet_i].append([])
                        self.plot_data["Delaunay_triangulation"][leaflet_i].append([])
                        self.plot_data["Prot_Flattened_points"][leaflet_i].append([])
                        self.plot_data["final_clusters"][leaflet_i].append([])

                    self.print_term("    ", "Starting protein nr", protein_nr)

                    leaf_grid_points = [(val["x"], val["y"]) for grid_point_nr, val in leaflet["grid"].items()]

                    prod_beads_for_plot = [
                        (val["x"], val["y"], val["z"])
                        for bead_nr, val in protein["beads_centered"].items()
                        if (leaflet["center"][2] < val["z"] < + leaflet["center"][2] + leaflet["lipid_dimensions"]["lipid_height"]
                            and leaflet["HG_direction"] == "up")
                        or (leaflet["center"][2] > val["z"] > - leaflet["center"][2] - leaflet["lipid_dimensions"]["lipid_height"]
                            and leaflet["HG_direction"] == "down")
                    ]
                    prod_beads_in_memb = [(x, y) for x, y, z in prod_beads_for_plot]

                    if self.plots_requested:
                        self.plot_data["Prot_Flattened_points"][leaflet_i][protein_i].append(prod_beads_for_plot)
        #                 plot_data["Prot_ConvexHull"][leaflet_i].append(scipy.spatial.ConvexHull(prod_beads_in_memb))
                        CH = scipy.spatial.ConvexHull(prod_beads_in_memb)
                        CH_vertices = CH.vertices
                        self.plot_data["Prot_ConvexHull"][leaflet_i][protein_i].append([(x, y) for x, y in np.asarray(prod_beads_in_memb)[CH_vertices]])

                    if len(prod_beads_in_memb) == 0:
                        self.print_term("    ", "No overlap found")
                    else:
                        self.print_term("    ", "Leaflet-protein overlap found")
                        prot_bead_grid_point_in_radius = [
                            any([
                                (px - lx) ** 2 + (py - ly) ** 2 < leaflet["lipid_dimensions"]["lipid_radius"] ** 2
                                for px, py in prod_beads_in_memb
                            ]) for lx, ly in leaf_grid_points
                        ]
                        if self.plots_requested:
                            self.plot_data["Grid_ProtBead_Overlap"][leaflet_i][protein_i].append(prot_bead_grid_point_in_radius)
                            ### Delaunay triangulation
                            delaunay_triangulation = scipy.spatial.Delaunay(prod_beads_in_memb, furthest_site = False)
                            self.plot_data["Delaunay_triangulation"][leaflet_i][protein_i].append(delaunay_triangulation)

                        ######################################################################
                        ### REMOVING GRID POINTS BASED ON PROT_BEAD TO GRID_POINT DISTANCE ### ALWAYS DONE
                        ######################################################################
                        for i, (grid_point_key, TF) in enumerate(list(zip(leaflet["grid"].keys(), prot_bead_grid_point_in_radius))):
                            if TF:
                                self.LEAFLETS[leaflet_nr]["grid"][grid_point_key]["Overlap"] = True
                                self.LEAFLETS[leaflet_nr]["grid"][grid_point_key]["InvalidPoint"] = True
                            if self.plots_requested:
                                if self.plot_data["Delaunay_triangulation"][leaflet_i][protein_i][-1].find_simplex(leaf_grid_points[i], tol = 0) >= 0:
                                    self.LEAFLETS[leaflet_nr]["grid"][grid_point_key]["Inside_ConvexHull"] = True

                        ######################################
                        ### CLUSTERING PROTEIN BEAD POINTS ###
                        ######################################
                        ### Clustering done by diameter of lipids
                        prot_beads = prod_beads_in_memb[:]
                        final_clusters = []
                        while len(prot_beads) > 0:
                            changes = True
                            cluster = []
                            while changes:
                                start_len = len(cluster)
                                for i, point in reversed(list(enumerate(prot_beads))):
                                    if cluster == []:
                                        cluster.append(point)
                                        del prot_beads[i]
                                    else:
                                        dist_check = any([math.dist(point, cluster_point) <= leaflet["lipid_dimensions"]["lipid_radius"] * 2 for cluster_point in cluster])
                                        if dist_check:
                                            cluster.append(point)
                                            del prot_beads[i]
                                end_len = len(cluster)
                                if start_len == end_len:
                                    changes = False

                            final_clusters.append(cluster)
                        if self.plots_requested:
                            self.plot_data["final_clusters"][leaflet_i][protein_i].append(final_clusters)

                        ###############################
                        ### CONCAVE HULL ALPHASHAPE ###
                        ###############################
                        ConcaveHulls_Polygon = []
                        ConcaveHulls_Polygon_Buffered = []
                        ConcaveHulls_Points = []
                        alpha = 1 / (leaflet["lipid_dimensions"]["lipid_radius"] * protein["alpha_mult"])
                        for cluster in final_clusters:
                            if len(cluster) > 2:
                                ### alpha = radius^-1
                                ALPHASHAPE = alphashape.alphashape(points = cluster, alpha = alpha)
                                self.print_term(type(ALPHASHAPE), debug = True)

                                ### Alphashape output can be a bit unpredictable so need to check all possibilities
                                if ALPHASHAPE.geom_type == "Polygon":
                                    self.print_term("Polygon", debug = True)
                                    ConcaveHulls_Polygon.append([ALPHASHAPE])
                                elif ALPHASHAPE.geom_type == "MultiPolygon":
                                    self.print_term("MultiPolygon", debug = True)
                                    ConcaveHulls_Polygon.append(list(ALPHASHAPE.geoms))
                                elif ConcaveHulls_Polygon[-1].geom_type == "MultiPolygon":
                                    self.print_term("MultiPolygon", debug = True)
                                    ConcaveHulls_Polygon.append(list(ALPHASHAPE[-1].geoms))

                                ### Iterate over polygons
                                for poly in ConcaveHulls_Polygon[-1]:
                                    ConcaveHulls_Points.append(list(poly.exterior.coords))
                                    ### Buffering to add a bit of space around protein concave hull
                                    poly_buffered = poly.buffer(distance = protein["buffer"]) # [Å]
                                    ConcaveHulls_Polygon_Buffered.append([poly_buffered])
                                    self.LEAFLETS[leaflet_nr]["prot_area"] += poly_buffered.area

                                    for grid_point_key in leaflet["grid"].keys():
                                        point = (leaflet["grid"][grid_point_key]["x"], leaflet["grid"][grid_point_key]["y"])
                                        shapely_point = shapely.geometry.Point(point)
                                        ### Only marks grid-point if lipids "inside == False"
                                        if poly.contains(shapely_point) and not protein["lipids_inside"]:
                                            self.LEAFLETS[leaflet_nr]["grid"][grid_point_key]["Inside_Poly"] = True
                                            self.LEAFLETS[leaflet_nr]["grid"][grid_point_key]["InvalidPoint"] = True

                        if self.plots_requested:
                            self.plot_data["Prot_ConcaveHull_Polygon"][leaflet_i][protein_i].append(ConcaveHulls_Polygon)
                            self.plot_data["Prot_ConcaveHull_Polygon_Buffered"][leaflet_i][protein_i].append(ConcaveHulls_Polygon_Buffered)
                            self.plot_data["Prot_ConcaveHulls_Points"][leaflet_i][protein_i].append(ConcaveHulls_Points)
            self.print_term("------------------------------ OVERLAP CHECKING COMPLETE", "\n")
    
    #################################################
    ### LEAFLET READJUSTMENTS AFTER OVERLAP CHECK ###
    #################################################
    def leaflet_adjuster(self):
        if len(self.LEAFLETS) != 0 and len(self.PROTEINS) != 0:
            self.print_term("------------------------------ ADJUSTING LEAFLETS DUE TO OVERLAP")
            for leaflet_i, (leaflet_nr, leaflet) in enumerate(self.LEAFLETS.items()):
                if leaflet_i != 0:
                    self.print_term()
                self.print_term("Starting leaflet nr", leaflet_nr)

                if leaflet["prot_area"] == 0:
                    continue
                else:
                    self.print_term("    ", "Leaflet area:", round(leaflet["leaflet_area"], 3))
                    self.print_term("    ", "Protein area:", round(leaflet["prot_area"], 3))
                    self.print_term("    ", "Area remaning:", round(leaflet["leaflet_area"] - leaflet["prot_area"], 3))

                    old_max_lipids_possible_decimal = leaflet["leaflet_area"] / leaflet["apl"]
                    new_max_lipids_possible_decimal = (leaflet["leaflet_area"] - leaflet["prot_area"]) / leaflet["apl"]

                    grid_points_NoOverlap = 0
                    for grid_point_nr in leaflet["grid"].keys():
                        if not leaflet["grid"][grid_point_nr]["InvalidPoint"]:
                            grid_points_NoOverlap += 1

                    ### Rounding max number of possible lipids according to requested rounding method
                    if leaflet["lip_round_func"][1] == int: # int rounding
                        old_max_lipids_possible = int(old_max_lipids_possible_decimal + 0.5)
                        new_max_lipids_possible = int(new_max_lipids_possible_decimal + 0.5)
                        func_for_printer = str(round(new_max_lipids_possible_decimal, 3)) + "+0.5"
                    else: # round(), math.floor() or math.ceil() rounding
                        old_max_lipids_possible = self.LEAFLETS[leaflet_nr]["lip_round_func"][1](old_max_lipids_possible_decimal)
                        new_max_lipids_possible = self.LEAFLETS[leaflet_nr]["lip_round_func"][1](new_max_lipids_possible_decimal)
                        func_for_printer = str(round(new_max_lipids_possible_decimal, 3))

                    self.LEAFLETS[leaflet_nr]["old_max_lipids_possible"] = old_max_lipids_possible
                    self.LEAFLETS[leaflet_nr]["new_max_lipids_possible"] = new_max_lipids_possible
                    self.LEAFLETS[leaflet_nr]["grid_points_NoOverlap"]   = grid_points_NoOverlap

                    self.print_term("    ", "Old maximum possible lipids:", round(old_max_lipids_possible_decimal, 3), "==>", old_max_lipids_possible)
                    self.print_term("    ", "New maximum possible lipids:", round(new_max_lipids_possible_decimal, 3), "==>", new_max_lipids_possible)
                    self.print_term("    ", "Rounding method:", "'" + str(leaflet["lip_round_func"][0]) + "(" + func_for_printer + ")'")
                    self.print_term("    ", "Grid points available after protein overlap:", grid_points_NoOverlap)

                    grid_nr_diff = new_max_lipids_possible - grid_points_NoOverlap
                    ###############################
                    ###### ADDING NEW POINTS ######
                    ###############################
                    if grid_nr_diff > 0:
                        ### Ensures that longest dimension is adjusted (more wiggle room for change along that axis)
                        if leaflet["x"] >= leaflet["y"]:
                            ax1 = "x"
                            ax2 = "y"
                        elif leaflet["x"] <leaflet["y"]:
                            ax1 = "y"
                            ax2 = "x"
                        
                        self.print_term("    ", "Squeezing in", grid_nr_diff, "extra grid points as adjustment")
                        rows_without_prot = []
                        current_row = 0
                        prot_in_row = False
                        point_counter = 0
                        rows = sorted(list(set([row for grid_nr, row, col in leaflet["grid"].keys()])))

                        ########################################################
                        ### FIGURING OUT WHICH ROWS/COLUMNS TO ADD POINTS TO ###
                        ########################################################
                        ### Ignores rows containing protein overlaps
                        for row in rows:
                            TF = [
                                leaflet["grid"][(grid_nr, row2, col)]["InvalidPoint"]
                                for grid_nr, row2, col in leaflet["grid"].keys()
                                if row == row2
                            ]
                            if not any(TF):
                                points_in_row = max(sorted(list(set([
                                    ### Use internal ax2 position due to centering done later
                                    (col + 1, vals["internal_" + ax2], vals[ax2 + "dim"])
                                    for (grid_nr, row2, col), vals
                                    in leaflet["grid"].items()
                                    if row == row2
                                ]))))
                                rows_without_prot.append((row,) + points_in_row)

                        ### Adding the new points ### Currently not spread evenly across leaflet rows
                        rows_sorted = sorted(rows_without_prot, key = lambda p: p[1])
                        for extra_points in range(grid_nr_diff):
                            extra_points -= extra_points // len(rows_sorted) * (len(rows_sorted))
                            r, p, ax2pos, dim2 = rows_sorted[extra_points]
                            rows_sorted[extra_points] = (r, p + 1, ax2pos, dim2)

                        ######################
                        ### POINT CREATION ###
                        ######################
                        grid = []
                        for d2_nr, dim1points, d2pos, dim2 in rows_sorted:
                            dim1 = leaflet[ax1] / dim1points
                            for d1_nr in range(dim1points):
                                d1pos = dim1 * d1_nr + dim1 * 0.5
                                row = d2_nr
                                col = d1_nr
                                if ax1 == "x":
                                    grid.append(((d1pos, d2pos), (dim1, dim2), (row, col)))
                                if ax1 == "y":
                                    grid.append(((d2pos, d1pos), (dim2, dim1), (row, col)))

                        ######################################
                        ### CONVERT TO NESTED DICTIONARIES ###
                        ######################################
                        new_grid_dict = {}
                        for new_point_nr, ((xpos, ypos), (xdim, ydim), (row, col)) in enumerate(grid):
                            new_grid_dict[new_point_nr, row, col] = {
                                    "internal_x": xpos,
                                    "internal_y": ypos,
                                    "x": xpos,
                                    "y": ypos,
                                    "xdim": xdim,
                                    "ydim": ydim,
                                    "Overlap": False,
                                    "Inside_Poly": False,
                                    "InvalidPoint": False,
                                    "Removed": False,
                                    ### For plotting:
                                    "Inside_ConvexHull": False,
                                }

                        ### Fixes points according to leaflet parameters
                        new_grid_dict = self.fix_points(new_grid_dict, leaflet)

                        ### Creates final adjusted dictionary containing only valid grid points and dictionary with the non-valid grid points
                        max_point_nr = max(point_nr for point_nr, row, point in new_grid_dict.keys())
                        old_remaining_grid = {
                            (new_point_nr, row, point):val
                            for new_point_nr, ((point_nr, row, point), val) in enumerate(leaflet["grid"].items(), max_point_nr)
                            if row not in list(zip(*rows_sorted))[0] and val["InvalidPoint"] == False
                        }
                        new_grid_dict.update(old_remaining_grid)

                        max_point_nr = max(point_nr for point_nr, row, point in new_grid_dict.keys())
                        old_overlapping_grid = {
                            (new_point_nr, row, point):val
                            for new_point_nr, ((point_nr, row, point), val) in enumerate(leaflet["grid"].items(), max_point_nr)
                            if row not in list(zip(*rows_sorted))[0] and val["InvalidPoint"] == True
                        }

                        self.LEAFLETS[leaflet_nr]["grid_adjusted"] = new_grid_dict
                        self.LEAFLETS[leaflet_nr]["grid_overlapping"] = old_overlapping_grid
                        self.LEAFLETS[leaflet_nr]["grid_removed"] = {}
                    
                    #############################
                    ###### REMOVING POINTS ######
                    #############################
                    elif grid_nr_diff < 0:
                        self.print_term("    ", "Removing", abs(grid_nr_diff), "random grid points as adjustment")

                        ### Creates adjusted dictionary containing only valid grid points and dictionary with the non-valid grid points due to overlapping
                        new_grid_dict = {
                            (new_point_nr, row, point):val
                            for new_point_nr, ((point_nr, row, point), val) in enumerate(leaflet["grid"].items())
                            if val["InvalidPoint"] == False
                        }
                        
                        max_point_nr = max(point_nr for point_nr, row, point in new_grid_dict.keys())
                        old_overlapping_grid = {
                            (new_point_nr, row, point):val
                            for new_point_nr, ((point_nr, row, point), val) in enumerate(leaflet["grid"].items(), max_point_nr)
                            if val["InvalidPoint"] == True
                        }

                        ### Removing grid points
                        removed_grid_dict = {}
                        for _ in range(abs(grid_nr_diff)):
                            grid_key = random.choice(list(new_grid_dict.keys()))
                            removed_grid_dict[grid_key] = new_grid_dict[grid_key]
                            del new_grid_dict[grid_key]

                        self.LEAFLETS[leaflet_nr]["grid_adjusted"] = new_grid_dict
                        self.LEAFLETS[leaflet_nr]["grid_overlapping"] = old_overlapping_grid
                        self.LEAFLETS[leaflet_nr]["grid_removed"] = removed_grid_dict
            self.print_term("------------------------------ LEAFLET ADJUSTMENTS COMPLETE", "\n")
    
    ######################
    ### LIPID INSERTER ###
    ######################
    def lipid_inserter(self):
        self.SYS_lipids_dict = {}
        if len(self.LEAFLETS) != 0:
            self.print_term("------------------------------ CALCULATING LIPID RATIOS")
            for leaflet_i, (leaflet_nr, leaflet) in enumerate(self.LEAFLETS.items()):
                if leaflet_i != 0:
                    self.print_term()
                self.print_term("Starting leaflet nr", leaflet_nr)

                ### Initial rounded estimations
                lipids = [(lipid_name, lipid_vals["ratio"]) for lipid_name, lipid_vals in leaflet["lipids"].items()]
                lipid_names = [name for name, ratio in lipids]
                lipid_ratios_decimal = [ratio for name, ratio in lipids]
                lipids_tot = sum(lipid_ratios_decimal)
                lipid_ratios = [round(int(i / lipids_tot * leaflet["new_max_lipids_possible"]), 3) for i in lipid_ratios_decimal]

                original_ratios_decimal = [round(i / sum(lipid_ratios) * 100, 3) for i in lipid_ratios]
                original_ratios = lipid_ratios[:]

                expected_ratios_decimal = [round(ratio / sum(lipid_ratios_decimal) * 100, 3) for ratio in lipid_ratios_decimal]

                self.print_term("    ", "Lipid optimization 'lipid_optim' setting:", leaflet["lipid_optim"])

                ### Just take integer lipid values and remove excess grid points
                if leaflet["lipid_optim"] == "abs_val":
                    lipid_ratios = [ratio for name, ratio in lipids]
                    ### Error out if values are not integers
                    assert all([type(val) == int for val in lipid_ratios]), "Not all supplied lipid values are integers. You can't have half a lipid, no matter how much you want to."
                    ### Error out if more lipids than available grid points
                    assert sum(lipid_ratios) <= leaflet["new_max_lipids_possible"], "You have specified too many lipids. Sorry they are too fat and won't fit in, try reducing number of lipids or the APL."
                    grid_points_to_remove = leaflet["new_max_lipids_possible"] - sum(lipid_ratios)

                elif leaflet["lipid_optim"] in ["force_fill", "fill", "avg_optimal"]:
                    ###############################
                    ### OPTIMIZING LIPID RATIOS ###
                    ###############################
                    grid_points_to_remove = 0

                    iters = [lipid_ratios[:]]

                    while sum(lipid_ratios) < leaflet["new_max_lipids_possible"]:
                        diffs = [(round(lipids[i][1], 3) - round(lip_val / (sum(lipid_ratios) / lipids_tot), 3), i) for i, lip_val in enumerate(lipid_ratios)]
                        biggest_diff, biggest_diff_i = max(diffs, key=lambda d: d[0])

                        ### In case ideal ratio is found, only continue if force_fill
                        if biggest_diff == float(0) and leaflet["lipid_optim"] == "force_fill":
                            lipid_ratios[lipid_ratios.index(min(lipid_ratios))] += 1
                            iters.append(lipid_ratios)
                        elif biggest_diff == float(0):
                            break
                        else:
                            lipid_ratios[biggest_diff_i] += 1
                            iters.append(lipid_ratios[:])

                    ### Filling in lipids
                    if leaflet["lipid_optim"] in ["force_fill", "fill"]:
                        if leaflet["lipid_optim"] == "fill":
                            self.print_term("    ", "    ", "Filling leaflet until perfect ratio between lipids is achieved or leaflet is full")
                        elif leaflet["lipid_optim"] == "force_fill":
                            self.print_term("    ", "    ", "Forcefully filling leaflet until all grid points have a lipid")
                        lipid_ratios = iters[-1]
                        grid_points_to_remove += leaflet["new_max_lipids_possible"] - sum(lipid_ratios)

                    ### Optimizing lipids
                    elif leaflet["lipid_optim"] == "avg_optimal": # lipid_ratios_decimal
                        self.print_term("    ", "    ", "Optimizing based on the average deviation from expected ratios")
                        
                        iters_decimal = [[round(iteration[i] / sum(iteration) * 100, 3) for i in range(len(iteration))] for iteration in iters]
                        iters_abs_diff = [
                            [abs(iteration[i] - expected_ratios_decimal[i]) for i in range(len(iteration))]
                            for iteration in iters_decimal
                        ]

                        iters_avg = [(np.mean(iters), i) for i, iters in enumerate(iters_abs_diff)]
                        iters_avg_optimal = min(iters_avg, key=lambda i: i[0])
                        lipid_ratios = iters[iters_avg_optimal[1]]
                        grid_points_to_remove += leaflet["new_max_lipids_possible"] - sum(lipid_ratios)

                ### No optimization or filling
                elif leaflet["lipid_optim"] == "no":
                    grid_points_to_remove = leaflet["new_max_lipids_possible"] - sum(lipid_ratios)

                if grid_points_to_remove:
                    self.print_term("    ", "    ", "Removing", "'" + str(grid_points_to_remove) + "'","random grid points to match number of lipids")
                    for _ in range(grid_points_to_remove):
                        grid_key = random.choice(list(leaflet["grid_adjusted"].keys()))
                        self.LEAFLETS[leaflet_nr]["grid_removed"][grid_key] = leaflet["grid_adjusted"][grid_key]
                        del self.LEAFLETS[leaflet_nr]["grid_adjusted"][grid_key]
                
                self.LEAFLETS[leaflet_nr]["lipid_ratios"] = lipid_ratios
                self.LEAFLETS[leaflet_nr]["lipid_names"] = lipid_names

                ###############################################
                ### PRINTING LEAFLET SPECIFIC LIPID DETAILS ###
                ###############################################
                if self.extra_info:
                    ### Printing mean, min and max deviation from wanted ratios
                    ratios_final_decimal = [round(lipid_ratios[i] / sum(lipid_ratios) * 100, 3) for i in range(len(lipid_ratios))]
                    ratios_final_abs_diff = [abs(ratios_final_decimal[i] - expected_ratios_decimal[i]) for i in range(len(ratios_final_decimal))]
                    ratios_final_avg = round(np.mean(ratios_final_abs_diff), 5) # Don't multiply by 100 as they are already percentages
                    ratios_final_min = round(min(ratios_final_abs_diff), 5) # Don't multiply by 100 as they are already percentages
                    ratios_final_max = round(max(ratios_final_abs_diff), 5) # Don't multiply by 100 as they are already percentages

                    self.print_term("    ", "Final deviations from expected ratios:", "Difference in %-values")
                    self.print_term("    ", "    ", "Maximum:", ratios_final_max)
                    self.print_term("    ", "    ", "Average:", ratios_final_avg)
                    self.print_term("    ", "    ", "Minimum:", ratios_final_min)

                    final_ratios_decimal = [round(ratio / sum(lipid_ratios) * 100, 3) for ratio in lipid_ratios]
                    final_lipid_data = [lipid_names, lipid_ratios, final_ratios_decimal, expected_ratios_decimal, original_ratios_decimal, original_ratios]

                    headers = ["Lipid name", "Final lipids", "Final %", "Expected %", "Starting %", "Starting lipids"]
                    lipids_for_printer = [[head] + data for head, data in zip(headers, final_lipid_data)]
                    max_lengths = [len(str(max(data, key = lambda d: len(str(d))))) for data in lipids_for_printer]

                    self.print_term("    ", "Leaflet specific lipid data")
                    for i, (L0, L1, L2, L3, L4, L5) in enumerate(list(zip(*lipids_for_printer))):
                        self.print_term(
                            "    ", "    ",
                            '{0: <{L}}'.format(L0, L = max_lengths[0]), ":",
                            '{0: <{L}}'.format(L1, L = max_lengths[1]), ":",
                            '{0: <{L}}'.format(L2, L = max_lengths[2]), ":",
                            '{0: <{L}}'.format(L3, L = max_lengths[3]), ":",
                            '{0: <{L}}'.format(L4, L = max_lengths[4]), ":",
                            '{0: <{L}}'.format(L5, L = max_lengths[5]),
                        )
                        if i != 0:
                            if L0 not in self.SYS_lipids_dict.keys():
                                self.SYS_lipids_dict[L0] = L1
                            else:
                                self.SYS_lipids_dict[L0] += L1
            for leaflet_i, (leaflet_nr, leaflet) in enumerate(self.LEAFLETS.items()):
                #########################################
                ### INSERTING LIPIDS INTO GRID POINTS ###
                #########################################
                randomized_list = random.sample([i for i in range(len(leaflet["grid_adjusted"].keys()))], k = len(leaflet["grid_adjusted"]))
                lipids_x_times = [l_name for l_name, l_nr in zip(leaflet["lipid_names"], leaflet["lipid_ratios"]) for _ in range(l_nr)]
                randomized_lipids = sorted([(rand, l_name) for rand, l_name in zip(randomized_list, lipids_x_times)], key = lambda L: L[0])
                random_lipids_in_grids = [(rand, l_name, grid_key) for (rand, l_name), grid_key in zip(randomized_lipids, leaflet["grid_adjusted"].keys())]

                if leaflet["HG_direction"] == "up":
                    sign = +1
                if leaflet["HG_direction"] == "down":
                    sign = -1
                
                for rand, l_name, grid_key in random_lipids_in_grids:
                    grid_point_x = leaflet["grid_adjusted"][grid_key]["x"]
                    grid_point_y = leaflet["grid_adjusted"][grid_key]["y"]
                    grid_point_z = leaflet["grid_adjusted"][grid_key]["z"] + (sign * leaflet["lipid_dimensions"]["zUnderLipCen"])
                    xs_rotated   = leaflet["lipids"][l_name]["x"]
                    new_x, new_y, new_z = [], [], []
                    random_rotaion_angle = random.uniform(0, 360)
                    for x, y, z in zip(leaflet["lipids"][l_name]["x"], leaflet["lipids"][l_name]["y"], leaflet["lipids"][l_name]["z"]):
                        nx, ny, nz = self.rotate_point(x, y, z, 0, 0, random_rotaion_angle)
                        new_x.append(nx)
                        new_y.append(ny)
                        new_z.append(nz)

                    lipid_dict = {
                        "name": l_name,
                        "beads": leaflet["lipids"][l_name]["beads"],
                        "x": [x + grid_point_x + random.uniform(-leaflet["kickx"], leaflet["kickx"]) for x in new_x],
                        "y": [y + grid_point_y + random.uniform(-leaflet["kicky"], leaflet["kicky"]) for y in new_y],
                        "z": [sign * z + grid_point_z + random.uniform(-leaflet["kickz"], leaflet["kickz"]) for z in new_z],
                    }
                    self.LEAFLETS[leaflet_nr]["grid_adjusted"][grid_key]["lipid"] = lipid_dict.copy()
                    self.LEAFLETS[leaflet_nr]["leaf_lipid_count"] = [(n, str(c)) for n, c in zip(leaflet["lipid_names"], leaflet["lipid_ratios"])]
                    self.system_charge += self.LEAFLETS[leaflet_nr]["lipids"][l_name]["charge_sum"]
            #########################################
            ### PRINTING SYSTEMWIDE LIPID DETAILS ###
            #########################################
            if self.extra_info:
                SYS_headers = ["Lipid name", "Total lipids", "Total %"]
                SYS_lipid_names = [name for name in self.SYS_lipids_dict.keys()]
                SYS_lipid_vals = [val for val in self.SYS_lipids_dict.values()]
                SYS_tot_lipids = sum(SYS_lipid_vals)
                SYS_lipid_percentages = [round(val / SYS_tot_lipids * 100, 3) for val in SYS_lipid_vals]

                SYS_lipids_for_printer = list(zip(SYS_lipid_names, SYS_lipid_vals, SYS_lipid_percentages))
                SYS_lipids_for_printer = [tuple(SYS_headers)] + SYS_lipids_for_printer
                SYS_max_lengths = [len(str(max(data, key = lambda d: len(str(d))))) for data in zip(*SYS_lipids_for_printer)]

                self.print_term("\nLipid data for whole system")
                for L0, L1, L2 in SYS_lipids_for_printer:
                    self.print_term(
                        "    ",
                        "    ",
                        '{0: <{L}}'.format(L0, L = SYS_max_lengths[0]), ":",
                        '{0: <{L}}'.format(L1, L = SYS_max_lengths[1]), ":",
                        '{0: <{L}}'.format(L2, L = SYS_max_lengths[2]),
                    )
            self.print_term("------------------------------ LIPID RATIO CALCULATIONS COMPLETE", "\n")
    
    ################
    ### SOLVATER ###
    ################
    def solvater(self):
        solv_beads_for_cell_checker = []
        if len(self.SOLVATIONS_cmds) != 0:
            self.print_term("------------------------------ SOLVATING SYSTEM")
            for solvation_i, (solvation_nr, solvation) in enumerate(self.SOLVATIONS.items()):
                if solvation_i != 0:
                    self.print_term("")
                self.print_term("Starting solvation nr", solvation_nr)

                #########################################
                ### CALCULATION FREE VOLUME OF SYSTEM ###
                #########################################
                self.print_term("    ", "Calculating box volume: (all values in [nm^3])")
                self.print_term("    ", "    ", "Bead radius used for volume calculations 'bead_radius':", solvation["bead_radius"], "[Å]")

                bead_radius = solvation["bead_radius"]
                ### Leaflet volume based on beads exclusively
                non_solv_beads = []
                prot_beads_for_cell_checker = []
                lipid_beads_for_cell_checker = []
                leafs_volume = 0
                if len(self.LEAFLETS) != None:
                    '''
                    Estimates a volume for all lipids, and finds the lipid bead positions
                    '''
                    for leaflet_i, (leaflet_nr, leaflet) in enumerate(self.LEAFLETS.items()):
                        for grid_point in leaflet["grid_adjusted"].values():
                            xs, ys, zs = grid_point["lipid"]["x"], grid_point["lipid"]["y"], grid_point["lipid"]["z"]
                            lipid_beads = list(zip(xs, ys, zs))
                            leafs_volume += (4/3 * math.pi * (bead_radius ** 3)) * len(lipid_beads) * 10**-27
                            non_solv_beads.extend(lipid_beads)
        #                     if leaflet["HG_direction"] == "up":
        #                         gridz = abs(max(zs)) + abs(min(zs))
        #                     if leaflet["HG_direction"] == "down":
        #                         gridz = abs(max(zs)) + abs(min(zs))
        #                     lipid_grid_points_for_alg.append((grid_point["x"], grid_point["y"], grid_point["z"], gridz, LEAFS[leaf_nr]["HG_direction"]))
                            gridz = abs(max(zs)) + abs(min(zs))
                            lipid_beads_for_cell_checker.extend(lipid_beads)

                ### Protein volume based on beads exclusively
                prots_volume = 0
                if len(self.PROTEINS) != None:
                    '''
                    Estimates a volume for all proteins, and finds the protein bead positions
                    '''
                    for protein_i, (protein_nr, protein) in enumerate(self.PROTEINS.items()):
                        prot_points = [
                            (protein["beads_centered"][key]["x"], protein["beads_centered"][key]["y"], protein["beads_centered"][key]["z"])
                            for key in protein["beads_centered"].keys()
                        ]
                        prots_volume += (4/3 * math.pi * (bead_radius ** 3)) * len(prot_points) * 10**-27
                        prot_beads_for_cell_checker.extend(prot_points)

                ### Solvent volume from previous solvation commands
                solvs_volume = (4/3 * math.pi * (bead_radius ** 3)) * len(solv_beads_for_cell_checker) * 10**-27

                N_A = 6.02214076 * 10**23
                box_volume = (self.pbc_box[0] * self.pbc_box[1] * self.pbc_box[2]) * 10**-27
                non_free_volume = leafs_volume + prots_volume + prots_volume
                box_free_volume = box_volume - non_free_volume
                self.print_term("    ", "    ", "Box volume:            ", round(box_volume * 10**24, 3))
                self.print_term("    ", "    ", "Lipid volume:          ", round(leafs_volume * 10**24, 3))
                self.print_term("    ", "    ", "Protein volume:        ", round(prots_volume * 10**24, 3))
                self.print_term("    ", "    ", "(Prior) Solvent volume:", round(solvs_volume * 10**24, 3))
                self.print_term("    ", "    ", "Excluded volume:       ", round(non_free_volume * 10**24, 3))
                self.print_term("    ", "    ", "Free volume:           ", round(box_free_volume * 10**24, 3))

                ###########################
                ### CALCULATING SOLVENT ###
                ###########################
                solv_volume = 0
                sol_charges = 0
                sol_ratios = []

                ### List used to find the maximum solvent/ion size
                solvent_molecules = []

                ### Solvent
                if solvation["solvent"] != []:
                    for solv_key, solv_vals in solvation["solvent"].items():
                        '''
                        Calculates the number of solvent particles based on molarity, density and solvent count (CG-to-AA mapping).
                        Calculates the solvent volume to precisely estimate the amount of volume that the solvent should occupy.
                        Solvent volume can be used for ion-centration calculations (Used by default)
                        '''
                        if solvation["count"]: ### Treats molarity as absolute number of molecules
                            count = solv_vals["molarity"]
                        else: ### Treats molarity as molarity
                            ### Using the free volume for concentration calculations (default)
                            if solvation["solvfreevol"] == True:
                                count = int(N_A * box_free_volume * solv_vals["molarity"] / solv_vals["solvcount"])
                            ### Using box volume for concentration calculations
                            elif solvation["solvfreevol"] == False:
                                count = int(N_A * box_volume * solv_vals["molarity"] / solv_vals["solvcount"])

                        ### Adds charge to solvent charge if it is specified in defines
                        if "charge" in solv_vals.keys():
                            sol_charges += count * solv_vals["charge"]

                        ### Checks if molar mass and density is specified for the solvent in defines
                        ### Only prints warning if solvent volume is to be used for ions
                        if "molar_mass" not in solv_vals.keys() and solvation["ionsvol"] == "solv":
                            self.print_term("WARNING: Chosen solvent [" + solv_key + "] is missing a 'molar_mass' entry in the solvent defines. Please add it if you wish to use solvent volume for ion concentration.", warn = True)
                        if "density" not in solv_vals.keys() and solvation["ionsvol"] == "solv":
                            self.print_term("WARNING: Chosen solvent [" + solv_key + "] is missing a 'density' entry in the solvent defines. Please add it if you wish to use solvent volume for ion concentration.", warn = True)

                        solvation["solvent"][solv_key]["count"] = count
                        ### Adds volume for solvent type to solvent volume if molar mass and density are both specified
                        if all([i in solv_vals.keys() for i in ["molar_mass", "density"]]):
                            solv_volume += (count * solv_vals["solvcount"] * solv_vals["molar_mass"]) / (N_A * solv_vals["density"]) * 10**-3

                        solvent_molecules.append(list(zip(solv_vals["x"], solv_vals["y"], solv_vals["z"])))

                self.print_term("    ", "    ", "Solvent volume:        ", round(solv_volume * 10**24, 3))
                self.print_term("    ", "    ",)

                ########################
                ### CALCULATING IONS ###
                ########################
                pos_charges = 0
                pos_ratios = []

                if solvation["pos_ions"] != []:
                    ### ### Calculating number of particles for each ion and their contribution to the total charge
                    ### Positive ions
                    for pos_key, pos_vals in solvation["pos_ions"].items():
                        if solvation["count"]: ### Treats molarity as absolute number of molecules
                            count = pos_vals["molarity"]
                        else: ### Treats molarity as molarity
                            ### Using the solvent volume for concentration calculations (default)
                            if solvation["ionsvol"] == "solv":
                                count = int(N_A * solv_volume * pos_vals["molarity"] / pos_vals["solvcount"])
                            ### Using the free volume for concentration calculations
                            elif solvation["ionsvol"] == "free":
                                count = int(N_A * box_free_volume * pos_vals["molarity"] / pos_vals["solvcount"])
                            ### Using the box volume for concentration calculations
                            elif solvation["ionsvol"] == "box":
                                count = int(N_A * box_volume * pos_vals["molarity"] / pos_vals["solvcount"])

                        solvation["pos_ions"][pos_key]["count"] = count
                        pos_charges += count * pos_vals["charge"]
                        pos_ratios.append((pos_key, pos_vals["charge"], count)) # Used for ion optimizer/neutralization

                        solvent_molecules.append(list(zip(pos_vals["x"], pos_vals["y"], pos_vals["z"])))

                neg_charges = 0
                neg_ratios = []

                if solvation["neg_ions"] != []:
                    ### Negative ions
                    for neg_key, neg_vals in solvation["neg_ions"].items():
                        if solvation["count"]: ### Treats molarity as absolute number of molecules
                            count = neg_vals["molarity"]
                        else: ### Treats molarity as molarity
                            ### Using the solvent volume for concentration calculations (default)
                            if solvation["ionsvol"] == "solv":
                                count = int(N_A * solv_volume * neg_vals["molarity"] / neg_vals["solvcount"])
                            ### Using the free volume for concentration calculations
                            elif solvation["ionsvol"] == "free":
                                count = int(N_A * box_free_volume * neg_vals["molarity"] / neg_vals["solvcount"])
                            ### Using the box volume for concentration calculations
                            elif solvation["ionsvol"] == "box":
                                count = int(N_A * box_volume * neg_vals["molarity"] / neg_vals["solvcount"])

                        solvation["neg_ions"][neg_key]["count"] = count
                        neg_charges += solvation["neg_ions"][neg_key]["count"] * neg_vals["charge"]
                        neg_ratios.append((neg_key, neg_vals["charge"], count)) # Used for ion optimizer/neutralization

                        solvent_molecules.append(list(zip(neg_vals["x"], neg_vals["y"], neg_vals["z"])))

                ### Summing positive and negative charges
                tot_charges = pos_charges + neg_charges
                if solvation["sys_charge"]:
                    ### Taking into account system charge if requested (default)
                    tot_charges += self.system_charge

                ### If both positive and negative ions
                def ions_optimizer(ions, charge_change):
                    ions_abscharge = [(name, abs(charge), ratio) for name, charge, ratio in ions]
                    charged_ratios = [ratio * charge for name, charge, ratio in ions_abscharge]
                    weighted_ratios = [w_ratio / sum(charged_ratios) for w_ratio in charged_ratios]
                    int_ratios = [int(w_ratio * charge_change / ions_abscharge[i][1]) for i, w_ratio in enumerate(weighted_ratios)]
                    charge_change -= sum([ratio * ions_abscharge[i][1] for i, ratio in enumerate(int_ratios)])
                    change = True
                    while charge_change != 0 and change: 
                        change = False
                        for i, ratio in enumerate(int_ratios):
                            if charge_change - ions_abscharge[i][1] >= 0:
                                int_ratios[i] += 1
                                charge_change -= ions_abscharge[i][1]
                                change = True
                    return int_ratios

                ### ### "Neutralizes" system towards the requested charge-value
                if str(solvation["charge"]) != "False" and not solvation["count"]: # Strings required, else if charge is 0 then, 0 == False would be true
                    ### Adds extra negative ions
                    if tot_charges > solvation["charge"]:
                        if solvation["neg_ions"] != []:
                            neg_charges = 0 
                            charge_diff = abs(solvation["charge"] - tot_charges)
                            extra_ions = ions_optimizer(neg_ratios, charge_diff)

                            for (name, vals), count in zip(solvation["neg_ions"].items(), extra_ions):
                                solvation["neg_ions"][name]["count"] += count
                                neg_charges += vals["charge"] * vals["count"]

                        else:
                            self.print_term("WARNING: I cannot neutralize with NEGATIVE ions when none are given")

                    ### Adds extra positive ions
                    elif tot_charges < solvation["charge"]:
                        if solvation["pos_ions"] != []:
                            pos_charges = 0
                            charge_diff = abs(solvation["charge"] - tot_charges)
                            extra_ions = ions_optimizer(pos_ratios, charge_diff)

                            for (name, vals), count in zip(solvation["pos_ions"].items(), extra_ions):
                                solvation["pos_ions"][name]["count"] += count
                                pos_charges += vals["charge"] * vals["count"]

                        else:
                            self.print_term("WARNING: I cannot neutralize with POSITIVE ions when none are given")

                ### Adding charges from ions and solvent to system charge
                self.system_charge += sol_charges + pos_charges + neg_charges

                ### Randomly shuffle solvent particles for random distribution in box
                collected_solvent = []
                solvent_sizes = []
                for stype in ["solvent", "pos_ions", "neg_ions"]:
                    for key, vals in solvation[stype].items():
                        sname = key
                        scount = vals["count"]
                        maxx, maxy, maxz, minx, miny, minz = [
                            func([
                                val + (solvation["WR"] * sign)
                                for val in solvation[stype][key][ax]
                            ])
                            for func, sign in [(max, +1), (min, -1)] for ax in ["x", "y", "z"]
                        ]
                        ssize = max([maxi - mini for maxi, mini in [(maxx, minx), (maxy, miny), (maxz, minz)]])
                        solvent_sizes.append(ssize)
                        collected_solvent.extend([[sname, stype, ssize]] * scount)
                
                #########################################################################
                ### DEFINING FUNCTIONS USED FOR SOLVENT CELL DIMENSIONALITY REDUCTION ###
                #########################################################################
                def cell_checker(cell, prot_beads, lipid_beads, solv_beads):
                    '''
                    Checks 3 different conditions for a solvent cell
                    1: If the cell is entirely contained within the hydrophobic area of a leaflet.
                        If True:
                            Signal that the cell should be deleted

                    2: If the cell is too close to any bead (protein, lipid and/or solvent).
                        If True:
                            Keep track of which beads are inside the cell.
                            Return the marked beads and signal that subcells should be spawned.

                    3: If the cell is partially contained within the hydrophobic area of a leaflet.
                        If True:
                            Return all original beads for the cell and signal that subcells should be spawned.

                    max_mol_size includes both 'solvation["buffer"]' and 'solvation["kick"]'
                    '''
                    ### Abbreviating some dictionary calls
                    ccx, ccy, ccz = cell["xc"], cell["yc"], cell["zc"]
                    cdx, cdy, cdz = cell["xdim"], cell["ydim"], cell["zdim"]
                    
                    marked_prot_beads = []
                    marked_lipid_beads = []
                    marked_solv_beads = []
                    internal = False
                    proximity = False
                    hydrophobic = False
                    
#                     solvent_buffer = solvation["bead_radius"] + solvation["gridres"]
                    solvent_buffer = solvation["gridres"] + solvation["buffer"]
                    
                    #################################################################
                    ### CHECK IF WHOLE SOLVENT CELL CONTAINED IN HYDROPHOBIC AREA ###
                    #################################################################
                    if len(self.LEAFLETS) != 0:
                        '''
                        pos: the "positive" (coordinate) side of a hydrophobic box in a given dimension (x/y/z)
                        neg: the "negative" (coordinate) side of a hydrophobic box in a given dimension (x/y/z)
                        '''
                        for leaflet_i, (leaflet_nr, leaflet) in enumerate(self.LEAFLETS.items()):
                            lcx, lcy, lcz = leaflet["center"] # Center of leaflet on given axis
                            llx, lly = leaflet["x"], leaflet["y"] # Length of leaflet in given axis
                            lhydrophob = leaflet["lipid_dimensions"]["zUnderLipCen"] # Height of hydrophobic volume

                            xhydr_pos = ccx + cdx/2 < lcx + llx/2 + solvent_buffer
                            xhydr_neg = ccx - cdx/2 > lcx - llx/2 - solvent_buffer
                            yhydr_pos = ccy + cdy/2 < lcy + lly/2 + solvent_buffer
                            yhydr_neg = ccy - cdy/2 > lcy - lly/2 - solvent_buffer
                            if leaflet["HG_direction"] == "up":
                                zhydr_pos = ccz + cdz/2 < lcz + lhydrophob + solvent_buffer
                                zhydr_neg = ccz - cdz/2 > lcz - solvent_buffer
                            if leaflet["HG_direction"] == "down":
                                zhydr_pos = ccz + cdz/2 < lcz + solvent_buffer
                                zhydr_neg = ccz - cdz/2 > lcz - lhydrophob - solvent_buffer
                            if all([xhydr_pos, xhydr_neg, yhydr_pos, yhydr_neg, zhydr_pos, zhydr_neg]):
                                return "hydrophobic", [], [], []

                    ##############################################
                    ### CHECK FOR BEADS INSIDE OF SOLVENT CELL ###
                    ##############################################
                    ### Protein beads
                    for bcx, bcy, bcz in prot_beads:
                        ### Checks if any non-solvent bead is within the solvent cell
                        xinternal = math.dist((ccx,), (bcx,)) <= cdx/2 + solvent_buffer
                        yinternal = math.dist((ccy,), (bcy,)) <= cdy/2 + solvent_buffer
                        zinternal = math.dist((ccz,), (bcz,)) <= cdz/2 + solvent_buffer
                        if all([xinternal, yinternal, zinternal]):
                            marked_prot_beads.append((bcx, bcy, bcz))
                            internal = True

                    ### Lipid beads
                    for bcx, bcy, bcz in lipid_beads:
                        ### Checks if any non-solvent bead is within the solvent cell
                        xinternal = math.dist((ccx,), (bcx,)) <= cdx/2 + solvent_buffer
                        yinternal = math.dist((ccy,), (bcy,)) <= cdy/2 + solvent_buffer
                        zinternal = math.dist((ccz,), (bcz,)) <= cdz/2 + solvent_buffer

                        ### Old code using leaflet grid point data instead of lipid bead data
                        ### Keep it around in case it is needed later
        #                 if grid_point[4] == "up":
        #                     zinternal_head_top = grid_point[2] + grid_point[3] + solvation["buffer"] > cell["zc"] + cell["zdim"] / 2
        #                     zinternal_head_bot = grid_point[2] + grid_point[3] + solvation["buffer"] > cell["zc"] - cell["zdim"] / 2
        #                     zinternal_tail_top = grid_point[2] - solvation["buffer"] < cell["zc"] + cell["zdim"] / 2
        #                     zinternal_tail_bot = grid_point[2] - solvation["buffer"] < cell["zc"] - cell["zdim"] / 2
        #                     zinternal_sum = sum([zinternal_head_top, zinternal_head_bot, zinternal_tail_top, zinternal_tail_bot])
        #                     if zinternal_sum >= 3 or zinternal_head_top + zinternal_tail_bot == 2:
        #                         zinternal = True
        #                     else:
        #                         zinternal = False
        #                 elif grid_point[4] == "down":
        #                     zinternal_head_top = grid_point[2] - grid_point[3] - solvation["buffer"] < cell["zc"] + cell["zdim"] / 2
        #                     zinternal_head_bot = grid_point[2] - grid_point[3] - solvation["buffer"] < cell["zc"] - cell["zdim"] / 2
        #                     zinternal_tail_top = grid_point[2] + solvation["buffer"] > cell["zc"] + cell["zdim"] / 2
        #                     zinternal_tail_bot = grid_point[2] + solvation["buffer"] > cell["zc"] - cell["zdim"] / 2
        #                     zinternal_sum = sum([zinternal_head_top, zinternal_head_bot, zinternal_tail_top, zinternal_tail_bot])
        #                     if zinternal_sum >= 3 or zinternal_tail_top + zinternal_head_bot == 2:
        #                         zinternal = True
        #                     else:
        #                         zinternal = False

                        if all([xinternal, yinternal, zinternal]):
                            marked_lipid_beads.append((bcx, bcy, bcz))
                            internal = True

                    ### Solvent beads (solvent created using earlier solvent commands)
                    for bcx, bcy, bcz in solv_beads:
                        ### Checks if any non-(current)solvent bead is within the solvent cell
                        xinternal = math.dist((ccx,), (bcx,)) <= cdx/2 + solvent_buffer
                        yinternal = math.dist((ccy,), (bcy,)) <= cdy/2 + solvent_buffer
                        zinternal = math.dist((ccz,), (bcz,)) <= cdz/2 + solvent_buffer
                        if all([xinternal, yinternal, zinternal]):
                            marked_solv_beads.append((bcx, bcy, bcz))
                            internal = True

                    if internal:
                        return "internal", marked_prot_beads, marked_lipid_beads, marked_solv_beads

                    ####################################################################################################
                    ### Check if partially within hydrophobic area without containing any internal non-solvent beads ###
                    ####################################################################################################
                    if len(self.LEAFLETS) != 0:
                        '''
                        Finds solvent cells that are only partially contained in a hydrophobic box.
                        Used to find edge cases where cells are partially contained in a hydrophobic box without overlapping with a particle.

                        pos: the "positive" (coordinate) side of a hydrophobic box in a given dimension (x/y/z)
                        neg: the "negative" (coordinate) side of a hydrophobic box in a given dimension (x/y/z)

                        top: the "positive" (coordinate) side of a solvent cell in a given dimension (x/y/z)
                        bot: the "negative" (coordinate) side of a solvent cell in a given dimension (x/y/z)
                        '''
                        for leaflet_i, (leaflet_nr, leaflet) in enumerate(self.LEAFLETS.items()):
                            lcx, lcy, lcz = leaflet["center"] # Center of leaflet on given axis
                            llx, lly = leaflet["x"], leaflet["y"] # Length of leaflet in given axis
                            lhydrophob = leaflet["lipid_dimensions"]["zUnderLipCen"] # Height of hydrophobic volume
                            
                            xhydr_pos_top = ccx + cdx/2 < lcx + llx / 2 + solvent_buffer
                            xhydr_pos_bot = ccx + cdx/2 > lcx - llx / 2 - solvent_buffer
                            xhydr_neg_top = ccx - cdx/2 < lcx + llx / 2 + solvent_buffer
                            xhydr_neg_bot = ccx - cdx/2 > lcx - llx / 2 - solvent_buffer

                            yhydr_pos_top = ccy + cdy/2 < lcy + lly / 2 + solvent_buffer
                            yhydr_pos_bot = ccy + cdy/2 > lcy - lly / 2 - solvent_buffer
                            yhydr_neg_top = ccy - cdy/2 < lcy + lly / 2 + solvent_buffer
                            yhydr_neg_bot = ccy - cdy/2 > lcy - lly / 2 - solvent_buffer

                            if leaflet["HG_direction"] == "up":
                                zhydr_pos_top = ccz + cdz/2 < lcz + lhydrophob + solvent_buffer
                                zhydr_pos_bot = ccz + cdz/2 > lcz - solvent_buffer
                                zhydr_neg_top = ccz - cdz/2 < lcz + lhydrophob + solvent_buffer
                                zhydr_neg_bot = ccz - cdz/2 > lcz - solvent_buffer
                            if leaflet["HG_direction"] == "down":
                                zhydr_pos_top = ccz + cdz/2 < lcz + solvent_buffer
                                zhydr_pos_bot = ccz + cdz/2 > lcz - lhydrophob - solvent_buffer
                                zhydr_neg_top = ccz - cdz/2 < lcz + solvent_buffer
                                zhydr_neg_bot = ccz - cdz/2 > lcz - lhydrophob - solvent_buffer

                            xhydr_pos = xhydr_pos_top and xhydr_pos_bot
                            xhydr_neg = xhydr_neg_top and xhydr_neg_bot
                            yhydr_pos = yhydr_pos_top and yhydr_pos_bot
                            yhydr_neg = yhydr_neg_top and yhydr_neg_bot
                            zhydr_pos = zhydr_pos_top and zhydr_pos_bot
                            zhydr_neg = zhydr_neg_top and zhydr_neg_bot

                            partially_contained = (xhydr_pos or xhydr_neg) and (yhydr_pos or yhydr_neg) and (zhydr_pos or zhydr_neg)
                            if partially_contained:
                                return "internal", [], [], []

                    ### If no hits from above checks, then assume the cell is free
                    return "free", [], [], []

                def subcell_spawner_ax(cell, ax):
                    '''
                    Divides a specific axis into two subcells
                    '''
                    new_axn = list(filter(lambda x: x != 0, [math.floor(cell[ax + "n"] / 2), math.ceil(cell[ax + "n"] / 2)]))
                    if len(new_axn) == 1:
                        new_axdim = [cell[ax + "dim"]]
                        new_axc = [cell[ax + "c"]]
                    else:
                        new_axdim = [cell[ax + "dim"] / cell[ax + "n"] * n for n in new_axn]
                        new_axc = [cell[ax + "c"] - sign * cell[ax + "dim"] / 2 + sign * axdim / 2 for axdim, sign in zip(new_axdim, [+1, -1])]
                    return new_axc, new_axdim, new_axn

                def subcell_spawner(cell, marked_prot_beads, marked_lipid_beads, marked_solv_beads):
                    '''
                    Spawns smaller subcell(s) from a primary cell.
                    Attempts to divide cell into two boxes in each dimension.
                    '''

                    new_axcs = []
                    new_axdims = []
                    new_axns = []

                    for ax in ["x", "y", "z"]:
                        output = subcell_spawner_ax(cell, ax)
                        new_axcs.append(output[0])
                        new_axdims.append(output[1])
                        new_axns.append(output[2])

                    subcells = []
                    for xc, xdim, xn in zip(new_axcs[0], new_axdims[0], new_axns[0]):
                        for yc, ydim, yn in zip(new_axcs[1], new_axdims[1], new_axns[1]):
                            for zc, zdim, zn in zip(new_axcs[2], new_axdims[2], new_axns[2]):
                                subcells.append({
                                    "xc": xc, "xdim": xdim, "xn": xn,
                                    "yc": yc, "ydim": ydim, "yn": yn,
                                    "zc": zc, "zdim": zdim, "zn": zn,
                                })

                    return [(subcell, marked_prot_beads, marked_lipid_beads, marked_solv_beads) for subcell in subcells]

                ##############################################
                ### RUNNING SOLVENT OPTIMIZATION ALGORITHM ###
                ##############################################

                ### Finds the maximum size of molecules used as solvent/ions
                ### Also includes buffer/kick size to prevent edge overlap cases
                max_mol_size = max([max([math.dist(bead1, bead2) for bead1 in molecule for bead2 in molecule]) for molecule in solvent_molecules]) #+ solvation["kick"] + solvation["buffer"]
                
                ### if the maximum molecule size is bigger than the designated grid resolution then change the gridres
                if max_mol_size*1.1 >= solvation["gridres"]:
                    gridres = max_mol_size * 1.1 # 10% larger than largest molecule
                    self.print_term("NOTE: Requested solvent is too large for the grid resolution. Adjusting grid resolution to prevent solvent molecule overlaps.", warn = True)
                    self.print_term("Original grid resolution was:", solvation["gridres"], warn = True)
                    self.print_term("New grid resolution is:      ", gridres, warn = True)
                    self.print_term("", warn = True)
                else:
                    gridres = solvation["gridres"]
                
                if solvation["kick"]*1.1 >= gridres/2:
                    self.print_term("NOTE: Kick is too large for grid resolution. Adjusting grid resolution to prevent solvent molecule overlaps.", warn = True)
                    self.print_term("Original grid resolution was:", gridres, warn = True)
                    self.print_term("Original kick was:           ", solvation["kick"], warn = True)
                    gridres = gridres + (solvation["kick"] - gridres/2)*1.1 # 10% extra
                    self.print_term("New grid resolution is:      ", gridres, warn = True)
                    self.print_term("", warn = True)
                

                start_cell = {
                    "xc": 0,
                    "yc": 0,
                    "zc": 0,
                    "xdim": self.pbcx,
                    "ydim": self.pbcy,
                    "zdim": self.pbcz,
                    "xn": int(self.pbcx / gridres),
                    "yn": int(self.pbcy / gridres),
                    "zn": int(self.pbcz / gridres),
                }
                self.print_term("    ", "Beads used for solvent cell overlaps:")
                self.print_term("    ", "    ", "Number of lipid beads:             ", len(lipid_beads_for_cell_checker))
                self.print_term("    ", "    ", "Number of protein beads:           ", len(prot_beads_for_cell_checker))
                self.print_term("    ", "    ", "Number of (previous) solvent beads:", len(solv_beads_for_cell_checker))
                self.print_term("    ", "    ", "Total number of beads:             ", len(lipid_beads_for_cell_checker) + len(prot_beads_for_cell_checker) + len(solv_beads_for_cell_checker))
                self.print_term("")
                self.print_term("    ", "Solvent cell calculations:")
                self.print_term("    ", "    ", "Buffer space 'buffer':                         ", solvation["buffer"], "[Å]")
                self.print_term("    ", "    ", "Grid resolution 'gridres':                     ", gridres, "[Å]")
                self.print_term("    ", "    ", "Minimum implicit solvent 'min_cellsize':       ", solvation["min_cellsize"])
                self.print_term("    ", "    ", "Max potential solvent particles in x dimension:", start_cell["xn"])
                self.print_term("    ", "    ", "Max potential solvent particles in y dimension:", start_cell["yn"])
                self.print_term("    ", "    ", "Max potential solvent particles in z dimension:", start_cell["zn"])
                self.print_term("    ", "    ", "Max potential total solvent:                   ", start_cell["xn"] * start_cell["yn"] * start_cell["zn"])
                self.print_term("    ", "    ", "    ", "Iteratively reducing solvent cell dimensions:")

                free_cells = []
                hydr_cells = []
                internal_cells = [(start_cell, prot_beads_for_cell_checker, lipid_beads_for_cell_checker, solv_beads_for_cell_checker)]

                iterations = 0
                max_iterations = 1
                max_n = max(start_cell["xn"], start_cell["yn"], start_cell["zn"])
                ### Counting the maximum possible number of iterations
                while max_n > 1:
                    max_iterations += 1
                    max_n = int(max_n / 2 + 0.5)


                self.print_term("Iteration nr:", iterations, "    internal:", len(internal_cells), "    free:", len(free_cells), "    hydrophobic:", len(hydr_cells), debug = True)
                self.print_term("xn", "yn", "zn", "  :  ", "xdim", "ydim", "zdim", "  :  ", "xc", "yc", "zc", debug = True,)
                self.print_term(
                    internal_cells[0][0]["xn"]  , internal_cells[0][0]["yn"]  , internal_cells[0][0]["zn"]  , "  :  ",
                    internal_cells[0][0]["xdim"], internal_cells[0][0]["ydim"], internal_cells[0][0]["zdim"], "  :  ",
                    internal_cells[0][0]["xc"]  , internal_cells[0][0]["yc"]  , internal_cells[0][0]["zc"]  ,
                    debug = True,
                )
                while len(internal_cells) > 0:
                    '''
                    Continuesly goes through all remaning non-free cells.
                    Checks them for overlaps and splits them if needed.
                    '''
                    iterations += 1
                    self.print_term("    ", "    ", "    ", "    ", "Iteration:", str(iterations), "/", str(max_iterations))
                    new_internal_cells = []
                    self.print_term("----------------------------------------------------------", debug = True)
                    for cell_i, (cell, prot_beads, lipid_beads, solv_beads) in enumerate(internal_cells):
                        check, marked_prot_beads, marked_lipid_beads, marked_solv_beads = cell_checker(cell, prot_beads, lipid_beads, solv_beads)
                        if check == "free":
                            free_cells.append(internal_cells[cell_i][0])
                        if check == "hydrophobic":
                            hydr_cells.append(internal_cells[cell_i][0])
                        if check == "internal":
                            if not all([axn <= solvation["min_cellsize"] for axn in [cell["xn"], cell["yn"], cell["zn"]]]):
                                '''Cells containing overlaps but are indivisible in all axes, should be removed'''
                                new_internal_cells.extend(subcell_spawner(cell, marked_prot_beads, marked_lipid_beads, marked_solv_beads))

                    ### For the next iteration
                    internal_cells = new_internal_cells

                    self.print_term("Iteration nr:", iterations, "    internal:", len(internal_cells), "    free:", len(free_cells), "    hydrophobic:", len(hydr_cells), debug = True)
                    if internal_cells != []:
                        self.print_term("xn", "yn", "zn", "  :  ", "xdim", "ydim", "zdim", "  :  ", "xc", "yc", "zc", debug = True,)
                        for i in [0, -1]:
                            self.print_term(
                                internal_cells[i][0]["xn"]  , internal_cells[i][0]["yn"]  , internal_cells[i][0]["zn"]  , "  :  ",
                                internal_cells[i][0]["xdim"], internal_cells[i][0]["ydim"], internal_cells[i][0]["zdim"], "  :  ",
                                internal_cells[i][0]["xc"]  , internal_cells[i][0]["yc"]  , internal_cells[i][0]["zc"]  ,
                                debug = True,
                            )

                if self.debug_prints:
                    counts = {}
                    for cell in sorted(free_cells, key = lambda cell: cell["xn"] * cell["yn"] * cell["zn"]):
                        prod = cell["xn"] * cell["yn"] * cell["zn"]
                        if prod not in counts.keys():
                            counts[prod] = 1
                        else:
                            counts[prod] += 1
                    for key, vals in counts.items():
                        self.print_term("xn*yn*zn product and number of cells", key, "    ", vals, debug = True)

                self.print_term("Lengths:", debug = True)
                self.print_term("    ", "Internal cells:", len(internal_cells), debug = True)
                self.print_term("    ", "Free cells:", len(free_cells), debug = True)
                self.print_term("    ", "Hydrophobic cells:", len(hydr_cells), debug = True)

                solv_grid_3D = [
                    [
                        (
                            cell["xc"] - cell["xdim"] / 2 + (cell["xdim"] / cell["xn"]) * (xn - 0.5),
                            cell["yc"] - cell["ydim"] / 2 + (cell["ydim"] / cell["yn"]) * (yn - 0.5),
                            cell["zc"] - cell["zdim"] / 2 + (cell["zdim"] / cell["zn"]) * (zn - 0.5),
                        )
                        for xn in range(1, cell["xn"] + 1) for yn in range(1, cell["yn"] + 1) for zn in range(1, cell["zn"] + 1)
                    ] for cell in free_cells
                ]
                solv_grid_3D = [item for sublist in solv_grid_3D for item in sublist]

                self.print_term("    ", "Final number of 3D grid points available for solvent placement:", len(solv_grid_3D))

                ### Assert that there are enough grid points for solvent
                assert len(collected_solvent) <= len(solv_grid_3D), "Not enough grid points for solvent. Use a lower 'gridres' value."

                ################################
                ### INSERT SOLVENT MOLECULES ###
                ################################

                def position_fixer(coords, ax):
                    '''
                    Checking if any beads are outside pbc and moves them if they are
                    If all beads are outside: Moved to other side of pbc
                    If at least 1 bead is inside: All beads are pushed slightly so that they are inside the pbc
                    '''
                    diffs = [coord - (np.sign(coord) * self.pbc_box[ax] / 2) for coord in coords]
                    checks = [np.sign(coords[i]) == np.sign(diffs[i]) for i in range(len(coords))]
                    if all(checks):
                        coords = [coord - (np.sign(coord) * self.pbc_box[ax]) for coord in coords]
                    elif any(checks):
                        coords = [coord - (np.sign(coord) * max(abs(min(diffs)), abs(max(diffs)))) for coord in coords]
                    return coords

                ### Shuffles the grid to ensure randomness
                random.shuffle(solv_grid_3D)
                ### random.sample extracts k random elements from the list without duplicates
                random_grid_points = random.sample(solv_grid_3D, k = len(collected_solvent))

                grid_solvated = []
                counter = 0
                generated_spots = []
                self.print_term("    ", "Inserting", len(collected_solvent), "solvent molecules into random grid points:")
                for (sname, stype, ssize), (gx, gy, gz) in zip(collected_solvent, random_grid_points):
                    '''
                    Finds a 3D-grid point for the solvent and places it there.
                    '''
                    counter += 1
                    if counter % 25000 == 0:
                        self.print_term("    ", "    ", "Currently at solvent number:", counter)
                    sdata = solvation[stype][sname]
                    sbeads = sdata["beads"]
                    scharge = sdata["charge"]
                    sx, sy, sz = sdata["x"].copy(), sdata["y"].copy(), sdata["z"].copy()
                    rx, ry, rz = random.uniform(0, 360), random.uniform(0, 360), random.uniform(0, 360)
                    for j, (x, y, z) in enumerate(zip(sx, sy, sz)):
                        sx[j], sy[j], sz[j] = self.rotate_point(x, y, z, rx, ry, rz)
                        kx = random.uniform(-solvation["kick"], solvation["kick"])
                        ky = random.uniform(-solvation["kick"], solvation["kick"])
                        kz = random.uniform(-solvation["kick"], solvation["kick"])
                        sx[j], sy[j], sz[j] = sx[j] + kx + gx, sy[j] + ky + gy, sz[j] + kz + gz

                    sx = position_fixer(sx, 0)
                    sy = position_fixer(sy, 1)
                    sz = position_fixer(sz, 2)

                    ### Order that beads should be written in for structure file and topology file
                    if stype == "solvent":
                        order = 1
                    if stype == "pos_ions":
                        order = 2
                    if stype == "neg_ions":
                        order = 3

                    grid_solvated.append({
                        "order": order,
                        "name": sname,
                        "type": stype,
                        "beads": sbeads,
                        "charge": scharge,
                        "coords": list(zip(sx, sy, sz)),
                    })
                ### Orders the solvent key:value pairs to ensure topology-structure file match
                self.SOLVATIONS[solvation_nr]["grid"] = sorted(grid_solvated, key=lambda g: (g["order"], g["name"]))

                ### Remembers the beads for subsequent solvation commands
                solv_beads_for_cell_checker.extend([(x, y, z) for vals in solvation["grid"] for x, y, z in vals["coords"]])

                self.SOLVATIONS[solvation_nr]["solv_count"] = {}
                for grid_point_3D in self.SOLVATIONS[solvation_nr]["grid"]:
                    '''
                    Counts all solvent for this specific solvation command
                    '''
                    if (grid_point_3D["name"], grid_point_3D["type"], grid_point_3D["charge"]) not in solvation["solv_count"].keys():
                        self.SOLVATIONS[solvation_nr]["solv_count"][(grid_point_3D["name"], grid_point_3D["type"], grid_point_3D["charge"])] = 1
                    else:
                        self.SOLVATIONS[solvation_nr]["solv_count"][(grid_point_3D["name"], grid_point_3D["type"], grid_point_3D["charge"])] += 1

            SYS_solv_count = {}
            for solvation_nr, solvation in self.SOLVATIONS.items():
                '''
                Counts all solvent in the system
                '''
                for (sname, stype, scharge), count in solvation["solv_count"].items():
                    if (sname, stype, scharge) not in SYS_solv_count.keys():
                        SYS_solv_count[(sname, stype, scharge)] = count
                    else:
                        SYS_solv_count[(sname, stype, scharge)] += count

            if self.extra_info:
                '''
                Prints system-wide solvent information to the terminal
                '''
                self.print_term()
                SYS_headers = ["Name", "Totals", "Total %", "Molarity(box)", "Molarity(free)", "Molarity(solvent)", "Charge"]
                SYS_names, SYS_charges, SYS_counts = list(zip(*[(sname, scharge, scount) for (sname, stype, scharge), scount in SYS_solv_count.items()]))
                SYS_tot = sum(SYS_counts)
                SYS_percentages = [round(count / SYS_tot * 100, 3) for count in SYS_counts]
                SYS_box_molarity = [round(scount / (N_A * box_volume), 3) for scount in SYS_counts]
                SYS_free_molarity = [round(scount / (N_A * box_free_volume), 3) for scount in SYS_counts]
                try:
                    SYS_solv_molarity = [round(scount / (N_A * solv_volume), 3) for scount in SYS_counts]
                except:
                    SYS_solv_molarity = ["NaN" for scount in SYS_counts]
                SYS_printer = [tuple(SYS_headers)] + list(zip(SYS_names, SYS_counts, SYS_percentages, SYS_box_molarity, SYS_free_molarity, SYS_solv_molarity, SYS_charges))
                SYS_max_lengths = [len(str(max(data, key = lambda d: len(str(d))))) for data in zip(*SYS_printer)]
                self.print_term("Solvent data for whole system")
                for L0, L1, L2, L3, L4, L5, L6 in SYS_printer:
                    self.print_term(
                        "    ",
                        "    ",
                        '{0: <{L}}'.format(L0, L = SYS_max_lengths[0]), ":",
                        '{0: <{L}}'.format(L1, L = SYS_max_lengths[1]), ":",
                        '{0: <{L}}'.format(L2, L = SYS_max_lengths[2]), ":",
                        '{0: <{L}}'.format(L3, L = SYS_max_lengths[3]), ":",
                        '{0: <{L}}'.format(L4, L = SYS_max_lengths[4]), ":",
                        '{0: <{L}}'.format(L5, L = SYS_max_lengths[5]), ":",
                        '{0: <{L}}'.format(L6, L = SYS_max_lengths[6]),
                    )
            self.print_term("------------------------------ SOLVATION COMPLETE", "\n")
    
    ###############
    ### PLOTTER ###
    ###############
    def plotter(self):
        if self.plots_requested and len(list(self.LEAFLETS.keys())) > 0:
            self.print_term("-------------------------------------------------------------")
            self.print_term("Plotting lipid-protein data for each individual leaflet into:", self.PLOT_cmd)

            nrows = 1
            ncols = len(list(self.LEAFLETS.keys()))

            # fig, axes = plt.subplots(nrows = nrows, ncols = ncols, figsize=(ncols * 12, nrows * 12), dpi = 120, squeeze = False)
            fig, axes = plt.subplots(nrows = nrows, ncols = ncols, figsize=(ncols * self.pbc_box[0] / 10, nrows * self.pbc_box[1] / 10), dpi = 120, squeeze = False)

            def alpha_normalizer(values, des_low, des_upp, act_low, act_upp):
                '''
                Adapted from:
                https://stackoverflow.com/questions/48109228/normalizing-data-to-certain-range-of-values
                '''
                return [des_low + (x - act_low) * (des_upp - des_low) / (act_upp - act_low) for x in values]

            for col, (leaflet_nr, leaflet) in enumerate(self.LEAFLETS.items()):
                axes[0, col].set(
                    xticks = np.arange(- self.pbc_box[0] / 2, self.pbc_box[0] / 2 + 1, 5),
                    yticks = np.arange(- self.pbc_box[1] / 2, self.pbc_box[1] / 2 + 1, 5),
                    xlim = (- self.pbc_box[0] / 2, self.pbc_box[0] / 2),
                    ylim = (- self.pbc_box[1] / 2, self.pbc_box[1] / 2),
                )

                if len(self.PROTEINS) != 0:
                    for protein_i in range(len(self.plot_data["Prot_ConvexHull"][col])):
                        ### Convex hull lines
                        for prot_CH in self.plot_data["Prot_ConvexHull"][col][protein_i]:
                            x, y = list(zip(*prot_CH))
                            x, y = x + (x[0],), y + (y[0],)
                            axes[0, col].plot(
                                x, y, linestyle = "-", color = "black"
                            )

                        ### Concave hull lines
                        for prot in self.plot_data["Prot_ConcaveHulls_Points"][col][protein_i]:
                            for cluster in prot:
                                CHx, CHy = list(zip(*cluster))
                                axes[0, col].plot(CHx, CHy, linestyle = "-", color = "firebrick")# marker = "-", markersize=6, markeredgecolor="black", markerfacecolor="black")

                        ## Concave hull Polygon area and buffer
                        for prot in self.plot_data["Prot_ConcaveHull_Polygon_Buffered"][col][protein_i]:
                            for cluster in prot:
                                for poly in cluster:
                                    shapely.plotting.plot_polygon(poly, ax=axes[0, col], add_points=False, color= "red", alpha=0.1)

                        ### Protein flattened
                        for Protein in self.plot_data["Prot_Flattened_points"][col][protein_i]:
                            prot_zmax = max([z for x, y, z in Protein])
                            prot_zmin = min([z for x, y, z in Protein])
                            for x, y, z in Protein:
                                alp = alpha_normalizer([z], 0.5, 1.0, prot_zmin, prot_zmax)[0]
                                axes[0, col].plot(
                                    x, y, marker = "o", markersize=6,
                                    markeredgecolor="darkgoldenrod", markeredgewidth=0.0,
                                    markerfacecolor="darkgoldenrod",
                                    alpha = alp
                                )

                ### Lipid grid points
                lipid_radius = leaflet["lipid_dimensions"]["lipid_radius"]
                ### After protein adjustment
                vals_adjusted = [(p["x"], p["y"], p["z"], p["xdim"], p["ydim"], p["Overlap"], p["Inside_Poly"], p["Inside_ConvexHull"], p["Removed"]) for key, p in leaflet["grid_adjusted"].items()]
                vals_overlap = [(p["x"], p["y"], p["z"], p["xdim"], p["ydim"], p["Overlap"], p["Inside_Poly"], p["Inside_ConvexHull"], p["Removed"]) for key, p in leaflet["grid_overlapping"].items()]
                vals_removed = [(p["x"], p["y"], p["z"], p["xdim"], p["ydim"], p["Overlap"], p["Inside_Poly"], p["Inside_ConvexHull"], p["Removed"]) for key, p in leaflet["grid_removed"].items()]
                vals = vals_adjusted + vals_overlap + vals_removed

                for grid_point_nr, (x, y, z, xl, yl, Overlap, Inside_Poly, Inside_ConvexHull, Removed) in enumerate(vals):
                    lw = 2
                    if Overlap and Inside_Poly:
                        color = "red"
                    elif Overlap and not Inside_Poly:
                        color = "orange"
                    elif Inside_Poly and not Overlap:
                        color = "purple"
                    elif Inside_ConvexHull and not Inside_Poly:
                        color = "green"
                    elif Removed:
                        color = "darkcyan"
                    else:
                        color = "grey"

                    axes[0, col].plot(x, y, marker = ".", markersize=6, markeredgecolor=color, markerfacecolor=color)

                    axes[0, col].add_patch(
                        patches.Circle(
                            (x, y),
                            radius = lipid_radius,
                            edgecolor = color,
                            fill = False,
                            linewidth = lw,
                        )
                    )

                    axes[0, col].add_patch(
                        patches.Rectangle(
                            (x - 0.5 * xl, y - 0.5 * yl),
                            xl,
                            yl,
                            edgecolor = "black",
                            fill = False,
                        )
                    )

                ### Lipid data
                lipids = [(p["lipid"]["x"], p["lipid"]["y"], p["lipid"]["z"], p["internal_z"]) for key, p in leaflet["grid_adjusted"].items()]

                ### Plotting beads in lipids
                for x, y, z, internal_z in lipids:
                    alphas = alpha_normalizer(z, 0.5, 1.0, min(z), max(z))
                    for xp, yp, alp in zip(x, y, alphas):
                        axes[0, col].plot(
                            xp, yp, marker = "o",markersize=3,
                            markeredgecolor="blue", markeredgewidth=0.0,
                            markerfacecolor="blue",
                            alpha = alp
                        )
                ### Normal plotting
                axes[0, col].set(
                    xticks = np.arange(- self.pbc_box[0] / 2, self.pbc_box[0] / 2 + 1, 5),
                    yticks = np.arange(- self.pbc_box[1] / 2, self.pbc_box[1] / 2 + 1, 5),
                    xlim = (- self.pbc_box[0] / 2, self.pbc_box[0] / 2),
                    ylim = (- self.pbc_box[1] / 2, self.pbc_box[1] / 2),
                )
                axes[0, col].set_title(
                    label = "leaf nr: " + str(leaflet_nr),
                    fontsize = 40,
                )
            if self.backup:
                self.backupper(self.PLOT_cmd)
            plt.savefig(self.PLOT_cmd)
            plt.close()
            self.print_term("-------------------------------------------------------------", "\n")
            
    ###############
    ### PICKLER ###
    ###############
    def pickler(self):
        if self.PICKLE_cmd:
            self.print_term("-------------------")
            self.print_term("Pickling data into:", self.PICKLE_cmd)
            self.print_term("-------------------", "\n")
            pickled_data = {
                "system_charge": self.system_charge,
                "pbc_box": self.pbc_box,
            }
            if len(self.LEAFLETS_cmds) != 0:
                pickled_data["leaflets_cmds"] = self.LEAFLETS_cmds
                pickled_data["leaflets"] = self.LEAFLETS
            if len(self.PROTEINS_cmds) != 0:
                pickled_data["proteins_cmds"] = self.PROTEINS_cmds
                pickled_data["proteins"] = self.PROTEINS
            if len(self.SOLVATIONS_cmds) != 0:
                pickled_data["solvations_cmds"] = self.SOLVATIONS_cmds
                pickled_data["solvations"] = self.SOLVATIONS
            if self.PLOT_cmd:
                pickled_data["plot_data"] = self.plot_data
            if self.backup:
                self.backupper(self.PICKLE_cmd)
            with open(self.PICKLE_cmd, 'wb') as f:
                pickle.dump(pickled_data, f)

#####################################################################
########################## HERE BE PARSING ##########################
#####################################################################

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawTextHelpFormatter,
    add_help = False
)

parser.add_argument("-h", "--help", dest = "help")

### Leaflet commands
parser.add_argument("--memb", "-memb", "-membrane", dest = "membrane_cmds", action="append", type=str, default = [], nargs="+")

### Protein commands
parser.add_argument("--prot", "-prot", "-protein", dest = "protein_cmds", action="append", type=str, default = [], nargs="+")

### Solvent commands
parser.add_argument("--solv", "-solv", "-solvation", dest = "solvation_cmds", action="append", type=str, default = [], nargs="+")

### Topology commands
parser.add_argument("--itp_input", "-itp_in", dest = "itp_input_cmds", action="append", type=str, default = [], nargs="+")

### Import commands
parser.add_argument("--solute_input", "-solute_in", dest = "solute_input_cmds", action="append", type=str, default = [], nargs="+")

### Plotting command
parser.add_argument("-plot", dest = "plot_cmd", default = False)

### Pickle commands
parser.add_argument("-pickle", dest = "pickle_cmd", default = False)

### Whether to backup files i they would be overwritten
parser.add_argument("-backup", dest = "backup", default = True)

### Random seed
parser.add_argument("-rand", dest = "randseed", default = False)

### System force field
parser.add_argument("-params", dest = "sys_params", default = "default")

### System name
parser.add_argument("-sn", "-system_name", dest = "system_name", default = "PLACEHOLDER_TITLE")

### pbc box size [nm]
parser.add_argument("-box", "-pbc", dest = "pbc_box", action="extend", type=str, default = [], nargs="+")

### x/y/z size of box [nm]
parser.add_argument("-x", dest = "pbcx", type=str, default = False)
parser.add_argument("-y", dest = "pbcy", type=str, default = False)
parser.add_argument("-z", dest = "pbcz", type=str, default = False)

####################
### OUTPUT FILES ###
####################
### Output pdb/gro file
parser.add_argument("--output_struct", "-out", "-out_sys", "-o", dest = "out_system_file_name", default = "output.pdb")

### Output topology file
parser.add_argument("--output_topol", "-top_out", "-t", dest = "out_topol_file_name", default = "topol.top")

### Output imported file
# parser.add_argument("--output_import", "-imp_out", dest = "output_imported", default = False)

#############################
### PRINTING AND LOG FILE ###
#############################
### Log file
parser.add_argument("--output_log", "-log", "-out_log", dest = "out_log_file_name", default = False)

### Prints
parser.add_argument("--print_quiet", "-quiet", dest = "quiet", default = False)
parser.add_argument("--print_debug", "-debug", dest = "debug", default = False)
parser.add_argument("--print_extra", "-extra", dest = "extra", default = True)
parser.add_argument("--print_warnings", "-warn", dest = "warnings", default = True)

### ### Parser for handling '-f' when importing module to Jupyter
parser.add_argument("-f", dest = "debug_flag_for_jupyter")

args = parser.parse_args()

parse_membrane_cmds      = [" ".join(i) for i in args.membrane_cmds]
parse_protein_cmds      = [" ".join(i) for i in args.protein_cmds]
parse_solvation_cmds   = [" ".join(i) for i in args.solvation_cmds]
parse_itp_input_cmds   = [" ".join(i) for i in args.itp_input_cmds]
parse_solute_input_cmds = [" ".join(i) for i in args.solute_input_cmds]

parse_plot_cmd   = args.plot_cmd
parse_pickle_cmd = args.pickle_cmd
parse_backup     = bool(ast.literal_eval(str(args.backup)))
parse_randseed   = args.randseed
if parse_randseed:
    random.seed(int(parse_randseed))

parse_sys_params  = args.sys_params
parse_pbc_box = args.pbc_box
if parse_pbc_box:
    parse_pbc_box = [ast.literal_eval(str(i)) for i in parse_pbc_box]

parse_pbcx = ast.literal_eval(str(args.pbcx))
parse_pbcy = ast.literal_eval(str(args.pbcy))
parse_pbcz = ast.literal_eval(str(args.pbcz))
if parse_pbcx and parse_pbcy and parse_pbcz:
    parse_pbc_box = [parse_pbcx, parse_pbcy, parse_pbcz]

parse_out_system_file_name = args.out_system_file_name
parse_out_topol_file_name  = args.out_topol_file_name
# parse_output_imported  = args.output_imported
parse_system_name      = args.system_name

parse_out_log_file_name = args.out_log_file_name

parse_quiet        = bool(ast.literal_eval(str(args.quiet)))
parse_debug_prints = bool(ast.literal_eval(str(args.debug)))
parse_extra_info   = bool(ast.literal_eval(str(args.extra)))
parse_warnings     = bool(ast.literal_eval(str(args.warnings)))

if any([i != [] for i in [parse_membrane_cmds, parse_protein_cmds, parse_solvation_cmds]]):
    ENSANE(
        box = parse_pbc_box,
        
        membrane  = parse_membrane_cmds,
        protein   = parse_protein_cmds,
        solvation = parse_solvation_cmds,

        itp_input    = parse_itp_input_cmds,
        solute_input = parse_solute_input_cmds,

        plot       = parse_plot_cmd,
        pickle     = parse_pickle_cmd,
        backup     = parse_backup,
        sys_params = parse_sys_params,

        out_sys     = parse_out_system_file_name,
        out_top     = parse_out_topol_file_name,
        out_log     = parse_out_log_file_name,
        system_name = parse_system_name,

#         output_imported = parse_output_imported,

        quiet = parse_quiet,
        debug = parse_debug_prints,
        extra = parse_extra_info,
        warn  = parse_warnings,

        run = True,
    )

##############################
### HELP FOR THOSE IN NEED ###
##############################

# if args.help != None or len(sys.argv) == 1:
#     print_helper()
#     sys.exit()

#####################################################################
########################## YOU HAVE PARSED ##########################
#####################################################################


