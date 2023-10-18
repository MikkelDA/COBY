import time
import_tic = time.time()

import ast
import itertools
import math
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
import operator
from operator import itemgetter

import matplotlib.pyplot as plt
from matplotlib import patches

import pickle

import_toc = time.time()
import_time = round(import_toc - import_tic, 4)
print("Time spent importing packages:", import_time)

lipid_defs = {}
solvent_defs = {}
ion_defs = {}
prot_defs = {}

### Diacyl glycerols
lipid_type, params = "lipid", "default"
lipid_defs[(lipid_type, params)] = {}
lipid_defs[(lipid_type, params)]["x"] = (    0, .5,  0,  0, .5,  0,  0, .5,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1)
lipid_defs[(lipid_type, params)]["y"] = (    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0)
lipid_defs[(lipid_type, params)]["z"] = (   10,  9,  9,  8,  8,  7,  6,  6,  5,  4,  3,  2,  1,  0,  5,  4,  3,  2,  1,  0)
lipid_defs[(lipid_type, params)]["center"] = 7 # PO4
lipid_defs[(lipid_type, params)]["bd"] = (0.25, 0.25, 0.3)
lipid_defs[(lipid_type, params)]["charges"] = (("NC3", 1), ("NH3", 1), ("PO4", -1))
lipid_defs[(lipid_type, params)]["lipids"] = {      # 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20
### Phospholipids
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

lipid_type, params = "INOSITOLLIPIDS", "default"
lipid_defs[(lipid_type, params)] = {}
lipid_defs[(lipid_type, params)]["x"] = (   .5,  .5,   0,   0,   1, .5,  0,  0,   .5,   0,   0,   0,   0,   0,   0,   1,   1,   1,   1,   1,   1)
lipid_defs[(lipid_type, params)]["y"] = (    0,   0,   0,   0,   0,  0,  0,  0,    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0)
lipid_defs[(lipid_type, params)]["z"] = (    8,   9,   9,   7,  10, 10, 10,  6,    6,   5,   4,   3,   2,   1,   0,   5,   4,   3,   2,   1,   0)
lipid_defs[(lipid_type, params)]["center"] = 7 # CP
lipid_defs[(lipid_type, params)]["bd"] = (0.25, 0.25, 0.3)
lipid_defs[(lipid_type, params)]["charges"] = (("NC3", 1), ("PO4", -1))
lipid_defs[(lipid_type, params)]["lipids"] = {      # 1     2    3    4    5   6   7   8    9    10    11    12    13    14   15    16    17    18    19   20 
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
lipid_type, params = "GLYCOLIPIDS", "default"
lipid_defs[(lipid_type, params)] = {}
lipid_defs[(lipid_type, params)]["x"] = (    0,  .5,   0,   0,  .5,  0,  0, .5,  0,    0,   .5,    0,    0,    0,   0,    0,    0,    0,   .5,   0,   0,   0,   0,   0,   1,   1,   1,   1,   1)
lipid_defs[(lipid_type, params)]["y"] = (    0,   0,   0,   0,   0,  0,  0,  0,  0,    0,    0,    0,   .5,    1,   1,    1,    1,    0,    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0)
lipid_defs[(lipid_type, params)]["z"] = (    6,   7,   7,   8,   9,  9, 10, 11, 11,   12,   13,   13,   10,    9,  10,    8,   11,    5,    5,   4,   3,   2,   1,   0,   4,   3,   2,   1,   0)
lipid_defs[(lipid_type, params)]["center"] = 6
lipid_defs[(lipid_type, params)]["bd"] = (0.25, 0.25, 0.3)
lipid_defs[(lipid_type, params)]["lipids"] = {      # 1     2    3    4    5   6   7   8   9    10    11    12    13    14   15    16    17    18    19   20   21   22   23   24   25   26   27   28   29
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

lipid_type, params = "QUINONES", "default"
lipid_defs[(lipid_type, params)] = {}
lipid_defs[(lipid_type, params)]["x"] = (    0,  .5,   0,    0,   0,   0,   0,   0,   0,    0,    0,    0)
lipid_defs[(lipid_type, params)]["y"] = (    0,   0,   0,    0,   0,   0,   0,   0,   0,    0,    0,    0)
lipid_defs[(lipid_type, params)]["z"] = (    6,   7,   7,   5.5,  5,  4.5,  4,  3.5, 2.5,   2,  1.5,    1)
lipid_defs[(lipid_type, params)]["center"] = 6 # PLQ3
lipid_defs[(lipid_type, params)]["bd"] = (0.25, 0.25, 0.3)
lipid_defs[(lipid_type, params)]["lipids"] = {      # 1     2    3    4    5    6    7    8    9    10    11    12
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
lipid_defs[(lipid_type, params)] = {}
lipid_defs[(lipid_type, params)]["x"] = (   0.5,   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1)
lipid_defs[(lipid_type, params)]["y"] = (     1,   0,  0,  1,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,   0,  0,  1,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1)
lipid_defs[(lipid_type, params)]["z"] = (     8,   7,  6,  6,  5,  4,  3,  2,  1,  0,  5,  4,  3,  2,  1,  0,   7,  6,  6,  5,  4,  3,  2,  1,  0,  5,  4,  3,  2,  1,  0)
lipid_defs[(lipid_type, params)]["center"] = 7 # PO4
lipid_defs[(lipid_type, params)]["bd"] = (0.25, 0.25, 0.3)
lipid_defs[(lipid_type, params)]["lipids"] = {      #  1    2   3   4   5   6   7   8   9  10  11  12  13  14  15  16   17  18  19  20  21  22  23  24  25  26  27  28  29  30  31
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

lipid_type, params = "MYCOLIC ACIDS", "default"
lipid_defs[(lipid_type, params)] = {}
lipid_defs[(lipid_type, params)]["x"] = (      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,    0,    1,    1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,  1,   1,   1)
lipid_defs[(lipid_type, params)]["y"] = (      0,   0,   0,   0,   0,   0,   0,   0,   1,   1,   1,   1,   1,   1,   1,    1,    1,    1,    1,   1,   1,   1,   1,   1,   0,   0,   0,   0,   0,  0,   0,   0)
lipid_defs[(lipid_type, params)]["z"] = (      7,   6,   5,   4,   3,   2,   1,   0,   0,   1,   2,   3,   4,   5,   6,    7,    7,    6,    5,   4,   3,   2,   1,   0,   1,   0,   2,   3,   4,  5,   6,   7)
lipid_defs[(lipid_type, params)]["center"] = 7
lipid_defs[(lipid_type, params)]["bd"] = (0.25, 0.25, 0.3)
lipid_defs[(lipid_type, params)]["lipids"] = {        # 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15    16    17    18    19   20   21   22   23   24   25   26   27   28   29   30   31   32
    ("AMA", "beads"):   ("  -    -    -  C1A  C2A  C3A  C4A  C5A  M1A  C1B  C2B  C3B  C4B    -    -     -     -     -   M1B  C1C  C2C  C3C    -    -  COH  OOH  C1D  C2D  C3D  C4D  C5D  C6D"),
    ("AMA.w", "beads"): ("  -    -    -  C1A  C2A  C3A  C4A  C5A  M1A  C1B  C2B  C3B  C4B    -    -     -     -     -   M1B  C1C  C2C  C3C    -    -  COH  OOH  C1D  C2D  C3D  C4D  C5D  C6D"),
    ("KMA", "beads"):   ("  -    -    -  C1A  C2A  C3A  C4A  C5A  M1A  C1B  C2B  C3B  C4B    -    -     -     -     -   M1B  C1C  C2C  C3C    -    -  COH  OOH  C1D  C2D  C3D  C4D  C5D  C6D"),
    ("MMA", "beads"):   ("  -    -    -  C1A  C2A  C3A  C4A  C5A  M1A  C1B  C2B  C3B  C4B    -    -     -     -     -   M1B  C1C  C2C  C3C    -    -  COH  OOH  C1D  C2D  C3D  C4D  C5D  C6D"),
}

### Sterols # These are from the cholesterol paper (https://doi.org/10.1021/acs.jctc.3c00547)
lipid_type, params = "sterol", "default"
lipid_defs[(lipid_type, params)] = {}
lipid_defs[(lipid_type, params)]["x"] = (       0,   1,   0,   0,  1,   0, 0.5, 0.5,   0,  0)
lipid_defs[(lipid_type, params)]["y"] = (       0,   0,   0,   0,  0,   0,   0,   0,   0,  0)
lipid_defs[(lipid_type, params)]["z"] = (     5.3, 4.5, 3.9, 3.3,  3, 2.6, 4.5, 2.6, 1.4,  0)
lipid_defs[(lipid_type, params)]["center"] = 4.9 # Between ROH and R1
lipid_defs[(lipid_type, params)]["bd"] = (0.25, 0.25, 0.3)
lipid_defs[(lipid_type, params)]["lipids"] = {
    ("CHOL", "beads"): (" ROH  R1  R2  R3  R4  -  R5  R6  C1  C2 "),
}

### Hopanoids
lipid_type, params = "Hopanoids", "default"
lipid_defs[(lipid_type, params)] = {}
lipid_defs[(lipid_type, params)]["x"] = (     0,  0,  0,  0, 0.5,-0.5,   0,   0, 0.5, 0.5,   0,   0,   0,   0,  0,  0,  0,  0)
lipid_defs[(lipid_type, params)]["y"] = (     0,  0,  0,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  0)
lipid_defs[(lipid_type, params)]["z"] = (     0,  0,  0,  0, 0.5, 1.4, 2.6,   3, 3.3, 3.9, 4.5, 5.0, 5.5, 6.0,  0,  0,  0,  0) 
lipid_defs[(lipid_type, params)]["center"] = 4.9
lipid_defs[(lipid_type, params)]["bd"] = (0.25, 0.25, 0.3)
lipid_defs[(lipid_type, params)]["lipids"] = {
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
params = "default"
solvent_defs[params] = {
    "W" : {"beads": "W",  "x": (0,), "y": (0,), "z": (0,), "solvcount": 4, "density": 0.99669, "molar_mass": 18.01528},
    "SW": {"beads": "SW", "x": (0,), "y": (0,), "z": (0,), "solvcount": 3, "density": 0.99669, "molar_mass": 18.01528},
    "TW": {"beads": "TW", "x": (0,), "y": (0,), "z": (0,), "solvcount": 2, "density": 0.99669, "molar_mass": 18.01528},
}
### Amino acids
solvent_defs[params].update({
    "GLY": {"beads": ("BB",), "solvcount": 1, "x": (0,), "y": (0,), "z": (0,)},
    "ALA": {"beads": ("BB",), "solvcount": 1, "x": (0,), "y": (0,), "z": (0,)},
    
    "ASN": {"beads": ("BB", "SC1"), "solvcount": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
    "ASP": {"beads": ("BB", "SC1"), "solvcount": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0), "charge": -1},
    "GLU": {"beads": ("BB", "SC1"), "solvcount": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0), "charge": -1},
    "GLN": {"beads": ("BB", "SC1"), "solvcount": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
    "LEU": {"beads": ("BB", "SC1"), "solvcount": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
    "ILE": {"beads": ("BB", "SC1"), "solvcount": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
    "VAL": {"beads": ("BB", "SC1"), "solvcount": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
    "SER": {"beads": ("BB", "SC1"), "solvcount": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
    "THR": {"beads": ("BB", "SC1"), "solvcount": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
    "CYS": {"beads": ("BB", "SC1"), "solvcount": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
    "LYS": {"beads": ("BB", "SC1"), "solvcount": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0), "charge": 1},
    "PRO": {"beads": ("BB", "SC1"), "solvcount": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
    "HYP": {"beads": ("BB", "SC1"), "solvcount": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
    
    "ARG": {"beads": ("BB", "SC1", "SC2"), "solvcount": 1, "x": (0.25, 0, -0.25), "y": (0, 0, 0.125), "z": (0, 0, 0), "charge": 1},
    "PHE": {"beads": ("BB", "SC1", "SC2", "SC3"), "solvcount": 1, "x": (0.25, 0, -0.25, -0.25), "y": (0, 0, 0.125, -0.125), "z": (0, 0, 0, 0)},
    "TYR": {"beads": ("BB", "SC1", "SC2", "SC3"), "solvcount": 1, "x": (0.25, 0, -0.25, -0.25), "y": (0, 0, 0.125, -0.125), "z": (0, 0, 0, 0)},
    "TRP": {"beads": ("BB", "SC1", "SC2", "SC3", "SC4"), "solvcount": 1, "x": (0.25, 0.25, 0, 0, -0.25), "y": (0.125, 0, -0.125, 0.125, 0), "z": (0, 0, 0, 0, 0)},
})

### Example of a multi-residue solvent molecule
# solvent_defs[params].update({
#     "LIG1": {
#         "residues": [
#             {
#                 "resname": "RES1",
#                 "beads": ("R1A1", "R1A2", "R1A3", "R1A4", "R1A5"),
#                 "x":     (     1,      2,      3,      4,      5),
#                 "y":     (     0,      0,    0.5,      0,   -0.5),
#                 "z":     (     0,      1,      0,    0.5,      0),
#             },
#             {
#                 "resname": "RES2",
#                 "beads": ("R2A1", "R2A2", "R2A3", "R2A4", "R2A5"),
#                 "x":     (     1,      2,      3,      4,      5),
#                 "y":     (    -1,      0,      0,   -0.5,      0),
#                 "z":     (     1,      2,      1,    1.5,      1),
#             },
#         ],
#         "solvcount": 1,
#         "charge": 0,
#     },
# })

############
### IONS ###
############
params = "default"
ion_defs[params] = {}
ion_defs[params]["positive"] = {
    "NA":  {"beads": "NA",  "charge": 1, "x": (0,), "y": (0,), "z": (0,)},
    "TMA": {"beads": "TMA", "charge": 1, "x": (0,), "y": (0,), "z": (0,)},
    "CA":  {"beads": "CA",  "charge": 2, "x": (0,), "y": (0,), "z": (0,)},
}

ion_defs[params]["negative"] = {
    "CL":   {"beads": "CL",  "charge": -1, "x": (0,), "y": (0,), "z": (0,)},
    "BR":   {"beads": "BR",  "charge": -1, "x": (0,), "y": (0,), "z": (0,)},
    "IOD":  {"beads": "ID",  "charge": -1, "x": (0,), "y": (0,), "z": (0,)},
    "ACE":  {"beads": "CL",  "charge": -1, "x": (0,), "y": (0,), "z": (0,)},
    "BF4":  {"beads": "BF4", "charge": -1, "x": (0,), "y": (0,), "z": (0,)},
    "PF6":  {"beads": "PF6", "charge": -1, "x": (0,), "y": (0,), "z": (0,)},
    "SCN":  {"beads": "SCN", "charge": -1, "x": (0,), "y": (0,), "z": (0,)},
    "CLO4": {"beads": "CLO", "charge": -1, "x": (0,), "y": (0,), "z": (0,)},
    "NO3":  {"beads": "NO3", "charge": -1, "x": (0,), "y": (0,), "z": (0,)},
}

###########################
### PROTEIN CHARGE DATA ###
###########################
params = "default"
prot_defs[params] = {}
prot_defs[params]["charges"] = {
    "ARG": 1, "LYS": 1, "ASP": -1, "GLU": -1,
    "GLY": 0, "ALA": 0, "ASN": 0,
    "GLN": 0, "LEU": 0, "ILE": 0,
    "VAL": 0, "SER": 0, "THR": 0,
    "CYS": 0, "PRO": 0, "HYP": 0,
    "PHE": 0, "TYR": 0, "TRP": 0,
    "ASH": 0, # Neutral version of ASP
    "MET": 0,
    "HSD": 0, "HSP": 0,
}
########################################################################################################
########################################### THE ACTUAL CLASSES #########################################
########################################################################################################

### General tools
def flatten(matrix):
    ### https://realpython.com/python-flatten-list/
    ### Fastest one shown
    flat_list = []
    for row in matrix:
        flat_list += row
    return flat_list

def center_coords(self, coords):
    '''
    Calculate lipid-spicific x/y-center
    '''
    coord_diff = (max(coords) + min(coords)) / 2
    centered_coords = [coord - coord_diff for coord in coords]
    return centered_coords

class ATOM:
    def __init__(self, bead, beadnr, x, y, z, resname, resnumber):
        self.bead = bead
        self.beadnr = beadnr
        self.x = x
        self.y = y
        self.z = z
        self.resname = resname
        self.resnr = resnumber
    
    def movex(self, x):
        self.x = x
    def movey(self, y):
        self.y = y
    def movez(self, z):
        self.z = z
        
    def move_atom(self, x, y, z):
        self.movex(x)
        self.movey(y)
        self.movez(z)
    
    def get_tuple(self):
        return (self.bead, self.beadnr, self.x, self.y, self.z, self.resname, self.resnr)
    
class RESIDUE:
    def __init__(self, resname, resnumber):
        self.resname = resname
        self.resnr = resnumber
        self.beads = []
        
    def add_bead_to_res(self, bead, beadnr, x, y, z):
        self.beads.append(ATOM(bead, beadnr, x, y, z, self.resname, self.resnr))
            
    def add_beads_to_res(self, beads = False, beadnrs = False, xs = False, ys = False, zs = False):
        assert beads and beadnrs and xs and ys and zs, "Lacking data for either 'beads', 'beadnr', 'xs', 'ys' or 'zs'"
        for bead, beadnr, x, y, z in zip(beads, beadnrs, xs, ys, zs):
            self.beads.append(ATOM(bead, beadnr, x, y, z, self.resname, self.resnr))
    
    def add_bead_data_to_res(self, bead_data):
        assert len(bead_data) == 5, "Length of list is " + str(len(bead_data)) + ". Lacking data for either 'beads', 'xs', 'ys' or 'zs'"
        for bead, beadnr, x, y, z in zip(bead_data):
            self.beads.append(ATOM(bead, beadnr, x, y, z, self.resname, self.resnr))
    
    def get_coords_res(self, AXs = "xyz"):
        AXs.lower()
        AXsList = []
        if AXs == "all":
            AXs = "xyz"
        for AX in AXs:
            if AX in ["x"]:
                AXsList.append([bead.x for bead in self.beads])
            if AX in ["y"]:
                AXsList.append([bead.y for bead in self.beads])
            if AX in ["z"]:
                AXsList.append([bead.z for bead in self.beads])
        if len(AXsList) == 1:
            AXsList = AXsList[0]
        return AXsList
    
    def get_center_point(self, centering = "ax", AXs = "xyz"):
        coords_AXs_res = self.get_coords_res(AXs)
        if len(AXs) == 1:
            coords_AXs_res = [coords_AXs_res]
        
        centers = []
        for ci, coords in enumerate(coords_AXs_res):
            if centering in ["cog", "axis", "ax"]:
                centers.append((max(coords) + min(coords)) / 2)
                
            elif centering == "mean":
                centers.append(np.mean(coords))
                
            elif centering == "beadnr":
                bead_nrs = []
                for bead_nr in target:
                    if "-" in bead_nr:
                        bead_nr1, bead_nr2 = bead_nr.split("-")
                        bead_nrs.extend(list(range(int(bead_nr1), int(bead_nr2)+1)))
                    else:
                        bead_nrs.append(int(bead_nr))
                bead_ax_vals = [coords[bead_nr] for bead_nr in bead_nrs]
                centers.append(np.mean(bead_ax_vals))
                
            elif centering == "vals":
                centers.append(target[ci])
        return tuple(centers)

class MOLECULE:
    def __init__(self, charge = False, molname = False):
        self.residues = []
        self.center = False
        self.n_residues = 0
        self.resnames = []
        if charge is False:
            self.charge = 0
            self.charge_has_been_set = False
        else:
            self.charge = charge
            self.charge_has_been_set = True
        self.molname = molname
        self.last_res_n = 0
    
    def add_res(self, resname, resnumber = False):
        if not resnumber:
            resnumber = self.n_residues
#         self.residues[resnumber] = RESIDUE(resname)
        self.residues.append(RESIDUE(resname, resnumber))
        self.resnames.append(resname)
        self.last_res_n = self.n_residues
        self.n_residues += 1
    
    def add_bead(self, bead, beadnr, x, y, z, resnumber = False, addcharge = 0):
        if not resnumber:
            resnumber = self.last_res_n
        self.residues[resnumber].add_bead_to_res(bead, beadnr, x, y, z)
        if addcharge != 0:
            self.charge_add(addcharge)
    
    def add_res_and_beads(self, resname, beads = False, beadnrs = False, xs = False, ys = False, zs = False, bead_data = False, resnumber = False):
        if not resnumber:
            resnumber = self.n_residues
        self.residues.append(RESIDUE(resname, resnumber))
        self.resnames.append(resname)
        self.last_res_n = self.n_residues
        self.n_residues += 1
        if bead_data:
            self.residues[self.last_res_n].add_bead_data_to_res(bead_data)
        else:
            self.residues[self.last_res_n].add_beads_to_res(beads, beadnrs, xs, ys, zs)
    
    def charge_add(self, charge):
        self.charge += charge
        self.charge_has_been_set = True
    
    def charge_set(self, charge):
        self.charge = charge
        self.charge_has_been_set = True
    
    def get_coords(self, AXs = "xyz"):
        AXs.lower()
        AXsList = []
        if AXs == "all":
            AXs = "xyz"
        for AX in AXs:
            AXsList.append(flatten([res.get_coords_res(AX) for res in self.residues]))
#         if len(AXsList) == 1:
#             AXsList = AXsList[0]
        return AXsList
    
    def get_beads(self, AXs = "all"):
        AXsList = self.get_coords(AXs)
        assert len(AXsList) > 1, "Length of AXsList must be greater than 1 to create beads"
        return list(zip(*AXsList))
    
    def get_centered_coords(self, centering = "ax", AXs = "xyz", target = False):
        coords_AXs = self.get_coords(AXs)
        centered_coords = []
        
        centers = self.get_center_point(centering = centering, AXs = AXs, target = target)
        
        for coords, center in zip(coords_AXs, centers):
            centered_coords.append([coord - center for coord in coords])

        return tuple(centered_coords)
    
    def get_center_point(self, centering = "ax", AXs = "xyz", target = False):
        coords_AXs = self.get_coords(AXs)
        
        centers = []
        for ci, coords in enumerate(coords_AXs):
            if centering in ["cog", "axis", "ax"]:
                centers.append((max(coords) + min(coords)) / 2)
                
            elif centering == "mean":
                centers.append(np.mean(coords))
                
            elif centering == "beadnr":
                bead_nrs = []
                for bead_nr in target:
                    if "-" in bead_nr:
                        bead_nr1, bead_nr2 = bead_nr.split("-")
                        bead_nrs.extend(list(range(int(bead_nr1), int(bead_nr2)+1)))
                    else:
                        bead_nrs.append(int(bead_nr))
                bead_ax_vals = [coords[bead_nr] for bead_nr in bead_nrs]
                centers.append(np.mean(bead_ax_vals))
                
            elif centering == "resnr":
                res_nrs = []
                ### Find all residues
                for res_nr in target:
                    if "-" in res_nr:
                        res_nr1, res_nr2 = res_nr.split("-")
                        res_nrs.extend(list(range(int(res_nr1), int(res_nr2)+1)))
                    else:
                        res_nrs.append(int(res_nr))
                ### Find center of each individual residue
                res_centers = [
                    self.residues[res_nr].get_center_point(centering = "mean", AXs = AXs[ci])
                    for res_nr in res_nrs
                ]
                ### Find center of residues
                centers.append(np.mean(res_centers))
                
            elif centering == "vals":
                centers.append(target[ci])
        return tuple(centers)
    
    def set_coords_to_center(self, centering = "ax", target = False):
        xsc, ysc, zsc = self.get_centered_coords(centering, "xyz", target)
        i = 0
        for ri, res in enumerate(self.residues):
            for bi, bead in enumerate(res.beads):
                self.residues[ri].beads[bi].move_atom(xsc[i], ysc[i], zsc[i])
                i += 1
    
    def move_coords(self, translation = [0, 0, 0]):
        tx, ty, tz = translation
        i = 0
        for ri, res in enumerate(self.residues):
            for bi, bead in enumerate(res.beads):
                self.residues[ri].beads[bi].move_atom(bead.x+tx, bead.y+ty, bead.z+tz)
                i += 1
    
    def move_atom(self, ri, bi, translation = [0, 0, 0]):
        tx, ty, tz = translation
        bead = self.residues[ri].beads[i]
        self.residues[ri].beads[bi].move_atom(bead.x+tx, bead.y+ty, bead.z+tz)
    
    def rotate_coords(self, rotation = [0, 0, 0]):
        x_deg, y_deg, z_deg = rotation
        i = 0
        for ri, res in enumerate(self.residues):
            for bi, bead in enumerate(res.beads):
                
                ### Get positions to limit calls to bead
                x = bead.x
                y = bead.y
                z = bead.z
                
                ### First convert angles to radians
                x_rad = math.radians(x_deg)
                y_rad = math.radians(y_deg)
                z_rad = math.radians(z_deg)

                ### Calculate rotation matrices (https://en.wikipedia.org/wiki/Rotation_matrix)
                rm_x = [
                    [1, 0, 0],
                    [0, math.cos(x_rad), -math.sin(x_rad)],
                    [0, math.sin(x_rad), math.cos(x_rad)]
                ]
                rm_y = [
                    [math.cos(y_rad), 0, math.sin(y_rad)],
                    [0, 1, 0],
                    [-math.sin(y_rad)  , 0, math.cos(y_rad)]
                ]
                rm_z = [
                    [math.cos(z_rad), -math.sin(z_rad), 0],
                    [math.sin(z_rad) , math.cos(z_rad), 0],
                    [0, 0, 1]
                ]

                ### Combine rotation matrices
                rm_xy  = [
                    [sum(r * c for r, c in zip(row, col)) for col in zip(*rm_y)]
                    for row in rm_x
                ]
                rm_xyz = [
                    [sum(r * c for r, c in zip(row, col)) for col in zip(*rm_z)]
                    for row in rm_xy
                ]

                ### New positions
                new_x = rm_xyz[0][0] * x + rm_xyz[0][1] * y + rm_xyz[0][2] * z
                new_y = rm_xyz[1][0] * x + rm_xyz[1][1] * y + rm_xyz[1][2] * z
                new_z = rm_xyz[2][0] * x + rm_xyz[2][1] * y + rm_xyz[2][2] * z
                
                self.residues[ri].beads[bi].move_atom(new_x, new_y, new_z)
                i += 1
    
    def get_max_length(self):
        beads = self.get_beads()
        max_length = max([math.dist(bead1, bead2) for bead1 in beads for bead2 in beads])
        return max_length

    def get_res_beads_info(self, output_type="ATOM"):
        res_beads_info = []
        for res in self.residues:
            for bead in res.beads:
                if output_type == "ATOM":
                    res_beads_info.append(bead)
                elif output_type in ["tuple", "ziptuple"]:
                    res_beads_info.append(bead.get_tuple())
        if output_type == ["ziptuple"]:
            res_beads_info = zip(*res_beads_info)
        return res_beads_info
    
class LIPID(MOLECULE):
    def __init__(self, charge = 0, hydr_z = 0, molname = False):
        super().__init__(charge = charge, molname = molname)
        self.hydr_z = hydr_z ### Hydrophobic z-height delimiter
        self.ratio = 0
    
    def get_radius(self):
        beads  = self.get_beads("xy")
        center = self.get_center_point("ax", "xy")
        radius = max([math.dist(bead, center) for bead in beads])
        return radius
    
    def set_xy_to_center(self, centering = "ax"):
        xsc, ysc = self.get_centered_coords(centering, "xy")
        i = 0
#         new_ress = []
        for ri, res in enumerate(self.residues):
#             new_ress.append(RESIDUE(res.resname, res.resnr))
            for bi, bead in enumerate(res.beads):
                self.residues[ri].beads[bi].move_atom(xsc[i], ysc[i], bead.z)
#                 new_ress[ri].add_bead_to_res(bead.bead, bead.beadnr, xsc[i], ysc[i], bead.z)
                i += 1
    
    def ratio_add(self, ratio):
        self.ratio = ratio

class SOLVENT(MOLECULE):
    def __init__(self, charge = 0, solvcount = 1, molname = False):
        super().__init__(charge = charge, molname = molname)
        self.solvcount = solvcount ### AA-to-CG convertion
        self.molarity = False
        self.density = False
        self.molar_mass = False
        
    def solvcount_set(self, solvcount):
        self.solvcount = solvcount

    def molarity_set(self, molarity):
        self.molarity = molarity
    
    def density_set(self, density):
        self.density = density
    
    def molar_mass_set(self, molar_mass):
        self.molar_mass = molar_mass
    
    def count_set(self, count):
        self.count = count
    
class PROTEIN(MOLECULE):
    def __init__(self, charge = 0, molname = False):
        super().__init__(charge = charge, molname = molname)

class CGSB:
    
    def __init__(self, **kwargs):
        try:
            self.lipid_defs = lipid_defs.copy()
        except:
            print_term("WARNING: No lipid definitions found 'lipid_defs'", warn=True)
            self.lipid_defs = {}
        try:
            self.solvent_defs = solvent_defs.copy()
        except:
            print_term("WARNING: No solvent definitions found 'solvent_defs'", warn=True)
            self.solvent_defs = {}
        try:
            self.ion_defs = ion_defs.copy()
        except:
            print_term("WARNING: No ions definitions found 'ion_defs'", warn=True)
            self.ion_defs = {}
        try:
            self.prot_defs = prot_defs.copy()
        except:
            print_term("WARNING: No protein definitions found 'prot_defs'", warn=True)
            self.prot_defs = {}
        
        self.lipid_dict = {}
        self.solvent_dict = {}
        self.ion_dict = {}
        self.prot_dict = {}
        
        self.PROTEINS = {}
        self.PROTEINS_cmds = []
        
        self.MEMBRANES = {}
        self.MEMBRANES_cmds = []
        
        self.SOLVATIONS = {}
        self.SOLVATIONS_cmds = []
        ### Floodings are just solvations with some different default settings
        self.FLOODINGS_cmds = []
        
        self.itp_moltypes = {}
        self.ITP_INPUT_cmds = []
        
        self.SOLUTE_INPUT_cmds = []
        
        self.PLOT_cmd = []
        self.plot_data = {}
        self.plots_requested = False
        
        self.plot_grid = False
        
        self.PICKLE_cmd = False
        
        self.randseed = round(time.time())
        
        self.sys_params = "default"
        self.prot_params = False # Only used if specifically set
        self.lipid_params = False # Only used if specifically set
        self.solv_params = False # Only used if specifically set
        
        self.system_charge = 0
        self.system_name = "PLACEHOLDER_TITLE"
        
#         self.output_system_file_name     = "output"
        self.output_system_pdb_file_name = False
        self.output_system_gro_file_name = False
        self.output_topol_file_name      = "topol.top"
        
        self.LOG_FILE = []
        self.output_log_file_name = False
        
        self.pbc_set  = []
        self.pbc_type = "rectangular"
        self.backup   = True
        self.pickle   = False
        
        self.debug_prints = False
        self.extra_info = True
        self.subleaflet_extra_info = True
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
                if type(cmd) not in [list, tuple]:
                    cmd = [cmd]
                for subcmd in cmd:
                    self.PROTEINS_cmds.extend([subcmd])
                
            if any(key.startswith(i) for i in ["membrane", "memb"]):
                if type(cmd) not in [list, tuple]:
                    cmd = [cmd]
                for subcmd in cmd:
                    self.MEMBRANES_cmds.extend([subcmd])
                
            if any(key.startswith(i) for i in ["solvation", "solv"]):
                if type(cmd) not in [list, tuple]:
                    cmd = [cmd]
                for subcmd in cmd:
                    self.SOLVATIONS_cmds.extend([subcmd])
                
            if any(key.startswith(i) for i in ["flooding", "flood"]):
                if type(cmd) not in [list, tuple]:
                    cmd = [cmd]
                for subcmd in cmd:
                    if "flooding:" in subcmd:
                        subcmd_split = subcmd.split()
                        subcmd = " ".join([string for string in subcmd_split if not string.startswith("flooding:")])
                    subcmd = " ".join(["flooding:True", subcmd, "count:True", "solv_molarity:1", "salt_molarity:1"])
                    self.FLOODINGS_cmds.extend([subcmd])
            
            ### Box type
            if key in ["pbc_type", "box_type"]:
                if cmd in ["rect", "rectangular"]:
                    self.pbc_type = "rectangular"
                elif cmd in ["hex", "hexagonal"]:
                    ### Hexagonal prism (hexagonal XY-plane, and straight Z-axis)
                    ### Y-values ignored as only one side-length is allowed. X-value used instead.
                    ### Z-value used for height
                    ### Not actually a hexagon, but instead a parallelepiped constituting a third of the hexagon
                    self.pbc_type = "hexagonal"
                elif cmd in ["optimal", "dodecahedron"]:
                    ### Rhombic dodecahedron (hexagonal XY-plane, Z-angled)
                    ### Y-values ignored as only one side-length is allowed. X-value used instead.
                    ### Z-value unused as it is calculated from x-value due to angles.
                    ### Not actually a rhombic dodecahedron, but instead a parallelepiped constituting a portion of dodecahedron
                    self.pbc_type = "dodecahedron"
                else:
                    assert False, "Incorrect pbc type given: " + str(key, cmd)
            
            ### Box size
            if key in ["pbc", "box"]:
                ### Evaluated after loop due to dependence on box type
                momentary_pbc = cmd
            
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
            if key in ["out_sys", "out", "o"]:
                if not any([cmd.endswith(i) for i in ["pdb", "gro"]]):
                    self.output_system_pdb_file_name = cmd + ".pdb"
                    self.output_system_gro_file_name = cmd + ".gro"
                elif cmd.endswith("pdb"):
                    self.output_system_pdb_file_name = cmd
                elif cmd.endswith("gro"):
                    self.output_system_gro_file_name = cmd
                else:
                    assert False, "Unknown file extension used for 'output_system': " + pdb
                    
            if key in ["out_sys_pdb", "out_pdb", "o_pdb"]:
                if cmd.endswith("pdb"):
                    self.output_system_pdb_file_name = cmd + ".pdb"
                else:
                    self.output_system_pdb_file_name = cmd + ".pdb"
                    
            if key in ["out_sys_gro", "out_gro", "o_gro"]:
                if cmd.endswith("gro"):
                    self.output_system_gro_file_name = cmd + ".gro"
                else:
                    self.output_system_gro_file_name = cmd + ".gro"
            
            if key in ["output_topol", "out_top", "top_out", "top_output"]:
                self.output_topol_file_name = cmd
                
            if key in ["output_log", "out_log", "log_out", "log"]:
                self.output_log_file_name = cmd
                
#             if key in ["imp_o", "output_imported"]:
#                 self.output_imported = cmd
            
            ### Extra functionalities
#             if key in ["plot"]:
#                 self.PLOT_cmd = cmd
#                 if self.PLOT_cmd:
#                     self.plots_requested = True
                    
            if key in ["plot_grid"]:
                self.plot_grid = True
                    
            if key in ["rand", "randseed"]:
                if not cmd is False:
                    if type(cmd) == int:
                        self.randseed = cmd
                    else:
                        isnumber, isint = self.is_number(cmd)
                        if isnumber:
                            if isint:
                                self.randseed = cmd
                            else:
                                self.randseed = round(cmd)
                    
            if key in ["pickle"]:
                self.PICKLE_cmd = cmd
                
            if key in ["backup"]:
                self.backup = cmd

            if key in ["params", "sys_params"]:
                self.sys_params = cmd
                
            if key in ["prot_params"]:
                self.prot_params = cmd
                
            if key in ["lipid_params"]:
                self.lipid_params = cmd
                
            if key in ["solv_params"]:
                self.solv_params = cmd
            
            ### Printer settings
            if key in ["quiet"]:
                self.quiet = cmd
                
            if key in ["debug"]:
                self.debug_prints = cmd
                
            if key in ["extra"]:
                self.extra_info = cmd
            
            if key in ["sl_extra", "subleaflet_extra"]:
                self.subleaflet_extra_info = cmd
                
            if key in ["warn"]:
                self.warnings = cmd
                
            ### Run the program
            if key in ["run"]:
                if type(cmd) == str:
                    cmd = ast.literal_eval(cmd)
                self.RUN = cmd
        
        if len(self.FLOODINGS_cmds) > 0:
            self.SOLVATIONS_cmds = self.FLOODINGS_cmds + self.SOLVATIONS_cmds
        
        ### PBC type settings:
        cmd = momentary_pbc
        if self.pbc_type == "rectangular":
            if len(cmd) == 2:
                self.pbcx, self.pbcy = cmd[0]*10
                self.pbcz = cmd[1]*10
                self.pbc_box = [self.pbcx, self.pbcy, self.pbcz]
            elif len(cmd) == 3:
                self.pbcx, self.pbcy, self.pbcz = [i*10 for i in cmd]
            else:
                assert False, "Incorrect pbc box dimensions for pbc type '"+self.pbc_type+"': " + str(cmd)
            self.gro_box_vectors = [
                float(self.pbcx/10),     # vector: vax or v1x
                float(self.pbcy/10),     # vector: vby or v2y
                float(self.pbcz/10),     # vector: vcz or v3z
                float(0),                # vector: vay or v1y
                float(0),                # vector: vaz or v1z
                float((self.pbcx/10)/2), # vector: vbx or v2x
                float(0),                # vector: vbz or v2z
                float(0),                # vector: vcx or v3x
                float(0),                # vector: vcy or v3y
            ]
            self.pdb_box_dimension = [
                float(self.pbcx), # Axis length:  x
                float(self.pbcy), # Axis length:  y
                float(self.pbcz), # Axis length:  z
                float(90),        # Corner angle: alpha
                float(90),        # Corner angle: beta
                float(90),        # Corner angle: gamma
            ]

        ### Following math taken from insane.py and https://manual.gromacs.org/current/reference-manual/algorithms/periodic-boundary-conditions.html
        elif self.pbc_type == "hexagonal":
            ### Not actually a hexagon, but instead a parallelepiped constituting a third of it
            if len(cmd) == 2:
                self.pbcx = cmd[0]*10
                self.pbcy = math.sqrt(3)*self.pbcx/2
                self.pbcz = cmd[1]*10
            else:
                assert False, "Incorrect pbc box dimensions for pbc type '"+self.pbc_type+"': " + str(cmd)
            self.gro_box_vectors = [
                float(self.pbcx/10),     # vector: vax or v1x
                float(self.pbcy/10),     # vector: vby or v2y
                float(self.pbcz/10),     # vector: vcz or v3z
                float(0),                # vector: vay or v1y
                float(0),                # vector: vaz or v1z
                float((self.pbcx/10)/2), # vector: vbx or v2x
                float(0),                # vector: vbz or v2z
                float(0),                # vector: vcx or v3x
                float(0),                # vector: vcy or v3y
            ]
            self.pdb_box_dimension = [
                float(self.pbcx), # Axis length:  x
                float(self.pbcx), # Axis length:  y
                float(self.pbcz), # Axis length:  z
                float(90),        # Corner angle: alpha
                float(90),        # Corner angle: beta
                float(60),        # Corner angle: gamma
            ]
            
        elif self.pbc_type == "dodecahedron":
            ### Not actually a dodecahedron, but instead a parallelepiped constituting a third of it
            if len(cmd) == 1:
                self.pbcx = cmd[0]
                self.pbcy = math.sqrt(3)*self.pbcx/2
                self.pbcz = math.sqrt(6)*self.pbcx/3
            else:
                assert False, "Incorrect pbc box dimensions for pbc type '"+self.pbc_type+"': " + str(cmd)
            self.gro_box_vectors = [
                float(self.pbcx/10),                  # vector: vax or v1x
                float(self.pbcy/10),                  # vector: vby or v2y
                float(self.pbcz/10),                  # vector: vcz or v3z
                float(0),                             # vector: vay or v1y
                float(0),                             # vector: vaz or v1z
                float((self.pbcx/10)/2),              # vector: vbx or v2x
                float(0),                             # vector: vbz or v2z
                float((self.pbcx/10)/2),              # vector: vcx or v3x
                float(math.sqrt(3)*(self.pbcx/10)/6), # vector: vcy or v3y
            ]
            self.pdb_box_dimension = [
                float(self.pbcx), # Axis length:  x
                float(self.pbcx), # Axis length:  y
                float(self.pbcx), # Axis length:  z
                float(60),        # Corner angle: alpha
                float(60),        # Corner angle: beta
                float(60),        # Corner angle: gamma
            ]
            
        self.pbc_box = [self.pbcx, self.pbcy, self.pbcz]
        
        ### Setting randseed
        self.print_term("\n" + "Setting random seed to:", self.randseed)
        random.seed(self.randseed)
        np.random.seed(self.randseed)
        
        if not self.output_system_pdb_file_name and not self.output_system_gro_file_name:
            self.output_system_pdb_file_name = "output.pdb"
            self.output_system_gro_file_name = "output.gro"
        
        if self.RUN:
            self.run()
        
    def run(self):
        '''
        Runs the entire system creation process
        '''
        ### Initial checks
#         assert len(self.pbc_box) > 0, "Box dimensions not set. Please do so using 'box=[x,y,z]'"
#         assert len(self.pbc_box) == 3, "Box dimensions improperly defined. 3 dimensions must be given"
#         assert all([self.is_number(ax)[0] for ax in self.pbc_box]), "Not all box values are numbers"
        
        assert any([len(cmd) > 0 for cmd in [self.PROTEINS_cmds, self.MEMBRANES_cmds, self.SOLVATIONS_cmds]]), (
            "Running requires at least one command of one of the following types: 'protein', 'membrane' or 'solvation'"
        )
        
        ### Topology
        self.itp_read_initiater()
        
        self.print_term("------------------------------ PREPROCESSING", spaces=0)
        preprocessing_tic = time.time()
        ### Definition preprocessing
        self.import_structures_handler()
        self.lipid_defs_preprocessor()
        self.solvent_defs_preprocessor()
        self.ion_defs_preprocessor()
        
        ### Command preprocessing
        self.prot_preprocessor()
        self.memb_preprocessor()
        self.solv_preprocessor()
        preprocessing_toc = time.time()
        preprocessing_time = round(preprocessing_toc - preprocessing_tic, 4)
        self.print_term("------------------------------ PREPROCESSING COMPLETE", "(time spent: "+str(preprocessing_time)+" [s])", "\n", spaces=0)
        
        ### Run the program
        self.prot_placer()
        self.subleaflet_poly_maker()
        self.holed_subleaflet_bbox_maker()
        self.lipid_calculator()
        self.planar_grid_maker()
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
    def print_term(self, *string, spaces=0, space="    ", sep=" ", end="\n", debug = False, extra = False, warn = False):
        '''
        Specialized printing function.
        Allows for easy customization of requirements for when a print statement should be printed
        '''
        print_true = False
        if type(string) in [tuple, list]:
            string = space*spaces + sep.join([str(i) for i in string])

        ### Check if string should be printed based on settings
        if debug:
            if self.debug_prints:
                print_true = True
                if not string.startswith("DEBUG:"):
                    string = "DEBUG: " + string
        elif extra:
            if self.extra_info:
                print_true = True
        elif warn:
            if self.warnings:
                print_true = True
                if not string.startswith("WARNING:"):
                    string = "WARNING: " + string
        else:
            ### If no special case, then assume it should be printed
            print_true = True

        ### Print to terminal if not quiet
        if print_true and not self.quiet:
            print(string, end=end)
        
        ### Appends string to log file, for later writing
        if print_true and self.output_log_file_name:
            self.LOG_FILE.append(string + end)
    
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
                if line_nr == 0:
                    continue
                elif line_nr == 1:
                    tot_atoms = int(ast.literal_eval(line.strip()))
                    continue
                if tot_atoms == 0:
                    break
                else:
                    if len(line):
                        tot_atoms -= 1
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
            "dihedraltypes": ["funct", "angle", "fc", "a0", "a1", "a2", "a3", "a4"],
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
                    try:
                        self.itp_defs[cur_cmd][string[1]] = {}
                        for i, st in enumerate(string[2:]):
                            self.itp_defs[cur_cmd][string[1]][defs_dict_names[cur_cmd][i]] = st
                    except:
                        continue

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
        PCL = [6, 5, 1 + 4, 1, 3 + 1, 1, 4, 1 + 3, 8, 8, 8, 6, 6, 10, 2, 2] 
        if len(r_name) > 4:
            r_name = r_name[:4]
        if len(r_name) < 4:
            r_name = r_name + " "
        if len(a_name) > 5:
            a_name = a_name[:5]
        if len(a_name) < 4:
            a_name = " " + a_name
        string = '{ATOM:<{L0}}{a_nr:>{L1}}{a_name:<{L2}}{aLoc:>{L3}}{r_name:>{L4}}{chain:^{L5}}{r_nr:>{L6}}{aChar:>{L7}}{x:>{L8}.3f}{y:>{L9}.3f}{z:>{L10}.3f}{oc:>{L11}.2f}{T:>{L12}.2f}{JUMP:>{L13}}{E:>{L14}}{C:>{L15}}'.format(
            ATOM = ATOM, L0 = PCL[0],
            a_nr = a_nr, L1 = PCL[1], # int
            a_name = a_name, L2 = PCL[2],
            aLoc = aLoc, L3 = PCL[3],
#             aLoc = "", L3 = 0,
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
        if len(a_name) > 5:
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
        
        out_stringlist = []
        files_string = ""
        if self.output_system_pdb_file_name:
            out_stringlist.append("PDB")
        if self.output_system_gro_file_name:
            out_stringlist.append("GRO")
        if len(out_stringlist) > 1:
            files_string = "files"
        else:
            files_string = "file"
        self.print_term("Writing structure", files_string, "(" + "/".join(out_stringlist) + ")")
            
        output_system_pdb_file_lines = []
        output_system_gro_file_lines = []
        
        ###############################
        ### Beginning of file lines ###
        ###############################
        if self.output_system_pdb_file_name:
            a, b, c, alpha, beta, gamma = self.pdb_box_dimension
            output_system_pdb_file_lines = [
                "TITLE     " + self.system_name,
                "REMARK    " + "PLACEHOLDER_REMARK",
                '{Rname:<{RnameL}}{a:>{aL}.3f}{b:>{bL}.3f}{c:>{cL}.3f}{alpha:>{alphaL}.2f}{beta:>{betaL}.2f}{gamma:>{gammaL}.2f} {sGroup:<{sGroupL}}{z:>{zL}}'.format(
                    Rname = "CRYST1", RnameL = 6,
                    a = a,            aL = 9,
                    b = b,            bL = 9,
                    c = c,            cL = 9,
                    alpha = alpha,    alphaL = 7,
                    beta = beta,      betaL = 7,
                    gamma = gamma,    gammaL = 7,
                    sGroup = "P 1",   sGroupL = 11,
                    z = 1,            zL = 4,
                ),
                "MODEL        1",
        #         "".join([str(i) for i in range(1, 10)]) + "".join([str(i) for _ in range(8) for i in range(0, 10)]),
            ]
        if self.output_system_gro_file_name:
            
            output_system_gro_file_lines = [
                self.system_name,
                "PLACEHOLDER_ATOM_COUNT",
        #         "".join([str(i) for i in range(1, 10)]) + "".join([str(i) for _ in range(8) for i in range(0, 10)]),
            ]

        atom_count = 0
        self.molecules_for_top = []

        old_res_nr = 0
        atom_nr    = 0
        res_nr     = 0
        
        #####################
        ### Protein lines ###
        #####################
        if len(self.PROTEINS) != 0:
            for protein_nr, protein in self.PROTEINS.items():
                for mol_name in protein["mol_names"]:
                    self.molecules_for_top.append((mol_name, str(1)))
                current_prot_res = 0
                
                original_bead_info = protein["beads"].keys()
                bead_vals = [bead.get_tuple() for bead in protein["protein"].get_res_beads_info()]
                for (i, atom, res), (a_name, beadnr, x, y, z, r_name, resnr) in zip(original_bead_info, bead_vals):
#                 for (i, atom, res), bead_vals in protein["beads_centered"].items():
                    if current_prot_res != res:
                        current_prot_res = res
                        res_nr += 1
                        if res_nr >= 10000:
                            res_nr -= 10000 * (res_nr // 10000)
                    atom_nr += 1
                    atom_count += 1
                    if atom_nr >= 100000:
                        atom_nr -= 100000 * (atom_nr // 100000)
                    x += self.pbc_box[0] / 2
                    y += self.pbc_box[1] / 2
                    z += self.pbc_box[2] / 2
                    if self.output_system_pdb_file_name:
                        output_system_pdb_file_lines.append(self.pdb_atom_writer("ATOM", atom_nr, a_name, " ", r_name, "A", res_nr, " ", float(x), float(y), float(z), float(1), float(0), " ", " ", " "))
                    if self.output_system_gro_file_name: ### gro coordinates are in [nm] not [Å]
                        output_system_gro_file_lines.append(self.gro_atom_writer(res_nr, r_name, a_name, atom_nr, x / 10, y / 10, z / 10, " ", " ", " "))

        ######################
        ### Membrane lines ###
        ######################
        if len(self.MEMBRANES) != 0:
            for memb_key, memb_dict in self.MEMBRANES.items():
                for leaflet_key, leaflet in memb_dict["leaflets"].items():
                    self.molecules_for_top.extend(sorted(leaflet["leaf_lipid_count"], key=lambda m: m[0]))
                    
                    ### Ordering lipids per leaflet for concise topology creation
#                     ordered_lipids = sorted(
#                         [(grid_point, grid_vals) for (grid_point, grid_vals) in leaflet["grid_adjusted"].items()],
#                         key=lambda gp: gp[1]["lipid"]["name"]
#                     )
                    ordered_lipids = sorted(
                        [grid_vals for grid_vals in leaflet["grid_lipids"]],
                        key=lambda gp: gp["lipid"]["name"]
                    )
            
                    for i, grid_vals in enumerate(ordered_lipids):
                        res_nr += 1
                        if res_nr >= 10000:
                            res_nr -= 10000 * (res_nr // 10000)
                        r_name = grid_vals["lipid"]["name"]

                        for x, y, z, bead_name in zip(grid_vals["lipid"]["x"], grid_vals["lipid"]["y"], grid_vals["lipid"]["z"], grid_vals["lipid"]["beads"]):
                            atom_nr += 1
                            atom_count += 1
                            if atom_nr >= 100000:
                                atom_nr -= 100000 * (atom_nr // 100000)
                            x += self.pbc_box[0] / 2
                            y += self.pbc_box[1] / 2
                            z += self.pbc_box[2] / 2
                            a_name = bead_name
                            string = ""
                            if self.output_system_pdb_file_name:
                                output_system_pdb_file_lines.append(self.pdb_atom_writer("ATOM", atom_nr, a_name, " ", r_name, "A", res_nr, " ", float(x), float(y), float(z), float(1), float(0), " ", " ", " "))
                            if self.output_system_gro_file_name: ### gro coordinates are in [nm] not [Å]
                                output_system_gro_file_lines.append(self.gro_atom_writer(res_nr, r_name, a_name, atom_nr, x / 10, y / 10, z / 10, " ", " ", " "))

        #########################
        ### Solvent/ion lines ###
        #########################
        if len(self.SOLVATIONS) != 0:
            for solvation_nr, solvation in self.SOLVATIONS.items():
                solvent_count = [(key_name, count) for (key_name, key_type, key_charge), count in solvation["solv_count"].items()]
                self.molecules_for_top.extend(solvent_count)

                for i, grid_point_3D in enumerate(solvation["grid"]):
                    res_nr += 1
                    if res_nr >= 10000:
                        res_nr -= 10000 * (res_nr // 10000)

                    for (x, y, z), r_name, bead_name in zip(grid_point_3D["coords"], grid_point_3D["resnames"], grid_point_3D["beads"]):
                        atom_nr += 1
                        atom_count += 1
                        if atom_nr >= 100000:
                            atom_nr -= 100000 * (atom_nr // 100000)
                        x += self.pbc_box[0] / 2
                        y += self.pbc_box[1] / 2
                        z += self.pbc_box[2] / 2
                        a_name = bead_name
                        string = ""
                        if self.output_system_pdb_file_name:
                            output_system_pdb_file_lines.append(self.pdb_atom_writer("ATOM", atom_nr, a_name, " ", r_name, "A", res_nr, " ", float(x), float(y), float(z), float(1), float(0), " ", " ", " "))
                        if self.output_system_gro_file_name: ### gro coordinates are in [nm] not [Å]
                            output_system_gro_file_lines.append(self.gro_atom_writer(res_nr, r_name, a_name, atom_nr, x / 10, y / 10, z / 10, " ", " ", " "))
        
        #########################
        ### End of file lines ###
        #########################
        if self.output_system_pdb_file_name:
            output_system_pdb_file_lines.append("TER")
            output_system_pdb_file_lines.append("END")
        
        if self.output_system_gro_file_name:
            v1x, v2y, v3z, v1y, v1z, v2x, v2z, v3x, v3y = self.gro_box_vectors
            output_system_gro_file_lines[1] = " " + str(atom_count)
            output_system_gro_file_lines.append( ### gro vectors are in [nm] not [Å]
                '{v1x:>{v1xL}.5f}{v2y:>{v2yL}.5f}{v3z:>{v3zL}.5f}{v1y:>{v1yL}.5f}{v1z:>{v1zL}.5f}{v2x:>{v2xL}.5f}{v2z:>{v2zL}.5f}{v3x:>{v3xL}.5f}{v3y:>{v3yL}.5f}'.format(
                    v1x = v1x, v1xL = 10,
                    v2y = v2y, v2yL = 10,
                    v3z = v3z, v3zL = 10,
                    v1y = v1y, v1yL = 10,
                    v1z = v1z, v1zL = 10,
                    v2x = v2x, v2xL = 10,
                    v2z = v2z, v2zL = 10,
                    v3x = v3x, v3xL = 10,
                    v3y = v3y, v3yL = 10,
                )
            )
            output_system_gro_file_lines.append("")
        
        if self.backup:
            if self.output_system_pdb_file_name:
                self.backupper(self.output_system_pdb_file_name)
            if self.output_system_gro_file_name:
                self.backupper(self.output_system_gro_file_name)

        if self.output_system_pdb_file_name:
            new_file = open(self.output_system_pdb_file_name, "w")
            for line in output_system_pdb_file_lines:
                new_file.write(line + "\n")
            new_file.close()
            self.print_term("---- PDB file written:", self.output_system_pdb_file_name)
        if self.output_system_gro_file_name:
            new_file = open(self.output_system_gro_file_name, "w")
            for line in output_system_gro_file_lines:
                new_file.write(line + "\n")
            new_file.close()
            self.print_term("---- GRO file written:", self.output_system_gro_file_name)
        self.print_term("")

    def topol_file_writer(self):
        if self.output_topol_file_name:
            self.print_term("Writing topology file:", self.output_topol_file_name)
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
            self.print_term("---- Topology file written", "\n")
    
    def log_file_writer(self):
        if self.output_log_file_name:
            self.print_term("Writing log file file:", self.output_log_file_name)
            if self.backup:
                self.backupper(self.output_log_file_name)
            new_file = open(self.output_log_file_name, "w")
            for line in self.LOG_FILE:
                new_file.write(line)
            new_file.close()
            self.print_term("---- Log file written", "\n")
    
    ###################################################
    ### Lipid/solvent/ion definitions preprocessors ###
    ###################################################
    def lipid_defs_preprocessor(self):
        '''
        Preprocesses lipid defintions, by rearranging the data into a class
        '''
        self.print_term("Preprocessing lipid definitions")
        ### Loop over generel lipid types (e.g. phospholipids, sterols, etc.)
        for (lipid_type, params), lipid_type_dict in self.lipid_defs.items():
            if params not in self.lipid_dict.keys():
                self.lipid_dict[params] = {}
            x, y, z = lipid_type_dict["x"], lipid_type_dict["y"], lipid_type_dict["z"]

            ### Centering x/y/z-values for current general lipid type
            if "center" in lipid_type_dict.keys():
                hydr_z = lipid_type_dict["center"]
            else:
                hydr_z = 0

            ### bead distance adjustment values for current general lipid type
            if "bd" in lipid_type_dict.keys():
                bdx, bdy, bdz = lipid_type_dict["bd"]
            else:
                bdx, bdy, bdz = 0.25, 0.25, 0.3
            
            if "charges" in lipid_type_dict.keys():
                lipid_type_charges = lipid_type_dict["charges"]
            else:
                lipid_type_charges = 0
            
            lipid_names = []
            
            for (lipid_name, data_specifier), lipid_details in lipid_type_dict["lipids"].items():
                lipid_names.append(lipid_name)
                if lipid_name not in self.lipid_dict[params].keys():
                    self.lipid_dict[params][lipid_name] = LIPID(hydr_z = hydr_z, molname = lipid_name)
                if data_specifier == "beads":
                    beads = list(filter(None, lipid_details.split(" ")))
                    ### Remove coordinate and bead indexes if no bead is assigned
                    ### x and z are centered using "set_xy_to_center()" method
                    ### z centering according to hydrophobic height "hydr_z"
                    lx = [xi * bdx * 10 for xi, bead in zip(x, beads) if bead != "-"]
                    ly = [yi * bdy * 10 for yi, bead in zip(y, beads) if bead != "-"]
                    lz = [(zi - hydr_z) * bdz * 10 for zi, bead in zip(z, beads) if bead != "-"]
                    beads = [bead for bead in beads if bead != "-"]
                    beadnrs = list([i for i in range(len(beads))])
                    
                    if not self.lipid_dict[params][lipid_name].charge_has_been_set:
                        self.lipid_dict[params][lipid_name].charge_set(0)
                        if lipid_type_charges != 0:
                            charge_beads, charge_vals = zip(*lipid_type_charges)
                            for bead in beads:
                                if bead in charge_beads:
                                    charge_index = charge_beads.index(bead)
                                    self.lipid_dict[params][lipid_name].charge_add(charge_vals[charge_index])
                    
                    self.lipid_dict[params][lipid_name].add_res_and_beads(
                        lipid_name,
                        beads = beads,
                        beadnrs = beadnrs,
                        xs = lx,
                        ys = ly,
                        zs = lz,
                    )
                    self.lipid_dict[params][lipid_name].set_xy_to_center()
                    
                if data_specifier == "charges":
                    self.lipid_dict[params][lipid_name].set_charge(lipid_details)
            
        tot_lipids = sum([len(vals) for vals in self.lipid_dict.values()])
        
        self.print_term("Number of lipids preprocessed:", tot_lipids, "\n", spaces=1)
    
    def solute_beads_checker(self, res_dict):
        xs = res_dict["x"]
        ys = res_dict["y"]
        zs = res_dict["z"]
        beads = res_dict["beads"]
        if type(xs) != tuple and type(xs) != list:
            xs = (xs,)
        if type(ys) != tuple and type(ys) != list:
            ys = (ys,)
        if type(zs) != tuple and type(zs) != list:
            zs = (zs,)
        if type(beads) != tuple and type(beads) != list:
            beads = (beads,)
        
        assert len(beads) == len(xs) == len(ys) == len(zs), (
            "Number of beads, x-values, y-values and z-values must be the same within a residue." + "\n"
            "len(beads) " + str(len(beads)) + ", values: " + str(beads) + "\n"
            "len(xs)    " + str(len(xs)) + ", values: " + str(xs) + "\n"
            "len(ys)    " + str(len(ys)) + ", values: " + str(ys) + "\n"
            "len(zs)    " + str(len(zs)) + ", values: " + str(zs) + "\n"
        )
        
        lx = [xi * 10 for xi in xs]
        ly = [yi * 10 for yi in ys]
        lz = [zi * 10 for zi in zs]

        return beads, lx, ly, lz
    
    def solute_defs_checker(self, mol_class, cur_name, cur_dict):
        if "solvcount" in cur_dict.keys():
            mol_class.solvcount_set(cur_dict["solvcount"])
        
        if "charge" in cur_dict.keys():
            charge = 0
            if type(cur_dict["charge"]) in [int, float]:
                charge = cur_dict["charge"]
            elif type(cur_dict["charge"]) == str:
                isnumber, isinteger = self.is_number(cur_dict["charge"])
                if isnumber:
                    charge = ast.literal_eval(cur_dict["charge"])
                elif cur_dict["charge"] in self.itp_moltypes.keys():
                    charge = self.itp_moltypes[cur_dict["charge"]]["charge_sum"]
                else:
                    self.print_term("WARNING: moleculetype '" + cur_dict["charge"] + "' for solute name '" + cur_name + "' not found in topology")
            mol_class.charge_set(charge)

        if "density" in cur_dict.keys():
            mol_class.density_set(cur_dict["density"])

        if "molar_mass" in cur_dict.keys():
            mol_class.molar_mass_set(cur_dict["molar_mass"])
        
        if "residues" in cur_dict.keys():
            molbeadnrs = []
            for res_dict in cur_dict["residues"]:
                resname = res_dict["resname"]
                beads, xs, ys, zs = self.solute_beads_checker(res_dict)
                if molbeadnrs == []:
                    beadnrs = list([i for i in range(len(beads))])
                    molbeadnrs = beadnrs[:]
                else:
                    beadnrs = list([i for i in range(max(molbeadnrs)+1, len(molbeadnrs)-1+len(beads))])

                mol_class.add_res_and_beads(
                    resname = resname, beads = beads, beadnrs = beadnrs,
                    xs = xs, ys = ys, zs = zs
                )
        else:
            resname = cur_name
            beads, xs, ys, zs = self.solute_beads_checker(cur_dict)
            beadnrs = list([i for i in range(len(beads))])
            mol_class.add_res_and_beads(
                resname = resname, beads = beads, beadnrs = beadnrs,
                xs = xs, ys = ys, zs = zs
            )
        mol_class.set_coords_to_center()
    
    def solvent_defs_preprocessor(self):
        '''
        Preprocesses solvent defintions by rearranging the data into a class
        '''
        self.print_term("Preprocessing solvent definitions")
        ### Loop over generel lipid types (e.g. phospholipids, sterols, etc.)
        for params, solvent_type_dict in self.solvent_defs.items():
            if params not in self.solvent_dict.keys():
                self.solvent_dict[params] = {}
            ### Loop over indexes
            for cur_name, cur_dict in solvent_type_dict.items():
                if cur_name not in self.solvent_dict[params].keys():
                    self.solvent_dict[params][cur_name] = SOLVENT(charge = 0, solvcount = 1, molname = cur_name)
                else:
                    continue
                
                self.solute_defs_checker(self.solvent_dict[params][cur_name], cur_name, cur_dict)
                
        tot_solvents = sum([len(vals) for vals in self.solvent_dict.values()])
        self.print_term("Number of solvents preprocessed:", tot_solvents, "\n", spaces=1)

    def ion_defs_preprocessor(self):
        '''
        Preprocesses ion defintions by rearranging the data into a class
        '''
        self.print_term("Preprocessing ion definitions")
        ### Loop over generel lipid types (e.g. phospholipids, sterols, etc.)
        for params, charge_type_dict in self.ion_defs.items():
            if params not in self.ion_dict.keys():
                self.ion_dict[params] = {}
            for charge_name, ion_type_dict in charge_type_dict.items():
                if charge_name not in self.ion_dict[params].keys():
                    self.ion_dict[params][charge_name] = {}
                ### Loop over indexes
                for cur_name, cur_dict in ion_type_dict.items():
                    if cur_name not in self.ion_dict[params][charge_name].keys():
                        self.ion_dict[params][charge_name][cur_name] = SOLVENT(charge = 0, solvcount = 1, molname = cur_name)
                    else:
                        continue
                    
                    self.solute_defs_checker(self.ion_dict[params][charge_name][cur_name], cur_name, cur_dict)
                    
        tot_pos_ions = sum([len(vals["positive"]) for vals in self.ion_dict.values()])
        self.print_term("Number of positive ions preprocessed:", tot_pos_ions, spaces=1)
        tot_neg_ions = sum([len(vals["negative"]) for vals in self.ion_dict.values()])
        self.print_term("Number of negative ions preprocessed:", tot_neg_ions, spaces=1)
    
    ##############################################
    ### Protein/membrane/solvent preprocessors ###
    ##############################################
    def prot_preprocessor(self):
        '''
        Preprocesses protein commands for later ease of use
        '''
        self.PROTEINS = {}
        if len(self.PROTEINS_cmds) != 0:
            self.print_term("\nPreprocessing protein requests")
            for cmd_nr, prot_cmd in enumerate(self.PROTEINS_cmds, 1):
                self.print_term("Starting protein:", cmd_nr, spaces=1)

                ### Defaults
                prot_dict = {
                    "tx": 0,
                    "ty": 0,
                    "tz": 0,
                    "rx": 0,
                    "ry": 0,
                    "rz": 0,
                    "cen_method": ("cog",), # "cog/axis/ax" (center of geometry), "mean" (mean of all points), "bead:INT", "res:INT" or "point:x:y:z"
                    "lipids_inside": False, # [bool]
                    "pbc_check": True,
                    "buffer": 1.32, # [Å] default = (vdw of regular beads) / 2
                    "mol_names": [],
                    "charge": "top", # int/float, "top" or "auto"
                    "tot_charge": 0,
                }

                ### ### Check protein command
                for cmd in prot_cmd.split():
                    sub_cmd = cmd.split(":")

                    ### Read pdb/gro file
                    if any([extension in sub_cmd[0].lower() for extension in [".pdb", ".gro"]]) or sub_cmd[0] in ["f", "file", "prot_file", "protein_file"]:
                        ### Find file name
                        if sub_cmd[0] in ["f", "file", "prot_file", "protein_file"]:
                            file_name = sub_cmd[1]
                        else:
                            file_name = sub_cmd[0]
                        
                        ### First check if file ends with .pdb or .gro
                        if file_name.endswith(".pdb"):
                            prot_dict["beads"] = self.pdb_reader(file_name)
                        elif file_name.endswith(".gro"):
                            prot_dict["beads"] = self.gro_reader(file_name)
                        ### If neither then check if it contains .pdb or .gro (e.g. #prot_file.pdb.1#)
                        elif ".pdb" in file_name.lower():
                            prot_dict["beads"] = self.pdb_reader(file_name)
                        elif ".gro" in file_name.lower():
                            prot_dict["beads"] = self.gro_reader(file_name)
                        ### Finally assume .pdb format if neither .pdb nor .gro is found
                        else:
                            prot_dict["beads"] = self.pdb_reader(file_name)

                    
                    
                    ### Center method "cog/axis/ax" (center of geometry), "mean", "bead:INT", "res:INT" or "point:x:y:z"
                    elif sub_cmd[0].lower() == "cen_method":
                        if sub_cmd[1].lower() in ["cog", "axis", "ax", "mean"]:
                            prot_dict["cen_method"] = (sub_cmd[1].lower(),)
                        elif sub_cmd[1].lower() in ["bead", "res"]:
                            prot_dict["cen_method"] = (sub_cmd[1].lower(), sub_cmd[2:])
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
                            params = self.prot_params or self.sys_params
                            if values["res_name"] in self.prot_defs[params]["charges"]:
                                charge += self.prot_defs[params]["charges"][values["res_name"]]
                            else:
                                self.print_term(values["res_name"], "residue name not in 'prot_defs' charges for system parameters:", params)
                            counted_residues.append(res_nr)
                    return charge
                
                ### Post-preprocessing (topology and charge determination)
                if isinstance(prot_dict["charge"], (float, int)):
                    prot_dict["tot_charge"] += prot_dict["charge"]
                    
                elif prot_dict["charge"] == "auto":
                    prot_dict["tot_charge"] += prot_charge_finder(prot_dict["beads"])
                    
                elif prot_dict["charge"] == "top":
                    if prot_dict["mol_names"] == []:
                        self.print_term("No protein names given.  Will estimate charges from residue names. Use 'mol_name' to assign protein names (name used in topology files)", spaces=2)
                        prot_dict["tot_charge"] = prot_charge_finder(prot_dict["beads"])
                    else:
                        for mol_name in prot_dict["mol_names"]:
                            if mol_name not in self.itp_moltypes.keys():
                                self.print_term("A protein name could not be found in your topology file(s): " +  mol_name, warn=True, spaces=2)
                            else:
                                prot_dict["tot_charge"] += self.itp_moltypes[mol_name]["charge_sum"]

                prot_dict["mol_names"] = prot_dict["mol_names"] or "_".join(["PROT", str(cmd_nr)])
                
                ### Converting data into protein class
                prot_dict["protein"] = PROTEIN(charge = prot_dict["tot_charge"], molname = "_".join(["PROT", str(cmd_nr)]))
                cur_res = False
                for key, vals in prot_dict["beads"].items():
                    if vals["res_nr"] != cur_res:
                        prot_dict["protein"].add_res(vals["res_name"], vals["res_nr"])
                        cur_res = vals["res_nr"]
                    
                    prot_dict["protein"].add_bead(
                        vals["atom_name"],
                        vals["atom_nr"],
                        vals["x"],
                        vals["y"],
                        vals["z"],
                    )
                
                self.PROTEINS[cmd_nr] = prot_dict.copy()
            self.print_term("Number of molecule insertions preprocessed:", len(self.PROTEINS), spaces=1)

    def memb_preprocessor(self):
        '''
        Preprocesses membrane commands for later ease of use
        '''
        self.MEMBRANES = {}
        ### ### Preprocessing membranes
        if len(self.MEMBRANES_cmds) != 0:
            self.print_term("\nPreprocessing membrane requests")
            for cmd_nr, leaf_cmd in enumerate(self.MEMBRANES_cmds, 1):
                self.print_term("Starting membrane command:", cmd_nr, spaces=1)
                ### "membrane" settings always apply to the whole membrane instead of individual leaflets
                ### "default" settings can be overwritten by individual leaflets
                settings_dict = {
                    "upper_leaf": {},
                    "lower_leaf": {},
                    "membrane"  : {
#                         "shape": "rectangle", # "rectangle"
                        
                        "x": self.pbcx / 10, # [nm] converted to [Å]
                        "y": self.pbcy / 10, # [nm] converted to [Å]
                        "center": [0, 0, 0], # [nm] converted to [Å]
                        
                        "gridsplits": ("auto", 500), # Å
                        
                        "pbc_check": True, # [bool]
                        
                        "lipid_distribution": "random",
                        
                        "optimize": "v5", # False/"no"/0, "limited", True/"full"/1
                        "optim_maxsteps": 100,
                        "optim_push_tol": 0.2,
                        "optim_push_mult": 1.0,
                        
                    },
                    "default": {
                        "lipids": {}, # empty
                        "lipids_preprocessing": [], # empty

                        "kickx": 0.025, # [nm] converted to [Å]
                        "kicky": 0.025, # [nm] converted to [Å]
                        "kickz": 0.025, # [nm] converted to [Å]

                        "apl": 0.6, # [nm^2] converted to [Å^2]
                        
                        "plane_buffer": 0.099, # 0.066, # 0.132, # [nm] converted to [Å], default = (vdw of regular beads) / 4
                        "height_buffer": 0.099, # 0.066, # 0.132, # [nm] converted to [Å], default = (vdw of regular beads) / 4

                        "prot_buffer": 0.132, # [nm] converted to [Å], default = (vdw of regular beads) / 2
                        "alpha_mult": 1.0,
                        
                        "lip_round_func": ("int", int), # int, round, math.floor or math.ceil

                        "lipid_optim": 'avg_optimal', # str: 'abs_val', 'force_fill', 'fill', 'avg_optimal', 'no'
                        "params": False, # False or str
                        "charge": "top", # "lib" or "top"
                    }
                }

                
                ### ### Membrane mono/bilayer names
                ### Mono as explicitly upper added by request
                monolayer_upper_designation = ["u", "up", "upper", "m", "mo", "mono", "monolayer", "mono_upper"]
                monolayer_lower_designation = ["d", "do", "down", "l", "lo", "lower", "mono_lower", "mono_down"]
                bilayer_designation = ["b", "bi", "bilayer", "memb", "membrane"]
                
                layer_definition = "bilayer"
                ### Some values must always be for the membrane, as they make little sense to do per-leaflet
                dict_target      = "membrane"
                                
                ### ### Check leaf command
                ### Liberal use of ast.literal_eval() below to convert strings to int/float
                ### ast.literal_eval checks if code is a valid python datatype before interpreting it
                for cmd in leaf_cmd.split():
                    sub_cmd = cmd.split(":")
                    
                    ### Bilayer/monolayer defition
                    if sub_cmd[0].lower() in ["type"]:
                        if sub_cmd[1].lower() in bilayer_designation:
                            layer_definition = "bilayer"
                        elif sub_cmd[1].lower() in monolayer_upper_designation:
                            layer_definition = "upper"
                        elif sub_cmd[1].lower() in monolayer_lower_designation:
                            layer_definition = "lower"
                    
                    ### Sets whether following sub commands are for a specific leaflet or the whole membrane
                    elif sub_cmd[0].lower() in ["leaflet", "leaf", "side"]:
                        if sub_cmd[1].lower() in bilayer_designation:
                            dict_target = "membrane"
                        elif sub_cmd[1].lower() in monolayer_upper_designation:
                            dict_target = "upper_leaf"
                        elif sub_cmd[1].lower() in monolayer_lower_designation:
                            dict_target = "lower_leaf"

                    ### Leaflet shape
                    elif sub_cmd[0].lower() == "gridsplits":
                        ELSE = False
                        if len(sub_cmd[1:]) == 1:
                            val = sub_cmd[1]
                            isnumber, isint = self.is_number(val)
                            if isnumber:
                                if val == "0":
                                    val = "1"
                                settings_dict["membrane"]["gridsplits"] = (int(ast.literal_eval(val)), int(ast.literal_eval(val)))
                            elif val == "False":
                                settings_dict["membrane"]["gridsplits"] = (1, 1)
                            elif val.lower() == "auto":
                                settings_dict["membrane"]["gridsplits"] = ("auto", 500)
                            else:
                                ELSE = True
                        elif len(sub_cmd[1:]) == 2:
                            val1, val2 = sub_cmd[1:]
                            isnumber1, isint1 = self.is_number(val1)
                            isnumber2, isint2 = self.is_number(val2)
                            if val1 == "auto" and isnumber2:
                                settings_dict["membrane"]["gridsplits"] = ("auto", ast.literal_eval(val2))
                            elif isnumber1 or isnumber2:
                                ### If values are numbers or "False"
                                if val1 in ["False", "0"]:
                                    val1 = "1"
                                if val2 in ["False", "0"]:
                                    val2 = "1"
                                settings_dict["membrane"]["gridsplits"] = (ast.literal_eval(val1), ast.literal_eval(val2))
                            else:
                                ELSE = True
                        else:
                            ELSE = True
                        if ELSE:
                            self.print_term(
                                "WARNING:",
                                "Subcommand:", "'" + sub_cmd[0] + "'", "was given with invalid value(s): ", "'" + sub_cmd[1:] + "'", "\n",
                                "For subcommand gridsplits:values; values must be one of: 'number/False', 'number1/False:number2/False', 'auto', 'auto:number'",
                                warn=True
                            )

                    ### Center [nm]
                    elif any(sub_cmd[0].lower() == cen for cen in ["c", "cen", "center"]):
                        if len(sub_cmd[1:]) == 3:
                            settings_dict["membrane"]["center"] = [ast.literal_eval(i) for i in sub_cmd[1:]]
                        elif len(sub_cmd[1:]) == 2:
                            settings_dict["membrane"]["center"] = [ast.literal_eval(i) for i in sub_cmd[1:]] + [0]
                    
                    elif sub_cmd[0].lower() == "cx":
                        settings_dict["membrane"]["center"][0] = ast.literal_eval(sub_cmd[1])
                    
                    elif sub_cmd[0].lower() == "cy":
                        settings_dict["membrane"]["center"][1] = ast.literal_eval(sub_cmd[1])
                    
                    elif sub_cmd[0].lower() == "cz":
                        settings_dict["membrane"]["center"][2] = ast.literal_eval(sub_cmd[1])

                    ### Area per lipid [nm^2] converted to [Å^2]
                    elif sub_cmd[0].lower() == "apl":
                        settings_dict[dict_target]["apl"] = ast.literal_eval(sub_cmd[1])

                    ### Wiggle room. Minimum distance from a bead to the edge of a bin [Å]
                    elif sub_cmd[0].lower() in ["plane_wr", "plane_buffer"]:
                        settings_dict[dict_target]["plane_buffer"] = ast.literal_eval(sub_cmd[1])

                    elif sub_cmd[0].lower() in ["height_wr", "height_buffer"]:
                        settings_dict[dict_target]["height_buffer"] = ast.literal_eval(sub_cmd[1])

                    ### Integer multiplier for radius used in alphashape function [multiplier]
                    elif sub_cmd[0].lower() == "prot_buffer":
                        prot_dict[dict_target] = ast.literal_eval(sub_cmd[1])

                    ### Integer multiplier for radius used in alphashape function [multiplier]
                    elif sub_cmd[0].lower() == "alpha_mult":
                        prot_dict[dict_target] = ast.literal_eval(sub_cmd[1])

                    ### Rounding method for rounding number of lipids in leaflet [function]
                    elif sub_cmd[0].lower() == "round": # "int", "round", "min" or "max"
                        rounding_types = {"int": ("int", int), "round": ("round", round), "floor": ("math.floor", math.floor), "ceil": ("math.ceil", math.ceil)}
                        settings_dict[dict_target]["lip_round_func"] = rounding_types[sub_cmd[1]] 

                    elif sub_cmd[0].lower() in ["kick", "kickx", "kicky", "kickz"]:
                        ### Random kick to beads x/y/z positions [nm] #[Å]
                        if sub_cmd[0].lower() == "kick":
                            settings_dict[dict_target]["kickx"] = ast.literal_eval(sub_cmd[1])
                            settings_dict[dict_target]["kicky"] = ast.literal_eval(sub_cmd[1])
                            settings_dict[dict_target]["kickz"] = ast.literal_eval(sub_cmd[1])
                        else:
                            settings_dict[dict_target][sub_cmd[0].lower()] = ast.literal_eval(sub_cmd[1])

                    ### Bead distance scaling for z is 0.3 by default and 0.25 by default for x and y [multiplier]
                    elif sub_cmd[0].lower() in ["bdx", "bdy", "bdz"]:
                        settings_dict[dict_target][sub_cmd[0].lower()] = ast.literal_eval(sub_cmd[1])

                    ### ### Rectangle specific subcommands:
                    ### xy dimensions of leaflet [nm]. Converted to [Å]
                    elif sub_cmd[0].lower() in ["x", "y"]:
                        settings_dict["membrane"][sub_cmd[0].lower()] = ast.literal_eval(sub_cmd[1])

                    ### True/False whether to check for pbc crossing
                    elif sub_cmd[0].lower() == "pbc_check":
                        settings_dict["membrane"]["pbc_check"] = ast.literal_eval(sub_cmd[1])

                    ### Lipid optimization method
                    elif sub_cmd[0].lower() == "lipid_optim":
                        valid_lipid_optim = sub_cmd[1] in ['avg_optimal', 'abs_val', 'force_fill', 'fill', 'no']
                        assert valid_lipid_optim, "Invalid lipid selection method: '" + str(sub_cmd[1]) + "'"
                        settings_dict[dict_target]["lipid_optim"] = sub_cmd[1]

                    ### Designates which force field to collect lipids from
                    elif sub_cmd[0].lower() == "params":
                        settings_dict[dict_target]["params"] = ast.literal_eval(sub_cmd[1])

                    ### Designates how lipid charges should be collected
                    elif sub_cmd[0].lower() == "charge": # "lib" or "top"
                        settings_dict[dict_target]["charge"] = sub_cmd[1]
                    
                    elif sub_cmd[0].lower() in ["ld", "lipid_dist", "lipid_distribution"]:
                        val = sub_cmd[1]
                        if val in ["e", "even", "evenly"]:
                            settings_dict[dict_target]["lipid_distribution"] = "evenly"
                        elif val in ["r", "ran", "rand", "random"]:
                            settings_dict[dict_target]["lipid_distribution"] = "random"
                    
                    #############################
                    ### OPTIMIZATION SETTINGS ###
                    #############################
                    elif sub_cmd[0].lower() in ["optim", "optimize", "optimization", "minim", "minimize", "minimization"]:
                        if sub_cmd[1] in ["False", "0", "no"]:
                            settings_dict[dict_target]["optimize"] = False
                        else:
                            settings_dict[dict_target]["optimize"] = sub_cmd[1].lower()
#                         if sub_cmd[1] in ["True", "1", "full"]:
#                             settings_dict[dict_target]["optimize"] = "full"
#                         elif sub_cmd[1] in ["limited"]:
#                             settings_dict[dict_target]["optimize"] = "limited"
#                         elif sub_cmd[1] in ["False", "0", "no"]:
#                             settings_dict[dict_target]["optimize"] = False
                    
                    elif sub_cmd[0].lower() in ["optim_maxsteps", "minim_maxsteps"]:
                        isnumber, isint = self.is_number(sub_cmd[1])
                        if isnumber:
                            settings_dict[dict_target]["optim_maxsteps"] = int(ast.literal_eval(sub_cmd[1]))
                    
                    elif sub_cmd[0].lower() in ["optim_push_tol", "minim_push_tol"]:
                        isnumber, isint = self.is_number(sub_cmd[1])
                        if isnumber:
                            settings_dict[dict_target]["optim_push_tol"] = ast.literal_eval(sub_cmd[1])
                    
                    elif sub_cmd[0].lower() in ["optim_push_mult", "minim_push_mult"]:
                        isnumber, isint = self.is_number(sub_cmd[1])
                        if isnumber:
                            settings_dict[dict_target]["optim_push_mult"] = ast.literal_eval(sub_cmd[1])
                    
                    ### Processes only lipids in lipid dictionary
                    elif any([sub_cmd[0] in self.lipid_dict[params].keys() for params in self.lipid_dict.keys()]):
                        if "lipids_preprocessing" not in settings_dict[dict_target].keys():
                            settings_dict[dict_target]["lipids_preprocessing"] = []
                        settings_dict[dict_target]["lipids_preprocessing"].append(sub_cmd)
                    ### Errors out if unknown subcommand used, and prints the subcommand to terminal
                    else:
                        assert False, "Unknown subcommand given to '-membrane'. The subcommand is: '" + str(sub_cmd) + "'"
                
                if layer_definition == "bilayer" and len(settings_dict["default"]["lipids_preprocessing"]) == 0:
                    if len(settings_dict["upper_leaf"]) > 0 and len(settings_dict["lower_leaf"]) == 0:
                        layer_definition = "upper"
                    if len(settings_dict["lower_leaf"]) > 0 and len(settings_dict["upper_leaf"]) == 0:
                        layer_definition = "lower"
                
                memb_dict = {"leaflets": {}, "membrane_type": layer_definition}
                
                if layer_definition == "upper" or layer_definition == "bilayer":
                    memb_dict["leaflets"]["upper_leaf"] = settings_dict["upper_leaf"]
                    memb_dict["leaflets"]["upper_leaf"].update({
                        "HG_direction": "up",
                        "leaflet_type": "upper",
                    })
                
                if layer_definition == "lower" or layer_definition == "bilayer":
                    memb_dict["leaflets"]["lower_leaf"] = settings_dict["lower_leaf"]
                    memb_dict["leaflets"]["lower_leaf"].update({
                        "HG_direction": "down",
                        "leaflet_type": "lower",
                    })
                
                ### Adding membrane-wide settings to specific leaflets if they are not already given for the specific leaflet
                for key, vals in settings_dict["membrane"].items():
                    for leaflet in memb_dict["leaflets"].values():
                        if key not in leaflet:
                            try:
                                leaflet[key] = vals.copy()
                            except:
                                leaflet[key] = vals
                
                ### Adding default settings for leaflets if none were given
                for key, vals in settings_dict["default"].items():
                    for leaflet in memb_dict["leaflets"].values():
                        if key not in leaflet:
                            try:
                                leaflet[key] = vals.copy()
                            except:
                                leaflet[key] = vals
                
                ### fixing a couple of values such that they are in angstrom
                for leaflet_type, leaflet in memb_dict["leaflets"].items():
                    leaflet["x"] *= 10
                    leaflet["y"] *= 10
                    leaflet["center"] = [val*10 for val in leaflet["center"]]
                    leaflet["apl"] *= 100
                    leaflet["prot_buffer"] *= 10
                    leaflet["plane_buffer"] *= 10
                    leaflet["height_buffer"] *= 10
                    leaflet["kickx"] *= 10
                    leaflet["kicky"] *= 10
                    leaflet["kickz"] *= 10
                    
                ################################
                ### Lipid data incorporation ###
                ################################
                
                ### Reconfigures lipid-specific data for leaflet according to subcommands
                for leaflet_key, leaflet in memb_dict["leaflets"].items():
                    assert len(leaflet["lipids_preprocessing"]) > 0, "No lipids given to '-membrane' flag. Please specify at least one lipid if you use it."
                    for l_name, *rest in leaflet["lipids_preprocessing"]:
                        """
                        Checks if parameters specified for lipid.
                        Order: lipid specific parameters, leaflet specific parameters, general system parameters (defaults to "default")
                        """
                        if len(rest) > 1:
                            l_ratio, l_params = rest
                        elif len(rest) == 0:
                            l_ratio  = 1 
                            l_params = leaflet["params"] or self.lipid_params or self.sys_params # (sys_params defaults to "default")
                        else:
                            l_ratio = rest[0]
                            l_params = leaflet["params"] or self.lipid_params or self.sys_params # (sys_params defaults to "default")

                        leaflet["lipids"][l_name] = self.lipid_dict[l_params][l_name]

                        if type(l_ratio) == str:
                            l_ratio = ast.literal_eval(l_ratio)
                        leaflet["lipids"][l_name].ratio_add(l_ratio)
                            
                        if leaflet["charge"] == "top" and l_name in self.itp_moltypes.keys():
                            leaflet["lipids"][l_name].charge_set(self.itp_moltypes[l_name]["charge_sum"])
                        
                    maxx, maxy, maxz, minx, miny, minz = [
                        func([
                            val + (leaflet["kick" + ax] + leaflet[wr + "_buffer"]) * sign
                            for lipid in leaflet["lipids"].keys()
                            for val in leaflet["lipids"][lipid].get_coords(ax)[0]
                        ])
                        for func, sign in [(max, +1), (min, -1)] for ax, wr in [("x", "plane"), ("y", "plane"), ("z", "height")]
                    ]
                    
                    ### Update leaflet dictionary. Some are duplicated here, but potentially overwritten later
                    leaflet["lipid_dimensions"] = {
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
                
                memb_dict.update({
                    "maxz"           : max([leaflet["lipid_dimensions"]["maxz"]         for leaflet in memb_dict["leaflets"].values()]),
                    "minz"           : min([leaflet["lipid_dimensions"]["minz"]         for leaflet in memb_dict["leaflets"].values()]),
                    "membrane_height": sum([leaflet["lipid_dimensions"]["lipid_height"] for leaflet in memb_dict["leaflets"].values()]),
                })
                    
                ### Upper leaflet monolayer
                if layer_definition in monolayer_upper_designation:
                    self.MEMBRANES[cmd_nr] = memb_dict
                ### Lower leaflet monolayer
                elif layer_definition in monolayer_lower_designation:
                    self.MEMBRANES[cmd_nr] = memb_dict
                ### Bilayer
                elif layer_definition in bilayer_designation:
                    self.MEMBRANES[cmd_nr] = memb_dict
            
            self.print_term("Number of membranes preprocessed:", len(self.MEMBRANES), spaces=1)

    def solv_preprocessor(self):
        '''
        Preprocesses solvation commands for later ease of use
        '''
        self.SOLVATIONS = {}
        
        if len(self.SOLVATIONS_cmds) != 0:
            self.print_term("\nPreprocessing solvent requests")
            for cmd_nr, solvate_cmd in enumerate(self.SOLVATIONS_cmds, 1):
                self.print_term("Starting Solvent command:", cmd_nr, spaces=1)
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
                    
                    ### ### Flooding specific
                    "flooding": False,
                    
                    ### ### General
                    "count": False, # [bool] Uses specific molarity value as absolute number of molecules instead of molarity ratio. Will be rounded using int(val+0.5)
                    "kick": 2.64/5*2, # [Å]
                    "bdx": 1.0, # [multiplier]
                    "bdy": 1.0, # [multiplier]
                    "bdz": 1.0, # [multiplier]
                    "params": False, # False or str
                    "bead_radius": 2.64, # [Å]
                    "gridres": 2.64, # [Å] # 1.32
                    "WR": 2.64, # [Å]
                    "buffer": 2.0, # [Å]
                    "protein_extra_buffer": 2,
                    "lipid_extra_buffer": 0,
                    "solute_extra_buffer": 0,
                    "min_cellsize": 1, # Integer value
                    
                    ### ### Choose solvation algorithm
                    "algorithm": "v2", # (v2: Fast), (v1: old, same result but slower, kept for now just in case)
                }
                
                ### Added "default" solvation command per user request
                if solvate_cmd == "":
                    solvate_cmd = "solv:W pos:NA neg:CL"
                
                for cmd in solvate_cmd.split():
                    sub_cmd = cmd.split(":")
                    ########################
                    ### SOLVENT SPECIFIC ###
                    ########################
                    ### Solvents
                    if sub_cmd[0].lower() in ["solv", "mol"]:
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
                        solv_dict["ionsvol"] = sub_cmd[1].lower()

                    ###############
                    ### GENERAL ###
                    ###############
                    ### Whether to use ratios as absolute number of molecules. True/False
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
                        
                    ### Bead scaling distances
                    ### Bead distance scaling for z is 0.3 by default and 0.25 by default for x and y [multiplier]
                    elif sub_cmd[0].lower() in ["bdx", "bdy", "bdz"]:
                        leaf_dict[sub_cmd[0].lower()] = ast.literal_eval(sub_cmd[1])

                    ### Designates force field
                    elif sub_cmd[0].lower() == "params":
                        solv_dict["params"] = sub_cmd[1]

                    ### Gerenal buffer for solvent placement
                    elif sub_cmd[0].lower() == "buffer":
                        solv_dict["buffer"] = ast.literal_eval(sub_cmd[1])

                    ### Extra buffer given to proteins
                    elif sub_cmd[0].lower() == "protein_extra_buffer":
                        solv_dict["protein_extra_buffer"] = ast.literal_eval(sub_cmd[1])

                    ### Extra buffer given to lipids
                    elif sub_cmd[0].lower() == "lipid_extra_buffer":
                        solv_dict["lipid_extra_buffer"] = ast.literal_eval(sub_cmd[1])

                    ### Extra buffer given to prior solvents (solutes)
                    elif sub_cmd[0].lower() == "solute_extra_buffer":
                        solv_dict["solute_extra_buffer"] = ast.literal_eval(sub_cmd[1])

                    ### Buffer for solvent placement from leaflet hydrophobic area
                    elif sub_cmd[0].lower() == "min_cellsize":
                        solv_dict["min_cellsize"] = int(ast.literal_eval(sub_cmd[1]) + 0.5)

                    ### Chooses solvation algorithm
                    elif sub_cmd[0].lower() in ["alg", "algorithm"]:
                        solv_dict["algorithm"] = sub_cmd[1]
                    
                    #########################
                    ### FLOODING SPECIFIC ###
                    #########################
                    ### Used to check if command is made as flooding or not
                    elif sub_cmd[0].lower() == "flooding":
                        solv_dict["flooding"] = ast.literal_eval(sub_cmd[1])
                        
                    elif solv_dict["flooding"] and sub_cmd[0] in [key for params in solvent_defs.keys() for key in solvent_defs[params].keys()]:
                        solv_dict["solv_preprocessing"].append(sub_cmd[0:])
                    
                    ### Errors out if unknown subcommand used, and prints the subcommand to terminal
                    else:
                        assert False, "Unknown subcommand given to '-solvate'. The subcommand is: '" + str(cmd) + "'"
                
                ######################################
                ### SOLVENT/ION DATA INCORPORATION ###
                ######################################
                assert solv_dict["solv_preprocessing"] != {}, "No solvent given to '-solvate' flag. Please specify at least one non-ionic solvent if you use it."

                params = solv_dict["params"] or self.solv_params or self.sys_params # (sys_params defaults to "default")
                
                def values_checker(subcmd_rest):
                    rest_params = False
                    rest_molarity = False
                    if len(subcmd_rest) > 0:
                        for val in subcmd_rest:
                            isnumber, isinteger = self.is_number(rest[-1])
                            if isnumber:
                                rest_molarity = ast.literal_eval(val)
                            else:
                                rest_params = val
                    return rest_params, rest_molarity 
                    
                ### ### SOLVENT DATA
                for subcommand_values in solv_dict["solv_preprocessing"]:
                    name = subcommand_values[0]
                    rest = subcommand_values[1:]
                    solv_type = "solvent"
                    rest_params, rest_molarity = values_checker(rest)
                    if rest_params:
                        params = rest_params
                    else:
                        params = solv_dict["params"] or self.solv_params or self.sys_params # (sys_params defaults to "default")
                    ### Adding name:dict combo to solvent dict
                    solv_dict[solv_type][name] = self.solvent_dict[params][name]
                    
                    ### If count has been set, then convert to integer value and treat as absolute number of molecules
                    if solv_dict["count"] and rest_molarity:
                        solv_dict[solv_type][name].molarity_set(int(ast.literal_eval(rest[-1]) + 0.5))
                    elif solv_dict["count"]:
                        solv_dict[solv_type][name].molarity_set(int(solv_dict["solv_molarity"] + 0.5))
                    else:
                        solv_dict[solv_type][name].molarity_set(solv_dict["solv_molarity"])
                
                ion_types_to_be_processed = []
                if solv_dict["pos_ions_preprocessing"]:
                    ion_types_to_be_processed.extend([("pos_ions", ions) for ions in solv_dict["pos_ions_preprocessing"]])
                if solv_dict["neg_ions_preprocessing"]:
                    ion_types_to_be_processed.extend([("neg_ions", ions) for ions in solv_dict["neg_ions_preprocessing"]])
                    
                ### ### Ion processing
                for solv_type, subcommand_values in ion_types_to_be_processed:
                    name = subcommand_values[0]
                    rest = subcommand_values[1:]
                    rest_params, rest_molarity = values_checker(rest)
                    if rest_params:
                        params = rest_params
                    else:
                        params = solv_dict["params"] or self.solv_params or self.sys_params # (sys_params defaults to "default")
                    ### Adding name:dict combo to solvent dict
                    if solv_type == "pos_ions":
                        solv_dict[solv_type][name] = self.ion_dict[params]["positive"][name]
                    elif solv_type == "neg_ions":
                        solv_dict[solv_type][name] = self.ion_dict[params]["negative"][name]
                    
                    rest_params, rest_molarity = values_checker(rest)
                    ### If count has been set, then convert to integer value and treat as absolute number of molecules
                    if solv_dict["count"] and rest_molarity:
                        solv_dict[solv_type][name].molarity_set(int(ast.literal_eval(rest[-1]) + 0.5))
                    elif solv_dict["count"]:
                        solv_dict[solv_type][name].molarity_set(int(solv_dict["salt_molarity"] + 0.5))
                    else:
                        solv_dict[solv_type][name].molarity_set(solv_dict["salt_molarity"])

                self.SOLVATIONS[cmd_nr] = solv_dict
                
            self.print_term("Number of solvent commands preprocessed:", len(self.SOLVATIONS), spaces=1)
    
    def itp_read_initiater(self):
        if len(self.ITP_INPUT_cmds) != 0:
            self.print_term("\nLoading topology file(s)", spaces=0)
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

            self.print_term("Finished loading topologies. Number of moleculetypes found:", len(self.itp_moltypes), spaces=1)
    
    def import_structures_handler(self):
        if len(self.SOLUTE_INPUT_cmds) != 0:
            for imp_struc in self.SOLUTE_INPUT_cmds:
                ### ### V1: Single subcommand:
                ### mol:name:n_residues:charge
                ### name can be either "resname" or a string of up to 4 characters
                ### n_residues must be integer value
                ### charge must be either int/float value or string of topology moleculetype
                ### each "mol" call is new molecule
                ### if no charge given, assume top using name
                ### if no n_residues given along with no charge then assume 1 residue
                ### ### Default:
                ### If no "mol" subcommands given. Assume all are single-residue molecules
                ### Residue name used for "name" and "moleculetype"
                ### If no moleculetype found for "name", then assume 0
                ### Examples: "ligands.pdb   mol:RHO:1:rhodamine   mol:PEP1:8:peptide_1   mol:resname:1:2"

                ### ### V2: multiple subcommands:
                ### names: list of names
                ### n_residues/nres: list of number of residues
                ### charges: list of charge designations
                ### Same number of values must be given to each of above subcommands
                ### ### Default:
                ### names: required unless n_residues 1 for all residues, for which resname will be used
                ### n_residues/nres: assumed 1 residue for all molecules
                ### charges: assumed "top" using names. Defaults to charge of 0 if no top found
                
                structures = [] # empty
                params = "default" # str
                
                mol_import_settings = []
                names_list = []
                nress_list = []
                charges_list = []
                ### ### settings dictionary:
                ### name: reference name for solvent selection
                ### nres: number of residues in molecule
                ### charge: whether to use name for topology or given value
                
                def name_checker(name):
                    ### Limiting name to 4 may not be necessary as residue names are used later
#                     assert len(name) <=4 or name == "resname", (
#                         "name can at most be 4 characters long or 'resname'" + "\n"
#                         "name: " + str(name)
#                     )
                    return name
                
                def nres_checker(nres):
                    nres = ast.literal_eval(nres)
                    return nres
                
                def charge_checker(charge):
                        isnum, isint = self.is_number(charge)
                        if isnum:
                            charge = ast.literal_eval(charge)
                        return charge
                
                ### Check for use of V1 or V2 system. User must choose one version for a given command.
                V1_check = [1 for i in imp_struc.split() if i.startswith("mol") and not any([i.endswith(j) for j in ["pdb", "gro"]])]
                V2_check = [1 for i in imp_struc.split() if any([i.startswith(j) for j in ["names", "nres", "n_residues", "charges"]]) and not any([i.endswith(j) for j in ["pdb", "gro"]])]
                assert not (sum(V1_check) > 0 and sum(V2_check) > 0), "Both V1 and V2 import commands used. Stick to one."
                
                for cmd in imp_struc.split():
                    sub_cmd = cmd.split(":")
                    ### Read pdb/gro file
                    if any([extension in sub_cmd[0].lower() for extension in [".pdb", ".gro"]]) or sub_cmd[0] in ["f", "file", "solute_file"]:
                        
                        ### Find file name
                        if sub_cmd[0] in ["f", "file", "solute_file"]:
                            file_name = sub_cmd[1]
                        else:
                            file_name = sub_cmd[0]
                        
                        ### First check if file ends with .pdb or .gro
                        if file_name.endswith(".pdb"):
                            structures.append(self.pdb_reader(file_name))
                        elif file_name.endswith(".gro"):
                            structures.append(self.gro_reader(file_name))
                        ### If neither then check if it contains .pdb or .gro (e.g. #prot_file.pdb.1#)
                        elif ".pdb" in file_name.lower():
                            structures.append(self.pdb_reader(file_name))
                        elif ".gro" in file_name.lower():
                            structures.append(self.gro_reader(file_name))
                        ### Finally assume .pdb format if neither .pdb nor .gro is found
                        else:
                            structures.append(self.pdb_reader(file_name))
                    
                    elif sub_cmd[0] == "params":
                        params = sub_cmd[1]
                    
                    elif sub_cmd[0] == "mol":
                        assert len(sub_cmd[1:]) == 3, "mol subcommand does not contain 3 values: " + str(sub_cmd)
                        name, nres, charge = sub_cmd[1:]
                        name   = name_checker(name)
                        nres   = nres_checker(nres)
                        charge = charge_checker(charge)
                        mol_import_settings.append({"name": name, "nres": nres, "charge": charge})
                        
                    elif sub_cmd[0] in ["name", "names"]:
                        names = sub_cmd[1:]
                        for name in names:
                            names_list.append(name_checker(name))
                        
                    elif sub_cmd[0] in ["n_residues", "residues", "nres"]:
                        nress = sub_cmd[1:]
                        for nres in nress:
                            nress_list.append(nres_checker(nres))
                        
                    elif sub_cmd[0] in ["charge", "charges"]:
                        charges = sub_cmd[1:]
                        for charge in charges:
                            charges_list.append(charge_checker(charge))
                
                assert len(names_list) == len(nress_list) == len(charges_list), "The number of names, n_residues and charges must be equal: names:{names}, n_residues:{nres}, charges:{charges}".format(names=names_list, nres=nress_list, charges=charges_list)
                for name, nres, charge in zip(names_list, nress_list, charges_list):
                    mol_import_settings.append({"name": name, "nres": nres, "charge": charge})
                
                assert len(structures) > 0, "No structure files given (.pdb and .gro are accepted file formats)"
                
                ### Residue sorting from across different structure files
                structures_residues_combined = []
                for struc_i, structure in enumerate(structures):
                    resnumber = "unset"
                    for key, values in structure.items():
                        if resnumber != key[2]:
                            if resnumber != "unset":
                                structures_residues_combined.append(reslist)
                            resnumber = key[2]
                            reslist   = [values]
                        else:
                            reslist.append(values)
                    structures_residues_combined.append(reslist)
                
                ### Combining settings and their residues from the structure file(s)
                molecule_residues = []
                current_res = 0
                for si, settings in enumerate(mol_import_settings):
                    molecule_residues.append(settings)
                    molecule_residues[si]["residues"] = []
                    nres = settings["nres"]
                    for res in structures_residues_combined[current_res:current_res+nres]:
                        molecule_residues[si]["residues"].append(res)
                    current_res += nres
                
                if params not in self.solvent_defs.keys():
                    self.solvent_defs[params] = {}
                if params not in self.ion_defs.keys():
                    self.ion_defs[params] = {"positive": {}, "negative": {}}
                
                ### Adding solute data to solvent and ion definition dictionaries
                for mi, molecule in enumerate(molecule_residues):
                    molname, nres, charge, residues = itemgetter("name", "nres", "charge", "residues")(molecule)
                    mol_dict = {
                        "residues": [],
                        "charge": charge,
                    }
                    for residue in residues:
                        beads, resnames, xs, ys, zs = [], [], [], [], []
                        for atom in residue:
                            bead, resname, x, y, z = itemgetter("atom_name", "res_name", "x", "y", "z")(atom)
                            beads.append(bead)
                            if resname not in resnames:
                                resnames.append(resname)
                            xs.append(x/10)
                            ys.append(y/10)
                            zs.append(z/10)
                        assert len(resnames) == 1, (
                            "Make sure only one residue name is used for each residue" + "\n"
                            "" + str(resnames)
                        )
                        mol_dict["residues"].append({
                            "resname": resnames[0],
                            "beads": tuple(beads),
                            "x": tuple(xs),
                            "y": tuple(ys),
                            "z": tuple(zs),
                        })
                    self.solvent_defs[params][molname]         = mol_dict.copy()
                    self.ion_defs[params]["positive"][molname] = mol_dict.copy()
                    self.ion_defs[params]["negative"][molname] = mol_dict.copy()
    
    #########################
    ### GENERAL FUNCTIONS ###
    #########################
    def rotate_point(self, posx, posy, posz, x_ang_deg, y_ang_deg, z_ang_deg):
        '''
        Rotates a point in 3D space using its original position and angles in degrees.
        '''
        ### Convert angles to radians
        x_ang = math.radians(x_ang_deg)
        y_ang = math.radians(y_ang_deg)
        z_ang = math.radians(z_ang_deg)

        ### Calculate rotation matrices
        rx = [[1               , 0               , 0               ],
              [0               , math.cos(x_ang) , -math.sin(x_ang)],
              [0               , math.sin(x_ang) , math.cos(x_ang) ]]
        
        ry = [[math.cos(y_ang) , 0               , math.sin(y_ang) ],
              [0               , 1               , 0               ],
              [-math.sin(y_ang), 0               , math.cos(y_ang) ]]
        
        rz = [[math.cos(z_ang) , -math.sin(z_ang), 0               ],
              [math.sin(z_ang) , math.cos(z_ang) , 0               ],
              [0               , 0               , 1               ]]
        
#         rx = [[1, 0, 0]                             , [0, math.cos(x_ang) , -math.sin(x_ang)], [0, math.sin(x_ang), math.cos(x_ang)]]
#         ry = [[math.cos(y_ang), 0, math.sin(y_ang)] , [0, 1, 0]                              , [-math.sin(y_ang)  , 0, math.cos(y_ang)]]
#         rz = [[math.cos(z_ang), -math.sin(z_ang), 0], [math.sin(z_ang) , math.cos(z_ang), 0] , [0, 0, 1]]

        ### Combine rotation matrices
        rxy  = [[sum(a * b for a, b in zip(row, col)) for col in zip(*ry)] for row in rx]
        rxyz = [[sum(a * b for a, b in zip(row, col)) for col in zip(*rz)] for row in rxy]

        ### Calculate new position of point
        new_posx, new_posy, new_posz = [rxyz[ax][0] * posx + rxyz[ax][1] * posy + rxyz[ax][2] * posz for ax in [0, 1, 2]]

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
            self.print_term("File " + output_file_name + " already exists. Backing it up", spaces=0)
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
    
    ########################
    ### PROTEIN INSERTER ###
    ########################
    def prot_placer(self):
        '''
        Places all proteins into the systems internal coordinate system
        Checks if all atoms/beads are within the pbc and moves them if they are outside
        '''
        if len(self.PROTEINS) != 0:
            prot_placer_tic = time.time()
            self.print_term("------------------------------ PROTEIN PLACEMENT", spaces=0)
            for protein_i, (protein_nr, protein) in enumerate(self.PROTEINS.items()):
                if protein_i != 0:
                    self.print_term("")
                self.print_term("Starting protein nr", protein_nr, spaces=0)
                
                #################
                ### CENTERING ###
                #################
                ### Centered on the mean coordinate of all beads (center of geometry)
                if protein["cen_method"][0] in ["cog", "axis", "ax"]:
                    centering = "cog"
                    target = False

                ### Centered on the mean of largest/smallest x/y/z coordinate
                if protein["cen_method"][0] == "mean":
                    centering = "mean"
                    target = False

                ### Centered on a single bead
                if protein["cen_method"][0] == "bead":
                    centering = "beadnr"
                    target = protein["cen_method"][1]

                ### Centered on the mean position of all beads in a single residue
                if protein["cen_method"][0] == "res":
                    centering = "resnr"
                    target = protein["cen_method"][1]

                ### Centered on the specific x/y/z coordinates
                if protein["cen_method"][0] == "point":
                    centering = "vals"
                    target = protein["cen_method"][1:]
                
                xcen, ycen, zcen = self.PROTEINS[protein_nr]["protein"].get_center_point(centering = centering, target = target)
                self.print_term("Centering protein using", "'" + " ".join([str(i) for i in protein["cen_method"]])+"'", "at x/y/z:", round(xcen, 3), round(ycen, 3), round(zcen, 3), "(Input file coordinate system [Å])", spaces=1)
                self.PROTEINS[protein_nr]["protein"].set_coords_to_center(centering = centering, target = target)

                #################
                ### ROTATIONS ###
                #################
                x_deg, y_deg, z_deg = protein["rx"], protein["ry"], protein["rz"]
                
                if any([ang != 0 for ang in [x_deg, y_deg, z_deg]]):
                    self.PROTEINS[protein_nr]["protein"].rotate_coords(rotation = [x_deg, y_deg, z_deg])

                ####################
                ### TRANSLATIONS ###
                ####################
                tx, ty, tz = protein["tx"], protein["ty"], protein["tz"]
                if any([ax != 0 for ax in [tx, ty, tz]]):
                    self.PROTEINS[protein_nr]["protein"].move_coords(translation = [tx, ty, tz])

                #################################################
                ### CHECKS IF COORDINATES ARE OUTSIDE THE BOX ###
                #################################################
                if protein["pbc_check"]:
                    beads = self.PROTEINS[protein_nr]["protein"].get_res_beads_info()
                    
                    errors_count = 0
                    for bead in beads:
                        bead_coords = [bead.x, bead.y, bead.z]
                        checked_beads, error = self.coord_checker(bead_coords, self.pbc_box, error_count = True)
                        if bead_coords != checked_beads:
                            bi = bead.bead
                            ri = bead.resnr
                            self.PROTEINS[protein_nr]["protein"].residues[ri].beads[bi].move_atom(checked_beads)
                            
                        if error > 0:
                            errors_count += 1
                    if errors_count > 0:
                        self.print_term("WARNING:", str(errors_count), "Beads are outside pbc. Moved to other side. Expect potential problems from this.", warn = True, spaces=2)
                        self.print_term("Please move the protein such that it fits within the pbc.", warn = True, spaces=2)

                xcen_new, ycen_new, zcen_new = 0, 0, 0
                xcen_new, ycen_new, zcen_new = xcen_new + tx, ycen_new + ty, zcen_new + tz
                self.print_term("New protein center at x/y/z:", round(xcen_new, 3), round(ycen_new, 3), round(zcen_new, 3), "(Internal coordinate system [Å])", spaces=1)

                self.system_charge += self.PROTEINS[protein_nr]["protein"].charge
                
                self.print_term("Finished placing protein nr", protein_nr, spaces=1)
                
            prot_placer_toc = time.time()
            prot_placer_time = round(prot_placer_toc - prot_placer_tic, 4)
            self.print_term("------------------------------ PLACEMENT COMPLETE", "(time spent: "+str(prot_placer_time)+" [s])", "\n", spaces=0)
    
    ######################
    ### POLYGON MAKERS ###
    ######################
    def subleaflet_poly_maker(self):
        if len(self.MEMBRANES) != 0:
            subleaflet_poly_maker_tic = time.time()
            self.print_term("------------------------------ CREATING LEAFLET BOUNDARY BOXES", spaces=0)
            for memb_i, (memb_key, memb_dict) in enumerate(self.MEMBRANES.items()):
                if memb_i != 0:
                    self.print_term("")
                self.print_term("Starting membrane nr", memb_key, spaces=0)
                
                for leaflet_i, (leaflet_key, leaflet) in enumerate(memb_dict["leaflets"].items()):
                    leaflet["subleaflets"] = {}
                    self.print_term("Starting", leaflet["leaflet_type"], "leaflet", spaces=1)
                    self.print_term("Base parameters:", "x=" + str(leaflet["x"] / 10) + "nm", "y=" + str(leaflet["y"] / 10) + "nm", "APL=" + str(leaflet["apl"] / 100) + "nm^2", spaces=2)
                    
                    if leaflet["gridsplits"][0] == "auto":
                        xsplits = leaflet["x"] // leaflet["gridsplits"][1] + 1
                        ysplits = leaflet["y"] // leaflet["gridsplits"][1] + 1
                    else:
                        xsplits, ysplits = leaflet["gridsplits"]
                    if xsplits > 1:
                        self.print_term("x-axis split into", xsplits, "subleaflets due to axis length or manual designation", spaces=2)
                    if ysplits > 1:
                        self.print_term("y-axis split into", ysplits, "subleaflets due to axis length or manual designation", spaces=2)
                    if xsplits * ysplits > 1:
                        self.print_term("A total of", xsplits*ysplits, "subleaflets have been made for the current leaflet", spaces=2)
                    
                    xmin, xmax = leaflet["center"][0] - leaflet["x"]/2, leaflet["center"][0] + leaflet["x"]/2
                    ymin, ymax = leaflet["center"][1] - leaflet["y"]/2, leaflet["center"][1] + leaflet["y"]/2
                    
                    xpoint_vals = np.linspace(start=xmin, stop=xmax, num=round(xsplits+1), endpoint=True)
                    ypoint_vals = np.linspace(start=ymin, stop=ymax, num=round(ysplits+1), endpoint=True)

                    for xi in range(round(xsplits)):
                        for yi in range(round(ysplits)):
                            xmin = xpoint_vals[xi]
                            xmax = xpoint_vals[xi+1]
                            ymin = ypoint_vals[yi]
                            ymax = ypoint_vals[yi+1]
                            
                            ### m = minus, p = plus
                            point_min_min = (xmin, ymin) # lower left
                            point_max_min = (xmax, ymin) # lower right
                            point_max_max = (xmax, ymax) # upper right
                            point_min_max = (xmin, ymax) # upper left
                            
                            box_points = [point_min_min, point_max_min, point_max_max, point_min_max]
                            box_Points = [shapely.Point(point) for point in box_points]
                            box_Poly   = shapely.Polygon(box_Points)
                            
                            leaflet["subleaflets"][(xi, yi)] = {
                                "xmin": xmin,
                                "xmax": xmax,
                                "ymin": ymin,
                                "ymax": ymax,
                                "box_points": box_points,
                                "box_poly":   box_Poly,
                            }
            
            subleaflet_poly_maker_toc = time.time()
            subleaflet_poly_maker_time = round(subleaflet_poly_maker_toc - subleaflet_poly_maker_tic, 4)
            self.print_term("------------------------------ LEAFLET BOUNDARY BOXES CREATED", "(time spent: "+str(subleaflet_poly_maker_time)+" [s])", "\n", spaces=0)
            
    def holed_subleaflet_bbox_maker(self):
        if len(self.MEMBRANES) != 0:
            holed_subleaflet_bbox_maker_tic = time.time()
            self.print_term("------------------------------ CREATING HOLED BOUNDARY BOXES", spaces=0)
            for memb_i, (memb_key, memb_dict) in enumerate(self.MEMBRANES.items()):
                if memb_i != 0:
                    self.print_term("")
                self.print_term("Starting membrane nr", memb_key, spaces=0)
                
                ### Defnining some default values for union, intersection and holed_bbox
                for leaflet_i, (leaflet_key, leaflet) in enumerate(memb_dict["leaflets"].items()):
                    leaflet["protein_poly"] = False
                    leaflet["union"] = False
                    for (slxi, slyi), subleaflet in leaflet["subleaflets"].items():
                        subleaflet["intersection"] = False
                        subleaflet["holed_bbox"]   = subleaflet["box_poly"]
                
                ######################################
                ### HANDLING PROTEIN RELATED HOLES ###
                ######################################
                ### Getting all the protein bead positions
                if len(self.PROTEINS) != 0:
                    self.print_term("Finding protein beads inside leaflets", spaces=1)

                    for leaflet_i, (leaflet_key, leaflet) in enumerate(memb_dict["leaflets"].items()):
                        self.print_term("Starting", leaflet["leaflet_type"], "leaflet", spaces=2)
                        polygon_list = []

                        ####################################################
                        ### Finding protein points contained in leaflets ###
                        ####################################################
                        prod_beads_in_memb = [
                            (beadx, beady)
                            for protein_i, (protein_nr, protein) in enumerate(self.PROTEINS.items())
                            for beadx, beady, beadz in protein["protein"].get_beads("xyz")
                            if (
                                leaflet["center"][2] - protein["buffer"] < beadz < leaflet["center"][2] + leaflet["lipid_dimensions"]["lipid_height"] + protein["buffer"]
                                and leaflet["HG_direction"] == "up"
                            )
                            or (
                                leaflet["center"][2] + protein["buffer"] > beadz > leaflet["center"][2] - leaflet["lipid_dimensions"]["lipid_height"] - protein["buffer"]
                                and leaflet["HG_direction"] == "down"
                            )
                        ]

                        if len(prod_beads_in_memb) == 0:
                            self.print_term("No leaflet-protein overlap found", spaces=3)
                            leaflet["prot_points"] = False
                        else:
                            self.print_term("Leaflet-protein overlap found for", str(len(prod_beads_in_memb)), "protein beads", spaces=3)
                            leaflet["prot_points"] = prod_beads_in_memb
                            #################################
                            ### CONCAVE HULL / ALPHASHAPE ###
                            #################################
                            
                            alpha = 1 / (leaflet["lipid_dimensions"]["lipid_radius"]*2 * leaflet["alpha_mult"])
                            ### alpha = radius^-1
                            ALPHASHAPE = alphashape.alphashape(points = prod_beads_in_memb, alpha = alpha)
                            self.print_term(type(ALPHASHAPE), debug = True)
                            
                            leaflet["alphashape_1"] = ALPHASHAPE
                            
                            ### Alphashape output can be a bit unpredictable so need to check all possibilities
                            if ALPHASHAPE.geom_type == "Polygon":
                                self.print_term("Polygon", debug = True)
                                ConcaveHulls_Polygon = [ALPHASHAPE]
                            elif ALPHASHAPE.geom_type == "MultiPolygon":
                                self.print_term("MultiPolygon", debug = True)
                                ConcaveHulls_Polygon = list(ALPHASHAPE.geoms)
                            elif ConcaveHulls_Polygon[-1].geom_type == "MultiPolygon":
                                self.print_term("list of MultiPolygon", debug = True)
                                ConcaveHulls_Polygon = list(ALPHASHAPE[-1].geoms)
                                
                            leaflet["ConcaveHulls_Polygon_1"] = ConcaveHulls_Polygon

                            ConcaveHulls_Polygon_Buffered = []
                            ### Iterate over polygons
                            for poly in ConcaveHulls_Polygon:
                                ### Buffering to add a bit of space around protein concave hull
                                poly_buffered = poly.buffer(distance = leaflet["prot_buffer"]) # [Å]
                                ConcaveHulls_Polygon_Buffered.append(poly_buffered)
                            
                            leaflet["protein_poly"] = ConcaveHulls_Polygon_Buffered
                            
                            ############################################
                            ### SECOND ALPHASHAPE TO REMOVE CREVICES ###
                            ############################################
                            buffered_poly_points = []
                            for poly in ConcaveHulls_Polygon_Buffered:
                                try:
                                    buffered_poly_points.extend(list(poly.coords))
                                except:
                                    buffered_poly_points.extend(list(poly.exterior.coords))
                            
                            ALPHASHAPE = alphashape.alphashape(points = buffered_poly_points, alpha = alpha/2)
                            
                            leaflet["alphashape_2"] = ALPHASHAPE
                            
                            ### Alphashape output can be a bit unpredictable so need to check all possibilities
                            if ALPHASHAPE.geom_type == "Polygon":
                                self.print_term("Polygon", debug = True)
                                New_ConcaveHulls_Polygon = [ALPHASHAPE]
                            elif ALPHASHAPE.geom_type == "MultiPolygon":
                                self.print_term("MultiPolygon", debug = True)
                                New_ConcaveHulls_Polygon = list(ALPHASHAPE.geoms)
                            elif ConcaveHulls_Polygon[-1].geom_type == "MultiPolygon":
                                self.print_term("list of MultiPolygon", debug = True)
                                New_ConcaveHulls_Polygon = list(ALPHASHAPE[-1].geoms)
                            
                            leaflet["ConcaveHulls_Polygon_2"] = New_ConcaveHulls_Polygon
                            
                            leaflet["protein_poly"] = ConcaveHulls_Polygon_Buffered + New_ConcaveHulls_Polygon
                
                ####################################################
                ### COMBINING SHAPELY SHAPES FOR ALL SUBLEAFLETS ###
                ####################################################
                self.print_term("")
                self.print_term("Calculating holed boundary box", spaces=1)
                for leaflet_i, (leaflet_key, leaflet) in enumerate(memb_dict["leaflets"].items()):
                    self.print_term("Starting", leaflet["leaflet_type"], "leaflet", spaces=2)
                    
                    ### All the various polygons to be removed from the bbox
                    ### Unique to each leaflet but the same for each subleaflet
                    ### To be expanded later with manually defined holes
                    protein_poly    = leaflet["protein_poly"]
                    poly_check_list = []

                    for val in [protein_poly]:
                        if val:
                            poly_check_list.extend(val)
                    
                    ### Combining polygons into single shape
                    union = shapely.unary_union(poly_check_list)
                    leaflet["union"] = union
                    
                    for (slxi, slyi), subleaflet in leaflet["subleaflets"].items():
                        ### The subleaflet box bbox polygon
                        box_poly     = subleaflet["box_poly"]
                        intersection = shapely.intersection(union, box_poly)
                        holed_bbox   = shapely.difference(box_poly, intersection)
                        
                        subleaflet["intersection"] = intersection
                        subleaflet["holed_bbox"]   = holed_bbox
            
            holed_subleaflet_bbox_maker_toc = time.time()
            holed_subleaflet_bbox_maker_time = round(holed_subleaflet_bbox_maker_toc - holed_subleaflet_bbox_maker_tic, 4)
            self.print_term("------------------------------ HOLED BOUNDARY BOXES CREATED", "(time spent: "+str(holed_subleaflet_bbox_maker_time)+" [s])", "\n", spaces=0)
    
    ########################
    ### LIPID CALCULATOR ###
    ########################
    def lipid_calculator(self):
        self.SYS_lipids_dict = {}
        if len(self.MEMBRANES) != 0:
            lipid_calculator_tic = time.time()
            self.print_term("------------------------------ CALCULATING LIPID RATIOS", spaces=0)
            printer_spacing = 1
            for memb_i, (memb_key, memb_dict) in enumerate(self.MEMBRANES.items()):
                if memb_i != 0:
                    self.print_term("")
                self.print_term("Starting membrane nr", memb_key, spaces=0)
                
                ### Checks if any leaflets contain multiple subleaflets for printer spacing
                for leaflet_i, (leaflet_key, leaflet) in enumerate(memb_dict["leaflets"].items()):
                    if len(leaflet["subleaflets"]) > 1:
                        printer_spacing = 2
                
                for leaflet_i, (leaflet_key, leaflet) in enumerate(memb_dict["leaflets"].items()):
                    if leaflet_i != 0:
                        self.print_term("")
                    self.print_term("Starting", leaflet["leaflet_type"], "leaflet", spaces=1)
                    
                    ### Lipid optimization printing
                    ### Printing here so it isn't spammed if multiple subleafs
                    self.print_term("Lipid optimization 'lipid_optim' setting:", leaflet["lipid_optim"], spaces=1+printer_spacing)
                    if leaflet["lipid_optim"] in ["force_fill", "fill"]:
                        if leaflet["lipid_optim"] == "fill":
                            self.print_term("Filling leaflet until perfect ratio between lipids is achieved or leaflet is full", spaces=2+printer_spacing)
                        elif leaflet["lipid_optim"] == "force_fill":
                            self.print_term("Forcefully filling leaflet until all grid points have a lipid", spaces=2+printer_spacing)
                    elif leaflet["lipid_optim"] == "avg_optimal":
                        self.print_term("Optimizing based on the average deviation from expected ratios", spaces=2+printer_spacing)

                    apl_sqrt = math.sqrt(leaflet["apl"])
                    
                    leaflet["leaf_lipid_count_dict"] = {}
                    leaflet["lipid_names"] = [lipid_name for lipid_name in leaflet["lipids"].keys()]
                    leaflet["lipid_ratios"] = [0 for _ in leaflet["lipid_names"]]
                    
                    for subleaflet_i, ((slxi, slyi), subleaflet) in enumerate(leaflet["subleaflets"].items()):
                        if len(leaflet["subleaflets"]) > 1:
                            if self.subleaflet_extra_info:
                                if subleaflet_i != 0:
                                    self.print_term("")
                                self.print_term("Starting subleaflet", str(subleaflet_i+1)+"/"+str(len(leaflet["subleaflets"])), spaces=2)
                            printer_leaf_name = "Subleaflet"
                        else:
                            printer_leaf_name = "Leaflet"
                        
                        ### Initial area per lipid calculations and max potential lipids in area
                        subleaflet_area = subleaflet["holed_bbox"].area
                        max_lipids_possible_decimal = subleaflet_area / leaflet["apl"]

                        ### Rounding max number of possible lipids according to requested rounding method
                        if leaflet["lip_round_func"][1] == int: # int rounding
                            max_lipids_possible = int(max_lipids_possible_decimal + 0.5)
                            func_for_printer = str(round(max_lipids_possible_decimal, 3)) + "+0.5"
                        else: # round(), math.floor() or math.ceil() rounding
                            max_lipids_possible = leaflet["lip_round_func"][1](max_lipids_possible_decimal)
                            func_for_printer = str(round(max_lipids_possible_decimal, 3))

                        ### Initial rounded estimations
                        lipids = [(lipid_name, lipid_vals.ratio) for lipid_name, lipid_vals in leaflet["lipids"].items()]
                        lipid_names = [name for name, ratio in lipids]
                        lipid_ratios_decimal = [ratio for name, ratio in lipids]
                        lipids_tot = sum(lipid_ratios_decimal)
                        lipid_ratios = [round(int(i / lipids_tot * max_lipids_possible), 3) for i in lipid_ratios_decimal]

                        original_ratios_decimal = [round(i / sum(lipid_ratios) * 100, 3) for i in lipid_ratios]
                        original_ratios = lipid_ratios[:]

                        expected_ratios_decimal = [round(ratio / sum(lipid_ratios_decimal) * 100, 3) for ratio in lipid_ratios_decimal]

                        

                        ### Just take integer lipid values and remove excess grid points
                        if leaflet["lipid_optim"] == "abs_val":
                            lipid_ratios = [ratio for name, ratio in lipids]
                            ### Error out if values are not integers
                            assert all([type(val) == int for val in lipid_ratios]), "Not all supplied lipid values are integers. You can't have half a lipid, no matter how much you want to."
                            ### Error out if more lipids than available grid points
                            assert sum(lipid_ratios) <= max_lipids_possible, "You have specified too many lipids. Sorry they are too fat and won't fit in, try reducing number of lipids or the APL."

                        elif leaflet["lipid_optim"] in ["force_fill", "fill", "avg_optimal"]:
                            ###############################
                            ### OPTIMIZING LIPID RATIOS ###
                            ###############################
                            iters = [lipid_ratios[:]]

                            while sum(lipid_ratios) < max_lipids_possible:
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
#                                 if leaflet["lipid_optim"] == "fill":
#                                     self.print_term("Filling leaflet until perfect ratio between lipids is achieved or leaflet is full", spaces=2+printer_spacing)
#                                 elif leaflet["lipid_optim"] == "force_fill":
#                                     self.print_term("Forcefully filling leaflet until all grid points have a lipid", spaces=2+printer_spacing)
                                lipid_ratios = iters[-1]

                            ### Optimizing lipids
                            elif leaflet["lipid_optim"] == "avg_optimal": # lipid_ratios_decimal
#                                 self.print_term("Optimizing based on the average deviation from expected ratios", spaces=2+printer_spacing)

                                iters_decimal = [[round(iteration[i] / sum(iteration) * 100, 3) for i in range(len(iteration))] for iteration in iters]
                                iters_abs_diff = [
                                    [abs(iteration[i] - expected_ratios_decimal[i]) for i in range(len(iteration))]
                                    for iteration in iters_decimal
                                ]

                                iters_avg = [(np.mean(iters), i) for i, iters in enumerate(iters_abs_diff)]
                                iters_avg_optimal = min(iters_avg, key=lambda i: i[0])
                                lipid_ratios = iters[iters_avg_optimal[1]]

                        subleaflet["lipid_names"]  = lipid_names
                        subleaflet["lipid_ratios"] = lipid_ratios
                        
                        for name, count in zip(subleaflet["lipid_names"], subleaflet["lipid_ratios"]):
                            if name not in leaflet["leaf_lipid_count_dict"]:
                                leaflet["leaf_lipid_count_dict"][name] = count
                            else:
                                leaflet["leaf_lipid_count_dict"][name] += count
                        
                        ##################################################
                        ### PRINTING SUBLEAFLET SPECIFIC LIPID DETAILS ###
                        ##################################################
                        if self.extra_info:
                            ### Printing mean, min and max deviation from wanted ratios
                            ratios_final_decimal = [round(lipid_ratios[i] / sum(lipid_ratios) * 100, 3) for i in range(len(lipid_ratios))]
                            ratios_final_abs_diff = [abs(ratios_final_decimal[i] - expected_ratios_decimal[i]) for i in range(len(ratios_final_decimal))]
                            ratios_final_avg = round(np.mean(ratios_final_abs_diff), 5) # Don't multiply by 100 as they are already percentages
                            ratios_final_min = round(min(ratios_final_abs_diff), 5) # Don't multiply by 100 as they are already percentages
                            ratios_final_max = round(max(ratios_final_abs_diff), 5) # Don't multiply by 100 as they are already percentages

                            final_ratios_decimal = [round(ratio / sum(lipid_ratios) * 100, 3) for ratio in lipid_ratios]
                            final_lipid_data = [lipid_names, lipid_ratios, final_ratios_decimal, expected_ratios_decimal, original_ratios_decimal, original_ratios]

                            headers = ["Lipid name", "Final lipids", "Final %", "Expected %", "Starting %", "Starting lipids"]
                            lipids_for_printer = [[head] + data for head, data in zip(headers, final_lipid_data)]
                            max_lengths = [len(str(max(data, key = lambda d: len(str(d))))) for data in lipids_for_printer]
                            
                            for L0, L1, L2, L3, L4, L5 in list(zip(*lipids_for_printer))[1:]:
                                if L0 not in self.SYS_lipids_dict.keys():
                                    self.SYS_lipids_dict[L0] = L1
                                else:
                                    self.SYS_lipids_dict[L0] += L1
                            
                            if printer_spacing == 1 or (printer_spacing == 2 and self.subleaflet_extra_info):
                                self.print_term("Final deviations from expected ratios:", "Difference in %-values", spaces=1+printer_spacing)
                                self.print_term("Maximum:", ratios_final_max, spaces=2+printer_spacing)
                                self.print_term("Average:", ratios_final_avg, spaces=2+printer_spacing)
                                self.print_term("Minimum:", ratios_final_min, spaces=2+printer_spacing)

                                self.print_term(
                                    printer_leaf_name, "specific lipid data",
                                    "(Max lipids: "+str(max_lipids_possible)+", Final lipids: "+str(sum(lipid_ratios))+")",
                                    spaces=1+printer_spacing
                                    
                                )
                                for i, (L0, L1, L2, L3, L4, L5) in enumerate(list(zip(*lipids_for_printer))):
                                    self.print_term(
                                        '{0: <{L}}'.format(L0, L = max_lengths[0]), ":",
                                        '{0: <{L}}'.format(L1, L = max_lengths[1]), ":",
                                        '{0: <{L}}'.format(L2, L = max_lengths[2]), ":",
                                        '{0: <{L}}'.format(L3, L = max_lengths[3]), ":",
                                        '{0: <{L}}'.format(L4, L = max_lengths[4]), ":",
                                        '{0: <{L}}'.format(L5, L = max_lengths[5]),
                                        spaces=2+printer_spacing
                                    )
                    
                    leaflet["leaf_lipid_count"] = [(n, c) for n, c in leaflet["leaf_lipid_count_dict"].items()]

                    ###############################################
                    ### PRINTING LEAFLET SPECIFIC LIPID DETAILS ###
                    ###############################################
                    if self.extra_info and len(leaflet["subleaflets"]) > 1:
                        self.print_term()
                        LEAF_headers = ["Lipid name", "Total lipids", "Total %"]
                        LEAF_lipid_names, LEAF_lipid_vals = zip(*leaflet["leaf_lipid_count"])
                        LEAF_tot_lipids = sum(LEAF_lipid_vals)
                        LEAF_lipid_percentages = [round(val / LEAF_tot_lipids * 100, 3) for val in LEAF_lipid_vals]

                        LEAF_lipids_for_printer = list(zip(LEAF_lipid_names, LEAF_lipid_vals, LEAF_lipid_percentages))
                        LEAF_lipids_for_printer = [tuple(LEAF_headers)] + LEAF_lipids_for_printer
                        LEAF_max_lengths = [len(str(max(data, key = lambda d: len(str(d))))) for data in zip(*LEAF_lipids_for_printer)]

                        self.print_term("Leaflet specific lipid data (Combined subleafs)", spaces=2)
                        for L0, L1, L2 in LEAF_lipids_for_printer:
                            self.print_term(
                                '{0: <{L}}'.format(L0, L = LEAF_max_lengths[0]), ":",
                                '{0: <{L}}'.format(L1, L = LEAF_max_lengths[1]), ":",
                                '{0: <{L}}'.format(L2, L = LEAF_max_lengths[2]),
                                spaces=2+printer_spacing
                            )

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

                self.print_term("\nLipid data for whole system", spaces=0)
                for L0, L1, L2 in SYS_lipids_for_printer:
                    self.print_term(
                        '{0: <{L}}'.format(L0, L = SYS_max_lengths[0]), ":",
                        '{0: <{L}}'.format(L1, L = SYS_max_lengths[1]), ":",
                        '{0: <{L}}'.format(L2, L = SYS_max_lengths[2]),
                        spaces=2+printer_spacing
                    )
            
            lipid_calculator_toc = time.time()
            lipid_calculator_time = round(lipid_calculator_toc - lipid_calculator_tic, 4)
            self.print_term("------------------------------ LIPID RATIO CALCULATIONS COMPLETE", "(time spent: "+str(lipid_calculator_time)+" [s])", "\n", spaces=0)
    
    #######################################
    ### PLANAR GRID MAKER AND OPTIMIZER ###
    #######################################
    def planar_grid_maker(self):
        self.SYS_lipids_dict = {}
        self.GRID_PLOTTING = {}
        if len(self.MEMBRANES) != 0:
            grid_making_tic = time.time()
            self.print_term("------------------------------ CREATING LIPID GRID", spaces=0)
            for memb_i, (memb_key, memb_dict) in enumerate(self.MEMBRANES.items()):
                if memb_i != 0:
                    self.print_term("")
                self.print_term("Starting membrane nr", memb_key, spaces=0)
                for leaflet_i, (leaflet_key, leaflet) in enumerate(memb_dict["leaflets"].items()):
                    if leaflet_i != 0:
                        self.print_term("")
                    self.print_term("Starting", leaflet["leaflet_type"], "leaflet", spaces=1)
                    
                    for sli, ((slxi, slyi), subleaflet) in enumerate(leaflet["subleaflets"].items()):
                        
                        lipid_names_nlipids_radii = [
                            (
                                name,
                                ratio,
                                leaflet["lipids"][name].get_radius() + max([leaflet["kickx"], leaflet["kicky"]]) + leaflet["plane_buffer"]
                            )
                            for name, ratio in zip(subleaflet["lipid_names"], subleaflet["lipid_ratios"])
                        ]
                        lipids = []
                        
                        leaflet_area = subleaflet["holed_bbox"].area
                        lipids_circle_area = sum([(math.pi*(radius**2))*ratio for name, ratio, radius in lipid_names_nlipids_radii])
                        lipids_square_area = sum([((radius*2)**2)*ratio for name, ratio, radius in lipid_names_nlipids_radii])
                        lipids_mean_area = (lipids_circle_area + lipids_square_area)/2
                        occupation_modifier = (leaflet_area-lipids_square_area)/(leaflet_area*2) # *2 to half the modifier
    
                        self.print_term("leaflet_area       ", leaflet_area, debug=True)
                        self.print_term("lipids_circle_area ", lipids_circle_area, debug=True)
                        self.print_term("lipids_square_area ", lipids_square_area, debug=True)
                        self.print_term("lipids_mean_area   ", lipids_mean_area, debug=True)
                        self.print_term("occupation_modifier", occupation_modifier, debug=True)
                        
                        if occupation_modifier < 0.1:
                            self.print_term(
                                "WARNING:",
                                "Chosen apl ("+str(leaflet["apl"]/100)+" [nm^2]) and buffer ("+str(leaflet["plane_buffer"]/10)+" [nm]) cause lipids to be very closely packed.",
                                "Optimization may be slow and some lipid overlaps may not be avoidable.",
                                "Consider increasing the apl using the subcommand 'apl' or decreasing the plane buffer using the subcommand 'plane_buffer'",
                                
                                warn=True
                            )
                        
                        
                        for name, nlipids, radius in lipid_names_nlipids_radii:
                            lipid = []
                            for i in range(nlipids):
                                lipid.append((name, radius))

                            lipids.append(lipid)
                        
                        if leaflet["lipid_distribution"] == "evenly":
                            lipids = self.n_list_mixer(*lipids)
                        elif leaflet["lipid_distribution"] == "random":
                            lipids = flatten(lipids)
                            random.shuffle(lipids)
                        
                        subleaflet["lipids"] = lipids
                        
                        if leaflet["HG_direction"] == "up":
                            sign = +1
                        if leaflet["HG_direction"] == "down":
                            sign = -1
                        
                        grid_points, grid_points_no_random = self.make_rect_grid(leaflet, subleaflet)
                        grid_points_arr = np.asarray(grid_points)
                        
                        grid_z_value = leaflet["center"][2] + (sign * leaflet["lipid_dimensions"]["zUnderLipCen"])
                        z_values_arr = np.asarray([[grid_z_value] for _ in range(len(grid_points_arr))])
                        
                        if leaflet["optimize"]:
                            if len(leaflet["subleaflets"]) > 1 and self.subleaflet_extra_info:
                                self.print_term("Starting optimization for subleaflet nr", sli+1, spaces=2)
                            else:
                                self.print_term("Starting optimization", spaces=2)
                            
                            optimized_grid_points_arr, POINT_STEPS, step, steps_time, mean_steps_time, max_push, bin_size = self.plane_grid_point_optimizer(
                                grid_points_arr,
                                lipid_sizes = list(zip(*lipids))[-1],
                                polygon = subleaflet["holed_bbox"],
                                
                                xcenter = np.mean([subleaflet["xmax"], subleaflet["xmin"]]),
                                ycenter = np.mean([subleaflet["ymax"], subleaflet["ymin"]]),
                                
                                xlen = subleaflet["xmax"] - subleaflet["xmin"],
                                ylen = subleaflet["ymax"] - subleaflet["ymin"],
                                
                                maxsteps = leaflet["optim_maxsteps"],
                                push_tolerance = leaflet["optim_push_tol"],
                                push_mult = leaflet["optim_push_mult"],
                                buffer = leaflet["plane_buffer"],
                                
                                occupation_modifier = occupation_modifier,
                                optimize = leaflet["optimize"],
                            )
                            subleaflet["grid_points"] = np.hstack([optimized_grid_points_arr, z_values_arr])
                            
                            if self.plot_grid:
                                self.GRID_PLOTTING[(memb_key, leaflet_key, slxi, slyi)] = {
                                    ### Inputs
                                    "lipids"      : lipids,
                                    "bbox_polygon": subleaflet["holed_bbox"],
                                    "prot_points" : leaflet["prot_points"],
                                    
                                    "alphashape_1" : leaflet["alphashape_1"],
                                    "ConcaveHulls_Polygon_1" : leaflet["ConcaveHulls_Polygon_1"],
                                    "protein_poly" : leaflet["protein_poly"],
                                    "alphashape_2" : leaflet["alphashape_2"],
                                    "ConcaveHulls_Polygon_2" : leaflet["ConcaveHulls_Polygon_2"],
                                    "protein_poly" : leaflet["protein_poly"],
                                    
                                    "xdims"       : (subleaflet["xmin"], subleaflet["xmax"]),
                                    "ydims"       : (subleaflet["ymin"], subleaflet["ymax"]),
                                    "push_mult"   : leaflet["optim_push_mult"],
                                    "apl"         : leaflet["apl"],
                                    "bin_size"    : bin_size,

                                    ### Outputs
                                    "grid_points"          : grid_points_arr,
                                    "grid_points_no_random": grid_points_no_random,
                                    "points_arr"           : optimized_grid_points_arr,
                                    "POINT_STEPS"          : POINT_STEPS,
                                    "step"                 : step,
                                    "steps_time"           : steps_time,
                                    "mean_steps_time"      : mean_steps_time,
                                    "max_push"             : max_push,
                                }
                        else:
                            subleaflet["grid_points"] = np.hstack([grid_points_arr, z_values_arr])
            grid_making_toc = time.time()
            grid_making_time = round(grid_making_toc - grid_making_tic, 4)
            self.print_term("------------------------------ LIPID GRID CREATED", "(time spent: "+str(grid_making_time)+" [s])", "\n", spaces=0)
    
    def make_rect_grid(self, leaflet, subleaflet):
        xmin, xmax, ymin, ymax = itemgetter("xmin", "xmax", "ymin", "ymax")(subleaflet)
        
        sidelen      = math.sqrt(leaflet["apl"])
        edge_buffer  = leaflet["lipid_dimensions"]["lipid_radius"]*1.1
        bbox_polygon = subleaflet["holed_bbox"]
        lipids       = subleaflet["lipids"]
        
        xmin_edge = xmin + edge_buffer
        xmax_edge = xmax - edge_buffer
        ymin_edge = ymin + edge_buffer
        ymax_edge = ymax - edge_buffer
        
        ### Making pointer for when number of points along x/y-axis should be expanded
        xpoints_ratio = (xmax_edge-xmin_edge)/sidelen
        ypoints_ratio = (ymax_edge-ymin_edge)/sidelen
        xsum = round(xpoints_ratio/sum([xpoints_ratio, ypoints_ratio])*1000)
        ysum = round(ypoints_ratio/sum([xpoints_ratio, ypoints_ratio])*1000)
        xy_pointer = self.n_list_mixer(["x"]*xsum, ["y"]*ysum)
        
        xpoints = 0
        ypoints = 0
        xy_pointer_counter = 0
        while xpoints * ypoints < len(lipids):
            if xy_pointer[xy_pointer_counter] == "x":
                xpoints += 1
            elif xy_pointer[xy_pointer_counter] == "y":
                ypoints += 1
            xy_pointer_counter += 1
            if xy_pointer_counter == 1000:
                xy_pointer_counter = 0
        
        enough_points = False
        while enough_points == False:
            xlinespace = np.linspace(start=xmin_edge, stop=xmax_edge, num=xpoints, endpoint=True)
            ylinespace = np.linspace(start=ymin_edge, stop=ymax_edge, num=ypoints, endpoint=True)
            np.round(
                np.linspace(start=0, stop=xpoints*ypoints, num=len(lipids), endpoint=False)
            ).astype(int)
            grid_points_list = []

            linepoints = []
            for xval in xlinespace:
                for yval in ylinespace:
                    point = (xval, yval)
                    valid_point = bbox_polygon.contains(shapely.Point(point))

                    if valid_point:
                        linepoints.append(point)
                    elif linepoints != []:
                        grid_points_list.append(linepoints)
                        linepoints = []
                if linepoints != []:
                    grid_points_list.append(linepoints)
                    linepoints = []
            if linepoints != []:
                grid_points_list.append(linepoints)
                linepoints = []

            ngridpoints = sum([len(l) for l in grid_points_list])

            if ngridpoints < len(lipids):
                if xy_pointer[xy_pointer_counter] == "x":
                    xpoints += 1
                elif xy_pointer[xy_pointer_counter] == "y":
                    ypoints += 1
                xy_pointer_counter += 1
                if xy_pointer_counter == 1000:
                    xy_pointer_counter = 0
            else:
                enough_points = True
        
        ### Finding the points that should be removed due to number of lipids
        ### No duplicate index values as "len(lipids)" is always equal to or smaller than "ngridpoints"
        point_is = np.round(
            np.linspace(start=0, stop=ngridpoints, num=len(lipids), endpoint=False)
        ).astype(int)

        ### Removing points due to number of lipids
        cur_i = 0
        grid_points_instructions = []
        for linepoints_old in grid_points_list:
            npoints_in_line = 0
            xvals, yvals = zip(*linepoints_old)
            xval = np.mean(xvals)
            ymin, ymax = min(yvals), max(yvals)

            for point in linepoints_old:
                if cur_i in point_is:
                    npoints_in_line += 1
                cur_i += 1
            if npoints_in_line != 0:
                grid_points_instructions.append((xval, ymin, ymax, npoints_in_line))

        ### Respacing points on a line due to removal of points
        grid_points = []
        grid_points_no_random = []
        rand_force = 10
        for xval, ymin, ymax, npoints_in_line in grid_points_instructions:
            new_yvals = np.linspace(start=ymin, stop=ymax, num=npoints_in_line)
            for yval in new_yvals:
                grid_points.append(
                    (xval+random.uniform(-sidelen/rand_force, sidelen/rand_force),
                     yval+random.uniform(-sidelen/rand_force, sidelen/rand_force))
                )
                grid_points_no_random.append((xval, yval))

#         return grid_points
        return grid_points, grid_points_no_random
    
    def plane_grid_point_optimizer(self, grid_points, lipid_sizes, polygon, xcenter, ycenter, xlen, ylen, maxsteps, push_tolerance, push_mult, buffer, occupation_modifier, optimize):

        def push_func(dist, power, mult):
            return dist**-power*mult

        def get_vector(A, B):
            ### Get vector from point A to point B
            return (round(B[0] - A[0], 3), round(B[1] - A[1], 3))
        
        def get_vector_len_fastsingle(v):
            '''
            Fastest single vector calculator
            https://stackoverflow.com/questions/7370801/how-do-i-measure-elapsed-time-in-python
            '''
            return math.sqrt(v[0]**2 + v[1]**2)

        def rand_sign():
            if random.random() < 0.5:
                return 1
            else:
                return -1

        def update_neighborlist(bins_arr):
            neighborlist = np.zeros((points_arr.shape[0], points_arr.shape[0]), dtype=int)
            nst_tic = time.time()
            for bix, binx in enumerate(bins_arr):
                for biy, biny in enumerate(binx):
                    points_for_nstlist = biny
                    for pj1, pi1 in enumerate(points_for_nstlist):
                        for pj2, pi2 in enumerate(points_for_nstlist[pj1+1:], pj1+1):
                            point1 = points_arr[pi1]
                            point2 = points_arr[pi2]
                            dist = scipy.spatial.distance.cdist([point1], [point2])[0]
#                             if dist < bin_size*3:
                            if dist < largest_lipid*4:
                                neighborlist[pi1, pi2] = 1
            nst_toc = time.time()
            self.print_term("    ", "    ", "    ", "    ", "Time spent calculating neighborlist:", round(nst_toc-nst_tic, 4), debug=True)

            return neighborlist

        def coord_to_indices(pos, dim, center, real_gridres):
            posi_dec = (pos+dim/2-center)/real_gridres
            posi_int = int(posi_dec)
            if posi_int < 0:
                posi_int = 0
            return posi_int

        def binning_func():
            ### Creating bins array
            new_bins_arr = np.empty((xnbins, ynbins), dtype=object)
            ### np.fill([]) makes all array indeces point to the same list so can't be used
            ### For loop is very fast so doesn't matter for time
            bin_tic = time.time()
            for xi in range(xnbins):
                for yi in range(ynbins):
                    new_bins_arr[xi][yi] = []

#             indeces = [-2, -1, 0, 1, 2]
            indeces = [-1, 0, 1]
            for pi, (px, py) in enumerate(points_arr):
                pix = coord_to_indices(px, xlen, xcenter, xbinlen)
                piy = coord_to_indices(py, ylen, ycenter, ybinlen)
                ### Surrounding bins
                for xi in indeces:
                    for yi in indeces:
                        if xnbins > pix+xi >= 0 and ynbins > piy+yi >= 0:
                            new_bins_arr[pix+xi][piy+yi].append(pi)
            bin_toc = time.time()
            self.print_term("    ", "    ", "    ", "    ", "Time spent calculating bins:        ", round(bin_toc-bin_tic, 4), debug=True)

            return new_bins_arr
        
        points_arr = np.array(grid_points)
        
        largest_lipid  = max(lipid_sizes)
        smallest_lipid = min(lipid_sizes)
        bin_size       = largest_lipid*2
        
        self.print_term("len(grid_points)", len(grid_points), debug=True)
        self.print_term("xlen            ", xlen,             debug=True)
        self.print_term("ylen            ", ylen,             debug=True)
        self.print_term("bin_size        ", bin_size,         debug=True)
        
        ### Calculating number of bins along each axis, by rounding down
        xnbins = int(xlen / bin_size)
        ynbins = int(ylen / bin_size)
        ### Finding the actual lengths of the bins
        xbinlen = xlen / xnbins
        ybinlen = ylen / ynbins
        
        dists_traveled = np.zeros((points_arr.shape[0],))
        
        if self.plot_grid:
            POINT_STEPS = []
            POINT_STEPS.append(points_arr.copy())
        else:
            POINT_STEPS = False

        start = 1
        end   = maxsteps

        steps_tic = time.time()
        
        ### Making initial bins and neighborlists
        self.print_term("STEP:", 0, debug=True, spaces=2)
        self.print_term("Updating bins and neigborlist", debug=True, spaces=3)
        dists_traveled = np.zeros((points_arr.shape[0],))
        bins_arr       = binning_func()
        neighborlist   = update_neighborlist(bins_arr)
        
        step_modifier_limit = 15
        
        bounce_counter = np.zeros((len(points_arr)))
        max_push = 0
        
        self.print_term("CURRENT STEP:", end=" ", spaces=3)
        
        for si, step in enumerate(np.arange(1, maxsteps+1)):
            '''CODE LEGEND

            'dists_traveled':
                Numpy array with shape (nlipids).
                Remembers how much all individual lipids have moved.
                Forces a bin and neighborlist update if largest movement is larger than the bin size.
                Resets to zero-array upon neighborlist updates.

            'bins_arr':
                Numpy array with shape (xbins, ybins).
                Each index is a list of points that are present within the bin or in the neighbouring bins.
                Used for neighborlist updates.

            'bin_pointers': Currently unused
                List with length (nlipids) with each index being a tuple containing (xpointer, ypointer).
                Contains the bin x/y indeces that points are in.

            'neighborlist':
                Numpy array with shape (nlipids, nlipids).
                Contains the neighbor indeces for each grid point.

            'points_arr':
                Numpy array with shape (nlipids, 2).
                Contains the x/y coordinates of all grid points.

            'push_arr':
                Numpy array with shape (nlipids, 2).
                Contains the combined vector push applied during a single time step.
                Modifies 'points_arr' at the end of a time step and then resets to only containing zeros.

            'pushes':
                Numpy array with shape (nlipids).
            '''
            
            if step >= step_modifier_limit:
                step_modifier = step_modifier_limit**-1
            else:
                step_modifier = 1-(step/step_modifier_limit)
            
            self.print_term(step, end=" ", spaces=0)
            
            push_arr = np.empty((len(points_arr)), dtype=object)
            for i in range(len(push_arr)):
                push_arr[i] = []

            pushes   = np.zeros((len(points_arr)))
            max_dists_traveled = max(dists_traveled)
            
            self.print_term("STEP:", step, "AT TIME:", round(time.time() - steps_tic, 4), debug=True, spaces=2)
            self.print_term("Max push:         ", round(max_push, 4), debug=True, spaces=3)
            self.print_term("Max dist traveled:", round(max_dists_traveled, 4), debug=True, spaces=3)
            
            if max_dists_traveled >= bin_size/4:
                self.print_term("Updating bins and neigborlist", debug=True, spaces=3)
                dists_traveled = np.zeros((points_arr.shape[0],))
                bins_arr       = binning_func()
                neighborlist   = update_neighborlist(bins_arr)
            
            for pi1, point1 in enumerate(points_arr):
                neighbor_is = np.nonzero(neighborlist[pi1])[0]
                for pi2 in neighbor_is:
                    if pi2 <= pi1:
                        continue
                    point2 = points_arr[pi2]
                    dist = scipy.spatial.distance.cdist([point1], [point2])[0][0]
                    combined_lipid_size = lipid_sizes[pi1] + lipid_sizes[pi2]
                    
                    if optimize in ["v1", "limited"] and dist < combined_lipid_size:
                        vector = np.array(get_vector(point1, point2))
                        if dist < combined_lipid_size/2:
                            push = combined_lipid_size-dist
                        else:
                            push = push_func(dist, power=2, mult=push_mult)

                        vector_push = vector*push
                        vector_push_len = get_vector_len_fastsingle(vector_push)
                        
                        ### Pseudo-normalizes the force behind the pushes
                        ### Modifies pushes according to the radius of each lipid
                        ### Smaller radius causes more of the push to go to the lipid
                        max_size = max([lipid_sizes[pi2], lipid_sizes[pi1]])
                        pi1_mult = 1/(max_size/lipid_sizes[pi2])
                        pi2_mult = 1/(max_size/lipid_sizes[pi1])

#                         print(lipid_size_sum, max_size, pi1_mult, pi2_mult)
#                         print()

                        push_arr[pi1] -= vector_push*pi1_mult
                        push_arr[pi2] += vector_push*pi2_mult

                        pushes[pi1] += vector_push_len*pi1_mult #abs(push*pi1_mult)
                        pushes[pi2] += vector_push_len*pi2_mult #abs(push*pi2_mult)
                        
                        dists_traveled[pi1] += vector_push_len*pi1_mult # abs(push*pi1_mult)
                        dists_traveled[pi2] += vector_push_len*pi2_mult # abs(push*pi2_mult)
                    
                    elif optimize == "v2":
                    
                        ideal_dist = combined_lipid_size*(1+occupation_modifier)
                        ideal_dist_diff = ideal_dist - combined_lipid_size
                        ideal_dist_upper = ideal_dist + ideal_dist_diff
                    
                        if dist < ideal_dist_upper:
                            vector = np.array(get_vector(point1, point2))

                            ### Difference in distance from ideal distance 
                            dist_diff = dist-combined_lipid_size
                            if dist_diff <= ideal_dist_diff:
                                dist_modifier = 1
                            else:
                                dist_modifier = abs(ideal_dist*dist_diff)**-1

                            push = (ideal_dist-dist)*occupation_modifier*step_modifier*dist_modifier

                            vector_push = vector*push
                            vector_push_len = get_vector_len_fastsingle(vector_push)

                            ### Pseudo-normalizes the force behind the pushes
                            max_size = max([lipid_sizes[pi2], lipid_sizes[pi1]])
                            pi1_mult = 1/(max_size/lipid_sizes[pi2])
                            pi2_mult = 1/(max_size/lipid_sizes[pi1])

                            push_arr[pi1] -= vector_push*pi1_mult
                            push_arr[pi2] += vector_push*pi2_mult

                            pushes[pi1] += vector_push_len*pi1_mult # abs(push*pi1_mult)
                            pushes[pi2] += vector_push_len*pi2_mult # abs(push*pi2_mult)

                            dists_traveled[pi1] += vector_push_len*pi1_mult # abs(push*pi1_mult)
                            dists_traveled[pi2] += vector_push_len*pi2_mult # abs(push*pi2_mult)
                    
                    elif optimize == "v3":
                    
                        ideal_dist = combined_lipid_size*(1+occupation_modifier)
                        ideal_dist_diff = ideal_dist - combined_lipid_size
                        ideal_dist_lower = combined_lipid_size
                        ideal_dist_upper = ideal_dist + ideal_dist_diff
                        
                        if dist <= ideal_dist_upper:
                            vector = np.array(get_vector(point1, point2))
                            vector_len = get_vector_len_fastsingle(vector)
                            vector /= vector_len

                            ### Difference in distance from ideal distance 
                            dist_diff = dist-combined_lipid_size
                            if dist <= ideal_dist_lower:
                                ### Push very hard if way too close
                                push = (combined_lipid_size-dist)#/2
                                if dist > combined_lipid_size-buffer:
                                    ### Modulated push with the step modifier if within buffer space
                                    push *= step_modifier
                            else:
                                dist_modifier = abs(ideal_dist-dist_diff)**-1
                                push = (ideal_dist-dist)*occupation_modifier*step_modifier*dist_modifier

                            vector_push = vector*push
                            vector_push_len = get_vector_len_fastsingle(vector_push)

                            ### Pseudo-normalizes the force behind the pushes
                            max_size = max([lipid_sizes[pi2], lipid_sizes[pi1]])
                            pi1_mult = 1/(max_size/lipid_sizes[pi2])
                            pi2_mult = 1/(max_size/lipid_sizes[pi1])

                            push_arr[pi1] -= vector_push*pi1_mult
                            push_arr[pi2] += vector_push*pi2_mult

                            pushes[pi1] += vector_push_len*pi1_mult
                            pushes[pi2] += vector_push_len*pi2_mult

                            dists_traveled[pi1] += vector_push_len*pi1_mult
                            dists_traveled[pi2] += vector_push_len*pi2_mult
                            
                    elif optimize == "v4":
                    
                        ideal_dist = combined_lipid_size*(1+occupation_modifier)
                        ideal_dist_diff = ideal_dist - combined_lipid_size
                        ideal_dist_lower = combined_lipid_size
                        ideal_dist_upper = ideal_dist + ideal_dist_diff
                        
                        if dist <= ideal_dist_upper:
                            vector = np.array(get_vector(point1, point2))
                            vector_len = get_vector_len_fastsingle(vector)
                            vector /= vector_len

                            ### Difference in distance from ideal distance 
#                             dist_diff = dist-combined_lipid_size
                            dist_diff = dist-ideal_dist
                            if dist <= ideal_dist_lower:
                                ### Push very hard if way too close
                                push = (combined_lipid_size-dist)#/2
                                if dist > combined_lipid_size-buffer:
                                    ### Modulated push with the step modifier if within buffer space
                                    push *= step_modifier
                            elif step < step_modifier_limit:
#                             elif ideal_dist_lower < dist <= ideal_dist_upper:
                                push = abs(ideal_dist-dist)*step_modifier
    
                                if dist > ideal_dist_upper:
                                    dist_modifier = dist_diff/ideal_dist
                                    print("dist_modifier", dist_modifier)
                                    push *= dist_modifier
                                
#                                 dist_modifier = abs(ideal_dist-dist_diff)**-1
                                
# #                                 dist_modifier = abs(ideal_dist-dist_diff)**-1*2
# #                                 push = (ideal_dist-dist)*occupation_modifier*step_modifier*dist_modifier
#                                 print("step_modifier", round(step_modifier, 3))
#                                 push = (ideal_dist-dist)*step_modifier
#                             else:
#                                 dist_modifier = abs(ideal_dist-dist_diff)**-1
# #                                 push = (ideal_dist-dist)*occupation_modifier*step_modifier*dist_modifier
#                                 print("step_modifier, dist_modifier", step_modifier, round(dist_modifier, 3), round(ideal_dist-dist_diff, 3))
#                                 if dist_modifier > 1:
#                                     print(ideal_dist, dist, dist_diff, abs(ideal_dist-dist_diff))
#                                     print("\n\n\n\n\n\n\n")
#                                 push = (ideal_dist-dist)*step_modifier*dist_modifier

                            vector_push = vector*push
                            vector_push_len = get_vector_len_fastsingle(vector_push)

                            ### Pseudo-normalizes the force behind the pushes
                            max_size = max([lipid_sizes[pi2], lipid_sizes[pi1]])
                            pi1_mult = 1/(max_size/lipid_sizes[pi2])
                            pi2_mult = 1/(max_size/lipid_sizes[pi1])

                            push_arr[pi1].append(-vector_push*pi1_mult)
                            push_arr[pi2].append(+vector_push*pi2_mult)

                    elif optimize == "v5":
                    
                        ideal_dist = combined_lipid_size*(1+occupation_modifier)
                        ideal_dist_diff = ideal_dist - combined_lipid_size
                        ideal_dist_lower = combined_lipid_size
                        ideal_dist_upper = ideal_dist + ideal_dist_diff
                        
                        if dist <= ideal_dist_upper:
                            vector = np.array(get_vector(point1, point2))
                            vector_len = get_vector_len_fastsingle(vector)
                            vector /= vector_len
                            push = 0
                            
                            ### Difference in distance from ideal distance 
                            if dist <= ideal_dist_lower: # combined_lipid_size
                                ### Push very hard if way too close
                                push += (combined_lipid_size-dist)/2

                            if dist <= ideal_dist_upper:
                                ### Always add smaller extra push
                                push += abs(ideal_dist_upper-dist)*step_modifier
                            
                            vector_push = vector*push
                            vector_push_len = get_vector_len_fastsingle(vector_push)

                            ### Pseudo-normalizes the force behind the pushes
                            max_size = max([lipid_sizes[pi2], lipid_sizes[pi1]])
                            pi1_mult = 1/(max_size/lipid_sizes[pi2])
                            pi2_mult = 1/(max_size/lipid_sizes[pi1])

                            push_arr[pi1].append(-vector_push*pi1_mult)
                            push_arr[pi2].append(+vector_push*pi2_mult)

            ### Pushes grid points
            for pi, push in enumerate(push_arr):
                if push:
                    push_vector         = np.mean(np.array(push), axis=0)
                    push_vector_len     = get_vector_len_fastsingle(push_vector)
#                     print(push)
#                     print(np.array(push))
#                     print(np.mean(np.array(push), axis=0))
#                     print(push_vector, push_vector_len)
                    points_arr[pi]     += np.mean(np.array(push), axis=0)
                    pushes[pi]         += push_vector_len
                    dists_traveled[pi] += push_vector_len

            ### BBOX and protein distance checks
            all_contained = True
            for pi1, point1 in enumerate(points_arr):
                if bounce_counter[pi1] > 0:
                    bounce_counter[pi1] -= 0.1
                ### Shapely Point
                point_Point = shapely.Point(point1)
                ### polygon.boundary includes both the outer surface and the surface of holes
                dist = polygon.boundary.distance(point_Point)
                point_contained = polygon.contains(point_Point)

                eq_dist = lipid_sizes[pi1]*(1+occupation_modifier/2*step_modifier)

                if point_contained and dist < eq_dist:
                    poly_nearest, point_Point = shapely.ops.nearest_points(polygon.boundary, point_Point)
                    vector = np.array(get_vector((poly_nearest.x, poly_nearest.y), (point_Point.x, point_Point.y)))
                    
                    diff        = eq_dist - dist
                    push        = diff/dist
                    
                    if dist < lipid_sizes[pi1]:
                        bounce_counter[pi1] += 1
                        push *= bounce_counter[pi1]

                    vector_push = vector*push
                    vector_push_len = get_vector_len_fastsingle(vector_push)
                    
                    points_arr[pi1]     += vector_push
                    pushes[pi1]         += vector_push_len # abs(push)
                    dists_traveled[pi1] += vector_push_len # abs(push)
                    
                elif not point_contained:
                    all_contained = False
                    
                    poly_nearest, point_Point = shapely.ops.nearest_points(polygon.boundary, point_Point)
                    vector = np.array(get_vector((poly_nearest.x, poly_nearest.y), (point_Point.x, point_Point.y)))
                    
                    diff        = -(eq_dist + dist)
                    push        = diff/dist
                    vector_push = vector*push
                    vector_push_len = get_vector_len_fastsingle(vector_push)
                    
                    points_arr[pi1]     += vector_push
                    pushes[pi1]         += vector_push_len # abs(push)
                    dists_traveled[pi1] += vector_push_len # abs(push)
            
            if self.plot_grid:
                POINT_STEPS.append(points_arr.copy())

            ### Checks if anything was pushed
            ### Values are floats so can't use '=='
            max_push = np.max(pushes)
            if max_push < push_tolerance and all_contained:
                break

        steps_toc = time.time()
        steps_time = steps_toc - steps_tic
        mean_steps_time = steps_time/end
        self.print_term("")
        self.print_term("Last step:         ", step, spaces=3)
        self.print_term("Optimization time: ", round(steps_time, 4), "[s]", spaces=3)
        self.print_term("Mean step time:    ", round(steps_time/end, 4), "[s]", spaces=3)
#         self.print_term("    ", "    ", "    ", "Max push:          ", max_push)
        
        return points_arr, POINT_STEPS, step, steps_time, mean_steps_time, max_push, bin_size

    ######################
    ### LIPID INSERTER ###
    ######################
    def lipid_inserter(self):
        if len(self.MEMBRANES) != 0:
            lipid_inserter_tic = time.time()
            self.print_term("------------------------------ CREATING LIPIDS", spaces=0)
            for memb_key, memb_dict in self.MEMBRANES.items():
                self.print_term("Starting membrane nr", memb_key, spaces=0)
                for leaflet_i, (leaflet_key, leaflet) in enumerate(memb_dict["leaflets"].items()):
                    self.print_term("Starting", leaflet["leaflet_type"], "leaflet", spaces=1)
                        
                    if leaflet["HG_direction"] == "up":
                        sign = +1
                    if leaflet["HG_direction"] == "down":
                        sign = -1
                    leaflet["grid_lipids"] = []
                    
                    for (slxi, slyi), subleaflet in leaflet["subleaflets"].items():
                        lipids                    = subleaflet["lipids"]
                        grid_points               = subleaflet["grid_points"]
                        for (l_name, l_radius), (grid_point_x, grid_point_y, grid_point_z) in zip(lipids, grid_points):
                            new_x, new_y, new_z = [], [], []
                            leaflet["grid_lipids"].append({
                                "grid_point_x": grid_point_x,
                                "grid_point_y": grid_point_y,
                                "grid_point_z": grid_point_z,
                            })
                            random_rotaion_angle = random.uniform(0, 360)
                            for x, y, z in leaflet["lipids"][l_name].get_beads("xyz"):
                                nx, ny, nz = self.rotate_point(x, y, z, 0, 0, random_rotaion_angle)
                                new_x.append(nx)
                                new_y.append(ny)
                                new_z.append(nz)
                            
                            lipid_dict = {
                                "name": l_name,
                                "beads": [bead.bead for bead in leaflet["lipids"][l_name].get_res_beads_info()],
                                "x": [x      + grid_point_x + random.uniform(-leaflet["kickx"], leaflet["kickx"]) for x in new_x],
                                "y": [y      + grid_point_y + random.uniform(-leaflet["kicky"], leaflet["kicky"]) for y in new_y],
                                "z": [z*sign + grid_point_z + random.uniform(-leaflet["kickz"], leaflet["kickz"]) for z in new_z],
                            }
                            leaflet["grid_lipids"][-1].update({"lipid": lipid_dict})
                            self.system_charge += leaflet["lipids"][l_name].charge
            lipid_inserter_toc = time.time()
            lipid_inserter_time = round(lipid_inserter_toc - lipid_inserter_tic, 4)
            self.print_term("------------------------------ LIPID CREATION COMPLETE", "(time spent: "+str(lipid_inserter_time)+" [s])", "\n", spaces=0)

    ################
    ### SOLVATER ###
    ################
    def solvater(self):
        if len(self.SOLVATIONS_cmds) != 0:
            solvation_tic = time.time()
            solv_beads_for_cell_checker = []
            self.print_term("------------------------------ SOLVATING SYSTEM", spaces=0)
            for solvation_i, (solvation_nr, solvation) in enumerate(self.SOLVATIONS.items()):
                if solvation_i != 0:
                    self.print_term("")
                self.print_term("Starting solvation nr", solvation_nr, spaces=0)

                #########################################
                ### CALCULATION FREE VOLUME OF SYSTEM ###
                #########################################
                self.print_term("Calculating box volume: (all values in [nm^3])", spaces=1)
                self.print_term("Bead radius used for volume calculations 'bead_radius':", solvation["bead_radius"], "[Å]", spaces=2)

                bead_radius = solvation["bead_radius"]
                ### Leaflet volume based on beads exclusively
                non_solv_beads = []
                prot_beads_for_cell_checker = []
                lipid_beads_for_cell_checker = []
                leafs_volume = 0
                if len(self.MEMBRANES) > 0:
                    '''
                    Estimates a volume for all lipids, and finds the lipid bead positions
                    '''
                    for memb_key, memb_dict in self.MEMBRANES.items():
                        for leaflet_key, leaflet in memb_dict["leaflets"].items():
                            for grid_point in leaflet["grid_lipids"]:
                                xs, ys, zs = grid_point["lipid"]["x"], grid_point["lipid"]["y"], grid_point["lipid"]["z"]
                                lipid_beads = list(zip(xs, ys, zs))
                                leafs_volume += (4/3 * math.pi * (bead_radius ** 3)) * len(lipid_beads) * 10**-27
                                non_solv_beads.extend(lipid_beads)
                                gridz = abs(max(zs)) + abs(min(zs))
                                lipid_beads_for_cell_checker.extend(lipid_beads)

                ### Protein volume based on beads exclusively
                prots_volume = 0
                if len(self.PROTEINS) > 0:
                    '''
                    Estimates a volume for all proteins, and finds the protein bead positions
                    '''
                    for protein_i, (protein_nr, protein) in enumerate(self.PROTEINS.items()):
                        prot_points = protein["protein"].get_beads("xyz")
                        prots_volume += (4/3 * math.pi * (bead_radius ** 3)) * len(prot_points) * 10**-27
                        prot_beads_for_cell_checker.extend(prot_points)

                ### Solvent volume from previous solvation commands
                solvs_volume = (4/3 * math.pi * (bead_radius ** 3)) * len(solv_beads_for_cell_checker) * 10**-27

                N_A = 6.02214076 * 10**23
                box_volume = (self.pbc_box[0] * self.pbc_box[1] * self.pbc_box[2]) * 10**-27
                non_free_volume = leafs_volume + prots_volume + prots_volume
                box_free_volume = box_volume - non_free_volume
                self.print_term("Box volume:            ", round(box_volume * 10**24, 3), spaces=2)
                self.print_term("Lipid volume:          ", round(leafs_volume * 10**24, 3), spaces=2)
                self.print_term("Protein volume:        ", round(prots_volume * 10**24, 3), spaces=2)
                self.print_term("(Prior) Solvent volume:", round(solvs_volume * 10**24, 3), spaces=2)
                self.print_term("Excluded volume:       ", round(non_free_volume * 10**24, 3), spaces=2)
                self.print_term("Free volume:           ", round(box_free_volume * 10**24, 3), spaces=2)

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
                            count = solv_vals.molarity
                        else: ### Treats molarity as molarity
                            ### Using the free volume for concentration calculations (default)
                            if solvation["solvfreevol"] == True:
                                count = int(N_A * box_free_volume * solv_vals.molarity / solv_vals.solvcount)
                            ### Using box volume for concentration calculations
                            elif solvation["solvfreevol"] == False:
                                count = int(N_A * box_volume * solv_vals.molarity / solv_vals.solvcount)
                        
                        ### Adds charge to solvent charge if it is specified in defines (0 by default)
                        sol_charges += count * solv_vals.charge

                        ### Checks if molar mass and density is specified for the solvent in defines
                        ### Only prints warning if solvent volume is to be used for ions
                        if not solv_vals.molar_mass and solvation["ionsvol"] == "solv" and not solvation["count"]:
                            self.print_term("WARNING: Chosen solvent [" + solv_key + "] is missing a 'molar_mass' entry in the solvent defines. Please add it if you wish to use solvent volume for ion concentration.", warn = True)
                        if not solv_vals.density and solvation["ionsvol"] == "solv" and not solvation["count"]:
                            self.print_term("WARNING: Chosen solvent [" + solv_key + "] is missing a 'density' entry in the solvent defines. Please add it if you wish to use solvent volume for ion concentration.", warn = True)

                        solvation["solvent"][solv_key].count_set(count)
                        ### Adds volume for solvent type to solvent volume if molar mass and density are both specified
                        if solv_vals.molar_mass and solv_vals.density:
                            solv_volume += (count * solv_vals.solvcount * solv_vals.molar_mass) / (N_A * solv_vals.density) * 10**-3

                        solvent_molecules.append(list(solv_vals.get_beads("xyz")))

                self.print_term("Solvent volume:        ", round(solv_volume * 10**24, 3), spaces=2)
                self.print_term(spaces=2)

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
                            count = pos_vals.molarity
                        else: ### Treats molarity as molarity
                            ### Using the solvent volume for concentration calculations (default)
                            if solvation["ionsvol"] == "solv":
                                count = int(N_A * solv_volume * pos_vals.molarity / pos_vals.solvcount)
                            ### Using the free volume for concentration calculations
                            elif solvation["ionsvol"] == "free":
                                count = int(N_A * box_free_volume * pos_vals.molarity / pos_vals.solvcount)
                            ### Using the box volume for concentration calculations
                            elif solvation["ionsvol"] == "box":
                                count = int(N_A * box_volume * pos_vals.molarity / pos_vals.solvcount)

                        solvation["pos_ions"][pos_key].count_set(count)
                        pos_charges += count * pos_vals.charge
                        pos_ratios.append((pos_key, pos_vals.charge, count)) # Used for ion optimizer/neutralization

                        solvent_molecules.append(list(pos_vals.get_beads("xyz")))

                neg_charges = 0
                neg_ratios = []

                if solvation["neg_ions"] != []:
                    ### Negative ions
                    for neg_key, neg_vals in solvation["neg_ions"].items():
                        if solvation["count"]: ### Treats molarity as absolute number of molecules
                            count = neg_vals.molarity
                        else: ### Treats molarity as molarity
                            ### Using the solvent volume for concentration calculations (default)
                            if solvation["ionsvol"] == "solv":
                                count = int(N_A * solv_volume * neg_vals.molarity / neg_vals.solvcount)
                            ### Using the free volume for concentration calculations
                            elif solvation["ionsvol"] == "free":
                                count = int(N_A * box_free_volume * neg_vals.molarity / neg_vals.solvcount)
                            ### Using the box volume for concentration calculations
                            elif solvation["ionsvol"] == "box":
                                count = int(N_A * box_volume * neg_vals.molarity / neg_vals.solvcount)

                        solvation["neg_ions"][neg_key].count_set(count)
                        neg_charges += count * neg_vals.charge
                        neg_ratios.append((neg_key, neg_vals.charge, count)) # Used for ion optimizer/neutralization

                        solvent_molecules.append(list(neg_vals.get_beads("xyz")))

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
                                solvation["neg_ions"][name].count_set(vals.count + count)
                                neg_charges += vals.charge * vals.count

                        else:
                            self.print_term("WARNING: I cannot neutralize with NEGATIVE ions when none are given")

                    ### Adds extra positive ions
                    elif tot_charges < solvation["charge"]:
                        if solvation["pos_ions"] != []:
                            pos_charges = 0
                            charge_diff = abs(solvation["charge"] - tot_charges)
                            extra_ions = ions_optimizer(pos_ratios, charge_diff)

                            for (name, vals), count in zip(solvation["pos_ions"].items(), extra_ions):
                                solvation["pos_ions"][name].count_set(vals.count + count)
                                pos_charges += vals.charge * vals.count

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
                        scount = vals.count
                        maxx, maxy, maxz, minx, miny, minz = [
                            func([
                                val + (solvation["WR"] * sign)
                                for val in solvation[stype][key].get_coords(ax)[0]
                            ])
                            for func, sign in [(max, +1), (min, -1)] for ax in ["x", "y", "z"]
                        ]
                        ssize = max([maxi - mini for maxi, mini in [(maxx, minx), (maxy, miny), (maxz, minz)]])
                        solvent_sizes.append(ssize)
                        collected_solvent.extend([[sname, stype, ssize]] * scount)
                                
                ##############################################
                ### RUNNING SOLVENT OPTIMIZATION ALGORITHM ###
                ##############################################
                ### Finds the maximum size of molecules used as solvent/ions
                ### Also includes buffer/kick size to prevent edge overlap cases
                max_mol_size = max([max([math.dist(bead1, bead2) for bead1 in molecule for bead2 in molecule]) for molecule in solvent_molecules]) #+ solvation["kick"] + solvation["buffer"]
                    
                self.print_term("gridres first:", solvation["gridres"], debug=True)
                self.print_term("max_mol_size:", max_mol_size, debug=True)
                self.print_term("max_mol_size*1.2:", max_mol_size*1.2, debug=True)
                self.print_term("kick:", solvation["kick"], debug=True)
                self.print_term("kick*1.2:", solvation["kick"]*1.2, debug=True)
                
                ### if the maximum molecule size is bigger than the designated grid resolution then change the gridres
                if max_mol_size*1.2 >= solvation["gridres"]:
                    gridres = (max_mol_size + solvation["kick"]*2) * 1.2 # 20% larger than largest molecule
                    self.print_term("NOTE: Requested solvent is too large for the grid resolution. Adjusting grid resolution to prevent solvent molecule overlaps.", warn = True)
                    self.print_term("Original grid resolution was:", round(solvation["gridres"]/10, 4), "[nm]", warn = True)
                    self.print_term("New grid resolution is:      ", round(gridres/10, 4), "[nm]", "\n", warn = True)
                else:
                    gridres = solvation["gridres"]
                    
                self.print_term("gridres mid:", gridres, debug=True)
                
                if (max_mol_size+solvation["kick"])*1.2 >= gridres:
                    self.print_term("NOTE: Kick is too large for grid resolution. Adjusting grid resolution to prevent solvent molecule overlaps.", warn = True)
                    self.print_term("Original grid resolution was:", round(gridres/10, 4), "[nm]", warn = True)
                    self.print_term("Original kick was:           ", round(solvation["kick"]/10, 4), "[nm]", warn = True)
                    gridres = gridres + solvation["kick"]*2*1.2 # 20% extra
                    self.print_term("New grid resolution is:      ", round(gridres/10, 4), "[nm]", "\n", warn = True)
                
                self.print_term("gridres last:", gridres, debug=True)
                
                solvent_buffer = gridres + solvation["buffer"]

                #################################
                ### CHOOSES ALGORITHM VERSION ###
                #################################
                self.print_term("Calculating the number of available grid points", spaces=1)
                if solvation["algorithm"] == "v1":
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
                        if len(self.MEMBRANES) != 0:
                            '''
                            pos: the "positive" (coordinate) side of a hydrophobic box in a given dimension (x/y/z)
                            neg: the "negative" (coordinate) side of a hydrophobic box in a given dimension (x/y/z)
                            '''
                            for memb_key, memb_dict in self.MEMBRANES.items():
                                for leaflet_key, leaflet in memb_dict["leaflets"].items():
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
                        if len(self.MEMBRANES) != 0:
                            '''
                            Finds solvent cells that are only partially contained in a hydrophobic box.
                            Used to find edge cases where cells are partially contained in a hydrophobic box without overlapping with a particle.

                            pos: the "positive" (coordinate) side of a hydrophobic box in a given dimension (x/y/z)
                            neg: the "negative" (coordinate) side of a hydrophobic box in a given dimension (x/y/z)

                            top: the "positive" (coordinate) side of a solvent cell in a given dimension (x/y/z)
                            bot: the "negative" (coordinate) side of a solvent cell in a given dimension (x/y/z)
                            '''
                            for memb_key, memb_dict in self.MEMBRANES.items():
                                for leaflet_key, leaflet in memb_dict["leaflets"].items():
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

                elif solvation["algorithm"] == "v2":
                    ### Creates a 3D matrix indicating the points in space that are allowed to have solvent
                    grid_bool_matrix = np.ones((int(self.pbcx/gridres), int(self.pbcy/gridres), int(self.pbcz/gridres)))
                    xpoints, ypoints, zpoints = grid_bool_matrix.shape
                    ### Calculates actual coordinate ranges for each axis and calculates the "real" grid resolution
                    xcoords, xreal_gridres = np.linspace(-self.pbcx//2-gridres/2, self.pbcx//2+gridres/2, xpoints+2, retstep=True)#[1:-1]
                    ycoords, yreal_gridres = np.linspace(-self.pbcy//2-gridres/2, self.pbcy//2+gridres/2, ypoints+2, retstep=True)#[1:-1]
                    zcoords, zreal_gridres = np.linspace(-self.pbcz//2-gridres/2, self.pbcz//2+gridres/2, zpoints+2, retstep=True)#[1:-1]
                    ### Real grid resolution is different from given to to axis lengths not being perfectly divisible by given grid resolution
#                     xreal_gridres = abs(xcoords[0]) - abs(xcoords[1])
#                     yreal_gridres = abs(ycoords[0]) - abs(ycoords[1])
#                     zreal_gridres = abs(zcoords[0]) - abs(zcoords[1])
                    ### Removes the first and last points as they are the actual edges of the box
                    xcoords = xcoords[1:-1]
                    ycoords = ycoords[1:-1]
                    zcoords = zcoords[1:-1]

                    def coord_to_indices(pos, dim, real_gridres, min_buffer, max_buffer):
                        ### Buffer limit
                        bead_min = pos-min_buffer
                        bead_max = pos+max_buffer

                        ### Convert solvent buffer limit to decimal index
                        beadi_min_dec = (bead_min+dim/2)/real_gridres
                        beadi_max_dec = (bead_max+dim/2)/real_gridres

                        ### Round buffer decimal index to integer index
                        beadi_min_int = math.floor(beadi_min_dec+0.5)
                        beadi_max_int = math.ceil(beadi_max_dec-0.5)
                        
                        if beadi_min_int < 0:
                            beadi_min_int = 0

                        return beadi_min_int, beadi_max_int


                    ### ### Marks coordinates as "occupied" by setting boolean to False
                    ### Marks points located in hydrophobic volume
                    if len(self.MEMBRANES) != 0:
                        for memb_key, memb_dict in self.MEMBRANES.items():
                            for leaflet_key, leaflet in memb_dict["leaflets"].items():
                                lcx, lcy, lcz = leaflet["center"] # Center of leaflet on given axis
                                llx, lly = leaflet["x"], leaflet["y"] # Length of leaflet in given axis
                                lhydrophob = leaflet["lipid_dimensions"]["zUnderLipCen"] # Height of hydrophobic volume
                                if leaflet["HG_direction"] == "up":
                                    zminbuffer = solvent_buffer
                                    zmaxbuffer = lhydrophob + solvent_buffer
                                if leaflet["HG_direction"] == "down":
                                    zminbuffer = lhydrophob + solvent_buffer
                                    zmaxbuffer = solvent_buffer

                                xmin, xmax = coord_to_indices(lcx, self.pbcx, xreal_gridres, llx/2, llx/2)
                                ymin, ymax = coord_to_indices(lcy, self.pbcy, yreal_gridres, lly/2, lly/2)
                                zmin, zmax = coord_to_indices(lcz, self.pbcz, zreal_gridres, zminbuffer, zmaxbuffer)
                                grid_bool_matrix[xmin:xmax, ymin:ymax, zmin:zmax] = 0

                    ### Protein beads
                    protein_buffer = solvent_buffer+solvation["protein_extra_buffer"]
                    for xpos, ypos, zpos in prot_beads_for_cell_checker:
                        ### Checks if any non-solvent bead is within the solvent cell
                        xmin, xmax = coord_to_indices(xpos, self.pbcx, xreal_gridres, protein_buffer, protein_buffer)
                        ymin, ymax = coord_to_indices(ypos, self.pbcy, yreal_gridres, protein_buffer, protein_buffer)
                        zmin, zmax = coord_to_indices(zpos, self.pbcz, zreal_gridres, protein_buffer, protein_buffer)
                        grid_bool_matrix[xmin:xmax, ymin:ymax, zmin:zmax] = 0

                    ### Lipid beads
                    lipid_buffer = solvent_buffer+solvation["lipid_extra_buffer"]
                    for xpos, ypos, zpos in lipid_beads_for_cell_checker:
                        ### Checks if any non-solvent bead is within the solvent cell
                        xmin, xmax = coord_to_indices(xpos, self.pbcx, xreal_gridres, lipid_buffer, lipid_buffer)
                        ymin, ymax = coord_to_indices(ypos, self.pbcy, yreal_gridres, lipid_buffer, lipid_buffer)
                        zmin, zmax = coord_to_indices(zpos, self.pbcz, zreal_gridres, lipid_buffer, lipid_buffer)
                        grid_bool_matrix[xmin:xmax, ymin:ymax, zmin:zmax] = 0

                    ### Prior solvent beads
                    solute_buffer = solvent_buffer+solvation["solute_extra_buffer"]
                    for xpos, ypos, zpos in solv_beads_for_cell_checker:
                        ### Checks if any non-solvent bead is within the solvent cell
                        xmin, xmax = coord_to_indices(xpos, self.pbcx, xreal_gridres, solute_buffer, solute_buffer)
                        ymin, ymax = coord_to_indices(ypos, self.pbcy, yreal_gridres, solute_buffer, solute_buffer)
                        zmin, zmax = coord_to_indices(zpos, self.pbcz, zreal_gridres, solute_buffer, solute_buffer)
                        grid_bool_matrix[xmin:xmax, ymin:ymax, zmin:zmax] = 0

                    ### Gets all the indices that are not occupied
                    free_indices = grid_bool_matrix.nonzero()
                    free_xs, free_ys, free_zs = free_indices
                    free_xs_coords = xcoords[free_xs]
                    free_ys_coords = ycoords[free_ys]
                    free_zs_coords = zcoords[free_zs]
                    solv_grid_3D = list(zip(free_xs_coords, free_ys_coords, free_zs_coords))

                self.print_term("Final number of 3D grid points available for solvent placement:", len(solv_grid_3D), spaces=1)

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
                self.print_term("Inserting", len(collected_solvent), "solvent molecules into random grid points:", spaces=1)
                random_grid_points = random.sample(solv_grid_3D, k = len(collected_solvent))
                grid_solvated = []
#                 counter = 0
                generated_spots = []
                for counter, ((sname, stype, ssize), (gx, gy, gz)) in enumerate(zip(collected_solvent, random_grid_points), 1):
                    '''
                    Finds a 3D-grid point for the solvent and places it there.
                    '''
                    if counter == 1 or counter % 25000 == 0:
                        self.print_term("Currently at solvent number:", counter, spaces=2)
                    sdata = solvation[stype][sname]
                    sinfo = solvation[stype][sname].get_res_beads_info()
                    
                    sbeads     = [bead.bead for bead in sinfo]
                    sresnames  = [bead.resname for bead in sinfo]
                    scharge    = sdata.charge
                    sx, sy, sz = sdata.get_coords("xyz")
                    rx, ry, rz = random.uniform(0, 360), random.uniform(0, 360), random.uniform(0, 360)
                    for j, (x, y, z) in enumerate(zip(sx, sy, sz)):
                        sx[j], sy[j], sz[j] = self.rotate_point(x, y, z, rx, ry, rz)
                        kx = random.uniform(-solvation["kick"], solvation["kick"])
                        ky = random.uniform(-solvation["kick"], solvation["kick"])
                        kz = random.uniform(-solvation["kick"], solvation["kick"])
                        sx[j], sy[j], sz[j] = sx[j] + kx + gx, sy[j] + ky + gy, sz[j] + kz + gz
                    
                    ### Shouldn't be possible for beads to be placed outside box using the current algorithms
#                     sx = position_fixer(sx, 0)
#                     sy = position_fixer(sy, 1)
#                     sz = position_fixer(sz, 2)

                    ### Order that beads should be written in for structure file and topology file
                    if stype == "solvent":
                        order = 1
                    if stype == "pos_ions":
                        order = 2
                    if stype == "neg_ions":
                        order = 3

                    grid_solvated.append({
                        "order":    order,
                        "name":     sname,
                        "resnames": sresnames,
                        "type":     stype,
                        "beads":    sbeads,
                        "charge":   scharge,
                        "coords":   list(zip(sx, sy, sz)),
                    })
                self.print_term("Currently at solvent number:", counter, spaces=2)
                
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
                SYS_solv_molarity = [round(scount / (N_A * solv_volume), 3) for scount in SYS_counts]
                SYS_printer = [tuple(SYS_headers)] + list(zip(SYS_names, SYS_counts, SYS_percentages, SYS_box_molarity, SYS_free_molarity, SYS_solv_molarity, SYS_charges))
                SYS_max_lengths = [len(str(max(data, key = lambda d: len(str(d))))) for data in zip(*SYS_printer)]
                self.print_term("Solvent data for whole system")
                for L0, L1, L2, L3, L4, L5, L6 in SYS_printer:
                    self.print_term(
                        '{0: <{L}}'.format(L0, L = SYS_max_lengths[0]), ":",
                        '{0: <{L}}'.format(L1, L = SYS_max_lengths[1]), ":",
                        '{0: <{L}}'.format(L2, L = SYS_max_lengths[2]), ":",
                        '{0: <{L}}'.format(L3, L = SYS_max_lengths[3]), ":",
                        '{0: <{L}}'.format(L4, L = SYS_max_lengths[4]), ":",
                        '{0: <{L}}'.format(L5, L = SYS_max_lengths[5]), ":",
                        '{0: <{L}}'.format(L6, L = SYS_max_lengths[6]),
                        spaces=2,
                    )
            
            solvation_toc = time.time()
            solvation_time = round(solvation_toc - solvation_tic, 4)
            self.print_term("------------------------------ SOLVATION COMPLETE", "(time spent: "+str(solvation_time)+" [s])", "\n", spaces=0)
    
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
#                             x, y = x + (x[0],), y + (y[0],)
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
            
            ### ### To unpickle use the following and change "pickled_file_path" to your pickled file name
            ### with open(pickled_file_path, 'rb') as pickled_file:
            ###     unpickled_class = pickle.load(pickled_file)
            
            pickled_data = self
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
parser.add_argument("--membrane", "-memb", "-membrane", dest = "membrane_cmds", action="append", type=str, default = [], nargs="+")

### Protein commands
parser.add_argument("--protein", "-prot", "-protein", dest = "protein_cmds", action="append", type=str, default = [], nargs="+")

### Solvent commands
parser.add_argument("--solvation", "-solv", "-solvation", dest = "solvation_cmds", action="append", type=str, default = [], nargs="+")

### Solvent commands
parser.add_argument("--flooding", "-flood", "-flooding", dest = "flooding_cmds", action="append", type=str, default = [], nargs="+")

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

### System parameters
parser.add_argument("--sys_params", "-params", dest = "sys_params", default = "default")
parser.add_argument("-prot_params", dest = "prot_params", default = False)
parser.add_argument("-lipid_params", dest = "lipid_params", default = False)
parser.add_argument("-solv_params", dest = "solv_params", default = False)

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
parser.add_argument("--output_struct"    , "-out"    , "-out_sys"    , "-o"    , dest = "out_system_file_name"    , default = False)
### Output pdb file
parser.add_argument("--output_struct_pdb", "-out_pdb", "-out_sys_pdb", "-o_pdb", dest = "out_system_pdb_file_name", default = False)
### Output gro file
parser.add_argument("--output_struct_gro", "-out_gro", "-out_sys_gro", "-o_gro", dest = "out_system_gro_file_name", default = False)

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
parser.add_argument("--print_quiet",    "-quiet", dest = "quiet",    default = False)
parser.add_argument("--print_debug",    "-debug", dest = "debug",    default = False)
parser.add_argument("--print_extra",    "-extra", dest = "extra",    default = True)
parser.add_argument("--print_warnings", "-warn",  dest = "warnings", default = True)

### ### Parser for handling '-f' when importing module to Jupyter
parser.add_argument("-f", dest = "debug_flag_for_jupyter")

args = parser.parse_args()

parse_membrane_cmds      = [" ".join(i) for i in args.membrane_cmds]
parse_protein_cmds       = [" ".join(i) for i in args.protein_cmds]
parse_solvation_cmds     = [" ".join(i) for i in args.solvation_cmds]
parse_flooding_cmds      = [" ".join(i) for i in args.flooding_cmds]
parse_itp_input_cmds     = [" ".join(i) for i in args.itp_input_cmds]
parse_solute_input_cmds  = [" ".join(i) for i in args.solute_input_cmds]

parse_plot_cmd   = args.plot_cmd
parse_pickle_cmd = args.pickle_cmd
parse_backup     = bool(ast.literal_eval(str(args.backup)))
parse_randseed   = args.randseed

parse_sys_params   = args.sys_params
parse_prot_params  = args.prot_params
parse_lipid_params = args.lipid_params
parse_solv_params  = args.solv_params

parse_pbc_box = args.pbc_box
if parse_pbc_box:
    parse_pbc_box = [ast.literal_eval(str(i)) for i in parse_pbc_box]

parse_pbcx = ast.literal_eval(str(args.pbcx))
parse_pbcy = ast.literal_eval(str(args.pbcy))
parse_pbcz = ast.literal_eval(str(args.pbcz))
if parse_pbcx and parse_pbcy and parse_pbcz:
    parse_pbc_box = [parse_pbcx, parse_pbcy, parse_pbcz]

parse_out_system_file_name     = args.out_system_file_name
parse_out_system_pdb_file_name = args.out_system_pdb_file_name
parse_out_system_gro_file_name = args.out_system_gro_file_name
parse_out_topol_file_name      = args.out_topol_file_name
# parse_output_imported  = args.output_imported
parse_system_name      = args.system_name

parse_out_log_file_name = args.out_log_file_name

parse_quiet        = bool(ast.literal_eval(str(args.quiet)))
parse_debug_prints = bool(ast.literal_eval(str(args.debug)))
parse_extra_info   = bool(ast.literal_eval(str(args.extra)))
parse_warnings     = bool(ast.literal_eval(str(args.warnings)))

if any([i != [] for i in [parse_membrane_cmds, parse_protein_cmds, parse_solvation_cmds]]):
    CGSB(
        box = parse_pbc_box,
        
        membrane  = parse_membrane_cmds,
        protein   = parse_protein_cmds,
        solvation = parse_solvation_cmds,
        flooding  = parse_flooding_cmds,
        
        rand = parse_randseed,

        itp_input    = parse_itp_input_cmds,
        solute_input = parse_solute_input_cmds,

        plot         = parse_plot_cmd,
        pickle       = parse_pickle_cmd,
        backup       = parse_backup,
        sys_params   = parse_sys_params,
        prot_params  = parse_prot_params,
        lipid_params = parse_lipid_params,
        solv_params  = parse_solv_params,

        out_sys     = parse_out_system_file_name,
        out_sys_pdb = parse_out_system_pdb_file_name,
        out_sys_gro = parse_out_system_gro_file_name,
        
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

