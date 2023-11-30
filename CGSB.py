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
import shapely.affinity
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
    "W" : {"beads": "W",  "x": (0,), "y": (0,), "z": (0,), "mapping_ratio": 4, "density": 0.99669, "molar_mass": 18.01528},
    "SW": {"beads": "SW", "x": (0,), "y": (0,), "z": (0,), "mapping_ratio": 3, "density": 0.99669, "molar_mass": 18.01528},
    "TW": {"beads": "TW", "x": (0,), "y": (0,), "z": (0,), "mapping_ratio": 2, "density": 0.99669, "molar_mass": 18.01528},
}
### Amino acids
solvent_defs[params].update({
    ### ### 1 atom residues
    "GLY":  {"beads": ("BB",), "mapping_ratio": 1, "x": (0,), "y": (0,), "z": (0,)},
    
    ### ### 2 atom residues
    ### Residues without variants
    "ALA":  {"beads": ("BB", "SC1"), "mapping_ratio": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
    "CYS":  {"beads": ("BB", "SC1"), "mapping_ratio": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
    "VAL":  {"beads": ("BB", "SC1"), "mapping_ratio": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
    "LEU":  {"beads": ("BB", "SC1"), "mapping_ratio": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
    "ILE":  {"beads": ("BB", "SC1"), "mapping_ratio": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
    "MET":  {"beads": ("BB", "SC1"), "mapping_ratio": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
    "PRO":  {"beads": ("BB", "SC1"), "mapping_ratio": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
    "HYP":  {"beads": ("BB", "SC1"), "mapping_ratio": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
    "ASN":  {"beads": ("BB", "SC1"), "mapping_ratio": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
    "GLN":  {"beads": ("BB", "SC1"), "mapping_ratio": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
    "THR":  {"beads": ("BB", "SC1"), "mapping_ratio": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
    "SER":  {"beads": ("BB", "SC1"), "mapping_ratio": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)},
    ### ASP variants
    "ASP":  {"beads": ("BB", "SC1"), "mapping_ratio": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0), "charges": ("SC1", -1)},
    "ASPP": {"beads": ("BB", "SC1"), "mapping_ratio": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)}, # Neutral ASP
    "ASH":  {"beads": ("BB", "SC1"), "mapping_ratio": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)}, # Neutral ASP
    ### Glutamate variants
    "GLU":  {"beads": ("BB", "SC1"), "mapping_ratio": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0), "charges": ("SC1", -1)},
    "GLUP": {"beads": ("BB", "SC1"), "mapping_ratio": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)}, # Neutral GLU
    "GLH":  {"beads": ("BB", "SC1"), "mapping_ratio": 1, "x": (0.25, -0.25), "y": (0, 0), "z": (0, 0)}, # Neutral GLU
    
    ### ### 3 atom residues
    "ARG":  {"beads": ("BB", "SC1", "SC2"), "mapping_ratio": 1, "x": (0.25, 0, -0.25), "y": (0, 0, 0.125), "z": (0, 0, 0), "charges": ("SC2", 1)},
    ### Lysine variants
    "LYS":  {"beads": ("BB", "SC1", "SC2"), "mapping_ratio": 1, "x": (0.25, 0, -0.25), "y": (0, 0, 0.125), "z": (0, 0, 0), "charges": ("SC2", 1)},
    "LSN":  {"beads": ("BB", "SC1", "SC2"), "mapping_ratio": 1, "x": (0.25, 0, -0.25), "y": (0, 0, 0.125), "z": (0, 0, 0)}, # Neutral LYS
    "LYN":  {"beads": ("BB", "SC1", "SC2"), "mapping_ratio": 1, "x": (0.25, 0, -0.25), "y": (0, 0, 0.125), "z": (0, 0, 0)}, # Neutral LYS
    
    ### ### 4 atom residues
    "PHE":  {"beads": ("BB", "SC1", "SC2", "SC3"), "mapping_ratio": 1, "x": (0.25, 0, -0.25, -0.25), "y": (0, 0, 0.125, -0.125), "z": (0, 0, 0, 0)},
    ### Histidines
    "HIS":  {"beads": ("BB", "SC1", "SC2", "SC3"), "mapping_ratio": 1, "x": (0.25, 0, -0.25, -0.25), "y": (0, 0, 0.125, -0.125), "z": (0, 0, 0, 0)},
    "HIE":  {"beads": ("BB", "SC1", "SC2", "SC3"), "mapping_ratio": 1, "x": (0.25, 0, -0.25, -0.25), "y": (0, 0, 0.125, -0.125), "z": (0, 0, 0, 0)},
    "HSE":  {"beads": ("BB", "SC1", "SC2", "SC3"), "mapping_ratio": 1, "x": (0.25, 0, -0.25, -0.25), "y": (0, 0, 0.125, -0.125), "z": (0, 0, 0, 0)},
    "HSD":  {"beads": ("BB", "SC1", "SC2", "SC3"), "mapping_ratio": 1, "x": (0.25, 0, -0.25, -0.25), "y": (0, 0, 0.125, -0.125), "z": (0, 0, 0, 0)},
    "HID":  {"beads": ("BB", "SC1", "SC2", "SC3"), "mapping_ratio": 1, "x": (0.25, 0, -0.25, -0.25), "y": (0, 0, 0.125, -0.125), "z": (0, 0, 0, 0)},
    "HSP":  {"beads": ("BB", "SC1", "SC2", "SC3"), "mapping_ratio": 1, "x": (0.25, 0, -0.25, -0.25), "y": (0, 0, 0.125, -0.125), "z": (0, 0, 0, 0), "charges": (("SC2", 0.5), ("SC3", 0.5))},
    "HIP":  {"beads": ("BB", "SC1", "SC2", "SC3"), "mapping_ratio": 1, "x": (0.25, 0, -0.25, -0.25), "y": (0, 0, 0.125, -0.125), "z": (0, 0, 0, 0), "charges": (("SC2", 0.5), ("SC3", 0.5))},
    
    ### ### 5 atom residues
    "TYR":  {"beads": ("BB", "SC1", "SC2", "SC3", "SC4"), "mapping_ratio": 1, "x": (0.25, 0.25, 0, 0, -0.25), "y": (0.125, 0, -0.125, 0.125, 0), "z": (0, 0, 0, 0, 0)},
    
    ### ### 6 atom residues
    "TRP":  {"beads": ("BB", "SC1", "SC2", "SC3", "SC4", "SC5"), "mapping_ratio": 1, "x": (0.25, 0.25, 0, 0, -0.25, -0.25), "y": (0.125, 0, -0.125, 0.125, 0, 0.125), "z": (0, 0, 0, 0, 0, 0)},
})

### Example of a multi-residue solvent molecule
# solvent_defs[params].update({
#     "LIG1": {
#         "residues": [
#             {
#                 "resname": "RES1",
#                 "beads":   ("R1A1", "R1A2", "R1A3", "R1A4", "R1A5"),
#                 "x":       (     1,      2,      3,      4,      5),
#                 "y":       (     0,      0,    0.5,      0,   -0.5),
#                 "z":       (     0,      1,      0,    0.5,      0),
#                 "charges": (R1A2, 1),
#             },
#             {
#                 "resname": "RES2",
#                 "beads":   ("R2A1", "R2A2", "R2A3", "R2A4", "R2A5"),
#                 "x":       (     1,      2,      3,      4,      5),
#                 "y":       (    -1,      0,      0,   -0.5,      0),
#                 "z":       (     1,      2,      1,    1.5,      1),
#                 "charges": (R2A5, -1),
#             },
#         ],
#         "mapping_ratio": 1,
#         "charges": ((RES1, R1A2, 1), (RES2, R2A5, -1)),
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
    def __init__(self, bead, beadnr, x, y, z, resname, resnumber, charge=0):
        self.bead   = bead
        self.beadnr = beadnr
        self.x = x
        self.y = y
        self.z = z
        self.resname = resname
        self.resnr   = resnumber
        self.charge  = charge
        
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
    
    def set_charge(self, charge):
        self.charge = charge
    
    def get_tuple(self):
        return (self.bead, self.beadnr, self.x, self.y, self.z, self.resname, self.resnr, self.charge)
    
class RESIDUE:
    def __init__(self, resname, resnumber):
        self.resname = resname
        self.resnr = resnumber
        self.beads = []
        
    def add_bead_to_res(self, bead, beadnr, x, y, z, charge=0):
        self.beads.append(ATOM(bead, beadnr, x, y, z, self.resname, self.resnr, charge))
            
    def add_beads_to_res(self, beads = False, beadnrs = False, xs = False, ys = False, zs = False, charges = False):
        assert beads and beadnrs and xs and ys and zs, "Lacking data for either 'beads', 'beadnr', 'xs', 'ys' or 'zs'"
        if not charges:
            charges = [0 for _ in range(len(beads))]
        for bead, beadnr, x, y, z, charge in zip(beads, beadnrs, xs, ys, zs, charges):
            self.beads.append(ATOM(bead, beadnr, x, y, z, self.resname, self.resnr, charge))
    
    def add_bead_data_to_res(self, bead_data):
#         assert len(bead_data) == 5, "Length of list is " + str(len(bead_data)) + ". Lacking data for either 'beads', 'xs', 'ys' or 'zs'"
        if len(bead_data) == 5:
            ### Adding charge if not given
            bead_data.append([0 for _ in range(len(bead_data[0]))])
        for bead, beadnr, x, y, z, charge in zip(bead_data):
            self.beads.append(ATOM(bead, beadnr, x, y, z, self.resname, self.resnr, charge))
    
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
    def __init__(self, molname = False, moleculetype = False):
        self.residues = []
        self.center = False
        self.n_residues = 0
        self.resnames = []
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
    
    def add_bead(self, bead, beadnr, x, y, z, resnumber = False, charge=0):
        if not resnumber:
            resnumber = self.last_res_n
        self.residues[resnumber].add_bead_to_res(bead, beadnr, x, y, z, charge)
    
    def add_res_and_beads(self, resname, beads = False, beadnrs = False, xs = False, ys = False, zs = False, bead_data = False, resnumber = False, charges = False):
        if not resnumber:
            resnumber = self.n_residues
        self.residues.append(RESIDUE(resname, resnumber))
        self.resnames.append(resname)
        self.last_res_n = self.n_residues
        self.n_residues += 1
        if bead_data:
            self.residues[self.last_res_n].add_bead_data_to_res(bead_data)
        else:
            self.residues[self.last_res_n].add_beads_to_res(beads, beadnrs, xs, ys, zs, charges)
    
    def set_bead_charges(self, charges):
        i = 0
        for ri, res in enumerate(self.residues):
            for bi, bead in enumerate(res.beads):
                self.residues[ri].beads[bi].set_charge(charges[i])
                i += 1
    
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
#         print(AXsList)
#         print(list(zip(*AXsList)))
        ### Not sure why following assert was made
#         assert len(AXsList) > 1, "Length of AXsList must be greater than 1 to create beads"
        return list(zip(*AXsList))

    def get_mol_charge(self):
        charge = sum(self.get_bead_charges())
        return charge

    def get_bead_charges(self):
        charges = [bead.charge for res in self.residues for bead in res.beads]
        return charges

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
    
    def get_radius(self, AXs="xyz"):
        beads  = self.get_beads(AXs)
        center = self.get_center_point("ax", AXs)
        radius = max([math.dist(bead, center) for bead in beads])
        return radius
    
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
                    [1, 0,                0              ],
                    [0, math.cos(x_rad), -math.sin(x_rad)],
                    [0, math.sin(x_rad),  math.cos(x_rad)],
                ]
                rm_y = [
                    [math.cos(y_rad),  0, math.sin(y_rad)],
                    [0,                1, 0              ],
                    [-math.sin(y_rad), 0, math.cos(y_rad)],
                ]
                rm_z = [
                    [math.cos(z_rad), -math.sin(z_rad), 0],
                    [math.sin(z_rad),  math.cos(z_rad), 0],
                    [0,                0,               1],
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
    def __init__(self, hydr_z = 0, molname = False, moleculetype = False):
        super().__init__(molname = molname, moleculetype = False)
        self.hydr_z = hydr_z ### Hydrophobic z-height delimiter
        self.ratio = 0
    
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
    def __init__(self, mapping_ratio = 1, molname = False, moleculetype = False):
        super().__init__(molname = molname, moleculetype = False)
        self.mapping_ratio = mapping_ratio ### AA-to-CG convertion
        self.molarity = False
        self.density = False
        self.molar_mass = False
        
    def mapping_ratio_set(self, mapping_ratio):
        self.mapping_ratio = mapping_ratio

    def molarity_set(self, molarity):
        self.molarity = molarity
    
    def density_set(self, density):
        self.density = density
    
    def ratio_set(self, ratio):
        self.ratio = ratio
    
    def molar_mass_set(self, molar_mass):
        self.molar_mass = molar_mass
    
    def count_set(self, count):
        self.count = count
    
class PROTEIN(MOLECULE):
    def __init__(self, molname = False, moleculetype = False):
        super().__init__(molname = molname, moleculetype = False)

class CGSB:
    
    def __init__(self, run = True, terminal_run_kwargs = False, **kwargs):
        self.RUN = run
        
        if terminal_run_kwargs:
            kwargs.update(terminal_run_kwargs)
        
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
        
        ### Stakced membranes are special combinations of membranes and solvations
        self.STACKED_MEMBRANES_cmds = []
        
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
        
        self.pbcx = 0
        self.pbcy = 0
        self.pbcz = 0
        
        self.verbose = 6
        
        if self.RUN:
            self.run(kwargs)
    
    ##############################
    ### GIVE COMMANDS TO CLASS ###
    ##############################
    def flooding_to_solvation_converter(self, subcmd):
        if "flooding:" in subcmd:
            subcmd_split = subcmd.split()
            subcmd = " ".join([string for string in subcmd_split if not string.startswith("flooding:")])
        subcmd = " ".join(["flooding:True", subcmd, "count:True", "solv_molarity:1", "salt_molarity:1"])
        return subcmd
    
    def commands_handler(self, kwargs):
        momentary_pbc = []
        momentary_x = 0
        momentary_y = 0
        momentary_z = 0
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
#                     if "flooding:" in subcmd:
#                         subcmd_split = subcmd.split()
#                         subcmd = " ".join([string for string in subcmd_split if not string.startswith("flooding:")])
#                     subcmd = " ".join(["flooding:True", subcmd, "count:True", "solv_molarity:1", "salt_molarity:1"])
                    subcmd = self.flooding_to_solvation_converter(subcmd)
                    self.FLOODINGS_cmds.extend([subcmd])
            
            if any(key.startswith(i) for i in ["stacked_membranes", "stack_memb"]):
                if type(cmd) not in [list, tuple]:
                    cmd = [cmd]
                for subcmd in cmd:
                    self.STACKED_MEMBRANES_cmds.extend([subcmd])
            
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
                for i, val in enumerate(cmd):
                    if type(val) == str:
                        isnumber, isint = self.is_number(val)
                        if isint:
                            val = int(val)
                        else:
                            val = float(val)
                    momentary_pbc.append(val)
                    
            if key == "x":
                val = cmd
                if type(val) == str:
                    isnumber, isint = self.is_number(val)
                    if isint:
                        val = int(val)
                    else:
                        val = float(val)
                momentary_x = val
                
            if key == "y":
                val = cmd
                if type(val) == str:
                    isnumber, isint = self.is_number(val)
                    if isint:
                        val = int(val)
                    else:
                        val = float(val)
                momentary_y = val
                
            if key == "z":
                val = cmd
                if type(val) == str:
                    isnumber, isint = self.is_number(val)
                    if isint:
                        val = int(val)
                    else:
                        val = float(val)
                momentary_z = val
            
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
            if key in ["out_all", "o_all"]:
                ### Cuts the extension so that all files can be generated with proper extensions
                if any(cmd.endswith(string) for string in [".pdb", ".gro", ".top", ".log"]):
                    cmd = cmd[:-4]
                self.output_system_pdb_file_name = cmd + ".pdb"
                self.output_system_gro_file_name = cmd + ".gro"
                self.output_topol_file_name      = cmd + ".top"
                self.output_log_file_name        = cmd + ".log"
            
            if key in ["out_sys", "o_sys"]:
                if not any([cmd.endswith(i) for i in [".pdb", ".gro"]]):
                    self.output_system_pdb_file_name = cmd + ".pdb"
                    self.output_system_gro_file_name = cmd + ".gro"
                elif cmd.endswith(".pdb"):
                    self.output_system_pdb_file_name = cmd
                elif cmd.endswith(".gro"):
                    self.output_system_gro_file_name = cmd
                else:
                    assert False, "Unknown file extension used for 'output_system': " + cmd
                    
            if key in ["out_pdb", "o_pdb"]:
                if not cmd.endswith("gro"):
                    cmd = cmd + ".pdb"
                self.output_system_pdb_file_name = cmd + ".pdb"
                    
            if key in ["out_gro", "o_gro"]:
                if not cmd.endswith("gro"):
                    cmd = cmd + ".gro"
                self.output_system_gro_file_name = cmd + ".gro"
            
            if key in ["out_top", "o_top"]:
                if not cmd.endswith(".top"):
                    cmd = cmd + ".top"
                self.output_topol_file_name = cmd
                
            if key in ["out_log", "o_log"]:
                if not cmd.endswith(".log"):
                    cmd = cmd + ".log"
                self.output_log_file_name = cmd
            
            if key in ["plot_grid"]:
                if type(cmd) == str:
                    cmd = ast.literal_eval(cmd)
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
                if type(cmd) == str:
                    cmd = ast.literal_eval(cmd)
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
                if type(cmd) == str:
                    cmd = ast.literal_eval(cmd)
                self.quiet = cmd
                
            if key in ["debug"]:
                if type(cmd) == str:
                    cmd = ast.literal_eval(cmd)
                self.debug_prints = cmd
                
            if key in ["extra"]:
                if type(cmd) == str:
                    cmd = ast.literal_eval(cmd)
                self.extra_info = cmd
            
            if key in ["warn"]:
                if type(cmd) == str:
                    cmd = ast.literal_eval(cmd)
                self.warnings = cmd
                
            if key == "verbose":
                number = self.get_number_from_string(cmd)
                if number is False:
                    self.verbose = len(cmd)
                else:
                    self.verbose = number
                
            ### Run the program
            if key in ["run"]:
                if type(cmd) == str:
                    cmd = ast.literal_eval(cmd)
                self.RUN = cmd
        
        ### Setting randseed
        self.print_term("\n" + "Setting random seed to:", self.randseed, verbose=1)
        random.seed(self.randseed)
        np.random.seed(self.randseed)
        
        self.print_term("------------------------------ PREPROCESSING DEFINITIONS", spaces=0, verbose=1)
        ### Definition preprocessing
        preprocessing_tic = time.time()
        self.itp_read_initiater()
        self.import_structures_handler()
        self.lipid_defs_preprocessor()
        self.solvent_defs_preprocessor()
        self.ion_defs_preprocessor()
        preprocessing_toc = time.time()
        preprocessing_time = round(preprocessing_toc - preprocessing_tic, 4)
        self.print_term("------------------------------ DEFINITIONS PREPROCESSING COMPLETE", "(time spent: "+str(preprocessing_time)+" [s])", "\n", spaces=0, verbose=1)
        
        if len(self.STACKED_MEMBRANES_cmds) > 0:
            self.print_term("------------------------------ PROCESSING SPECIAL COMMANDS", spaces=0, verbose=1)
            momentary_z = self.stacked_membranes_preprocessor()
            self.print_term("------------------------------ SPECIAL COMMAND PROCESSING COMPLETE", "(time spent: "+str(preprocessing_time)+" [s])", "\n", spaces=0, verbose=1)
        
        if len(self.FLOODINGS_cmds) > 0:
            self.SOLVATIONS_cmds = self.FLOODINGS_cmds + self.SOLVATIONS_cmds
        
        ### Setting box size values to be used in PBC type settings
        if momentary_pbc:
            cmd = momentary_pbc
        else:
            cmd = [i for i in [momentary_x, momentary_y, momentary_z] if i]

        ### PBC type settings:
        if self.pbc_type == "rectangular":
            if len(cmd) == 1:
                self.pbcx = cmd[0]*10
                self.pbcy = cmd[0]*10
                self.pbcz = cmd[0]*10
            elif len(cmd) == 2:
                self.pbcx = cmd[0]*10
                self.pbcy = cmd[0]*10
                self.pbcz = cmd[1]*10
                self.pbc_box = [self.pbcx, self.pbcy, self.pbcz]
            elif len(cmd) == 3:
                self.pbcx, self.pbcy, self.pbcz = [i*10 for i in cmd]
            else:
                assert False, "Incorrect pbc box dimensions for pbc type '"+self.pbc_type+"': " + str(cmd)
            self.gro_box_vectors = [
                float(self.pbcx/10), # vector: vax or v1x
                float(self.pbcy/10), # vector: vby or v2y
                float(self.pbcz/10), # vector: vcz or v3z
                float(0),            # vector: vay or v1y
                float(0),            # vector: vaz or v1z
                float(0),            # vector: vbx or v2x
                float(0),            # vector: vbz or v2z
                float(0),            # vector: vcx or v3x
                float(0),            # vector: vcy or v3y
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
            ### Uses 'rhombic dodecahedron (xy-hexagon) "a" and "b" vectors' and 'cubic "c" vector'
            if len(cmd) == 1:
                self.pbcx = cmd[0]*10
                self.pbcy = math.sqrt(3)*self.pbcx/2
                self.pbcz = cmd[0]*10
            elif len(cmd) == 2:
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
            
        ### Following math taken from insane.py and https://manual.gromacs.org/current/reference-manual/algorithms/periodic-boundary-conditions.html
        elif self.pbc_type == "dodecahedron":
            ### Not actually a dodecahedron, but instead a parallelepiped constituting a third of it
            ### Uses 'rhombic dodecahedron (xy-hexagon) vectors'
            if len(cmd) == 1:
                if momentary_z:
                    ### Done when "stacked_membranes" is used as the box size is defined by the z-height
                    self.pbcz = cmd[0]*10
                    self.pbcx = self.pbcz/math.sqrt(6)*3
                    self.pbcy = math.sqrt(3)*self.pbcx/2
                else:
                    self.pbcx = cmd[0]*10
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
        
        if not self.output_system_pdb_file_name and not self.output_system_gro_file_name:
            self.output_system_pdb_file_name = "output.pdb"
            self.output_system_gro_file_name = "output.gro"
    
    def run(self, kwargs):
        '''
        Runs the entire system creation process
        '''
        
        CGSB_run_tic = time.time()
        self.commands_handler(kwargs)

        ### Initial checks
#         assert len(self.pbc_box) > 0, "Box dimensions not set. Please do so using 'box=[x,y,z]'"
#         assert len(self.pbc_box) == 3, "Box dimensions improperly defined. 3 dimensions must be given"
#         assert all([self.is_number(ax)[0] for ax in self.pbc_box]), "Not all box values are numbers"

        assert any([len(cmd) > 0 for cmd in [self.PROTEINS_cmds, self.MEMBRANES_cmds, self.SOLVATIONS_cmds]]), (
            "Running requires at least one command of one of the following types: 'protein', 'membrane' or 'solvation'"
        )

#         self.print_term("------------------------------ PREPROCESSING DEFINITIONS", spaces=0, verbose=1)
#         ### Definition preprocessing
#         preprocessing_tic = time.time()
#         self.itp_read_initiater()
#         self.import_structures_handler()
#         self.lipid_defs_preprocessor()
#         self.solvent_defs_preprocessor()
#         self.ion_defs_preprocessor()
#         preprocessing_time = round(preprocessing_toc - preprocessing_tic, 4)
#         self.print_term("------------------------------ DEFINITIONS PREPROCESSING COMPLETE", "(time spent: "+str(preprocessing_time)+" [s])", "\n", spaces=0, verbose=1)

        self.print_term("------------------------------ PREPROCESSING COMMANDS", spaces=0, verbose=1)
        ### Command preprocessing
        preprocessing_tic = time.time()
        self.prot_preprocessor()
        self.memb_preprocessor(return_self = "self", cmds_given = False)
        self.solv_preprocessor()
        preprocessing_toc = time.time()
        preprocessing_time = round(preprocessing_toc - preprocessing_tic, 4)
        self.print_term("------------------------------ PREPROCESSING COMPLETE", "(time spent: "+str(preprocessing_time)+" [s])", "\n", spaces=0, verbose=1)

        ### Run the program
        self.prot_placer()
        self.subleaflet_poly_maker()
        self.holed_subleaflet_bbox_maker()
        self.lipid_calculator()
        self.planar_grid_maker()
        self.lipid_inserter()
        self.solvater()

        self.print_term("--------------------", verbose=1)
        self.print_term("Final system charge:", self.system_charge, verbose=1)
        self.print_term("--------------------", "\n", verbose=1)

        self.pickler()

        ### Write the files
        self.system_file_writer()
        self.topol_file_writer()
        self.log_file_writer()

        self.print_term("My task is complete. Did i do a good job?", verbose=1)
        CGSB_run_toc  = time.time()
        CGSB_run_time = round(CGSB_run_toc - CGSB_run_tic, 4)
        self.print_term("Time spent running CGSB:", CGSB_run_time, verbose=1)
    
    #####################################
    ### Specialized printing function ###
    #####################################
    def print_term(self, *string, spaces=0, verbose=0, space="    ", sep=" ", end="\n", debug = False, extra = False, warn = False):
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
        if print_true and not self.quiet and self.verbose >= verbose:
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
        self.print_term("Writing structure", files_string, "(" + "/".join(out_stringlist) + ")", verbose=1)
            
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
            ]
        if self.output_system_gro_file_name:
            
            output_system_gro_file_lines = [
                self.system_name,
                "PLACEHOLDER_ATOM_COUNT",
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
                for (i, atom, res), (a_name, beadnr, x, y, z, r_name, resnr, charge) in zip(original_bead_info, bead_vals):
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
            self.print_term("---- PDB file written:", self.output_system_pdb_file_name, verbose=1)
        if self.output_system_gro_file_name:
            new_file = open(self.output_system_gro_file_name, "w")
            for line in output_system_gro_file_lines:
                new_file.write(line + "\n")
            new_file.close()
            self.print_term("---- GRO file written:", self.output_system_gro_file_name, verbose=1)
        self.print_term("", verbose=1)

    def topol_file_writer(self):
        if self.output_topol_file_name:
            self.print_term("Writing topology file:", self.output_topol_file_name, verbose=1)
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
            self.print_term("---- Topology file written", "\n", verbose=1)
    
    def log_file_writer(self):
        if self.output_log_file_name:
            self.print_term("Writing log file file:", self.output_log_file_name, verbose=1)
            if self.backup:
                self.backupper(self.output_log_file_name)
            new_file = open(self.output_log_file_name, "w")
            for line in self.LOG_FILE:
                new_file.write(line)
            new_file.close()
            self.print_term("---- Log file written", "\n", verbose=1)
    
    ###################################################
    ### Lipid/solvent/ion definitions preprocessors ###
    ###################################################
    def lipid_defs_preprocessor(self):
        '''
        Preprocesses lipid defintions, by rearranging the data into a class
        '''
        self.print_term("Preprocessing lipid definitions", verbose=2)
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
            
            ### Checks what specifiers have been given for each lipid name
            ### Used to check if lipid-specific charges have been given
            dict_keys, dict_vals = zip(*lipid_type_dict["lipids"].items())
            lipid_names, data_specifiers = zip(*dict_keys)
            lipid_name_data_specifiers = {lipid_name: {} for lipid_name in lipid_names}
            for (lipid_name, data_specifier), lipid_details in lipid_type_dict["lipids"].items():
                if data_specifier not in lipid_name_data_specifiers[lipid_name]:
                    lipid_name_data_specifiers[lipid_name][data_specifier] = lipid_details
            
            ### Removing duplicates from lipid_names list while keeping the order
            lipid_names = list(dict.fromkeys(lipid_names))
            
            for lipid_name in lipid_names:
        
                if lipid_name not in self.lipid_dict[params].keys():
                    self.lipid_dict[params][lipid_name] = LIPID(hydr_z = hydr_z, molname = lipid_name)
                
                ### Checks if lipid-specific charges are given and joins them with lipid_type charges
                if "charges" in lipid_name_data_specifiers[lipid_name]:
                    lipid_details = lipid_type_dict["lipids"][(lipid_name, "charges")]
                    specifc_charge_names = zip(*lipid_details)[0]
                    lipid_charge_dict = lipid_details + tuple([(key, val) for key, val in lipid_type_charges if key not in specifc_charge_names])
                else:
                    lipid_charge_dict = lipid_type_charges
            
                if "beads" in lipid_name_data_specifiers[lipid_name]:
                    lipid_details = lipid_type_dict["lipids"][(lipid_name, "beads")]
                    beads = list(filter(None, lipid_details.split(" ")))
                    ### Remove coordinate and bead indexes if no bead is assigned
                    ### x and z are centered in the plane using "set_xy_to_center()" method
                    ### z centering according to hydrophobic height "hydr_z"
                    lx = [xi * bdx * 10 for xi, bead in zip(x, beads) if bead != "-"]
                    ly = [yi * bdy * 10 for yi, bead in zip(y, beads) if bead != "-"]
                    lz = [(zi - hydr_z) * bdz * 10 for zi, bead in zip(z, beads) if bead != "-"]
                    beads = [bead for bead in beads if bead != "-"]
                    charges = []
                    beadnrs = list([i for i in range(len(beads))])
                    
                    if lipid_charge_dict:
                        charge_beads, charge_vals = zip(*lipid_type_charges)
                        for bead in beads:
                            if bead in charge_beads:
                                charge_index = charge_beads.index(bead)
                                charges.append(charge_vals[charge_index])
                            else:
                                charges.append(0)
                    else:
                        charges = [0 for _ in range(len(beads))]
                    
                    self.lipid_dict[params][lipid_name].add_res_and_beads(
                        lipid_name,
                        beads = beads,
                        beadnrs = beadnrs,
                        xs = lx,
                        ys = ly,
                        zs = lz,
                        charges = charges,
                    )
                    self.lipid_dict[params][lipid_name].set_xy_to_center()
                    
        tot_lipids = sum([len(vals) for vals in self.lipid_dict.values()])
        
        self.print_term("Number of lipids preprocessed:", tot_lipids, "\n", spaces=1, verbose=3)
    
    def solute_beads_checker(self, res_dict, cur_name, resname):
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
            "molecule name: " + cur_name + "\n"
            "name of problematic residue: " + resname + "\n"
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
        charge_dict = {}
        tot_charge = False
        moleculetype = False
        
        if "mapping_ratio" in cur_dict.keys():
            mol_class.mapping_ratio_set(cur_dict["mapping_ratio"])

        if "density" in cur_dict.keys():
            mol_class.density_set(cur_dict["density"])

        if "molar_mass" in cur_dict.keys():
            mol_class.molar_mass_set(cur_dict["molar_mass"])
        
        if "residues" not in cur_dict.keys():
            residues_dict = [{
                    "resname": cur_name,
                    "beads": cur_dict["beads"],
                    "x": cur_dict["x"],
                    "y": cur_dict["y"],
                    "z": cur_dict["z"],
            }]
            if "charges" in cur_dict:
                if type(cur_dict["charges"][0]) == str:
                    cur_dict["charges"] = (cur_dict["charges"],)
                residues_dict[0]["charges"] = tuple([(bead, charge) for bead, charge in cur_dict["charges"]])

        elif type(cur_dict["residues"]) == dict:
            residues_dict = [cur_dict["residues"]]
        elif isinstance(cur_dict["residues"], (list, tuple)):
            residues_dict = cur_dict["residues"]
        
        ### Molecule-wide charges
        if "charges" in cur_dict.keys() and "residues" in cur_dict.keys():
            if type(cur_dict["charges"][0]) == str:
                cur_dict["charges"] = (cur_dict["charges"],)
            if "residues" in cur_dict.keys():
                ress, beads, charges = zip(*cur_dict["charges"])
            else:
                ress = [cur_name]
                beads, charges = zip(*cur_dict["charges"])
            for res, bead, charge in zip(ress, beads, charges):
                if res not in res_charge_dict.keys():
                    charge_dict[res] = {}
                charge_dict[res][bead] = charge
        
        ### Residue-specific charges
        if "residues" in cur_dict.keys():
            for res_dict in residues_dict:
                resname = res_dict["resname"]
                if "charges" in res_dict:
                    if type(res_dict["charges"][0]) == str:
                        res_dict["charges"] = (res_dict["charges"],)
                    for bead, charge in res_dict["charges"]:
                        charge_dict[resname][bead] = charge
        
        ### Total charge.
        ### Spreads charge across all beads. Should only be used if no bead-specfic charges are given.
        ### Ideally only use this for single bead solutes/solvents where the bead is obvious
        if "charge" in cur_dict.keys():
            number = self.get_number_from_string(cur_dict["charge"])
            if number is not False:
                tot_charge = number
            else:
                moleculetype = cur_dict["charge"]
        
        ### Dictionaries without "residues" key have had the beads, and positions converted to "residues" key:val pair
        ### Thus "residues" will always be present
        molbeadnrs = []
        for res_dict in residues_dict:
            resname = res_dict["resname"]
            beads, xs, ys, zs = self.solute_beads_checker(res_dict, cur_name, resname)

            if molbeadnrs == []:
                beadnrs = list([i for i in range(len(beads))])
                molbeadnrs = beadnrs[:]
            else:
                beadnrs = list([i for i in range(max(molbeadnrs)+1, len(molbeadnrs)-1+len(beads))])
            
            if tot_charge is not False:
                charges = [tot_charge/len(beads) for _ in range(len(beads))]
            elif moleculetype is not False:
                ### If topology is given
                moleculetype = tot_charge
                ### Set charges to be 0 for now, since they will be obtained from topology later
                charges = [0 for _ in range(len(beads))]
            else:
                charges = []
                for bead in beads:
                    if resname in charge_dict:
                        if bead in charge_dict[resname]:
                            charges.append(charge_dict[resname][bead])
                        else:
                            charges.append(0)
                    else:
                        charges.append(0)

            mol_class.add_res_and_beads(
                resname = resname, beads = beads, beadnrs = beadnrs,
                xs = xs, ys = ys, zs = zs,
                charges = charges,
            )
        
        if moleculetype:
            ### If explicit topology name is given
            mol_class.moleculetype = moleculetype
        else:
            ### Else assume cur_name to be topology name
            mol_class.moleculetype = cur_name
        
        mol_class.set_coords_to_center()
    
    def solvent_defs_preprocessor(self):
        '''
        Preprocesses solvent defintions by rearranging the data into a class
        '''
        self.print_term("Preprocessing solvent definitions", verbose=2)
        for params, solvent_type_dict in self.solvent_defs.items():
            if params not in self.solvent_dict.keys():
                self.solvent_dict[params] = {}
            ### Loop over indexes
            for cur_name, cur_dict in solvent_type_dict.items():
                if cur_name not in self.solvent_dict[params].keys():
                    self.solvent_dict[params][cur_name] = SOLVENT(mapping_ratio = 1, molname = cur_name)
                else:
                    continue
                
                self.solute_defs_checker(self.solvent_dict[params][cur_name], cur_name, cur_dict)
                
        tot_solvents = sum([len(vals) for vals in self.solvent_dict.values()])
        self.print_term("Number of solvents preprocessed:", tot_solvents, "\n", spaces=1, verbose=3)

    def ion_defs_preprocessor(self):
        '''
        Preprocesses ion defintions by rearranging the data into a class
        '''
        self.print_term("Preprocessing ion definitions", verbose=2)
        ### Creates parameter libraries in the solvent defs for positive and negative ions
        ### In case people want to flood a system with a specific number of ions
        if "pos_ions" not in self.solvent_defs.keys():
            self.solvent_defs["pos_ions"] = {}
        if "neg_ions" not in self.solvent_defs.keys():
            self.solvent_defs["neg_ions"] = {}
            
        for params, charge_type_dict in self.ion_defs.items():
            if params not in self.ion_dict.keys():
                self.ion_dict[params] = {}
            for charge_name, ion_type_dict in charge_type_dict.items():
                if charge_name not in self.ion_dict[params].keys():
                    self.ion_dict[params][charge_name] = {}
                ### Loop over indexes
                for cur_name, cur_dict in ion_type_dict.items():
                    if cur_name not in self.ion_dict[params][charge_name].keys():
                        ion = SOLVENT(mapping_ratio = 1, molname = cur_name)
                        self.ion_dict[params][charge_name][cur_name] = ion
                        if params == "default":
                            if charge_name not in self.solvent_defs["pos_ions"].keys():
                                self.solvent_defs["pos_ions"][cur_name] = ion
                            if charge_name not in self.solvent_defs["neg_ions"].keys():
                                self.solvent_defs["neg_ions"][cur_name] = ion
                    else:
                        continue
                    
                    self.solute_defs_checker(self.ion_dict[params][charge_name][cur_name], cur_name, cur_dict)
                    
        tot_pos_ions = sum([len(vals["positive"]) for vals in self.ion_dict.values()])
        self.print_term("Number of positive ions preprocessed:", tot_pos_ions, spaces=1, verbose=3)
        tot_neg_ions = sum([len(vals["negative"]) for vals in self.ion_dict.values()])
        self.print_term("Number of negative ions preprocessed:", tot_neg_ions, spaces=1, verbose=3)
    
    ####################################
    ### Special system preprocessors ###
    ####################################
    def stacked_membranes_preprocessor(self):
        '''
        Preprocesses stacked membrane commands for later ease of use
        '''
        box_heights = []
        if len(self.STACKED_MEMBRANES_cmds) > 0:
            self.print_term("Preprocessing stacked membrane commands", verbose=2)
            for cmd_nr, sm_cmd in enumerate(self.STACKED_MEMBRANES_cmds, 1):
                self.print_term("Starting command:", cmd_nr, spaces=1, verbose=3)

                ### Defaults
                settings_dict = {
                    "dz": [5], # Can be set for each individual space
                    "dtype": ["surface"], # Can be set for each individual space
                    "dn": 0,
                    "membrane_commands": [],
                    "solvation_commands": [],
                }
                current_sub_cmd = False
                memb_i = -1
                solv_i = -1
                
                dn_lipid_sizes = []
                
                if type(sm_cmd) == str:
                    split_cmds = sm_cmd.split()
                elif isinstance(sm_cmd, list, tuple):
                    split_cmds = " ".join(sm_cmd).split()
                
                ### ### Check protein command
                for cmd in split_cmds:
                    sub_cmd = cmd.split(":")
                    
                    if not current_sub_cmd and sub_cmd[0] not in ["membrane", "memb", "solvation", "solv"]:
                        if sub_cmd[0].lower() == "dz":
                            settings_dict["dz"] = [self.get_number_from_string(dz) for dz in sub_cmd[1:]]

                        elif sub_cmd[0].lower() == "dtype":
                            if sub_cmd[1] in ["surface", "center", "mean_surface"]:
                                settings_dict["dtype"] = [d_type for d_type in sub_cmd[1:]]

                        elif sub_cmd[0].lower() == "dn":
                            settings_dict["dn"] = int(self.get_number_from_string(sub_cmd[1]))
                    
                    elif sub_cmd[0].lower() in ["membrane", "memb"] and len(sub_cmd) == 1:
                        settings_dict["membrane_commands"].append({"subcommands": []})
                        memb_i += 1
                        current_sub_cmd = "membrane"
                    
                    elif sub_cmd[0].lower() in ["solvation", "solv"] and len(sub_cmd) == 1:
                        settings_dict["solvation_commands"].append({"subcommands": []})
                        solv_i += 1
                        current_sub_cmd = "solvation"
                    
                    else:
                        if current_sub_cmd == "membrane":
                            if sub_cmd[0] == "dn":
                                settings_dict["membrane_commands"][memb_i]["dn"] = [
                                    int(self.get_number_from_string(dn)) for dn in sub_cmd[1:]
                                ]
                            else:
                                settings_dict["membrane_commands"][memb_i]["subcommands"].append(cmd)
                                
                        elif current_sub_cmd == "solvation":
                            if sub_cmd[0] == "dn":
                                settings_dict["solvation_commands"][solv_i]["dn"] = [
                                    int(self.get_number_from_string(dn)) for dn in sub_cmd[1:]
                                ]
                            else:
                                settings_dict["solvation_commands"][solv_i]["subcommands"].append(cmd)
                
                ### Sets command "dn" if no value has been explicitly given
                if settings_dict["dn"] == 0:
                    settings_dict["dn"] = max(len(settings_dict["membrane_commands"]), len(settings_dict["solvation_commands"]))
                
                ### Checks if mulitple "dz" values have been given
                if len(settings_dict["dz"]) == 1:
                    settings_dict["dz"] = [settings_dict["dz"][0] for dn in range(settings_dict["dn"])]
                elif len(settings_dict["dz"]) > 1 and len(settings_dict["dz"]) != settings_dict["dn"]:
                    assert False, (
                        "Multiple 'dz' values ("+str(len(settings_dict["dz"]))+") have been given but the amount does not equal 'dn' ("+str(settings_dict["dn"])+")",
                        "I you want to provide multiple 'dz' values, please ensure that you provide a number of values equal to 'dn'"
                    )
                
                ### converts all "dz" values from [nm] to [Å]
                settings_dict["dz"] = [dz*10 for dz in settings_dict["dz"]]
                
                ### Checks if mulitple "dtype" values have been given
                if len(settings_dict["dtype"]) == 1:
                    settings_dict["dtype"] = [settings_dict["dtype"][0] for dn in range(settings_dict["dn"])]
                elif len(settings_dict["dtype"]) > 1 and len(settings_dict["dtype"]) != settings_dict["dn"]:
                    assert False, (
                        "Multiple 'd_type' values ("+str(len(settings_dict["dtype"]))+") have been given but the amount does not equal 'dn' ("+str(settings_dict["dn"])+")",
                        "I you want to provide multiple 'd_type' values, please ensure that you provide a number of values equal to 'dn'"
                    )
                
                ### Checks largest "dn" value used for membranes
                ### Also runs the membrane commands through the memb_preprocessor method
                memb_max_dn = 0
                for memb_i, memb_dict in enumerate(settings_dict["membrane_commands"]):
                    if "dn" in memb_dict.keys():
                        memb_max_dn = max(memb_max_dn, max(memb_dict["dn"]))
                    memb_dict["preprocessed"] = self.memb_preprocessor(return_self = "return", cmds_given = [" ".join(memb_dict["subcommands"])])
                
                if memb_max_dn > settings_dict["dn"]:
                    self.print_term(
                        "Largest 'dn' value used for membranes ("+str(memb_max_dn)+") is larger than stacked_membranes 'dn' ("+str(settings_dict["dn"])+").",
                        "Will ignore membranes using 'dn' values larger than ("+str(settings_dict["dn"])+")",
                        warn=True
                    )
                
                ### Checks largest "dn" value used for solvations
                solv_max_dn = 0
                for solv_i, solv_dict in enumerate(settings_dict["solvation_commands"]):
                    if "dn" in solv_dict.keys():
                        solv_max_dn = max(solv_max_dn, max(solv_dict["dn"]))
                
                if solv_max_dn > settings_dict["dn"]:
                    self.print_term(
                        "Largest 'dn' value used for solvations ("+str(solv_max_dn)+") is larger than command 'dn' ("+str(settings_dict["dn"])+").",
                        "Will ignore solvations using 'dn' values larger than ("+str(settings_dict["dn"])+")",
                        warn=True
                    )
                
                ### Makes sure all "dn" values have been used
                dn_memb_pointers = {}
                dn_solv_pointers = {}
                for dn in range(1, settings_dict["dn"] + 1):
                    for memb_i, memb_dict in enumerate(settings_dict["membrane_commands"]):
                        if "dn" in memb_dict.keys():
                            if dn in memb_dict["dn"]:
                                if dn not in dn_memb_pointers.keys():
                                    dn_memb_pointers[dn] = memb_i
                                else:
                                    self.print_term(
                                        "One 'dn' value ("+str(dn)+") has been used for multiple membranes.",
                                        "Make sure to only use a 'dn' value once.",
                                        warn=True
                                    )
                                    
                    for solv_i, solv_dict in enumerate(settings_dict["solvation_commands"]):
                        if "dn" in solv_dict.keys():
                            if dn in solv_dict["dn"]:
                                if dn not in dn_solv_pointers.keys():
                                    dn_solv_pointers[dn] = solv_i
                                else:
                                    self.print_term(
                                        "One 'dn' value ("+str(dn)+") has been used for multiple solvations.",
                                        "Make sure to only use a 'dn' value once.",
                                        warn=True
                                    )
                
                ### Checks if number of given "dn" pointers are zero or matches "dn" setting for command
                if not (len(dn_memb_pointers) == 0 or len(dn_memb_pointers) == settings_dict["dn"]):
                    assert False, (
                        "Number of membrane 'dn' values ("+str(len(dn_memb_pointers))+") must be either 0 or equal the 'dn' value of the command ("+str(settings_dict["dn"])+")."
                    )
                
                if not (len(dn_solv_pointers) == 0 or len(dn_solv_pointers) == settings_dict["dn"]):
                    assert False, (
                        "Number of solvation 'dn' values ("+str(len(dn_solv_pointers))+") must be either 0 or equal the 'dn' value of the command ("+str(settings_dict["dn"])+")."
                    )
                
                ### Sets "dn" pointers if they have not been specified
                if len(dn_memb_pointers) == 0 and len(settings_dict["membrane_commands"]) > 0:
                    if len(settings_dict["membrane_commands"]) == settings_dict["dn"]:
                        ### If multiple membrane commands
                        for i in range(settings_dict["dn"]):
                            settings_dict["membrane_commands"][i]["dn"] = [i+1]
                            dn_memb_pointers[i+1] = i
                    else:
                        ### If single membrane command
                        settings_dict["membrane_commands"][0]["dn"] = [i for i in range(settings_dict["dn"])]
                        dn_memb_pointers = {i+1: 0 for i in range(settings_dict["dn"])}
                
                if len(dn_solv_pointers) == 0 and len(settings_dict["solvation_commands"]) > 0:
                    if len(settings_dict["solvation_commands"]) == settings_dict["dn"]:
                        ### If multiple solvation commands
                        for i in range(settings_dict["dn"]):
                            settings_dict["solvation_commands"][i]["dn"] = [i+1]
                            dn_solv_pointers[i+1] = i
                    else:
                        ### If single solvation command
                        settings_dict["solvation_commands"][0]["dn"] = [i for i in range(settings_dict["dn"])]
                        dn_solv_pointers = {i+1: 0 for i in range(settings_dict["dn"])}
                
                ### Sorting according to dn value just to be sure. Probably unnecessary
                dn_memb_pointers = dict(sorted(dn_memb_pointers.items(), key=lambda x: x[0]))
                dn_solv_pointers = dict(sorted(dn_solv_pointers.items(), key=lambda x: x[0]))
                
                ### Calculates intermembrane spacing
                membrane_centers = []
                solvation_delimiters = []
                cur_center = 0
                box_height = 0
                
                ### Creates membrane commands and finds solvent command z-delimiters
                for dni, dn in enumerate(dn_memb_pointers.keys()):
                    
                    dz    = settings_dict["dz"][dni]
                    dtype = settings_dict["dtype"][dni]
                    
                    if dn-1 in dn_memb_pointers:
                        upper_dn = dn-1
                    else:
                        upper_dn = max(dn_memb_pointers.keys())
                    lower_dn = dn
                    
                    upper_memb = settings_dict["membrane_commands"][dn_memb_pointers[upper_dn]]["preprocessed"]
                    lower_memb = settings_dict["membrane_commands"][dn_memb_pointers[lower_dn]]["preprocessed"]
                    
                    if dtype == "center":
                        upper_space = 0
                        lower_space = 0
                    elif dtype == "surface":
                        upper_space = abs(upper_memb["bead_maxz"])
                        lower_space = abs(lower_memb["bead_minz"])
                    elif dtype == "mean_surface":
                        upper_space = abs(upper_memb["bead_mean_maxz"])
                        lower_space = abs(lower_memb["bead_mean_minz"])
                    
                    ### (total distance, upper solvent space, lower solvent space)
                    upper_space = dz/2 + upper_space
                    lower_space = dz/2 + lower_space
                    total_space = upper_space + lower_space
                    
                    if dni == 0:
                        first_dn_lower_delimiters = (cur_center, cur_center + lower_space)
                        last_space = round(upper_space, 3)
                        cur_center += round(lower_space, 3)
                    else:
                        lower_delimiters = (cur_center, cur_center + lower_space)
                        cur_center += round(lower_space, 3)
                        upper_delimiters = (cur_center, cur_center + upper_space)
                        cur_center += round(upper_space, 3)
                        solvation_delimiters.append((lower_delimiters, upper_delimiters))
                    box_height += total_space
                    membrane_centers.append(cur_center)
                
                first_dn_upper_delimiters = (cur_center, cur_center + last_space)
                solvation_delimiters = [(first_dn_lower_delimiters, first_dn_upper_delimiters)] + solvation_delimiters
                
                membrane_commands = []
                for dni, (dn, memb_pointer) in enumerate(dn_memb_pointers.items()):
                    membrane_command = settings_dict["membrane_commands"][memb_pointer]["subcommands"].copy()
                    membrane_center = round((membrane_centers[dni] - box_height/2)/10, 3)
                    membrane_command.append("cz:" + str(membrane_center))
                    membrane_commands.append(" ".join(membrane_command))
                
                ### Creates solvation commands
                solvation_commands = []
                for dni, (dn, solv_pointer) in enumerate(dn_solv_pointers.items()):
                    main_solvation_command = settings_dict["solvation_commands"][solv_pointer]["subcommands"].copy()
                    
                    for bot_delimiter, top_delimiter in solvation_delimiters[dni]:
                        cz = round((np.mean([top_delimiter, bot_delimiter]) - box_height/2)/10, 3)
                        z = round((top_delimiter - bot_delimiter)/10, 3)
                        solvation_command = main_solvation_command.copy()
                        solvation_command.append("cz:" + str(cz))
                        solvation_command.append("z:" + str(z))
                        solvation_commands.append(" ".join(solvation_command))
                
                self.MEMBRANES_cmds.extend(membrane_commands)
                self.SOLVATIONS_cmds.extend(solvation_commands)
                box_heights.append(round(box_height/10, 3))
        
        return max(box_heights)
    
    ##############################################
    ### Protein/membrane/solvent preprocessors ###
    ##############################################
    def prot_preprocessor(self):
        '''
        Preprocesses protein commands for later ease of use
        '''
        if len(self.PROTEINS_cmds) != 0:
            self.print_term("Preprocessing protein requests", verbose=2)
            for cmd_nr, prot_cmd in enumerate(self.PROTEINS_cmds, 1):
                self.print_term("Starting protein:", cmd_nr, spaces=1, verbose=3)

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
                    ### Buffer used during solvations
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
                        prot_dict["lipids_inside"] = ast.literal_eval(str(sub_cmd[1]))

                    ### True/False whether to check for pbc crossing
                    elif sub_cmd[0].lower() == "pbc_check":
                        prot_dict["pbc_check"] = ast.literal_eval(str(sub_cmd[1]))

                    ### True/False whether to check for pbc crossing
                    elif sub_cmd[0].lower() == "buffer":
                        prot_dict["buffer"] = ast.literal_eval(str(sub_cmd[1]))

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
                            prot_dict["charge"] = ast.literal_eval(str(sub_cmd[1]))

                    ### xyz axis translations [nm]
                    elif sub_cmd[0].lower() in ["tx", "ty", "tz", "cx", "cy", "cz"]:
                        if sub_cmd[0].lower().startswith("c"):
                            sub_cmd[0] = "t" + sub_cmd[0].lower()[-1]
                        prot_dict[sub_cmd[0].lower()] = ast.literal_eval(str(sub_cmd[1])) * 10

                    ### xyz axis rotations [degrees]
                    elif sub_cmd[0].lower() in ["rx", "ry", "rz"]:
                        prot_dict[sub_cmd[0].lower()] = ast.literal_eval(str(sub_cmd[1]))

                    ### Random kick to beads x/y/z positions [nm] #[Å]
                    elif sub_cmd[0].lower() in ["kick", "kickx", "kicky", "kickz"]:
                        if sub_cmd[0].lower() == "kick":
                            prot_dict["kickx"] = ast.literal_eval(str(sub_cmd[1]))
                            prot_dict["kicky"] = ast.literal_eval(str(sub_cmd[1]))
                            prot_dict["kickz"] = ast.literal_eval(str(sub_cmd[1]))
                        else:
                            prot_dict[sub_cmd[0].lower()] = ast.literal_eval(str(sub_cmd[1]))

                    ### Errors out if unknown subcommand used, and prints the subcommand to console
                    else:
                        assert False, "Unknown subcommand given to '-prot'. The subcommand is: '" + str(cmd) + "'"
                
                def prot_charge_finder(beads):
                    '''
                    Finds the charge of a protein based on residue names and the 'prot_defs' dictionary
                    '''
                    protein_bead_charges = []
                    for (bead_i, bead_nr, res_nr), values in beads.items():
                        params   = self.prot_params or self.sys_params
                        res_name  = values["res_name"]
                        atom_name = values["atom_name"]
                        if res_name in self.prot_defs[params]["charges"]:
                            if atom_name in self.prot_defs[params]["charges"][res_name]:
                                protein_bead_charges.append(self.prot_defs[params]["charges"][res_name][atom_name])
                            else:
                                self.print_term(atom_name, "atom name not in 'prot_defs' charges for residue ", res_name, " for system parameters:", params, warn=True)
                                protein_bead_charges.append(0)
                        else:
                            self.print_term(res_name, "residue name not in 'prot_defs' charges for system parameters:", params, warn=True)
                            protein_bead_charges.append(0)
                    return protein_bead_charges
                
                ### Post-preprocessing (topology and charge determination)
                if isinstance(prot_dict["charge"], (float, int)):
                    ### Sets charge manually
                    ### Evenly distributes charges accross all beads as there is no way of knowing where they are located based on a single charge value
                    protein_bead_charges = [prot_dict["charge"]/len(prot_dict["beads"]) for _ in range(len(prot_dict["beads"]))]
                    
                elif prot_dict["charge"] == "auto":
                    ### Finds charge information from prot_defs charge dictionary
                    protein_bead_charges = prot_charge_finder(prot_dict["beads"])
                
                elif prot_dict["charge"] == "top":
                    ### Finds charge information in topology files
                    if prot_dict["mol_names"] == []:
                        ### Reverts to getting charge data from prot_defs charge dictionary if molecule names not given
                        self.print_term("No protein names given. Will estimate charges from residue names. Use 'mol_name' to assign protein names (name used in topology files)", spaces=2, warn=True)
                        protein_bead_charges = prot_charge_finder(prot_dict["beads"])
                    else:
                        ### Finds charge data from topology files
                        protein_bead_charges = []
                        erase_charges = False
                        for mol_name in prot_dict["mol_names"]:
                            if mol_name not in self.itp_moltypes.keys():
                                self.print_term("A molecule name could not be found in your topology file(s): " +  mol_name, warn=True, spaces=2)
                                self.print_term("Setting all molecule bead charges to 0", warn=True, spaces=2)
                                erase_charges = True
                            else:
                                protein_bead_charges.extend([
                                    atom["charge"] for atom in self.itp_moltypes[mol_name]["atoms"].values()
                                ])
                        if erase_charges:
                            protein_bead_charges = [0 for _ in range(len(prot_dict["beads"]))]
                
                prot_dict["mol_names"] = prot_dict["mol_names"] or ["_".join(["PROT", str(cmd_nr)])]
                
                ### ### Mapping charge values to int/float data type
                ### First ensure values are strings otherwise ast.literal_eval will error
                protein_bead_charges = list(map(str, protein_bead_charges))
                ### The map strings to int/float
                protein_bead_charges = list(map(ast.literal_eval, protein_bead_charges))
                
                ### Converting data into protein class
                prot_dict["protein"] = PROTEIN(molname = prot_dict["mol_names"])
                
                cur_res = False
                for (key, vals), charge in zip(prot_dict["beads"].items(), protein_bead_charges):
                    if vals["res_nr"] != cur_res:
                        prot_dict["protein"].add_res(vals["res_name"], vals["res_nr"])
                        cur_res = vals["res_nr"]
                    
                    prot_dict["protein"].add_bead(
                        vals["atom_name"],
                        vals["atom_nr"],
                        vals["x"],
                        vals["y"],
                        vals["z"],
                        charge=charge,
                    )
                
                self.system_charge += prot_dict["protein"].get_mol_charge()
                
                self.PROTEINS[cmd_nr] = prot_dict.copy()
            self.print_term("Number of molecule insertions preprocessed:", len(self.PROTEINS), spaces=1, verbose=3)

    def memb_preprocessor(self, return_self = "self", cmds_given = False):
        '''
        Preprocesses membrane commands for later ease of use
        '''
        ### ### Preprocessing membranes
        
        if return_self == "return":
            MEMBRANES = {}
        
        if cmds_given:
            MEMBRANES_cmds = cmds_given
        else:
            MEMBRANES_cmds = self.MEMBRANES_cmds
            
        if len(MEMBRANES_cmds) != 0:
            self.print_term("Preprocessing membrane requests", verbose=2)
            for cmd_nr, leaf_cmd in enumerate(MEMBRANES_cmds, 1):
                self.print_term("Starting membrane command:", cmd_nr, spaces=1, verbose=3)
                ### "membrane" settings always apply to the whole membrane instead of individual leaflets
                ### "default" settings can be overwritten by individual leaflets
                
                settings_dict = {
                    "upper_leaf": {},
                    "lower_leaf": {},
                    "membrane"  : {
                        "x": self.pbcx / 10, # [nm] converted to [Å]
                        "y": self.pbcy / 10, # [nm] converted to [Å]
                        "center": [0, 0, 0], # [nm] converted to [Å]
                        
                        "gridsplits": ("auto", 500), # Å
                        
                        "pbc_check": True, # [bool]
                        
                        "lipid_distribution": "random",
                        
                        "optimize": "v5", # False/"no"/0, "limited", True/"full"/1
                        "optim_maxsteps": 100,
                        "optim_push_tol": 0.5,
                        "optim_push_mult": 1.0,
                        
                        ### Readjusts the max/min x/y-values if cutouts/holes are made
                        "readjust_bbox": True,
                        ### Readjusts the max/min x/y-values if cutouts/holes are made
                        "split_bbox": True,
                        ### Allows empty space (holes/cutouts) to be solvated
                        "solvate_empty": True,
                    },
                    "default": {
                        "lipids": {}, # empty
                        "lipids_preprocessing": [], # empty

                        "kickx": 0.025, # [nm] converted to [Å]
                        "kicky": 0.025, # [nm] converted to [Å]
                        "kickz": 0.025, # [nm] converted to [Å]

                        "apl": 0.6, # [nm^2] converted to [Å^2]
                        
                        ### 0.264 = vdw of regular beads. 0.264/16*6 = 0.099, 0.264/4=0.066, # 0.264/2=0.132
                        "plane_buffer": 0.264/16*6, # [nm] converted to [Å]
                        "height_buffer": 0.264/16*6, # [nm] converted to [Å]

                        "prot_buffer": 0.132, # [nm] converted to [Å], default = (vdw of regular beads) / 2
                        "alpha_mult": 1.0,
                        
                        "lip_round_func": ("int", int), # int, round, math.floor, math.ceil

                        "lipid_optim": 'avg_optimal', # str: see beloww
                        ### 'avg_optimal': Ratios are optimised based on average lipid distance from optimal composition
                        ### 'abs_val':     Ratio-values are treated as absolute number of lipids
                        ### 'no':          No optimization.
                        ### 'insane':      'no' and initial lipid calculation is intentionally wrong like in insane
                        ### 'force_fill':  Fills up to the allowed number number of lipids
                        ### 'fill:         Same as 'force_fill' but stops if perfect lipid ratio has been achieved
                        "params": False, # False or str
                        "charge": "top", # "lib" or "top"
                        
                        "holes": {
                            "ellipses": [],
                            "polygons": [],
                        },
                    }
                }
                
                hole_types = ["circle", "ellipse", "square", "rectangle", "poly", "polygon"]
                
                ellipse_abbrevs = {}
                ellipse_abbrevs_key_vals = [
                    (("rot", "rotate", "rotation"), "rotation"),
                    (("b", "buf", "buffer"),        "buffer"),
                    (("inverse",),                  "inverse"),
                ]
                for keys, val in ellipse_abbrevs_key_vals:
                    ellipse_abbrevs.update(dict.fromkeys(keys, val))
                
                polygon_abbrevs = {}
                polygon_abbrevs_key_vals = [
                    (("rot", "rotate", "rotation"), "rotation"),
                    (("b", "buf", "buffer"),        "buffer"),
                    (("p", "point"),                "point"),
                    (("inverse",),                  "inverse"),
                ]
                for keys, val in polygon_abbrevs_key_vals:
                    polygon_abbrevs.update(dict.fromkeys(keys, val))
                
                last_hole_type = False
                
                ### ### Membrane mono/bilayer names
                ### Mono as explicitly upper added by request
                monolayer_upper_designation = ["u", "up", "upper", "m", "mo", "mono", "monolayer", "mono_upper"]
                monolayer_lower_designation = ["d", "do", "down", "l", "lo", "lower", "mono_lower", "mono_down"]
                bilayer_designation = ["b", "bi", "bilayer", "memb", "membrane", "both"]
                
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
                            dict_target      = "membrane"
                        elif sub_cmd[1].lower() in monolayer_upper_designation:
                            layer_definition = "upper"
                            dict_target      = "upper_leaf"
                        elif sub_cmd[1].lower() in monolayer_lower_designation:
                            layer_definition = "lower"
                            dict_target      = "lower_leaf"
                    
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
                    
                    ### Creates holes or cutouts in membrane
                    elif sub_cmd[0].lower() in ["h", "hole"] and sub_cmd[1].lower() in hole_types or sub_cmd[0].lower() in hole_types:
                        if sub_cmd[0].lower() in ["h", "hole"]:
                            sub_cmd = sub_cmd[1:]
                        sub_cmd[0] = sub_cmd[0].lower()
                        
                        if "holes" not in settings_dict[dict_target]:
                            settings_dict[dict_target]["holes"] = {
                                "ellipses": [],
                                "polygons": [],
                            }
                        
                        key_vals = []
                        if sub_cmd[-1] in ["inverse", "True", "true"]:
                            key_vals.append(("inverse", True))
                        
                        if sub_cmd[0] in ["circle", "ellipse"]:
                            settings_dict[dict_target]["holes"]["ellipses"].append({
                                "cx": 0,
                                "cy": 0,
                                "xradius": 1,
                                "yradius": 1,
                                "rotation": 0,
                                "buffer": 0,
                                "inverse": False,
                            })

                            ### Command styles:
                            ### circle:cx:cy:radius:inverse
                            ### circle:cx:cy:radius
                            ### circle:radius:inverse
                            ### circle:radius
                            if sub_cmd[0] == "circle":
                                if len(sub_cmd) in [5, 4]:
                                    key_vals.append(("cx", self.get_number_from_string(sub_cmd[1])))
                                    key_vals.append(("cy", self.get_number_from_string(sub_cmd[2])))
                                    key_vals.append(("xradius", self.get_number_from_string(sub_cmd[3])))
                                    key_vals.append(("yradius", self.get_number_from_string(sub_cmd[3])))
                                elif len(sub_cmd) in [3, 2]:
                                    key_vals.append(("xradius", self.get_number_from_string(sub_cmd[1])))
                                    key_vals.append(("yradius", self.get_number_from_string(sub_cmd[1])))
                            
                            if sub_cmd[0] == "ellipse":
                                ### Command styles:
                                ### ellipse:cx:cy:xradius:yradius:rotation:inverse
                                ### ellipse:cx:cy:xradius:yradius:rotation
                                ### ellipse:xradius:yradius:rotation:inverse
                                ### ellipse:xradius:yradius:rotation
                                if len(sub_cmd) in [5, 4]:
                                    key_vals.append(("cx", self.get_number_from_string(sub_cmd[1])))
                                    key_vals.append(("cy", self.get_number_from_string(sub_cmd[2])))
                                    key_vals.append(("xradius", self.get_number_from_string(sub_cmd[3])))
                                    key_vals.append(("yradius", self.get_number_from_string(sub_cmd[4])))
                                elif len(sub_cmd) in [3, 2]:
                                    key_vals.append(("xradius", self.get_number_from_string(sub_cmd[1])))
                                    key_vals.append(("yradius", self.get_number_from_string(sub_cmd[2])))
                            
                            for key, val in key_vals:
                                settings_dict[dict_target]["holes"]["ellipses"][-1][key] = val
                            
                            last_hole_type = "ellipse"
                        
                        elif sub_cmd[0] in ["square", "rectangle", "poly", "polygon"]:
                            settings_dict[dict_target]["holes"]["polygons"].append({
                                "points": [],
                                "rotation": 0,
                                "buffer": 0,
                                "inverse": False,
                            })
                            last_hole_type = "polygon"
                            
                            if sub_cmd[0] in ["square", "rectangle"]:
                                ### Command styles:
                                ### square:cx:cy:sidelen:inverse
                                ### square:cx:cy:sidelen
                                ### square:sidelen:inverse
                                ### square:sidelen

                                ### rectangle:cx:cy:xlen:ylen:inverse
                                ### rectangle:cx:cy:xlen:ylen
                                ### rectangle:xlen:ylen:inverse
                                ### rectangle:xlen:ylen
                                if sub_cmd[0] == "square":
                                    if len(sub_cmd) in [5, 4]:
                                        cx = self.get_number_from_string(sub_cmd[1])
                                        cy = self.get_number_from_string(sub_cmd[2])
                                        xlen = self.get_number_from_string(sub_cmd[3])
                                        ylen = self.get_number_from_string(sub_cmd[3])

                                    elif len(sub_cmd) in [3, 2]:
                                        cx, cy = 0, 0
                                        xlen = self.get_number_from_string(sub_cmd[1])
                                        ylen = self.get_number_from_string(sub_cmd[1])

                                elif sub_cmd[0] == "rectangle":
                                    if len(sub_cmd) in [6, 5]:
                                        cx = self.get_number_from_string(sub_cmd[1])
                                        cy = self.get_number_from_string(sub_cmd[2])
                                        xlen = self.get_number_from_string(sub_cmd[3])
                                        ylen = self.get_number_from_string(sub_cmd[4])

                                    elif len(sub_cmd) in [4, 3]:
                                        cx, cy = 0, 0
                                        xlen = self.get_number_from_string(sub_cmd[1])
                                        ylen = self.get_number_from_string(sub_cmd[2])
                            
                                points = [
                                    (cx-xlen/2, cy-ylen/2),
                                    (cx+xlen/2, cy-ylen/2),
                                    (cx+xlen/2, cy+ylen/2),
                                    (cx-xlen/2, cy+ylen/2),
                                ]

                                key_vals.append(("points", points))
                                
                            for key, val in key_vals:
                                settings_dict[dict_target]["holes"]["polygons"][-1][key] = val
                    
                    elif sub_cmd[0].lower() in ["h", "hole"] and sub_cmd[1] in list(ellipse_abbrevs.keys()) + list(polygon_abbrevs.keys()):
                        if last_hole_type == "ellipse":
                            sub_cmd[1] = ellipse_abbrevs[sub_cmd[1].lower()]

                            key_vals = []

                            if sub_cmd[1] in ["rotation", "buffer"]:
                                key_vals.append((sub_cmd[1], self.get_number_from_string(sub_cmd[2])))

                            elif sub_cmd[1] == "inverse":
                                if len(sub_cmd) == 2 or sub_cmd[-1].lower() in ["inverse", "true"]:
                                    key_vals.append(("inverse", True))
                                else:
                                    key_vals.append(("inverse", False))

                            for key, val in key_vals:
                                settings_dict[dict_target]["holes"]["ellipses"][-1][key] = val
                    
                        elif last_hole_type == "polygon":# or sub_cmd[0].lower() in ["s", "square", "r", "rectangle", "p", "poly", "polygon"]:
                            sub_cmd[1] = polygon_abbrevs[sub_cmd[1].lower()]

                            if sub_cmd[1] in ["rotation", "buffer"]:
                                key_vals.append((sub_cmd[1], self.get_number_from_string(sub_cmd[2])))

                            elif sub_cmd[1] == "inverse":
                                if len(sub_cmd) == 2 or sub_cmd[-1].lower() in ["inverse", "true"]:
                                    key_vals.append((sub_cmd[1], True))
                                else:
                                    key_vals.append((sub_cmd[1], False))

                            ### Points are a bit special as they must be appended to a list instead of overwriting a key:val pair
                            elif sub_cmd[1] == "point":
                                xval = self.get_number_from_string(sub_cmd[2])
                                yval = self.get_number_from_string(sub_cmd[3])
                                settings_dict[dict_target]["holes"]["polygons"][-1]["points"].append((xval, yval))

                            for key, val in key_vals:
                                settings_dict[dict_target]["holes"]["polygons"][-1][key] = val
                    
                    ### For use with "poly/polygon"
                    ### Second way to designate points for polygons
                    ### Essentially an abbreviation of "hole:point" / "h:p"
                    elif sub_cmd[0].lower() in ["point", "poly_point", "polygon_point"]:
                        ### Command style:
                        ### point:xval:yval
                        xval = self.get_number_from_string(sub_cmd[1])
                        yval = self.get_number_from_string(sub_cmd[2])
                        settings_dict[dict_target]["holes"]["polygons"][-1]["points"].append((xval, yval))
                        
                    ### Center [nm]
                    ### "c" removed from abbreviations due to it being an abbreviation for "circle"
#                     elif any(sub_cmd[0].lower() == cen for cen in ["c", "cen", "center"]):
                    elif any(sub_cmd[0].lower() == cen for cen in ["cen", "center"]):
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
                        valid_lipid_optim = sub_cmd[1] in ["avg_optimal", "abs_val", "force_fill", "fill", "no", "insane"]
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
                    
                    ### Readjusts the max/min x/y-values if cutouts/holes are made
                    elif sub_cmd[0].lower() == "readjust_bbox":
                        settings_dict["membrane"]["readjust_bbox"] = ast.literal_eval(sub_cmd[1])
                    
                    ### Splits subleaflet bbox into multiple subleaflets if they are not all connected
                    elif sub_cmd[0].lower() == "split_bbox":
                        settings_dict["membrane"]["split_bbox"] = ast.literal_eval(sub_cmd[1])
                    
                    ### Splits subleaflet bbox into multiple subleaflets if they are not all connected
                    elif sub_cmd[0].lower() == "solvate_empty":
                        settings_dict["membrane"]["solvate_empty"] = ast.literal_eval(sub_cmd[1])
                    
                    elif sub_cmd[0].lower() == "lipid":
                        if any([sub_cmd[1] in self.lipid_dict[params].keys() for params in self.lipid_dict.keys()]):
                            if "lipids_preprocessing" not in settings_dict[dict_target].keys():
                                settings_dict[dict_target]["lipids_preprocessing"] = []
                            settings_dict[dict_target]["lipids_preprocessing"].append(sub_cmd[1])
                    
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
                        "HG_sign":      +1,
                        "leaflet_type": "upper",
                    })
                
                if layer_definition == "lower" or layer_definition == "bilayer":
                    memb_dict["leaflets"]["lower_leaf"] = settings_dict["lower_leaf"]
                    memb_dict["leaflets"]["lower_leaf"].update({
                        "HG_direction": "down",
                        "HG_sign":      -1,
                        "leaflet_type": "lower",
                    })
                
                ### Adding membrane-wide settings to specific leaflets if they are not already given for the specific leaflet
                for key, vals in settings_dict["membrane"].items():
                    for leaflet in memb_dict["leaflets"].values():
                        if key not in leaflet:
                            try:
                                leaflet[key] = copy.deepcopy(vals)
                            except:
                                leaflet[key] = vals
                
                ### Adding default settings for leaflets if none were given
                for key, vals in settings_dict["default"].items():
                    for leaflet in memb_dict["leaflets"].values():
                        if key not in leaflet:
                            try:
                                leaflet[key] = copy.deepcopy(vals)
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
                    for ellipse in leaflet["holes"]["ellipses"]:
                        ellipse["cx"] *= 10
                        ellipse["cy"] *= 10
                        ellipse["xradius"] *= 10
                        ellipse["yradius"] *= 10
                        ellipse["buffer"] *= 10
                    for polygon in leaflet["holes"]["polygons"]:
                        polygon["buffer"] *= 10
                        for pi, (xval, yval) in enumerate(polygon["points"]):
                            polygon["points"][pi] = (xval*10, yval*10)
                    
                ################################
                ### Lipid data incorporation ###
                ################################
                
                ### Reconfigures lipid-specific data for leaflet according to subcommands
                ### Adding zero to "memb_directional_heights" to ensure "min" or "max" is defaults to zero
                memb_directional_heights = [0]
                memb_mean_directional_heights = [0]
                for leaflet_key, leaflet in memb_dict["leaflets"].items():
                    tot_ratio = 0
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

                        leaflet["lipids"][l_name] = copy.deepcopy(self.lipid_dict[l_params][l_name])

                        if type(l_ratio) == str:
                            l_ratio = ast.literal_eval(l_ratio)
                        leaflet["lipids"][l_name].ratio_add(l_ratio)
                        tot_ratio += leaflet["lipids"][l_name].ratio
                            
                        ### Finds charge data from topology files
                        if len(self.ITP_INPUT_cmds) > 0:
                            if leaflet["charge"] == "top" and l_name in self.itp_moltypes.keys():
                                bead_charges = list(map(self.get_number_from_string, [
                                    atom["charge"] for atom in self.itp_moltypes[l_name]["atoms"].values()
                                ]))
                                leaflet["lipids"][l_name].set_bead_charges(bead_charges)
                            else:
                                self.print_term("Lipid ({lipid_name}) could not be found in the topology".format(lipid_name=l_name), warn=True)
                        
                    maxx, maxy, maxz, minx, miny, minz = [
                        func([
                            val + (leaflet["kick" + ax] + leaflet[wr + "_buffer"]) * sign
                            for lipid in leaflet["lipids"].keys()
                            for val in leaflet["lipids"][lipid].get_coords(ax)[0]
                        ])
                        for func, sign in [(max, +1), (min, -1)] for ax, wr in [("x", "plane"), ("y", "plane"), ("z", "height")]
                    ]
                    
                    lipid_radii = [
                            lipid_class.get_radius(AXs="xy")
                            for lipid_class in leaflet["lipids"].values()
                    ]
                    
                    lipid_heights = [
                            lipid_class.get_radius(AXs="z") * 2
                            for lipid_class in leaflet["lipids"].values()
                    ]
                    
                    memb_directional_heights.append(max(lipid_heights) * leaflet["HG_sign"])
                    
                    lipid_heights_ratios = [
                            lipid_class.get_radius(AXs="z") * 2 * lipid_class.ratio / tot_ratio
                            for lipid_class in leaflet["lipids"].values()
                    ]
                    
                    memb_mean_directional_heights.append(np.sum(lipid_heights_ratios) * leaflet["HG_sign"])
                    
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
                        "lipid_height": max(lipid_heights),
                        "lipid_radius": max(lipid_radii),
                    }
                    
                memb_dict.update({
#                     "maxz"           : max([leaflet["lipid_dimensions"]["maxz"]         for leaflet in memb_dict["leaflets"].values()]),
#                     "minz"           : min([leaflet["lipid_dimensions"]["minz"]         for leaflet in memb_dict["leaflets"].values()]),
                    "bead_maxz"      : max(memb_directional_heights),
                    "bead_minz"      : min(memb_directional_heights),
                    "bead_mean_maxz" : max(memb_mean_directional_heights),
                    "bead_mean_minz" : min(memb_mean_directional_heights),
                    "membrane_height": sum([leaflet["lipid_dimensions"]["lipid_height"] for leaflet in memb_dict["leaflets"].values()]),
                })
                
                if return_self == "self":
                    self.MEMBRANES[cmd_nr] = memb_dict
#                     ### Upper leaflet monolayer
#                     if layer_definition in monolayer_upper_designation:
#                         self.MEMBRANES[cmd_nr] = memb_dict
#                     ### Lower leaflet monolayer
#                     elif layer_definition in monolayer_lower_designation:
#                         self.MEMBRANES[cmd_nr] = memb_dict
#                     ### Bilayer
#                     elif layer_definition in bilayer_designation:
#                         self.MEMBRANES[cmd_nr] = memb_dict
                elif return_self == "return":
                    MEMBRANES[cmd_nr] = memb_dict
        
            if return_self == "return":
                return MEMBRANES[cmd_nr]
            
            self.print_term("Number of membranes preprocessed:", len(self.MEMBRANES), spaces=1, verbose=3)

    def solv_preprocessor(self):
        '''
        Preprocesses solvation commands for later ease of use
        '''
        if len(self.SOLVATIONS_cmds) != 0:
            self.print_term("Preprocessing solvent requests", verbose=2)
            for cmd_nr, solvate_cmd in enumerate(self.SOLVATIONS_cmds, 1):
                self.print_term("Starting Solvent command:", cmd_nr, spaces=1, verbose=3)
                ### Defaults
                solv_dict = {
                    ### ### Solvent
                    "solvent": {}, # empty
                    "solv_preprocessing": [], # empty
                    "solv_molarity": 55.56, # [int] or [float] # mol/L
                    "solvfreevol": True,

                    ### ### Ions
                    "neg_ions": {}, # empty
                    "pos_ions": {}, # empty
                    "neg_ions_preprocessing": [], # empty
                    "pos_ions_preprocessing": [], # empty
                    "charge": 0, # int or float # Target charge for system
                    "sys_charge": True, # [bool]
                    ### starting mol/liter concentration of both negative and positive ions (each will have the concentration) 
                    "salt_molarity": 0.15, # int or float # mol/L
                    "ionsvol": "solv", # "box", "free", "solv"
                    "salt_method": "add",
                    
                    ### ### Flooding specific
                    "flooding": False,
                    
                    ### ### General
                    ### Solvent box center:
                    "x": self.pbcx / 10, # [nm] converted to [Å]
                    "y": self.pbcy / 10, # [nm] converted to [Å]
                    "z": self.pbcz / 10, # [nm] converted to [Å]
                    "center": [0, 0, 0], # [nm] converted to [Å]
                    
                    "charges_from_top": True, # bool
                    
                    ### "count" and "solv_per_lipid" are mutually exclusive
                    ### "solv_per_lipid" takes priority if both are given
                    "solv_per_lipid": False, # Number of solv particles per lipid contained in solvent box
                    "count": False, # [bool] Uses specific molarity value as absolute number of molecules instead of molarity ratio. Will be rounded using int(val+0.5)
#                     "kick": 0.264/5*2, # [nm] converted to [Å]
                    "kick": 0.264/4, # [nm] converted to [Å]
                    "mapping": True,# bool, Whether AA-to-CG mapping should be considered
                    
                    "bdx": 1.0, # [multiplier]
                    "bdy": 1.0, # [multiplier]
                    "bdz": 1.0, # [multiplier]
                    
                    "params": False, # False or str
                    "bead_radius": 0.264, # [nm] converted to [Å] # Used for volume calculations
                    "gridres": 0.264, # [nm] converted to [Å] # 1.32
                    
                    "WR": 0.264, # [nm] converted to [Å]
                    "buffer": 0.2, # [nm] converted to [Å]
                    "protein_extra_buffer": 0.2, # [nm] converted to [Å]
                    "lipid_extra_buffer": 0, # [nm] converted to [Å]
                    "solute_extra_buffer": 0, # [nm] converted to [Å]
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

                    ### Target charge [int/float] or False for no "neutralization"
                    elif sub_cmd[0].lower() == "charge":
                        solv_dict["charge"] = ast.literal_eval(sub_cmd[1])

                    ### Whether to consider system charge for charge calculations
                    elif sub_cmd[0].lower() == "sys_charge":
                        if type(sub_cmd[1]) == bool:
                            solv_dict["sys_charge"] = sub_cmd[1]
                        else:
                            solv_dict["sys_charge"] = ast.literal_eval(sub_cmd[1])

                    ### Salt concentration
                    elif sub_cmd[0].lower() in ["salt_molarity", "salt"]:
                        solv_dict["salt_molarity"] = ast.literal_eval(sub_cmd[1])

                    ### Whether to consider free volume (excluding lipid and protein) or purely the box volume
                    elif sub_cmd[0].lower() == "ionsvol":
                        solv_dict["ionsvol"] = sub_cmd[1].lower()

                    ### Whether to consider system charge for charge calculations
                    elif sub_cmd[0].lower() == "salt_method":
                        valid_settings = ["add", "remove", "mean"]
                        if sub_cmd[1] in valid_settings:
                            solv_dict["salt_method"] = sub_cmd[1]
                        else:
                            self.print_term(
                                "Subcommand", "'salt_method'", "has been given with an invalid setting",
                                "\nSetting used was:", "'" + sub_cmd[1] + "'",
                                "\nValid settings are:", " ".join(["'"+s+"'" for s in valid_settings]),
                                warn=True
                            )

                    ###############
                    ### GENERAL ###
                    ###############
                    ### Solvation box x-length
                    elif sub_cmd[0].lower() == "x":
                        solv_dict["x"] = ast.literal_eval(sub_cmd[1])
                    ### Solvation box y-length
                    elif sub_cmd[0].lower() == "y":
                        solv_dict["y"] = ast.literal_eval(sub_cmd[1])
                    ### Solvation box z-length
                    elif sub_cmd[0].lower() == "z":
                        solv_dict["z"] = ast.literal_eval(sub_cmd[1])
                        
                    ### Solvation box x center
                    elif sub_cmd[0].lower() == "cx":
                        solv_dict["center"][0] = ast.literal_eval(sub_cmd[1])
                    ### Solvation box y center
                    elif sub_cmd[0].lower() == "cy":
                        solv_dict["center"][1] = ast.literal_eval(sub_cmd[1])
                    ### Solvation box z center
                    elif sub_cmd[0].lower() == "cz":
                        solv_dict["center"][2] = ast.literal_eval(sub_cmd[1])
                        
                    ### Solvation box center
                    elif sub_cmd[0].lower() == "center":
                        solv_dict["center"] = [ast.literal_eval(val) for val in sub_cmd[1:]]
                    
                    ### Whether charges should be obtained from topology
                    elif sub_cmd[0].lower() in ["charges_from_top", "top"]:
                        val = sub_cmd[1]
                        if type(val) == str:
                            val = ast.literal_eval(sub_cmd[1])
                        solv_dict["charges_from_top"] = val
                        
                    ### Whether to use ratios as absolute number of molecules. True/False
                    elif sub_cmd[0].lower() == "solv_per_lipid":
                        solv_dict["solv_per_lipid"] = self.get_number_from_string(sub_cmd[1])
                        
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
                    
                    ### Whether charges should be obtained from topology
                    elif sub_cmd[0].lower() == "mapping":
                        val = sub_cmd[1]
                        if type(val) == str:
                            val = ast.literal_eval(sub_cmd[1])
                        solv_dict["mapping"] = val
                    
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
                    
                    #########################
                    ### FLOODING SPECIFIC ###
                    #########################
                    ### Used to check if command is made as flooding or not
                    elif sub_cmd[0].lower() == "flooding":
                        solv_dict["flooding"] = ast.literal_eval(sub_cmd[1])
                        
                    elif solv_dict["flooding"] and sub_cmd[0] in [key for params in self.solvent_defs.keys() for key in self.solvent_defs[params].keys()]:
                        solv_dict["solv_preprocessing"].append(sub_cmd[0:])
                    
                    ### Errors out if unknown subcommand used, and prints the subcommand to terminal
                    else:
                        assert False, "Unknown subcommand given to '-solvate'. The subcommand is: '" + str(cmd) + "'"
                
                ### fixing a couple of values such that they are in angstrom
                solv_dict.update({
                    "x": solv_dict["x"]*10,
                    "y": solv_dict["y"]*10,
                    "z": solv_dict["z"]*10,
                    "center": [val*10 for val in solv_dict["center"]],
                    
                    "kick":        solv_dict["kick"]*10,
                    "bead_radius": solv_dict["bead_radius"]*10,
                    "gridres":     solv_dict["gridres"]*10,
                    "WR":          solv_dict["WR"]*10,
                    
                    "buffer":               solv_dict["buffer"]*10,
                    "protein_extra_buffer": solv_dict["protein_extra_buffer"]*10,
                    "lipid_extra_buffer":   solv_dict["lipid_extra_buffer"]*10,
                    "solute_extra_buffer":  solv_dict["solute_extra_buffer"]*10,
                })
                
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
                            isnumber, isinteger = self.is_number(val)
                            if isnumber:
                                rest_molarity = ast.literal_eval(val)
                            else:
                                rest_params = val
                    return rest_params, rest_molarity
                
                solv_dict["solv_tot_ratio"] = 0
                solv_dict["pos_tot_ratio"] = 0
                solv_dict["neg_tot_ratio"] = 0
                
                ### ### SOLVENT DATA
                for subcommand_values in solv_dict["solv_preprocessing"]:
                    name = subcommand_values[0]
                    rest = subcommand_values[1:]
                    solv_type = "solvent"
                    rest_params, rest_ratio = values_checker(rest)
                    if rest_params:
                        params = rest_params
                    else:
                        params = solv_dict["params"] or self.solv_params or self.sys_params # (sys_params defaults to "default")
                    ### Adding name:dict combo to solvent dict
                    solv_dict[solv_type][name] = copy.deepcopy(self.solvent_dict[params][name])
                    
                    ### Finds charge data from topology files
                    if len(self.ITP_INPUT_cmds) > 0:
                        if solv_dict["charges_from_top"] and name in self.itp_moltypes.keys():
                            bead_charges = list(map(self.get_number_from_string, [
                                atom["charge"] for atom in self.itp_moltypes[name]["atoms"].values()
                            ]))
                            solv_dict[solv_type][name].set_bead_charges(bead_charges)
                        else:
                            self.print_term("Solvent ({name}) could not be found in the topology".format(name=name), warn=True)
                    
                    ### If solv_dict["solv_per_lipid"] has been set, then ignore count command
                    if solv_dict["solv_per_lipid"] and solv_dict["count"]:
                        solv_dict["count"] = False
                    
                    ### If count has been set, then convert to integer value and treat as absolute number of molecules
                    if solv_dict["count"]:
                        if rest_ratio:
                            solv_dict[solv_type][name].molarity_set(int(rest_ratio + 0.5))
                            solv_dict[solv_type][name].ratio_set(0)
                        else:
                            solv_dict[solv_type][name].molarity_set(int(solv_dict["solv_molarity"] + 0.5))
                            solv_dict[solv_type][name].ratio_set(0)
                    else:
                        if rest_ratio:
                            solv_dict[solv_type][name].molarity_set(solv_dict["solv_molarity"])
                            solv_dict[solv_type][name].ratio_set(rest_ratio)
                            solv_dict["solv_tot_ratio"] += rest_ratio
                        else:
                            solv_dict[solv_type][name].molarity_set(solv_dict["solv_molarity"])
                            solv_dict[solv_type][name].ratio_set(1)
                            solv_dict["solv_tot_ratio"] += 1
                    
                    ### If mapping ratio should be ignored, set to 1 for all solvents
                    if not solv_dict["mapping"]:
                        solv_dict[solv_type][name].mapping_ratio_set(1)
                
                ion_types_to_be_processed = []
                if solv_dict["pos_ions_preprocessing"]:
                    ion_types_to_be_processed.extend([("pos_ions", ions) for ions in solv_dict["pos_ions_preprocessing"]])
                if solv_dict["neg_ions_preprocessing"]:
                    ion_types_to_be_processed.extend([("neg_ions", ions) for ions in solv_dict["neg_ions_preprocessing"]])
                    
                ### ### Ion processing
                for solv_type, subcommand_values in ion_types_to_be_processed:
                    name = subcommand_values[0]
                    rest = subcommand_values[1:]
                    rest_params, rest_ratio = values_checker(rest)
                    if rest_params:
                        params = rest_params
                    else:
                        params = solv_dict["params"] or self.solv_params or self.sys_params # (sys_params defaults to "default")
                    
                    ### Adding name:dict combo to solvent dict
                    if solv_type == "pos_ions":
                        solv_dict[solv_type][name] = copy.deepcopy(self.ion_dict[params]["positive"][name])
                        solv_ratio_type = "pos_tot_ratio"
                    elif solv_type == "neg_ions":
                        solv_dict[solv_type][name] = copy.deepcopy(self.ion_dict[params]["negative"][name])
                        solv_ratio_type = "neg_tot_ratio"
                    
                    ### Finds charge data from topology files
                    if len(self.ITP_INPUT_cmds) > 0:
                        if solv_dict["charges_from_top"] and name in self.itp_moltypes.keys():
                            bead_charges = list(map(self.get_number_from_string, [
                                atom["charge"] for atom in self.itp_moltypes[name]["atoms"].values()
                            ]))
                            solv_dict[solv_type][name].set_bead_charges(bead_charges)
                        else:
                            self.print_term("Solvent ({name}) could not be found in the topology".format(name=name), warn=True)
                    
                    ### If count has been set, then convert to integer value and treat as absolute number of molecules
                    if solv_dict["count"]:
                        if rest_ratio:
                            solv_dict[solv_type][name].molarity_set(int(rest_ratio + 0.5))
                            solv_dict[solv_type][name].ratio_set(0)
                        else:
                            solv_dict[solv_type][name].molarity_set(int(solv_dict["salt_molarity"] + 0.5))
                            solv_dict[solv_type][name].ratio_set(0)
                    else:
                        if rest_ratio:
                            solv_dict[solv_type][name].molarity_set(solv_dict["salt_molarity"])
                            solv_dict[solv_type][name].ratio_set(rest_ratio)
                            solv_dict[solv_ratio_type] += rest_ratio
                        else:
                            solv_dict[solv_type][name].molarity_set(solv_dict["salt_molarity"])
                            solv_dict[solv_type][name].ratio_set(1)
                            solv_dict[solv_ratio_type] += 1
                    
                    ### If mapping ratio should be ignored, set to 1 for all solvents
                    if not solv_dict["mapping"]:
                        solv_dict[solv_type][name].mapping_ratio_set(1)

                self.SOLVATIONS[cmd_nr] = solv_dict.copy()
                
            self.print_term("Number of solvent commands preprocessed:", len(self.SOLVATIONS), spaces=1, verbose=3)
    
    def itp_read_initiater(self):
        if len(self.ITP_INPUT_cmds) != 0:
            self.print_term("Loading topology file(s)", spaces=0, verbose=2)
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

            self.print_term("Finished loading topologies. Number of moleculetypes found:", len(self.itp_moltypes), "\n", spaces=1, verbose=3)
    
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
        ##################### CURRENTLY NO LONGER USED

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
        
        if type(s) == str:
            if s.startswith("-"):
                s = s[1:]
        else:
            s = str(s)
        
        try:
            test = float(s)
            number = True
            if s.isdigit():
                integer = True
            else:
                integer = False
            return number, integer
        
        except ValueError:
            return False, False
    
    def get_number_from_string(self, s):
        isnumber, isint = self.is_number(s)
        if isnumber and isint:
            return int(s)
        elif isnumber:
            return float(s)
        else:
            return False
    
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
            self.print_term("------------------------------ PROTEIN PLACEMENT", spaces=0, verbose=1)
            for protein_i, (protein_nr, protein) in enumerate(self.PROTEINS.items()):
                if protein_i != 0:
                    self.print_term("", verbose=3)
                self.print_term("Starting protein nr", protein_nr, spaces=0, verbose=2)
                
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
                self.print_term(
                    "Centering protein using", "'" + " ".join([str(i) for i in protein["cen_method"]])+"'",
                    "at x/y/z:", round(xcen, 3), round(ycen, 3), round(zcen, 3), "(Input file coordinate system [Å])",
                    spaces=1,
                    verbose=3
                )
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
                self.print_term("New protein center at x/y/z:", round(xcen_new, 3), round(ycen_new, 3), round(zcen_new, 3), "(Internal coordinate system [Å])", spaces=1, verbose=3)
                
                self.print_term("Finished placing protein nr", protein_nr, spaces=1, verbose=3)
                
            prot_placer_toc = time.time()
            prot_placer_time = round(prot_placer_toc - prot_placer_tic, 4)
            self.print_term("------------------------------ PLACEMENT COMPLETE", "(time spent: "+str(prot_placer_time)+" [s])", "\n", spaces=0, verbose=1)
    
    ######################
    ### POLYGON MAKERS ###
    ######################
    def subleaflet_poly_maker(self):
        if len(self.MEMBRANES) != 0:
            subleaflet_poly_maker_tic = time.time()
            self.print_term("------------------------------ CREATING LEAFLET BOUNDARY BOXES", spaces=0, verbose=1)
            for memb_i, (memb_key, memb_dict) in enumerate(self.MEMBRANES.items()):
                if memb_i != 0:
                    self.print_term("", verbose=2)
                self.print_term("Starting membrane nr", memb_key, spaces=0, verbose=2)
                
                for leaflet_i, (leaflet_key, leaflet) in enumerate(memb_dict["leaflets"].items()):
                    leaflet["subleaflets"] = {}
                    self.print_term("Starting", leaflet["leaflet_type"], "leaflet", spaces=1, verbose=3)
                    self.print_term("Base parameters:", "x=" + str(leaflet["x"] / 10) + "nm", "y=" + str(leaflet["y"] / 10) + "nm", "APL=" + str(leaflet["apl"] / 100) + "nm^2", spaces=2, verbose=4)
                    
                    if leaflet["gridsplits"][0] == "auto":
                        xsplits = leaflet["x"] // leaflet["gridsplits"][1] + 1
                        ysplits = leaflet["y"] // leaflet["gridsplits"][1] + 1
                    else:
                        xsplits, ysplits = leaflet["gridsplits"]
                    if xsplits > 1:
                        self.print_term("x-axis split into", xsplits, "subleaflets due to axis length or manual designation", spaces=2, verbose=4)
                    if ysplits > 1:
                        self.print_term("y-axis split into", ysplits, "subleaflets due to axis length or manual designation", spaces=2, verbose=4)
                    if xsplits * ysplits > 1:
                        self.print_term("A total of", xsplits*ysplits, "subleaflets have been made for the current leaflet", spaces=2, verbose=4)
                    
                    xmin, xmax = leaflet["center"][0] - leaflet["x"]/2, leaflet["center"][0] + leaflet["x"]/2
                    ymin, ymax = leaflet["center"][1] - leaflet["y"]/2, leaflet["center"][1] + leaflet["y"]/2
                    
                    xpoint_vals = np.linspace(start=xmin, stop=xmax, num=round(xsplits+1), endpoint=True)
                    ypoint_vals = np.linspace(start=ymin, stop=ymax, num=round(ysplits+1), endpoint=True)
                    
                    sli = 0
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
                            
                            leaflet["subleaflets"][(xi, yi, sli)] = {
                                "xmin":          xmin, # Cutouts change these values
                                "xmax":          xmax, # Cutouts change these values
                                "ymin":          ymin, # Cutouts change these values
                                "ymax":          ymax, # Cutouts change these values
                                "xmin_original": xmin, # Cutouts DO NOT change these values
                                "xmax_original": xmax, # Cutouts DO NOT change these values
                                "ymin_original": ymin, # Cutouts DO NOT change these values
                                "ymax_original": ymax, # Cutouts DO NOT change these values
                                "box_points":    box_points,
                                "box_poly":      box_Poly,
                            }
                            sli += 1
            
            subleaflet_poly_maker_toc = time.time()
            subleaflet_poly_maker_time = round(subleaflet_poly_maker_toc - subleaflet_poly_maker_tic, 4)
            self.print_term("------------------------------ LEAFLET BOUNDARY BOXES CREATED", "(time spent: "+str(subleaflet_poly_maker_time)+" [s])", "\n", spaces=0, verbose=1)
            
    def holed_subleaflet_bbox_maker(self):
        if len(self.MEMBRANES) != 0:
            holed_subleaflet_bbox_maker_tic = time.time()
            self.print_term("------------------------------ CREATING HOLED BOUNDARY BOXES", spaces=0, verbose=1)
            for memb_i, (memb_key, memb_dict) in enumerate(self.MEMBRANES.items()):
                if memb_i != 0:
                    self.print_term("", verbose=3)
                self.print_term("Starting membrane nr", memb_key, spaces=0, verbose=2)
                
                ### Defining some default values for union, intersection and holed_bbox
                for leaflet_i, (leaflet_key, leaflet) in enumerate(memb_dict["leaflets"].items()):
                    leaflet["protein_poly"]           = False
                    leaflet["remove_union"]           = False
                    leaflet["require_union"]          = False
                    leaflet["prot_points"]            = False
                    leaflet["prot_Points_buffered"]   = False
                    leaflet["alphashape_1"]           = False
                    leaflet["ConcaveHulls_Polygon_1"] = False
                    for (slxi, slyi, sli), subleaflet in leaflet["subleaflets"].items():
                        subleaflet["intersection"] = False
                        subleaflet["holed_bbox"]   = subleaflet["box_poly"]
                
                ######################################
                ### HANDLING PROTEIN RELATED HOLES ###
                ######################################
                ### Getting all the protein bead positions
                if len(self.PROTEINS) != 0:
                    self.print_term("Finding protein beads inside leaflets", spaces=1, verbose=3)

                    for leaflet_i, (leaflet_key, leaflet) in enumerate(memb_dict["leaflets"].items()):
                        self.print_term("Starting", leaflet["leaflet_type"], "leaflet", spaces=2, verbose=4)
                        polygon_list = []

                        ####################################################
                        ### Finding protein points contained in leaflets ###
                        ####################################################
                        prot_beads_in_memb = [
                            (beadx, beady)
                            for protein_i, (protein_nr, protein) in enumerate(self.PROTEINS.items())
                            for beadx, beady, beadz in protein["protein"].get_beads("xyz")
                            if (
                                leaflet["HG_direction"] == "up" and (
                                    leaflet["center"][2] - protein["buffer"] < beadz < leaflet["center"][2] + leaflet["lipid_dimensions"]["lipid_height"] + protein["buffer"]
                                )
                            )
                            or (
                                leaflet["HG_direction"] == "down" and (
                                    leaflet["center"][2] + protein["buffer"] > beadz > leaflet["center"][2] - leaflet["lipid_dimensions"]["lipid_height"] - protein["buffer"]
                                )
                            )
                        ]

                        if len(prot_beads_in_memb) == 0:
                            self.print_term("No leaflet-protein overlap found", spaces=3, verbose=5)
                        else:
                            self.print_term("Leaflet-protein overlap found for", str(len(prot_beads_in_memb)), "protein beads", spaces=3, verbose=5)
                            leaflet["prot_points"] = prot_beads_in_memb
                            #################################
                            ### CONCAVE HULL / ALPHASHAPE ###
                            #################################
                            Points = [shapely.Point(point) for point in prot_beads_in_memb]
                            
                            Points_buffered = []
                            for Point in Points:
                                Points_buffered.append(
                                    Point.buffer(
                                        distance = leaflet["prot_buffer"], # [Å]
                                        resolution = 3, # "quad_segs" in the docs
                                    )
                                )
                            
                            leaflet["prot_Points_buffered"] = Points_buffered
                            
                            Points_buffered_union = shapely.ops.unary_union(Points_buffered)
                            
                            if Points_buffered_union.geom_type == "Polygon":
                                self.print_term("Polygon", debug = True)
                                Points_buffered_union_polygons = [Points_buffered_union]
                            elif Points_buffered_union.geom_type == "MultiPolygon":
                                self.print_term("MultiPolygon", debug = True)
                                Points_buffered_union_polygons = list(Points_buffered_union.geoms)
                            elif Points_buffered_union[-1].geom_type == "MultiPolygon":
                                self.print_term("list of MultiPolygon", debug = True)
                                Points_buffered_union_polygons = list(Points_buffered_union[-1].geoms)
                            
                            poly_points = []
                            for poly in list(Points_buffered_union_polygons):
                                xs, ys = poly.exterior.coords.xy
                                poly_points.extend(list(zip(xs, ys)))
                            
                            leaflet["prot_poly_points"] = poly_points
    
                            alpha = 1 / (leaflet["lipid_dimensions"]["lipid_radius"] * 2 * leaflet["alpha_mult"] + leaflet["prot_buffer"])
                            ### alpha = radius^-1
                            ALPHASHAPE = alphashape.alphashape(points = poly_points, alpha = alpha)
                            self.print_term(type(ALPHASHAPE), debug = True)

                            leaflet["alphashape_1"] = ALPHASHAPE
                            
                            ### Alphashape output can be a bit unpredictable so need to check all possibilities
                            if ALPHASHAPE.geom_type == "Polygon":
                                self.print_term("Polygon", debug = True)
                                ConcaveHulls_Polygon = [ALPHASHAPE]
                            elif ALPHASHAPE.geom_type == "MultiPolygon":
                                self.print_term("MultiPolygon", debug = True)
                                ConcaveHulls_Polygon = list(ALPHASHAPE.geoms)
                            elif ALPHASHAPE[-1].geom_type == "MultiPolygon":
                                self.print_term("list of MultiPolygon", debug = True)
                                ConcaveHulls_Polygon = list(ALPHASHAPE[-1].geoms)
                            
                            leaflet["ConcaveHulls_Polygon_1"] = ConcaveHulls_Polygon
                
                ####################################################
                ### COMBINING SHAPELY SHAPES FOR ALL SUBLEAFLETS ###
                ####################################################
                self.print_term("", verbose=4) # verbose=4 to avoid blank line for verbose < 4
                self.print_term("Calculating holed boundary box", spaces=1, verbose=3)
                for leaflet_i, (leaflet_key, leaflet) in enumerate(memb_dict["leaflets"].items()):
                    self.print_term("Starting", leaflet["leaflet_type"], "leaflet", spaces=2, verbose=4)
                    
                    ### All the various polygons to be removed from the bbox
                    ### Unique to each leaflet but the same for each subleaflet
                    ### To be expanded later with manually defined holes
                    poly_remove_list = []
                    poly_require_list = []
                    
                    if "ConcaveHulls_Polygon_1" in leaflet.keys():
                        protein_poly = leaflet["ConcaveHulls_Polygon_1"]
                        for val in [protein_poly]:
                            if val:
                                poly_remove_list.extend(val)
                    
                    #################################
                    ### HANDLING CIRCLES/ELLIPSES ###
                    #################################
                    if leaflet["holes"]["ellipses"]:
                        for settings in leaflet["holes"]["ellipses"]:
                            circle = shapely.geometry.Point((settings["cx"], settings["cy"])).buffer(1)
                            circle = shapely.affinity.scale(circle, settings["xradius"], settings["yradius"])
                            if settings["rotation"]:
                                circle = shapely.affinity.rotate(circle, settings["rotation"])
                            if settings["buffer"]:
                                circle = circle.buffer(settings["buffer"])
                            if settings["inverse"]:
                                poly_require_list.append(circle)
                            else:
                                poly_remove_list.append(circle)
                    
                    #########################
                    ### HANDLING POLYGONS ###
                    #########################
                    if leaflet["holes"]["polygons"]:
                        for settings in leaflet["holes"]["polygons"]:
                            if len(settings["points"]) > 2:
                                ### Polygon must contain at least two points
                                polygon = shapely.Polygon(settings["points"])
                            elif len(settings["points"]) == 2 and settings["buffer"] > 0:
                                ### or be a line with a buffer
                                polygon = shapely.LineString(settings["points"])
                            elif len(settings["points"]) == 1 and settings["buffer"] > 0:
                                ### or be a point with a buffer (e.g. just a circle)
                                polygon = shapely.Point(settings["points"])
                            else:
                                continue
                            if settings["rotation"]:
                                polygon = shapely.affinity.rotate(polygon, settings["rotation"])
                            if settings["buffer"]:
                                polygon = polygon.buffer(settings["buffer"])
                            if settings["inverse"]:
                                poly_require_list.append(polygon)
                            else:
                                poly_remove_list.append(polygon)
                    
                    ### Combining polygons into single shape
                    remove_union = shapely.unary_union(poly_remove_list)
                    leaflet["remove_union"] = remove_union
                    
                    if poly_require_list:
                        require_union = shapely.unary_union(poly_require_list)
                        leaflet["require_union"] = require_union
                    
                    for (slxi, slyi, sli), subleaflet in leaflet["subleaflets"].items():
                        ### The subleaflet box bbox polygon
                        box_poly = subleaflet["box_poly"]
                        if leaflet["require_union"]:
                            box_poly = shapely.intersection(leaflet["require_union"], box_poly)
                        intersection = shapely.intersection(leaflet["remove_union"], box_poly)
                        holed_bbox   = shapely.difference(box_poly, intersection)
                        
                        subleaflet["intersection"] = intersection
                        subleaflet["holed_bbox"]   = holed_bbox
                        
                        ### Readjusts the dimensions of the subleaflet
                        if leaflet["readjust_bbox"]:
                            new_xmin, new_ymin, new_xmax, new_ymax = subleaflet["holed_bbox"].bounds
                            subleaflet["xmin"] = new_xmin
                            subleaflet["ymin"] = new_ymin
                            subleaflet["xmax"] = new_xmax
                            subleaflet["ymax"] = new_ymax
                        
                    ### Splits subleaflet into multiple subleaflets if the individual parts are not touching
                    if leaflet["split_bbox"]:
                        try:
                            ### Following code can easily encounter problems
                            ### Using "try" until all potential errors have been handled 
                            new_subleaflets = {}
                            new_sli = 0
                            for (slxi, slyi, sli), subleaflet in leaflet["subleaflets"].items():
                                if subleaflet["holed_bbox"].geom_type == "Polygon":
                                    self.print_term("Polygon", debug = True)
                                    geoms = [subleaflet["holed_bbox"]]
                                elif subleaflet["holed_bbox"].geom_type == "MultiPolygon":
                                    self.print_term("MultiPolygon", debug = True)
                                    geoms = list(subleaflet["holed_bbox"].geoms)
                                elif subleaflet["holed_bbox"][-1].geom_type == "MultiPolygon":
                                    self.print_term("list of MultiPolygon", debug = True)
                                    geoms = list(subleaflet["holed_bbox"][-1].geoms)

                                for geom in geoms:
                                    new_subleaflet = copy.deepcopy(subleaflet)
                                    new_subleaflet["original_holed_bbox"] = new_subleaflet["holed_bbox"]
                                    new_subleaflet["holed_bbox"] = geom
                                    if leaflet["readjust_bbox"]:
                                        new_xmin, new_ymin, new_xmax, new_ymax = new_subleaflet["holed_bbox"].bounds
                                        new_subleaflet["xmin"] = new_xmin
                                        new_subleaflet["ymin"] = new_ymin
                                        new_subleaflet["xmax"] = new_xmax
                                        new_subleaflet["ymax"] = new_ymax
                                    new_subleaflets[(slxi, slyi, new_sli)] = new_subleaflet
                                    new_sli += 1
                            leaflet["subleaflets"] = new_subleaflets
                        except:
                            ### If error encounted then don't do it
                            self.print_term(
                                "CGSB tried to split the subleaflet bbox but encountered an error.\n",
                                "This bit of code is prone to errors due to many potential inputs.\n",
                                "Will continue running without splitting the bbox.",
                                warn=True
                            )
            
            holed_subleaflet_bbox_maker_toc = time.time()
            holed_subleaflet_bbox_maker_time = round(holed_subleaflet_bbox_maker_toc - holed_subleaflet_bbox_maker_tic, 4)
            self.print_term("------------------------------ HOLED BOUNDARY BOXES CREATED", "(time spent: "+str(holed_subleaflet_bbox_maker_time)+" [s])", "\n", spaces=0, verbose=1)
    
    ########################
    ### LIPID CALCULATOR ###
    ########################
    def lipid_calculator(self):
        self.SYS_lipids_dict = {}
        if len(self.MEMBRANES) != 0:
            lipid_calculator_tic = time.time()
            self.print_term("------------------------------ CALCULATING LIPID RATIOS", spaces=0, verbose=1)
            printer_spacing = 1
            for memb_i, (memb_key, memb_dict) in enumerate(self.MEMBRANES.items()):
                membrane_lipid_count_dict = {}
                if memb_i != 0:
                    self.print_term("", verbose=3)
                self.print_term("Starting membrane nr", memb_key, spaces=0, verbose=2)
                
                ### Checks if any leaflets contain multiple subleaflets for printer spacing
                for leaflet_i, (leaflet_key, leaflet) in enumerate(memb_dict["leaflets"].items()):
                    if len(leaflet["subleaflets"]) > 1:
                        printer_spacing = 2
                
                for leaflet_i, (leaflet_key, leaflet) in enumerate(memb_dict["leaflets"].items()):
                    if leaflet_i != 0:
                        self.print_term("", verbose=4)
                    self.print_term("Starting", leaflet["leaflet_type"], "leaflet", spaces=1, verbose=3)
                    
                    ### Lipid optimization printing
                    ### Printing here so it isn't spammed if multiple subleafs
                    self.print_term("Lipid optimization 'lipid_optim' setting:", leaflet["lipid_optim"], spaces=2, verbose=4)
                    if leaflet["lipid_optim"] in ["force_fill", "fill"]:
                        if leaflet["lipid_optim"] == "fill":
                            self.print_term("Filling leaflet until perfect ratio between lipids is achieved or leaflet is full", spaces=3, verbose=5)
                        elif leaflet["lipid_optim"] == "force_fill":
                            self.print_term("Forcefully filling leaflet until all grid points have a lipid", spaces=3, verbose=5)
                    elif leaflet["lipid_optim"] == "avg_optimal":
                        self.print_term("Optimizing based on the average deviation from expected ratios", spaces=3, verbose=5)
                    
                    apl_sqrt = math.sqrt(leaflet["apl"])
                    
                    leaflet["leaf_lipid_count_dict"] = {}
                    leaflet["lipid_names"] = [lipid_name for lipid_name in leaflet["lipids"].keys()]
                    leaflet["lipid_ratios"] = [0 for _ in leaflet["lipid_names"]]
                    
                    if len(leaflet["subleaflets"]) > 1:
                        printer_leaf_name = "Subleaflet"
                    else:
                        printer_leaf_name = "Leaflet"
                        
                    for subleaflet_i, ((slxi, slyi, sli), subleaflet) in enumerate(leaflet["subleaflets"].items()):
                        if subleaflet_i != 0:
                            self.print_term("", verbose=5)
                        if len(leaflet["subleaflets"]) > 1:
                            self.print_term("Starting subleaflet", str(subleaflet_i+1)+"/"+str(len(leaflet["subleaflets"])), spaces=2, verbose=5)
                        
                        ### Initial area per lipid calculations and max potential lipids in area
                        
                        ### Rounding max number of possible lipids according to requested rounding method
                        if leaflet["lipid_optim"] == "insane":
                            nlipidsx = abs(subleaflet["xmax"] - subleaflet["xmin"]) / math.sqrt(leaflet["apl"])
                            nlipidsy = abs(subleaflet["ymax"] - subleaflet["ymin"]) / math.sqrt(leaflet["apl"])
                            if "intersection" in subleaflet:
                                area_ratio = math.sqrt(subleaflet["holed_bbox"].area / subleaflet["box_poly"].area)
                                ### Using 'round' here instead of 'int' because insane-protein compatability part is very iffy
                                ### Definitely not the same as insane, but too much has to be changed for it to be possible to make systems identical to insane systems
                                nlipidsx = round(nlipidsx * area_ratio)
                                nlipidsy = round(nlipidsy * area_ratio)
                            max_lipids_possible = int(nlipidsx) * int(nlipidsy)
                        else:
                            subleaflet_area = subleaflet["holed_bbox"].area
                            max_lipids_possible_decimal = subleaflet_area / leaflet["apl"]
                            if leaflet["lip_round_func"][1] == int: # int rounding
                                max_lipids_possible = int(max_lipids_possible_decimal + 0.5)
                            else: # round(), math.floor() or math.ceil() rounding
                                max_lipids_possible = leaflet["lip_round_func"][1](max_lipids_possible_decimal)
                        
                        self.print_term("Maximum number of lipids allowed:", max_lipids_possible, spaces=3, verbose=5)
                        self.print_term("", verbose=4)
                        
                        ### Initial rounded estimations
                        lipids = [(lipid_name, lipid_vals.ratio) for lipid_name, lipid_vals in leaflet["lipids"].items()]
                        lipid_names = [name for name, ratio in lipids]
                        lipid_ratios_decimal = [ratio for name, ratio in lipids]
                        lipids_tot = sum(lipid_ratios_decimal)
                        lipid_ratios = [round(int(i / lipids_tot * max_lipids_possible), 3) for i in lipid_ratios_decimal]
                        lipid_ratios_TEST = [(i, lipids_tot, max_lipids_possible) for i in lipid_ratios_decimal]

                        original_ratios = lipid_ratios[:]
                        original_ratios_decimal = [round(i / sum(original_ratios) * 100, 3) for i in original_ratios]

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
                                lipid_ratios = iters[-1]

                            ### Optimizing lipids
                            elif leaflet["lipid_optim"] == "avg_optimal": # lipid_ratios_decimal
                                iters_decimal = [[round(iteration[i] / sum(iteration) * 100, 3) for i in range(len(iteration))] for iteration in iters]
                                iters_abs_diff = [
                                    [abs(iteration[i] - expected_ratios_decimal[i]) for i in range(len(iteration))]
                                    for iteration in iters_decimal
                                ]

                                iters_avg = [(np.mean(iters), i) for i, iters in enumerate(iters_abs_diff)]
                                iters_avg_optimal = min(iters_avg, key=lambda i: i[0])
                                lipid_ratios = iters[iters_avg_optimal[1]]
                        
                        elif leaflet["lipid_optim"] in ["no", "insane"]:
                            ### E.g. do nothing. Just here to show that the options are understood.
                            pass
                        
                        subleaflet["lipid_names"]  = lipid_names
                        subleaflet["lipid_ratios"] = lipid_ratios
                        
                        for name, count in zip(subleaflet["lipid_names"], subleaflet["lipid_ratios"]):
                            if name not in leaflet["leaf_lipid_count_dict"]:
                                leaflet["leaf_lipid_count_dict"][name] = count
                            else:
                                leaflet["leaf_lipid_count_dict"][name] += count
                        
                        if self.extra_info:
                            if (len(leaflet["subleaflets"]) > 1 and self.verbose >= 5) or (len(leaflet["subleaflets"]) == 1 and self.verbose >= 4):
                                ### Printing mean, min and max deviation from wanted ratios
                                ratios_final_decimal = [round(lipid_ratios[i] / sum(lipid_ratios) * 100, 3) for i in range(len(lipid_ratios))]
                                ratios_final_abs_diff = [abs(ratios_final_decimal[i] - expected_ratios_decimal[i]) for i in range(len(ratios_final_decimal))]
                                ratios_final_avg = round(np.mean(ratios_final_abs_diff), 5) # Don't multiply by 100 as they are already percentages
                                ratios_final_min = round(min(ratios_final_abs_diff), 5) # Don't multiply by 100 as they are already percentages
                                ratios_final_max = round(max(ratios_final_abs_diff), 5) # Don't multiply by 100 as they are already percentages

                                self.print_term("Final deviations from expected ratios:", "Difference in %-values",
                                                spaces=1+printer_spacing,                             )
                                self.print_term("Maximum:", ratios_final_max, spaces=2+printer_spacing)
                                self.print_term("Average:", ratios_final_avg, spaces=2+printer_spacing)
                                self.print_term("Minimum:", ratios_final_min, spaces=2+printer_spacing)
                        
                        ##################################################
                        ### PRINTING SUBLEAFLET SPECIFIC LIPID DETAILS ###
                        ##################################################
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

                        if self.extra_info and len(leaflet["subleaflets"]) > 1 and self.verbose >= 5:
                            self.print_term(
                                printer_leaf_name, "specific lipid data",
                                "(Max lipids: "+str(max_lipids_possible)+", Final lipids: "+str(sum(lipid_ratios))+")",
                                spaces=1+printer_spacing,
                                verbose=4,

                            )
                            for i, (L0, L1, L2, L3, L4, L5) in enumerate(list(zip(*lipids_for_printer))):
                                self.print_term(
                                    '{0: <{L}}'.format(L0, L = max_lengths[0]), ":",
                                    '{0: <{L}}'.format(L1, L = max_lengths[1]), ":",
                                    '{0: <{L}}'.format(L2, L = max_lengths[2]), ":",
                                    '{0: <{L}}'.format(L3, L = max_lengths[3]), ":",
                                    '{0: <{L}}'.format(L4, L = max_lengths[4]), ":",
                                    '{0: <{L}}'.format(L5, L = max_lengths[5]),
                                    spaces=2+printer_spacing,
                                    verbose=4,
                                )
                    
                    leaflet["leaf_lipid_count"] = []
                    for n, c in leaflet["leaf_lipid_count_dict"].items():
                        leaflet["leaf_lipid_count"].append((n, c))
                        if n not in membrane_lipid_count_dict:
                            membrane_lipid_count_dict[n] = c
                        else:
                            membrane_lipid_count_dict[n] += c
                    
                    ###############################################
                    ### PRINTING LEAFLET SPECIFIC LIPID DETAILS ###
                    ###############################################
                    if self.extra_info:
                        self.print_term("", verbose=4)
                        LEAF_headers = ["Lipid name", "Total lipids", "Total %"]
                        LEAF_lipid_names, LEAF_lipid_vals = zip(*leaflet["leaf_lipid_count"])
                        LEAF_tot_lipids = sum(LEAF_lipid_vals)
                        LEAF_lipid_percentages = [round(val / LEAF_tot_lipids * 100, 3) for val in LEAF_lipid_vals]

                        LEAF_lipids_for_printer = list(zip(LEAF_lipid_names, LEAF_lipid_vals, LEAF_lipid_percentages))
                        LEAF_lipids_for_printer = [tuple(LEAF_headers)] + LEAF_lipids_for_printer
                        LEAF_max_lengths = [len(str(max(data, key = lambda d: len(str(d))))) for data in zip(*LEAF_lipids_for_printer)]

                        self.print_term("Leaflet specific lipid data (Combined subleafs)", spaces=2, verbose=4)
                        for L0, L1, L2 in LEAF_lipids_for_printer:
                            self.print_term(
                                '{0: <{L}}'.format(L0, L = LEAF_max_lengths[0]), ":",
                                '{0: <{L}}'.format(L1, L = LEAF_max_lengths[1]), ":",
                                '{0: <{L}}'.format(L2, L = LEAF_max_lengths[2]),
                                spaces=2+printer_spacing,
                                verbose=4,
                            )

                ################################################
                ### PRINTING MEMBRANE SPECIFIC LIPID DETAILS ###
                ################################################
                if self.extra_info and len(memb_dict["leaflets"]) > 1:
                    self.print_term("", verbose=3)
                    membrane_lipid_count_tuples = [(n, c) for n, c in membrane_lipid_count_dict.items()]
                    MEMB_headers = ["Lipid name", "Total lipids", "Total %"]
                    MEMB_lipid_names, MEMB_lipid_vals = zip(*membrane_lipid_count_tuples)
                    MEMB_tot_lipids = sum(MEMB_lipid_vals)
                    MEMB_lipid_percentages = [round(val / MEMB_tot_lipids * 100, 3) for val in MEMB_lipid_vals]

                    MEMB_lipids_for_printer = list(zip(MEMB_lipid_names, MEMB_lipid_vals, MEMB_lipid_percentages))
                    MEMB_lipids_for_printer = [tuple(MEMB_headers)] + MEMB_lipids_for_printer
                    MEMB_max_lengths = [len(str(max(data, key = lambda d: len(str(d))))) for data in zip(*MEMB_lipids_for_printer)]

                    self.print_term("Membrane-wide lipid data", spaces=1, verbose=3)
                    for L0, L1, L2 in MEMB_lipids_for_printer:
                        self.print_term(
                            '{0: <{L}}'.format(L0, L = MEMB_max_lengths[0]), ":",
                            '{0: <{L}}'.format(L1, L = MEMB_max_lengths[1]), ":",
                            '{0: <{L}}'.format(L2, L = MEMB_max_lengths[2]),
                            spaces=2+printer_spacing,
                            verbose=3,
                        )
            
            #########################################
            ### PRINTING SYSTEMWIDE LIPID DETAILS ###
            #########################################
            if self.extra_info and len(self.MEMBRANES) > 1:
                self.print_term("", verbose=2)
                SYS_headers = ["Lipid name", "Total lipids", "Total %"]
                SYS_lipid_names = [name for name in self.SYS_lipids_dict.keys()]
                SYS_lipid_vals = [val for val in self.SYS_lipids_dict.values()]
                SYS_tot_lipids = sum(SYS_lipid_vals)
                SYS_lipid_percentages = [round(val / SYS_tot_lipids * 100, 3) for val in SYS_lipid_vals]

                SYS_lipids_for_printer = list(zip(SYS_lipid_names, SYS_lipid_vals, SYS_lipid_percentages))
                SYS_lipids_for_printer = [tuple(SYS_headers)] + SYS_lipids_for_printer
                SYS_max_lengths = [len(str(max(data, key = lambda d: len(str(d))))) for data in zip(*SYS_lipids_for_printer)]
                
                self.print_term("\nLipid data for whole system", spaces=0, verbose=3)
                for L0, L1, L2 in SYS_lipids_for_printer:
                    self.print_term(
                        '{0: <{L}}'.format(L0, L = SYS_max_lengths[0]), ":",
                        '{0: <{L}}'.format(L1, L = SYS_max_lengths[1]), ":",
                        '{0: <{L}}'.format(L2, L = SYS_max_lengths[2]),
                        spaces=2+printer_spacing,
                        verbose=3,
                    )
            
            lipid_calculator_toc = time.time()
            lipid_calculator_time = round(lipid_calculator_toc - lipid_calculator_tic, 4)
            self.print_term("------------------------------ LIPID RATIO CALCULATIONS COMPLETE", "(time spent: "+str(lipid_calculator_time)+" [s])", "\n", spaces=0, verbose=1)
    
    #######################################
    ### PLANAR GRID MAKER AND OPTIMIZER ###
    #######################################
    def planar_grid_maker(self):
        self.SYS_lipids_dict = {}
        self.GRID_PLOTTING = {}
        if len(self.MEMBRANES) != 0:
            grid_making_tic = time.time()
            self.print_term("------------------------------ CREATING LIPID GRID", spaces=0, verbose=1)
            for memb_i, (memb_key, memb_dict) in enumerate(self.MEMBRANES.items()):
                if memb_i != 0:
                    self.print_term("", verbose=3)
                self.print_term("Starting membrane nr", memb_key, spaces=0, verbose=2)
                for leaflet_i, (leaflet_key, leaflet) in enumerate(memb_dict["leaflets"].items()):
                    if leaflet_i != 0:
                        self.print_term("", verbose=4)
                    self.print_term("Starting", leaflet["leaflet_type"], "leaflet", spaces=1, verbose=3)
    
                    
                    if (self.verbose == 4 and len(leaflet["subleaflets"]) > 1) or (self.verbose > 4 and len(leaflet["subleaflets"]) == 1):
                        self.print_term("Starting optimization", spaces=2)
                    
                    for sli, ((slxi, slyi, sli), subleaflet) in enumerate(leaflet["subleaflets"].items()):
                        if self.verbose > 4 and len(leaflet["subleaflets"]) > 1:
                            self.print_term("Starting optimization for subleaflet nr", sli+1, spaces=2, verbose=5)
                        if self.plot_grid:
                            self.GRID_PLOTTING[(memb_key, leaflet_key, slxi, slyi, sli)] = {}
                        lipid_names_nlipids_radii = [
                            (
                                name,
                                ratio,
                                leaflet["lipids"][name].get_radius(AXs="xy") + max([leaflet["kickx"], leaflet["kicky"]]) + leaflet["plane_buffer"]
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
                        
                        if occupation_modifier < 0.2:
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
                        
                        grid_points, grid_points_no_random, dict_for_plotting = self.make_rect_grid(leaflet, subleaflet, occupation_modifier)
                        grid_points_arr = np.asarray(grid_points)
                        if self.plot_grid:
                            for key, vals in dict_for_plotting.items():
                                self.GRID_PLOTTING[(memb_key, leaflet_key, slxi, slyi, sli)][key] = vals
                        
                        grid_z_value = leaflet["center"][2] + (sign * leaflet["lipid_dimensions"]["zUnderLipCen"])
                        z_values_arr = np.asarray([[grid_z_value] for _ in range(len(grid_points_arr))])
                        
                        if leaflet["optimize"]:
                            
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
                                for key in ["prot_points", "prot_Points_buffered", "prot_Points_buffered_union", "prot_poly_points", "alphashape_1", "ConcaveHulls_Polygon_1", "protein_poly", "remove_union", "require_union"]:
                                    if key in leaflet and leaflet[key]:
                                        self.GRID_PLOTTING[(memb_key, leaflet_key, slxi, slyi, sli)].update({
                                            key: leaflet[key]
                                        })
                                for key in ["original_holed_bbox"]:
                                    if key in subleaflet and subleaflet[key]:
                                        self.GRID_PLOTTING[(memb_key, leaflet_key, slxi, slyi, sli)].update({
                                            key: subleaflet[key]
                                        })
                                    
                                self.GRID_PLOTTING[(memb_key, leaflet_key, slxi, slyi, sli)].update({
                                    ### Inputs
                                    "lipids"      : lipids,
                                    "bbox_polygon": subleaflet["holed_bbox"],
                                    
                                    "xdims"          : (subleaflet["xmin"], subleaflet["xmax"]),
                                    "ydims"          : (subleaflet["ymin"], subleaflet["ymax"]),
                                    "xdims_original" : (subleaflet["xmin_original"], subleaflet["xmax_original"]),
                                    "ydims_original" : (subleaflet["ymin_original"], subleaflet["ymax_original"]),
                                    "push_mult"      : leaflet["optim_push_mult"],
                                    "apl"            : leaflet["apl"],
                                    "bin_size"       : bin_size,

                                    ### Outputs
                                    "grid_points"          : grid_points_arr,
                                    "grid_points_no_random": grid_points_no_random,
                                    "points_arr"           : optimized_grid_points_arr,
                                    "POINT_STEPS"          : POINT_STEPS,
                                    "step"                 : step,
                                    "steps_time"           : steps_time,
                                    "mean_steps_time"      : mean_steps_time,
                                    "max_push"             : max_push,
                                })
                        else:
                            subleaflet["grid_points"] = np.hstack([grid_points_arr, z_values_arr])
            grid_making_toc = time.time()
            grid_making_time = round(grid_making_toc - grid_making_tic, 4)
            self.print_term("------------------------------ LIPID GRID CREATED", "(time spent: "+str(grid_making_time)+" [s])", "\n", spaces=0, verbose=1)
    
    def make_rect_grid(self, leaflet, subleaflet, occupation_modifier):
        xmin, xmax, ymin, ymax = itemgetter("xmin", "xmax", "ymin", "ymax")(subleaflet)
        
        sidelen      = math.sqrt(leaflet["apl"])
        ### "edge_buffer" uses "occupation_modifier*2" while "mean_lipid_radius" uses only "occupation_modifier"
        ### Done to ensure lipids are allowed to be placed right near the border
        edge_buffer  = (leaflet["lipid_dimensions"]["lipid_radius"] + max([leaflet["kickx"], leaflet["kicky"]]) + leaflet["plane_buffer"]) * (1+occupation_modifier*2)
        mean_lipid_radius = (leaflet["lipid_dimensions"]["lipid_radius"] + max([leaflet["kickx"], leaflet["kicky"]]) + leaflet["plane_buffer"]) * (1+occupation_modifier)
        bbox_polygon = subleaflet["holed_bbox"]
        lipids       = subleaflet["lipids"]
        
        xmin_edge = xmin + edge_buffer
        xmax_edge = xmax - edge_buffer
        ymin_edge = ymin + edge_buffer
        ymax_edge = ymax - edge_buffer
        
        ### Making pointer for when number of points along x/y-axis should be expanded
        xlines_ideal = (xmax_edge-xmin_edge)/sidelen
        ylines_ideal = (ymax_edge-ymin_edge)/sidelen
        
        ratio_tot_ideal = xlines_ideal + ylines_ideal
        xlines_ratio_ideal = xlines_ideal/ratio_tot_ideal
        ylines_ratio_ideal = ylines_ideal/ratio_tot_ideal
        
        xlines = int(xlines_ideal)
        ylines = int(ylines_ideal)
        
        while xlines * ylines < len(lipids):
            new_ratio_tot = xlines + ylines
            if xlines / new_ratio_tot < xlines_ratio_ideal:
                xlines += 1
            elif ylines / new_ratio_tot < ylines_ratio_ideal:
                ylines += 1
            else:
                xlines += 1
        
        bbox_polygon_BufferedForLineStrings = shapely.buffer(bbox_polygon, -mean_lipid_radius)
        dict_for_plotting = {"bbox_polygon_BufferedForLineStrings": bbox_polygon_BufferedForLineStrings}
        
        enough_points = False
        while enough_points == False:
            xlinespace, xspace = np.linspace(start=xmin_edge, stop=xmax_edge, num=xlines, endpoint=True, retstep=True)
            ylinespace, yspace = np.linspace(start=ymin_edge, stop=ymax_edge, num=ylines, endpoint=True, retstep=True)
            
            LineStrings            = []
            LineStringsOverlapping = []
            LineStringsContained   = []
            
            c = 0
            for xval in xlinespace:
                top_point = (xval, ylinespace[-1])
                bot_point = (xval, ylinespace[0])
                LineString = shapely.LineString([top_point, bot_point])
                LineStrings.append(LineString)
                
                ### Finding parts of LineStrings that are contained within legal area and at least a certain distance from BBOX
                LineStringContained = shapely.intersection(LineString, bbox_polygon_BufferedForLineStrings)
                ### If LineStrings are cut into multiple smaller LineStrings
                if LineStringContained.geom_type == "MultiLineString":
                    momentary_LineStrings = [LineString for LineString in LineStringContained.geoms]
                else:
                    momentary_LineStrings = [LineStringContained]
                ### Removes "empty" LineStrings
                momentary_LineStrings = [LineString for LineString in momentary_LineStrings if not LineString.is_empty]
                LineStringsOverlapping.extend(momentary_LineStrings)
                
                ### Checks if any endpoints on LineStrings are too close to each other
                ### Cuts a portion of each if they are
                new_momentary_LineStrings = []
                for LS, LineString in enumerate(momentary_LineStrings):
                    curr_xs, curr_ys = LineString.xy
                    curr_xmax, curr_xmin = max(curr_xs), min(curr_xs)
                    curr_ymax, curr_ymin = max(curr_ys), min(curr_ys)
                    curr_length = curr_ymax - curr_ymin
                    if len(new_momentary_LineStrings) > 0:
                        y_diff = last_ymin - curr_ymax
                        if y_diff <= yspace:
                            tot_LS_lengths  = curr_length + last_length
                            curr_LS_portion = curr_length / tot_LS_lengths
                            last_LS_portion = last_length / tot_LS_lengths

                            yspace_diff = yspace - y_diff
                            curr_cut = curr_LS_portion * yspace_diff
                            last_cut = last_LS_portion * yspace_diff
                            if last_cut >= last_length:
                                ### Remove last LineString if it is too short
                                new_momentary_LineStrings = new_momentary_LineStrings[:-1]
                            else:
                                ### Else cut a part of it off and add it to list
                                last_ymin += last_cut
                                last_top_point, last_bot_point = (last_xmax, last_ymax), (last_xmin, last_ymin)
                                new_last_LS = shapely.LineString([last_top_point, last_bot_point])
                                new_momentary_LineStrings[-1] = new_last_LS
                            if curr_cut >= curr_length:
                                ### If current LineString too short, then just don't do anything
                                pass
                            else:
                                ### Else cut a part of it off and add it to list
                                curr_ymax -= curr_cut
                                curr_top_point, curr_bot_point = (curr_xmax, curr_ymax), (curr_xmin, curr_ymin)
                                new_curr_LS = shapely.LineString([curr_top_point, curr_bot_point])
                                new_momentary_LineStrings.append(new_curr_LS)
                        else:
                            new_momentary_LineStrings.append(LineString)
                    else:
                        new_momentary_LineStrings.append(LineString)
                    
                    ### Could happen that all LineStrings have been removed due to being too short.
                    if len(new_momentary_LineStrings) > 0:
                        last_xs, last_ys = new_momentary_LineStrings[-1].xy
                        last_xmax, last_xmin = max(last_xs), min(last_xs)
                        last_ymax, last_ymin = max(last_ys), min(last_ys)
                        last_length = last_ymax - last_ymin
                
                LineStringsContained.extend(new_momentary_LineStrings)
            
            dict_for_plotting["LineStrings"]            = LineStrings
            dict_for_plotting["LineStringsOverlapping"] = LineStringsOverlapping
            dict_for_plotting["LineStringsContained"]   = LineStringsContained
            
            ngridpoints = 0
            lines_info = []
            for LineString in LineStringsContained:
                line_length = LineString.length
                npoints_on_line = round(line_length // yspace) + 1
                ngridpoints += npoints_on_line
                
                if npoints_on_line > 1:
                    real_yspace = line_length / (npoints_on_line-1)
                else:
                    real_yspace = 0
                
                xs, ys = LineString.xy
                xval = xs[0] # both x-values are the same
                ymin, ymax = min(ys), max(ys)
                
                lines_info.append({
                    "LineString": LineString,
                    "xval": xval,
                    "ymin": ymin,
                    "ymax": ymax,
                    "line_length": line_length,
                    "npoints_on_line": npoints_on_line,
                    "real_yspace": real_yspace,
                })

            if ngridpoints < len(lipids):
                new_ratio_tot = xlines + ylines
                if xlines / new_ratio_tot < xlines_ratio_ideal:
                    xlines += 1
                elif ylines / new_ratio_tot < ylines_ratio_ideal:
                    ylines += 1
                else:
                    xlines += 1
            else:
                enough_points = True
        
        ### Finding the points that should be removed due to number of lipids
        ### No duplicate index values as "len(lipids)" is always equal to or smaller than "ngridpoints"
        while ngridpoints > len(lipids):
            lines = [line for line in lines_info if line["real_yspace"] != 0]
            
            line_smallest_yspace = min(lines, key=lambda line: line["real_yspace"])
            
            line_smallest_yspace["npoints_on_line"] -= 1
            
            if line_smallest_yspace["npoints_on_line"] > 1:
                line_smallest_yspace["real_yspace"] = line_smallest_yspace["line_length"] / (line_smallest_yspace["npoints_on_line"]-1)
            else:
                line_smallest_yspace["real_yspace"] = 0
            
            ngridpoints = sum([line["npoints_on_line"] for line in lines_info])
                
        ### Creating the points
        grid_points = []
        grid_points_no_random = []
        rand_force = 1
        for line in lines_info:
            xval = line["xval"]
            if line["real_yspace"] == 0:
                yvals = [np.mean([line["ymin"], line["ymax"]])]
            else:
                yvals = np.linspace(start=line["ymin"], stop=line["ymax"], num=line["npoints_on_line"], endpoint=True)
            
            for yval in yvals:
                point = (xval, yval)
                randx = random.uniform(-leaflet["kickx"]/rand_force, leaflet["kickx"]/rand_force)
                randy = random.uniform(-leaflet["kicky"]/rand_force, leaflet["kicky"]/rand_force)
                grid_points.append((xval+randx, yval+randy))
                grid_points_no_random.append((xval, yval))

#         return grid_points
        return grid_points, grid_points_no_random, dict_for_plotting
    
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
        
        self.print_term("CURRENT STEP:", end=" ", spaces=3, verbose=5)
        
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
            
            self.print_term(step, end=" ", spaces=0, verbose=5)
            
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
            
            ### Lipid-lipid pushes
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

                eq_dist = lipid_sizes[pi1]*(1+occupation_modifier*2)
#                 eq_dist = lipid_sizes[pi1]*(1+occupation_modifier)
                
                ### Lipid contained in legal areas but close to edge
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
                
                ### Lipid outside legal areas
                ### Push it just enough for it to be contained in legal areas
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
        self.print_term("", verbose=5)
        self.print_term("Last step:         ", step, spaces=3, verbose=5)
        self.print_term("Optimization time: ", round(steps_time, 4), "[s]", spaces=3, verbose=5)
        self.print_term("Mean step time:    ", round(steps_time/end, 4), "[s]", spaces=3, verbose=5)
        
        return points_arr, POINT_STEPS, step, steps_time, mean_steps_time, max_push, bin_size

    ######################
    ### LIPID INSERTER ###
    ######################
    def lipid_inserter(self):
        if len(self.MEMBRANES) != 0:
            lipid_inserter_tic = time.time()
            self.print_term("------------------------------ CREATING LIPIDS", spaces=0, verbose=1)
            for memb_key, memb_dict in self.MEMBRANES.items():
                self.print_term("Starting membrane nr", memb_key, spaces=0, verbose=2)
                for leaflet_i, (leaflet_key, leaflet) in enumerate(memb_dict["leaflets"].items()):
                    self.print_term("Starting", leaflet["leaflet_type"], "leaflet", spaces=1, verbose=3)
                        
                    if leaflet["HG_direction"] == "up":
                        sign = +1
                    if leaflet["HG_direction"] == "down":
                        sign = -1
                    leaflet["grid_lipids"] = []
                    
                    for (slxi, slyi, sli), subleaflet in leaflet["subleaflets"].items():
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
                                "charges": [bead.charge for bead in leaflet["lipids"][l_name].get_res_beads_info()],
                            }
                            leaflet["grid_lipids"][-1].update({"lipid": lipid_dict})
                            self.system_charge += leaflet["lipids"][l_name].get_mol_charge()
            lipid_inserter_toc = time.time()
            lipid_inserter_time = round(lipid_inserter_toc - lipid_inserter_tic, 4)
            self.print_term("------------------------------ LIPID CREATION COMPLETE", "(time spent: "+str(lipid_inserter_time)+" [s])", "\n", spaces=0, verbose=1)

    ################
    ### SOLVATER ###
    ################
    def get_solute_volume(self, cxmin, cxmax, cymin, cymax, czmin, czmax):
        solvent_beads_for_cell_checker = []
        solvent_box_solute_charge = 0
        if len(self.SOLVATIONS) > 0:
            '''
            Finds the solvents beads and estimates their volume
            '''
            for solv_key, solv_dict in self.SOLVATIONS.items():
                if "grid" in solv_dict:
                    for grid_point in solv_dict["grid"]:
                        beads, charges = itemgetter("coords", "bead_charges")(grid_point)
                        for x, y, z in beads:
                            if cxmin < x <= cxmax and cymin < y <= cymax and czmin < z <= czmax:
                                solvent_beads_for_cell_checker.append((x, y, z))
                                solvent_box_solute_charge += sum(charges)
        return solvent_beads_for_cell_checker, solvent_box_solute_charge

    def get_lipid_volume(self, cxmin, cxmax, cymin, cymax, czmin, czmax):
        lipid_beads_for_cell_checker = []
        solvent_box_lipids_charge = 0
        if len(self.MEMBRANES) > 0:
            '''
            Finds the lipids beads and estimates their volume
            '''
            for memb_key, memb_dict in self.MEMBRANES.items():
                for leaflet_key, leaflet in memb_dict["leaflets"].items():
                    for grid_point in leaflet["grid_lipids"]:
                        xs, ys, zs, charges = itemgetter("x", "y", "z", "charges")(grid_point["lipid"])
                        for x, y, z, charge in zip(xs, ys, zs, charges):
                            if cxmin < x <= cxmax and cymin < y <= cymax and czmin < z <= czmax:
                                lipid_beads_for_cell_checker.append((x, y, z))
                                solvent_box_lipids_charge += charge
        return lipid_beads_for_cell_checker, solvent_box_lipids_charge

    def get_protein_volume(self, cxmin, cxmax, cymin, cymax, czmin, czmax):
        prot_beads_for_cell_checker = []
        solvent_box_proteins_charge = 0
        if len(self.PROTEINS) > 0:
            '''
            Finds the proteins beads and estimates their volume
            '''
            for protein_i, (protein_nr, protein) in enumerate(self.PROTEINS.items()):
                beads = protein["protein"].get_beads("xyz")
                charges = protein["protein"].get_bead_charges()
                for (x, y, z), charge in zip(beads, charges):
                    if cxmin < x <= cxmax and cymin < y <= cymax and czmin < z <= czmax:
                        prot_beads_for_cell_checker.append((x, y, z))
                        solvent_box_proteins_charge += charge
        return prot_beads_for_cell_checker, solvent_box_proteins_charge
    
    def solvater(self):
        if len(self.SOLVATIONS_cmds) != 0:
            solvation_tic = time.time()
            self.print_term("------------------------------ SOLVATING SYSTEM", spaces=0, verbose=1)
            solv_beads_for_cell_checker = []
            
            for solvation_i, (solvation_nr, solvation) in enumerate(self.SOLVATIONS.items()):
                if solvation_i != 0:
                    self.print_term("", verbose=3)
                self.print_term("Starting solvation nr", solvation_nr, spaces=0, verbose=2)

                #########################################
                ### CALCULATION FREE VOLUME OF SYSTEM ###
                #########################################
                self.print_term("Calculating box volume: (all values in [nm^3])", spaces=1, verbose=4)
                self.print_term("Bead radius used for volume calculations 'bead_radius':", solvation["bead_radius"]/10, "[nm]", spaces=2, verbose=4)

                bead_radius = solvation["bead_radius"]
                bead_volume = (4/3 * math.pi * (bead_radius ** 3)) * 10**-27
                
                cx, cy, cz = solvation["center"]
                xlen, ylen, zlen = itemgetter("x", "y", "z")(solvation)
                cxmin, cxmax = cx-xlen/2, cx+xlen/2
                cymin, cymax = cy-ylen/2, cy+ylen/2
                czmin, czmax = cz-zlen/2, cz+zlen/2
                
                self.print_term("cx, cy, cz:", cx, cy, cz, debug=True)
                self.print_term("xlen, ylen, zlen:", xlen, ylen, zlen, debug=True)
                self.print_term("cxmin, cxmax:", cxmin, cxmax, debug=True)
                self.print_term("cymin, cymax:", cymin, cymax, debug=True)
                self.print_term("czmin, czmax:", czmin, czmax, debug=True)
                
                ### Solvent/solute volume
                solvent_beads_for_cell_checker, solvent_box_solute_charge = self.get_solute_volume(
                    cxmin, cxmax, cymin, cymax, czmin, czmax
                )
                solvs_volume = bead_volume * len(solvent_beads_for_cell_checker)
                
                ### Lipid volume
                lipid_beads_for_cell_checker, solvent_box_lipids_charge = self.get_lipid_volume(
                    cxmin, cxmax, cymin, cymax, czmin, czmax
                )
                leafs_volume = bead_volume * len(lipid_beads_for_cell_checker)

                ### Protein volume
                prot_beads_for_cell_checker, solvent_box_proteins_charge = self.get_protein_volume(
                    cxmin, cxmax, cymin, cymax, czmin, czmax
                )
                prots_volume = bead_volume * len(prot_beads_for_cell_checker)
                
                n_lipids = 0
                if solvation["solv_per_lipid"]:
                    if len(self.MEMBRANES) > 0:
                        ### Counts the number of lipids inside the solvent box
                        for memb_key, memb_dict in self.MEMBRANES.items():
                            for leaflet_key, leaflet in memb_dict["leaflets"].items():
                                for grid_point in leaflet["grid_lipids"]:
                                    xs, ys, zs = itemgetter("x", "y", "z")(grid_point["lipid"])
                                    beads_total = len(xs)
                                    beads_inside = 0
                                    for x, y, z in zip(xs, ys, zs):
                                        if cxmin < x <= cxmax and cymin < y <= cymax and czmin < z <= czmax:
                                            beads_inside += 1
                                    if beads_inside > beads_total/2:
                                        n_lipids += 1
                
                solvent_box_charge = solvent_box_solute_charge + solvent_box_lipids_charge + solvent_box_proteins_charge
                
                N_A = 6.02214076 * 10**23
                box_volume = (xlen * ylen * zlen) * 10**-27
                non_free_volume = leafs_volume + prots_volume + prots_volume
                box_free_volume = box_volume - non_free_volume
                self.print_term("Solvent box volume:", round(box_volume      * 10**24, 3), spaces=2, verbose=4)
                self.print_term("Excluded volume:   ", round(non_free_volume * 10**24, 3), spaces=2, verbose=4)
                
                self.print_term("Solute volume: ", round(solvs_volume    * 10**24, 3), spaces=3, verbose=5)
                self.print_term("Lipid volume:  ", round(leafs_volume    * 10**24, 3), spaces=3, verbose=5)
                self.print_term("Protein volume:", round(prots_volume    * 10**24, 3), spaces=3, verbose=5)
                
                self.print_term("Free volume:       ", round(box_free_volume * 10**24, 3), spaces=2, verbose=4)
                
                ###########################
                ### CALCULATING SOLVENT ###
                ###########################

                ### List used to find the maximum solvent/ion size
                solvent_radii = []
                                
                solvent_radii = [
                    vals.get_radius()
                    for dict_pointer in ["solvent", "pos_ions", "neg_ions"]
                    for key, vals in solvation[dict_pointer].items()
                ]

                def get_solvent_ratios(solvent_type_dict, tot_ratio):
                    ratios = [
                        round(1 / tot_ratio * vals.ratio, 4)
                        for key, vals in solvent_type_dict.items()
                    ]
                    return ratios
                
                if not solvation["count"]:
                    sol_ratios = get_solvent_ratios(solvation["solvent"], solvation["solv_tot_ratio"])
                    pos_ratios = get_solvent_ratios(solvation["pos_ions"], solvation["pos_tot_ratio"])
                    neg_ratios = get_solvent_ratios(solvation["neg_ions"], solvation["neg_tot_ratio"])
                else:
                    sol_ratios = [1 for _ in solvation["solvent"].keys()]
                    pos_ratios = [1 for _ in solvation["pos_ions"].keys()]
                    neg_ratios = [1 for _ in solvation["neg_ions"].keys()]

#                 def solvent_count_calculator(solvent_type_dict, ratios, count_bool, solv_per_lipid, n_lipids, volume):
                def solvent_count_calculator(solvent_type, ratios, volume):
                    counts = []
                    for (key, vals), ratio in zip(solvation[solvent_type].items(), ratios):
                        if solvation["count"]:
                            ### Treats molarity as absolute number of molecules
                            counts.append(vals.molarity)
                        elif solvation["solv_per_lipid"] and solvent_type == "solvent":
                            counts.append(round(solvation["solv_per_lipid"] * n_lipids * ratio))
                        else:
                            ### Treats molarity as molarity
                            counts.append(int(N_A * volume * vals.molarity / vals.mapping_ratio * ratio))
                    return counts
                
                ### Using the free volume for concentration calculations (default)
                if solvation["solvfreevol"] == True:
                    volume_for_solv = box_free_volume
                ### Using box volume for concentration calculations
                elif solvation["solvfreevol"] == False:
                    volume_for_solv = box_volume
                
                sol_counts = solvent_count_calculator(
                    "solvent",
                    sol_ratios,
                    volume_for_solv,
                )
                
                solv_volume = 0
                sol_charges = 0

                for (key, vals), count in zip(solvation["solvent"].items(), sol_counts):
                    vals.count_set(count)
                    sol_charges += vals.count * vals.get_mol_charge()
                    if not vals.molar_mass and solvation["ionsvol"] == "solv" and not solvation["count"]:
                        self.print_term("WARNING: Chosen solvent [" + key + "] is missing a 'molar_mass' entry in the solvent defines. Please add it if you wish to use solvent volume for ion concentration.", warn = True)
                    if not vals.density and solvation["ionsvol"] == "solv" and not solvation["count"]:
                        self.print_term("WARNING: Chosen solvent [" + key + "] is missing a 'density' entry in the solvent defines. Please add it if you wish to use solvent volume for ion concentration.", warn = True)
                    if vals.molar_mass and vals.density:
                        solv_volume += (count * vals.mapping_ratio * vals.molar_mass) / (N_A * vals.density) * 10**-3
                
                self.print_term("Solvent volume:    ", round(solv_volume * 10**24, 3), spaces=2, verbose=4)

                def ions_neutraliser(pos_ions, neg_ions, target_charge, current_charge, direction):
                    charge_difference = target_charge - current_charge
                    if charge_difference > 0:
                        ions = pos_ions
                    elif charge_difference < 0:
                        ions = neg_ions
                    else:
                        return pos_ions, neg_ions
                    
                    if direction == "add":
                        sign = +1
                    if direction == "remove":
                        sign = -1
                        
                    keys = list(ions.keys())
                    keys_sorted_after_ratio = list(sorted(keys, key=lambda key: ions[key]["ratio"], reverse = True))
                    tot_ions = sum([vals["count"] for vals in ions.values()])
                    counter = 0
                    while round(current_charge) != round(target_charge):
                        counter += 1
                        ### Find which ion is furthest from ideal ratio
                        if tot_ions == 0:
                            ### No ions yet, just sort them from highest to lowest ratio
                            furthest_dist_keys = keys_sorted_after_ratio[::sign]
                        else:
                            ### Calculate how far each ion is from its ideal ratio
                            dists_from_ratio = []
                            for key, vals in list(ions.items()):
                                dists_from_ratio.append((key, vals["count"]/tot_ions - vals["ratio"]))
                            ### Sort them from most underrepresented to most overrepresented
                            ### Secondarily sort them according to their ideal ratios in case the are equally underrepresented
                            dists_from_ratio_sorted = sorted(
                                dists_from_ratio,
                                key=lambda x: (x[1], keys_sorted_after_ratio.index(x[0])),
                                reverse = False
                            )
                            ### Get just the keys
                            furthest_dist_keys = [key for key, dist_from_ratio in dists_from_ratio_sorted][::sign]
                        
                        ### Checks from most underrepresented to most overrepresented if adding another ion
                        ### would cause the total added charge to become greater that then difference that is being corrected
                        for key in furthest_dist_keys:
                            if abs(target_charge - current_charge) >= abs(ions[key]["charge"]):
                                furthest_dist_key = key
                                break
                        else:
                            ### If an ion should be added but no ion can be added without adding
                            ### more charge than is needed for neutralization
                            ### then just return as is
                            return pos_ions, neg_ions
                        
                        tot_ions += sign
                        current_charge += ions[furthest_dist_key]["charge"]
                        ions[furthest_dist_key]["count"] += sign
                        
                    return pos_ions, neg_ions
    
                ### Using the solvent volume for concentration calculations (default)
                if solvation["ionsvol"] == "solv" and solv_volume > 0 and not solvation["count"]:
                    volume_for_ions = solv_volume
                ### Using the free volume for concentration calculations
                elif solvation["ionsvol"] == "free" or (solvation["ionsvol"] == "solv" and solvation["count"]):
                    volume_for_ions = box_free_volume
                ### Using the box volume for concentration calculations
                elif solvation["ionsvol"] == "box":
                    volume_for_ions = box_volume

                ### First ensures that the ion concentration is equal to salt_molarity
                pos_counts = solvent_count_calculator(
                    "pos_ions",
                    pos_ratios,
                    volume_for_ions
                )
                neg_counts = solvent_count_calculator(
                    "neg_ions",
                    neg_ratios,
                    volume_for_ions
                )
                
                pos_charges = 0
                neg_charges = 0
                
                ions_dict = {}
                zipped1 = zip(["pos_ions", "neg_ions"], [pos_counts, neg_counts], [pos_ratios, neg_ratios])
                for ion_type, ion_counts, ion_ratios in zipped1:
                    ions_dict[ion_type] = {}
                    zipped2 = zip(solvation[ion_type].items(), ion_counts, ion_ratios)
                    for (key, vals), count, ratio in zipped2:
                        vals.count_set(count)
                        if ion_type == "pos_ions":
                            pos_charges += vals.count * vals.get_mol_charge()
                        elif ion_type == "neg_ions":
                            neg_charges += vals.count * vals.get_mol_charge()
                        ions_dict[ion_type][key] = {
                            "charge": vals.get_mol_charge(),
                            "ratio": ratio,
                            "count": vals.count,
                        }
                
                current_charge = solvent_box_charge + sol_charges + pos_charges + neg_charges
                
                pos_charges = 0
                neg_charges = 0
                
                if solvation["salt_method"] == "add":
                    '''
                    Adds extra ions to neutralize the solvent box
                    '''
                    pos_ions, neg_ions = ions_neutraliser(
                        ions_dict["pos_ions"],
                        ions_dict["neg_ions"],
                        solvation["charge"],
                        current_charge,
                        direction = "add",
                    )
                    zipped = zip(["pos_ions", "neg_ions"],
                                 [pos_ions, neg_ions])
                    for ion_type, ions_dicts in zipped:
                        for (key, vals), (ions_dict_key, ions_dict_vals) in zip(solvation[ion_type].items(), ions_dicts.items()):
                            vals.count_set(ions_dict_vals["count"])
                            if ion_type == "pos_ions":
                                pos_charges += vals.count * vals.get_mol_charge()
                            elif ion_type == "neg_ions":
                                neg_charges += vals.count * vals.get_mol_charge()

                elif solvation["salt_method"] == "remove":
                    '''
                    Removes excess ions to neutralize the solvent box
                    '''
                    pos_ions, neg_ions = ions_neutraliser(
                        ions_dict["pos_ions"],
                        ions_dict["neg_ions"],
                        solvation["charge"],
                        -current_charge,
                        direction = "remove",
                    )
                    zipped = zip(["pos_ions", "neg_ions"],
                                 [pos_ions, neg_ions])
                    for ion_type, ions_dicts in zipped:
                        for (key, vals), (ions_dict_key, ions_dict_vals) in zip(solvation[ion_type].items(), ions_dicts.items()):
                            vals.count_set(ions_dict_vals["count"])
                            if ion_type == "pos_ions":
                                pos_charges += vals.count * vals.get_mol_charge()
                            elif ion_type == "neg_ions":
                                neg_charges += vals.count * vals.get_mol_charge()
                    
                elif solvation["salt_method"] == "mean":
                    '''
                    Adds ions with either a positive or negative charge and removes the other type
                    Ends up with ion concentrations in between "add" and "remove"
                    Does an extra ion addition at the end to ensure no errors have been made when
                        adding/removing ions.
                    '''
                    ### Need to deepcopy, otherwise the nested dictionaries will be modified 
                    ions_dict_add = copy.deepcopy(ions_dict)
                    pos_ions_add, neg_ions_add = ions_neutraliser(
                        ions_dict_add["pos_ions"],
                        ions_dict_add["neg_ions"],
                        solvation["charge"],
                        +math.ceil(current_charge/2),
                        direction = "add",
                    )
                    ions_dict_remove = copy.deepcopy(ions_dict)
                    pos_ions_remove, neg_ions_remove = ions_neutraliser(
                        ions_dict_remove["pos_ions"],
                        ions_dict_remove["neg_ions"],
                        solvation["charge"],
                        -math.floor(current_charge/2),
                        direction = "remove",
                    )
                    
                    zipped1 = zip(["pos_ions", "neg_ions"],
                                  [pos_counts, neg_counts],
                                  [pos_ions_add, neg_ions_add],
                                  [pos_ions_remove, neg_ions_remove])
                    for ion_type, ion_counts, ions_dicts_add, ions_dicts_remove in zipped1:
                        zipped2 = zip(solvation[ion_type].items(),
                                      ion_counts,
                                      ions_dicts_add.values(),
                                      ions_dicts_remove.values())
                        for (key, vals), count, ions_dict_add_vals, ions_dict_remove_vals in zipped2:
                            add_diff = abs(count - ions_dict_add_vals["count"])
                            rem_diff = abs(count - ions_dict_remove_vals["count"])
                            vals.count_set(count + add_diff - rem_diff)
                            if ion_type == "pos_ions":
                                pos_charges += vals.count * vals.get_mol_charge()
                            elif ion_type == "neg_ions":
                                neg_charges += vals.count * vals.get_mol_charge()

                    ### Extra ion addtion to prevent errors (e.g. non-neutralized systems)
                    ions_dict = {}
                    for ion_type, ion_ratios in zip(["pos_ions", "neg_ions"], [pos_ratios, neg_ratios]):
                        ions_dict[ion_type] = {}
                        for (key, vals), ratio in zip(solvation[ion_type].items(), ion_ratios):
                            ions_dict[ion_type][key] = {
                                "charge": vals.get_mol_charge(),
                                "ratio": ratio,
                                "count": vals.count,
                            }

                    current_charge = solvent_box_charge + sol_charges + pos_charges + neg_charges

                    pos_ions, neg_ions = ions_neutraliser(
                        ions_dict["pos_ions"],
                        ions_dict["neg_ions"],
                        solvation["charge"],
                        current_charge,
                        direction = "add",
                    )

                    pos_charges = 0
                    neg_charges = 0

                    zipped = zip(["pos_ions", "neg_ions"], [pos_ions, neg_ions])
                    for ion_type, ions_dicts in zipped:
                        for (key, vals), (ions_dict_key, ions_dict_vals) in zip(solvation[ion_type].items(), ions_dicts.items()):
                            vals.count_set(ions_dict_vals["count"])
                            if ion_type == "pos_ions":
                                pos_charges += vals.count * vals.get_mol_charge()
                            elif ion_type == "neg_ions":
                                neg_charges += vals.count * vals.get_mol_charge()
                
                ### Adding charges from ions and solvent to system charge
                solvent_box_solvent_charge = sol_charges + pos_charges + neg_charges
                self.system_charge += solvent_box_solvent_charge
                
                self.print_term("", verbose=4)
                self.print_term("Solvent box charge information:", spaces=1, verbose=4)
                self.print_term("Solvent box charge before solvent insertion:", round(solvent_box_charge, 1), spaces=2, verbose=4)
                self.print_term("Prior solvent beads: ", round(solvent_box_solute_charge, 1), spaces=3, verbose=5)
                self.print_term("Lipid beads:         ", round(solvent_box_lipids_charge, 1), spaces=3, verbose=5)
                self.print_term("Protein beads:       ", round(solvent_box_proteins_charge, 1), spaces=3, verbose=5)
                solvent_box_charge += solvent_box_solvent_charge
                self.print_term("Solvent box charge after solvent insertion: ", round(solvent_box_charge, 1), spaces=2, verbose=4)
                self.print_term("New solvent beads:", round(solvent_box_solvent_charge, 1), spaces=3, verbose=5)
                self.print_term("Solvent       (solv):", round(sol_charges, 1), spaces=4, verbose=5)
                self.print_term("Positive ions (pos): ", round(pos_charges, 1), spaces=4, verbose=5)
                self.print_term("Negative ions (neg): ", round(neg_charges, 1), spaces=4, verbose=5)

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
                max_mol_size = max(solvent_radii)
                
                self.print_term("gridres first:", solvation["gridres"], debug=True)
                self.print_term("max_mol_size:", max_mol_size, debug=True)
                self.print_term("max_mol_size*1.2:", max_mol_size*1.2, debug=True)
                self.print_term("kick:", solvation["kick"], debug=True)
                self.print_term("kick*1.2:", solvation["kick"]*1.2, debug=True)
                
                gridres = solvation["gridres"]
                
                ### if the maximum molecule size is bigger than the designated grid resolution then change the gridres
                if max_mol_size*1.2 >= gridres:
                    gridres = (max_mol_size + solvation["kick"]*2) * 1.2 # 20% larger than largest molecule
                    self.print_term("NOTE: Requested solvent is too large for the grid resolution. Adjusting grid resolution to prevent solvent molecule overlaps.", warn = True)
                    self.print_term("Original grid resolution was:", round(solvation["gridres"]/10, 4), "[nm]", warn = True)
                    self.print_term("New grid resolution is:      ", round(gridres/10, 4), "[nm]", "\n", warn = True)
                    
                self.print_term("gridres mid:", gridres, debug=True)
                
                if (max_mol_size+solvation["kick"])*1.2 >= gridres:
                    self.print_term("NOTE: Kick is too large for grid resolution. Adjusting grid resolution to prevent solvent molecule overlaps.", warn = True)
                    self.print_term("Original grid resolution was:", round(gridres/10, 4), "[nm]", warn = True)
                    self.print_term("Original kick was:           ", round(solvation["kick"]/10, 4), "[nm]", warn = True)
                    gridres = gridres + solvation["kick"]*2*1.2 # 20% extra
                    self.print_term("New grid resolution is:      ", round(gridres/10, 4), "[nm]", "\n", warn = True)
                
                self.print_term("gridres last:", gridres, debug=True)
                
#                 solvent_buffer = gridres + solvation["buffer"]
                solvent_buffer = solvation["buffer"]

                #################################
                ### CHOOSES ALGORITHM VERSION ###
                #################################
                self.print_term("", verbose=4)
                self.print_term("Calculating the number of available grid points", spaces=1, verbose=4)

                def coord_to_indices(pos, dim, center, real_gridres, max_int, min_buffer, max_buffer):
                    ### Buffer limit
                    bead_min = pos-min_buffer
                    bead_max = pos+max_buffer

                    ### Convert solvent buffer limit to decimal index
                    beadi_min_dec = (bead_min+dim/2-center)/real_gridres
                    beadi_max_dec = (bead_max+dim/2-center)/real_gridres

                    ### Round buffer decimal index to integer index
                    beadi_min_int = math.floor(beadi_min_dec+0.5)
                    beadi_max_int = math.ceil(beadi_max_dec-0.5)

                    if beadi_min_int < 0:
                        beadi_min_int = 0
                    if beadi_max_int < 0:
                        beadi_max_int = 0

                    return beadi_min_int, beadi_max_int
                
                xpoints, ypoints, zpoints = int(xlen/gridres), int(ylen/gridres), int(zlen/gridres)
                ### Calculates actual coordinate ranges for each axis and calculates the "real" grid resolution
                xcoords, xreal_gridres = np.linspace(cxmin-gridres/2, cxmax+gridres/2, xpoints+2, retstep=True)
                ycoords, yreal_gridres = np.linspace(cymin-gridres/2, cymax+gridres/2, ypoints+2, retstep=True)
                zcoords, zreal_gridres = np.linspace(czmin-gridres/2, czmax+gridres/2, zpoints+2, retstep=True)
                ### Removes the first and last points as they are the actual edges of the box
                xcoords = xcoords[1:-1]
                ycoords = ycoords[1:-1]
                zcoords = zcoords[1:-1]
                
                solute_buffer = solvent_buffer+solvation["solute_extra_buffer"]
                lipid_buffer = solvent_buffer+solvation["lipid_extra_buffer"]
                protein_buffer = solvent_buffer+solvation["protein_extra_buffer"]
                
                self.print_term("xpoints, xreal_gridres", xpoints, round(xreal_gridres, 4), debug=True)
                self.print_term("ypoints, yreal_gridres", ypoints, round(yreal_gridres, 4), debug=True)
                self.print_term("zpoints, zreal_gridres", zpoints, round(zreal_gridres, 4), debug=True)
                self.print_term("points total", xpoints*ypoints*zpoints, debug=True)
                self.print_term("solvent_buffer     ", solvent_buffer, debug=True)
                self.print_term('solvation["buffer"]', solvent_buffer, debug=True)
                self.print_term("solute_buffer      ", solute_buffer, debug=True)
                self.print_term("lipid_buffer       ", lipid_buffer, debug=True)
                self.print_term("protein_buffer     ", protein_buffer, debug=True)
                
                ### Creates a 3D matrix indicating the points in space that are allowed to have solvent
                ### Starts out being completely filled with ones indicating allowed space
                grid_bool_matrix = np.ones((xpoints, ypoints, zpoints))
                
                ### ### Marks coordinates as "occupied" by setting boolean to False
                ### Marks points located in hydrophobic volume
                if len(self.MEMBRANES) != 0:
                    for memb_key, memb_dict in self.MEMBRANES.items():
                        for leaflet_key, leaflet in memb_dict["leaflets"].items():
                            lcz = leaflet["center"][2]
                            ### Z-related stuff is always the same for all subleaflets
                            lhydrophob = leaflet["lipid_dimensions"]["zUnderLipCen"] # Height of hydrophobic volume
                            if leaflet["HG_direction"] == "up":
                                zminbuffer = solvent_buffer
                                zmaxbuffer = lhydrophob + solvent_buffer
                            if leaflet["HG_direction"] == "down":
                                zminbuffer = lhydrophob + solvent_buffer
                                zmaxbuffer = solvent_buffer
                            zmin, zmax = coord_to_indices(lcz, zlen, cz, zreal_gridres, zpoints, zminbuffer, zmaxbuffer)
                            
                            ### Advanced (and slow) solvent-in-membrane prevention method
                            ### Creates a series of points in 2D and marks matrix coordinates one point at a time
                            ### Allows hole to be solvated
                            if leaflet["solvate_empty"] and any([leaflet["remove_union"], leaflet["require_union"]]):
                                ### The actual leaflet
                                for (slxi, slyi, sli), subleaflet in leaflet["subleaflets"].items():
                                    xmin = subleaflet["xmin"] + solvent_buffer/2
                                    xmax = subleaflet["xmax"] - solvent_buffer/2
                                    ymin = subleaflet["ymin"] + solvent_buffer/2
                                    ymax = subleaflet["ymax"] - solvent_buffer/2
                                    xvals = np.arange(xmin, xmax+solvent_buffer/2, solvent_buffer/2)
                                    yvals = np.arange(ymin, ymax+solvent_buffer/2, solvent_buffer/2)

                                    for lcx in xvals:
                                        for lcy in yvals:
                                            point = (lcx, lcy)
                                            if subleaflet["holed_bbox"].contains(shapely.Point(point)):
                                                xmin, xmax = coord_to_indices(lcx, xlen, cx, xreal_gridres, xpoints, solvent_buffer, solvent_buffer)
                                                ymin, ymax = coord_to_indices(lcy, ylen, cy, yreal_gridres, ypoints, solvent_buffer, solvent_buffer)
                                                grid_bool_matrix[xmin:xmax, ymin:ymax, zmin:zmax] = 0
                                
                                ### Protein inserted into the leaflet
                                if leaflet["prot_points"]:
                                    for lcx, lcy in leaflet["prot_points"]:
                                        xmin, xmax = coord_to_indices(lcx, xlen, cx, xreal_gridres, xpoints, solvent_buffer*2, solvent_buffer*2)
                                        ymin, ymax = coord_to_indices(lcy, ylen, cy, yreal_gridres, ypoints, solvent_buffer*2, solvent_buffer*2)
                                        grid_bool_matrix[xmin:xmax, ymin:ymax, zmin:zmax] = 0
                            
                            ### Simple (and fast) solvent-in-membrane prevention method
                            ### Simply marks all matrix coordinates overlapping with the membrane
                            ### Prevents holes from being solvated
                            else:
                                lcx, lcy, lcz = leaflet["center"] # Center of leaflet on given axis
                                llx, lly = leaflet["x"], leaflet["y"] # Length of leaflet in given axis
                                xmin, xmax = coord_to_indices(lcx, xlen, cx, xreal_gridres, xpoints, llx/2, llx/2)
                                ymin, ymax = coord_to_indices(lcy, ylen, cy, yreal_gridres, ypoints, lly/2, lly/2)
                                grid_bool_matrix[xmin:xmax, ymin:ymax, zmin:zmax] = 0

                ### Prior solvent beads
                for xpos, ypos, zpos in solv_beads_for_cell_checker:
                    ### Checks if any non-solvent bead is within the solvent cell
                    xmin, xmax = coord_to_indices(xpos, xlen, cx, xreal_gridres, xpoints, solute_buffer, solute_buffer)
                    ymin, ymax = coord_to_indices(ypos, ylen, cy, yreal_gridres, ypoints, solute_buffer, solute_buffer)
                    zmin, zmax = coord_to_indices(zpos, zlen, cz, zreal_gridres, zpoints, solute_buffer, solute_buffer)
                    grid_bool_matrix[xmin:xmax, ymin:ymax, zmin:zmax] = 0

                ### Lipid beads
                for xpos, ypos, zpos in lipid_beads_for_cell_checker:
                    ### Checks if any non-solvent bead is within the solvent cell
                    xmin, xmax = coord_to_indices(xpos, xlen, cx, xreal_gridres, xpoints, lipid_buffer, lipid_buffer)
                    ymin, ymax = coord_to_indices(ypos, ylen, cy, yreal_gridres, ypoints, lipid_buffer, lipid_buffer)
                    zmin, zmax = coord_to_indices(zpos, zlen, cz, zreal_gridres, zpoints, lipid_buffer, lipid_buffer)
                    grid_bool_matrix[xmin:xmax, ymin:ymax, zmin:zmax] = 0

                ### Protein beads
                for xpos, ypos, zpos in prot_beads_for_cell_checker:
                    ### Checks if any non-solvent bead is within the solvent cell
                    xmin, xmax = coord_to_indices(xpos, xlen, cx, xreal_gridres, xpoints, protein_buffer, protein_buffer)
                    ymin, ymax = coord_to_indices(ypos, ylen, cy, yreal_gridres, ypoints, protein_buffer, protein_buffer)
                    zmin, zmax = coord_to_indices(zpos, zlen, cz, zreal_gridres, zpoints, protein_buffer, protein_buffer)
                    grid_bool_matrix[xmin:xmax, ymin:ymax, zmin:zmax] = 0


                ### Gets all the indices that are not occupied
                free_indices = grid_bool_matrix.nonzero()
                free_xs, free_ys, free_zs = free_indices
                free_xs_coords = xcoords[free_xs]
                free_ys_coords = ycoords[free_ys]
                free_zs_coords = zcoords[free_zs]
                
                solv_grid_3D = list(zip(free_xs_coords, free_ys_coords, free_zs_coords))

                self.print_term("Final number of 3D grid points available for solvent placement:", len(solv_grid_3D), spaces=2, verbose=4)

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
                
                ### random.sample extracts k random elements from the list without duplicates
                self.print_term("Inserting", len(collected_solvent), "solvent molecules into random grid points:", spaces=2, verbose=4)
                random_grid_points = random.sample(solv_grid_3D, k = len(collected_solvent))
                
                grid_solvated = []
                generated_spots = []
                for counter, ((sname, stype, ssize), (gx, gy, gz)) in enumerate(zip(collected_solvent, random_grid_points), 1):
                    '''
                    Finds a 3D-grid point for the solvent and places it there.
                    '''
                    if counter == 1 or counter % 25000 == 0:
                        self.print_term("Currently at solvent number:", counter, spaces=3, verbose=4)
                    sdata = solvation[stype][sname]
                    sinfo = solvation[stype][sname].get_res_beads_info()
                    
                    sbeads        = [bead.bead for bead in sinfo]
                    sresnames     = [bead.resname for bead in sinfo]
                    sbeadcharges  = [bead.charge for bead in sinfo]
                    scharge       = sdata.get_mol_charge()
                    sx, sy, sz    = sdata.get_coords("xyz")
                    rx, ry, rz    = random.uniform(0, 360), random.uniform(0, 360), random.uniform(0, 360)
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
                        "order":        order,
                        "name":         sname,
                        "resnames":     sresnames,
                        "type":         stype,
                        "beads":        sbeads,
                        "bead_charges": sbeadcharges,
                        "charge":       scharge,
                        "coords":       list(zip(sx, sy, sz)),
                    })
                self.print_term("Currently at solvent number:", counter, spaces=3, verbose=4)
                
                ### Orders the solvent key:value pairs to ensure topology-structure file match
                solvation["grid"] = sorted(grid_solvated, key=lambda g: (g["order"], g["name"]))

#                 ### Remembers the beads for subsequent solvation commands
#                 solv_beads_for_cell_checker.extend([(x, y, z) for vals in solvation["grid"] for x, y, z in vals["coords"]])

                solvation["solv_count"] = {}
                for grid_point_3D in solvation["grid"]:
                    '''
                    Counts all solvent for this specific solvation command
                    '''
                    if (grid_point_3D["name"], grid_point_3D["type"], grid_point_3D["charge"]) not in solvation["solv_count"].keys():
                        solvation["solv_count"][(grid_point_3D["name"], grid_point_3D["type"], grid_point_3D["charge"])] = 1
                    else:
                        solvation["solv_count"][(grid_point_3D["name"], grid_point_3D["type"], grid_point_3D["charge"])] += 1
                
                
                if solvation_i == 0:
                    SYS_solv_count = {}
                    SYS_solv_volume = solv_volume
                else:
                    SYS_solv_volume += solv_volume
                
                CMD_solv_count = {}
                for (sname, stype, scharge), count in solvation["solv_count"].items():
                    '''
                    Counts all solvent in the solvation command
                    '''
                    if (sname, stype, scharge) not in CMD_solv_count.keys():
                        CMD_solv_count[(sname, stype, scharge)] = count
                    else:
                        CMD_solv_count[(sname, stype, scharge)] += count
                    
                    if (sname, stype, scharge) not in SYS_solv_count.keys():
                        SYS_solv_count[(sname, stype, scharge)] = count
                    else:
                        SYS_solv_count[(sname, stype, scharge)] += count

                CMD_headers = ["Name", "Totals", "Total %", "Molarity(box)", "Molarity(free)", "Molarity(solvent)", "Charge"]
                CMD_names, CMD_charges, CMD_counts = list(zip(*[(sname, scharge, scount) for (sname, stype, scharge), scount in CMD_solv_count.items()]))
                CMD_tot = sum(CMD_counts)
                CMD_percentages = [round(count / CMD_tot * 100, 3) for count in CMD_counts]
                CMD_box_molarity = [round(scount / (N_A * box_volume), 3) for scount in CMD_counts]
                CMD_free_molarity = [round(scount / (N_A * box_free_volume), 3) for scount in CMD_counts]
                if solv_volume > 0:
                    CMD_solv_molarity = [round(scount / (N_A * solv_volume), 3) for scount in CMD_counts]
                else:
                    CMD_solv_molarity = ["NaN" for scount in CMD_counts]
                CMD_printer = [tuple(CMD_headers)] + list(zip(CMD_names, CMD_counts, CMD_percentages, CMD_box_molarity, CMD_free_molarity, CMD_solv_molarity, CMD_charges))
                CMD_max_lengths = [len(str(max(data, key = lambda d: len(str(d))))) for data in zip(*CMD_printer)]
                if self.extra_info:
                    '''
                    Prints command-specific solvent information to the terminal
                    '''
                    self.print_term("", verbose=4)
                    self.print_term("Solvent data for the the specific solvation command", spaces=1, verbose=4)
                    for string in CMD_printer:
                        for_printer = []
                        for si, substring in enumerate(string):
                            for_printer.append('{0: <{L}}'.format(substring, L = CMD_max_lengths[si]))
                            for_printer.append(":")
                        self.print_term(
                            *for_printer,
                            spaces=2,
                            verbose=4,
                        )

            if self.extra_info:
                '''
                Prints system-wide solvent information to the terminal
                '''
                self.print_term("", verbose=3)
                SYS_headers = ["Name", "Totals", "Total %", "Molarity(box)", "Molarity(free)", "Molarity(solvent)", "Charge"]
                SYS_names, SYS_charges, SYS_counts = list(zip(*[(sname, scharge, scount) for (sname, stype, scharge), scount in SYS_solv_count.items()]))
                SYS_tot = sum(SYS_counts)
                SYS_percentages = [round(count / SYS_tot * 100, 3) for count in SYS_counts]
                ### Lipid volume
                lipid_beads_for_cell_checker, solvent_box_lipids_charge = self.get_lipid_volume(
                    -self.pbcx/2, self.pbcx/2, -self.pbcy/2, self.pbcy/2, -self.pbcz/2, self.pbcz/2
                )
                SYS_leafs_volume = bead_volume * len(lipid_beads_for_cell_checker)

                ### Protein volume
                prot_beads_for_cell_checker, solvent_box_proteins_charge = self.get_protein_volume(
                    -self.pbcx/2, self.pbcx/2, -self.pbcy/2, self.pbcy/2, -self.pbcz/2, self.pbcz/2
                )
                SYS_prots_volume = bead_volume * len(prot_beads_for_cell_checker)
                SYS_box_volume = (self.pbcx * self.pbcy * self.pbcz) * 10**-27
                SYS_free_volume = SYS_box_volume - leafs_volume - prots_volume
                
                SYS_box_molarity = [round(scount / (N_A * SYS_box_volume), 3) for scount in SYS_counts]
                SYS_free_molarity = [round(scount / (N_A * SYS_free_volume), 3) for scount in SYS_counts]
                if SYS_solv_volume > 0:
                    SYS_solv_molarity = [round(scount / (N_A * SYS_solv_volume), 3) for scount in SYS_counts]
                else:
                    SYS_solv_molarity = ["NaN" for scount in SYS_counts]
                SYS_printer = [tuple(SYS_headers)] + list(zip(SYS_names, SYS_counts, SYS_percentages, SYS_box_molarity, SYS_free_molarity, SYS_solv_molarity, SYS_charges))
                SYS_max_lengths = [len(str(max(data, key = lambda d: len(str(d))))) for data in zip(*SYS_printer)]
                self.print_term("Solvent data for whole system", verbose=3)
                for string in SYS_printer:
                    for_printer = []
                    for si, substring in enumerate(string):
                        for_printer.append('{0: <{L}}'.format(substring, L = SYS_max_lengths[si]))
                        for_printer.append(":")
                    self.print_term(
                        *for_printer,
                        spaces=2,
                        verbose=3,
                    )
            
            solvation_toc = time.time()
            solvation_time = round(solvation_toc - solvation_tic, 3)
            self.print_term("------------------------------ SOLVATION COMPLETE", "(time spent: "+str(solvation_time)+" [s])", "\n", spaces=0, verbose=1)
    
    ###############
    ### PICKLER ###
    ###############
    def pickler(self):
        if self.PICKLE_cmd:
            self.print_term("-------------------", verbose=1)
            self.print_term("Pickling data into:", self.PICKLE_cmd, verbose=1)
            self.print_term("-------------------", "\n", verbose=1)
            
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

### Custom Action classes for to check if arguments have been given.
given_arguments = set()

class IsStored_ActionStore(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        given_arguments.add(self.dest)
        setattr(namespace, self.dest + '_set', True)
        setattr(namespace, self.dest, values)

class IsStored_ActionAppend(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        given_arguments.add(self.dest)
        setattr(namespace, self.dest + '_set', True)
        items = getattr(namespace, self.dest, None)
        items = argparse._copy_items(items)
        items.append(values)
        setattr(namespace, self.dest, items)

class IsStored_ActionExtend(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        given_arguments.add(self.dest)
        setattr(namespace, self.dest + '_set', True)
        items = getattr(namespace, self.dest, None)
        items = argparse._copy_items(items)
        items.extend(values)
        setattr(namespace, self.dest, items)

### Does not work
# class IsStored_ActionCount(argparse.Action):
#     def __call__(self, parser, namespace, values, option_string=None):
#         given_arguments.add(self.dest)
#         setattr(namespace, self.dest + '_set', True)
#         items = getattr(namespace, self.dest, None)
#         items = argparse._copy_items(items)
#         items += 1
#         setattr(namespace, self.dest, items)

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawTextHelpFormatter,
    add_help = False,
)

parser.add_argument("--help", "-h", dest = "help", action=IsStored_ActionStore)

#######################
### SYSTEM CREATION ###
#######################
### Leaflet commands
parser.add_argument("--membrane",  "-membrane",  "-memb", dest = "membrane_cmds", action=IsStored_ActionAppend, type=str, default = [], nargs="+")

### Protein commands
parser.add_argument("--protein",   "-protein",   "-prot", dest = "protein_cmds", action=IsStored_ActionAppend, type=str, default = [], nargs="+")

### Solvent commands
parser.add_argument("--solvation", "-solvation", "-solv", dest = "solvation_cmds", action=IsStored_ActionAppend, type=str, default = [], nargs="+")

### Solvent commands
parser.add_argument("--flooding",  "-flooding",  "-flood", dest = "flooding_cmds", action=IsStored_ActionAppend, type=str, default = [], nargs="+")

###############################
### SPECIAL SYSTEM CREATION ###
###############################
### Solvent commands
parser.add_argument("--stacked_membranes", "-stack_memb", "-stacked_membranes", dest = "stacked_membranes_cmds", action=IsStored_ActionAppend, type=str, default = [], nargs="+")

############
### MISC ###
############
### Topology commands
parser.add_argument("--itp_input", "-itp_input", "-itp_in", dest = "itp_input_cmds", action=IsStored_ActionAppend, type=str, default = [], nargs="+")

### Import commands
parser.add_argument("--solute_input", "-solute_in", "-sol_in", dest = "solute_input_cmds", action=IsStored_ActionAppend, type=str, default = [], nargs="+")

### Plotting command
parser.add_argument("--plot_grid", "-plot_grid", "-plot", dest = "plot_grid_cmd", action=IsStored_ActionStore)

### Pickle commands
parser.add_argument("--pickle", "-pickle", dest = "pickle_cmd", action=IsStored_ActionStore)

### Whether to backup files i they would be overwritten
parser.add_argument("--backup", "-backup", dest = "backup_cmd", action=IsStored_ActionStore)

### Random seed
parser.add_argument("--randseed", "-randseed", "-rand", dest = "randseed_cmd", action=IsStored_ActionStore)

### System parameters
parser.add_argument("--sys_params",   "-sys_params",   "-sysp", dest = "sys_params",   action=IsStored_ActionStore)
parser.add_argument("--prot_params",  "-prot_params",  "-pp",   dest = "prot_params",  action=IsStored_ActionStore)
parser.add_argument("--lipid_params", "-lipid_params", "-lp",   dest = "lipid_params", action=IsStored_ActionStore)
parser.add_argument("--solv_params",  "-solv_params",  "-sp",   dest = "solv_params",  action=IsStored_ActionStore)

### System name
parser.add_argument("--system_name", "-system_name", "-sn", dest = "system_name", action=IsStored_ActionStore)

#########################
### BOX SIZE AND TYPE ###
#########################
### pbc box size [nm]
parser.add_argument("--box", "-box", dest = "pbc_box", action=IsStored_ActionExtend, type=str, default = [], nargs="+")
parser.add_argument("--pbc", "-pbc", dest = "pbc_box", action=IsStored_ActionExtend, type=str, default = [], nargs="+")

### x/y/z size of box [nm]
parser.add_argument("--x", "-x", dest = "pbcx", type=str, action=IsStored_ActionStore)
parser.add_argument("--y", "-y", dest = "pbcy", type=str, action=IsStored_ActionStore)
parser.add_argument("--z", "-z", dest = "pbcz", type=str, action=IsStored_ActionStore)

### pbc box type
parser.add_argument("--box_type", "-box_type", dest = "pbc_box_type", type=str, action=IsStored_ActionStore)
parser.add_argument("--pbc_type", "-pbc_type", dest = "pbc_box_type", type=str, action=IsStored_ActionStore)

####################
### OUTPUT FILES ###
####################
### Output pdb/gro/top/log file
parser.add_argument("--out_all", "-out_all", "-o_all", dest = "out_all_file_name", action=IsStored_ActionStore)

### Output pdb/gro file
parser.add_argument("--out_sys", "-out_sys", "-o_sys", dest = "out_sys_file_name", action=IsStored_ActionStore)
### Output pdb file
parser.add_argument("--out_pdb", "-out_pdb", "-o_pdb", dest = "out_pdb_file_name", action=IsStored_ActionStore)
### Output gro file
parser.add_argument("--out_gro", "-out_gro", "-o_gro", dest = "out_gro_file_name", action=IsStored_ActionStore)

### Output topology file
parser.add_argument("--out_top", "-out_top", "-t",     dest = "out_top_file_name", action=IsStored_ActionStore)

### Log file
parser.add_argument("--out_log", "-out_log", "-log",   dest = "out_log_file_name", action=IsStored_ActionStore)

################
### PRINTING ###
################
### Prints
parser.add_argument("--print_quiet",    "-quiet",    dest = "quiet",    action=IsStored_ActionStore)
parser.add_argument("--print_debug",    "-debug",    dest = "debug",    action=IsStored_ActionStore)
parser.add_argument("--print_extra",    "-extra",    dest = "extra",    action=IsStored_ActionStore)
parser.add_argument("--print_warnings", "-warn",     dest = "warnings", action=IsStored_ActionStore)
parser.add_argument("--verbose",        "-verbose", "-v",  dest = "verbose",  action=IsStored_ActionStore)

### ### Parser for handling '-f' when importing module to Jupyter
parser.add_argument("-fff", "-f", dest = "debug_flag_for_ipython")
### unknown variable includes the weird -f flag that jupyter puts in so no need for added argument.
### Argument still needed otherwise jupyter will throw the following error:
### "ipykernel_launcher.py: error: ambiguous option: -f could match -flood, -flooding"

args, unknown = parser.parse_known_args()

##############################
### HELP FOR THOSE IN NEED ###
##############################
if "help" in given_arguments or len(sys.argv) == 1:
    parser.print_helper()
    sys.exit()

parser_kwargs = {}

parse_membrane_cmds          = [" ".join(i) for i in args.membrane_cmds]
parse_protein_cmds           = [" ".join(i) for i in args.protein_cmds]
parse_solvation_cmds         = [" ".join(i) for i in args.solvation_cmds]
parse_flooding_cmds          = [" ".join(i) for i in args.flooding_cmds]
parse_stacked_membranes_cmds = [" ".join(i) for i in args.stacked_membranes_cmds]
parse_itp_input_cmds         = [" ".join(i) for i in args.itp_input_cmds]
parse_solute_input_cmds      = [" ".join(i) for i in args.solute_input_cmds]

for CGSB_cmd, parse, arg_name in [
    ("membrane",          parse_membrane_cmds,          "membrane_cmds"),
    ("protein",           parse_protein_cmds,           "protein_cmds"),
    ("solvation",         parse_solvation_cmds,         "solvation_cmds"),
    ("flooding",          parse_flooding_cmds,          "flooding_cmds"),
    ("stacked_membranes", parse_stacked_membranes_cmds, "stacked_membranes_cmds"),
    ("itp_input",         parse_itp_input_cmds,         "itp_input_cmds"),
    ("solute_input",      parse_solute_input_cmds,      "solute_input_cmds"),
    
    ("plot_grid", args.plot_grid_cmd, "plot_grid_cmd"),
    ("pickle",    args.pickle_cmd,    "pickle_cmd"),
    ("backup",    args.backup_cmd,    "backup_cmd"),
    ("randseed",  args.randseed_cmd,  "randseed_cmd"),
    
    ("sys_params",   args.sys_params,   "sys_params"),
    ("prot_params",  args.prot_params,  "prot_params"),
    ("lipid_params", args.lipid_params, "lipid_params"),
    ("solv_params",  args.solv_params,  "solv_params"),
    
    ("box",      args.pbc_box,      "pbc_box"),
    ("x",        args.pbcx,         "pbcx"),
    ("y",        args.pbcy,         "pbcy"),
    ("z",        args.pbcz,         "pbcz"),
    ("box_type", args.pbc_box_type, "pbc_box_type"),
    
    ("out_all", args.out_all_file_name, "out_all_file_name"),
    ("out_sys", args.out_sys_file_name, "out_sys_file_name"),
    ("out_pdb", args.out_pdb_file_name, "out_pdb_file_name"),
    ("out_gro", args.out_gro_file_name, "out_gro_file_name"),
    ("out_top", args.out_top_file_name, "out_top_file_name"),
    ("out_log", args.out_log_file_name, "out_log_file_name"),
    
    ("sn", args.system_name, "system_name"),
    
    ("quiet",    args.quiet,        "quiet"),
    ("debug",    args.debug,        "debug"),
    ("extra",    args.extra,        "extra"),
    ("warn",     args.warnings,     "warnings"),
    ("verbose",  args.verbose,      "verbose"),
]:
    if arg_name in given_arguments:
        parser_kwargs[CGSB_cmd] = parse

if parser_kwargs:
    if ("verbose" in parser_kwargs and int(parser_kwargs["verbose"]) > 0) or "verbose" not in parser_kwargs:
        print("Time spent importing packages:", import_time)
    CGSB(
        run = True,
        terminal_run_kwargs = parser_kwargs,
    )

#####################################################################
########################## YOU HAVE PARSED ##########################
#####################################################################
