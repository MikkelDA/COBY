# CGSB - Coarse Grained System Builder
CGSB (Coarse Grained System Builder) can be used to create coarse-grained systems in Martini 3

The program is written completely from scratch only reusing the initial lipid, protein and solvent definitions (though the format has been changed) from the original python2 insane script by Tsjerk.

Note that the program is currently under development and thus may not be feature complete and may contain redundant/outcommented code.

CGSB imports the following packages:
- time
- ast
- itertools
- math
- numpy as np
- random
- copy
- scipy.spatial
- alphashape
- shapely
- shapely.plotting
- inspect
- sys
- os
- argparse
- operator
- matplotlib
- pickle

Most are likely already installed but the following might require installations

[alphashape](https://pypi.org/project/alphashape/ ): pip install alphashape

[shapely](https://pypi.org/project/shapely/): pip install shapely


