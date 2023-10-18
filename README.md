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

Installation procedure:

Version 1:

conda env create -f environment.yml

conda activate CGSB

python -m ipykernel install --user --name=CGSB

Version 2:

conda create --name CGSB python==3.9

conda activate CGSB

pip install numpy==1.21.5 scipy==1.7.3 alphashape matplotlib ipykernel ipywidgets

python -m ipykernel install --user --name=CGSB
