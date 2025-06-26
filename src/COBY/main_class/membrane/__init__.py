'''
DOCSTRING TO BE WRITTEN
'''

from COBY.main_class.membrane.step1_polygon_makers.__init__ import *
from COBY.main_class.membrane.step2_lipid_calculator.__init__ import *
from COBY.main_class.membrane.step3_planar_grid_maker_and_optimizer.__init__ import *
from COBY.main_class.membrane.step4_lipid_inserter.__init__ import *
from COBY.main_class.membrane.grid_plotter.__init__ import *

class membrane(
    polygon_makers,
    lipid_calculator,
    planar_grid_maker_and_optimizer,
    lipid_inserter,
    grid_plotter,
):
    def __init__(self):
        pass

