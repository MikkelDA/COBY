'''
DOCSTRING TO BE WRITTEN
'''

from COBY.main_class.definition_preprocessors.lipid_scaffolds_preprocessor import *
from COBY.main_class.definition_preprocessors.lipid_defs_preprocessor import *
from COBY.main_class.definition_preprocessors.molecule_defs_checker import *
from COBY.main_class.definition_preprocessors.molecule_beads_checker import *
from COBY.main_class.definition_preprocessors.solvent_defs_preprocessor import *
from COBY.main_class.definition_preprocessors.ion_defs_preprocessor import *

class definition_preprocessors(
    lipid_scaffolds_preprocessor,
    lipid_defs_preprocessor,
    molecule_defs_checker,
    molecule_beads_checker,
    solvent_defs_preprocessor,
    ion_defs_preprocessor,
):
    def __init__(self):
        pass

