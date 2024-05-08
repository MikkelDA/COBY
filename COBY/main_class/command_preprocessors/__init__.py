'''
DOCSTRING TO BE WRITTEN
'''

from COBY.main_class.command_preprocessors.prot_preprocessor import *
from COBY.main_class.command_preprocessors.memb_preprocessor import *
from COBY.main_class.command_preprocessors.solv_preprocessor import *

class command_preprocessors(
    prot_preprocessor,
    memb_preprocessor,
    solv_preprocessor,
):
    def __init__(self):
        pass

