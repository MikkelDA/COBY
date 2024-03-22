'''
DOCSTRING TO BE WRITTEN
'''

from CGSB.main_class.topology_handlers.itp_reader import *
from CGSB.main_class.topology_handlers.topol_file_writer import *
from CGSB.main_class.topology_handlers.itp_read_initiater import *

class topology_handlers(
    itp_reader,
    topol_file_writer,
    itp_read_initiater,
):
    def __init__(self):
        pass

