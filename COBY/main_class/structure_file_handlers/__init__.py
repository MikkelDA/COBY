'''
DOCSTRING TO BE WRITTEN
'''

from COBY.main_class.structure_file_handlers.pdb_reader import *
from COBY.main_class.structure_file_handlers.gro_reader import *
from COBY.main_class.structure_file_handlers.cif_reader import *

from COBY.main_class.structure_file_handlers.pdb_atom_writer import *
from COBY.main_class.structure_file_handlers.gro_atom_writer import *
from COBY.main_class.structure_file_handlers.cif_atom_writer import *

from COBY.main_class.structure_file_handlers.system_file_writer import *
from COBY.main_class.structure_file_handlers.molecule_importer import *
from COBY.main_class.structure_file_handlers.structure_importer import *

class structure_file_handlers(
    pdb_reader,
    gro_reader,
    cif_reader,

    pdb_atom_writer,
    gro_atom_writer,
    cif_atom_writer,

    system_file_writer,
    molecule_importer,
    structure_importer,
):
    def __init__(self):
        pass
