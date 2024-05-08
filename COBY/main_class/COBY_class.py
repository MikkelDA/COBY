import time

import ast
import math
import numpy as np
import random
import copy

import importlib

### Placeholders in case none are defined in "COBY.molecule_definitions.__init__"
lipid_scaffolds = {}
lipid_defs      = {}
solvent_defs    = {}
ion_defs        = {}
prot_defs       = {}

from COBY.molecule_definitions.__init__ import *
from COBY.structure_classes.__init__ import *
from COBY.general_functions.__init__ import *


from COBY.main_class.structure_file_handlers.__init__ import *
from COBY.main_class.topology_handlers.__init__ import *
from COBY.main_class.general_tools.__init__ import *
from COBY.main_class.definition_preprocessors.__init__ import *
from COBY.main_class.command_preprocessors.__init__ import *
from COBY.main_class.special_preprocessors.__init__ import *

from COBY.main_class.protein_inserter.__init__ import *
from COBY.main_class.polygon_makers.__init__ import *
from COBY.main_class.lipid_calculator.__init__ import *
from COBY.main_class.planar_grid_maker_and_optimizer.__init__ import *
from COBY.main_class.lipid_inserter.__init__ import *
from COBY.main_class.solvater.__init__ import *

class COBY(
    structure_file_handlers,
    topology_handlers,
    general_tools,
    definition_preprocessors,
    command_preprocessors,
    special_preprocessors,

    protein_inserter,
    polygon_makers,
    lipid_calculator,
    planar_grid_maker_and_optimizer,
    lipid_inserter,
    solvater,
):
    
    def __init__(self, run = True, terminal_run_kwargs = False, **kwargs):
        self.COBY_run_tic = time.time()

        self.RUN = run
        
        self.debug_prints          = False
        self.extra_info            = True
        self.subleaflet_extra_info = True
        self.warnings              = True
        self.quiet                 = False
        self.verbose               = 1
        
        if terminal_run_kwargs:
            kwargs.update(terminal_run_kwargs)
        
        
        self.PROTEINS      = {}
        self.PROTEINS_cmds = []
        
        self.MEMBRANES      = {}
        self.MEMBRANES_cmds = []
        
        self.SOLVATIONS      = {}
        self.SOLVATIONS_cmds = []
        
        ### Floodings are just solvations with some different default settings
        self.FLOODINGS_cmds = []
        
        ### Stakced membranes are special combinations of membranes and solvations
        self.STACKED_MEMBRANES_cmds = []
        
        # self.SOLUTE_IMPORT_cmds = []
        # self.LIPID_IMPORT_cmds = []
        self.MOLECULE_IMPORT_cmds = []
        
        self.PLOT_cmd        = []
        self.plot_data       = {}
        self.plots_requested = False
        
        self.plot_grid = False
        
        self.PICKLE_cmd = False
        
        self.randseed = round(time.time())
        
        self.sys_params      = "default"
        self.prot_params     = False # Only used if specifically set
        self.lipid_params    = False # Only used if specifically set
        self.solv_params     = False # Only used if specifically set
        self.topology_params = False # Only used if specifically set
        
        self.itp_defs = {
            "atomtypes":       {},
            "bondtypes":       {},
            "pairtypes":       {},
            "angletypes":      {},
            "dihedraltypes":   {},
            "constrainttypes": {},
        }
        self.itp_moleculetypes      = {}
        self.ITP_INPUT_cmds         = []
        self.TOP_include_statements = []

        self.system_charge = 0
        self.system_name = "PLACEHOLDER_TITLE"
        
#         self.output_system_file_name     = "output"
        self.output_system_pdb_file_name = False
        self.output_system_gro_file_name = False
        self.output_topol_file_name      = "topol.top"
        
        self.LOG_FILE                     = []
        self.output_log_file_name         = False
        self.terminalupdate_string_length = 80

        ### Adds given commands to log file
        self.LOG_FILE.append("The following COBY arguments will be processed:" + "\n")
        self.LOG_FILE.append("COBY(" + "\n")
        for key, val in kwargs.copy().items():
            if type(val) == str:
                val = "\"" + val + "\""
                self.LOG_FILE.append("    " + str(key) + " = " + str(val) + "," + "\n")
            elif type(val) in [list, tuple]:
                self.LOG_FILE.append("    " + str(key) + " = " + "[" + "\n")
                for subval in val:
                    if type(subval) == str:
                        subval = "\"" + subval + "\""
                    self.LOG_FILE.append("    " + "    " + str(subval) + "," + "\n")
                self.LOG_FILE.append("    " + "]," + "\n")
        self.LOG_FILE.append(")" + "\n")
            
        self.pbc_set  = []
        self.pbc_type = "rectangular"
        self.backup   = True
        self.pickle   = False
        
        self.gro_unitcell = []
        self.pdb_unitcell = []

        self.pbcx = 0
        self.pbcy = 0
        self.pbcz = 0
        
        try:
            self.lipid_scaffolds = copy.deepcopy(lipid_scaffolds)
            if len(self.lipid_scaffolds) == 0:
                self.print_term("WARNING: No lipid scaffolds found", warn=True)
        except:
            self.print_term("WARNING: No lipid scaffolds found in 'lipid_scaffolds'", warn=True)
            self.lipid_scaffolds = {}
        
        try:
            self.lipid_defs = copy.deepcopy(lipid_defs)
            if len(self.lipid_defs) == 0:
                self.print_term("WARNING: No lipid definitions found", warn=True)
        except:
            self.print_term("WARNING: No lipid definitions found in 'lipid_defs'", warn=True)
            self.lipid_defs = {}
        
        try:
            self.solvent_defs = copy.deepcopy(solvent_defs)
            if len(self.solvent_defs) == 0:
                self.print_term("WARNING: No solvent definitions found", warn=True)
        except:
            self.print_term("WARNING: No solvent definitions found in 'solvent_defs'", warn=True)
            self.solvent_defs = {}
        
        try:
            self.ion_defs = copy.deepcopy(ion_defs)
            if len(self.ion_defs) == 0:
                self.print_term("WARNING: No ion definitions found", warn=True)
        except:
            self.print_term("WARNING: No ions definitions found in 'ion_defs'", warn=True)
            self.ion_defs = {}
        
        try:
            self.prot_defs = copy.deepcopy(prot_defs)
            if len(self.prot_defs) == 0:
                self.print_term("WARNING: No protein charge definitions found", warn=True)
        except:
            self.print_term("WARNING: No protein definitions found in 'prot_defs'", warn=True)
            self.prot_defs = {}
        
        self.lipid_defs_imported   = {}
        self.solvent_defs_imported = {}
        self.ion_defs_imported     = {}


        self.lipid_dict   = {}
        self.solvent_dict = {}
        self.ion_dict     = {}
        self.prot_dict    = {}

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
        momentary_x   = 0
        momentary_y   = 0
        momentary_z   = 0
        sm_z_value    = False # stacked membranze z value
        for key, cmd in kwargs.items():
            ### Molecule definitions and scaffolds
            if key.startswith("molecule_definition"):
                if type(cmd) not in [list, tuple]:
                    cmd = [cmd]
                for subcmd in cmd:
                    assert subcmd.endswith(".py"), "Molecule definitions / lipid scaffolds file must be python file: '" + subcmd + "'"
                    
                    spec = importlib.util.spec_from_file_location(
                        name="defs_module",  # note that ".test" is not a valid module name
                        location=subcmd,
                    )
                    defs_module = importlib.util.module_from_spec(spec)
                    spec.loader.exec_module(defs_module)
                    
                    if hasattr(defs_module, "lipid_scaffolds"):
                        self.lipid_scaffolds.update(defs_module.lipid_scaffolds)
                    if hasattr(defs_module, "lipid_defs"):
                        self.lipid_defs.update(defs_module.lipid_defs)
                    if hasattr(defs_module, "solvent_defs"):
                        self.solvent_defs.update(defs_module.solvent_defs)
                    if hasattr(defs_module, "ion_defs"):
                        self.ion_defs.update(defs_module.ion_defs)
                    if hasattr(defs_module, "prot_defs"):
                        self.prot_defs.update(defs_module.prot_defs)

            ### General system inputs
            elif any(key.startswith(i) for i in ["protein", "prot"]):
                if type(cmd) not in [list, tuple]:
                    cmd = [cmd]
                for subcmd in cmd:
                    self.PROTEINS_cmds.extend([subcmd])
                
            elif any(key.startswith(i) for i in ["membrane", "memb"]):
                if type(cmd) not in [list, tuple]:
                    cmd = [cmd]
                for subcmd in cmd:
                    self.MEMBRANES_cmds.extend([subcmd])
                
            elif any(key.startswith(i) for i in ["solvation", "solv"]):
                if type(cmd) not in [list, tuple]:
                    cmd = [cmd]
                for subcmd in cmd:
                    self.SOLVATIONS_cmds.extend([subcmd])
                
            elif any(key.startswith(i) for i in ["flooding", "flood"]):
                if type(cmd) not in [list, tuple]:
                    cmd = [cmd]
                for subcmd in cmd:
                    subcmd = self.flooding_to_solvation_converter(subcmd)
                    self.FLOODINGS_cmds.extend([subcmd])
            
            elif any(key.startswith(i) for i in ["stacked_membranes", "stack_memb"]):
                if type(cmd) not in [list, tuple]:
                    cmd = [cmd]
                for subcmd in cmd:
                    self.STACKED_MEMBRANES_cmds.extend([subcmd])
            
            ### Box type
            elif key in ["pbc_type", "box_type"]:
                if cmd in ["rectangular"]:
                    self.pbc_type = "rectangular"

                elif cmd in ["hexagonal", "hexagonal_prism"]:
                    ### Hexagonal prism (hexagonal XY-plane, and straight Z-axis)
                    ### Y-values ignored as only one side-length is allowed. X-value used instead.
                    ### Z-value used for height
                    ### Not actually a hexagon, but instead a parallelepiped constituting a third of the hexagon
                    self.pbc_type = "hexagonal"

                elif cmd in ["skewed_hexagonal", "skewed_hexagonal_prism"]:
                    ### Skewed hexagonal prism /  (hexagonal XY-plane, and angled Z-axis)
                    ### Y-values ignored as only one side-length is allowed. X-value used instead.
                    ### Z-value used for height
                    ### Not actually a hexagon, but instead a parallelepiped constituting a third of the hexagon
                    self.pbc_type = "skewed_hexagonal"

                elif cmd in ["dodecahedron", "rhombic_dodecahedron"]:
                    ### Rhombic dodecahedron (hexagonal XY-plane, Z-angled)
                    ### Y-values ignored as only one side-length is allowed. X-value used instead.
                    ### Z-value unused as it is calculated from x-value due to angles.
                    ### Not actually a rhombic dodecahedron, but instead a parallelepiped constituting a portion of dodecahedron
                    self.pbc_type = "dodecahedron"

                else:
                    assert False, "Incorrect pbc type given: " + str(key, cmd)
            
            elif key == "pdb_unitcell":
                assert len(cmd) in [3, 6], "number of values given to 'pdb_unitcell' must be either 3 or 6 (pdb_unitcell: {})".format(cmd)
                assert all([self.is_number(i)[0] for i in cmd]), "All values given to 'pdb_unitcell' must be numbers"
                self.pdb_unitcell = [self.get_number_from_string(i) for i in cmd]
            
            elif key == "gro_unitcell":
                assert len(cmd) in [3, 9], "Number of values given to 'gro_unitcell' must be either 3 or 9 (gro_unitcell: {})".format(cmd)
                assert all([self.is_number(i)[0] for i in cmd]), "All values given to 'gro_unitcell' must be numbers"
                self.gro_unitcell = [self.get_number_from_string(i) for i in cmd]

            ### Box size
            elif key in ["pbc", "box"]:
                ### Evaluated after loop due to dependence on box type
                for i, val in enumerate(cmd):
                    if type(val) == str:
                        isnumber, isint = self.is_number(val)
                        if isint:
                            val = int(val)
                        else:
                            val = float(val)
                    momentary_pbc.append(val)
                    
            elif key == "x":
                val = cmd
                if type(val) == str:
                    isnumber, isint = self.is_number(val)
                    if isint:
                        val = int(val)
                    else:
                        val = float(val)
                momentary_x = val
                
            elif key == "y":
                val = cmd
                if type(val) == str:
                    isnumber, isint = self.is_number(val)
                    if isint:
                        val = int(val)
                    else:
                        val = float(val)
                momentary_y = val
                
            elif key == "z":
                val = cmd
                if type(val) == str:
                    isnumber, isint = self.is_number(val)
                    if isint:
                        val = int(val)
                    else:
                        val = float(val)
                momentary_z = val
            
            ### Imports
            elif key in ["itp_input", "itp_in"]:
                if type(cmd) != list:
                    cmd = [cmd]
                for subcmd in cmd:
                    self.ITP_INPUT_cmds.extend([subcmd])
            
            elif key in ["molecule_import"]:
                ### Puts individual string inside list
                if type(cmd) not in [list, tuple]:
                    cmd = [cmd]
                ### Converts tuple to list
                if type(cmd) == tuple:
                    cmd = list(cmd)
                for subcmd in cmd:
                    self.MOLECULE_IMPORT_cmds.extend([subcmd])
            
            ### Outputs
            elif key in ["out_all", "o_all"]:
                ### Cuts the extension if present so that all files can be generated with proper extensions
                if any(cmd.lower().endswith(string) for string in [".pdb", ".gro", ".top", ".log"]):
                    cmd = cmd[:-4]
                self.output_system_pdb_file_name = cmd + ".pdb"
                self.output_system_gro_file_name = cmd + ".gro"
                self.output_topol_file_name      = cmd + ".top"
                self.output_log_file_name        = cmd + ".log"
            
            elif key in ["out_sys", "o_sys"]:
                if not any([cmd.lower().endswith(i) for i in [".pdb", ".gro"]]):
                    self.output_system_pdb_file_name = cmd + ".pdb"
                    self.output_system_gro_file_name = cmd + ".gro"
                elif cmd.lower().endswith(".pdb"):
                    self.output_system_pdb_file_name = cmd
                elif cmd.lower().endswith(".gro"):
                    self.output_system_gro_file_name = cmd
                else:
                    assert False, "Unknown file extension used for 'output_system': " + cmd
                    
            elif key in ["out_pdb", "o_pdb"]:
                if not cmd.lower().endswith(".gro"):
                    cmd = cmd + ".pdb"
                self.output_system_pdb_file_name = cmd
                    
            elif key in ["out_gro", "o_gro"]:
                if not cmd.lower().endswith(".gro"):
                    cmd = cmd + ".gro"
                self.output_system_gro_file_name = cmd
            
            elif key in ["out_top", "o_top"]:
                if not cmd.lower().endswith(".top"):
                    cmd = cmd + ".top"
                self.output_topol_file_name = cmd
                
            elif key in ["out_log", "o_log"]:
                if not cmd.lower().endswith(".log"):
                    cmd = cmd + ".log"
                self.output_log_file_name = cmd
            
            elif key in ["system_name", "sn"]:
                self.system_name = cmd
            
            elif key in ["plot_grid"]:
                if type(cmd) == str:
                    cmd = ast.literal_eval(cmd)
                self.plot_grid = True
                    
            elif key in ["rand", "randseed"]:
                cmd = self.get_number_from_string(cmd)
                assert type(cmd) in [int, float], "Value given to rand/randseed must be a number (string-numbers are allowed): " + str(cmd)
                self.randseed = round(cmd)
                    
            elif key in ["pickle"]:
                self.PICKLE_cmd = cmd
                
            elif key in ["backup"]:
                assert cmd in ["False", "True", "0", "1", False, True, 0, 1], "Value given to 'backup' must be False/True/0/1 (strings allowed): " + cmd
                if type(cmd) == str:
                    cmd = ast.literal_eval(cmd)
                self.backup = cmd

            elif key in ["params", "sys_params"]:
                self.sys_params = cmd
                
            elif key in ["prot_params"]:
                self.prot_params = cmd
                
            elif key in ["lipid_params"]:
                self.lipid_params = cmd
                
            elif key in ["solv_params"]:
                self.solv_params = cmd
            
            ### Printer settings
            elif key in ["quiet"]:
                assert cmd in ["False", "True", "0", "1", False, True, 0, 1], "Value given to 'quiet' must be False/True/0/1 (strings allowed): " + cmd
                if type(cmd) == str:
                    cmd = ast.literal_eval(cmd)
                self.quiet = cmd
                
            elif key in ["debug"]:
                assert cmd in ["False", "True", "0", "1", False, True, 0, 1], "Value given to 'debug' must be False/True/0/1 (strings allowed): " + cmd
                if type(cmd) == str:
                    cmd = ast.literal_eval(cmd)
                self.debug_prints = cmd
                
            elif key in ["extra"]:
                assert cmd in ["False", "True", "0", "1", False, True, 0, 1], "Value given to 'extra' must be False/True/0/1 (strings allowed): " + cmd
                if type(cmd) == str:
                    cmd = ast.literal_eval(cmd)
                self.extra_info = cmd
            
            elif key in ["warn"]:
                assert cmd in ["False", "True", "0", "1", False, True, 0, 1], "Value given to 'warn' must be False/True/0/1 (strings allowed): " + cmd
                if type(cmd) == str:
                    cmd = ast.literal_eval(cmd)
                self.warnings = cmd
                
            elif key == "verbose":
                number = self.get_number_from_string(cmd)
                if number is False:
                    self.verbose = len(cmd)
                else:
                    self.verbose = number
                
            ### Run the program
            elif key in ["run"]:
                if type(cmd) == str:
                    cmd = ast.literal_eval(cmd)
                self.RUN = cmd
            
            else:
                assert False, "Invalid argument given to COBY: " + str((key, cmd))
        
        ### Setting randseed
        self.print_term("\n" + "Setting random seed to:", self.randseed, verbose=1)
        random.seed(self.randseed)
        np.random.seed(self.randseed)
        
        string = " ".join(["", "PREPROCESSING DEFINITIONS", ""])
        self.print_term("{string:-^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), spaces=0, verbose=1)

        ### Definition preprocessing
        preprocessing_tic = time.time()
        self.itp_read_initiater()
        self.molecule_importer()
        self.lipid_scaffolds_preprocessor()
        self.lipid_defs_preprocessor()
        ### Preprocess ions first as they add the "default" parameter libraries of "neg_ions" and "pos_ions" to solvent defs
        self.ion_defs_preprocessor()
        self.solvent_defs_preprocessor()
        preprocessing_toc = time.time()
        preprocessing_time = round(preprocessing_toc - preprocessing_tic, 4)
        string = " ".join(["", "DEFINITIONS PREPROCESSING COMPLETE", ""])
        self.print_term("{string:-^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), spaces=0, verbose=1)
        string = " ".join(["", "(Time spent:", str(preprocessing_time), "[s])", ""])
        self.print_term("{string:^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), "\n", spaces=0, verbose=1)
        
        assert len(self.STACKED_MEMBRANES_cmds) in [0, 1], "Only a single 'stacked_membranes' argument is allowed."

        if len(self.STACKED_MEMBRANES_cmds) > 0:
            string = " ".join(["", "PRE-PREPROCESSING STACKED MEMBRANE COMMANDS", ""])
            self.print_term("{string:-^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), verbose=1)
            sm_preprocessing_tic = time.time()

            assert len(self.gro_unitcell) == 0 and len(self.pdb_unitcell) == 0, "'stacked_membranes' argument is not compatable with custom unit cells."

            if self.pbc_type == "rectangular":
                assert len(momentary_pbc) in [0, 2] or sum([bool(val) for val in [momentary_x, momentary_y, momentary_z]]) in [0, 2], "\n".join([
                    "'stacked_membrane' calculates one box length value, meaning that you should give one less value than usually to the 'box'/'pbc' argument for the 'rectangular' pbc type.",
                    "    "+"Rectangular: Give 2 values to the 'box'/'pbc' argument. Values sets the x/y-lengths of the box.",
                    "    "+"Rectangular (cubic): Do not use the 'box'/'pbc' argument.",
                ])

            if self.pbc_type == "hexagonal":
                assert len(momentary_pbc) in [0, 1] or sum([bool(val) for val in [momentary_x, momentary_y, momentary_z]]) in [0, 1], "\n".join([
                    "'stacked_membrane' calculates one box length value, meaning that you should give one less value than usually to the 'box'/'pbc' argument for the 'hexagonal' pbc type.",
                    "    "+"Hexagonal: Give 1 value to the 'box'/'pbc' argument. Value sets the x/y-lengths of the box.",
                    "    "+"Hexagonal (z-length is the same as the x/y-plane): Give 0 values to the 'box'/'pbc' argument.",
                ])

            if self.pbc_type == "skewed_hexagonal":
                assert len(momentary_pbc) == 0 or sum([bool(val) for val in [momentary_x, momentary_y, momentary_z]]) == 0, "\n".join([
                    "'stacked_membrane' calculates one box length value, meaning that you should give one less value than usually to the 'box'/'pbc' argument for the 'skewed_hexagonal' pbc type.",
                    "    "+"Skewed hexagonal: Give 0 values to the 'box'/'pbc' argument.",
                ])

            assert self.pbc_type != "dodecahedron", "\n".join([
                    "The pbc type 'dodecahedron' is not compatible with the 'stacked_membrane' argument. Consider using the 'skewed_hexagonal' pbc type instead",
                ])

            sm_z_value = self.stacked_membranes_preprocessor()

            sm_preprocessing_toc = time.time()
            sm_preprocessing_time = round(sm_preprocessing_toc - sm_preprocessing_tic, 4)
            string = " ".join(["", "STACKED MEMBRANE PRE-PREPROCESSING COMPLETE", ""])
            self.print_term("{string:-^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), spaces=0, verbose=1)
            string = " ".join(["", "(Time spent:", str(sm_preprocessing_time), "[s])", ""])
            self.print_term("{string:^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), "\n", spaces=0, verbose=1)
        
        if len(self.FLOODINGS_cmds) > 0:
            self.SOLVATIONS_cmds = self.FLOODINGS_cmds + self.SOLVATIONS_cmds
        
        ### Setting box size values to be used in PBC type settings
        if self.pdb_unitcell or self.gro_unitcell:
            assert not (self.pdb_unitcell and self.gro_unitcell), "Only one of the 'pdb_unitcell' and 'gro_unitcell' arguments may be used in a single call to COBY"
            assert len(momentary_pbc) == 0 and sum([bool(val) for val in [momentary_x, momentary_y, momentary_z]]) == 0, "Box dimensions may not be given if custom unit cells are designated."

            if self.gro_unitcell:
                if len(self.gro_unitcell) == 3:
                    self.pbcx, self.pbcy, self.pbcz = [i for i in self.gro_unitcell]
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
                elif len(self.gro_unitcell) == 9:
                    self.pbcx, self.pbcy, self.pbcz = [i*10 for i in self.gro_unitcell[:3]]
                    self.gro_box_vectors = [float(i) for i in self.gro_unitcell]
                
                ### ### Converting gro unit cell vectors to pdb unit cell lengths and angles
                ### Create vectors
                vector1 = np.array([self.gro_box_vectors[0], self.gro_box_vectors[3], self.gro_box_vectors[4]])
                vector2 = np.array([self.gro_box_vectors[5], self.gro_box_vectors[1], self.gro_box_vectors[6]])
                vector3 = np.array([self.gro_box_vectors[7], self.gro_box_vectors[8], self.gro_box_vectors[2]])
                
                assert vector1[1] == vector1[2] == vector2[2] == 0, "Gromacs only supports boxes with v1(y)=v1(z)=v2(z)=0, https://manual.gromacs.org/current/reference-manual/file-formats.html#gro"

                ### Calculate angles
                ### https://math.stackexchange.com/questions/361412/finding-the-angle-between-three-points
                ### https://en.wikipedia.org/wiki/Lattice_constant#/media/File:UnitCell.png
                alpha = math.degrees(np.arccos(np.dot(vector2, vector3) / (np.linalg.norm(vector2) * np.linalg.norm(vector3))))
                beta  = math.degrees(np.arccos(np.dot(vector1, vector3) / (np.linalg.norm(vector1) * np.linalg.norm(vector3))))
                gamma = math.degrees(np.arccos(np.dot(vector1, vector2) / (np.linalg.norm(vector1) * np.linalg.norm(vector2))))

                ### Calculate lengths
                a = np.linalg.norm(vector1)
                b = np.linalg.norm(vector2)
                c = np.linalg.norm(vector3)

                self.pdb_box_dimension = [
                    round(a*10, 3),
                    round(b*10, 3),
                    round(c*10, 3),
                    round(alpha, 2),
                    round(beta, 2),
                    round(gamma, 2),
                ]
            
            elif self.pdb_unitcell:
                if len(self.pdb_unitcell) == 3:
                    # self.pbcx, self.pbcy, self.pbcz = [i*10 for i in self.pdb_unitcell]
                    self.pdb_box_dimension = [
                        float(self.pdb_unitcell[0]*10), # Axis length:  x
                        float(self.pdb_unitcell[1]*10), # Axis length:  y
                        float(self.pdb_unitcell[2]*10), # Axis length:  z
                        float(90),        # Corner angle: alpha
                        float(90),        # Corner angle: beta
                        float(90),        # Corner angle: gamma
                    ]
                elif len(self.pdb_unitcell) == 6:
                    # self.pbcx, self.pbcy, self.pbcz = [i*10 for i in self.pdb_unitcell[:3]]
                    self.pdb_box_dimension = [float(i)*10 for i in self.pdb_unitcell[:3]]
                    self.pdb_box_dimension.extend([float(i) for i in self.pdb_unitcell[3:]])
                
                ### ### Converting pdb unit cell lengths and angles to gro unit cell vectors
                ### Get lengths
                a = self.pdb_box_dimension[0]
                b = self.pdb_box_dimension[1]
                c = self.pdb_box_dimension[2]

                ### Get angles
                alpha = self.pdb_box_dimension[3]
                beta  = self.pdb_box_dimension[4]
                gamma = self.pdb_box_dimension[5]

                ### Convert degrees to radians
                alpha_rad = math.radians(alpha)
                beta_rad  = math.radians(beta)
                gamma_rad = math.radians(gamma)

                ### First vector aligned with x-axis
                vector1 = [
                    a,
                    0,
                    0,
                ]

                ### Second vector in x/y-plane
                vector2 = [
                    b * math.cos(gamma_rad),
                    b * math.sin(gamma_rad),
                    0,
                ]

                ### Third vector exists in all dimensions
                vector3 = [
                    c * math.cos(beta_rad),
                    c * (math.cos(alpha_rad) - math.cos(beta_rad) * math.cos(gamma_rad)) / math.sin(gamma_rad),
                    c * math.sqrt(1 - math.cos(beta_rad)**2 - ((math.cos(alpha_rad) - math.cos(beta_rad) * math.cos(gamma_rad)) / math.sin(gamma_rad))**2),
                ]
                
                assert vector1[1] == vector1[2] == vector2[2] == 0, "Gromacs only supports boxes with v1(y)=v1(z)=v2(z)=0, https://manual.gromacs.org/current/reference-manual/file-formats.html#gro"

                self.gro_box_vectors = [
                    round(vector1[0]/10, 5),
                    round(vector2[1]/10, 5),
                    round(vector3[2]/10, 5),
                    round(vector1[1]/10, 5),
                    round(vector1[2]/10, 5),
                    round(vector2[0]/10, 5),
                    round(vector2[2]/10, 5),
                    round(vector3[0]/10, 5),
                    round(vector3[1]/10, 5),
                ]

                self.pbcx = vector1[0]
                self.pbcy = vector2[1]
                self.pbcz = vector3[2]

        else:
            if len(momentary_pbc) > 0:
                if sm_z_value:
                    pbc = momentary_pbc + [sm_z_value]
                else:
                    pbc = momentary_pbc
            else:
                if sm_z_value:
                    pbc = [i for i in [momentary_x, momentary_y, sm_z_value] if i != 0]
                else:
                    pbc = [i for i in [momentary_x, momentary_y, momentary_z] if i != 0]

            assert len(pbc) != 0, "No box dimensions have been given. Please do so using the 'box'/'pbc' or 'x'/'y'/'z' arguments."
            
            ### PBC type settings:
            if self.pbc_type == "rectangular":
                assert len(pbc) in [1, 3], "Either one or three values must given as box size when using the 'rectangular' 'pbc_type'. Values given: {values}".format(str(pbc))
                if len(pbc) == 1:
                    ### Cubic
                    self.pbcx = pbc[0]*10
                    self.pbcy = pbc[0]*10
                    self.pbcz = pbc[0]*10
                elif len(pbc) == 3:
                    self.pbcx, self.pbcy, self.pbcz = [i*10 for i in pbc]
                else:
                    assert False, "Incorrect pbc box dimensions for pbc type '"+self.pbc_type+"': " + str(pbc)
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
                assert len(pbc) in [1, 2], "Either one or two values must given as box size when using the 'hexagonal' 'pbc_type'. Values given: {values}".format(str(pbc))
                ### Not actually a hexagon, but instead a parallelepiped constituting a third of it
                ### Uses 'rhombic dodecahedron (xy-hexagon) "a" and "b" vectors' and 'cubic "c" vector'
                self.pbcx = pbc[0]*10
                self.pbcy = math.sqrt(3)*self.pbcx/2
                if len(pbc) == 1:
                    self.pbcz = pbc[0]*10
                elif len(pbc) == 2:
                    self.pbcz = pbc[1]*10
                else:
                    assert False, "Incorrect pbc box dimensions for pbc type '"+self.pbc_type+"': " + str(pbc)
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
            elif self.pbc_type == "skewed_hexagonal":
                assert len(pbc) == 1, "Exactly one value must given as box size when using the 'skewed_hexagonal' 'pbc_type'. Values given: {values}".format(str(pbc))
                ### Not actually a dodecahedron, but instead a parallelepiped constituting a third of it
                ### Uses 'rhombic dodecahedron (xy-hexagon) vectors'
                if len(pbc) == 1:
                    if sm_z_value:
                        ### Done when "stacked_membranes" is used as the box size is defined by the z-height
                        self.pbcz = pbc[0]*10
                        self.pbcx = self.pbcz/math.sqrt(6)*3
                        self.pbcy = math.sqrt(3)*self.pbcx/2
                    else:
                        self.pbcx = pbc[0]*10
                        self.pbcy = math.sqrt(3)*self.pbcx/2
                        self.pbcz = math.sqrt(6)*self.pbcx/3
                else:
                    assert False, "Incorrect pbc box dimensions for pbc type '"+self.pbc_type+"': " + str(pbc)
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
            
            ### Following math taken from insane.py and https://manual.gromacs.org/current/reference-manual/algorithms/periodic-boundary-conditions.html
            elif self.pbc_type == "dodecahedron":
                assert len(self.MEMBRANES_cmds) == 0, "\n".join([
                    "'dodecahedron' box type should not be used with systems containing membranes",
                    "If you are looking for 'optimal' from insane.py, then use 'skewed_hexagonal'",
                    "If you are certain you wish to use it with a membrane then input the unitcell parameters manually using the arguments 'pdb_unitcell' or 'gro_unitcell'",
                ])
                assert len(pbc) == 1, "Exactly one value must given as box size when using the 'dodecahedron' 'pbc_type'. Values given: {values}".format(str(pbc))
                ### Not actually a hexagon, but instead a parallelepiped constituting a third of it
                ### Uses 'rhombic dodecahedron (xy-hexagon) "a" and "b" vectors' and 'cubic "c" vector'
                if len(pbc) == 1:
                    if sm_z_value:
                        ### Done when "stacked_membranes" is used as the box size is defined by the z-height
                        self.pbcz = pbc[0]*10
                        self.pbcx = self.pbcz/math.sqrt(2)*2
                        self.pbcy = self.pbcz/math.sqrt(2)*2
                    else:
                        self.pbcx = pbc[0]*10
                        self.pbcy = pbc[0]*10
                        self.pbcz = math.sqrt(2)*self.pbcx/2
                else:
                    assert False, "Incorrect pbc box dimensions for pbc type '"+self.pbc_type+"': " + str(pbc)
                self.gro_box_vectors = [
                    float(self.pbcx/10),                       # vector: vax or v1x
                    float(self.pbcy/10),                       # vector: vby or v2y
                    float(self.pbcz/10),                       # vector: vcz or v3z
                    
                    float(0),                                  # vector: vay or v1y
                    float(0),                                  # vector: vaz or v1z
                    float(0),                                  # vector: vbx or v2x
                    
                    float(0),                                  # vector: vbz or v2z
                    float((self.pbcx/10)/2),                   # vector: vcx or v3x
                    float((self.pbcx/10)/2),                   # vector: vcy or v3y
                ]
                self.pdb_box_dimension = [
                    float(self.pbcx), # Axis length:  x
                    float(self.pbcx), # Axis length:  y
                    float(self.pbcx), # Axis length:  z
                    float(60),        # Corner angle: alpha
                    float(60),        # Corner angle: beta
                    float(90),        # Corner angle: gamma
                ]
            
        self.pbc_box = [self.pbcx, self.pbcy, self.pbcz]
        
        if not self.output_system_pdb_file_name and not self.output_system_gro_file_name:
            self.output_system_pdb_file_name = "output.pdb"
            self.output_system_gro_file_name = "output.gro"
    
    def run(self, kwargs):
        '''
        Runs the entire system creation process
        '''
        
        self.commands_handler(kwargs)

        assert any([len(cmd) > 0 for cmd in [self.PROTEINS_cmds, self.MEMBRANES_cmds, self.SOLVATIONS_cmds]]), (
            "Running requires at least one argument of one of the following types: 'protein', 'membrane', 'solvation' or 'flooding'"
        )
        
        string = " ".join(["", "PREPROCESSING COMMANDS", ""])
        self.print_term("{string:-^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), spaces=0, verbose=1)
        ### Command preprocessing
        preprocessing_tic = time.time()
        self.prot_preprocessor()
        self.memb_preprocessor(return_self = "self", cmds_given = False)
        self.solv_preprocessor()
        preprocessing_toc = time.time()
        preprocessing_time = round(preprocessing_toc - preprocessing_tic, 4)
        string = " ".join(["", "COMMAND PREPROCESSING COMPLETE", ""])
        self.print_term("{string:-^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), spaces=0, verbose=1)
        string = " ".join(["", "(Time spent:", str(preprocessing_time), "[s])", ""])
        self.print_term("{string:^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), "\n", spaces=0, verbose=1)

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
        COBY_run_toc  = time.time()
        COBY_run_time = round(COBY_run_toc - self.COBY_run_tic, 4)
        self.print_term("Time spent running COBY:", COBY_run_time, verbose=1)
