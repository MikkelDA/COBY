import time

import ast
import math
import numpy as np
import random
import copy

### Placeholders in case none are defined in "COBY.molecule_definitions.__init__" and "COBY.fragment_definitions.__init__"
lipid_scaffolds = {}
lipid_defs      = {}
solvent_defs    = {}
pos_ion_defs    = {}
neg_ion_defs    = {}
fragment_defs   = {}

from COBY.molecule_definitions.__init__ import *
from COBY.fragment_definitions.__init__ import *
from COBY.structure_classes.__init__ import *
from COBY.general_functions.__init__ import *

from COBY.main_class.structure_file_handlers.__init__ import *
from COBY.main_class.topology_handlers.__init__ import *
from COBY.main_class.general_tools.__init__ import *
from COBY.main_class.definition_preprocessors.__init__ import *
from COBY.main_class.molecule_fragment_builder.__init__ import *


class Crafter(
    structure_file_handlers,
    topology_handlers,
    general_tools,
    definition_preprocessors,
    molecule_fragment_builder,
):

    def __init__(self, run = True, terminal_run_kwargs = False, **kwargs):
        self.COBY_run_tic = time.time()
        self.PROGRAM = "Crafter"

        self.RUN = run
        
        self.debug_prints = False
        self.debug_keys   = []
        self.extra_info   = True
        self.warnings     = True
        self.quiet        = False
        self.verbose      = 1

        if terminal_run_kwargs:
            kwargs.update(terminal_run_kwargs)

        self.MOLECULE_FRAGMENT_BUILDER_cmds = []
        self.MOLECULE_cmds = []
        self.MOLECULES = []

        self.output_system_pdb_file_name = False
        self.output_system_gro_file_name = False
        
        self.LOG_FILE                     = []
        self.output_log_file_name         = False
        self.terminalupdate_string_length = 80

        ### Adds given commands to log file
        self.LOG_FILE.append("The following MoleculeBuilder arguments will be processed:" + "\n")
        self.LOG_FILE.append("MoleculeBuilder(" + "\n")
        for key, val in kwargs.copy().items():
            if type(val) in [str, int, float, bool]:
                if type(val) == str:
                    val = "\"" + val + "\""
                if type(val) in [int, float, bool]:
                    val = str(val)
                self.LOG_FILE.append("    " + str(key) + " = " + str(val) + "," + "\n")
            elif type(val) in [list, tuple]:
                self.LOG_FILE.append("    " + str(key) + " = " + "[" + "\n")
                for subval in val:
                    if type(subval) == str:
                        subval = "\"" + subval + "\""
                    self.LOG_FILE.append("    " + "    " + str(subval) + "," + "\n")
                self.LOG_FILE.append("    " + "]," + "\n")
        self.LOG_FILE.append(")" + "\n")

        self.backup = True
        self.pickle = False

        self.randseed = round(time.time())
        
        self.itp_defs = {
            "atomtypes":       {},
            "bondtypes":       {},
            "pairtypes":       {},
            "angletypes":      {},
            "dihedraltypes":   {},
            "constrainttypes": {},
        }
        self.itp_defs_all_defnames  = set()
        self.itp_moleculetypes      = {}
        self.ITP_INPUT_cmds         = []
        self.TOP_include_statements = []

        try:
            self.lipid_scaffolds = copy.deepcopy(lipid_scaffolds)
        except:
            self.print_term("WARNING: No lipid scaffolds found", warn=True)
            self.lipid_scaffolds = {}
        
        try:
            self.lipid_defs = copy.deepcopy(lipid_defs)
        except:
            self.print_term("WARNING: No lipid definitions found", warn=True)
            self.lipid_defs = {}
        
        try:
            self.solvent_defs = copy.deepcopy(solvent_defs)
        except:
            self.print_term("WARNING: No solvent definitions found", warn=True)
            self.solvent_defs = {}
        
        try:
            self.pos_ion_defs = copy.deepcopy(pos_ion_defs)
        except:
            self.print_term("WARNING: No positive ion definitions found", warn=True)
            self.pos_ion_defs = {}

        try:
            self.neg_ion_defs = copy.deepcopy(neg_ion_defs)
        except:
            self.print_term("WARNING: No negative ion definitions found", warn=True)
            self.neg_ion_defs = {}
        
        try:
            self.fragment_defs = copy.deepcopy(fragment_defs)
        except:
            self.print_term("WARNING: No fragment definitions found", warn=True)
            self.fragment_defs = {}

        self.lipid_dict   = {}
        self.solvent_dict = {}
        self.pos_ion_dict = {}
        self.neg_ion_dict = {}

        self.lipid_defs_built   = {}
        self.solvent_defs_built = {}
        self.pos_ion_defs_built     = {}
        self.neg_ion_defs_built     = {}

        if self.RUN:
            self.run(kwargs)

    ##############################
    ### GIVE COMMANDS TO CLASS ###
    ##############################
    def commands_handler(self, kwargs):
        for key, cmd in kwargs.items():
            ### Molecule definitions and scaffolds
            if key in ["import_library"]:
                if type(cmd) not in [list, tuple]:
                    cmd = [cmd]
                for subcmd in cmd:
                    assert subcmd.endswith(".py"), "Molecule definitions / lipid scaffolds / fragment definitions file must be a python file: '" + subcmd + "'"
                    
                    self.import_library(subcmd)

            ### Importing lipid fragment builder arguments from topology files
            elif key in ["itp_input", "itp_in"]:
                if type(cmd) != list:
                    cmd = [cmd]
                for subcmd in cmd:
                    self.ITP_INPUT_cmds.extend([subcmd])

            ### Molecule fragment builder
            elif key in ["molecule_builder"]:
                ### Puts individual string inside list
                if type(cmd) not in [list, tuple]:
                    cmd = [cmd]
                ### Converts tuple to list
                if type(cmd) == tuple:
                    cmd = list(cmd)
                for subcmd in cmd:
                    self.MOLECULE_FRAGMENT_BUILDER_cmds.extend([subcmd])
            
            ### Molecule builder
            elif key in ["molecule"]:
                ### Puts individual string inside list
                if type(cmd) not in [list, tuple]:
                    cmd = [cmd]
                ### Converts tuple to list
                if type(cmd) == tuple:
                    cmd = list(cmd)
                for subcmd in cmd:
                    self.MOLECULE_cmds.extend([subcmd])

            ### Outputs
            elif key in ["out_all", "o_all"]:
                ### Cuts the extension if present so that all files can be generated with proper extensions
                if any(cmd.lower().endswith(string) for string in [".pdb", ".gro", ".log"]):
                    cmd = cmd[:-4]
                self.output_system_pdb_file_name = cmd + ".pdb"
                self.output_system_gro_file_name = cmd + ".gro"
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
                if not cmd.lower().endswith(".pdb"):
                    cmd = cmd + ".pdb"
                self.output_system_pdb_file_name = cmd
                    
            elif key in ["out_gro", "o_gro"]:
                if not cmd.lower().endswith(".gro"):
                    cmd = cmd + ".gro"
                self.output_system_gro_file_name = cmd

            elif key in ["out_log", "o_log"]:
                if not cmd.lower().endswith(".log"):
                    cmd = cmd + ".log"
                self.output_log_file_name = cmd

            elif key in ["backup"]:
                assert cmd in ["False", "True", "0", "1", False, True, 0, 1], "Value given to 'backup' must be False/True/0/1 (strings allowed): " + cmd
                if type(cmd) == str:
                    cmd = ast.literal_eval(cmd)
                self.backup = cmd
            
            ### Misc
            elif key in ["rand", "randseed"]:
                cmd = self.get_number_from_string(cmd)
                assert type(cmd) in [int, float], "Value given to rand/randseed must be a number (string-numbers are allowed): " + str(cmd)
                self.randseed = round(cmd)
                    
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
                
            elif key in ["debug_keys"]:
                if type(cmd) == str:
                    self.debug_keys.append(cmd)
                elif type(cmd) in [list, tuple]:
                    for subcmd in cmd:
                        self.debug_keys.append(subcmd)
                
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
                assert False, "Invalid argument given to LipidBuilder: " + str((key, cmd))

        ### Setting randseed
        self.print_term("Setting random seed to:", self.randseed, verbose=1)
        random.seed(self.randseed)
        np.random.seed(self.randseed)

        ################################
        ### DEFINITION PREPROCESSING ###
        ################################
        string = " ".join(["", "PREPROCESSING DEFINITIONS", ""])
        self.print_term("{string:-^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), spaces=0, verbose=1)
        preprocessing_tic = time.time()

        self.itp_read_initiater()

        self.lipid_scaffolds_preprocessor()

        ### Fragment builder before defs as it adds to "self.lipid_defs_built", "self.solvent_defs_built", "self.pos_ion_defs_built" and "self.neg_ion_defs_built" dicts
        self.molecule_fragment_builder()

        self.lipid_defs_preprocessor()

        ### Preprocess ions first as they add the "default" parameter libraries of "neg_ions" and "pos_ions" to solvent defs
        self.ion_defs_preprocessor()
        self.solvent_defs_preprocessor()

        preprocessing_toc = time.time()
        preprocessing_time = round(preprocessing_toc - preprocessing_tic, 4)
        string1 = " ".join(["", "DEFINITIONS PREPROCESSING COMPLETE", ""])
        string2 = " ".join(["", "(Time spent:", str(preprocessing_time), "[s])", ""])
        self.print_term("{string:-^{string_length}}".format(string=string1, string_length=self.terminalupdate_string_length), spaces=0, verbose=1)
        self.print_term("{string:^{string_length}}".format(string=string2, string_length=self.terminalupdate_string_length), spaces=0, verbose=1)
        self.print_term("", spaces=0, verbose=1)

    def run(self, kwargs):
        self.commands_handler(kwargs)

        assert len(self.MOLECULE_cmds) + len(self.MOLECULE_FRAGMENT_BUILDER_cmds) > 0, (
            "Running 'MoleculeBuilder' requires at least one argument of one of the following types: 'molecule', 'lipid_builder'"
        )

        ########################
        ### MOLECULE BUILDER ###
        ########################
        string = " ".join(["", "CRAFING MOLCULE SYSTEMS", ""])
        self.print_term("{string:-^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), spaces=0, verbose=1)
        molcrafter_tic = time.time()
        self.molecule_crafter()
        molcrafter_toc = time.time()
        molcrafter_time = round(molcrafter_toc - molcrafter_tic, 4)
        string1 = " ".join(["", "MOLECULE SYSTEM CRAFTING COMPLETE", ""])
        string2 = " ".join(["", "(Time spent:", str(molcrafter_time), "[s])", ""])
        self.print_term("{string:-^{string_length}}".format(string=string1, string_length=self.terminalupdate_string_length), spaces=0, verbose=1)
        self.print_term("{string:^{string_length}}".format(string=string2, string_length=self.terminalupdate_string_length), spaces=0, verbose=1)
        self.print_term("", spaces=0, verbose=1)

        ####################
        ### FILE WRITERS ###
        ####################
        string = " ".join(["", "WRITING FILES", ""])
        self.print_term("{string:-^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), spaces=0, verbose=1)
        filewriting_tic = time.time()

        self.molecule_system_writer()
        self.log_file_writer()
        
        filewriting_toc = time.time()
        filewriting_time = round(filewriting_toc - filewriting_tic, 4)
        string1 = " ".join(["", "FILE WRITING COMPLETE", ""])
        string2 = " ".join(["", "(Time spent:", str(filewriting_time), "[s])", ""])
        self.print_term("{string:-^{string_length}}".format(string=string1, string_length=self.terminalupdate_string_length), spaces=0, verbose=1)
        self.print_term("{string:^{string_length}}".format(string=string2, string_length=self.terminalupdate_string_length), spaces=0, verbose=1)
        self.print_term("", spaces=0, verbose=1)

        ######################
        ### END OF PROGRAM ###
        ######################
        self.print_term("My task is complete. Did i do a good job?", verbose=1)
        COBY_run_toc  = time.time()
        COBY_run_time = round(COBY_run_toc - self.COBY_run_tic, 4)
        self.print_term("Time spent running COBY:", COBY_run_time, verbose=1)

    def molecule_crafter(self):
        '''
        Returns each molecule as a separate pdb/gro file with only one copy of the given molecule
        '''

        ###################
        ### FILE NAMING ###
        ###################
        self.print_term("Construction system names", spaces=0, verbose=2)

        ### Create initial dicts for each argument given to molecule builder
        cmd_dicts = []
        for cmd in self.MOLECULE_cmds:
            cmd_dict = {}
            cmd_dict["params"] = "default"
            cmd_dict.update({key: val for key, val in [subcmd.split(":", maxsplit=1) for subcmd in cmd.split()]})

            ### Check that values are properly given
            assert "moltype" in cmd_dict.keys() and cmd_dict["moltype"] in ["lipid", "solvent", "solute", "pos_ion", "neg_ion"], "\n".join([
                "The subargument 'moltype' for the argument '{cmd}' was not specified or the 'moltype' was not valid.".format(cmd=cmd),
                "The molecule type of the molecule to be created must be specified using the subargument 'moltype'.",
                "The accepted molecule types are:",
                "    ", "'lipid', 'solvent', 'solute', 'pos_ion' and 'neg_ion'",
            ])
            assert "name" in cmd_dict.keys(), "A molecule name must be specified using the subargument 'name'"

            cmd_dict["filename"] = cmd_dict["name"]
            cmd_dicts.append(copy.deepcopy(cmd_dict))
        
        def listofdicts_counter(l, name):
            '''
            Counts the number of occurences of a given filename in a list of dicts
            '''
            counter = 0
            for d in l:
                if name == d["filename"]:
                    counter += 1
            return counter

        ### First check if duplicate names (from multiple parameter libraries or molecule types) and add parameter library name to filename
        cmd_dicts = [
            {"moltype": d["moltype"], "name": d["name"], "params": d["params"], "filename": d["filename"] + "_" + d["params"]}
            if listofdicts_counter(cmd_dicts, d["filename"]) > 1
            else d
            for d in cmd_dicts
        ]
        ### Secondly check if there are still duplicate names (from multiple molecule types) and add molecule type to filename
        cmd_dicts = [
            {"moltype": d["moltype"], "name": d["name"], "params": d["params"], "filename": d["filename"] + "_" + d["moltype"]}
            if listofdicts_counter(cmd_dicts, d["filename"]) > 1
            else d
            for d in cmd_dicts
        ]

        ### Removes any remaining duplicates as it is pointless to process them multiple times
        cmd_dicts_new = []
        filenames = []
        for d in cmd_dicts:
            if d["filename"] not in filenames:
                filenames.append(d["filename"])
                cmd_dicts_new.append(copy.deepcopy(d))
        cmd_dicts = cmd_dicts_new
        
        ###############################################
        ### PROCESS MOLECULE BUILDER ARGUMENT DICTS ###
        ###############################################
        self.print_term("Building molecules", spaces=0, verbose=2)

        for cmd_nr, cmd_dict in enumerate(cmd_dicts, 1):
            self.print_term("Starting molecule number {nr} with name '{name}' of type '{moltype}' from parameter library '{params}'".format(nr=cmd_nr, name=cmd_dict["name"], moltype=cmd_dict["moltype"], params=cmd_dict["params"]), spaces=1, verbose=2)
            
            ###########################
            ### PREPROCESS ARGUMENT ###
            ###########################
            ### Designate explicit variables to avoid unnecessary dictionary lookups (also looks less cluttered)
            moltype, params, name, molecule_suffix = cmd_dict["moltype"], cmd_dict["params"], cmd_dict["name"], cmd_dict["filename"]

            ### Obtain molecule from corresponding molecule type dictionary
            if moltype == "lipid":
                assert params in self.lipid_dict.keys(), "\n".join([
                    "The lipid parameter library '{params}' was not found".format(params=params),
                    "The available lipid parameter libraries are:"
                    "    ", " ".join(list(self.lipid_dict.keys()))
                ])
                assert name in self.lipid_dict[params].keys(), "\n".join([
                    "The lipid '{name}' was not found in the lipid parameter library '{params}'".format(name=name, params=params),
                    "The available lipids in the lipid parameter library '{params}' are:".format(params=params),
                    "    ", " ".join(list(self.lipid_dict[params].keys()))
                ])
                molecule = copy.deepcopy(self.lipid_dict[params][name])

            elif moltype in ["solvent", "solute"]:
                assert params in self.solvent_dict.keys(), "\n".join([
                    "The solvent/solute parameter library '{params}' was not found".format(params=params),
                    "The available solvent/solute parameter libraries are:"
                    "    ", " ".join(list(self.solvent_dict.keys()))
                ])
                assert name in self.solvent_dict[params].keys(), "\n".join([
                    "The solvent/solute '{name}' was not found in the solvent/solute parameter library '{params}'".format(name=name, params=params),
                    "The available solvents/solutes in the solvent/solute parameter library '{params}' are:".format(params=params),
                    "    ", " ".join(list(self.solvent_dict[params].keys()))
                ])
                molecule = copy.deepcopy(self.solvent_dict[params][name])

            elif moltype in ["pos_ion"]:
                assert params in self.pos_ion_dict.keys(), "\n".join([
                    "The ion parameter library '{params}' was not found".format(params=params),
                    "The available ion parameter libraries are:"
                    "    ", " ".join(list(self.pos_ion_dict.keys()))
                ])
                assert name in self.pos_ion_dict[params].keys(), "\n".join([
                    "The positive ion '{name}' was not found in the ion parameter library '{params}'".format(name=name, params=params),
                    "The available positive ions in the ion parameter library '{params}' are:".format(params=params),
                    "    ", " ".join(list(self.pos_ion_dict[params].keys()))
                ])
                molecule = copy.deepcopy(self.pos_ion_dict[params][name])

            elif moltype in ["neg_ion"]:
                assert params in self.neg_ion_dict.keys(), "\n".join([
                    "The ion parameter library '{params}' was not found".format(params=params),
                    "The available ion parameter libraries are:"
                    "    ", " ".join(list(self.neg_ion_dict.keys()))
                ])
                assert name in self.neg_ion_dict[params].keys(), "\n".join([
                    "The negative ion '{name}' was not found in the ion parameter library '{params}'".format(name=name, params=params),
                    "The available negative ions in the ion parameter library '{params}' are:".format(params=params),
                    "    ", " ".join(list(self.neg_ion_dict[params].keys()))
                ])
                molecule = copy.deepcopy(self.neg_ion_dict[params][name])

            ################################################
            ### CENTER MOLECULE AND OBTAIN SYSTEM VALUES ###
            ################################################

            ### Centering molecule with a 0.5 nm buffer on each side
            bead_xs, bead_ys, bead_zs = molecule.get_coords()
            edgebuffer = 5 # [Å]
            xmin, ymin, zmin = min(bead_xs), min(bead_ys), min(bead_zs)
            beads_centered = [(bead_name, beadnr, x - xmin + edgebuffer, y - ymin + edgebuffer, z - zmin + edgebuffer, resname, resnr, charge) for bead_name, beadnr, x, y, z, resname, resnr, charge in molecule.get_res_beads_info(output_type="tuple")]

            ### Get PDB and GRO box parameters
            bead_name, beadnrs, xs, ys, zs, resnames, resnrs, charges = list(zip(*beads_centered))
            xlen, ylen, zlen = max(xs) + edgebuffer, max(ys) + edgebuffer, max(zs) + edgebuffer
            gro_box_vectors = [
                float(xlen/10), # vector: vax or v1x
                float(ylen/10), # vector: vby or v2y
                float(zlen/10), # vector: vcz or v3z
                float(0),       # vector: vay or v1y
                float(0),       # vector: vaz or v1z
                float(0),       # vector: vbx or v2x
                float(0),       # vector: vbz or v2z
                float(0),       # vector: vcx or v3x
                float(0),       # vector: vcy or v3y
            ]
            pdb_box_dimension = [
                float(xlen), # Axis length:  x
                float(ylen), # Axis length:  y
                float(zlen), # Axis length:  z
                float(90),   # Corner angle: alpha
                float(90),   # Corner angle: beta
                float(90),   # Corner angle: gamma
            ]

            ### Get system names
            system_name = molecule_suffix

            if not self.output_system_pdb_file_name and not self.output_system_gro_file_name:
                # molecule_specific_output_system_pdb_file_name = "molecule" + "_" + molecule_suffix + ".pdb"
                # molecule_specific_output_system_gro_file_name = "molecule" + "_" + molecule_suffix + ".gro"
                molecule_specific_output_system_pdb_file_name = molecule_suffix + ".pdb"
                molecule_specific_output_system_gro_file_name = molecule_suffix + ".gro"
            else:
                if self.output_system_pdb_file_name:
                    # molecule_specific_output_system_pdb_file_name = self.output_system_pdb_file_name.rstrip(".pdb") + "_" + molecule_suffix + ".pdb"
                    molecule_specific_output_system_pdb_file_name = self.output_system_pdb_file_name.rstrip(".pdb").rstrip(".") + molecule_suffix + ".pdb"
                else:
                    molecule_specific_output_system_pdb_file_name = False
                
                if self.output_system_gro_file_name:
                    # molecule_specific_output_system_gro_file_name = self.output_system_gro_file_name.rstrip(".gro") + "_" + molecule_suffix + ".gro"
                    molecule_specific_output_system_gro_file_name = self.output_system_gro_file_name.rstrip(".gro").rstrip(".") + molecule_suffix + ".gro"
                else:
                    molecule_specific_output_system_gro_file_name = False
            
            self.MOLECULES.append((beads_centered, gro_box_vectors, pdb_box_dimension, molecule_specific_output_system_pdb_file_name, molecule_specific_output_system_gro_file_name, system_name))
            
        self.print_term("Number of molecules built:", len(cmd_dicts), spaces=0, verbose=2)
        
    def molecule_system_writer(self):

        #############################################
        ### WRITING STRUCTURE FILES FOR MOLECULES ###
        #############################################
        self.print_term("Writing structure files (PDB/GRO)", spaces=0, verbose=1)

        for beads_centered, gro_box_vectors, pdb_box_dimension, molecule_specific_output_system_pdb_file_name, molecule_specific_output_system_gro_file_name, system_name in self.MOLECULES:
            output_system_pdb_file_lines = []
            output_system_gro_file_lines = []
            ###############################
            ### BEGINNING OF FILE LINES ###
            ###############################
            if molecule_specific_output_system_pdb_file_name:
                a, b, c, alpha, beta, gamma = pdb_box_dimension
                output_system_pdb_file_lines = [
                    "TITLE     " + system_name,
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
            if molecule_specific_output_system_gro_file_name:
                output_system_gro_file_lines = [
                    system_name,
                    "PLACEHOLDER_ATOM_COUNT",
                ]
            
            ################################
            ### MOLECULE BEAD/ATOM LINES ###
            ################################

            atom_count       = 0
            atom_nr          = 0
            res_nr           = 0
            last_bead_res_nr = 0

            for bead_name, bead_nr, bead_x, bead_y, bead_z, bead_res_name, bead_res_nr, bead_charge in beads_centered:
                if last_bead_res_nr != bead_res_nr:
                    res_nr += 1
                last_bead_res_nr = bead_res_nr
                if res_nr >= 10000:
                    res_nr -= 10000 * (res_nr // 10000)
                res_name = bead_res_name

                atom_nr += 1
                atom_count += 1
                if atom_nr >= 100000:
                    atom_nr -= 100000 * (atom_nr // 100000)

                x, y, z = bead_x, bead_y, bead_z

                atom_name = bead_name
                if molecule_specific_output_system_pdb_file_name:
                    output_system_pdb_file_lines.append(self.pdb_atom_writer("ATOM", atom_nr, atom_name, " ", res_name, "A", res_nr, " ", float(x), float(y), float(z), float(1), float(0), " ", " ", " "))
                if molecule_specific_output_system_gro_file_name: ### gro coordinates are in [nm] not [Å]
                    output_system_gro_file_lines.append(self.gro_atom_writer(res_nr, res_name, atom_name, atom_nr, x / 10, y / 10, z / 10, " ", " ", " "))

            #########################
            ### END OF FILE LINES ###
            #########################
            if molecule_specific_output_system_pdb_file_name:
                output_system_pdb_file_lines.append("TER")
                output_system_pdb_file_lines.append("END")
            
            if molecule_specific_output_system_gro_file_name:
                v1x, v2y, v3z, v1y, v1z, v2x, v2z, v3x, v3y = gro_box_vectors
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
                # output_system_gro_file_lines.append("") # Apparantly should not be present

            ##############################
            ### BACKUP AND WRITE FILES ###
            ##############################
            ### PDB
            if molecule_specific_output_system_pdb_file_name:
                if self.backup:
                    if molecule_specific_output_system_pdb_file_name:
                        self.backupper(molecule_specific_output_system_pdb_file_name)
                new_file = open(molecule_specific_output_system_pdb_file_name, "w")
                for line in output_system_pdb_file_lines:
                    new_file.write(line + "\n")
                new_file.close()
                self.print_term("PDB file written:", molecule_specific_output_system_pdb_file_name, spaces=1, verbose=1)
            
            ### GRO
            if molecule_specific_output_system_gro_file_name:
                if self.backup:
                    if molecule_specific_output_system_gro_file_name:
                        self.backupper(molecule_specific_output_system_gro_file_name)
                new_file = open(molecule_specific_output_system_gro_file_name, "w")
                for line in output_system_gro_file_lines:
                    new_file.write(line + "\n")
                new_file.close()
                self.print_term("GRO file written:", molecule_specific_output_system_gro_file_name, spaces=1, verbose=1)

        self.print_term("", verbose=1)
