import time

import ast
import math
import numpy as np
import random
import copy
import os

### Placeholders in case none are defined in "COBY.molecule_definitions.__init__" and "COBY.fragment_definitions.__init__"
lipid_scaffolds = {}
lipid_defs      = {}
solvent_defs    = {}
pos_ion_defs    = {}
neg_ion_defs    = {}
prot_defs       = {}
fragment_defs   = {}

lipid_metadata    = {}
solvent_metadata  = {}
pos_ion_metadata  = {}
neg_ion_metadata  = {}
prot_metadata     = {}
fragment_metadata = {}

from COBY.molecule_definitions.__init__ import *
from COBY.fragment_definitions.__init__ import *
from COBY.structure_classes.__init__ import *
from COBY.general_functions.__init__ import *

from COBY.main_class.structure_file_handlers.__init__ import *
from COBY.main_class.topology_handlers.__init__ import *
from COBY.main_class.general_tools.__init__ import *
from COBY.main_class.definition_preprocessors.__init__ import *
from COBY.main_class.molecule_fragment_builder.__init__ import *

class Library(
    structure_file_handlers,
    topology_handlers,
    general_tools,
    definition_preprocessors,
    molecule_fragment_builder,
):
    
    def __init__(self, run = True, terminal_run_kwargs = False, **kwargs):
        self.COBY_run_tic = time.time()
        self.PROGRAM = "COBY"

        self.RUN = run
        
        self.debug_prints = False
        self.debug_keys   = []
        self.extra_info   = True
        self.warnings     = True
        self.quiet        = False
        self.verbose      = 1
        
        if terminal_run_kwargs:
            kwargs.update(terminal_run_kwargs)
        
        self.MOLECULE_IMPORT_cmds = []

        self.MOLECULE_FRAGMENT_BUILDER_cmds = []
        
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
        self.itp_defs_all_defnames  = set()
        self.itp_moleculetypes      = {}
        self.ITP_INPUT_cmds         = []
        self.TOP_include_statements = []

        self.system_charge = 0
        self.system_name = "PLACEHOLDER_TITLE"
        
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
        
        try:
            self.prot_defs = copy.deepcopy(prot_defs)
        except:
            self.print_term("WARNING: No protein charge definitions found", warn=True)
            self.prot_defs = {}

        try:
            self.lipid_metadata = copy.deepcopy(lipid_metadata)
        except:
            self.print_term("WARNING: No lipid metadata found", warn=True)
            self.lipid_metadata = {}

        try:
            self.solvent_metadata = copy.deepcopy(solvent_metadata)
        except:
            self.print_term("WARNING: No solvent metadata found", warn=True)
            self.solvent_metadata = {}

        try:
            self.pos_ion_metadata = copy.deepcopy(pos_ion_metadata)
        except:
            self.print_term("WARNING: No positive ion metadata found", warn=True)
            self.pos_ion_metadata = {}

        try:
            self.neg_ion_metadata = copy.deepcopy(neg_ion_metadata)
        except:
            self.print_term("WARNING: No negative ion metadata found", warn=True)
            self.neg_ion_metadata = {}

        try:
            self.prot_metadata = copy.deepcopy(prot_metadata)
        except:
            self.print_term("WARNING: No fragment metadata found", warn=True)
            self.prot_metadata = {}

        try:
            self.fragment_metadata = copy.deepcopy(fragment_metadata)
        except:
            self.print_term("WARNING: No fragment metadata found", warn=True)
            self.fragment_metadata = {}

        
        self.lipid_dict   = {}
        self.solvent_dict = {}
        self.pos_ion_dict = {}
        self.neg_ion_dict = {}
        self.prot_dict    = {}

        self.lipid_defs_built   = {}
        self.solvent_defs_built = {}
        self.pos_ion_defs_built     = {}
        self.neg_ion_defs_built     = {}

        self.lipid_defs_imported   = {}
        self.solvent_defs_imported = {}
        self.pos_ion_defs_imported     = {}
        self.neg_ion_defs_imported     = {}
        
        if self.RUN:
            self.run(kwargs)
    
    ##############################
    ### GIVE COMMANDS TO CLASS ###
    ##############################
    def commands_handler(self, kwargs):
        invalid_args_given = False

        for key, cmd in kwargs.items():
            ### Molecule definitions and scaffolds
            if key in ["import_library"]:
                if type(cmd) not in [list, tuple]:
                    cmd = [cmd]
                for subcmd in cmd:
                    assert subcmd.endswith(".py"), "Molecule definitions / lipid scaffolds / fragment definitions file must be a python file: '" + subcmd + "'"
                    if subcmd.startswith("file:"):
                        assert subcmd != "file:", "'file' subargument given to 'import_library' argument but no destination was given alongside 'file:'."
                        subcmd = subcmd[subcmd.index(":")+1:]
                    
                    self.import_library(subcmd)
            
            ### Importing charges and lipid fragment builder arguments from topology files
            elif key in ["itp_input", "itp_in"]:
                if type(cmd) != list:
                    cmd = [cmd]
                for subcmd in cmd:
                    self.ITP_INPUT_cmds.extend([subcmd])
            
            ### Importing molecules from pdb/gro files
            elif key in ["molecule_import"]:
                ### Puts individual string inside list
                if type(cmd) not in [list, tuple]:
                    cmd = [cmd]
                ### Converts tuple to list
                if type(cmd) == tuple:
                    cmd = list(cmd)
                for subcmd in cmd:
                    self.MOLECULE_IMPORT_cmds.extend([subcmd])
            
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
            
            ### Outputs
            elif key in ["out_all", "o_all"]:
                ### Cuts the extension if present so that all files can be generated with proper extensions
                if any(cmd.lower().endswith(string) for string in [".log"]):
                    cmd = cmd[:-4]
                self.output_log_file_name = os.path.join(cmd + ".log")
            
            elif key in ["out_log", "o_log"]:
                if not cmd.lower().endswith(".log"):
                    cmd = cmd + ".log"
                self.output_log_file_name = os.path.join(cmd)
                
            elif key in ["backup"]:
                assert cmd in ["False", "True", "0", "1", False, True, 0, 1], "Value given to 'backup' must be False/True/0/1 (strings allowed): " + cmd
                if type(cmd) == str:
                    cmd = ast.literal_eval(cmd)
                self.backup = cmd
            
            elif key in ["rand", "randseed"]:
                cmd = self.get_number_from_string(cmd)
                assert type(cmd) in [int, float], "Value given to rand/randseed must be a number (string-numbers are allowed): " + str(cmd)
                self.randseed = round(cmd)
                    
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
                
            elif key in ["debug_keys"]:
                if type(cmd) == str:
                    self.debug_keys.append(cmd)
                elif type(cmd) in [list, tuple]:
                    for subcmd in cmd:
                        self.debug_keys.append(subcmd)
                
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
            elif key == "run":
                if type(cmd) == str:
                    cmd = ast.literal_eval(cmd)
                self.RUN = cmd
            
            else:
                if not invalid_args_given:
                    self.print_term("This is the COBY.Library() method. Arguments given to COBY.Library() that are not related to molecule definition file importing, molecule structure importing, topology importing or the molecule fragment builder will not stop the program, but will also not be processed.", warn=True)
                    invalid_args_given = True
                if type(cmd) == str:
                    self.print_term("Invalid argument detected:", '"'+key, "=", "'"+cmd+"'"+'"', spaces=1, warn=True)
                else:
                    self.print_term("Invalid argument detected:", '"'+key, "=", str(cmd)+'"', spaces=1, warn=True)
        
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
        
        self.molecule_importer()

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
        '''
        Runs the entire system creation process
        '''
        
        self.commands_handler(kwargs)

        self.ILR_layer0_main()

        ####################
        ### FILE WRITING ###
        ####################
        ### Only log files can be written by COBY.Library, so checks if terminal printing should be done:
        if self.output_log_file_name:
            string = " ".join(["", "WRITING FILES", ""])
            self.print_term("{string:-^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), spaces=0, verbose=1)
            filewriting_tic = time.time()

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

    def ILR_invalid_answer(self, val):
        self.print_term("I did not understand your answer: '{}'.".format(val), verbose=0)
        self.print_term("", verbose=0)

    def print_mol_data(self, moltype, parameter_library, mol_name):
        if moltype == "lipid":
            argtypes = ["membrane"]
            subargtypes = ["lipid"]
            mol_dict = self.lipid_dict
            moltype_longname = "lipid"
        
        elif moltype == "solvent":
            argtypes = ["solvation", "flooding"]
            subargtypes = ["solv", "solute"]
            mol_dict = self.solvent_dict
            moltype_longname = "solvent"
        
        elif moltype == "pos_ion":
            argtypes = ["solvation"]
            subargtypes = ["pos"]
            mol_dict = self.pos_ion_dict
            moltype_longname = "positive ion"
        
        elif moltype == "neg_ion":
            argtypes = ["solvation"]
            subargtypes = ["neg"]
            mol_dict = self.neg_ion_dict
            moltype_longname = "negative ion"
        else:
            assert False, "Invalid 'moltype' recieved by 'print_mol_data': {moltype}".format(moltype=moltype)

        self.print_term("Below information is for the {moltype_longname} '{mol_name}' from the parameter library '{parameter_library}'".format(moltype_longname=moltype_longname, mol_name=mol_name, parameter_library=parameter_library), verbose=0)
        self.print_term("General data:",                                                                        verbose=0, spaces=1)
        self.print_term("Number of residues:", mol_dict[parameter_library][mol_name].n_residues,                verbose=0, spaces=2)
        self.print_term("Number of beads:", len(mol_dict[parameter_library][mol_name].get_bead_charges()),      verbose=0, spaces=2)
        self.print_term("Total charge:", mol_dict[parameter_library][mol_name].get_mol_charge(),                verbose=0, spaces=2)
        self.print_term("Tags:", " ".join(["'"+tag+"'" for tag in mol_dict[parameter_library][mol_name].tags]), verbose=0, spaces=2)

        self.print_term("Residues:",                                                                          verbose=0, spaces=1)
        for residue in mol_dict[parameter_library][mol_name].residues:
            self.print_term("Residue number "+str(residue.resnr)+":",                                         verbose=0, spaces=2)
            self.print_term("Residue name:", residue.resname,                                                 verbose=0, spaces=3)
            self.print_term("Number of beads in residue:", len(residue.beads),                                verbose=0, spaces=3)
            self.print_term("Total charge of residue:", round(sum(bead.charge for bead in residue.beads), 3), verbose=0, spaces=3)
            self.print_term("Beads:",                                                                         verbose=0, spaces=3)
            table_values  = {"numbers": [], "names": [], "xs": [], "ys": [], "zs": [], "charges": []}
            for bead in residue.beads:
                table_values["numbers"].append(str(bead.beadnr))
                table_values["names"].append(bead.bead)
                table_values["xs"].append(str(round(bead.x, 3)))
                table_values["ys"].append(str(round(bead.y, 3)))
                table_values["zs"].append(str(round(bead.z, 3)))
                table_values["charges"].append(str(round(bead.charge, 3)))
            
            ### Centers all numbers on the decimal position.
            for key in ["numbers", "xs", "ys", "zs", "charges"]:
                ### Adds "+" to all positive values to help with alignment if any value is negative. Plusses are removed later and replaced with empty spaces " ".
                if any(val.startswith("-") for val in table_values[key]):
                    table_values[key] = [val if val.startswith("-") else "+"+val for val in table_values[key]]

                ### Splitting values based on decimal position
                split_vals = [val.split(".") for val in table_values[key]]

                ### Obtaining number of digits to the left and right of the decimal position
                max_left   = max(len(split_val[0]) for split_val in split_vals)
                max_right  = max(len(split_val[1]) if len(split_val) > 1 else 0 for split_val in split_vals)

                ### Adds blank spaces (with rjust and ljust) to each value to align them according to decimal position
                for i, split_val in enumerate(split_vals):
                    left_digits = split_val[0]
                    right_digits = split_val[1] if len(split_val) > 1 else ''
                    table_values[key][i] = left_digits.rjust(max_left) + "." + right_digits.ljust(max_right)
                table_values[key] = [val.replace("+", " ") for val in table_values[key]]
                table_values[key] = [val.replace(".", " ") if len(val.split(".")[1]) == 0 else val for val in table_values[key]]

            columns = [["Number"], ["Name"], ["X"], ["Y"], ["Z"], ["Charge"]]
            for beadi in range(len(table_values["names"])):
                for vali, val in enumerate([table_values["numbers"][beadi], table_values["names"][beadi], table_values["xs"][beadi], table_values["ys"][beadi], table_values["zs"][beadi], table_values["charges"][beadi]]):
                    columns[vali].append(str(val))

            max_column_lengths = [max([len(val) for val in col]) for col in columns]

            tot_length = sum(max_column_lengths) + len(" : ")*5 + len(" ")*2
            for rowi in range(len(columns[0])):
                string = '{nr:^{L0}} : {name:^{L1}} : {x:^{L2}} : {y:^{L3}} : {z:^{L4}} : {charge:^{L5}}'.format(
                    nr     = columns[0][rowi], L0 = max_column_lengths[0],
                    name   = columns[1][rowi], L1 = max_column_lengths[1],
                    x      = columns[2][rowi], L2 = max_column_lengths[2],
                    y      = columns[3][rowi], L3 = max_column_lengths[3],
                    z      = columns[4][rowi], L4 = max_column_lengths[4],
                    charge = columns[5][rowi], L5 = max_column_lengths[5],
                )

                self.print_term(string, verbose=0, spaces=4)


        for argtype, subargtype in zip(argtypes, subargtypes):
            self.print_term("", verbose=0)
            self.print_term("This {moltype_longname} can be accessed in '{argtype}' arguments via the following subargument:".format(moltype_longname=moltype_longname, argtype=argtype), verbose=0, spaces=1)
            self.print_term("'{subargtype}:params:{parameter_library}:name:{mol_name}'".format(subargtype=subargtype, mol_name=mol_name, parameter_library=parameter_library),            verbose=0, spaces=2)

        self.print_term("", verbose=0)

    def print_mols_with_tags(self, moltype, parameter_libraries, tagsarg):
        if moltype == "lipid":
            mol_dict = self.lipid_dict
            moltype_longname = "lipid"
        
        elif moltype == "solvent":
            mol_dict = self.solvent_dict
            moltype_longname = "solvent"
        
        elif moltype == "pos_ion":
            mol_dict = self.pos_ion_dict
            moltype_longname = "positive ion"
        
        elif moltype == "neg_ion":
            mol_dict = self.neg_ion_dict
            moltype_longname = "negative ion"
        else:
            assert False, "Invalid 'moltype' recieved by 'print_mol_data': {moltype}".format(moltype=moltype)

        for parameter_library in parameter_libraries:
            found_mols_with_tags = []
            func_str = False
            func = False
            given_tags = []
            breaked_loop = False
            for subval in tagsarg.split(":") + ["end"]:
                if func_str is not False and subval.lower() in ["tag", "any", "all", "end"]:
                    if len(given_tags) == 0:
                        self.print_term("No tags have ben given to the '{}' subargument. Please try again.".format(str(func_str)), verbose=0)
                        self.print_term("", verbose=0)
                        breaked_loop = True
                        break

                    if func_str == "tag":
                        if len(given_tags) > 1:
                            self.print_term("More than 1 tag has been given to the 'tag' subargument. Please try again with only 1 tag.".format(str(func_str)), verbose=0)
                            self.print_term("", verbose=0)
                            breaked_loop = True
                            break

                        found_mols_with_tags.append(sorted(list(set([
                            key
                            for key, val in mol_dict[parameter_library].items()
                            if given_tags[0] in val.tags and len(val.tags) > 0
                        ]))))
                    else:
                        found_mols_with_tags.append(sorted(list(set([
                            key
                            for key, val in mol_dict[parameter_library].items()
                            if func([tag in val.tags for tag in given_tags]) and len(val.tags) > 0
                        ]))))
                    func_str = False
                    func = False
                    given_tags = []
                
                if subval.lower() == "tag":
                    func_str = "tag"
                elif subval.lower() == "any":
                    func_str = "any"
                    func = any
                elif subval.lower() == "all":
                    func_str = "all"
                    func = all
                elif subval.lower() == "end":
                    break
                else:
                    given_tags.append(subval)

            if not breaked_loop:
                mols_with_tags_in_params = {}
                for l in found_mols_with_tags:
                    for m in l:
                        if m not in mols_with_tags_in_params:
                            mols_with_tags_in_params[m] = 1
                        else:
                            mols_with_tags_in_params[m] += 1
                mols_with_tags_in_params = [key for key, val in mols_with_tags_in_params.items() if val == len(found_mols_with_tags)]

                self.print_term("Parameter library name: '{parameter_library}'".format(parameter_library=parameter_library), verbose=0)
                if len(mols_with_tags_in_params) > 0:
                    self.print_term(" ".join(mols_with_tags_in_params), verbose=0, spaces=1)
                self.print_term("", verbose=0)

    def print_prot_charge_data(self, parameter_library, res_name):
        self.print_term("Below information is for the protein residue '{res_name}' from the parameter library '{parameter_library}'".format(res_name=res_name, parameter_library=parameter_library), verbose=0)
        self.print_term("General data:",                                                                        verbose=0, spaces=1)
        self.print_term("Number of beads:", len(self.prot_defs[parameter_library]["charges"][res_name].keys()), verbose=0, spaces=2)
        self.print_term("Total charge:", sum(self.prot_defs[parameter_library]["charges"][res_name].values()),  verbose=0, spaces=2)
        self.print_term("Beads:",                                                                               verbose=0, spaces=1)
        table_values  = {"numbers": [], "names": [], "charges": []}
        for i, (name, charge) in enumerate(self.prot_defs[parameter_library]["charges"][res_name].items()):
            table_values["numbers"].append(str(i))
            table_values["names"].append(name)
            table_values["charges"].append(str(round(charge, 3)))
            
        ### Centers all numbers on the decimal position.
        for key in ["numbers", "charges"]:
            ### Adds "+" to all positive values to help with alignment if any value is negative. Plusses are removed later and replaced with empty spaces " ".
            if any(val.startswith("-") for val in table_values[key]):
                table_values[key] = [val if val.startswith("-") else "+"+val for val in table_values[key]]

            ### Splitting values based on decimal position
            split_vals = [val.split(".") for val in table_values[key]]

            ### Obtaining number of digits to the left and right of the decimal position
            max_left   = max(len(split_val[0]) for split_val in split_vals)
            max_right  = max(len(split_val[1]) if len(split_val) > 1 else 0 for split_val in split_vals)

            ### Adds blank spaces (with rjust and ljust) to each value to align them according to decimal position
            for i, split_val in enumerate(split_vals):
                left_digits = split_val[0]
                right_digits = split_val[1] if len(split_val) > 1 else ''
                table_values[key][i] = left_digits.rjust(max_left) + "." + right_digits.ljust(max_right)
            table_values[key] = [val.replace("+", " ") for val in table_values[key]]
            table_values[key] = [val.replace(".", " ") if len(val.split(".")[1]) == 0 else val for val in table_values[key]]

        columns = [["Number"], ["Name"], ["Charge"]]
        for beadi in range(len(table_values["names"])):
            for vali, val in enumerate([table_values["numbers"][beadi], table_values["names"][beadi], table_values["charges"][beadi]]):
                columns[vali].append(str(val))

        max_column_lengths = [max([len(val) for val in col]) for col in columns]

        tot_length = sum(max_column_lengths) + len(" : ")*2 + len(" ")*2
        for rowi in range(len(columns[0])):
            string = '{nr:^{L0}} : {name:^{L1}} : {charge:^{L2}}'.format(
                nr     = columns[0][rowi], L0 = max_column_lengths[0],
                name   = columns[1][rowi], L1 = max_column_lengths[1],
                charge = columns[2][rowi], L2 = max_column_lengths[2],
            )

            self.print_term(string, verbose=0, spaces=2)

        self.print_term("", verbose=0)

    ### ILR = interactive library roamer
    def ILR_layer0_main(self):
        ILR_restart_layer = self.ILR_layer0_main
        self.print_term("-"*self.terminalupdate_string_length, verbose=0)

        val = self.print_term(
            "\n".join([
                "What library would you like to investigate? Your options are shown below:",
                "    "+"Quit:                                      'q' or 'quit' ",
                "    "+"For the lipid library:                     'lipid(s)'",
                "    "+"For the solvent/solute library:            'solvent(s)'",
                "    "+"For the positive ion library:              'pos_ion(s)'",
                "    "+"For the negative ion library:              'neg_ion(s)'",
                "    "+"For the protein charge library:            'protein(s)'",
                "    "+"For the molecule fragment builder library: 'fragment(s)'",
                "",
            ]),
            inp=True
        )
        self.print_term("", verbose=0)
        val = val.lstrip(" ").rstrip(" ")

        if val.lower() in ["q", "quit"]:
            pass
        elif val.lower() in ["lipid", "lipids"]:
            self.ILR_layer1_lipid_solvent_ions("lipid")
        elif val.lower() in ["solvent", "solvents"]:
            self.ILR_layer1_lipid_solvent_ions("solvent")
        elif val.lower() in ["pos_ion", "pos_ions"]:
            self.ILR_layer1_lipid_solvent_ions("pos_ion")
        elif val.lower() in ["neg_ion", "neg_ions"]:
            self.ILR_layer1_lipid_solvent_ions("neg_ion")
        elif val.lower() in ["protein", "proteins"]:
            self.ILR_layer1_protein()
        elif val.lower() in ["fragment", "fragments"]:
            self.ILR_layer1_MFB()
        else:
            self.ILR_invalid_answer(val)
            ILR_restart_layer()

    def ILR_layer1_lipid_solvent_ions(self, moltype):
        ILR_restart_layer = self.ILR_layer1_lipid_solvent_ions
        restart_dict={"moltype": moltype}
        self.print_term("-"*self.terminalupdate_string_length, verbose=0)

        if moltype == "lipid":
            mol_dict = self.lipid_dict
            moltype_longname = "lipid"
        
        elif moltype == "solvent":
            mol_dict = self.solvent_dict
            moltype_longname = "solvent"
        
        elif moltype == "pos_ion":
            mol_dict = self.pos_ion_dict
            moltype_longname = "positive ion"
        
        elif moltype == "neg_ion":
            mol_dict = self.neg_ion_dict
            moltype_longname = "negative ion"
        else:
            assert False, "Invalid 'moltype' recieved by 'ILR_layer1': {moltype}".format(moltype=moltype)
        
        parameter_libraries = sorted(list(mol_dict.keys()))
        all_tags_in_params = sorted(list(set([
            tag
            for parameter_library in parameter_libraries
            for mol in mol_dict[parameter_library].values()
            for tag in mol.tags
        ])))

        longest_string = "Print all {moltype_longname} names in a parameter library:".format(moltype_longname=moltype_longname)
        longest_string_len = len(longest_string) + 1

        val = self.print_term(
            "\n".join([
                "You are in the general {moltype_longname} library. The names of all available {moltype_longname} parameter libraries are shown below:".format(moltype_longname=moltype_longname),
                "    "+" ".join(["'"+parameter_library+"'" for parameter_library in parameter_libraries]),
                "",
                "What would you like to do? Your options are shown below:",
                "    "+"Quit:".ljust(longest_string_len)                                                                                                + "'q' or 'quit'",
                "    "+"Return to previous question:".ljust(longest_string_len)                                                                         + "'r' or 'return'",
                "    "+"Examine parameter library:".ljust(longest_string_len)                                                                           + "'param(s):[parameter library name]' or '[parameter library name]'",
                "    "+"Print all {moltype_longname} names in a parameter library:".format(moltype_longname=moltype_longname).ljust(longest_string_len) + "'printparam(s):[parameter library name]' or 'pp:[parameter library name]'",
                "    "+""+"Print names of {moltype_longname}s with given tag(s):".format(moltype_longname=moltype_longname).ljust(longest_string_len),
                "    "+"    "+"- Only this specific tag:".ljust(longest_string_len)                                                                     + "'tag:[tag1]' (only 1 tag)",
                "    "+"    "+"- With any of these tags:".ljust(longest_string_len)                                                                     + "'any:[tag1]:[tag2]:[etc.]' (any number of tags)",
                "    "+"    "+"- With all of these tags:".ljust(longest_string_len)                                                                     + "'all:[tag1]:[tag2]:[etc.]' (any number of tags)",
                "    "+"    "+"- Combining the above:".ljust(longest_string_len)                                                                        + "'tag:[tag1]:any:[tag2]:[etc.]:all:[tag3]:[etc.]'",
                "    "+"Print all tags:".ljust(longest_string_len)                                                                                      + "'printtag(s)' or 'pt",
                "    "+"Examine {moltype_longname} from given parameter library:".format(moltype_longname=moltype_longname).ljust(longest_string_len)   + "'{moltype}(s):[parameter library name]:[{moltype_longname} name]'".format(moltype=moltype, moltype_longname=moltype_longname),
                "",
            ]),
            inp=True
        )
        self.print_term("", verbose=0)
        val = val.lstrip(" ").rstrip(" ")

        if val.lower() in ["q", "quit"]:
            pass
        
        elif val in ["r", "return"]:
            self.ILR_layer0_main()
        
        elif val.lower().startswith(("param:", "params:")) or val in parameter_libraries:
            parameter_library = False
            if val.lower().startswith(("param:", "params:")):
                if len(val.split(":")) == 2:
                    parameter_library = val.split(":")[1]
                else:
                    self.print_term("You must specify exactly one parameter library: '{val}'.".format(val=val), verbose=0)
                    self.print_term("", verbose=0)
            else:
                parameter_library = val
            
            if parameter_library is not False:
                if parameter_library in parameter_libraries:
                    self.ILR_layer2_lipid_solvent_ions(moltype=moltype, parameter_library=parameter_library)
                else:
                    self.print_term("Parameter library not found. You must specify an existing parameter library: '{val}'.".format(val=val), verbose=0)
                    self.print_term("", verbose=0)
                    ILR_restart_layer(**restart_dict)
            else:
                ILR_restart_layer(**restart_dict)
        
        elif val.lower().startswith(("printparam:", "printparams:", "pp:")):
            if len(val.split(":")) == 2:
                parameter_library = val.split(":")[1]
                if parameter_library in parameter_libraries:
                    self.print_term("Following {moltype_longname}s are present in the '{parameter_library}' parameter library:".format(moltype_longname=moltype_longname, parameter_library=parameter_library), verbose=0)
                    self.print_term(" ".join(list(mol_dict[parameter_library].keys())), verbose=0, spaces=1)
                    self.print_term("", verbose=0)
                else:
                    self.print_term("Parameter library not found. You must specify an existing parameter library: '{val}'.".format(val=val), verbose=0)
                    self.print_term("", verbose=0)
            else:
                self.print_term("You must specify exactly one parameter library: '{val}'.".format(val=val), verbose=0)
                self.print_term("", verbose=0)
            ILR_restart_layer(**restart_dict)

        elif val.lower().startswith(("tag:", "any:" ,"all:")):
            self.print_mols_with_tags(moltype=moltype, parameter_libraries=parameter_libraries, tagsarg=val)
            ILR_restart_layer(**restart_dict)

        elif val.lower() in ("printtag", "printtags", "pt"):
            self.print_term("Following tags are present across all {moltype_longname} parameter libraries:".format(moltype_longname=moltype_longname), verbose=0)
            self.print_term(" ".join(["'"+tag+"'" for tag in all_tags_in_params]), verbose=0, spaces=1)
            self.print_term("", verbose=0)
            ILR_restart_layer(**restart_dict)
        
        elif val.lower().startswith((moltype+":", moltype+"s:")):
            if len(val.split(":")) == 3:
                parameter_library, mol_name = val.split(":")[1:]
                if parameter_library in parameter_libraries:
                    if mol_name in mol_dict[parameter_library].keys():
                        self.print_mol_data(moltype, parameter_library, mol_name)
                    else:
                        self.print_term("{moltype_longname} name not found in parameter library.".format(moltype_longname=moltype_longname).capitalize(), "You must specify an existing {moltype_longname} in the parameter library: '{val}'.".format(moltype_longname=moltype_longname, val=val), verbose=0)
                        self.print_term("", verbose=0)
                else:
                    self.print_term("Parameter library not found. You must specify an existing parameter library: '{val}'.".format(val=val), verbose=0)
                    self.print_term("", verbose=0)
            else:
                self.print_term("You must specify both a parameter library and a {moltype_longname} name: '{val}'.".format(moltype_longname=moltype_longname, val=val), verbose=0)
                self.print_term("", verbose=0)
            ILR_restart_layer(**restart_dict)
        
        else:
            self.ILR_invalid_answer(val)
            ILR_restart_layer(**restart_dict)

    def ILR_layer2_lipid_solvent_ions(self, moltype, parameter_library):
        ILR_restart_layer = self.ILR_layer2_lipid_solvent_ions
        restart_dict={"moltype": moltype, "parameter_library": parameter_library}
        self.print_term("-"*self.terminalupdate_string_length, verbose=0)

        if moltype == "lipid":
            mol_dict = self.lipid_dict
            moltype_longname = "lipid"
            metadata_dict = self.lipid_metadata

        elif moltype == "solvent":
            mol_dict = self.solvent_dict
            moltype_longname = "solvent"
            metadata_dict = self.solvent_metadata
        
        elif moltype == "pos_ion":
            mol_dict = self.pos_ion_dict
            moltype_longname = "positive ion"
            metadata_dict = self.pos_ion_metadata
        
        elif moltype == "neg_ion":
            mol_dict = self.neg_ion_dict
            moltype_longname = "negative ion"
            metadata_dict = self.neg_ion_metadata
        
        else:
            assert False, "Invalid 'moltype' recieved by 'ILR_layer1': {moltype}".format(moltype=moltype)

        metadata_str = ""
        if parameter_library in metadata_dict.keys():
            metadata_str_list = []
            metadata_str_list.append("")
            metadata_str_list.append("This parameter library contains the following metadata:")
            for descriptor, lines in metadata_dict[parameter_library].items():
                descriptor = str(descriptor)
                if descriptor[-1] != ":":
                    descriptor = descriptor + ":"
                metadata_str_list.append("    " + descriptor)
                if type(lines) not in [list, tuple]:
                    lines = [lines]
                for line in lines:
                    metadata_str_list.append("    " + "    " + line)
            metadata_str_list.append("")
            metadata_str = "\n".join(metadata_str_list)

        tags_in_parameter_library = sorted(list(set([
            tag
            for mol in mol_dict[parameter_library].values()
            for tag in mol.tags
        ])))

        mol_names_in_parameter_library = sorted([mol for mol in mol_dict[parameter_library].keys()])

        longest_string = "Print names of all {moltype_longname}s in the parameter library:".format(moltype_longname=moltype_longname)
        longest_string_len = len(longest_string) + 1

        val = self.print_term(
            "\n".join([
                "You are in the {moltype_longname} parameter library '{parameter_library}':".format(moltype_longname=moltype_longname, parameter_library=parameter_library),
                metadata_str,
                "What would you like to do? Your options are shown below:",
                "    "+"Quit:".ljust(longest_string_len)                                                                                                      + "'q' or 'quit'",
                "    "+"Return to previous question:".ljust(longest_string_len)                                                                               + "'r' or 'return'",
                "    "+"Print names of {moltype_longname} with given tag(s):".format(moltype_longname=moltype_longname).ljust(longest_string_len),
                "    "+"    "+"- Only this specific tag:".ljust(longest_string_len)                                                                           + "'tag:[tag1]' (only 1 tag)",
                "    "+"    "+"- With any of these tags:".ljust(longest_string_len)                                                                           + "'any:[tag1]:[tag2]:[etc.]' (any number of tags)",
                "    "+"    "+"- With all of these tags:".ljust(longest_string_len)                                                                           + "'all:[tag1]:[tag2]:[etc.]' (any number of tags)",
                "    "+"    "+"- Combining the above:".ljust(longest_string_len)                                                                              + "'tag:[tag1]:any:[tag2]:[etc.]:all:[tag3]:[etc.]'",
                "    "+"Print all tags used in the parameter library:".ljust(longest_string_len)                                                              + "'printtag(s)' or 'pt'",
                "    "+"Print names of all {moltype_longname}s in the parameter library:".format(moltype_longname=moltype_longname).ljust(longest_string_len) + "'print{moltype}(s)' or 'p{first_letter}'".format(moltype=moltype, first_letter=moltype[0]),
                "    "+"Examine a specific {moltype_longname} in the parameter library:".format(moltype_longname=moltype_longname).ljust(longest_string_len)  + "'{moltype}(s):[{moltype_longname} name]' or '[{moltype_longname} name]'".format(moltype=moltype, moltype_longname=moltype_longname),
                "",
            ]),
            inp=True
        )
        self.print_term("", verbose=0)
        val = val.lstrip(" ").rstrip(" ")
        
        if val.lower() in ["q", "quit"]:
            pass
        
        elif val in ["r", "return"]:
            self.ILR_layer1_lipid_solvent_ions(moltype)
        
        elif val.lower() in ("print{moltype}".format(moltype=moltype), "print{moltype}s".format(moltype=moltype), "p{first_letter}".format(first_letter=moltype[0])):
            self.print_term("Following {moltype_longname}s are present in the '{parameter_library}' parameter library:".format(moltype_longname=moltype_longname, parameter_library=parameter_library), verbose=0)
            self.print_term(" ".join(mol_names_in_parameter_library), verbose=0, spaces=1)
            self.print_term("", verbose=0)
            ILR_restart_layer(**restart_dict)
        
        elif val.lower().startswith(("tag:", "any:" ,"all:")):
            self.print_mols_with_tags(moltype=moltype, parameter_libraries=parameter_libraries, tagsarg=val)
            ILR_restart_layer(**restart_dict)
        
        elif val.lower() in ("printtag", "printtags", "pt"):
            self.print_term("Following tags are present in the '{parameter_library}' parameter library:".format(parameter_library=parameter_library), verbose=0)
            self.print_term("    "+" ".join(["'"+tag+"'" for tag in tags_in_parameter_library]), verbose=0)
            self.print_term("", verbose=0)
            ILR_restart_layer(**restart_dict)
        
        elif val.startswith((moltype+":", moltype+"s:")) or val in mol_names_in_parameter_library:
            mol_name = False
            if val.startswith((moltype+":", moltype+"s:")):
                if len(val.split(":")) == 2:
                    mol_name = val.split(":")[1]
                else:
                    self.print_term("You must specify exactly one {moltype_longname} name: '{val}'.".format(moltype_longname=moltype_longname, val=val), verbose=0)
                    self.print_term("", verbose=0)
            else:
                mol_name = val
            
            if mol_name is not False:
                if mol_name in mol_dict[parameter_library].keys():
                    self.print_mol_data(moltype, parameter_library, mol_name)
                else:
                    self.print_term("{moltype_longname} name not found in parameter library.".format(moltype_longname=moltype_longname).capitalize(), "You must specify a {moltype_longname} existing in the parameter library: '{val}'.".format(moltype_longname=moltype_longname, val=val), verbose=0)
                    self.print_term("", verbose=0)
            ILR_restart_layer(**restart_dict)
        
        else:
            self.ILR_invalid_answer(val)
            ILR_restart_layer(**restart_dict)

    def ILR_layer1_protein(self):
        ILR_restart_layer = self.ILR_layer1_protein
        self.print_term("-"*self.terminalupdate_string_length, verbose=0)

        parameter_libraries = sorted(list(self.prot_defs.keys()))

        longest_string = "Print all residue names in a parameter library:"
        longest_string_len = len(longest_string) + 1

        val = self.print_term(
            "\n".join([
                "You are in the general protein residue charge library. The names of all available protein charge parameter libraries are shown below:",
                "    "+" ".join(["'"+parameter_library+"'" for parameter_library in parameter_libraries]),
                "",
                "What would you like to do? Your options are shown below:",
                "    "+"Quit:".ljust(longest_string_len)                                           + "'q' or 'quit'",
                "    "+"Return to previous question:".ljust(longest_string_len)                    + "'r' or 'return'",
                "    "+"Examine parameter library:".ljust(longest_string_len)                      + "'param(s):[parameter library name]' or '[parameter library name]'",
                "    "+"Print all residue names in a parameter library:".ljust(longest_string_len) + "'printres(idue)(s):[parameter library name]' or 'pr:[parameter library name]'",
                "    "+"Examine residue from given parameter library:".ljust(longest_string_len)   + "'res(idue)(s):[parameter library name]:[residue name]'",
                "",
            ]),
            inp=True
        )
        self.print_term("", verbose=0)
        val = val.lstrip(" ").rstrip(" ")

        if val.lower() in ["q", "quit"]:
            pass
        
        elif val in ["r", "return"]:
            self.ILR_layer0_main()
        
        elif val.lower().startswith(("param:", "params:")) or val in parameter_libraries:
            parameter_library = False
            if val.lower().startswith(("param:", "params:")):
                if len(val.split(":")) == 2:
                    parameter_library = val.split(":")[1]
                else:
                    self.print_term("You must specify exactly one parameter library: '{val}'.".format(val=val), verbose=0)
                    self.print_term("", verbose=0)
            else:
                parameter_library = val
            
            if parameter_library is not False:
                if parameter_library in parameter_libraries:
                    self.ILR_layer2_protein(parameter_library=parameter_library)
                else:
                    self.print_term("Parameter library not found. You must specify an existing parameter library: '{val}'.".format(val=val), verbose=0)
                    self.print_term("", verbose=0)
                    ILR_restart_layer()
            else:
                ILR_restart_layer()
        
        elif val.lower().startswith(("printres:", "printress:", "printresidue:", "printresidues:", "pr:")):
            if len(val.split(":")) == 2:
                parameter_library = val.split(":")[1]
                if parameter_library in parameter_libraries:
                    self.print_term("Following residues are present in the '{parameter_library}' parameter library.".format(parameter_library=parameter_library), verbose=0)
                    self.print_term(" ".join(list(self.prot_defs[parameter_library]["charges"].keys())), verbose=0, spaces=1)
                    self.print_term("", verbose=0)
                else:
                    self.print_term("Parameter library not found. You must specify an existing parameter library: '{val}'.".format(val=val), verbose=0)
                    self.print_term("", verbose=0)
            else:
                self.print_term("You must specify exactly one parameter library: '{val}'.".format(val=val), verbose=0)
                self.print_term("", verbose=0)
            ILR_restart_layer()

        elif val.lower().startswith(("res:", "ress:", "residue:", "residues:")):
            if len(val.split(":")) == 3:
                parameter_library, res_name = val.split(":")[1:]
                if parameter_library in parameter_libraries:
                    if res_name in self.prot_defs[parameter_library]["charges"].keys():
                        self.print_prot_charge_data(parameter_library, res_name)
                    else:
                        self.print_term("Protein residue name not found in parameter library. You must specify an existing protein residue in the parameter library: '{val}'.".format(val=val), verbose=0)
                        self.print_term("", verbose=0)
                else:
                    self.print_term("Parameter library not found. You must specify an existing parameter library: '{val}'.".format(val=val), verbose=0)
                    self.print_term("", verbose=0)
            else:
                self.print_term("You must specify both a parameter library and a protein residue name: '{val}'.".format(val=val), verbose=0)
                self.print_term("", verbose=0)
            ILR_restart_layer()
        
        else:
            self.ILR_invalid_answer(val)
            ILR_restart_layer()

    def ILR_layer2_protein(self, parameter_library):
        ILR_restart_layer = self.ILR_layer2_protein
        restart_dict={"parameter_library": parameter_library}
        self.print_term("-"*self.terminalupdate_string_length, verbose=0)

        metadata_dict = self.prot_metadata
        metadata_str = ""
        if parameter_library in metadata_dict.keys():
            metadata_str_list = []
            metadata_str_list.append("")
            metadata_str_list.append("This parameter library contains the following metadata:")
            for descriptor, lines in metadata_dict[parameter_library].items():
                descriptor = str(descriptor)
                if descriptor[-1] != ":":
                    descriptor = descriptor + ":"
                metadata_str_list.append("    " + descriptor)
                if type(lines) not in [list, tuple]:
                    lines = [lines]
                for line in lines:
                    metadata_str_list.append("    " + "    " + line)
            metadata_str_list.append("")
            metadata_str = "\n".join(metadata_str_list)

        res_names_in_parameter_library = sorted([res for res in self.prot_defs[parameter_library]["charges"].keys()])

        longest_string = "Print names of all protein residues in the parameter library:"
        longest_string_len = len(longest_string) + 1

        val = self.print_term(
            "\n".join([
                "You are in the protein residue charge parameter library '{parameter_library}':".format(parameter_library=parameter_library),
                metadata_str,
                "What would you like to do? Your options are shown below:",
                "    "+"Quit:".ljust(longest_string_len)                                                         + "'q' or 'quit'",
                "    "+"Return to previous question:".ljust(longest_string_len)                                  + "'r' or 'return'",
                "    "+"Print names of all protein residues in the parameter library:".ljust(longest_string_len) + "'printres(idues)(s)' or 'pr'",
                "    "+"Examine a specific protein residue in the parameter library:".ljust(longest_string_len)  + "'res(idues)(s):[protein residue name]' or '[protein residue name]'",
                "",
            ]),
            inp=True
        )
        self.print_term("", verbose=0)
        val = val.lstrip(" ").rstrip(" ")
        
        if val.lower() in ["q", "quit"]:
            pass
        
        elif val in ["r", "return"]:
            self.ILR_layer1_protein()
        
        elif val.lower() in ("printres", "printress", "printresidue", "printresidues"):
            self.print_term("Following residues are present in the '{parameter_library}' parameter library.".format(parameter_library=parameter_library), verbose=0)
            self.print_term(" ".join(res_names_in_parameter_library), verbose=0, spaces=1)
            self.print_term("", verbose=0)
            ILR_restart_layer(**restart_dict)
        
        elif val.startswith(("res:", "ress:", "residue:", "residues:")) or val in res_names_in_parameter_library:
            res_name = False
            if val.startswith(("res:", "residue:", "residues:")):
                if len(val.split(":")) == 2:
                    res_name = val.split(":")[1]
                else:
                    self.print_term("You must specify exactly one protein residue name: '{val}'.".format(val=val), verbose=0)
                    self.print_term("", verbose=0)
            else:
                res_name = val
            
            if res_name is not False:
                if res_name in self.prot_defs[parameter_library]["charges"].keys():
                    self.print_prot_charge_data(parameter_library, res_name)
                else:
                    self.print_term("Protein residue name not found in parameter library. You must specify a protein residue existing in the parameter library: '{val}'.".format(val=val), verbose=0)
                    self.print_term("", verbose=0)
            ILR_restart_layer(**restart_dict)
        
        else:
            self.ILR_invalid_answer(val)
            ILR_restart_layer(**restart_dict)

    def ILR_layer1_MFB(self):
        ILR_restart_layer = self.ILR_layer1_MFB
        self.print_term("-"*self.terminalupdate_string_length, verbose=0)

        moltypes = sorted(list(self.fragment_defs.keys()))

        longest_string = "Examine molecule type (moltype):"
        longest_string_len = len(longest_string) + 1

        val = self.print_term(
            "\n".join([
                "You are in the general molecule fragment builder library. The names of all available fragment molecule types (moltype) are shown below:",
                "    "+" ".join(["'"+moltype+"'" for moltype in moltypes]),
                "",
                "What would you like to do? Your options are shown below:",
                "    "+"Quit:".ljust(longest_string_len)                                              + "'q' or 'quit'",
                "    "+"Return to previous question:".ljust(longest_string_len)                       + "'r' or 'return'",
                "    "+"Examine molecule type (moltype):".ljust(longest_string_len)                   + "'moltype(s):[moltype]' or '[moltype]'",
                "",
            ]),
            inp=True
        )
        self.print_term("", verbose=0)
        val = val.lstrip(" ").rstrip(" ")

        if val.lower() in ["q", "quit"]:
            pass
        
        elif val in ["r", "return"]:
            self.ILR_layer0_main()
        
        elif val.lower().startswith(("moltype:", "moltypes:")) or val in moltypes:
            moltype = False
            if val.lower().startswith(("moltype:", "moltypes:")):
                if len(val.split(":")) == 2:
                    moltype = val.split(":")[1]
                else:
                    self.print_term("You must specify exactly one molecule type (moltype): '{val}'.".format(val=val), verbose=0)
                    self.print_term("", verbose=0)
            else:
                moltype = val
            
            if moltype is not False:
                if moltype in moltypes:
                    self.ILR_layer2_MFB(moltype=moltype)
                else:
                    self.print_term("Molecule type (moltype) not found. You must specify an existing molecule type: '{val}'.".format(val=val), verbose=0)
                    self.print_term("", verbose=0)
                    ILR_restart_layer()
            else:
                ILR_restart_layer()
        
        else:
            self.ILR_invalid_answer(val)
            ILR_restart_layer()

    def ILR_layer2_MFB(self, moltype):
        ILR_restart_layer = self.ILR_layer2_MFB
        restart_dict={"moltype": moltype}
        self.print_term("-"*self.terminalupdate_string_length, verbose=0)

        metadata_dict = self.fragment_metadata
        metadata_str = ""
        if moltype in metadata_dict.keys():
            metadata_str_list = []
            metadata_str_list.append("")
            metadata_str_list.append("This molecule type (moltype) contains the following metadata:")
            for descriptor, lines in metadata_dict[moltype].items():
                descriptor = str(descriptor)
                if descriptor[-1] != ":":
                    descriptor = descriptor + ":"
                metadata_str_list.append("    " + descriptor)
                if type(lines) not in [list, tuple]:
                    lines = [lines]
                for line in lines:
                    metadata_str_list.append("    " + "    " + line)
            metadata_str_list.append("")
            metadata_str = "\n".join(metadata_str_list)
        
        moltype_info_str_list = []
        if "accepted_parts" in self.fragment_defs[moltype]:
            moltype_info_str_list.append("The moltype '{moltype}' accepts the following part types:".format(moltype=moltype))
            moltype_info_str_list.append("    "+" ".join(["'"+val+"'" for val in self.fragment_defs[moltype]["accepted_parts"]]))
        else:
            moltype_info_str_list.append("The key 'accepted_parts' was not found in the dictionary for the moltype '{moltype}'. Please add it if you wish to use this moltype.".format(moltype=moltype))
        moltype_info_str_list.append("")
        
        if "order" in self.fragment_defs[moltype]:
            moltype_info_str_list.append("The moltype '{moltype}' will place the particles from the part types in the following order:".format(moltype=moltype))
            moltype_info_str_list.append("    "+" ".join(["'"+val+"'" for val in self.fragment_defs[moltype]["order"]]))
        else:
            moltype_info_str_list.append("The key 'order' was not found in the dictionary for the moltype '{moltype}'. Please add it if you wish to use this moltype.".format(moltype=moltype))
        moltype_info_str_list.append("")
        
        if "parts" in self.fragment_defs[moltype]:
            moltype_info_str_list.append("The moltype '{moltype}' contains the following parts for the given part types:".format(moltype=moltype))
            for parttype, parts in self.fragment_defs[moltype]["parts"].items():
                parts_str_list = []
                parts_str_list.append(parttype + ":")
                if "function" in parts.keys():
                    parts_str_list.append("'Part type built via function'")
                else:
                    for part in parts:
                        parts_str_list.append("'{part}'".format(part=part))
                moltype_info_str_list.append("    "+" ".join(parts_str_list))

        else:
            moltype_info_str_list.append("The key 'parts' was not found in the dictionary for the moltype '{moltype}'. Please add it if you wish to use this moltype.".format(moltype=moltype))
        moltype_info_str_list.append("")

        if "default_parts" in self.fragment_defs[moltype]:
            moltype_info_str_list.append("The moltype '{moltype}' contains the following default parts:".format(moltype=moltype))
            for parttype, part in self.fragment_defs[moltype]["default_parts"].items():
                moltype_info_str_list.append("    "+"{parttype}: {part}".format(parttype=parttype, part=part))
            moltype_info_str_list.append("")

        if "ascii" in self.fragment_defs[moltype]:
            moltype_info_str_list.append("The moltype '{moltype}' has the following ascii art detailing its structure:".format(moltype=moltype))
            for line in self.fragment_defs[moltype]["ascii"]:
                moltype_info_str_list.append("    "+line)
            moltype_info_str_list.append("")

        moltype_info_str = "\n".join(moltype_info_str_list)

        longest_string = "Print details and possible parts for a for a specific part type:"
        longest_string_len = len(longest_string) + 1

        val = self.print_term(
            "\n".join([
                "You are in the library entry for the molecule fragment builder molecule type (moltype) '{moltype}':".format(moltype=moltype),
                metadata_str,
                moltype_info_str,
                "What would you like to do? Your options are shown below:",
                "    "+"Quit:".ljust(longest_string_len)                                                            + "'q' or 'quit'",
                "    "+"Return to previous question:".ljust(longest_string_len)                                     + "'r' or 'return'",
                "    "+"Print details and possible parts for a for a specific part type:".ljust(longest_string_len) + "'printpart(s):[parttype]' / 'pp(s):[parttype]' or [parttype]",
                "",
            ]),
            inp=True
        )
        self.print_term("", verbose=0)
        val = val.lstrip(" ").rstrip(" ")
        
        if val.lower() in ["q", "quit"]:
            pass
        
        elif val in ["r", "return"]:
            self.ILR_layer1_MFB()
        
        elif val.startswith(("printpart:", "pp:", "printparts:", "pps:")) or val in self.fragment_defs[moltype]["parts"].keys():
            parttype = val
            if val.startswith(("printpart:", "pp:", "printparts:", "pps:")) and len(val.split(":")) > 1:
                parttype = val.split(":")[1]
            if parttype in self.fragment_defs[moltype]["parts"].keys():
                if "function" in self.fragment_defs[moltype]["parts"][parttype].keys():
                    self.print_term("    "+"This part is built via a function called '{function}' - Details for the building function for part type '{parttype}':".format(function=self.fragment_defs[moltype]["parts"][parttype]["function"].__name__, parttype=parttype), verbose=0)
                    if "kwargs" in self.fragment_defs[moltype]["parts"][parttype].keys():
                        self.print_term("    "+"    "+"This function uses the following kwargs:", verbose=0)
                        self.recursive_dict_printer(iterable=self.fragment_defs[moltype]["parts"][parttype]["kwargs"], current_spacing="    "+"    "+"    ", spacing="    ")
                else:
                    self.print_term("    "+"This part is selected from a list of prebuilt parts - Available parts for part type '{parttype}':".format(parttype=parttype), verbose=0)
                    for part_name, part_vals in self.fragment_defs[moltype]["parts"][parttype].items():
                        self.print_term("    "+"    "+"Part '{part_name}':".format(part_name=part_name), verbose=0)
                        self.print_term("    "+"    "+"    "+"It contains the following beads", verbose=0)
                        for bead in part_vals["beads"]:
                            self.print_term("    "+"    "+"    "+"    "+"'name': {name}, 'charge': {charge}, 'x': {x}, 'y': {y}, 'z': {z}, 'resname': {resname}, 'resnr': {resnr}".format(name=bead["name"], charge=bead["charge"], x=bead["x"], y=bead["y"], z=bead["z"], resname=bead["resname"], resnr=bead["resnr"]), verbose=0)
                        if "join_to" in part_vals.keys():
                            self.print_term("    "+"    "+"    "+"It joins to the following part at the given position", verbose=0)
                            self.print_term("    "+"    "+"    "+"    "+"{join_to_part}: {join_to_position}".format(join_to_part=part_vals["join_to"][0], join_to_position=part_vals["join_to"][0]), verbose=0)
                        if "join_from" in part_vals.keys():
                            self.print_term("    "+"    "+"    "+"It is joined to by the following part at the given position", verbose=0)
                            for frompart, position in part_vals["join_from"].items():
                                self.print_term("    "+"    "+"    "+"    "+"'join_from': {join_from}".format(join_from=part_vals["join_from"]), verbose=0)
                self.print_term("", verbose=0)

            else:
                self.print_term("Part type '{parttype}' not found. You must specify a valid part type. Valid part types:".format(parttype=parttype), verbose=0)
                self.print_term("    "+" ".join(["'"+parttype+"'" for parttype in self.fragment_defs[moltype]["parts"].keys()]), verbose=0)
                self.print_term("", verbose=0)

            ILR_restart_layer(**restart_dict)
        
        else:
            self.ILR_invalid_answer(val)
            ILR_restart_layer(**restart_dict)
