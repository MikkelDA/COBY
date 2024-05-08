'''
    If you want to run this program from the command line,
    then the command must involve the -m flag before the program.
    ex:
    python -m CPPM_calculator #Note that '.py' must not be included
    Flags should follow the above.
    ex:
    python -m CPPM_calculator -N normfile.txt
'''
import argparse
import sys
import ast

### Imports the package part
from COBY.__init__ import COBY

#####################################################################
########################## HERE BE PARSING ##########################
#####################################################################

if __name__ == "__main__":
    ### Custom Action classes to check if arguments have been given.
    given_arguments = set()

    class IsStored_ActionStore(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            given_arguments.add(self.dest)
            setattr(namespace, self.dest + '_set', True)
            setattr(namespace, self.dest, values)

    class IsStored_ActionAppend(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            given_arguments.add(self.dest)
            setattr(namespace, self.dest + '_set', True)
            items = getattr(namespace, self.dest, None)
            items = argparse._copy_items(items)
            items.append(values)
            setattr(namespace, self.dest, items)

    class IsStored_ActionExtend(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            given_arguments.add(self.dest)
            setattr(namespace, self.dest + '_set', True)
            items = getattr(namespace, self.dest, None)
            items = argparse._copy_items(items)
            items.extend(values)
            setattr(namespace, self.dest, items)

    ### Does not work
    # class IsStored_ActionCount(argparse.Action):
    #     def __call__(self, parser, namespace, values, option_string=None):
    #         given_arguments.add(self.dest)
    #         setattr(namespace, self.dest + '_set', True)
    #         items = getattr(namespace, self.dest, None)
    #         items = argparse._copy_items(items)
    #         items += 1
    #         setattr(namespace, self.dest, items)

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        add_help = False,
    )

    parser.add_argument("--help", "-h", dest = "help", action=IsStored_ActionStore)

    #######################
    ### SYSTEM CREATION ###
    #######################
    ### Leaflet commands
    parser.add_argument("--membrane",  "-membrane",  "-memb", dest = "membrane_cmds", action=IsStored_ActionAppend, type=str, default = [], nargs="+")

    ### Protein commands
    parser.add_argument("--protein",   "-protein",   "-prot", dest = "protein_cmds", action=IsStored_ActionAppend, type=str, default = [], nargs="+")

    ### Solvent commands
    parser.add_argument("--solvation", "-solvation", "-solv", dest = "solvation_cmds", action=IsStored_ActionAppend, type=str, default = [], nargs="+")

    ### Solvent commands
    parser.add_argument("--flooding",  "-flooding",  "-flood", dest = "flooding_cmds", action=IsStored_ActionAppend, type=str, default = [], nargs="+")

    ###############################
    ### SPECIAL SYSTEM CREATION ###
    ###############################
    ### Stacked membrane commands
    parser.add_argument("--stacked_membranes", "-stack_memb", "-stacked_membranes", dest = "stacked_membranes_cmds", action=IsStored_ActionAppend, type=str, default = [], nargs="+")

    ############
    ### MISC ###
    ############
    ### Topology commands
    parser.add_argument("--itp_input", "-itp_input", "-itp_in", dest = "itp_input_cmds", action=IsStored_ActionAppend, type=str, default = [], nargs="+")

    ### Molecule structure import commands
    parser.add_argument("--molecule_import", "-molecule_import", dest = "molecule_import_cmds", action=IsStored_ActionAppend, type=str, default = [], nargs="+")

    ### Molecule definition import commands
    parser.add_argument("--molecule_definition", "-molecule_definition", dest = "molecule_definition_cmds", action=IsStored_ActionAppend, type=str, default = [], nargs="+")

    ### Plotting command. Developer feature to test different algorithms.
    parser.add_argument("--plot_grid", "-plot_grid", "-plot", dest = "plot_grid_cmd", action=IsStored_ActionStore)

    ### Pickle commands
    parser.add_argument("--pickle", "-pickle", dest = "pickle_cmd", action=IsStored_ActionStore)

    ### Whether to backup files if they would be overwritten
    parser.add_argument("--backup", "-backup", dest = "backup_cmd", action=IsStored_ActionStore)

    ### Random seed
    parser.add_argument("--randseed", "-randseed", "-rand", dest = "randseed_cmd", action=IsStored_ActionStore)

    ### System parameters
    parser.add_argument("--sys_params",   "-sys_params",   "-sysp", dest = "sys_params",   action=IsStored_ActionStore)
    parser.add_argument("--prot_params",  "-prot_params",  "-pp",   dest = "prot_params",  action=IsStored_ActionStore)
    parser.add_argument("--lipid_params", "-lipid_params", "-lp",   dest = "lipid_params", action=IsStored_ActionStore)
    parser.add_argument("--solv_params",  "-solv_params",  "-sp",   dest = "solv_params",  action=IsStored_ActionStore)

    ### System name
    parser.add_argument("--system_name", "-system_name", "-sn", dest = "system_name", action=IsStored_ActionStore)

    #########################
    ### BOX SIZE AND TYPE ###
    #########################
    ### pbc box size [nm]
    parser.add_argument("--box", "-box", dest = "pbc_box", action=IsStored_ActionExtend, type=str, default = [], nargs="+")
    parser.add_argument("--pbc", "-pbc", dest = "pbc_box", action=IsStored_ActionExtend, type=str, default = [], nargs="+")
    
    ### x/y/z size of box [nm]
    parser.add_argument("--x", "-x", dest = "pbcx", type=str, action=IsStored_ActionStore)
    parser.add_argument("--y", "-y", dest = "pbcy", type=str, action=IsStored_ActionStore)
    parser.add_argument("--z", "-z", dest = "pbcz", type=str, action=IsStored_ActionStore)

    ### pbc box type
    parser.add_argument("--box_type", "-box_type", dest = "pbc_box_type", type=str, action=IsStored_ActionStore)
    parser.add_argument("--pbc_type", "-pbc_type", dest = "pbc_box_type", type=str, action=IsStored_ActionStore)

    ### Manual unit cell designation
    parser.add_argument("--pdb_unitcell", "-pdb_unitcell", dest = "pdb_unitcell", action=IsStored_ActionExtend, type=str, default = [], nargs="+")
    parser.add_argument("--gro_unitcell", "-gro_unitcell", dest = "gro_unitcell", action=IsStored_ActionExtend, type=str, default = [], nargs="+")

    ####################
    ### OUTPUT FILES ###
    ####################
    ### Output pdb/gro/top/log file
    parser.add_argument("--out_all", "-out_all", "-o_all", dest = "out_all_file_name", action=IsStored_ActionStore)

    ### Output pdb/gro file
    parser.add_argument("--out_sys", "-out_sys", "-o_sys", dest = "out_sys_file_name", action=IsStored_ActionStore)
    ### Output pdb file
    parser.add_argument("--out_pdb", "-out_pdb", "-o_pdb", dest = "out_pdb_file_name", action=IsStored_ActionStore)
    ### Output gro file
    parser.add_argument("--out_gro", "-out_gro", "-o_gro", dest = "out_gro_file_name", action=IsStored_ActionStore)

    ### Output topology file
    parser.add_argument("--out_top", "-out_top", "-t",     dest = "out_top_file_name", action=IsStored_ActionStore)

    ### Log file
    parser.add_argument("--out_log", "-out_log", "-log",   dest = "out_log_file_name", action=IsStored_ActionStore)

    ################
    ### PRINTING ###
    ################
    ### Prints
    parser.add_argument("--print_quiet",    "-quiet",          dest = "quiet",    action=IsStored_ActionStore)
    parser.add_argument("--print_debug",    "-debug",          dest = "debug",    action=IsStored_ActionStore)
    parser.add_argument("--print_extra",    "-extra",          dest = "extra",    action=IsStored_ActionStore)
    parser.add_argument("--print_warnings", "-warn",           dest = "warnings", action=IsStored_ActionStore)
    parser.add_argument("--verbose",        "-verbose", "-v",  dest = "verbose",  action=IsStored_ActionStore)

    ### ### Parser for handling '-f' when importing module to Jupyter
    # parser.add_argument("-fff", "-f", dest = "debug_flag_for_ipython")
    ### unknown variable includes the weird -f flag that jupyter puts in so no need for added argument.
    ### Argument still needed otherwise jupyter will throw the following error:
    ### "ipykernel_launcher.py: error: ambiguous option: -f could match -flood, -flooding"

    args, unknown = parser.parse_known_args()

    ##############################
    ### HELP FOR THOSE IN NEED ###
    ##############################

    if __name__ == "__main__":
        if "help" in given_arguments or len(given_arguments) == 0:
            parser.print_help()
            sys.exit()

    parser_kwargs = {}

    parse_membrane_cmds            = [" ".join(i) for i in args.membrane_cmds]
    parse_protein_cmds             = [" ".join(i) for i in args.protein_cmds]
    parse_solvation_cmds           = [" ".join(i) for i in args.solvation_cmds]
    parse_flooding_cmds            = [" ".join(i) for i in args.flooding_cmds]
    parse_stacked_membranes_cmds   = [" ".join(i) for i in args.stacked_membranes_cmds]
    parse_itp_input_cmds           = [" ".join(i) for i in args.itp_input_cmds]
    parse_molecule_import_cmds     = [" ".join(i) for i in args.molecule_import_cmds]
    parse_molecule_definition_cmds = [" ".join(i) for i in args.molecule_definition_cmds]

    for COBY_cmd, parse, arg_name in [
        ("membrane",            parse_membrane_cmds,            "membrane_cmds"),
        ("protein",             parse_protein_cmds,             "protein_cmds"),
        ("solvation",           parse_solvation_cmds,           "solvation_cmds"),
        ("flooding",            parse_flooding_cmds,            "flooding_cmds"),
        ("stacked_membranes",   parse_stacked_membranes_cmds,   "stacked_membranes_cmds"),
        ("itp_input",           parse_itp_input_cmds,           "itp_input_cmds"),
        ("molecule_import",     parse_molecule_import_cmds,     "molecule_import_cmds"),
        ("molecule_definition", parse_molecule_definition_cmds, "molecule_definition_cmds"),
        
        ("plot_grid", args.plot_grid_cmd, "plot_grid_cmd"),
        ("pickle",    args.pickle_cmd,    "pickle_cmd"),
        ("backup",    args.backup_cmd,    "backup_cmd"),
        ("randseed",  args.randseed_cmd,  "randseed_cmd"),
        
        ("sys_params",   args.sys_params,   "sys_params"),
        ("prot_params",  args.prot_params,  "prot_params"),
        ("lipid_params", args.lipid_params, "lipid_params"),
        ("solv_params",  args.solv_params,  "solv_params"),
        
        ("box",          args.pbc_box,      "pbc_box"),
        ("x",            args.pbcx,         "pbcx"),
        ("y",            args.pbcy,         "pbcy"),
        ("z",            args.pbcz,         "pbcz"),
        ("box_type",     args.pbc_box_type, "pbc_box_type"),
        ("pdb_unitcell", args.pdb_unitcell, "pdb_unitcell"),
        ("gro_unitcell", args.gro_unitcell, "gro_unitcell"),
        
        ("out_all", args.out_all_file_name, "out_all_file_name"),
        ("out_sys", args.out_sys_file_name, "out_sys_file_name"),
        ("out_pdb", args.out_pdb_file_name, "out_pdb_file_name"),
        ("out_gro", args.out_gro_file_name, "out_gro_file_name"),
        ("out_top", args.out_top_file_name, "out_top_file_name"),
        ("out_log", args.out_log_file_name, "out_log_file_name"),
        
        ("sn", args.system_name, "system_name"),
        
        ("quiet",    args.quiet,        "quiet"),
        ("debug",    args.debug,        "debug"),
        ("extra",    args.extra,        "extra"),
        ("warn",     args.warnings,     "warnings"),
        ("verbose",  args.verbose,      "verbose"),
    ]:
        if arg_name in given_arguments:
            parser_kwargs[COBY_cmd] = parse

    if parser_kwargs:
        COBY(
            run = True,
            terminal_run_kwargs = parser_kwargs,
        )

    #####################################################################
    ########################## YOU HAVE PARSED ##########################
    #####################################################################
