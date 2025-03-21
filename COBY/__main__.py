'''
    If you want to run this program from the argument line,
    then the argument must involve the -m flag before the program.
    ex:
    python -m CPPM_calculator #Note that '.py' must not be included
    Flags should follow the above.
    ex:
    python -m CPPM_calculator -N normfile.txt
'''
import argparse
import sys

### Imports the package part
from COBY.__init__ import *
from COBY.version import __version__, version_changes_str

#####################################################################
########################## HERE BE PARSING ##########################
#####################################################################

if __name__ == "__main__":
    ### Custom Action classes to check if arguments have been given.
    given_arguments = set()

    args_for_COBY = []

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

    class IsStored_ActionVersion(argparse._VersionAction):
        def __call__(self, parser, namespace, values, option_string=None):
            version = self.version
            formatter = parser._get_formatter()
            formatter.add_text(version)
            parser._print_message(formatter.format_help(), sys.stdout)
            parser.exit()

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

    ### Initial parser argument to check which subprogram is being called
    ### COBY is the default program
    parser.add_argument("--program", "-program", dest = "program", action=IsStored_ActionStore, type=str, default = "COBY")

    parser.add_argument("--help", "-h", dest = "help", action=IsStored_ActionStore)

    args, unknown = parser.parse_known_args()

    assert args.program in ["COBY", "Crafter", "Library"], "\n".join([
        "Invalid program given to '--program'. The valid programs are:",
        "    "+"COBY: (default) Runs the main COBY system builder",
        "    "+"Crafter: Runs the individual molecule system builder",
        "    "+"Library: Runs the interactive parameter library explorer",
    ])

    #########################################################################
    ########################## UNIVERSAL ARGUMENTS ##########################
    #########################################################################

    ################
    ### PRINTING ###
    ################
    ### Prints
    parser.add_argument("--print_quiet",      "-quiet",         dest = "quiet",      action=IsStored_ActionStore)
    parser.add_argument("--print_debug",      "-debug",         dest = "debug",      action=IsStored_ActionStore)
    parser.add_argument("--print_debug_keys", "-debug_keys",    dest = "debug_keys", action=IsStored_ActionStore)
    parser.add_argument("--print_extra",      "-extra",         dest = "extra",      action=IsStored_ActionStore)
    parser.add_argument("--print_warnings",   "-warn",          dest = "warnings",   action=IsStored_ActionStore)
    parser.add_argument("--verbose",          "-verbose", "-v", dest = "verbose",    action=IsStored_ActionStore)

    ###############
    ### VERSION ###
    ###############
    ### Program exits if one of the following arguments are used
    parser.add_argument("--version","-version",action=IsStored_ActionVersion, version="COBY version:" + " " + __version__)
    parser.add_argument("--version_changes", "--changes", "--changelog", "-version_changes",  "-changes", "-changelog", action=IsStored_ActionVersion, version="Changes from previous version (Current version: {cv}):\n".format(cv=__version__) + version_changes_str)

    ############
    ### MISC ###
    ############
    ### Topology arguments
    parser.add_argument("--itp_input", "-itp_input", "-itp_in", dest = "itp_input_args", action=IsStored_ActionAppend, type=str, default = [], nargs="+")

    ### Molecule definition import arguments
    parser.add_argument("--import_library", "-import_library", dest = "import_library_args", action=IsStored_ActionAppend, type=str, default = [], nargs="+")

    ### Lipid fragment builder arguments
    parser.add_argument("--molecule_builder", "-molecule_builder", dest = "molecule_builder_args", action=IsStored_ActionAppend, type=str, default = [], nargs="+")

    ### Random seed
    parser.add_argument("--randseed", "-randseed", "-rand", dest = "randseed_arg", action=IsStored_ActionStore)

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
    ### Output cif file
    parser.add_argument("--out_cif", "-out_cif", "-o_cif", dest = "out_cif_file_name", action=IsStored_ActionStore)

    ### Log file
    parser.add_argument("--out_log", "-out_log", "-log",   dest = "out_log_file_name", action=IsStored_ActionStore)

    ### Pickle arguments
    parser.add_argument("--pickle", "-pickle", dest = "pickle_arg", action=IsStored_ActionStore)

    ### Whether to backup files if they would be overwritten
    parser.add_argument("--backup", "-backup", dest = "backup_arg", action=IsStored_ActionStore)

    args, unknown = parser.parse_known_args()

    parse_itp_input_args         = [" ".join(i) for i in args.itp_input_args]
    parse_import_library_args    = [" ".join(i) for i in args.import_library_args]
    parse_molecule_builder_args  = [" ".join(i) for i in args.molecule_builder_args]

    args_for_COBY.extend([
        ("itp_input",         parse_itp_input_args,         "itp_input_args"),
        ("import_library",    parse_import_library_args,    "import_library_args"),
        ("molecule_builder",  parse_molecule_builder_args,  "molecule_builder_args"),
        
        ("pickle",    args.pickle_arg,    "pickle_arg"),
        ("backup",    args.backup_arg,    "backup_arg"),
        ("randseed",  args.randseed_arg,  "randseed_arg"),
        
        ("out_all", args.out_all_file_name, "out_all_file_name"),
        ("out_sys", args.out_sys_file_name, "out_sys_file_name"),
        ("out_pdb", args.out_pdb_file_name, "out_pdb_file_name"),
        ("out_gro", args.out_gro_file_name, "out_gro_file_name"),
        ("out_cif", args.out_cif_file_name, "out_cif_file_name"),
        ("out_log", args.out_log_file_name, "out_log_file_name"),
        
        ("quiet",      args.quiet,      "quiet"),
        ("debug",      args.debug,      "debug"),
        ("debug_keys", args.debug_keys, "debug_keys"),
        ("extra",      args.extra,      "extra"),
        ("warn",       args.warnings,   "warnings"),
        ("verbose",    args.verbose,    "verbose"),
    ])

    ### ### Parser for handling '-f' when importing module to Jupyter
    # parser.add_argument("-fff", "-f", dest = "debug_flag_for_ipython")
    ### unknown variable includes the weird -f flag that jupyter puts in so no need for added argument.
    ### Argument still needed otherwise jupyter will throw the following error:
    ### "ipykernel_launcher.py: error: ambiguous option: -f could match -flood, -flooding"

    #############################################################################
    ########################## COBY-SPECIFIC ARGUMENTS ##########################
    #############################################################################
    if args.program in ["COBY", "Library"]:

        #######################
        ### SYSTEM CREATION ###
        #######################
        ### Leaflet arguments
        parser.add_argument("--membrane",  "-membrane",  "-memb", dest = "membrane_args", action=IsStored_ActionAppend, type=str, default = [], nargs="+")

        ### Protein arguments
        parser.add_argument("--protein",   "-protein",   "-prot", dest = "protein_args", action=IsStored_ActionAppend, type=str, default = [], nargs="+")

        ### Solvent arguments
        parser.add_argument("--solvation", "-solvation", "-solv", dest = "solvation_args", action=IsStored_ActionAppend, type=str, default = [], nargs="+")

        ### Solvent arguments
        parser.add_argument("--flooding",  "-flooding",  "-flood", dest = "flooding_args", action=IsStored_ActionAppend, type=str, default = [], nargs="+")

        ################################
        ### STACKED MEMBRANE SYSTEMS ###
        ################################
        ### Stacked membrane arguments
        parser.add_argument("--stacked_membranes", "-stack_memb", "-stacked_membranes", dest = "stacked_membranes_args", action=IsStored_ActionAppend, type=str, default = [], nargs="+")

        ############
        ### MISC ###
        ############
        ### Molecule structure import arguments
        parser.add_argument("--molecule_import", "-molecule_import", dest = "molecule_import_args", action=IsStored_ActionAppend, type=str, default = [], nargs="+")

        ### Plotting argument. Developer feature to test different algorithms.
        parser.add_argument("--plot_grid", "-plot_grid", "-plot", dest = "plot_grid_arg", action=IsStored_ActionStore)

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
        parser.add_argument("--box", "--pbc", "-box", "-pbc", dest = "pbc_box", action=IsStored_ActionExtend, type=str, default = [], nargs="+")
        
        ### x/y/z size of box [nm]
        parser.add_argument("--x", "-x", dest = "pbcx", type=str, action=IsStored_ActionStore)
        parser.add_argument("--y", "-y", dest = "pbcy", type=str, action=IsStored_ActionStore)
        parser.add_argument("--z", "-z", dest = "pbcz", type=str, action=IsStored_ActionStore)

        ### pbc box type
        parser.add_argument("--box_type", "--pbc_type", "-box_type", "-pbc_type", dest = "pbc_box_type", type=str, action=IsStored_ActionStore)

        ### Manual unit cell designation
        parser.add_argument("--pdb_unitcell", "-pdb_unitcell", dest = "pdb_unitcell", action=IsStored_ActionExtend, type=str, default = [], nargs="+")
        parser.add_argument("--gro_unitcell", "-gro_unitcell", dest = "gro_unitcell", action=IsStored_ActionExtend, type=str, default = [], nargs="+")

        ####################
        ### OUTPUT FILES ###
        ####################
        ### Output topology file
        parser.add_argument("--out_top", "-out_top", "-t",     dest = "out_top_file_name", action=IsStored_ActionStore)

        args, unknown = parser.parse_known_args()
        
        parse_membrane_args          = [" ".join(i) for i in args.membrane_args]
        parse_protein_args           = [" ".join(i) for i in args.protein_args]
        parse_solvation_args         = [" ".join(i) for i in args.solvation_args]
        parse_flooding_args          = [" ".join(i) for i in args.flooding_args]
        parse_stacked_membranes_args = [" ".join(i) for i in args.stacked_membranes_args]
        parse_molecule_import_args   = [" ".join(i) for i in args.molecule_import_args]

        args_for_COBY.extend([
            ("membrane",          parse_membrane_args,          "membrane_args"),
            ("protein",           parse_protein_args,           "protein_args"),
            ("solvation",         parse_solvation_args,         "solvation_args"),
            ("flooding",          parse_flooding_args,          "flooding_args"),
            ("stacked_membranes", parse_stacked_membranes_args, "stacked_membranes_args"),
            ("molecule_import",   parse_molecule_import_args,   "molecule_import_args"),
            
            ("plot_grid", args.plot_grid_arg, "plot_grid_arg"),
            
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
            
            ("out_top", args.out_top_file_name, "out_top_file_name"),
            
            ("sn", args.system_name, "system_name"),
        ])

    ################################################################################
    ########################## CRAFTER-SPECIFIC ARGUMENTS ##########################
    ################################################################################
    if args.program in ["Crafter", "Library"]:
        #######################
        ### SYSTEM CREATION ###
        #######################
        ### Molecule system argument
        parser.add_argument("--molecule", "-molecule", "-mol", dest = "molecule_args", action=IsStored_ActionAppend, type=str, default = [], nargs="+")

        args, unknown = parser.parse_known_args()
        
        parse_molecule_args = [" ".join(i) for i in args.molecule_args]

        args_for_COBY.extend([
            ("molecule", parse_molecule_args, "molecule_args"),
        ])

    parser_kwargs = {}
    for COBY_arg, parse, arg_name in args_for_COBY:
        if arg_name in given_arguments:
            parser_kwargs[COBY_arg] = parse

    ##############################
    ### HELP FOR THOSE IN NEED ###
    ##############################
    if __name__ == "__main__":
        if "help" in given_arguments or len(given_arguments) == 0:
            parser.print_help()
            sys.exit()

    ############################
    ### RUNNING THE PROGRAMS ###
    ############################
    programs_dict = {
        "COBY": COBY,
        "Crafter": Crafter,
        "Library": Library,
    }

    if parser_kwargs:
        programs_dict[args.program](
            run = True,
            terminal_run_kwargs = parser_kwargs,
        )

#####################################################################
########################## YOU HAVE PARSED ##########################
#####################################################################
