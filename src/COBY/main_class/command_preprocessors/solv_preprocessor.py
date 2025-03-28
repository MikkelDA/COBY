import ast
import copy

class solv_preprocessor:

    def solv_preprocessor(self):
        '''
        Preprocesses solvation commands for later ease of use
        '''
        if len(self.SOLVATIONS_cmds) != 0:
            self.print_term("Preprocessing solvation requests", verbose=2)
            for cmd_nr, solvate_cmd in enumerate(self.SOLVATIONS_cmds, 1):
                self.print_term("Starting solvation argument:", cmd_nr, spaces=1, verbose=2)
                ### Defaults
                solv_dict = {
                    ### ### Solvent
                    "solvent": {}, # empty
                    "solvent_preprocessing": [], # empty
                    "solv_molarity": 55.56, # [int] or [float] # mol/L (molar)
                    "solvfreevol": True,

                    ### ### Ions
                    "neg_ions": {}, # empty
                    "pos_ions": {}, # empty
                    "neg_ions_preprocessing": [], # empty
                    "pos_ions_preprocessing": [], # empty
                    "target_charge": 0, # int or float # Target charge for system
                    # "sys_charge": True, # [bool] # Don't think it is used anymore
                    ### starting mol/liter concentration of both negative and positive ions (each will have the concentration) 
                    "salt_molarity": 0.15, # int or float # mol/L
                    "ionsvol": "solv", # "box", "free", "solv"
                    "salt_method": "add",
                    
                    ### ### Flooding specific
                    "flooding": False,
                    
                    ### ### General
                    ### Solvent box center:
                    "xlength": self.pbcx / 10, # [nm] converted to [Å]
                    "ylength": self.pbcy / 10, # [nm] converted to [Å]
                    "zlength": self.pbcz / 10, # [nm] converted to [Å]
                    "center": [0, 0, 0], # [nm] converted to [Å]
                    
                    "charge": "top", # "top" or "lib"
                    
                    ### "count" and "solv_per_lipid" are mutually exclusive
                    ### "solv_per_lipid" takes priority if both are given
                    "solv_per_lipid": False, # Number of solv particles per lipid contained in solvent box
                    "solv_per_lipid_cutoff": 0.5, # Number of solv particles per lipid contained in solvent box
                    "count": False, # [bool] Uses specific molarity value as absolute number of molecules instead of molarity ratio. Will be rounded using int(val+0.5)
                    "kick": 0.264/4, # [nm] converted to [Å]
                    "mapping": True,# bool, Whether AA-to-CG mapping should be considered
                    "ratio_method": "coarse", # "coarse", "atomistic"
                    
                    "params": False, # False or str
                    "bead_radius": 0.264, # [nm] converted to [Å] # Used for volume calculations
                    "gridres": 0.264, # [nm] converted to [Å] # 2.64
                    
                    "buffer": 0.2, # [nm] converted to [Å]
                    "protein_extra_buffer": 0.2, # [nm] converted to [Å]
                    "lipid_extra_buffer": 0, # [nm] converted to [Å]
                    "solute_extra_buffer": 0, # [nm] converted to [Å]

                    ### Can only be defined by the stacked membrane preprocessor
                    "extra_volume_box": [], # Adds extra volume for volume calculations

                    ### Can only be defined by the stacked membrane preprocessor
                    "extra_insertion_box": [], # Adds extra insertion space (positions) for insertion algorithm

                    ### Sets whether the hydrophobic volume of leaflets should be solvated
                    ### True values must be paired with True values in a "membrane / leaflet" argument
                    "solvate_hydrophobic_volume": False,
                }
                
                ### Added "default" solvation command per user request
                if solvate_cmd.startswith("default"):
                    solvate_cmd = solvate_cmd.lstrip("default")
                    solvate_cmd = " ".join(["solv:W pos:NA neg:CL", solvate_cmd])

                sbox = {}
                sbox["xmin"], sbox["xmax"], sbox["ymin"], sbox["ymax"], sbox["zmin"], sbox["zmax"] = False, False, False, False, False, False
                sbox["cx"], sbox["cz"], sbox["cy"], sbox["xlength"], sbox["ylength"], sbox["zlength"] = False, False, False, False, False, False

                for cmd in solvate_cmd.split():
                    sub_cmd = cmd.split(":")
                    
                    ########################
                    ### SOLVENT SPECIFIC ###
                    ########################
                    ### Solvents
                    if sub_cmd[0].lower() in ["solvent", "solv", "solute"]:
                        if solv_dict["flooding"]:
                            assert sub_cmd[0] == "solute", "'solv'/'solvent' subcommand used during solvation. Please use 'solute' instead"
                        else:
                            assert sub_cmd[0] in ["solvent", "solv"], "'solute' subcommand used during flooding. Please use 'solv'/'solvent' instead"
                        solv_dict["solvent_preprocessing"].append(sub_cmd[1:])

                    ### Solvent concentration
                    elif sub_cmd[0].lower() == "solv_molarity":
                        solv_dict["solv_molarity"] = ast.literal_eval(sub_cmd[1])

                    ### Whether to consider free volume (excluding lipid and protein) or purely the box volume
                    elif sub_cmd[0].lower() == "solvfreevol":
                        solv_dict["solvfreevol"] = ast.literal_eval(sub_cmd[1])

                    ####################
                    ### ION SPECIFIC ###
                    ####################
                    ### Negative ions
                    elif sub_cmd[0].lower() == "neg":
                        assert not solv_dict["flooding"], "subcommand 'neg' may not be used inside a flooding command"
                        solv_dict["neg_ions_preprocessing"].append(sub_cmd[1:])

                    ### Positive ions
                    elif sub_cmd[0].lower() == "pos":
                        assert not solv_dict["flooding"], "subcommand 'pos' may not be used inside a flooding command"
                        solv_dict["pos_ions_preprocessing"].append(sub_cmd[1:])

                    ### Target charge [int/float] or False for no "neutralization"
                    elif sub_cmd[0].lower() == "target_charge":
                        solv_dict["target_charge"] = ast.literal_eval(sub_cmd[1])

                    ### Whether to consider system charge for charge calculations
                    elif sub_cmd[0].lower() == "sys_charge":
                        if type(sub_cmd[1]) == bool:
                            solv_dict["sys_charge"] = sub_cmd[1]
                        else:
                            solv_dict["sys_charge"] = ast.literal_eval(sub_cmd[1])

                    ### Salt concentration
                    elif sub_cmd[0].lower() in ["salt_molarity"]:
                        solv_dict["salt_molarity"] = ast.literal_eval(sub_cmd[1])

                    ### Whether to consider free volume (excluding lipid and protein) or purely the box volume
                    elif sub_cmd[0].lower() == "ionsvol":
                        solv_dict["ionsvol"] = sub_cmd[1].lower()

                    ### Whether to consider system charge for charge calculations
                    elif sub_cmd[0].lower() == "salt_method":
                        valid_settings = ["add", "remove", "mean"]
                        assert sub_cmd[1] in ["add", "remove", "mean"], "\n".join([
                            "Subcommand 'salt_method' has been given with an invalid setting",
                            "Setting used was: " + str(sub_cmd[1]),
                            "Valid settings are: " + str(valid_settings),
                        ])
                        solv_dict["salt_method"] = sub_cmd[1]

                    ###############
                    ### GENERAL ###
                    ###############
                    ### Solvation box x-length
                    elif sub_cmd[0].lower() == "xlength":
                        # solv_dict["xlength"] = ast.literal_eval(sub_cmd[1])
                        sbox["xlength"] = ast.literal_eval(sub_cmd[1])
                    ### Solvation box y-length
                    elif sub_cmd[0].lower() == "ylength":
                        # solv_dict["ylength"] = ast.literal_eval(sub_cmd[1])
                        sbox["ylength"] = ast.literal_eval(sub_cmd[1])
                    ### Solvation box z-length
                    elif sub_cmd[0].lower() == "zlength":
                        # solv_dict["zlength"] = ast.literal_eval(sub_cmd[1])
                        sbox["zlength"] = ast.literal_eval(sub_cmd[1])
                        
                    ### Solvation box x center
                    elif sub_cmd[0].lower() == "cx":
                        # solv_dict["center"][0] = ast.literal_eval(sub_cmd[1])
                        sbox["cx"] = ast.literal_eval(sub_cmd[1])
                    ### Solvation box y center
                    elif sub_cmd[0].lower() == "cy":
                        # solv_dict["center"][1] = ast.literal_eval(sub_cmd[1])
                        sbox["cy"] = ast.literal_eval(sub_cmd[1])
                    ### Solvation box z center
                    elif sub_cmd[0].lower() == "cz":
                        # solv_dict["center"][2] = ast.literal_eval(sub_cmd[1])
                        sbox["cz"] = ast.literal_eval(sub_cmd[1])
                        
                    ### Solvation box center
                    elif sub_cmd[0].lower() == "center":
                        # solv_dict["center"] = [ast.literal_eval(val) for val in sub_cmd[1:]]
                        sbox["cx"], sbox["cy"], sbox["cz"] = [ast.literal_eval(val) for val in sub_cmd[1:]]
                    
                    ### Solvation box minima and maxima
                    elif sub_cmd[0].lower() == "xmin":
                        sbox["xmin"] = ast.literal_eval(sub_cmd[1])
                    elif sub_cmd[0].lower() == "xmax":
                        sbox["xmax"] = ast.literal_eval(sub_cmd[1])
                    elif sub_cmd[0].lower() == "ymin":
                        sbox["ymin"] = ast.literal_eval(sub_cmd[1])
                    elif sub_cmd[0].lower() == "ymax":
                        sbox["ymax"] = ast.literal_eval(sub_cmd[1])
                    elif sub_cmd[0].lower() == "zmin":
                        sbox["zmin"] = ast.literal_eval(sub_cmd[1])
                    elif sub_cmd[0].lower() == "zmax":
                        sbox["zmax"] = ast.literal_eval(sub_cmd[1])
                    
                    ### Whether charges should be obtained from topology
                    elif sub_cmd[0].lower()  == "charge":
                        assert sub_cmd[1] in ["topology", "top", "library", "lib"], "The 'charge' subcommand must be given either 'topology'/'top' or 'library'/'lib' as value. You have given: " + str(sub_cmd[1])
                        if sub_cmd[1] in ["topology", "top"]:
                            solv_dict["charge"] = "top"
                        elif sub_cmd[1] in ["library", "lib"]:
                            solv_dict["charge"] = "lib"
                        
                    ### Whether to use ratios as absolute number of molecules. True/False
                    elif sub_cmd[0].lower() == "solv_per_lipid":
                        solv_dict["solv_per_lipid"] = self.get_number_from_string(sub_cmd[1])
                        
                    ### Whether to use ratios as absolute number of molecules. [int/float]
                    elif sub_cmd[0].lower() == "solv_per_lipid_cutoff":
                        val = ast.literal_eval(sub_cmd[1])
                        assert 1 >= val >= 0, "The value for 'solv_per_lipid_cutoff' must be between 1 and 0. Your given value is '{cutoff}'.".format(cutoff=val)
                        solv_dict["solv_per_lipid_cutoff"] = val
                        
                    ### Whether to use ratios as absolute number of molecules. True/False
                    elif sub_cmd[0].lower() == "count":
                        solv_dict["count"] = ast.literal_eval(sub_cmd[1])

                    ### Radius of non-solvent beads for volume calculations
                    elif sub_cmd[0].lower() == "bead_radius":
                        solv_dict["bead_radius"] = ast.literal_eval(sub_cmd[1])

                    ### Grid resolution
                    elif sub_cmd[0].lower() == "gridres":
                        solv_dict["gridres"] = ast.literal_eval(sub_cmd[1])

                    ### Minimum wiggle room for beads
                    elif sub_cmd[0].lower() == "wr":
                        solv_dict["WR"] = ast.literal_eval(sub_cmd[1])

                    elif sub_cmd[0].lower() in ["kick"]:
                        ### Random kick to beads x/y/z positions [Å]
                        solv_dict["kick"] = ast.literal_eval(sub_cmd[1])
                    
                    ### Whether mapping should be used
                    elif sub_cmd[0].lower() == "mapping":
                        val = sub_cmd[1]
                        if type(val) == str:
                            val = ast.literal_eval(sub_cmd[1])
                        solv_dict["mapping"] = val
                    
                    ### Whether mapping should be used
                    elif sub_cmd[0].lower() == "ratio_method":
                        val = sub_cmd[1]
                        if val in ["coarse", "atomistic"]:
                            solv_dict["ratio_method"] = val
                    
                    ### Bead scaling distances
                    ### Bead distance scaling for z is 0.3 by default and 0.25 by default for x and y [multiplier]
                    elif sub_cmd[0].lower() in ["bdx", "bdy", "bdz"]:
                        leaf_dict[sub_cmd[0].lower()] = ast.literal_eval(sub_cmd[1])

                    ### Designates force field
                    elif sub_cmd[0].lower() == "params":
                        solv_dict["params"] = sub_cmd[1]

                    ### Gerenal buffer for solvent placement
                    elif sub_cmd[0].lower() == "buffer":
                        solv_dict["buffer"] = ast.literal_eval(sub_cmd[1])

                    ### Extra buffer given to proteins
                    elif sub_cmd[0].lower() == "protein_extra_buffer":
                        solv_dict["protein_extra_buffer"] = ast.literal_eval(sub_cmd[1])

                    ### Extra buffer given to lipids
                    elif sub_cmd[0].lower() == "lipid_extra_buffer":
                        solv_dict["lipid_extra_buffer"] = ast.literal_eval(sub_cmd[1])

                    ### Extra buffer given to prior solvents (solutes)
                    elif sub_cmd[0].lower() == "solute_extra_buffer":
                        solv_dict["solute_extra_buffer"] = ast.literal_eval(sub_cmd[1])

                    ### "extra_insertion_box" is used with stacked membrane arguments.
                    ### "extra_volume_box" is currently unsused. Previously used with stacked membrane arguments instead of "extra_insertion_box".
                    ### Most preprocessing is identical for the two so it is done jointly.
                    elif sub_cmd[0].lower() in ["extra_volume_box", "evb", "extra_insertion_box", "eib"]:
                        
                        if sub_cmd[0].lower() in ["extra_volume_box", "evb"]:
                            fn, abr = "extra_volume_box", "evb"
                        elif sub_cmd[0].lower() in ["extra_insertion_box", "eib"]:
                            fn, abr = "extra_insertion_box", "eib"
                        
                        sub_cmd_dict_given = {sub_cmd[i*2+1]: sub_cmd[i*2+2] for i in range(len(sub_cmd[1:]) // 2)}

                        for ax in ["x", "y", "z"]:
                            ### Following assert is a NAND (Not AND) gate as either one type or neither must be given 
                            assert not (all([i in sub_cmd_dict_given.keys() for i in [ax+"min", ax+"max"]]) and any([i in sub_cmd_dict_given.keys() for i in ["c"+ax, ax+"length"]])), "\n".join([
                                "Incorrect subsubarguments given to '{fn}' / '{abr}' x-dimensions must be specified by using either '{ax}min' and '{ax}max' OR 'c{ax}' or '{ax}length', not both".format(fn=fn, abr=abr, ax=ax),
                                "    "+"Your subargument:",
                                "    "+"    "+cmd,
                                "    "+"Correct example with '{ax}min' and '{ax}max' (subsubarguments for the other axes are not included):".format(ax=ax),
                                "    "+"    "+"{abr}:{ax}min:-2:{ax}max:8".format(abr=abr, ax=ax),
                                "    "+"Correct example with 'c{ax}' and '{ax}length' (subsubarguments for the other axes are not included):".format(ax=ax),
                                "    "+"    "+"{abr}:c{ax}:3:{ax}length:10".format(abr=abr, ax=ax),
                            ])
                            
                            ### Converting center/length to min/max
                            if "c"+ax in sub_cmd_dict_given.keys() or ax+"length" in sub_cmd_dict_given.keys():
                                if "c"+ax not in sub_cmd_dict_given.keys():
                                    cax = solv_dict["center"][("x", "y", "z").index(ax)]
                                if ax+"length" not in sub_cmd_dict_given.keys():
                                    axlength = solv_dict[ax+"length"]
                                sub_cmd_dict_given[ax+"min"] = cax - axlength/2
                                sub_cmd_dict_given[ax+"max"] = cax + axlength/2

                        sub_cmd_dict = {}
                        for axis_direction, default_val in [("xmin", -self.pbcx/2), ("xmax", self.pbcx/2), ("ymin", -self.pbcy/2), ("ymax", self.pbcy/2), ("zmin", -self.pbcz/2), ("zmax", self.pbcz/2)]:
                            if axis_direction not in sub_cmd:
                                sub_cmd_dict[axis_direction] = str(default_val/10)
                            else:
                                sub_cmd_dict[axis_direction] = sub_cmd_dict_given[axis_direction]
                        
                        if sub_cmd[0].lower() in ["extra_volume_box", "evb"]:
                            solv_dict["extra_volume_box"].append(tuple([ast.literal_eval(val) for val in sub_cmd_dict.values()]))
                        elif sub_cmd[0].lower() in ["extra_insertion_box", "eib"]:
                            solv_dict["extra_insertion_box"].append(tuple([ast.literal_eval(val) for val in sub_cmd_dict.values()]))

                    ### Used with stacked membrane arguments
                    elif sub_cmd[0].lower() == "solvate_hydrophobic_volume":
                        solv_dict["solvate_hydrophobic_volume"] = [ast.literal_eval(val) for val in sub_cmd[1:]]
                    
                    #########################
                    ### FLOODING SPECIFIC ###
                    #########################
                    ### Used to check if command is made as flooding or not
                    elif sub_cmd[0].lower() == "flooding":
                        solv_dict["flooding"] = ast.literal_eval(sub_cmd[1])
                    
                    ### Errors out if unknown subcommand used, and prints the subcommand to terminal
                    else:
                        assert False, "Unknown subcommand given to '-solvate'. The subcommand is: '" + str(cmd) + "'"
                
                for ax_i, ax in enumerate(["x", "y", "z"]):
                    ### NAND (Not AND) gate checking if any solvent box position and size settings have been given.
                    ### Values could be zero, so can't simply check if false-like, but must check if explicitly False
                    assert not ((sbox["c"+ax] is not False or sbox[ax+"length"] is not False) and (sbox[ax+"min"] is not False or sbox[ax+"max"] is not False)), "\n".join([
                        "Only one method for selecting the solvent box is allowed for each axis. You have used subarguments from both methods.",
                        "    "+"Available method 1 (center-length) subarguments: 'c{ax}' and '{ax}length'".format(ax=ax),
                        "    "+"    "+"Note that 'c{ax}' may be interpreted from the 'center' subargument if you used it.".format(ax=ax),
                        "    "+"Available method 2 (min-max) subarguments: '{ax}min' and '{ax}max'".format(ax=ax),
                        "    "+"Your subarguments (False, means that they were not given):",
                        "    "+"    "+"Method 1: '{cax}' and '{axlength}'".format(cax=sbox["c"+ax], axlength=sbox[ax+"length"]),
                        "    "+"    "+"Method 2: '{axmin}' and '{axmax}'".format(axmin=sbox[ax+"min"], axmax=sbox[ax+"max"]),
                    ])
                    ### If center-length method is used.
                    if sbox["c"+ax] is not False or sbox[ax+"length"] is not False:
                        if sbox["c"+ax] is not False:
                            solv_dict["center"][ax_i] = sbox["c"+ax]
                        if sbox[ax+"length"] is not False:
                            assert sbox[ax+"length"] > 0, "Axis length must be greater than zero. Your {ax}length: {axlength}".format(ax, sbox[ax+"length"])
                            solv_dict[ax+"length"] = sbox[ax+"length"]
                    
                    ### If min-max method is used.
                    elif sbox[ax+"min"] is not False or sbox[ax+"max"] is not False:
                        if sbox[ax+"min"] is not False:
                            axmin = sbox[ax+"min"]
                        else:
                            axmin = solv_dict["center"][ax_i] - solv_dict[ax+"length"]/2
                        if sbox[ax+"max"] is not False:
                            axmax = sbox[ax+"max"]
                        else:
                            axmax = solv_dict["center"][ax_i] + solv_dict[ax+"length"]/2
                        
                        assert axmax > axmin, "Axis maximum must be greater than axis minimum. Your {ax}min and {ax}max: {axmin} and {axmax} (note that system box values are used for {ax}min and {ax}max if not specified by you)".format(ax, ax, axmin, axmax, ax, ax)

                        solv_dict["center"][ax_i] = (axmax + axmin) / 2
                        solv_dict[ax+"length"] = axmax - axmin
                            
                    ### No values given for axis. Uses default values.
                    else:
                        pass
                    

                ### fixing a couple of values such that they are in angstrom
                ### 'extra_volume_box' and 'extra_insertion_box' values are fixed further down
                solv_dict.update({
                    "xlength": solv_dict["xlength"]*10,
                    "ylength": solv_dict["ylength"]*10,
                    "zlength": solv_dict["zlength"]*10,
                    "center": [val*10 for val in solv_dict["center"]],
                    
                    "kick":        solv_dict["kick"]*10,
                    "bead_radius": solv_dict["bead_radius"]*10,
                    "gridres":     solv_dict["gridres"]*10,
                    
                    "buffer":               solv_dict["buffer"]*10,
                    "protein_extra_buffer": solv_dict["protein_extra_buffer"]*10,
                    "lipid_extra_buffer":   solv_dict["lipid_extra_buffer"]*10,
                    "solute_extra_buffer":  solv_dict["solute_extra_buffer"]*10,
                })

                ### Ensures that no solvent boxes extends beyond the PBC in any direction
                def bounds_checker(c, length, pbc, err_str_ax, err_str_boxtype):
                    '''
                    Checks whether a box stretches outside the pbc.
                    '''
                    ### ### Positive direction, both upper and lower are calculated because they are used in the error messages.
                    upper = round(c + length/2, 3)
                    lower = round(c - length/2, 3)
                    if upper > round(pbc/2+0.001, 3): # Small addition to pbc to prevent activation by float errors
                        exceeding = upper - pbc/2
                        new_c = c - exceeding/2
                        new_length = length - exceeding
                        self.print_term("\n".join([
                            "{err_str_boxtype} exceeds {err_str_ax}-bounds of the pbc box in the positive {err_str_ax}-direction. Will correct to fit pbc.".format(err_str_boxtype=err_str_boxtype, err_str_ax=err_str_ax),

                            "    " + "Positive pbc {err_str_ax}-axis box bound: {val} nm".format(err_str_ax=err_str_ax, val=round(pbc/2/10, 4)),
                            "    " + "Positive {err_str_boxtype} {err_str_ax}-axis bound: {val} nm".format(err_str_boxtype=err_str_boxtype, err_str_ax=err_str_ax, val=round(upper/10, 4)),
                            "    " + "    " + "Exceeds pbc bounds by: {exceeding} nm".format(exceeding=round(exceeding/10, 4)),

                            "    " + "c{err_str_ax}:".format(err_str_ax=err_str_ax),
                            "    " + "    " + "Old: {val} nm".format(err_str_ax=err_str_ax, val=round(c/10, 4)),
                            "    " + "    " + "New: {val} nm".format(err_str_ax=err_str_ax, val=round(new_c/10, 4)),

                            "    " + "{err_str_ax}length:".format(err_str_ax=err_str_ax),
                            "    " + "    " + "Old: {val} nm".format(err_str_ax=err_str_ax, val=round(length/10, 4)),
                            "    " + "    " + "New: {val} nm".format(err_str_ax=err_str_ax, val=round(new_length/10, 4)),

                            "    " + "{err_str_ax}-axis bounds:".format(err_str_ax=err_str_ax),
                            "    " + "    " + "Old: {upper} nm to {lower} nm".format(err_str_ax=err_str_ax, upper=round(upper/10, 4), lower=round(lower/10, 4)),
                            "    " + "    " + "New: {upper} nm to {lower} nm".format(err_str_ax=err_str_ax, upper=round((new_c + new_length/2)/10, 4), lower=round((new_c - new_length/2)/10, 4)),
                        ]), warn=True)
                        c = new_c
                        length = new_length

                    ### ### Negative direction, both upper and lower are calculated because they are used in the error messages.
                    upper = round(c + length/2, 3)
                    lower = round(c - length/2, 3)
                    if lower < -round(pbc/2+0.001, 3): # Small addition to pbc to prevent activation by float errors
                        exceeding = lower + pbc/2
                        new_c = c - exceeding/2
                        new_length = length + exceeding
                        self.print_term("\n".join([
                            "{err_str_boxtype} exceeds {err_str_ax}-bounds of the pbc box in the negative {err_str_ax}-direction. Will correct to fit pbc.".format(err_str_boxtype=err_str_boxtype, err_str_ax=err_str_ax),

                            "    " + "Negative pbc {err_str_ax}-axis box boundary: {val} nm".format(err_str_ax=err_str_ax, val=round(-pbc/2/10, 4)),
                            "    " + "Negative {err_str_boxtype} {err_str_ax}-axis boundary: {val} nm".format(err_str_boxtype=err_str_boxtype, err_str_ax=err_str_ax, val=round(lower/10, 4)),
                            "    " + "    " + "Exceeds pbc bounds by: {exceeding} nm".format(exceeding=round(exceeding/10, 4)),

                            "    " + "c{err_str_ax}:".format(err_str_ax=err_str_ax),
                            "    " + "    " + "Old: {val} nm".format(err_str_ax=err_str_ax, val=round(c/10, 4)),
                            "    " + "    " + "New: {val} nm".format(err_str_ax=err_str_ax, val=round(new_c/10, 4)),

                            "    " + "{err_str_ax}length:".format(err_str_ax=err_str_ax),
                            "    " + "    " + "Old: {val} nm".format(err_str_ax=err_str_ax, val=round(length/10, 4)),
                            "    " + "    " + "New: {val} nm".format(err_str_ax=err_str_ax, val=round(new_length/10, 4)),

                            "    " + "{err_str_ax}-axis boundaries:".format(err_str_ax=err_str_ax),
                            "    " + "    " + "Old: {upper} nm to {lower} nm".format(err_str_ax=err_str_ax, upper=round(upper/10, 4), lower=round(lower/10, 4)),
                            "    " + "    " + "New: {upper} nm to {lower} nm".format(err_str_ax=err_str_ax, upper=round((new_c + new_length/2)/10, 4), lower=round((new_c - new_length/2)/10, 4)),

                        ]), warn=True)
                        c = new_c
                        length = new_length

                    return c, length
                    
                new_cx, new_xlength = bounds_checker(solv_dict["center"][0], solv_dict["xlength"], self.pbcx, "x", "Solvent box")
                new_cy, new_ylength = bounds_checker(solv_dict["center"][1], solv_dict["ylength"], self.pbcy, "y", "Solvent box")
                new_cz, new_zlength = bounds_checker(solv_dict["center"][2], solv_dict["zlength"], self.pbcz, "z", "Solvent box")
                solv_dict["center"] = (new_cx, new_cy, new_cz)
                solv_dict["xlength"] = new_xlength
                solv_dict["ylength"] = new_ylength
                solv_dict["zlength"] = new_zlength

                for key, err_str in [("extra_volume_box", "Extra volume box"), ("extra_insertion_box", "Extra insertion box")]:
                    for i, (cxmin_extra, cxmax_extra, cymin_extra, cymax_extra, czmin_extra, czmax_extra) in enumerate(solv_dict[key]):
                        cxmin_extra, cxmax_extra, cymin_extra, cymax_extra, czmin_extra, czmax_extra = cxmin_extra*10, cxmax_extra*10, cymin_extra*10, cymax_extra*10, czmin_extra*10, czmax_extra*10
                        cx_extra, cy_extra, cz_extra = (cxmin_extra + cxmax_extra)/2, (cymin_extra + cymax_extra)/2, (czmin_extra + czmax_extra)/2
                        xlength_extra = cxmax_extra - cxmin_extra
                        ylength_extra = cymax_extra - cymin_extra
                        zlength_extra = czmax_extra - czmin_extra
                        new_cx, new_xlength = bounds_checker(cx_extra, xlength_extra, self.pbcx, "x", err_str)
                        new_cy, new_ylength = bounds_checker(cy_extra, ylength_extra, self.pbcy, "y", err_str)
                        new_cz, new_zlength = bounds_checker(cz_extra, zlength_extra, self.pbcz, "z", err_str)
                        new_cxmin_extra = new_cx - new_xlength/2
                        new_cxmax_extra = new_cx + new_xlength/2
                        new_cymin_extra = new_cy - new_ylength/2
                        new_cymax_extra = new_cy + new_ylength/2
                        new_czmin_extra = new_cz - new_zlength/2
                        new_czmax_extra = new_cz + new_zlength/2
                        solv_dict[key][i] = (new_cxmin_extra, new_cxmax_extra, new_cymin_extra, new_cymax_extra, new_czmin_extra, new_czmax_extra)

                ### If solv_dict["solv_per_lipid"] has been set, then ignore count command
                if solv_dict["solv_per_lipid"] and solv_dict["count"]:
                    solv_dict["count"] = False

                ######################################
                ### SOLVENT/ION DATA INCORPORATION ###
                ######################################
                assert solv_dict["solvent_preprocessing"] != {}, "No solvent given to '-solvation' flag. Please specify at least one non-ionic solvent if you use it."

                params = solv_dict["params"] or self.solv_params or self.sys_params # (sys_params defaults to "default")
                
                solv_dict["solv_tot_ratio"] = 0
                solv_dict["pos_tot_ratio"]  = 0
                solv_dict["neg_tot_ratio"]  = 0

                types_to_be_processed = []
                if solv_dict["solvent_preprocessing"]:
                    types_to_be_processed.extend([("solvent", molecules) for molecules in solv_dict["solvent_preprocessing"]])
                
                ### Error out if positive OR negative ions have been set but not both. (XNOR gate)
                assert bool(solv_dict["pos_ions_preprocessing"]) == bool(solv_dict["neg_ions_preprocessing"]), "\n".join([
                    "You have given either positive or negative ions without specifying the other.",
                    "If you wish to neutralize using ions then you must provide both positive and negative ions.",
                    "Ions can be inserted without neutralization at a given concentration by using their automatically generated parameter libraries as shown in the examples below.",
                    "    " + "Positive ions: 'solvent:CL:params:pos_ions solv_molarity:0.5' (inserts CL ions at a concentration of 0.5 mol/L)",
                    "    " + "Negative ions: 'solvent:NA:params:neg_ions solv_molarity:0.5' (inserts NA ions at a concentration of 0.5 mol/L)",
                    "All positive and negative ions in the molecule definitions are added to the 'pos_ions' and 'neg_ions' parameter libraries automatically.",
                ])
                    
                if solv_dict["pos_ions_preprocessing"]:
                    types_to_be_processed.extend([("pos_ions", molecules) for molecules in solv_dict["pos_ions_preprocessing"]])
                if solv_dict["neg_ions_preprocessing"]:
                    types_to_be_processed.extend([("neg_ions", molecules) for molecules in solv_dict["neg_ions_preprocessing"]])
                
                ### ### Designated solvent, solute and ion preprocessing
                for solv_type, sub_cmd in types_to_be_processed:

                    name = False
                    ratio = 1
                    params = solv_dict["params"] or self.solv_params or self.sys_params
                    charge = solv_dict["charge"]

                    all_subcmd_names = ["params", "charge", "ratio", "name"]

                    i = 0
                    while i < len(sub_cmd):
                        if sub_cmd[i] == "params":
                            params = sub_cmd[i+1]
                            i += 2
                        elif sub_cmd[i] == "charge":
                            assert sub_cmd[i+1] in ["top", "topology", "lib", "library"], "Only 'top' and 'lib' are allowed values for solvent/solute/ion specific 'charge'" + "\n" + "Your subcommand: " + str(sub_cmd)
                            if sub_cmd[i+1] in ["top", "topology"]:
                                charge = "top"
                            elif sub_cmd[i+1] in ["lib", "library"]:
                                charge = "lib"
                            i += 2
                        elif sub_cmd[i] == "ratio":
                            ratio = ast.literal_eval(sub_cmd[i+1])
                            i += 2
                        elif sub_cmd[i] == "name":
                            name = sub_cmd[i+1]
                            i += 2
                        else: ### Assumes that it is the solvent/solute/ion name and potentially ratio
                            name = sub_cmd[i]
                            if len(sub_cmd[i:]) > 1 and sub_cmd[i+1] not in all_subcmd_names and self.get_number_from_string(sub_cmd[i+1]) is not False:
                                ratio = ast.literal_eval(sub_cmd[i+1])
                                i += 1
                            i += 1
                    
                    ### Checks if solvent/solute/ion name has been given and if it exists in the specified parameter library
                    assert name is not False, "A name has not been given for a solvent/solute/ion. Full solvent/solute/ion command: " + str(sub_cmd)

                    ### Adding name:dict combo to solvent dict
                    if solv_type == "solvent":
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
                        solv_dict[solv_type][name] = copy.deepcopy(self.solvent_dict[params][name])
                        solv_ratio_type = "solv_tot_ratio"
                        molarity = solv_dict["solv_molarity"]
                    
                    elif solv_type in ["pos_ions", "neg_ions"]:
                        assert params in self.pos_ion_dict.keys(), "\n".join([
                            "The ion parameter library '{params}' was not found".format(params=params),
                            "The available ion parameter libraries are:"
                            "    ", " ".join(list(self.pos_ion_dict.keys()))
                        ])
                        if solv_type == "pos_ions":
                            # assert "positive" in self.pos_ion_dict[params].keys(), "No positive ions found in the ion parameter library '{params}'".format(params=params)
                            assert name in self.pos_ion_dict[params].keys(), "\n".join([
                                "The positive ion '{name}' was not found in the ion parameter library '{params}'".format(name=name, params=params),
                                "The available positive ions in the ion parameter library '{params}' are:".format(params=params),
                                "    ", " ".join(list(self.pos_ion_dict[params].keys()))
                            ])
                            solv_dict[solv_type][name] = copy.deepcopy(self.pos_ion_dict[params][name])
                            solv_ratio_type = "pos_tot_ratio"
                            molarity = solv_dict["salt_molarity"]
                        
                        elif solv_type == "neg_ions":
                            # assert "negative" in self.neg_ion_dict[params].keys(), "No negative ions found in the ion parameter library '{params}'".format(params=params)
                            assert name in self.neg_ion_dict[params].keys(), "\n".join([
                                "The negative ion '{name}' was not found in the ion parameter library '{params}'".format(name=name, params=params),
                                "The available negative ions in the ion parameter library '{params}' are:".format(params=params),
                                "    ", " ".join(list(self.neg_ion_dict[params].keys()))
                            ])
                            solv_dict[solv_type][name] = copy.deepcopy(self.neg_ion_dict[params][name])
                            solv_ratio_type = "neg_tot_ratio"
                            molarity = solv_dict["salt_molarity"]
                    
                    if solv_dict[solv_type][name].moleculetype is not False:
                        moleculetype = solv_dict[solv_type][name].moleculetype
                    else:
                        moleculetype = False

                    ### Finds charge data from topology files
                    if len(self.ITP_INPUT_cmds) > 0 and charge == "top":
                        if solv_type == "solvent" and solv_dict["flooding"]:
                            assert_cmd_string = "solute"
                        elif solv_type == "solvent" and not solv_dict["flooding"]:
                            assert_cmd_string = "solv"
                        elif solv_type == "neg_ions":
                            assert_cmd_string = "neg"
                        elif solv_type == "pos_ions":
                            assert_cmd_string = "pos"
                        
                        assert moleculetype is not False, "\n".join([
                            "Charges are set to be obtained from topology, but this molecule does not have a moleculetype: {name}".format(name=name),
                            "If you have imported the molecule using 'molecule_import' without specifying a moleculetype, then please change solvation/flooding command to the following:",
                            "'" + assert_cmd_string + ":" + ":".join(sub_cmd) + ":charge:lib'",
                        ])
                        assert moleculetype in self.itp_moleculetypes.keys(), "\n".join([
                            "Solvent/ion/solute name '{moleculetype}' could not be found in the topology.".format(moleculetype=moleculetype),
                            "You can exempt all solvents/solutes/ions (in this solvation/flooding command) from being searched for in the topology by adding the following to the solvation/flooding command:",
                            "    "+"charge:lib",
                            "Alternatively, you can exempt a solvent/solute/ion from being searched for in the topology by adding the following to the solvent/solute/ion command:",
                            "    "+"charge:lib",
                            "Example:",
                            "    "+"neg:CL:charge:lib",
                        ])
                        bead_charges = list(map(
                            self.get_number_from_string,
                            self.itp_moleculetypes[moleculetype].topology_types_in_molecule["atoms"].get_value_type("charge")
                        ))
                        assert len(bead_charges) == len(solv_dict[solv_type][name].get_beads()), "\n".join([
                            "Mismatch found between number of beads in structure and number of beads in topology for '{name}' during preprocessing of solvation/flooding command".format(name=name)
                        ])
                        solv_dict[solv_type][name].set_bead_charges(bead_charges)
                    
                    ### If count has been set, then convert to integer value and treat as absolute number of molecules
                    if solv_dict["count"]:
                        solv_dict[solv_type][name].molarity_set(int(ratio + 0.5))
                        solv_dict[solv_type][name].ratio_set(0)
                    else:
                        solv_dict[solv_type][name].molarity_set(molarity)
                        solv_dict[solv_type][name].ratio_set(ratio)
                        solv_dict[solv_ratio_type] += ratio
                    
                    ### If mapping ratio should be ignored, set to 1 for all solvents
                    if not solv_dict["mapping"]:
                        solv_dict[solv_type][name].mapping_ratio_set(1)
                
                self.SOLVATIONS[cmd_nr] = solv_dict.copy()
                
            self.print_term("Number of solvent commands preprocessed:", len(self.SOLVATIONS), spaces=1, verbose=2)
    