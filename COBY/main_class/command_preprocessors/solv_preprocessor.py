import ast
import copy

class solv_preprocessor:

    def solv_preprocessor(self):
        '''
        Preprocesses solvation commands for later ease of use
        '''
        if len(self.SOLVATIONS_cmds) != 0:
            self.print_term("Preprocessing solvent requests", verbose=2)
            for cmd_nr, solvate_cmd in enumerate(self.SOLVATIONS_cmds, 1):
                self.print_term("Starting Solvent command:", cmd_nr, spaces=1, verbose=2)
                ### Defaults
                solv_dict = {
                    ### ### Solvent
                    "solvent": {}, # empty
                    "solvent_preprocessing": [], # empty
                    "solv_molarity": 55.56, # [int] or [float] # mol/L
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
                    "count": False, # [bool] Uses specific molarity value as absolute number of molecules instead of molarity ratio. Will be rounded using int(val+0.5)
                    "kick": 0.264/4, # [nm] converted to [Å]
                    "mapping": True,# bool, Whether AA-to-CG mapping should be considered
                    "ratio_method": "coarse", # "coarse", "atomistic"
                    
                    "bdx": 1.0, # [multiplier]
                    "bdy": 1.0, # [multiplier]
                    "bdz": 1.0, # [multiplier]
                    
                    "params": False, # False or str
                    "bead_radius": 0.264, # [nm] converted to [Å] # Used for volume calculations
                    "gridres": 0.264, # [nm] converted to [Å] # 1.32
                    
                    "WR": 0.264, # [nm] converted to [Å]
                    "buffer": 0.2, # [nm] converted to [Å]
                    "protein_extra_buffer": 0.2, # [nm] converted to [Å]
                    "lipid_extra_buffer": 0, # [nm] converted to [Å]
                    "solute_extra_buffer": 0, # [nm] converted to [Å]
                }
                
                ### Added "default" solvation command per user request
                if solvate_cmd.startswith("default"):
                    solvate_cmd = solvate_cmd.lstrip("default")
                    solvate_cmd = " ".join(["solv:W pos:NA neg:CL", solvate_cmd])
                
                for cmd in solvate_cmd.split():
                    sub_cmd = cmd.split(":")
                    
                    ########################
                    ### SOLVENT SPECIFIC ###
                    ########################
                    ### Solvents
                    if sub_cmd[0].lower() in ["solvent", "solv", "solute"]:
                        if solv_dict["flooding"]:
                            assert sub_cmd[0] == "solute", "'solv'/'solvent' subcommand used during solvation. Please used 'solute' instead"
                        else:
                            assert sub_cmd[0] in ["solvent", "solv"], "'solute' subcommand used during flooding. Please used 'solv'/'solvent' instead"
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
                        solv_dict["xlength"] = ast.literal_eval(sub_cmd[1])
                    ### Solvation box y-length
                    elif sub_cmd[0].lower() == "ylength":
                        solv_dict["ylength"] = ast.literal_eval(sub_cmd[1])
                    ### Solvation box z-length
                    elif sub_cmd[0].lower() == "zlength":
                        solv_dict["zlength"] = ast.literal_eval(sub_cmd[1])
                        
                    ### Solvation box x center
                    elif sub_cmd[0].lower() == "cx":
                        solv_dict["center"][0] = ast.literal_eval(sub_cmd[1])
                    ### Solvation box y center
                    elif sub_cmd[0].lower() == "cy":
                        solv_dict["center"][1] = ast.literal_eval(sub_cmd[1])
                    ### Solvation box z center
                    elif sub_cmd[0].lower() == "cz":
                        solv_dict["center"][2] = ast.literal_eval(sub_cmd[1])
                        
                    ### Solvation box center
                    elif sub_cmd[0].lower() == "center":
                        solv_dict["center"] = [ast.literal_eval(val) for val in sub_cmd[1:]]
                    
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
                    
                    #########################
                    ### FLOODING SPECIFIC ###
                    #########################
                    ### Used to check if command is made as flooding or not
                    elif sub_cmd[0].lower() == "flooding":
                        solv_dict["flooding"] = ast.literal_eval(sub_cmd[1])
                    
                    ### Errors out if unknown subcommand used, and prints the subcommand to terminal
                    else:
                        assert False, "Unknown subcommand given to '-solvate'. The subcommand is: '" + str(cmd) + "'"
                
                ### fixing a couple of values such that they are in angstrom
                solv_dict.update({
                    "xlength": solv_dict["xlength"]*10,
                    "ylength": solv_dict["ylength"]*10,
                    "zlength": solv_dict["zlength"]*10,
                    "center": [val*10 for val in solv_dict["center"]],
                    
                    "kick":        solv_dict["kick"]*10,
                    "bead_radius": solv_dict["bead_radius"]*10,
                    "gridres":     solv_dict["gridres"]*10,
                    "WR":          solv_dict["WR"]*10,
                    
                    "buffer":               solv_dict["buffer"]*10,
                    "protein_extra_buffer": solv_dict["protein_extra_buffer"]*10,
                    "lipid_extra_buffer":   solv_dict["lipid_extra_buffer"]*10,
                    "solute_extra_buffer":  solv_dict["solute_extra_buffer"]*10,
                })
                
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
                        assert name in self.solvent_dict[params].keys(), "Solvent/solute name '{name}' was not found in the parameter library '{params}'".format(name=name, params=params)
                        solv_dict[solv_type][name] = copy.deepcopy(self.solvent_dict[params][name])
                        solv_ratio_type = "solv_tot_ratio"
                        molarity = solv_dict["solv_molarity"]
                    
                    elif solv_type == "pos_ions":
                        assert name in self.ion_dict[params]["positive"].keys(), "Positive ion name '{name}' was not found in the parameter library '{params}'".format(name=name, params=params)
                        solv_dict[solv_type][name] = copy.deepcopy(self.ion_dict[params]["positive"][name])
                        solv_ratio_type = "pos_tot_ratio"
                        molarity = solv_dict["salt_molarity"]
                    
                    elif solv_type == "neg_ions":
                        assert name in self.ion_dict[params]["negative"].keys(), "Negative ion name '{name}' was not found in the parameter library '{params}'".format(name=name, params=params)
                        solv_dict[solv_type][name] = copy.deepcopy(self.ion_dict[params]["negative"][name])
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
    