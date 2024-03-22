import numpy as np
import sys
import time

class stacked_membranes_preprocessor:
    def stacked_membranes_preprocessor(self):
        '''
        Preprocesses stacked membrane commands for later ease of use
        '''
        box_heights = []
        if len(self.STACKED_MEMBRANES_cmds) > 0:
            self.print_term("Pre-prepreprocessing commands:", spaces=0, verbose=2)
            for cmd_nr, sm_cmd in enumerate(self.STACKED_MEMBRANES_cmds, 1):
                self.print_term("Starting command:", cmd_nr, spaces=1, verbose=2)

                ### Defaults
                settings_dict = {
                    "distance": [5], # Can be set for each individual space
                    "distance_type": ["surface"], # Can be set for each individual space
                    "number": 0,
                    "membrane_commands": [],
                    "solvation_commands": [],
                }
                current_sub_cmd = False
                memb_i = -1
                solv_i = -1
                
                if type(sm_cmd) == str:
                    split_cmds = sm_cmd.split()
                elif isinstance(sm_cmd, list, tuple):
                    split_cmds = " ".join(sm_cmd).split()
                
                arg_extensions = ("subargumetn", "subarg", "argument", "arg", "command", "cmd")
                
                ### ### Check protein command
                for cmd in split_cmds:
                    sub_cmd = cmd.split(":")
                    
                    if not current_sub_cmd and not (sub_cmd[0].startswith(("membrane", "solvation")) and sub_cmd[0].endswith(arg_extensions)):
                        if sub_cmd[0].lower() == "distance":
                            settings_dict["distance"] = [self.get_number_from_string(distance) for distance in sub_cmd[1:]]

                        elif sub_cmd[0].lower() == "distance_type":
                            if sub_cmd[1] in ["surface", "center", "mean_surface"]:
                                settings_dict["distance_type"] = [d_type for d_type in sub_cmd[1:]]

                        elif sub_cmd[0].lower() == "number":
                            settings_dict["number"] = int(self.get_number_from_string(sub_cmd[1]))
                    
                    elif sub_cmd[0].lower().startswith("membrane") and sub_cmd[0].lower().endswith(arg_extensions):
                        settings_dict["membrane_commands"].append({"subcommands": []})
                        memb_i += 1
                        current_sub_cmd = "membrane"
                        if len(sub_cmd) > 1:
                            cur_subcmd = False
                            positions = []
                            i = 1
                            while i < len(sub_cmd):
                                if sub_cmd[i] in ["position", "positions"]:
                                    cur_subcmd = "positions"
                                    positions.append(self.get_number_from_string(sub_cmd[i+1]))
                                    i += 2
                                else: ### Assumes another value to last given subcmd
                                    if cur_subcmd == "positions":
                                        positions.append(self.get_number_from_string(sub_cmd[i]))
                                        i += 1
                                    else:
                                        assert False, "Unreadable values given to 'membrane_argument': {sub_cmd}".format(sub_cmd)
                        else:
                            positions = False
                        settings_dict["membrane_commands"][memb_i]["positions"] = positions
                    
                    elif sub_cmd[0].lower().startswith("solvation") and sub_cmd[0].lower().endswith(arg_extensions):
                        settings_dict["solvation_commands"].append({"subcommands": []})
                        solv_i += 1
                        current_sub_cmd = "solvation"
                        if len(sub_cmd) > 1:
                            cur_subcmd = False
                            positions = []
                            i = 1
                            while i < len(sub_cmd):
                                if sub_cmd[i] in ["position", "positions"]:
                                    cur_subcmd = "positions"
                                    positions.append(self.get_number_from_string(sub_cmd[i+1]))
                                    i += 2
                                else: ### Assumes another value to last given subcmd
                                    if cur_subcmd == "positions":
                                        positions.append(self.get_number_from_string(sub_cmd[i]))
                                        i += 1
                                    else:
                                        assert False, "Unreadable values given to 'solvation_argument': {sub_cmd}".format(sub_cmd)
                        else:
                            positions = False
                        settings_dict["solvation_commands"][solv_i]["positions"] = positions
                    
                    else:
                        if current_sub_cmd == "membrane":
                            settings_dict["membrane_commands"][memb_i]["subcommands"].append(cmd)
                                
                        elif current_sub_cmd == "solvation":
                            settings_dict["solvation_commands"][solv_i]["subcommands"].append(cmd)

                ### Checks if subargument "number" has been set
                assert settings_dict["number"] > 0, "You have not set the number of membranes you wish to be stacked. Please do so with the subargument 'number:[int]'"
                
                ### Checks if mulitple "distance" values have been given
                if len(settings_dict["distance"]) == 1:
                    settings_dict["distance"] = [settings_dict["distance"][0] for _ in range(settings_dict["number"])]
                else:
                    ### Throws error if multiple distances but number of distances does not match the number of designated membranes/solvations
                    assert len(settings_dict["distance"]) > 1 and len(settings_dict["distance"]) == settings_dict["number"], "\n".join([
                        "Multiple 'distance' values ({distance}) have been given but the number of values does not equal 'number' ({number})".format(distance=str(settings_dict["distance"]), number=str(settings_dict["number"])),
                        "I you want to provide multiple 'distance' values, please ensure that you provide a number of values equal to 'number'",
                    ])
                
                ### Converts all "distance" values from [nm] to [Å]
                settings_dict["distance"] = [distance*10 for distance in settings_dict["distance"]]
                
                ### Checks if mulitple "distance_type" values have been given
                if len(settings_dict["distance_type"]) == 1:
                    settings_dict["distance_type"] = [settings_dict["distance_type"][0] for _ in range(settings_dict["number"])]
                else:
                    ### Throws error if multiple distance_types but number of distances does not match the number of designated membranes/solvations
                    assert len(settings_dict["distance_type"]) > 1 and len(settings_dict["distance_type"]) == settings_dict["number"], "\n".join([
                        "Multiple 'distance_type' values ({distance_types}) have been given but the number of values does not equal 'number' ({number})".format(distance_types=str(settings_dict["distance_type"]), number=str(settings_dict["number"])),
                        "I you want to provide multiple 'distance_type' values, please ensure that you provide a number of values equal to 'number'",
                    ])
                
                ### If only one membrane commands has been given without specifying position, then add all positions (assuming identical membranes)
                if len(settings_dict["membrane_commands"]) == 1:
                    if settings_dict["membrane_commands"][0]["positions"] is False:
                        settings_dict["membrane_commands"][0]["positions"] = [i for i in range(1, settings_dict["number"]+1)]
                
                ### If only one solvation commands has been given without specifying position, then add all positions (assuming identical solvations)
                if len(settings_dict["solvation_commands"]) == 1:
                    if settings_dict["solvation_commands"][0]["positions"] is False:
                        settings_dict["solvation_commands"][0]["positions"] = [i for i in range(1, settings_dict["number"]+1)]
                
                ### Checks that all positions have been used for membrane commands, that there are no duplicates and that each membrane command has had positions designated
                memb_positions = []
                for memb_i, memb_dict in enumerate(settings_dict["membrane_commands"]):
                    ### Checking if positions are designated
                    assert memb_dict["positions"] is not False, "\n".join([
                        "Multiple membranes have been given but their positions have not been specified by adding 'positions:[int]' to the 'membrane_argument' subargument.",
                        "Any number of values can be added to 'positions' as an shown below:",
                        "    "+"membrane_argument:positions:1:3",
                        "    "+"membrane_argument:positions:2",
                    ])
                    for position in memb_dict["positions"]:
                        ### Checks if position has not been used before
                        assert position not in memb_positions, "One membrane position has been given multiple times: {position}".format(position=position)
                        memb_positions.append(position)
                    ### Checks that subcommands have been given to a 'membrane_argument' subcommand
                    assert len(memb_dict["subcommands"]) > 0, "No membrane-specific subarguments given for membrane_argument in 'stacked_membranes'"
                    ### Preprocesses the membrane commands for later use (not needed for solvations)
                    memb_dict["preprocessed"] = self.memb_preprocessor(return_self = "return", cmds_given = [" ".join(memb_dict["subcommands"])], extra_verbose=5)
                ### Checks that all positions have been used
                assert len(memb_positions) == settings_dict["number"], "{number} membranes have been requested but not all membrane positions have been used. Used membrane positions are: {positions}".format(number=settings_dict["number"], positions=memb_positions)
                
                ### Checks that all positions have been used for solvation commands, that there are no duplicates and that each solvation command has had positions designated
                solv_positions = []
                for solv_i, solv_dict in enumerate(settings_dict["solvation_commands"]):
                    ### Checking if positions are designated
                    assert solv_dict["positions"] is not False, "\n".join([
                        "Multiple solvations have been given but their positions have not been specified by adding 'positions:[int]' to the 'solvation_argument' subargument.",
                        "Any number of values can be added to 'positions' as an shown below:",
                        "    "+"solvation_argument:positions:1:3",
                        "    "+"solvation_argument:positions:2",
                    ])
                    for position in solv_dict["positions"]:
                        ### Checks if position has not been used before
                        assert position not in solv_positions, "One solvation position has been given multiple times: {position}".format(position=position)
                        solv_positions.append(position)
                    ### Checks that subcommands have been given to a 'solvation_argument' subcommand
                    assert len(solv_dict["subcommands"]) > 0, "No solvation-specific subarguments given for solvation_argument in 'stacked_membranes'"
                ### Checks that all positions have been used
                assert len(solv_positions) == settings_dict["number"], "{number} solvations have been requested but not all membrane positions have been used. Used solvation positions are: {positions}".format(number=settings_dict["number"], positions=solv_positions)

                ### Creates position pointers
                position_memb_pointers = {}
                position_solv_pointers = {}
                for position in range(1, settings_dict["number"] + 1):
                    for memb_i, memb_dict in enumerate(settings_dict["membrane_commands"]):
                        if position in memb_dict["positions"]:
                            position_memb_pointers[position] = memb_i
                    for solv_i, solv_dict in enumerate(settings_dict["solvation_commands"]):
                        if position in solv_dict["positions"]:
                            position_solv_pointers[position] = solv_i
                
                ### Sorting according to number value just to be sure. Probably unnecessary
                position_memb_pointers = dict(sorted(position_memb_pointers.items(), key=lambda x: x[0]))
                position_solv_pointers = dict(sorted(position_solv_pointers.items(), key=lambda x: x[0]))
                
                ### Calculates intermembrane spacing
                membrane_centers = []
                solvation_delimiters = []
                cur_center = 0
                box_height = 0
                
                ### Creates membrane commands and finds solvent command z-delimiters
                for position_i, position in enumerate(position_memb_pointers.keys()):
                    
                    distance      = settings_dict["distance"][position_i]
                    distance_type = settings_dict["distance_type"][position_i]
                    
                    if position-1 in position_memb_pointers:
                        upper_position = position-1
                    else:
                        upper_position = max(position_memb_pointers.keys())
                    lower_position = position
                    
                    upper_memb = settings_dict["membrane_commands"][position_memb_pointers[upper_position]]["preprocessed"]
                    lower_memb = settings_dict["membrane_commands"][position_memb_pointers[lower_position]]["preprocessed"]
                    
                    if distance_type == "center":
                        upper_space = 0
                        lower_space = 0
                    elif distance_type == "surface":
                        upper_space = abs(upper_memb["bead_maxz"])
                        lower_space = abs(lower_memb["bead_minz"])
                    elif distance_type == "mean_surface":
                        upper_space = abs(upper_memb["bead_mean_maxz"])
                        lower_space = abs(lower_memb["bead_mean_minz"])
                    
                    ### (total distance, upper solvent space, lower solvent space)
                    upper_space = distance/2 + upper_space
                    lower_space = distance/2 + lower_space
                    total_space = upper_space + lower_space
                    
                    if position_i == 0:
                        first_position_lower_delimiters = (cur_center, cur_center + lower_space)
                        last_space = round(upper_space, 3)
                        cur_center += round(lower_space, 3)
                    else:
                        lower_delimiters = (cur_center, cur_center + lower_space)
                        cur_center += round(lower_space, 3)
                        upper_delimiters = (cur_center, cur_center + upper_space)
                        cur_center += round(upper_space, 3)
                        solvation_delimiters.append((lower_delimiters, upper_delimiters))
                    box_height += total_space
                    membrane_centers.append(cur_center)
                
                first_position_upper_delimiters = (cur_center, cur_center + last_space)
                solvation_delimiters = [(first_position_lower_delimiters, first_position_upper_delimiters)] + solvation_delimiters
                
                membrane_commands = []
                for position_i, (position, memb_pointer) in enumerate(position_memb_pointers.items()):
                    membrane_command = settings_dict["membrane_commands"][memb_pointer]["subcommands"].copy()
                    membrane_center = round((membrane_centers[position_i] - box_height/2)/10, 3)
                    membrane_command.append("cz:" + str(membrane_center))
                    membrane_commands.append(" ".join(membrane_command))
                
                ### Creates solvation commands
                solvation_commands = []
                for position_i, (position, solv_pointer) in enumerate(position_solv_pointers.items()):
                    main_solvation_command = settings_dict["solvation_commands"][solv_pointer]["subcommands"].copy()
                    
                    for bot_delimiter, top_delimiter in solvation_delimiters[position_i]:
                        cz = round((np.mean([top_delimiter, bot_delimiter]) - box_height/2)/10, 3)
                        z = round((top_delimiter - bot_delimiter)/10, 3)
                        solvation_command = main_solvation_command.copy()
                        solvation_command.append("cz:" + str(cz))
                        solvation_command.append("zlength:" + str(z))
                        solvation_commands.append(" ".join(solvation_command))
                
                self.MEMBRANES_cmds.extend(membrane_commands)
                self.SOLVATIONS_cmds.extend(solvation_commands)
                box_heights.append(round(box_height/10, 3))

        self.print_term("Number of membranes preprocessed:", len(self.STACKED_MEMBRANES_cmds), spaces=1, verbose=2)
        
        return max(box_heights)
