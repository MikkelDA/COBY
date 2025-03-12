import numpy as np
import sys
import time
import copy

class stacked_membranes_preprocessor:
    def stacked_membranes_preprocessor(self):
        '''
        Preprocesses stacked membrane arguments for later ease of use
        '''
        box_heights = []
        if len(self.STACKED_MEMBRANES_cmds) > 0:
            self.print_term("Prepreprocessing stacked membranes request:", spaces=0, verbose=2)
            for cmd_nr, sm_cmd in enumerate(self.STACKED_MEMBRANES_cmds, 1):
                self.print_term("Starting argument:", cmd_nr, spaces=1, verbose=2)

                ### Defaults
                settings_dict = {
                    "distance": [5], # Can be set for each individual space
                    "distance_type": ["surface"], # Can be set for each individual space
                    "number": 0, # sets the number of membranes/solvent spaces
                    ### Following sets how the solvent should be treated
                    ### ### split: Split solvent position 1 across the pbc (half-half)
                    ### ### bottom: All of position 1 is located at the bottom
                    ### ### top: All of position 1 is located at the top
                    "pbc": "split", # split, bottom, top
                    "membrane_arguments": {},
                    "solvation_arguments": {},
                }
                current_sub_cmd = False
                
                if type(sm_cmd) == str:
                    split_cmds = sm_cmd.split()
                elif isinstance(sm_cmd, list, tuple):
                    split_cmds = " ".join(sm_cmd).split()
                
                arg_extensions = ("subargument", "subarg", "argument", "arg", "subcommand", "subcmd", "command", "cmd")
                
                ### ### Check protein argument
                for cmd in split_cmds:
                    sub_cmd = cmd.split(":")
                    
                    if not current_sub_cmd and not (sub_cmd[0].startswith(("membrane", "solvation")) and sub_cmd[0].endswith(arg_extensions)):
                        if sub_cmd[0].lower() == "distance":
                            settings_dict["distance"] = [self.get_number_from_string(distance) for distance in sub_cmd[1:]]

                        elif sub_cmd[0].lower() == "distance_type":
                            if sub_cmd[1] in ["surface", "center"]:
                                settings_dict["distance_type"] = [d_type for d_type in sub_cmd[1:]]

                        elif sub_cmd[0].lower() == "number":
                            settings_dict["number"] = int(self.get_number_from_string(sub_cmd[1]))
                            settings_dict["membrane_arguments"]  = {i: [] for i in range(1, settings_dict["number"] + 1)}
                            settings_dict["solvation_arguments"] = {i: [] for i in range(1, settings_dict["number"] + 1)}

                        elif sub_cmd[0].lower() == "pbc":
                            assert sub_cmd[1] in ["split", "bottom", "bot", "top"], "The 'pbc' subargument must have a value of either 'split', 'bottom'/'bot' or 'top'"
                            settings_dict["pbc"] = sub_cmd[1]
                    
                    elif sub_cmd[0].lower().startswith(("membrane", "solvation")) and sub_cmd[0].lower().endswith(arg_extensions):
                        ### Appends previous arguments to their respective positions
                        if current_sub_cmd is not False:
                            if current_sub_cmd == "membrane":
                                for position in curr_positions:
                                    settings_dict["membrane_arguments"][position].append({"subarguments": copy.deepcopy(curr_subarguments)})
                            if current_sub_cmd == "solvation":
                                for position in curr_positions:
                                    settings_dict["solvation_arguments"][position].append({"subarguments": copy.deepcopy(curr_subarguments)})
                        curr_subarguments = []

                        ### Membrane specific check and error message
                        if sub_cmd[0].lower().startswith("membrane"):
                            current_sub_cmd = "membrane"
                            assert len(sub_cmd) > 1, "\n".join([
                                "You must supply a 'position(s)' subsubargument with a 'membrane' subargument like shown below",
                                "    "+"membrane_argument:positions:1:3",
                                "    "+"membrane_argument:positions:2",
                            ])
                        ### Solvation specific check and error message
                        elif sub_cmd[0].lower().startswith("solvation"):
                            current_sub_cmd = "solvation"
                            assert len(sub_cmd) > 1, "\n".join([
                                "You must supply a 'position(s)' subsubargument with a 'solvation' subargument like shown below",
                                "    "+"solvation_argument:positions:1:3",
                                "    "+"solvation_argument:positions:2",
                            ])
                        
                        ### Finds all the positions that this set of subarguments should be used in
                        curr_positions = []
                        if len(sub_cmd) > 1:
                            cur_subcmd = False
                            i = 1
                            while i < len(sub_cmd):
                                if sub_cmd[i] in ["position", "positions"]:
                                    cur_subcmd = "positions"
                                    curr_positions.append(self.get_number_from_string(sub_cmd[i+1]))
                                    i += 2
                                else: ### Assumes another value to last given subcmd
                                    if cur_subcmd == "positions":
                                        curr_positions.append(self.get_number_from_string(sub_cmd[i]))
                                        i += 1
                                    else:
                                        assert False, "Unreadable values given to '{subargtype}_argument': {sub_cmd}".format(subargtype=current_sub_cmd, sub_cmd=sub_cmd)
                        
                    else:
                        ### Adds the subargument to a list
                        curr_subarguments.append(cmd)
                
                ### Appends the last set of arguments to their respective positions
                if current_sub_cmd is not False:
                    if current_sub_cmd == "membrane":
                        for position in curr_positions:
                            settings_dict["membrane_arguments"][position].append({"subarguments": copy.deepcopy(curr_subarguments)})
                    if current_sub_cmd == "solvation":
                        for position in curr_positions:
                            settings_dict["solvation_arguments"][position].append({"subarguments": copy.deepcopy(curr_subarguments)})

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
                
                ### Checks that all positions have been used for membrane arguments and preprocesses membrane arguments
                memb_positions = []
                for memb_position_i, memb_position_ListOfDicts in settings_dict["membrane_arguments"].items():
                    if len(memb_position_ListOfDicts) > 0:
                        memb_positions.append(memb_position_i)
                    ### Checks that subarguments have been given to all 'membrane_argument' subarguments
                    for memb_position_dict in memb_position_ListOfDicts:
                        assert len(memb_position_dict["subarguments"]) > 0, "No membrane-specific subarguments given for membrane_argument in 'stacked_membranes'"
                        ### Preprocesses the membrane arguments for later use (not needed for solvations)
                        memb_position_dict["preprocessed"] = self.memb_preprocessor(return_self = "return", cmds_given = [" ".join(memb_position_dict["subarguments"])], extra_verbose=5)
                ### Checks that all positions have been used
                assert all([i in memb_positions for i in range(1, settings_dict["number"] + 1)]), "{number} membranes have been requested but not all membrane positions have been used. Used membrane positions are: {positions}".format(number=settings_dict["number"], positions=sorted(list(set(memb_positions))))
                self.print_term("memb_positions:", memb_positions, debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                
                ### Checks that all positions have been used for solvation arguments
                solv_positions = []
                for solv_position_i, solv_position_ListOfDicts in settings_dict["solvation_arguments"].items():
                    if len(solv_position_ListOfDicts) > 0:
                        solv_positions.append(solv_position_i)
                    ### Checks that subarguments have been given to all 'solvation_argument' subarguments
                    for solv_position_dict in solv_position_ListOfDicts:
                        assert len(solv_position_dict["subarguments"]) > 0, "No solvation-specific subarguments given for solvation_argument in 'stacked_membranes'"
                ### Checks that all positions have been used
                assert all([i in solv_positions for i in range(1, settings_dict["number"] + 1)]), "{number} solvations have been requested but not all solvation positions have been used. Used solvation positions are: {positions}".format(number=settings_dict["number"], positions=sorted(list(set(solv_positions))))
                self.print_term("solv_positions:", solv_positions, debug=True, debug_keys=["stacked_membranes", "sm", "SM"])

                ### Calculates intermembrane spacing
                membranes_zmax_zmin = []
                inter_leaflet_buffer = 1.5
                ### bead_vdw added to all memb_zmax and memb_zmin to prevent beads from being placed exactly on the pbc, which can lead to them being pushed outside of it due to kicks placed on all lipid beads.
                bead_vdw = 0.264/2*10 # vdw of regular beads) / 2 * (Angstrom to nanometer convertions)

                ### Finds membrane z-delimiters based
                for position_i, memb_position_ListOfDicts in settings_dict["membrane_arguments"].items():
                    membranes_zmax_zmin.append((0, 0)) # For initial comparison
                    for memb_position_dict in memb_position_ListOfDicts:
                        self.print_term('Membrane; position_i, subarguments:', position_i, memb_position_dict["subarguments"], debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                        
                        distance_type = settings_dict["distance_type"][position_i-1]

                        memb_zmax = bead_vdw
                        if "upper_leaf" in memb_position_dict["preprocessed"]["leaflets"]:
                            memb_zmax = abs(memb_position_dict["preprocessed"]["bead_maxz"] + inter_leaflet_buffer) + bead_vdw
                        
                        memb_zmin = bead_vdw
                        if "lower_leaf" in memb_position_dict["preprocessed"]["leaflets"]:
                            memb_zmin = abs(memb_position_dict["preprocessed"]["bead_minz"] - inter_leaflet_buffer) + bead_vdw
                        
                        membranes_zmax_zmin[-1] = (max(memb_zmax, membranes_zmax_zmin[-1][0]), max(memb_zmin, membranes_zmax_zmin[-1][1]))
                        
                self.print_term("membranes_zmax_zmin:", membranes_zmax_zmin, debug=True, debug_keys=["stacked_membranes", "sm", "SM"])

                ### Calculates initial z-positionings for solvations and membranes
                memb_cz = []
                solv_cz_zupperspacing_zlowerspacing = []
                curr_z = 0
                for position_i in range(settings_dict["number"]):
                    self.print_term("", debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                    self.print_term("position_i, len(memb_position_ListOfDicts), len(solv_position_ListOfDicts):", position_i, len(memb_position_ListOfDicts), len(solv_position_ListOfDicts), debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                    distance      = settings_dict["distance"][position_i]
                    distance_type = settings_dict["distance_type"][position_i]

                    above_memb_zmin = membranes_zmax_zmin[position_i][1]
                    below_memb_zmax = membranes_zmax_zmin[position_i-1][0]

                    self.print_term("distance_type:  ", distance_type,   spaces=1, debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                    self.print_term("distance:       ", distance,        spaces=1, debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                    self.print_term("above_memb_zmin:", above_memb_zmin, spaces=1, debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                    self.print_term("below_memb_zmax:", below_memb_zmax, spaces=1, debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                    self.print_term("curr_z:         ", curr_z,          spaces=1, debug=True, debug_keys=["stacked_membranes", "sm", "SM"])

                    leaflets_lengths_combined = below_memb_zmax + above_memb_zmin
                    ### Don't count 'previous lower membranes upper leaflet' the first time.
                    ### Must be done to avoid the 'previous lower membranes upper leaflet' from being counted during the first position as it is placed at the top of the system and not "below" the solvent box.
                    if position_i != 0:
                        curr_z += below_memb_zmax

                    if distance_type == "center":
                        assert distance > leaflets_lengths_combined, "'distance' must be greater than the combined length of the longest lipids in the neighboring leaflets when using 'distance_type' 'center'"
                        full_length = distance
                    
                    elif distance_type == "surface":
                        assert distance > 0, "'distance' must be greater than 0 when using 'distance_type' 'center'"
                        full_length = leaflets_lengths_combined + distance

                    solv_zspacing_no_memb_overlap = full_length - leaflets_lengths_combined
                    solv_zspacing_lowerbox        = solv_zspacing_no_memb_overlap/2 + below_memb_zmax
                    solv_zspacing_upperbox        = solv_zspacing_no_memb_overlap/2 + above_memb_zmin
                    solvbox_full_zmid             = solv_zspacing_no_memb_overlap/2 + curr_z

                    self.print_term("full_length:                  ", full_length,                   spaces=1, debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                    self.print_term("solv_zspacing_no_memb_overlap:", solv_zspacing_no_memb_overlap, spaces=1, debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                    self.print_term("solv_zspacing_lowerbox:       ", solv_zspacing_lowerbox,        spaces=1, debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                    self.print_term("solv_zspacing_upperbox:       ", solv_zspacing_upperbox,        spaces=1, debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                    self.print_term("solvbox_full_zmid:            ", solvbox_full_zmid,             spaces=1, debug=True, debug_keys=["stacked_membranes", "sm", "SM"])

                    curr_z += solv_zspacing_upperbox + solv_zspacing_no_memb_overlap/2
                    solv_cz_zupperspacing_zlowerspacing.append((solvbox_full_zmid, solv_zspacing_upperbox, solv_zspacing_lowerbox))
                    memb_cz.append(curr_z)
                    ### Count 'last lower membranes upper leaflet' from first position during the last position
                    ### Counted after memb_cz.append() to get total box height
                    if position_i == settings_dict["number"]-1:
                        curr_z += membranes_zmax_zmin[position_i][0]

                total_height = curr_z
                self.print_term("", debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                self.print_term("Total height:", total_height, debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                
                self.print_term("", debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                self.print_term("position_i, ('solv cz', 'solv upper zspacing', 'solv lower zspacing'), 'memb cz':", debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                for position_i, ((s_cz, s_u_spacing, s_l_spacing), m_cz) in enumerate(zip(solv_cz_zupperspacing_zlowerspacing, memb_cz)):
                    self.print_term(position_i, (round(s_cz, 3), round(s_u_spacing, 3), round(s_l_spacing, 3)), round(m_cz, 3), spaces=1, debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                
                ### Adjusting z-positioning of solvation boxes and membranes depending on 'pbc' setting.
                ### ### 'solvation_individual_grouped_zmax_zmin_extra' is a list of ((zmax, zmin, extra_info),  (zmax, zmin, extra_info))
                ### ### ### 'extra_info' can either be False or (zmax_extra, zmin_extra)
                ### ### memb_cz is simply overwritten instead of creating a new list
                solvation_individual_grouped_zmax_zmin_extra = []
                curr_z = 0
                self.print_term("", debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                self.print_term("pbc:", settings_dict["pbc"], debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                self.print_term("position_i, ('solv cz', 'solv upper zspacing', 'solv lower zspacing'), ('solv zmax position', 'solv zmin position'), 'memb cz':", debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                for position_i, ((solvbox_full_zmid, solv_zspacing_upperbox, solv_zspacing_lowerbox), memb_z) in enumerate(zip(solv_cz_zupperspacing_zlowerspacing, memb_cz)):
                    position_nr = position_i + 1
                    self.print_term(position_i, (round(solvbox_full_zmid, 3), round(solv_zspacing_upperbox, 3), round(solv_zspacing_lowerbox, 3)), (round(solvbox_full_zmid+solv_zspacing_upperbox, 3), round(solvbox_full_zmid-solv_zspacing_lowerbox, 3)), round(memb_z, 3), spaces=1, debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                    if settings_dict["pbc"] in ["bottom", "bot"]:
                        ### The generic calculations above are done on a "bottom" pbc system. So most relevant values are already correct.
                        ### Note: memb_z does not need to be modified for 'bottom' 'pbc' systems
                        if position_i == 0:
                            extra_z_spacing = solvbox_full_zmid
                            z_pos0_offset = solv_zspacing_lowerbox - solvbox_full_zmid
                            self.print_term("extra_z_spacing:", extra_z_spacing, spaces=2, debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                            self.print_term("z_pos0_offset:", z_pos0_offset, spaces=2, debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                            solvation_individual_grouped_zmax_zmin_extra.append(
                                (
                                    (
                                        round(extra_z_spacing, 3),
                                        0,
                                        (round(total_height, 3), round(total_height - z_pos0_offset, 3)),
                                    ),
                                    (
                                        round(extra_z_spacing + solv_zspacing_upperbox, 3),
                                        round(extra_z_spacing, 3),
                                        False,
                                    )
                                )
                            )
                        else:
                            solvation_individual_grouped_zmax_zmin_extra.append(
                                (
                                    (
                                        round(solvbox_full_zmid, 3),
                                        round(solvbox_full_zmid - solv_zspacing_lowerbox, 3),
                                        False,
                                    ),
                                    (
                                        round(solvbox_full_zmid + solv_zspacing_upperbox, 3),
                                        round(solvbox_full_zmid, 3),
                                        False,
                                    )
                                )
                            )
                        ### memb_cz is already correct so no need to change it for "bottom" systems
                    
                    elif settings_dict["pbc"] in ["top"]:
                        if position_i == 0:
                            extra_z_spacing = solvbox_full_zmid
                            z_pos0_offset = solv_zspacing_upperbox - solvbox_full_zmid
                            z_offset = solvbox_full_zmid+solv_zspacing_upperbox - z_pos0_offset
                            self.print_term("extra_z_spacing:", extra_z_spacing, spaces=2, debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                            self.print_term("z_pos0_offset:", z_pos0_offset, spaces=2, debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                            self.print_term("z_offset:", z_offset, spaces=2, debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                            solvation_individual_grouped_zmax_zmin_extra.append(
                                (
                                    (
                                        round(total_height - extra_z_spacing, 3),
                                        round(total_height - extra_z_spacing - solv_zspacing_lowerbox, 3),
                                        False,
                                    ),
                                    (
                                        round(total_height, 3),
                                        round(total_height - extra_z_spacing, 3),
                                        (round(z_pos0_offset, 3), 0),
                                    )
                                )
                            )
                        else:
                            solvation_individual_grouped_zmax_zmin_extra.append(
                                (
                                    (
                                        round(solvbox_full_zmid                          - z_offset, 3),
                                        round(solvbox_full_zmid - solv_zspacing_lowerbox - z_offset, 3),
                                        False,
                                    ),
                                    (
                                        round(solvbox_full_zmid + solv_zspacing_upperbox - z_offset, 3),
                                        round(solvbox_full_zmid                          - z_offset, 3),
                                        False,
                                    )
                                )
                            )
                        memb_cz[position_i] = round(memb_z - z_offset, 3)
                    
                    elif settings_dict["pbc"] in ["split"]:
                        if position_i == 0:
                            z_offset = solvbox_full_zmid
                            self.print_term("z_offset:", z_offset, spaces=2, debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                            solvation_individual_grouped_zmax_zmin_extra.append(
                                (
                                    (
                                        round(total_height, 3),
                                        round(total_height - solv_zspacing_lowerbox, 3),
                                        False,
                                    ),
                                    (
                                        round(solv_zspacing_upperbox, 3),
                                        round(0, 3),
                                        False,
                                    )
                                )
                            )
                        else:
                            solvation_individual_grouped_zmax_zmin_extra.append(
                                (
                                    (
                                        round(solvbox_full_zmid                          - z_offset, 3),
                                        round(solvbox_full_zmid - solv_zspacing_lowerbox - z_offset, 3),
                                        False,
                                    ),
                                    (
                                        round(solvbox_full_zmid + solv_zspacing_upperbox - z_offset, 3),
                                        round(solvbox_full_zmid                          - z_offset, 3),
                                        False,
                                    )
                                )
                            )
                        memb_cz[position_i] = round(memb_z - z_offset, 3)

                self.print_term("", debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                self.print_term("position_i, ('Lower box zmax', 'Lower box zmin', ''('LB extra zmax', 'LB extra zmin') OR False''), ('Upper box zmax', 'Upper box zmin', ''('UB extra zmax', 'UB extra zmin') OR False''), 'Memb cz'", debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                for position_i, (((l_zmax, l_zmin, l_extra), (u_zmax, u_zmin, u_extra)), mz) in enumerate(zip(solvation_individual_grouped_zmax_zmin_extra, memb_cz)):
                    if l_extra is False:
                        l_extra = "False"
                    else:
                        l_extra = (round(l_extra[0], 3), round(l_extra[1], 3))
                    if u_extra is False:
                        u_extra = "False"
                    else:
                        u_extra = (round(u_extra[0], 3), round(u_extra[1], 3))
                    self.print_term(position_i, (round(l_zmax, 3), round(l_zmin, 3), l_extra), (round(u_zmax, 3), round(u_zmin, 3), u_extra), round(mz, 3), spaces=1, debug=True, debug_keys=["stacked_membranes", "sm", "SM"])

                ### Creating solvation and membrane arguments
                solvation_arguments = []
                membrane_arguments = []

                half_height = total_height / 2

                self.print_term("", debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                self.print_term("Membrane and Solvation arguments: 'Membrane Or Solvation', subarguments", debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                for position_i in range(settings_dict["number"]):
                    position_nr = position_i + 1 # arguments are numbered 1,2,3,... not 0,1,2,...
                    self.print_term("position_i, position_nr:", position_i, position_nr, spaces=1, debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                    memb_z = memb_cz[position_i] - half_height
                    for memb_i, memb_position_dict in enumerate(settings_dict["membrane_arguments"][position_nr]):
                        self.print_term("Membrane:", memb_i, memb_position_dict["subarguments"], spaces=2, debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                        cz = round(memb_z/10, 4)
                        membrane_arguments.append(memb_position_dict["subarguments"].copy())
                        membrane_arguments[-1].append("cz:{cz}".format(cz=cz))
                        membrane_arguments[-1] = " ".join(membrane_arguments[-1])
                    
                    for solv_i, solv_position_dict in enumerate(settings_dict["solvation_arguments"][position_nr]):
                        self.print_term("Solvation:", solv_i, solv_position_dict["subarguments"], spaces=2, debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                        for solv_zmax, solv_zmin, extra in solvation_individual_grouped_zmax_zmin_extra[position_i]:
                            cz = round(((solv_zmax + solv_zmin) / 2 - half_height)/10, 4)
                            zlength = round((solv_zmax - solv_zmin)/10, 4)
                            solvation_arguments.append(solv_position_dict["subarguments"].copy())
                            solvation_arguments[-1].append("cz:{cz}".format(cz=cz))
                            solvation_arguments[-1].append("zlength:{zlength}".format(zlength=zlength))
                            if extra is not False:
                                extra_arg = "eib:zmin:{zmin}:zmax:{zmax}".format(
                                    zmin = round(extra[1] - total_height/2, 4)/10,
                                    zmax = round(extra[0] - total_height/2, 4)/10,
                                )
                                ### Following is done in case the requested solvation is only to be done in parts of the x/y "space"
                                for sub_arg in solv_position_dict["subarguments"]:
                                    if sub_arg.startswith("cx:"):
                                        extra_arg = extra_arg + ":cx:" + str(sub_arg)
                                    elif sub_arg.startswith("cy:"):
                                        extra_arg = extra_arg + ":cy:" + str(sub_arg)
                                    elif sub_arg.startswith("xlength:"):
                                        extra_arg = extra_arg + ":xlength:" + str(sub_arg)
                                    elif sub_arg.startswith("ylength:"):
                                        extra_arg = extra_arg + ":ylength:" + str(sub_arg)
                                solvation_arguments[-1].append(extra_arg)
                            solvation_arguments[-1] = " ".join(solvation_arguments[-1])

                self.print_term("", debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                self.print_term("Final Membrane arguments:", debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                for args in membrane_arguments:
                    self.print_term(args, spaces=1, debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                
                self.print_term("", debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                self.print_term("Final Solvation arguments:", debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                for args in solvation_arguments:
                    self.print_term(args, spaces=1, debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                
                self.print_term("", debug=True, debug_keys=["stacked_membranes", "sm", "SM"])
                self.print_term("Box height:", total_height/10, spaces=1, debug=True, debug_keys=["stacked_membranes", "sm", "SM"])

                self.MEMBRANES_cmds.extend(membrane_arguments)
                self.SOLVATIONS_cmds.extend(solvation_arguments)

        self.print_term("Number of membranes preprocessed:", len(self.STACKED_MEMBRANES_cmds), verbose=2)
        
        return total_height/10
