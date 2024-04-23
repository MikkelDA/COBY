import time
import math
import numpy as np
import random
import shapely
from operator import itemgetter
import copy

class solvater:
    def solvater(self):
        if len(self.SOLVATIONS_cmds) != 0:
            solvation_tic = time.time()
            string = " ".join(["", "SOLVATING SYSTEM", ""])
            self.print_term("{string:-^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), spaces=0, verbose=1)
            
            for solvation_i, (solvation_nr, solvation) in enumerate(self.SOLVATIONS.items()):
                if solvation_i != 0:
                    self.print_term("", verbose=2)
                self.print_term("Starting solvation nr", solvation_nr, spaces=0, verbose=2)

                #########################################
                ### CALCULATION FREE VOLUME OF SYSTEM ###
                #########################################
                self.print_term("Calculating box volume: (all values in [nm^3])", spaces=1, verbose=2)
                self.print_term("Bead radius used for volume calculations 'bead_radius':", solvation["bead_radius"]/10, "[nm]", spaces=2, verbose=2)

                bead_radius = solvation["bead_radius"]
                bead_volume = (4/3 * math.pi * (bead_radius ** 3)) * 10**-27
                
                cx, cy, cz = solvation["center"]
                xlen, ylen, zlen = itemgetter("xlength", "ylength", "zlength")(solvation)
                cxmin, cxmax = cx-xlen/2, cx+xlen/2
                cymin, cymax = cy-ylen/2, cy+ylen/2
                czmin, czmax = cz-zlen/2, cz+zlen/2
                
                self.print_term("cx, cy, cz:", cx, cy, cz, debug=True)
                self.print_term("xlen, ylen, zlen:", xlen, ylen, zlen, debug=True)
                self.print_term("cxmin, cxmax:", cxmin, cxmax, debug=True)
                self.print_term("cymin, cymax:", cymin, cymax, debug=True)
                self.print_term("czmin, czmax:", czmin, czmax, debug=True)
                
                ### Solvent/solute volume
                solvent_beads_for_cell_checker, solvent_box_solute_charge = self.get_solute_volume(
                    cxmin, cxmax, cymin, cymax, czmin, czmax
                )
                solvs_volume = bead_volume * len(solvent_beads_for_cell_checker)
                
                ### Lipid volume
                lipid_beads_for_cell_checker, solvent_box_lipids_charge = self.get_lipid_volume(
                    cxmin, cxmax, cymin, cymax, czmin, czmax
                )
                leafs_volume = bead_volume * len(lipid_beads_for_cell_checker)

                ### Protein volume
                prot_beads_for_cell_checker, solvent_box_proteins_charge = self.get_protein_volume(
                    cxmin, cxmax, cymin, cymax, czmin, czmax
                )
                prots_volume = bead_volume * len(prot_beads_for_cell_checker)
                
                n_lipids = 0
                if solvation["solv_per_lipid"]:
                    if len(self.MEMBRANES) > 0:
                        ### Counts the number of lipids inside the solvent box
                        for memb_key, memb_dict in self.MEMBRANES.items():
                            for leaflet_key, leaflet in memb_dict["leaflets"].items():
                                for grid_point in leaflet["grid_lipids"]:
                                    xs, ys, zs = itemgetter("x", "y", "z")(grid_point["lipid"])
                                    beads_total = len(xs)
                                    beads_inside = 0
                                    for x, y, z in zip(xs, ys, zs):
                                        if cxmin < x <= cxmax and cymin < y <= cymax and czmin < z <= czmax:
                                            beads_inside += 1
                                    if beads_inside > beads_total/2:
                                        n_lipids += 1
                
                solvent_box_charge = solvent_box_solute_charge + solvent_box_lipids_charge + solvent_box_proteins_charge
                
                N_A = 6.02214076 * 10**23
                box_volume = (xlen * ylen * zlen) * 10**-27
                non_free_volume = leafs_volume + prots_volume + prots_volume
                box_free_volume = box_volume - non_free_volume
                self.print_term("Solvent box volume:", round(box_volume      * 10**24, 3), spaces=2, verbose=2)
                self.print_term("Excluded volume:   ", round(non_free_volume * 10**24, 3), spaces=2, verbose=2)
                
                self.print_term("Solute volume: ", round(solvs_volume    * 10**24, 3), spaces=3, verbose=2)
                self.print_term("Lipid volume:  ", round(leafs_volume    * 10**24, 3), spaces=3, verbose=2)
                self.print_term("Protein volume:", round(prots_volume    * 10**24, 3), spaces=3, verbose=2)
                
                self.print_term("Free volume:       ", round(box_free_volume * 10**24, 3), spaces=2, verbose=2)
                
                ###########################
                ### CALCULATING SOLVENT ###
                ###########################
                solvent_box_solvent_charge = 0

                ### List used to find the maximum solvent/ion size
                solvent_radii = []
                
                solvent_radii = [
                    vals.get_radius()
                    for dict_pointer in ["solvent", "pos_ions", "neg_ions"]
                    for key, vals in solvation[dict_pointer].items()
                ]
                
                def get_solvent_ratios(solvent_type_dict, tot_ratio):
                    ratios = [
                        round(1 / tot_ratio * vals.ratio, 4)
                        for key, vals in solvent_type_dict.items()
                    ]
                    return ratios
                
                if not solvation["count"]:
                    sol_ratios = get_solvent_ratios(solvation["solvent"], solvation["solv_tot_ratio"])
                    pos_ratios = get_solvent_ratios(solvation["pos_ions"], solvation["pos_tot_ratio"])
                    neg_ratios = get_solvent_ratios(solvation["neg_ions"], solvation["neg_tot_ratio"])
                else:
                    sol_ratios = [1 for _ in solvation["solvent"].keys()]
                    pos_ratios = [1 for _ in solvation["pos_ions"].keys()]
                    neg_ratios = [1 for _ in solvation["neg_ions"].keys()]

                def solvent_count_calculator(solvent_type, ratios, volume):

                    if solvent_type == "solvent":
                        molarity = solvation["solv_molarity"]
                    else:
                        molarity = solvation["salt_molarity"]
                    
                    N_A = 6.02214076 * 10**23
                    number_of_particles = N_A * volume * molarity

                    if solvation["count"]:
                        counts = [vals.molarity for vals in solvation[solvent_type].values()]
                    elif solvation["solv_per_lipid"] and solvent_type == "solvent":
                        ### Int rounding
                        counts = [int((solvation["solv_per_lipid"] * n_lipids * ratio)+0.5) for ratio in ratios]

                    elif solvation["ratio_method"] == "coarse":
                        ### Treats molarity as molarity
                        ### Ratios are treated as ratios between number of inserted CG molecules
                        mappings_weighted_mean = np.sum([vals.mapping_ratio * ratio for vals, ratio in zip(solvation[solvent_type].values(), ratios)])
                        ### Int rounding
                        counts = [int((number_of_particles * ratio / mappings_weighted_mean)+0.5) for ratio in ratios]

                    elif solvation["ratio_method"] == "atomistic":
                        ### Treats molarity as molarity
                        ### Ratios are treated as ratios between number of atoms (taking mapping ratio of solvent into account)
                        ### Int rounding
                        counts = [int((number_of_particles * ratio / vals.mapping_ratio)+0.5) for vals, ratio in zip(solvation[solvent_type].values(), ratios)]

                    return counts
                
                ### Using the free volume for concentration calculations (default)
                if solvation["solvfreevol"] == True:
                    volume_for_solv = box_free_volume
                ### Using box volume for concentration calculations
                elif solvation["solvfreevol"] == False:
                    volume_for_solv = box_volume
                
                sol_counts = solvent_count_calculator(
                    "solvent",
                    sol_ratios,
                    volume_for_solv,
                )
                
                solv_volume = 0
                sol_charges = 0

                for (key, vals), count in zip(solvation["solvent"].items(), sol_counts):
                    vals.count_set(count)
                    sol_charges += vals.count * vals.get_mol_charge()
                    ### Only throws warnings if you are actually trying to insert ions
                    if solvation["pos_ions"] and solvation["neg_ions"]:
                        if not vals.molar_mass and solvation["ionsvol"] == "solv" and not solvation["count"]:
                            self.print_term("WARNING: Chosen solvent [" + key + "] is missing a 'molar_mass' entry in the solvent defines. Please add it if you wish to use solvent volume for ion concentration.", warn = True)
                        if not vals.density and solvation["ionsvol"] == "solv" and not solvation["count"]:
                            self.print_term("WARNING: Chosen solvent [" + key + "] is missing a 'density' entry in the solvent defines. Please add it if you wish to use solvent volume for ion concentration.", warn = True)
                        if vals.molar_mass and vals.density:
                            solv_volume += (count * vals.mapping_ratio * vals.molar_mass) / (N_A * vals.density) * 10**-3
                
                self.print_term("Solvent volume:    ", round(solv_volume * 10**24, 3), spaces=2, verbose=2)

                solvent_box_solvent_charge += sol_charges

                def ions_neutraliser(pos_ions, neg_ions, target_charge, current_charge, direction):
                    charge_difference = target_charge - current_charge
                    if charge_difference > 0:
                        ions = pos_ions
                    elif charge_difference < 0:
                        ions = neg_ions
                    else:
                        return pos_ions, neg_ions
                    
                    if direction == "add":
                        sign = +1
                    if direction == "remove":
                        sign = -1
                    
                    keys = list(ions.keys())
                    keys_sorted_after_ratio = list(sorted(keys, key=lambda key: ions[key]["ratio"], reverse = True))
                    tot_ions = sum([vals["count"] for vals in ions.values()])
                    counter = 0
                    while round(current_charge) != round(target_charge):
                        counter += 1
                        ### Find which ion is furthest from ideal ratio
                        if tot_ions == 0:
                            ### No ions yet, just sort them from highest to lowest ratio
                            furthest_dist_keys = keys_sorted_after_ratio[::sign]
                        else:
                            ### Calculate how far each ion is from its ideal ratio
                            dists_from_ratio = []
                            for key, vals in list(ions.items()):
                                dists_from_ratio.append((key, vals["count"]/tot_ions - vals["ratio"]))
                            ### Sort them from most underrepresented to most overrepresented
                            ### Secondarily sort them according to their ideal ratios in case the are equally underrepresented
                            dists_from_ratio_sorted = sorted(
                                dists_from_ratio,
                                key=lambda x: (x[1], keys_sorted_after_ratio.index(x[0])),
                                reverse = False
                            )
                            ### Get just the keys
                            furthest_dist_keys = [key for key, dist_from_ratio in dists_from_ratio_sorted][::sign]
                        
                        ### Checks from most underrepresented to most overrepresented if adding another ion
                        ### would cause the total added charge to become greater than the difference that is being corrected
                        for key in furthest_dist_keys:
                            if abs(target_charge - current_charge) >= abs(ions[key]["charge"]):
                                furthest_dist_key = key
                                break
                        else:
                            ### If an ion should be added but no ion can be added without adding
                            ### more charge than is needed for neutralization
                            ### then just return as is
                            return pos_ions, neg_ions
                        
                        tot_ions += sign
                        current_charge += ions[furthest_dist_key]["charge"]
                        ions[furthest_dist_key]["count"] += sign
                        
                    return pos_ions, neg_ions
                
                ### Using the solvent volume for concentration calculations (default)
                if solvation["ionsvol"] == "solv" and solv_volume > 0 and not solvation["count"]:
                    volume_for_ions = solv_volume
                ### Using the free volume for concentration calculations
                elif solvation["ionsvol"] == "free" or (solvation["ionsvol"] == "solv" and solvation["count"]):
                    volume_for_ions = box_free_volume
                ### Using the box volume for concentration calculations
                elif solvation["ionsvol"] == "box":
                    volume_for_ions = box_volume

                ### Only runs ion-specific calculations if both positive and negative ions are designated.
                if pos_ratios and neg_ratios:
                    ### First ensures that the ion concentration is equal to salt_molarity
                    pos_counts = solvent_count_calculator(
                        "pos_ions",
                        pos_ratios,
                        volume_for_ions
                    )
                
                    neg_counts = solvent_count_calculator(
                        "neg_ions",
                        neg_ratios,
                        volume_for_ions
                    )
                    
                    pos_charges = 0
                    neg_charges = 0
                    
                    ions_dict = {}
                    zipped1 = zip(["pos_ions", "neg_ions"], [pos_counts, neg_counts], [pos_ratios, neg_ratios])
                    for ion_type, ion_counts, ion_ratios in zipped1:
                        ions_dict[ion_type] = {}
                        zipped2 = zip(solvation[ion_type].items(), ion_counts, ion_ratios)
                        for (key, vals), count, ratio in zipped2:
                            vals.count_set(count)
                            if ion_type == "pos_ions":
                                pos_charges += vals.count * vals.get_mol_charge()
                            elif ion_type == "neg_ions":
                                neg_charges += vals.count * vals.get_mol_charge()
                            ions_dict[ion_type][key] = {
                                "charge": vals.get_mol_charge(),
                                "ratio": ratio,
                                "count": vals.count,
                            }
                    
                    current_charge = solvent_box_charge + sol_charges + pos_charges + neg_charges
                    
                    pos_charges = 0
                    neg_charges = 0
                    
                    if solvation["salt_method"] == "add":
                        '''
                        Adds extra ions to neutralize the solvent box
                        '''
                        pos_ions, neg_ions = ions_neutraliser(
                            ions_dict["pos_ions"],
                            ions_dict["neg_ions"],
                            solvation["target_charge"],
                            current_charge,
                            direction = "add",
                        )
                        zipped = zip(["pos_ions", "neg_ions"], [pos_ions, neg_ions])
                        for ion_type, ions_dicts in zipped:
                            for (key, vals), (ions_dict_key, ions_dict_vals) in zip(solvation[ion_type].items(), ions_dicts.items()):
                                vals.count_set(ions_dict_vals["count"])

                    elif solvation["salt_method"] == "remove":
                        '''
                        Removes excess ions to neutralize the solvent box
                        '''
                        pos_ions, neg_ions = ions_neutraliser(
                            ions_dict["pos_ions"],
                            ions_dict["neg_ions"],
                            solvation["target_charge"],
                            -current_charge,
                            direction = "remove",
                        )
                        zipped = zip(["pos_ions", "neg_ions"], [pos_ions, neg_ions])
                        for ion_type, ions_dicts in zipped:
                            for (key, vals), (ions_dict_key, ions_dict_vals) in zip(solvation[ion_type].items(), ions_dicts.items()):
                                vals.count_set(ions_dict_vals["count"])
                        
                    elif solvation["salt_method"] == "mean":
                        '''
                        Adds ions with either a positive or negative charge and removes the other type
                        Ends up with ion concentrations in between "add" and "remove"
                        Does an extra ion addition at the end to ensure no errors have been made when
                            adding/removing ions.
                        '''
                        ### Need to deepcopy, otherwise the nested dictionaries will be modified 
                        ions_dict_add = copy.deepcopy(ions_dict)
                        pos_ions_add, neg_ions_add = ions_neutraliser(
                            ions_dict_add["pos_ions"],
                            ions_dict_add["neg_ions"],
                            solvation["target_charge"],
                            +math.ceil(current_charge/2),
                            direction = "add",
                        )
                        ions_dict_remove = copy.deepcopy(ions_dict)
                        pos_ions_remove, neg_ions_remove = ions_neutraliser(
                            ions_dict_remove["pos_ions"],
                            ions_dict_remove["neg_ions"],
                            solvation["target_charge"],
                            -math.floor(current_charge/2),
                            direction = "remove",
                        )
                        
                        zipped1 = zip(["pos_ions", "neg_ions"],
                                    [pos_counts, neg_counts],
                                    [pos_ions_add, neg_ions_add],
                                    [pos_ions_remove, neg_ions_remove])
                        for ion_type, ion_counts, ions_dicts_add, ions_dicts_remove in zipped1:
                            zipped2 = zip(solvation[ion_type].items(),
                                        ion_counts,
                                        ions_dicts_add.values(),
                                        ions_dicts_remove.values())
                            for (key, vals), count, ions_dict_add_vals, ions_dict_remove_vals in zipped2:
                                add_diff = abs(count - ions_dict_add_vals["count"])
                                rem_diff = abs(count - ions_dict_remove_vals["count"])
                                vals.count_set(count + add_diff - rem_diff)

                    ### Extra ion addtion to prevent errors (e.g. non-neutralized systems)
                    ions_dict = {}
                    for ion_type, ion_ratios in zip(["pos_ions", "neg_ions"], [pos_ratios, neg_ratios]):
                        ions_dict[ion_type] = {}
                        for (key, vals), ratio in zip(solvation[ion_type].items(), ion_ratios):
                            
                            ### Checks if negative number of the ion should be inserted and set to zero if it is the case.
                            ### Incase "remove" or "mean" methods have removed more ions than should even be present.
                            if vals.count < 0:
                                vals.count_set(0)

                            ions_dict[ion_type][key] = {
                                "charge": vals.get_mol_charge(),
                                "ratio": ratio,
                                "count": vals.count,
                            }

                            if ion_type == "pos_ions":
                                pos_charges += vals.count * vals.get_mol_charge()
                            elif ion_type == "neg_ions":
                                neg_charges += vals.count * vals.get_mol_charge()
                    
                    current_charge = solvent_box_charge + sol_charges + pos_charges + neg_charges

                    pos_ions, neg_ions = ions_neutraliser(
                        ions_dict["pos_ions"],
                        ions_dict["neg_ions"],
                        solvation["target_charge"],
                        current_charge,
                        direction = "add",
                    )

                    pos_charges = 0
                    neg_charges = 0

                    for ion_type, ions_dicts in zip(["pos_ions", "neg_ions"], [pos_ions, neg_ions]):
                        for (key, vals), (ions_dict_key, ions_dict_vals) in zip(solvation[ion_type].items(), ions_dicts.items()):
                            vals.count_set(ions_dict_vals["count"])
                            if ion_type == "pos_ions":
                                pos_charges += vals.count * vals.get_mol_charge()
                            elif ion_type == "neg_ions":
                                neg_charges += vals.count * vals.get_mol_charge()
                    
                    solvent_box_solvent_charge += pos_charges + neg_charges

                ### Adding charges from ions and solvent to system charge
                self.system_charge += solvent_box_solvent_charge
                
                self.print_term("", verbose=2)
                self.print_term("Solvent box charge information:",                                            spaces=1, verbose=2)
                self.print_term("Solvent box charge before solvent insertion:", round(solvent_box_charge, 1), spaces=2, verbose=2)
                self.print_term("Prior solvent beads: ", round(solvent_box_solute_charge, 1),                 spaces=3, verbose=2)
                self.print_term("Lipid beads:         ", round(solvent_box_lipids_charge, 1),                 spaces=3, verbose=2)
                self.print_term("Protein beads:       ", round(solvent_box_proteins_charge, 1),               spaces=3, verbose=2)
                solvent_box_charge += solvent_box_solvent_charge
                self.print_term("Solvent box charge after solvent insertion: ", round(solvent_box_charge, 1), spaces=2, verbose=2)
                self.print_term("New solvent beads:", round(solvent_box_solvent_charge, 1),                   spaces=3, verbose=2)
                self.print_term("Solvent       (solv):", round(sol_charges, 1),                               spaces=4, verbose=2)
                if solvation["pos_ions"]:
                    self.print_term("Positive ions (pos): ", round(pos_charges, 1),                               spaces=4, verbose=2)
                if solvation["neg_ions"]:
                    self.print_term("Negative ions (neg): ", round(neg_charges, 1),                               spaces=4, verbose=2)

                ### Randomly shuffle solvent particles for random distribution in box
                ### Finding number of particles here to prevent extra calls to "get_res_beads_info()" later
                collected_solvent = []
                for stype in ["solvent", "pos_ions", "neg_ions"]:
                    for key, vals in solvation[stype].items():
                        sname = key
                        snbeads = len(vals.get_res_beads_info())
                        scount = vals.count
                        collected_solvent.extend([[sname, stype, snbeads]] * scount)
                
                ##############################################
                ### RUNNING SOLVENT OPTIMIZATION ALGORITHM ###
                ##############################################
                ### Finds the maximum size of molecules used as solvent/ions
                ### Also includes buffer/kick size to prevent edge overlap cases
                ### *2 to get actual size (diameter) instead of radius
                max_mol_size = max(solvent_radii)*2
                
                self.print_term("max_mol_size:    ", max_mol_size,          debug=True)
                self.print_term("max_mol_size*1.2:", max_mol_size*1.2,      debug=True)
                self.print_term("kick:            ", solvation["kick"],     debug=True)
                self.print_term("kick*2.2:        ", solvation["kick"]*2.2, debug=True)
                
                gridres = solvation["gridres"]
                
                self.print_term("gridres first:",    solvation["gridres"],  debug=True)
                
                ### if the largest molecule + 2*kick is bigger than the designated grid resolution then change the gridres
                if (max_mol_size+solvation["kick"]*2) >= gridres:
                    gridres = (max_mol_size + solvation["kick"]*2)*1.2 # 20% extra
                    self.print_term("Requested solvent is too large for the grid resolution. Adjusting grid resolution to prevent molecule overlaps.", warn = True)
                    self.print_term("Original grid resolution was:", round(solvation["gridres"]/10, 4), "[nm]", warn = True)
                    self.print_term("New grid resolution is:      ", round(gridres/10, 4), "[nm]", warn = True)
                    self.print_term("Calculation: (largest molecule size + 2 * particle kick value) * 1.2", "\n",  warn = True)
                
                self.print_term("gridres last:", gridres, debug=True)
                
                solvent_buffer = solvation["buffer"]

                #################################
                ### CHOOSES ALGORITHM VERSION ###
                #################################
                self.print_term("", verbose=2)
                self.print_term("Calculating the number of available grid points", spaces=1, verbose=2)

                def coord_to_indices(pos, dim, center, real_gridres, max_int, min_buffer, max_buffer):
                    ### Buffer limit
                    bead_min = pos-min_buffer
                    bead_max = pos+max_buffer

                    ### Convert solvent buffer limit to decimal index
                    beadi_min_dec = (bead_min+dim/2-center)/real_gridres
                    beadi_max_dec = (bead_max+dim/2-center)/real_gridres

                    ### Round buffer decimal index to integer index
                    beadi_min_int = int(beadi_min_dec)
                    beadi_max_int = int(beadi_max_dec)+1 # +1 due to last index not being counted
                    
                    if beadi_min_int < 0:
                        beadi_min_int = 0
                    if beadi_max_int < 0:
                        beadi_max_int = 0
                        
                    return beadi_min_int, beadi_max_int
                
                enough_grid_points = False
                while not enough_grid_points:
                    
                    xpoints, ypoints, zpoints = int(xlen/gridres), int(ylen/gridres), int(zlen/gridres)
                    ### Calculates actual coordinate ranges for each axis and calculates the "real" grid resolution
                    xcoords, xreal_gridres = np.linspace(cxmin-gridres/2, cxmax+gridres/2, xpoints+2, retstep=True)
                    ycoords, yreal_gridres = np.linspace(cymin-gridres/2, cymax+gridres/2, ypoints+2, retstep=True)
                    zcoords, zreal_gridres = np.linspace(czmin-gridres/2, czmax+gridres/2, zpoints+2, retstep=True)
                    ### Removes the first and last points as they are the actual edges of the box
                    xcoords = xcoords[1:-1]
                    ycoords = ycoords[1:-1]
                    zcoords = zcoords[1:-1]
                    
                    solute_buffer  = solvent_buffer+solvation["solute_extra_buffer"]
                    lipid_buffer   = solvent_buffer+solvation["lipid_extra_buffer"]
                    protein_buffer = solvent_buffer+solvation["protein_extra_buffer"]
                    
                    self.print_term("xpoints, xreal_gridres", xpoints, round(xreal_gridres, 4), debug=True)
                    self.print_term("ypoints, yreal_gridres", ypoints, round(yreal_gridres, 4), debug=True)
                    self.print_term("zpoints, zreal_gridres", zpoints, round(zreal_gridres, 4), debug=True)
                    self.print_term("points total          ", xpoints*ypoints*zpoints,          debug=True)
                    self.print_term("solvent_buffer        ", solvent_buffer,                   debug=True)
                    self.print_term('solvation["buffer"]   ', solvent_buffer,                   debug=True)
                    self.print_term("solute_buffer         ", solute_buffer,                    debug=True)
                    self.print_term("lipid_buffer          ", lipid_buffer,                     debug=True)
                    self.print_term("protein_buffer        ", protein_buffer,                   debug=True)
                    
                    ### Creates a 3D matrix indicating the points in space that are allowed to have solvent
                    ### Starts out being completely filled with ones indicating allowed space
                    grid_bool_matrix = np.ones((xpoints, ypoints, zpoints), dtype=np.int8)
                    
                    ### ### Marks coordinates as "occupied" by setting boolean to False
                    ### Marks points located in hydrophobic volume
                    if len(self.MEMBRANES) != 0:
                        for memb_key, memb_dict in self.MEMBRANES.items():
                            for leaflet_key, leaflet in memb_dict["leaflets"].items():
                                lcz = leaflet["center"][2]
                                ### Z-related stuff is always the same for all subleaflets
                                if leaflet["HG_direction"] == "up":
                                    zminbuffer = solvent_buffer
                                    zmaxbuffer = leaflet["lipid_dimensions"]["min_max_zs"] + solvent_buffer
                                if leaflet["HG_direction"] == "down":
                                    ### Also"min_max_zs" here as the lipid dimensions are done for a generalized upwards pointing lipid
                                    ### It is subtracted inside coord_to_indices()
                                    zminbuffer = leaflet["lipid_dimensions"]["min_max_zs"] + solvent_buffer
                                    zmaxbuffer = solvent_buffer
                                zmin, zmax = coord_to_indices(lcz, zlen, cz, zreal_gridres, zpoints, zminbuffer, zmaxbuffer)
                                
                                ### Advanced (and slower) solvent-in-membrane prevention method
                                ### Creates a series of points in 2D and marks matrix coordinates one point at a time
                                ### Allows holes to be solvated
                                if leaflet["solvate_hole"] and any([leaflet["remove_union_without_protein"], leaflet["require_union"], leaflet["inside_protein"]]):
                                    geoms = self.get_geoms_list(leaflet["holed_bbox_without_protein_removed"])
                                    for geom in geoms:
                                        poly_xmin, poly_ymin, poly_xmax, poly_ymax = geom.bounds
                                        lcx = np.mean([poly_xmax, poly_xmin])
                                        lcy = np.mean([poly_ymax, poly_ymin])
                                        llx = poly_xmax - poly_xmin
                                        lly = poly_ymax - poly_ymin

                                        xmin_i, xmax_i = coord_to_indices(lcx, xlen, cx, xreal_gridres, xpoints, llx/2, llx/2)
                                        ymin_i, ymax_i = coord_to_indices(lcy, ylen, cy, yreal_gridres, ypoints, lly/2, lly/2)
                                        xcoords_in_leaflet = xcoords[xmin_i:xmax_i+1]
                                        ycoords_in_leaflet = ycoords[ymin_i:ymax_i+1]

                                        points_to_check = [(round(x, 3), round(y, 3)) for x in xcoords_in_leaflet for y in ycoords_in_leaflet]
                                        Points_to_check = [shapely.Point(point) for point in points_to_check]
                                        
                                        geom_buffered = shapely.buffer(geom, solvent_buffer)
                                        points_contained_bools = geom_buffered.contains(Points_to_check)
                                        points_in_geom = [point for point, point_bool in zip(points_to_check, points_contained_bools) if point_bool]
                                        
                                        for lcx, lcy in points_in_geom:
                                            xmin, xmax = coord_to_indices(lcx, xlen, cx, xreal_gridres, xpoints, solvent_buffer, solvent_buffer)
                                            ymin, ymax = coord_to_indices(lcy, ylen, cy, yreal_gridres, ypoints, solvent_buffer, solvent_buffer)
                                            grid_bool_matrix[xmin:xmax, ymin:ymax, zmin:zmax] = 0

                                ### Simple (and fast) solvent-in-membrane prevention method
                                ### Simply marks all matrix coordinates overlapping with the membrane
                                ### Prevents holes from being solvated
                                else:
                                    lcx, lcy = leaflet["center"][:2] # Center of leaflet on given axis
                                    llx, lly = leaflet["xlength"], leaflet["ylength"] # Length of leaflet in given axis
                                    xmin, xmax = coord_to_indices(lcx, xlen, cx, xreal_gridres, xpoints, llx/2, llx/2)
                                    ymin, ymax = coord_to_indices(lcy, ylen, cy, yreal_gridres, ypoints, lly/2, lly/2)
                                    grid_bool_matrix[xmin:xmax, ymin:ymax, zmin:zmax] = 0

                    ### Prior solvent beads
                    for xpos, ypos, zpos in solvent_beads_for_cell_checker:
                        ### Checks if any non-solvent bead is within the solvent cell
                        xmin, xmax = coord_to_indices(xpos, xlen, cx, xreal_gridres, xpoints, solute_buffer, solute_buffer)
                        ymin, ymax = coord_to_indices(ypos, ylen, cy, yreal_gridres, ypoints, solute_buffer, solute_buffer)
                        zmin, zmax = coord_to_indices(zpos, zlen, cz, zreal_gridres, zpoints, solute_buffer, solute_buffer)
                        grid_bool_matrix[xmin:xmax, ymin:ymax, zmin:zmax] = 0

                    ### Lipid beads
                    for xpos, ypos, zpos in lipid_beads_for_cell_checker:
                        ### Checks if any lipid bead is within the solvent cell
                        xmin, xmax = coord_to_indices(xpos, xlen, cx, xreal_gridres, xpoints, lipid_buffer, lipid_buffer)
                        ymin, ymax = coord_to_indices(ypos, ylen, cy, yreal_gridres, ypoints, lipid_buffer, lipid_buffer)
                        zmin, zmax = coord_to_indices(zpos, zlen, cz, zreal_gridres, zpoints, lipid_buffer, lipid_buffer)
                        grid_bool_matrix[xmin:xmax, ymin:ymax, zmin:zmax] = 0

                    ### Protein beads
                    for xpos, ypos, zpos in prot_beads_for_cell_checker:
                        ### Checks if any protein bead is within the solvent cell
                        xmin, xmax = coord_to_indices(xpos, xlen, cx, xreal_gridres, xpoints, protein_buffer, protein_buffer)
                        ymin, ymax = coord_to_indices(ypos, ylen, cy, yreal_gridres, ypoints, protein_buffer, protein_buffer)
                        zmin, zmax = coord_to_indices(zpos, zlen, cz, zreal_gridres, zpoints, protein_buffer, protein_buffer)
                        
                        grid_bool_matrix[xmin:xmax, ymin:ymax, zmin:zmax] = 0
                    
                    ### Gets all the indices that are not occupied
                    free_indices = grid_bool_matrix.nonzero()
                    free_xs, free_ys, free_zs = free_indices
                    free_xs_coords = xcoords[free_xs]
                    free_ys_coords = ycoords[free_ys]
                    free_zs_coords = zcoords[free_zs]
                    
                    solv_grid_3D = list(zip(free_xs_coords, free_ys_coords, free_zs_coords))
                    
                    if len(solv_grid_3D) >= len(collected_solvent):
                        enough_grid_points = True
                    else:
                        self.print_term("Number of available grid points ({GP}) is less than number of solvent to be inserted ({NS}).".format(GP=len(solv_grid_3D), NS=len(collected_solvent)), warn=True)
                        self.print_term("Reducing grid size (resolution) by 5% and recalculating grid occupancy.", warn=True)
                        self.print_term("    ", "Old: {GR}".format(GR=gridres), warn=True)
                        gridres = gridres*0.95
                        self.print_term("    ", "New: {GR}".format(GR=gridres), warn=True)
                        
                        if gridres < max(solvent_radii)*2:
                            self.print_term("New grid resolution ({GR}) is very small compared to largest solvent molecule. This will likely result in overlapping molecules. Consider changing the concentration.".format(GR=gridres), warn=True)

                self.print_term("Final number of 3D grid points available for solvent placement:", len(solv_grid_3D), spaces=2, verbose=2)

                ################################
                ### INSERT SOLVENT MOLECULES ###
                ################################

                def position_fixer(coords, ax):
                    '''
                    Checking if any beads are outside pbc and moves them if they are
                    If all beads are outside: Moved to other side of pbc
                    If at least 1 bead is inside: All beads are pushed slightly so that they are inside the pbc
                    '''
                    diffs = [coord - (np.sign(coord) * self.pbc_box[ax] / 2) for coord in coords]
                    checks = [np.sign(coords[i]) == np.sign(diffs[i]) for i in range(len(coords))]
                    if all(checks):
                        coords = [coord - (np.sign(coord) * self.pbc_box[ax]) for coord in coords]
                    elif any(checks):
                        coords = [coord - (np.sign(coord) * max(abs(min(diffs)), abs(max(diffs)))) for coord in coords]
                    return coords
                
                self.print_term("Inserting", len(collected_solvent), "solvent molecules into random grid points:", spaces=2, verbose=2)
                self.print_term("number of grid points:", len(solv_grid_3D), debug=True)
                self.print_term("number of solvent:    ", len(collected_solvent), debug=True)

                assert len(solv_grid_3D) >= len(collected_solvent), "\n".join([
                    "Too few 3D-grid points available for solvent insertion.",
                    "3D Grid points:                     " + str(len(solv_grid_3D)),
                    "Solvent molecules (including ions): " + str(len(collected_solvent)),
                    "Either reduce solvent or salt concentration 'solv_molarity' / 'salt_molarity' or make grid resolution value smaller 'gridres'",
                    "If you are using the 'solvation' command but are trying to flood the system, then try using the 'flooding' command instead",
                ])
                
                ### random.sample extracts k random elements from the list without duplicates
                random_grid_points = random.sample(solv_grid_3D, k = len(collected_solvent))

                ### Calculating kicks first might be faster for large number of solvent
                kick = solvation["kick"]
                kxs = [random.uniform(-kick, kick) for _ in range(len(collected_solvent))]
                kys = [random.uniform(-kick, kick) for _ in range(len(collected_solvent))]
                kzs = [random.uniform(-kick, kick) for _ in range(len(collected_solvent))]
                kicks = zip(kxs, kys, kzs)

                grid_solvated = []
                for counter, ((sname, stype, snbeads), (gx, gy, gz), (kx, ky, kz)) in enumerate(zip(collected_solvent, random_grid_points, kicks), 1):
                    '''
                    Finds a 3D-grid point for the solvent and places it there.
                    '''
                    if counter == 1 or counter % 25000 == 0:
                        self.print_term("Currently at solvent number:", counter, spaces=3, verbose=2)
                        # print("Currently at solvent number:", counter)
                    sdata = solvation[stype][sname]
                    sinfo = sdata.get_res_beads_info()
                    
                    sbeads       = [bead.bead for bead in sinfo]
                    sresnames    = [bead.resname for bead in sinfo]
                    sbeadcharges = [bead.charge for bead in sinfo]
                    scharge      = sdata.get_mol_charge()
                    sx, sy, sz   = sdata.get_coords("xyz")

                    ### Only makes sense to rotate if molecule contains more than 1 bead
                    if snbeads > 1:
                        rx, ry, rz   = random.uniform(0, 360), random.uniform(0, 360), random.uniform(0, 360)
                        for j, (x, y, z) in enumerate(zip(sx, sy, sz)):
                            sx[j], sy[j], sz[j] = self.rotate_point(x, y, z, rx, ry, rz)
                    
                    for j, (x, y, z) in enumerate(zip(sx, sy, sz)):
                        sx[j], sy[j], sz[j] = sx[j] + kx + gx, sy[j] + ky + gy, sz[j] + kz + gz
                    
                    ### Shouldn't be possible for beads to be placed outside box using the current algorithms
#                     sx = position_fixer(sx, 0)
#                     sy = position_fixer(sy, 1)
#                     sz = position_fixer(sz, 2)

                    ### Order that beads should be written in structure and topology files
                    ### if/elif faster than dictionary lookup due to most molecules (roughly 95-99%) being "solvent"
                    if stype == "solvent":
                        order = 1
                    elif stype == "pos_ions":
                        order = 2
                    elif stype == "neg_ions":
                        order = 3

                    grid_solvated.append({
                        "order":        order,
                        "name":         sname,
                        "resnames":     sresnames,
                        "type":         stype,
                        "beads":        sbeads,
                        "bead_charges": sbeadcharges,
                        "charge":       scharge,
                        "coords":       list(zip(sx, sy, sz)),
                    })
                self.print_term("Currently at solvent number:", counter, spaces=3, verbose=2)

                ### Orders the solvent key:value pairs to ensure topology-structure file match
                solvation["grid"] = sorted(grid_solvated, key=lambda g: (g["order"], g["name"]))

                solvation["solv_count"] = {}
                for grid_point_3D in solvation["grid"]:
                    '''
                    Counts all solvent for this specific solvation command
                    '''
                    if (grid_point_3D["name"], grid_point_3D["type"], grid_point_3D["charge"]) not in solvation["solv_count"].keys():
                        solvation["solv_count"][(grid_point_3D["name"], grid_point_3D["type"], grid_point_3D["charge"])] = 1
                    else:
                        solvation["solv_count"][(grid_point_3D["name"], grid_point_3D["type"], grid_point_3D["charge"])] += 1
                
                
                if solvation_i == 0:
                    SYS_solv_count = {}
                    SYS_solv_volume = solv_volume
                else:
                    SYS_solv_volume += solv_volume
                
                CMD_solv_count = {}
                for (sname, stype, scharge), count in solvation["solv_count"].items():
                    '''
                    Counts all solvent in the solvation command
                    '''
                    if (sname, stype, scharge) not in CMD_solv_count.keys():
                        CMD_solv_count[(sname, stype, scharge)] = count
                    else:
                        CMD_solv_count[(sname, stype, scharge)] += count
                    
                    if (sname, stype, scharge) not in SYS_solv_count.keys():
                        SYS_solv_count[(sname, stype, scharge)] = count
                    else:
                        SYS_solv_count[(sname, stype, scharge)] += count

                CMD_headers = ["Name", "Totals", "Total %", "Molarity(box)", "Molarity(free)", "Molarity(solvent)", "Charge"]
                CMD_names, CMD_charges, CMD_counts = list(zip(*[(sname, scharge, scount) for (sname, stype, scharge), scount in CMD_solv_count.items()]))
                CMD_tot = sum(CMD_counts)
                CMD_percentages = [round(count / CMD_tot * 100, 3) for count in CMD_counts]
                CMD_box_molarity = [round(scount / (N_A * box_volume), 3) for scount in CMD_counts]
                CMD_free_molarity = [round(scount / (N_A * box_free_volume), 3) for scount in CMD_counts]
                if solv_volume > 0:
                    CMD_solv_molarity = [round(scount / (N_A * solv_volume), 3) for scount in CMD_counts]
                else:
                    CMD_solv_molarity = ["NaN" for scount in CMD_counts]
                CMD_printer = [tuple(CMD_headers)] + list(zip(CMD_names, CMD_counts, CMD_percentages, CMD_box_molarity, CMD_free_molarity, CMD_solv_molarity, CMD_charges))
                CMD_max_lengths = [len(str(max(data, key = lambda d: len(str(d))))) for data in zip(*CMD_printer)]
                if self.extra_info and len(self.SOLVATIONS) > 1:
                    '''
                    Prints command-specific solvent information to the terminal
                    '''
                    self.print_term("", verbose=2)
                    self.print_term("Solvent data for the the specific solvation command", spaces=1, verbose=2)
                    for string in CMD_printer:
                        for_printer = []
                        for si, substring in enumerate(string):
                            for_printer.append('{0: <{L}}'.format(substring, L = CMD_max_lengths[si]))
                            if si+1 < len(string):
                                for_printer.append(":")
                        self.print_term(
                            *for_printer,
                            spaces=2,
                            verbose=2,
                        )

            if self.extra_info:
                '''
                Prints system-wide solvent information to the terminal
                '''
                self.print_term("", verbose=2)
                SYS_headers = ["Name", "Totals", "Total %", "Molarity(box)", "Molarity(free)", "Molarity(solvent)", "Charge"]
                SYS_names, SYS_charges, SYS_counts = list(zip(*[(sname, scharge, scount) for (sname, stype, scharge), scount in SYS_solv_count.items()]))
                SYS_tot = sum(SYS_counts)
                SYS_percentages = [round(count / SYS_tot * 100, 3) for count in SYS_counts]
                ### Lipid volume
                lipid_beads_for_cell_checker, solvent_box_lipids_charge = self.get_lipid_volume(
                    -self.pbcx/2, self.pbcx/2, -self.pbcy/2, self.pbcy/2, -self.pbcz/2, self.pbcz/2
                )
                SYS_leafs_volume = bead_volume * len(lipid_beads_for_cell_checker)

                ### Protein volume
                prot_beads_for_cell_checker, solvent_box_proteins_charge = self.get_protein_volume(
                    -self.pbcx/2, self.pbcx/2, -self.pbcy/2, self.pbcy/2, -self.pbcz/2, self.pbcz/2
                )
                SYS_prots_volume = bead_volume * len(prot_beads_for_cell_checker)
                SYS_box_volume = (self.pbcx * self.pbcy * self.pbcz) * 10**-27
                SYS_free_volume = SYS_box_volume - leafs_volume - prots_volume
                
                SYS_box_molarity = [round(scount / (N_A * SYS_box_volume), 3) for scount in SYS_counts]
                SYS_free_molarity = [round(scount / (N_A * SYS_free_volume), 3) for scount in SYS_counts]
                if SYS_solv_volume > 0:
                    SYS_solv_molarity = [round(scount / (N_A * SYS_solv_volume), 3) for scount in SYS_counts]
                else:
                    SYS_solv_molarity = ["NaN" for scount in SYS_counts]
                SYS_printer = [tuple(SYS_headers)] + list(zip(SYS_names, SYS_counts, SYS_percentages, SYS_box_molarity, SYS_free_molarity, SYS_solv_molarity, SYS_charges))
                SYS_max_lengths = [len(str(max(data, key = lambda d: len(str(d))))) for data in zip(*SYS_printer)]
                self.print_term("Solvent data for whole system", verbose=2)
                for string in SYS_printer:
                    for_printer = []
                    for si, substring in enumerate(string):
                        for_printer.append('{0: <{L}}'.format(substring, L = SYS_max_lengths[si]))
                        if si+1 < len(string):
                            for_printer.append(":")
                    self.print_term(
                        *for_printer,
                        spaces=2,
                        verbose=2,
                    )
            
            solvation_toc = time.time()
            solvation_time = round(solvation_toc - solvation_tic, 3)
            string = " ".join(["", "SOLVATION COMPLETE", ""])
            self.print_term("{string:-^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), spaces=0, verbose=1)
            string = " ".join(["", "(Time spent:", str(solvation_time), "[s])", ""])
            self.print_term("{string:^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), "\n", spaces=0, verbose=1)
