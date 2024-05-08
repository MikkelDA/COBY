from COBY.structure_classes.LIPID_class import LIPID

class lipid_scaffolds_preprocessor:
    def lipid_scaffolds_preprocessor(self):
        '''
        Preprocesses lipid scaffolds, by rearranging the data into a class
        '''
        self.print_term("Preprocessing lipid scaffolds", verbose=2)
        ### Loop over generel lipid types (e.g. phospholipids, sterols, etc.)
        for (lipid_type, params), lipid_type_dict in self.lipid_scaffolds.items():
            if params not in self.lipid_dict.keys():
                self.lipid_dict[params] = {}
            x, y, z = lipid_type_dict["x"], lipid_type_dict["y"], lipid_type_dict["z"]

            ### bead distance adjustment values for current general lipid type
            if "bd" in lipid_type_dict.keys():
                bdx, bdy, bdz = lipid_type_dict["bd"]
            else:
                bdx, bdy, bdz = 0.25, 0.25, 0.3
            
            if "charges" in lipid_type_dict.keys() and len(lipid_type_dict["charges"]) > 0:
                ### Checks if nested list/tuple and converts to nested tuples
                if type(lipid_type_dict["charges"][0]) in [list, tuple]:
                    lipid_type_charges = tuple([tuple(details) for details in lipid_type_dict["charges"]])
                else:
                    ### Converts to nested tuples
                    if type(lipid_details) == list:
                        lipid_type_charges = (tuple(lipid_type_dict["charges"]),)
            else:
                lipid_type_charges = 0
            
            ### Checks what specifiers have been given for each lipid name
            ### Used to check if lipid-specific charges have been given
            dict_keys, dict_vals = zip(*lipid_type_dict["lipids"].items())
            lipid_names, data_specifiers = zip(*dict_keys)
            lipid_name_data_specifiers = {lipid_name: {} for lipid_name in lipid_names}
            for (lipid_name, data_specifier), lipid_details in lipid_type_dict["lipids"].items():
                if data_specifier not in lipid_name_data_specifiers[lipid_name]:
                    lipid_name_data_specifiers[lipid_name][data_specifier] = lipid_details
            
            ### Removing duplicates from lipid_names list while keeping the order
            lipid_names = list(dict.fromkeys(lipid_names))
            
            for lipid_name in lipid_names:
                
                ### Creates new LIPID instance if it does not already exist.
                if lipid_name not in self.lipid_dict[params].keys():
                    self.lipid_dict[params][lipid_name] = LIPID(molname = lipid_name, moleculetype = lipid_name)
                else:
                    self.print_term(
                        "The name \"{name}\" has been detected multiple times for lipids within the same parameter library.".format(name=lipid_name),
                        "If multiple lipids have the same name then only the first is loaded into the library.",
                        warn=True,
                    )
                    continue
                
                ### Checks if lipid-specific charges are given and joins them with lipid_type charges
                if "charges" in lipid_name_data_specifiers[lipid_name]:
                    lipid_details = lipid_type_dict["lipids"][(lipid_name, "charges")]
                    ### Checks if nested list/tuple and converts to nested tuples
                    if type(lipid_details[0]) in [list, tuple]:
                        lipid_details = tuple([tuple(details) for details in lipid_details])
                    else:
                        ### Converts to nested tuples
                        if type(lipid_details) == list:
                            lipid_details = (tuple(lipid_details),)
                    
                    if lipid_type_charges != 0:
                        specifc_charge_names = list(zip(*lipid_details))[0]
                        lipid_charge_dict = lipid_details + tuple([(key, val) for key, val in lipid_type_charges if key not in specifc_charge_names])
                    else:
                        lipid_charge_dict = lipid_details
                else:
                    lipid_charge_dict = lipid_type_charges
            
                if "beads" in lipid_name_data_specifiers[lipid_name]:
                    lipid_details = lipid_type_dict["lipids"][(lipid_name, "beads")]
                    beads = list(filter(None, lipid_details.split(" ")))
                    ### Remove coordinate and bead indexes if no bead is assigned
                    ### x and z are centered in the plane using set_coords_to_center(AXs="xy") method further down
                    lx = [xi * bdx * 10 for xi, bead in zip(x, beads) if bead != "-"]
                    ly = [yi * bdy * 10 for yi, bead in zip(y, beads) if bead != "-"]
                    minz = min([zi for zi, bead in zip(z, beads) if bead != "-"])
                    inter_leaflet_buffer = 1.5
                    lz = [(zi - minz) * bdz * 10 + inter_leaflet_buffer for zi, bead in zip(z, beads) if bead != "-"]
                    beads = [bead for bead in beads if bead != "-"]
                    beadnrs = list([i for i in range(len(beads))])
                    
                    if lipid_charge_dict:
                        charges = []
                        charge_beads, charge_vals = zip(*lipid_charge_dict)
                        for bead in beads:
                            if bead in charge_beads:
                                charge_index = charge_beads.index(bead)
                                charges.append(charge_vals[charge_index])
                            else:
                                charges.append(0)
                    else:
                        charges = [0 for _ in range(len(beads))]
                    
                    self.lipid_dict[params][lipid_name].add_res_and_beads(
                        lipid_name,
                        beads   = beads,
                        beadnrs = beadnrs,
                        xs      = lx,
                        ys      = ly,
                        zs      = lz,
                        charges = charges,
                    )
                    self.lipid_dict[params][lipid_name].set_coords_to_center(centering = "axis", AXs = "xy")
                    
        tot_lipids = sum([len(vals) for vals in self.lipid_dict.values()])
        
        self.print_term("Number of lipids preprocessed:", tot_lipids, "\n", spaces=1, verbose=2)
    