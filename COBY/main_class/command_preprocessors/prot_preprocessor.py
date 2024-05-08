import ast
import copy

from COBY.structure_classes.PROTEIN_class import PROTEIN

class prot_preprocessor:
    def prot_preprocessor(self):
        '''
        Preprocesses protein commands for later ease of use
        '''
        if len(self.PROTEINS_cmds) != 0:
            self.print_term("Preprocessing protein requests", verbose=2)
            for cmd_nr, prot_cmd in enumerate(self.PROTEINS_cmds, 1):
                self.print_term("Starting protein:", cmd_nr, spaces=1, verbose=2)

                ### Defaults
                prot_dict = {
                    "cx": 0,
                    "cy": 0,
                    "cz": 0,
                    "rx": 0,
                    "ry": 0,
                    "rz": 0,
                    
                    "cen_method": ("cog",), # "axis" (axial center), "mean" (mean of all points), "bead:INT", "res:INT" or "point:x:y:z"
                    "lipids_inside": False, # [bool]
                    
                    "pbc_check": True,
                    ### Buffer used during solvations
                    "buffer": 1.32, # [Å] default = (vdw of regular beads) / 2
                    
                    "moleculetypes": [],
                    "charge": "top", # int/float, "top" or "lib"

                    "membrane_border": False,
                }
                
                ### ### Check protein command
                for cmd in prot_cmd.split():
                    sub_cmd = cmd.split(":")

                    ### Read pdb/gro file
                    if sub_cmd[0].lower() == "file":
                        ### Find file name
                        file_name = sub_cmd[1]
                        
                        ### First check if file ends with .pdb or .gro
                        if file_name.endswith(".pdb"):
                            prot_dict["beads"] = self.pdb_reader(file_name)
                        elif file_name.endswith(".gro"):
                            prot_dict["beads"] = self.gro_reader(file_name)
                        ### If neither then check if it contains .pdb or .gro (e.g. #prot_file.pdb.1#)
                        elif ".pdb" in file_name.lower():
                            prot_dict["beads"] = self.pdb_reader(file_name)
                        elif ".gro" in file_name.lower():
                            prot_dict["beads"] = self.gro_reader(file_name)
                        ### Finally assume .pdb format if neither .pdb nor .gro is found
                        else:
                            prot_dict["beads"] = self.pdb_reader(file_name)
                            
                    ### Center methods "cog" / "mean_of_beads" (mean of all points) (default), "axis" / "mean_of_extremes" (axial center), "bead:INT", "res:INT" or "point:x:y:z"
                    elif sub_cmd[0].lower() in ["cen_method", "centering_method"]:
                        if sub_cmd[1].lower() in ["cog", "mean_of_beads"] + ["axis", "mean_of_extremes"]:
                            prot_dict["cen_method"] = (sub_cmd[1].lower(),)
                        elif any([sub_cmd[1].lower().startswith(value) for value in ["bead", "res"]]) and any([sub_cmd[1].lower().endswith(value) for value in ["bead", "res", "cog", "mean_of_beads", "axis", "mean_of_extremes"]]):
                            prot_dict["cen_method"] = (sub_cmd[1].lower(), sub_cmd[2:])
                        elif sub_cmd[1].lower() in ["point"]:
                            prot_dict["cen_method"] = (sub_cmd[1].lower(), ast.literal_eval(sub_cmd[2])*10, ast.literal_eval(sub_cmd[3])*10, ast.literal_eval(sub_cmd[4])*10)

                    ### True/False whether or not lipds may be placed within protein area
                    elif sub_cmd[0].lower() == "lipids_inside":
                        prot_dict["lipids_inside"] = ast.literal_eval(str(sub_cmd[1]))

                    ### True/False whether to check for pbc crossing
                    elif sub_cmd[0].lower() == "pbc_check":
                        prot_dict["pbc_check"] = ast.literal_eval(str(sub_cmd[1]))

                    ### True/False whether to check for pbc crossing
                    elif sub_cmd[0].lower() == "buffer":
                        prot_dict["buffer"] = ast.literal_eval(str(sub_cmd[1]))

                    ### Sets the topology name(s) present in the structure file
                    elif sub_cmd[0].lower() in ["moleculetypes", "moleculetype"]:
                        moli = 1
                        while moli < len(sub_cmd):
                            if len(sub_cmd[moli:]) > 1:
                                isnumber, isinteger = self.is_number(sub_cmd[moli+1])
                                if isnumber:
                                    prot_dict["moleculetypes"].extend([sub_cmd[moli]] * self.get_number_from_string(sub_cmd[moli+1]))
                                    moli += 2
                                else:
                                    prot_dict["moleculetypes"].append(sub_cmd[moli])
                                    moli += 1
                            else:
                                prot_dict["moleculetypes"].append(sub_cmd[moli])
                                moli += 1

                    ### Determine charge of protein int/float or "top"
                    elif sub_cmd[0].lower() == "charge":
                        assert prot_dict["charge"] in ["topology", "top", "library", "lib"] or self.get_number_from_string(sub_cmd[1]) is not False, "\n".join([
                            "The 'charge' subcommand must be given either 'topology'/'top', 'library'/'lib' or a [float]/[int] as value. You have given: " + str(sub_cmd[1])
                        ])
                        if sub_cmd[1] in ["topology", "top"]:
                            prot_dict["charge"] = "top"
                        elif sub_cmd[1] in ["library", "lib"]:
                            prot_dict["charge"] = "lib"
                        else:
                            prot_dict["charge"] = ast.literal_eval(sub_cmd[1])

                    ### xyz axis translations [nm]
                    elif sub_cmd[0].lower() == "center" and len(sub_cmd[1:]) == 3:
                        prot_dict["cx"], prot_dict["cy"], prot_dict["cz"] = [ast.literal_eval(str(val)) * 10 for val in sub_cmd[1:]]

                    ### xyz axis translations [nm]
                    elif sub_cmd[0].lower() in ["cx", "cy", "cz"]:
                        prot_dict[sub_cmd[0].lower()] = ast.literal_eval(str(sub_cmd[1])) * 10

                    ### xyz axis rotations [degrees]
                    elif sub_cmd[0].lower() in ["rx", "ry", "rz"]:
                        prot_dict[sub_cmd[0].lower()] = ast.literal_eval(str(sub_cmd[1]))

                    ### Sets whether protein should constitute a membrane border
                    elif sub_cmd[0].lower() == "membrane_border":
                        prot_dict["membrane_border"] = ast.literal_eval(sub_cmd[1])


                    ### Errors out if unknown subcommand used, and prints the subcommand to console
                    else:
                        assert False, "Unknown subcommand given to '-prot'. The subcommand is: '" + str(cmd) + "'"
                
                ### Post-preprocessing (topology and charge determination)
                if isinstance(prot_dict["charge"], (float, int)):
                    ### Sets charge manually
                    ### Evenly distributes charges accross all beads as there is no way of knowing where they are located based on a single charge value
                    protein_bead_charges = [prot_dict["charge"]/len(prot_dict["beads"]) for _ in range(len(prot_dict["beads"]))]
                    
                elif prot_dict["charge"] == "lib":
                    ### Finds charge information from prot_defs charge dictionary
                    protein_bead_charges = []
                    for (bead_i, bead_nr, res_nr), values in prot_dict["beads"].items():
                        params   = self.prot_params or self.sys_params
                        res_name  = values["res_name"]
                        atom_name = values["atom_name"]
                        assert res_name in self.prot_defs[params]["charges"], "Residue name '{resname}' not in 'prot_defs' charges for system parameters '{params}'.".format(resname=res_name, params=params)
                        assert atom_name in self.prot_defs[params]["charges"][res_name], "Residue name '{atomname}' not in 'prot_defs' charges for system parameters '{params}'.".format(atomname=atom_name, params=params)
                        protein_bead_charges.append(self.prot_defs[params]["charges"][res_name][atom_name])

                
                elif prot_dict["charge"] == "top":
                    ### Finds charge information in topology files
                    assert len(prot_dict["moleculetypes"]) > 0, "\n".join([
                        "No 'moleculetypes' given. At least one 'molecutype' must be given if charges are to be determined from topology.",
                        "If you want protein charges to be automatically determined, then use the 'charge:lib' subcommand.",
                        "If you want to set the total charge of the protein, then use the 'charge:[float]/[int]' subcommand.",
                    ])
                    ### Finds charge data from topology files
                    protein_bead_charges = []
                    erase_charges = False
                    for moleculetype in prot_dict["moleculetypes"]:
                        assert moleculetype in self.itp_moleculetypes.keys(), "Moleculetype '{moleculetype}' not in topology.".format(moleculetype=moleculetype)
                        bead_charges = list(map(
                            self.get_number_from_string,
                            self.itp_moleculetypes[moleculetype].topology_types_in_molecule["atoms"].get_value_type("charge")
                        ))
                        ### Assert for len(bead_charges) == n_atoms done further below
                        protein_bead_charges.extend(bead_charges)
                

                prot_dict["moleculetypes"] = prot_dict["moleculetypes"] or ["_".join(["PROT", str(cmd_nr)])]
                
                ### Mapping charge values to int/float data type
                protein_bead_charges = list(map(self.get_number_from_string, protein_bead_charges))

                assert all(charge is not False for charge in protein_bead_charges), "One of the charges found for '{moleculetypes}' is not a number".format(moleculetypes=prot_dict["moleculetypes"])
                
                ### Converting data into protein class
                prot_dict["protein"] = PROTEIN(molname = prot_dict["moleculetypes"])

                assert len(protein_bead_charges) == len(prot_dict["beads"]), "\n".join([
                    "Mismatch found between number of beads in structure and number of beads in topology for \"{name}\" during preprocessing of protein command".format(name=prot_dict["moleculetypes"]),
                    "Number of protein bead charges: {charges}".format(charges = len(protein_bead_charges)),
                    "Number of protein beads: {charges}".format(charges = len(prot_dict["beads"])),

                ])

                cur_res = False
                for (key, vals), charge in zip(prot_dict["beads"].items(), protein_bead_charges):
                    if vals["res_nr"] != cur_res:
                        prot_dict["protein"].add_res(vals["res_name"], vals["res_nr"])
                        cur_res = vals["res_nr"]
                    
                    prot_dict["protein"].add_bead(
                        vals["atom_name"],
                        vals["atom_nr"],
                        vals["x"],
                        vals["y"],
                        vals["z"],
                        charge=charge,
                    )
                
                self.system_charge += prot_dict["protein"].get_mol_charge()
                
                self.PROTEINS[cmd_nr] = prot_dict.copy()
            self.print_term("Number of molecule insertions preprocessed:", len(self.PROTEINS), spaces=1, verbose=2)
