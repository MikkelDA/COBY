import ast

class cif_reader:
    def cif_reader(self, cif_file_dest, return_box_info = False):
        '''
        Reads a pdbx/mmCIF file and returns it as a dictionary
        '''
        key_value_types = {
            "atom_nr"  : int,
            "atom_name": str,
            "res_name" : str,
            "res_nr"   : int,
            "x"        : float, # [Å]
            "y"        : float, # [Å]
            "z"        : float, # [Å]
        }
        cif_to_datatype = {
            "id"          : "atom_nr",
            "auth_atom_id": "atom_name",
            "auth_comp_id": "res_name",
            "auth_seq_id" : "res_nr",
            "Cartn_x"     : "x",
            "Cartn_y"     : "y",
            "Cartn_z"     : "z",
        }
        dict_order = []
        processed_file = {}
        atom_nr = 0
        box_info = { # Default values in case not specified in file
            "x":     0,  # [Å]
            "y":     0,  # [Å]
            "z":     0,  # [Å]
            "alpha": 90, # [angle] (in degrees)
            "beta":  90, # [angle] (in degrees)
            "gamma": 90, # [angle] (in degrees)
        }
        with open(cif_file_dest, "r") as input_file:
            for line_nr, line in enumerate(input_file):
                if line.startswith("_cell"):
                    line_split = line.split()
                    if line_split[0] == "_cell.length_a":
                        box_info["x"] = float(line_split[1]) # [Å]
                    elif line_split[0] == "_cell.length_b":
                        box_info["y"] = float(line_split[1]) # [Å]
                    elif line_split[0] == "_cell.length_c":
                        box_info["z"] = float(line_split[1]) # [Å]
                    elif line_split[0] == "_cell.angle_alpha":
                        box_info["alpha"] = float(line_split[1]) # [angle] (in degrees)
                    elif line_split[0] == "_cell.angle_beta":
                        box_info["beta"] = float(line_split[1]) # [angle] (in degrees)
                    elif line_split[0] == "_cell.angle_gamma":
                        box_info["gamma"] = float(line_split[1]) # [angle] (in degrees)

                ### Ignore everything not atom-related data
                elif line.startswith("_atom_site"):
                    data_type = line.split(".")[-1]
                    ### All needed data types
                    if data_type in cif_to_datatype.keys():
                        dict_order.append(cif_to_datatype[data_type])
                    
                    ### Unused data types
                    else:
                        dict_order.append("PLACEHOLDER")

                elif any([line.split(" ") == string for string in ["ATOM", "HETATM"]]):
                    ### Adds only needed data types to atom_dict and converts to needed value type
                    atom_dict = {key: key_value_types[key](val) for key, val in zip(dict_order, line.split()) if key in key_value_types.keys()}
                    processed_file[(atom_nr, atom_dict["atom_nr"], atom_dict["res_nr"])] = atom_dict
                    atom_nr += 1
        
        if return_box_info:
            return (processed_file, box_info)
        else:
            return processed_file

