import ast

class cif_reader:
    def cif_reader(self, cif_file_dest):
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
        with open(cif_file_dest, "r") as input_file:
            for line_nr, line in enumerate(input_file):
                ### Ignore everything not atom-related data
                if line.startswith("_atom_site"):
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
        return processed_file

