class pdb_reader:
    def pdb_reader(self, pdb_file_dest, return_box_info = False):
        '''
        Reads a pdb file and returns it as a dictionary
        '''
        # PCL = [6, 5, 4, 1, 3 + 1, 1, 4, 1 + 3, 8, 8, 8, 6, 6, 10, 2, 2]
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
        with open(pdb_file_dest, "r") as input_file:
            for line_nr, line in enumerate(input_file):
                if line.startswith("CRYST1"):
                    box_info = {
                        "x":     float(line[6:15]),  # [Å]
                        "y":     float(line[15:24]), # [Å]
                        "z":     float(line[24:33]), # [Å]
                        "alpha": float(line[33:40]), # [angle] (in degrees)
                        "beta":  float(line[40:47]), # [angle] (in degrees)
                        "gamma": float(line[47:54]), # [angle] (in degrees)
                    }
                
                if any([line.startswith(string) for string in ["ATOM", "HETATM"]]):
                    atom_dict = {}
                    key_indexes = [
                        ("atom_nr"  , 6, 11,  int),
                        ("atom_name", 11, 16, str),
                        ("res_name" , 16, 21, str),
                        ("res_nr"   , 22, 26, int),
                        ("x"        , 30, 38, float), # [Å]
                        ("y"        , 38, 46, float), # [Å]
                        ("z"        , 46, 54, float), # [Å]
                    ]
                    for key, i1, i2, func in key_indexes:
                        if len(line) >= i2:
                            atom_dict[key] = func(line[i1:i2].replace(" ",""))
                    processed_file[(atom_nr, atom_dict["atom_nr"], atom_dict["res_nr"])] = atom_dict
                    atom_nr += 1
        
        if return_box_info:
            return (processed_file, box_info)
        else:
            return processed_file
