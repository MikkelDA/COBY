class pdb_reader:
    def pdb_reader(self, pdb_file_dest):
        '''
        Reads a pdb file and returns it as a dictionary
        '''
        # PCL = [6, 5, 4, 1, 3 + 1, 1, 4, 1 + 3, 8, 8, 8, 6, 6, 10, 2, 2]
        processed_file = {}
        atom_nr = 0
        with open(pdb_file_dest, "r") as input_file:
            for line_nr, line in enumerate(input_file):
                if any([line.startswith(string) for string in ["ATOM", "HETATM"]]):
                    atom_dict = {}
                    key_indexes = [
                        ("atom_nr"  , 6, 11, int),
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
        return processed_file

