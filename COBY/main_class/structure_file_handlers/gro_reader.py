import ast

class gro_reader:
    def gro_reader(self, gro_file_dest):
        '''
        Reads a gro file and returns it as a dictionary
        '''
        # GCL = [5, 5, 5, 5, 8, 8, 8, 8, 8, 8] # gro atom column lengths
        processed_file = {}
        atom_nr = 0
        with open(gro_file_dest, "r") as input_file:
            file_length = len(input_file.readlines())
        with open(gro_file_dest, "r") as input_file:
            for line_nr, line in enumerate(input_file):
                if line_nr == 0:
                    continue
                elif line_nr == 1:
                    tot_atoms = int(ast.literal_eval(line.strip()))
                    continue
                if tot_atoms == 0:
                    break
                else:
                    if len(line):
                        tot_atoms -= 1
                        atom_dict = {}
                        key_indexes = [
                            ("res_nr"   , 0, 5, int),
                            ("res_name" , 5, 10, str),
                            ("atom_name", 10, 15, str),
                            ("atom_nr"  , 15, 20, int),
                            ("x"        , 20, 28, float), # [nm]
                            ("y"        , 28, 36, float), # [nm]
                            ("z"        , 36, 44, float), # [nm]
                        ]
                        for key, i1, i2, func in key_indexes:
                            if len(line) >= i2:
                                atom_dict[key] = func(line[i1:i2].replace(" ",""))
                                if key in ["x", "y", "z"]: # Convert [nm] to [Å]
                                    atom_dict[key] *= 10
                        processed_file[(atom_nr, atom_dict["atom_nr"], atom_dict["res_nr"])] = atom_dict
                        atom_nr += 1
        return processed_file


