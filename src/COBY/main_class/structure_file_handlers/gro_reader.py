import ast
import numpy as np
import math

class gro_reader:
    def gro_reader(self, gro_file_dest, return_box_info = False):
        '''
        Reads a gro file and returns it as a dictionary
        '''
        # GCL = [5, 5, 5, 5, 8, 8, 8, 8, 8, 8] # gro atom column lengths
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
        with open(gro_file_dest, "r") as input_file:
            for line_nr, line in enumerate(input_file):
                if line_nr == 0:
                    continue
                elif line_nr == 1:
                    tot_atoms = int(ast.literal_eval(line.strip()))
                    continue

                if tot_atoms == 0 and box_info:
                    break
                elif tot_atoms == 0 and not box_info:
                    gro_box_vectors = [float(val) for val in line.split()]

                    ### gro files allow for less than 9 values but all are needed for errorless conversion to pdb format so adding extra zeros.
                    if len(gro_box_vectors) < 9:
                        while len(gro_box_vectors) < 9:
                            gro_box_vectors.append(0)
                    
                    gro_box_vectors = [float(val) for val in gro_box_vectors]

                    ### ### Converting gro unit cell vectors to pdb unit cell lengths and angles
                    ### Create vectors
                    vector1 = np.array([gro_box_vectors[0], gro_box_vectors[3], gro_box_vectors[4]])
                    vector2 = np.array([gro_box_vectors[5], gro_box_vectors[1], gro_box_vectors[6]])
                    vector3 = np.array([gro_box_vectors[7], gro_box_vectors[8], gro_box_vectors[2]])
                    
                    ### Calculate angles
                    ### https://math.stackexchange.com/questions/361412/finding-the-angle-between-three-points
                    ### https://en.wikipedia.org/wiki/Lattice_constant#/media/File:UnitCell.png
                    alpha = math.degrees(np.arccos(np.dot(vector2, vector3) / (np.linalg.norm(vector2) * np.linalg.norm(vector3))))
                    beta  = math.degrees(np.arccos(np.dot(vector1, vector3) / (np.linalg.norm(vector1) * np.linalg.norm(vector3))))
                    gamma = math.degrees(np.arccos(np.dot(vector1, vector2) / (np.linalg.norm(vector1) * np.linalg.norm(vector2))))

                    ### Calculate lengths
                    a = np.linalg.norm(vector1)
                    b = np.linalg.norm(vector2)
                    c = np.linalg.norm(vector3)

                    box_info = {
                        "x":     a*10,  # From [nm] to [Å]
                        "y":     b*10,  # From [nm] to [Å]
                        "z":     c*10,  # From [nm] to [Å]
                        "alpha": alpha, # [angle] (in degrees)
                        "beta":  beta,  # [angle] (in degrees)
                        "gamma": gamma, # [angle] (in degrees)
                    }


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
        
        if return_box_info:
            return (processed_file, box_info)
        else:
            return processed_file


