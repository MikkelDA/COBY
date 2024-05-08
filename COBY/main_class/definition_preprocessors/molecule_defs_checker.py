from COBY.main_class.definition_preprocessors.molecule_beads_checker import molecule_beads_checker
from COBY.main_class.general_tools.rotation_matrix_from_vectors import rotation_matrix_from_vectors
import numpy as np
import math

class molecule_defs_checker(molecule_beads_checker, rotation_matrix_from_vectors):
    def molecule_defs_checker(self, mol_class, cur_name, cur_dict, type_of_molecule, structure_origin):
        charge_dict = {}
        tot_charge = False
        
        if structure_origin == "defs":
            if "moleculetype" in cur_dict:
                mol_class.moleculetype = cur_dict["moleculetype"]
            else:
                mol_class.moleculetype = cur_name
        
        elif structure_origin == "import":
            if "moleculetype" in cur_dict:
                mol_class.moleculetype = cur_dict["moleculetype"]

        self.print_term("", debug=True) ### Spacer for terminal and log file
        self.print_term("Current name", cur_name, debug=True)
        self.print_term("    ", "Current dictionary", cur_dict, debug=True)
        
        if "bd" in cur_dict:
            bd = cur_dict["bd"]
        else:
            bd = (1, 1, 1)
        
        if type_of_molecule == "solvent":
            if "mapping_ratio" in cur_dict.keys():
                mol_class.mapping_ratio_set(cur_dict["mapping_ratio"])

            if "density" in cur_dict.keys():
                mol_class.density_set(cur_dict["density"])

            if "molar_mass" in cur_dict.keys():
                mol_class.molar_mass_set(cur_dict["molar_mass"])
        
        
        ### Ensures molecule has the minimally required information
        assert ("residues" in cur_dict.keys()) or (all(val in cur_dict.keys() for val in ["beads", "x", "y", "z"])), "\n".join([
            " ".join([
                "A solute/solvent molecule must either contain a list of residues dictionaries each containing 'resname', 'beads', 'x', 'y' and 'z'",
                "or 'beads', 'x', 'y' and 'z' must be present within the molecule dictionary, in which case it will be treated as a single-residue molecule",
            ]),
            "    "+"Name of current molecule: " + cur_name,
            "    "+"'resname' is the name of the residue",
            "    "+"'beads' is a list or tuple of all bead names in the residue",
            "    "+"'x' is a list or tuple of the x-coordinates for all beads in the residue",
            "    "+"'y' is a list or tuple of the y-coordinates for all beads in the residue",
            "    "+"'z' is a list or tuple of the z-coordinates for all beads in the residue",
            "The i'th value in 'beads', 'x', 'y' and 'z' is connected to each other",
            # "Full solute/solvent dictionary: "+str([str(key)+": "+str(vals) for key, vals in cur_dict.items()]),
            "Full solute/solvent dictionary: "+str(cur_dict),
        ])

        ### Specifies residue for single-residue solutes that do not have it specified already
        ### Charges are handled later
        if "residues" not in cur_dict.keys():
            residues_dict = [{
                    "resname": cur_name,
                    "beads": cur_dict["beads"],
                    "x": cur_dict["x"],
                    "y": cur_dict["y"],
                    "z": cur_dict["z"],
            }]

        elif type(cur_dict["residues"]) == dict:
            residues_dict = [cur_dict["residues"]]
        
        elif isinstance(cur_dict["residues"], (list, tuple)):
            assert all(val in residue for val in ["resname", "beads", "x", "y", "z"] for residue in cur_dict["residues"]), "\n".join([
                "Residues for solute/solvent molecule are specified but do not contain all required values",
                "    "+"Name of current molecule: " + cur_name,
                "A residue must contain 'resname', 'beads', 'x', 'y' and 'z'",
                "    "+"'resname' is the name of the residue",
                "    "+"'beads' is a list or tuple of all bead names in the residue",
                "    "+"'x' is a list or tuple of the x-coordinates for all beads in the residue",
                "    "+"'y' is a list or tuple of the y-coordinates for all beads in the residue",
                "    "+"'z' is a list or tuple of the z-coordinates for all beads in the residue",
                "The i'th value in 'beads', 'x', 'y' and 'z' is connected to each other",
            ])
            residues_dict = cur_dict["residues"]
        
        ### Molecule-wide charges
        if "charges" in cur_dict.keys():

            ### Ensures tuple of tuples
            if type(cur_dict["charges"][0]) in [str, int, float]:
                cur_dict["charges"] = (cur_dict["charges"],)

            for charge_values in cur_dict["charges"]:
                
                ### ### At this point "residues" has been added to the dictionary even if it wasn't originally
                ### Multiple residues (Residues specification is required)
                if len(residues_dict) > 1:
                    assert len(charge_values) == 3, "\n".join([
                        "Bead charges specified outside residues when more than 1 residue is present must contain exactly three numbers",
                        "The numbers must be (residue nr [int], bead nr [int], charge [float] or [int])",
                        "Name of current molecule: " + cur_name,
                        "The problematic residue-bead-charge designation: " + str(charge_values)
                    ])

                    res_nr  = self.get_number_from_string(charge_values[0])
                    bead_nr = self.get_number_from_string(charge_values[1])
                    charge  = self.get_number_from_string(charge_values[2])

                    ### Set residue number to 0 as there is only one residue
                    assert all(val is not False for val in [res_nr, bead_nr, charge]), "\n".join([
                        "An incorrect value has been given for residue bead charges",
                        "The first value must be the residue number [int] (not the name)",
                        "The second value must be the bead number [int] (not the name)",
                        "The third value must be the charge value [float] or [int]",
                        "    "+"Name of current molecule: " + cur_name,
                        "    "+"The problematic residue-bead-charge designation: " + str(charge_values)
                    ])
                
                ### Single residue (Residues specification is optional)
                else:
                    assert len(charge_values) in [2, 3], "\n".join([
                        "Bead charges specified outside residues when exactly 1 residue is present must contain either two or three numbers",
                        "    "+"Two number version: The numbers must be (bead nr [int], charge [float] or [int])",
                        "    "+"Three number version: The numbers must be (residue nr [int], bead nr [int], charge [float] or [int])",
                        "Name of current molecule: " + cur_name,
                        "Name and number of current residue: " + res_dict["resname"] + ", " + str(res_nr),
                        "The problematic charge designation: " + str(charge_values)
                    ])

                    if len(charge_values) == 2:
                        res_nr  = 0
                        bead_nr = self.get_number_from_string(charge_values[0])
                        charge  = self.get_number_from_string(charge_values[1])

                        ### Set residue number to 0 as there is only one residue
                        assert all(val is not False for val in [res_nr, bead_nr, charge]), "\n".join([
                            "An incorrect value has been given for residue bead charges",
                            "The first value must be the residue number [int] (not the name)",
                            "The second value must be the bead number [int] (not the name)",
                            "The third value must be the charge value [float] or [int]",
                            "    "+"Name of current molecule: " + cur_name,
                            "    "+"The problematic residue-bead-charge designation: " + str(charge_values)
                        ])
                    
                    elif len(charge_values) == 3:
                        res_nr  = self.get_number_from_string(charge_values[0])
                        bead_nr = self.get_number_from_string(charge_values[1])
                        charge  = self.get_number_from_string(charge_values[2])

                        ### Set residue number to 0 as there is only one residue
                        assert all(val is not False for val in [res_nr, bead_nr, charge]), "\n".join([
                            "An incorrect value has been given for residue bead charges",
                            "The first value must be the residue number [int] (not the name)",
                            "The second value must be the bead number [int] (not the name)",
                            "The third value must be the charge value [float] or [int]",
                            "    "+"Name of current molecule: " + cur_name,
                            "    "+"The problematic residue-bead-charge designation: " + str(charge_values)
                        ])
                    
                if res_nr not in charge_dict.keys():
                    charge_dict[res_nr] = {}
                charge_dict[res_nr][bead_nr] = charge
        
        self.print_term("    ", "Residue dictionary", residues_dict, debug=True)

        ### Residue-specific charges
        ### No need to check for residues as they are created earlier in this script
        for res_nr, res_dict in enumerate(residues_dict):
            if "charges" in res_dict:
                if res_nr not in charge_dict.keys():
                    charge_dict[res_nr] = {}
                
                ### Ensures tuple of tuples
                if type(res_dict["charges"][0]) in [str, int, float]:
                    res_dict["charges"] = (res_dict["charges"],)

                for bead_nr_charge in res_dict["charges"]:
                    assert len(bead_nr_charge) == 2, "\n".join([
                        "Bead charges specified within residues must contain exactly two numbers",
                        "Name of current molecule: " + cur_name,
                        "Name and number of current residue: " + res_dict["resname"] + ", " + str(res_nr),
                        "The problematic bead-charge designation: " + str(bead_nr_charge)
                    ])
                    bead_nr = self.get_number_from_string(bead_nr_charge[0])
                    charge = self.get_number_from_string(bead_nr_charge[1])
                    assert all(val is not False for val in [bead_nr, charge]), "\n".join([
                        "An incorrect value has been given for within-residue charges",
                        "The first value must be the bead number [int] (not the name)",
                        "The second value must be the charge value [float] or [int]",
                        "    "+"Name of current molecule: " + cur_name,
                        "    "+"Name and number of current residue: " + res_dict["resname"] + ", " + str(res_nr),
                        "    "+"The problematic bead-charge designation: " + str(bead_nr_charge)
                    ])
                    charge_dict[res_nr][bead_nr] = charge

        ### Total charge.
        ### Spreads charge across all beads. Should only be used if no bead-specfic charges are given.
        ### Ideally only use this for single bead solutes/solvents where the bead is obvious
        if "charge" in cur_dict.keys():
            number = self.get_number_from_string(cur_dict["charge"])
            if number is not False:
                tot_charge = number
        
        self.print_term("    ", "Charge dictionary", charge_dict, debug=True)
        self.print_term("    ", "tot_charge       ", tot_charge, debug=True)
        self.print_term("    ", "moleculetype     ", mol_class.moleculetype, debug=True)

        ### Rotates a lipid such that it is vertically aligned based on the designated upwards and downwards pointing beads
        if type_of_molecule == "lipid" and "upbeads" in cur_dict.keys() and "downbeads" in cur_dict.keys():
            up_coords = []
            down_coords = []
            for beadname, coords_list in [("upbeads", up_coords), ("downbeads", down_coords)]:
                for res_nr, bead_nr in cur_dict[beadname]:
                    x = round(residues_dict[res_nr]["x"][bead_nr], 4)
                    y = round(residues_dict[res_nr]["y"][bead_nr], 4)
                    z = round(residues_dict[res_nr]["z"][bead_nr], 4)
                    coords_list.append((x, y, z))
            
            up_coords_array        = np.array(up_coords)
            down_coords_array      = np.array(down_coords)
            up_coords_mean         = np.mean(up_coords_array, axis=0)
            down_coords_mean       = np.mean(down_coords_array, axis=0)
            original_vector        = up_coords_mean - down_coords_mean
            original_vector_length = math.sqrt(original_vector[0]**2 + original_vector[1]**2+original_vector[2]**2)
            alignment_vector       = [0, 0, original_vector_length]
            rotation_matrix        = self.rotation_matrix_from_vectors(original_vector, alignment_vector)
            for res_nr, res_dict in enumerate(residues_dict, 0):
                xs, ys, zs = res_dict["x"], res_dict["y"], res_dict["z"]
                new_xs, new_ys, new_zs = [], [], []
                for x, y, z in zip(xs, ys, zs):
                    coord = np.array((x, y, z))
                    new_coord = np.dot(coord, rotation_matrix.T)
                    new_xs.append(new_coord[0])
                    new_ys.append(new_coord[1])
                    new_zs.append(new_coord[2])
                res_dict["x"], res_dict["y"], res_dict["z"] = new_xs, new_ys, new_zs
        
        ### Dictionaries without "residues" key have had the beads and positions converted to "residues" key:val pair
        ### Thus "residues" will always be present
        molbeadnrs = []
        total_number_of_beads = sum([len(res_dict["beads"]) for res_dict in residues_dict])
        for res_nr, res_dict in enumerate(residues_dict, 0):
            resname = res_dict["resname"]
            beads, xs, ys, zs = self.molecule_beads_checker(res_dict, cur_name, resname, bd, type_of_molecule)

            if molbeadnrs == []:
                beadnrs = [i for i in range(len(beads))]
                molbeadnrs = beadnrs[:]
            else:
                beadnrs = [i for i in range(max(molbeadnrs)+1, len(molbeadnrs)+len(beads))]
            
            ### Total charge distributed across beads
            if tot_charge is not False:
                charges = [tot_charge/total_number_of_beads for _ in range(len(beads))]
            
            ### Specific charges for specific beads specified in library or in import command
            elif len(charge_dict) > 0:
                charges = []
                for bead_nr, bead in enumerate(beads, 0):
                    if res_nr in charge_dict:
                        if bead_nr in charge_dict[res_nr]:
                            charges.append(charge_dict[res_nr][bead_nr])
                        else:
                            charges.append(0)
                    else:
                        charges.append(0)
            
            ### Finally set charges to zero if none of the above are triggered
            else:
                charges = [0 for _ in range(len(beads))]

            mol_class.add_res_and_beads(
                resname = resname, beads = beads, beadnrs = beadnrs,
                xs = xs, ys = ys, zs = zs,
                charges = charges,
            )

            self.print_term("    ", "resname", resname, debug=True)
            self.print_term("    ", "    ", "beads  ", beads, debug=True)
            self.print_term("    ", "    ", "beadnrs", beadnrs, debug=True)
            self.print_term("    ", "    ", "xs     ", xs, debug=True)
            self.print_term("    ", "    ", "ys     ", ys, debug=True)
            self.print_term("    ", "    ", "zs     ", zs, debug=True)
            self.print_term("    ", "    ", "charges", charges, debug=True)
        
        self.print_term("", debug=True) ### Spacer for terminal and log file

        
        if type_of_molecule == "solvent":
            mol_class.set_coords_to_center(centering = "axis")
        elif type_of_molecule == "lipid":
            mol_class.set_coords_to_center(centering = "axis", AXs = "xy")
    