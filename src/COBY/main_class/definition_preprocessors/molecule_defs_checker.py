from COBY.main_class.definition_preprocessors.molecule_beads_checker import molecule_beads_checker
from COBY.main_class.general_tools.rotation_matrix_from_vectors import rotation_matrix_from_vectors
from COBY.main_class.general_tools.tags_checker import tags_checker
import numpy as np
import math
import copy

class molecule_defs_checker(molecule_beads_checker, rotation_matrix_from_vectors, tags_checker):
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
        residues_list = []
        if "residues" in cur_dict.keys() and isinstance(cur_dict["residues"], (list, tuple)):
            assert all(val in residue for val in ["resname", "beads"] for residue in cur_dict["residues"]) or all(val in residue for val in ["resname", "names", "x", "y", "z"] for residue in cur_dict["residues"]), "\n".join([
                "Residues for molecule are specified but do not contain all required values or do not all use the same format.",
                "    "+"Name of current molecule: " + cur_name,
                "Residues can be written in two different ways:",
                "    "+"Version 1: A residue must contain 'resname' and 'beads'.",
                "    "+"    "+"'resname' is the name of the residue.",
                "    "+"    "+"'beads' is a list of dictionaries each containing information for a specific bead. The dictionaries must contain the following keys.",
                "    "+"    "+"    "+"'name' Bead name.",
                "    "+"    "+"    "+"'x' is the bead x-coordinate.",
                "    "+"    "+"    "+"'y' is the bead y-coordinate.",
                "    "+"    "+"    "+"'z' is the bead z-coordinate.",
                "    "+"    "+"    "+"'charge' is the charge of the bead (optional, assumed to be zero if not specified).",
                "    "+"Version 2: A residue must contain 'resname', 'names', 'x', 'y' and 'z' keys.",
                "    "+"    "+"'resname' is the name of the residue.",
                "    "+"    "+"'names' is a list or tuple of all bead names in the residue.",
                "    "+"    "+"'x' is a list or tuple of the x-coordinates for all names in the residue.",
                "    "+"    "+"'y' is a list or tuple of the y-coordinates for all names in the residue.",
                "    "+"    "+"'z' is a list or tuple of the z-coordinates for all names in the residue.",
                "    "+"    "+"The i'th value in 'names', 'x', 'y' and 'z' is connected to each other.",
            ])
            for residue_dict in cur_dict["residues"]:
                ### If separate lists/tuples data type
                if "names" in residue_dict:
                    residues_list.append(copy.deepcopy(residue_dict))
                
                ### If dictionary data type
                elif "beads" in residue_dict:
                    names, x, y, z, charges = zip(*[
                        (d["name"], d["x"], d["y"], d["z"], (di, d["charge"]))
                        if "charge" in d
                        else (d["name"], d["x"], d["y"], d["z"], (di, 0))
                        for di, d in enumerate(residue_dict["beads"])
                    ])
                    residue = {
                        "resname": residue_dict["resname"],
                        "names": names,
                        "x": x,
                        "y": y,
                        "z": z,
                        "charges": charges,
                        **{key:val for key, val in residue_dict.items() if key not in ["resname", "beads"]}
                    }
                    residues_list.append(copy.deepcopy(residue))
        
        else:
            assert all(val in cur_dict for val in ["beads"]) or all(val in cur_dict for val in ["names", "x", "y", "z"]), "\n".join([
                "Residues for molecule are not specified but the molecule does not contain all required values that a normal residue would need.",
                "    "+"Name of current molecule: " + cur_name,
                "Residues can be written in two different ways:",
                "    "+"Version 1: A residue must contain 'beads'.",
                "    "+"    "+"'beads' is a list of dictionaries each containing information for a specific bead. The dictionaries must contain the following keys.",
                "    "+"    "+"    "+"'name' Bead name.",
                "    "+"    "+"    "+"'x' is the bead x-coordinate.",
                "    "+"    "+"    "+"'y' is the bead y-coordinate.",
                "    "+"    "+"    "+"'z' is the bead z-coordinate.",
                "    "+"    "+"    "+"'charge' is the charge of the bead (optional, assumed to be zero if not specified).",
                "    "+"Version 2: A residue must contain 'names', 'x', 'y' and 'z' keys.",
                "    "+"    "+"'names' is a list or tuple of all bead names in the residue.",
                "    "+"    "+"'x' is a list or tuple of the x-coordinates for all names in the residue.",
                "    "+"    "+"'y' is a list or tuple of the y-coordinates for all names in the residue.",
                "    "+"    "+"'z' is a list or tuple of the z-coordinates for all names in the residue.",
                "    "+"    "+"The i'th value in 'names', 'x', 'y' and 'z' is connected to each other.",
            ])

            residue_dict = cur_dict

            ### If separate lists/tuples data type
            if "names" in residue_dict:
                residue = {
                    "resname": cur_name,
                    **{key:val for key, val in residue_dict.items()}
                }
                residues_list.append(copy.deepcopy(residue))
            
            ### If dictionary data type
            elif "beads" in residue_dict:
                names, x, y, z, charges = zip(*[
                    (d["name"], d["x"], d["y"], d["z"], (di, d["charge"]))
                    if "charge" in d
                    else (d["name"], d["x"], d["y"], d["z"], (di, 0))
                    for di, d in enumerate(residue_dict["beads"])
                ])
                residue = {
                    "resname": cur_name,
                    "names": names,
                    "x": x,
                    "y": y,
                    "z": z,
                    "charges": charges,
                    **{key:val for key, val in residue_dict.items() if key not in ["beads"]}
                }
                residues_list.append(copy.deepcopy(residue))
        
        ### Ensures that "names" are always a tuple of strings rather than an individual string
        for res_nr, res_dict in enumerate(residues_list):
            if type(res_dict["names"]) in [str, int, float]:
                residues_list[res_nr]["names"] = (res_dict["names"],)

        ### Ensures that "x" are always a tuple of strings rather than an individual string
        for res_nr, res_dict in enumerate(residues_list):
            if type(res_dict["x"]) in [str, int, float]:
                residues_list[res_nr]["x"] = (res_dict["x"],)

        ### Ensures that "y" are always a tuple of strings rather than an individual string
        for res_nr, res_dict in enumerate(residues_list):
            if type(res_dict["y"]) in [str, int, float]:
                residues_list[res_nr]["y"] = (res_dict["y"],)

        ### Ensures that "z" are always a tuple of strings rather than an individual string
        for res_nr, res_dict in enumerate(residues_list):
            if type(res_dict["z"]) in [str, int, float]:
                residues_list[res_nr]["z"] = (res_dict["z"],)

        ### Molecule-wide charges
        if "charges" in cur_dict.keys():

            ### Ensures tuple of tuples
            if type(cur_dict["charges"][0]) in [str, int, float]:
                cur_dict["charges"] = (cur_dict["charges"],)

            for charge_values in cur_dict["charges"]:
                
                ### ### At this point "residues" has been added to the dictionary even if it wasn't originally
                ### Multiple residues (Residues specification is required)
                if len(residues_list) > 1:
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
            
        self.print_term("    ", "Residue dictionary", residues_list, debug=True)

        ### Residue-specific charges
        ### No need to check for residues as they are created earlier in this script
        for res_nr, res_dict in enumerate(residues_list):
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
        ### Ideally only use this for single bead molecules where the bead is obvious
        if "charge" in cur_dict.keys():
            number = self.get_number_from_string(cur_dict["charge"])
            if number is not False:
                tot_charge = number
        
        self.print_term("    ", "Charge dictionary", charge_dict, debug=True)
        self.print_term("    ", "tot_charge       ", tot_charge, debug=True)
        self.print_term("    ", "moleculetype     ", mol_class.moleculetype, debug=True)

        ### Rotates a lipid such that it is vertically aligned based on the designated upwards and downwards pointing beads
        ### ### First version aligns based on line between middle of upbeads and middle of downbeads
        if "alignment" in cur_dict.keys():
            if cur_dict["alignment"] == "manual":
                up_coords = []
                down_coords = []
                for beadname, coords_list in [("upbeads", up_coords), ("downbeads", down_coords)]:
                    for res_nr, bead_nr in cur_dict[beadname]:
                        x = round(residues_list[res_nr]["x"][bead_nr], 4)
                        y = round(residues_list[res_nr]["y"][bead_nr], 4)
                        z = round(residues_list[res_nr]["z"][bead_nr], 4)
                        coords_list.append((x, y, z))
                
                up_coords_array        = np.array(up_coords)
                down_coords_array      = np.array(down_coords)
                up_coords_mean         = np.mean(up_coords_array, axis=0)
                down_coords_mean       = np.mean(down_coords_array, axis=0)
                original_vector        = up_coords_mean - down_coords_mean
                original_vector_length = math.sqrt(original_vector[0]**2 + original_vector[1]**2+original_vector[2]**2)
                alignment_vector       = [0, 0, original_vector_length]
                rotation_matrix        = self.rotation_matrix_from_vectors(original_vector, alignment_vector)
                for res_nr, res_dict in enumerate(residues_list, 0):
                    xs, ys, zs = res_dict["x"], res_dict["y"], res_dict["z"]
                    new_xs, new_ys, new_zs = [], [], []
                    for x, y, z in zip(xs, ys, zs):
                        coord = np.array((x, y, z))
                        new_coord = np.dot(coord, rotation_matrix.T)
                        new_xs.append(new_coord[0])
                        new_ys.append(new_coord[1])
                        new_zs.append(new_coord[2])
                    res_dict["x"], res_dict["y"], res_dict["z"] = new_xs, new_ys, new_zs
            
            ### ### Second version aligns based on principal axis using a principal component analysis.
            ### ### Uses "upbeads" or "downbeads" to determine if final structure should be rotated 180 degrees around the x-axis to ensure "upbeads" are pointing up / "downbeads" are pointing down.
            elif cur_dict["alignment"] == "principal":

                ### Code in the following three functions was specifically written with a lot of help from ChatGPT. Beware of potential errors.
                def rodrigues_rotation_formula(axis, theta):
                    """
                    Rodrigues' rotation formula to create a rotation matrix.
                    """
                    axis = np.asarray(axis, dtype=float)
                    axis = axis / np.linalg.norm(axis)
                    a    = np.cos(theta)
                    b    = np.sin(theta)

                    K = np.array([
                        [       0, -axis[2],  axis[1]],
                        [ axis[2],        0, -axis[0]],
                        [-axis[1],  axis[0],        0],
                    ])

                    R = a * np.eye(3) + b * K + (1 - a) * np.outer(axis, axis)
                    return R

                def principal_axis_to_z(points):
                    """
                    Given an (N,3) array of 3D points, compute a rotation matrix that aligns the
                    principal axis (direction of maximum variance) with the z-axis [0, 0, 1].
                    Returns:
                        R : (3,3) rotation matrix
                        pts_rot : (N,3) rotated points (after centering and rotation)
                        principal : (3,) principal axis (unit vector, in original coordinates)
                        var_explained : fraction of total variance along the principal axis
                        centroid : (3,) centroid that was subtracted
                    """
                    pts = np.asarray(points, dtype=float)
                    if pts.ndim != 2 or pts.shape[1] != 3:
                        raise ValueError("points must be an (N,3) array")
                    
                    # Center the points
                    centroid = pts.mean(axis=0)
                    X        = pts - centroid

                    ### Covariance and eigen decomposition
                    cov              = np.cov(X, rowvar=False, bias=False)
                    eigvals, eigvecs = np.linalg.eigh(cov)
                    idx              = np.argsort(eigvals)[::-1]
                    eigvals          = eigvals[idx]
                    eigvecs          = eigvecs[:, idx]
                    principal        = eigvecs[:, 0] / np.linalg.norm(eigvecs[:, 0])
                    var_explained    = eigvals[0] / np.sum(eigvals)

                    ### Compute rotation to z-axis
                    z     = np.array([0.0, 0.0, 1.0])
                    dot   = np.clip(np.dot(principal, z), -1.0, 1.0)
                    angle = math.acos(dot)

                    if math.isclose(angle, 0.0, abs_tol=1e-12):
                        R = np.eye(3)
                    elif math.isclose(angle, np.pi, rel_tol=1e-12, abs_tol=1e-12):
                        ### 180 degree rotation: choose arbitrary perpendicular axis
                        v     = np.array([1.0, 0.0, 0.0]) if abs(principal[0]) < abs(principal[1]) else np.array([0.0, 1.0, 0.0])
                        axis  = np.cross(principal, v)
                        axis /= np.linalg.norm(axis)
                        R     = rodrigues_rotation_formula(axis, np.pi)
                    else:
                        axis  = np.cross(principal, z)
                        axis /= np.linalg.norm(axis)
                        R     = rodrigues_rotation_formula(axis, angle)

                    pts_rot = (R @ X.T).T

                    return R, pts_rot, principal, var_explained, centroid

                def rotate_x_180(points):
                    """
                    Rotates an (N,3) array of 3D points 180 degrees around the x-axis.
                    Returns rotated points.
                    """
                    ### 180° rotation matrix around x-axis
                    R = np.array([
                        [1,  0,  0],
                        [0, -1,  0],
                        [0,  0, -1]
                    ])
                    return (R @ points.T).T

                ### First get particle points into a usable format
                points = []
                for res_nr, res_dict in enumerate(residues_list, 0):
                    xs, ys, zs = res_dict["x"], res_dict["y"], res_dict["z"]
                    points.extend(list(zip(xs, ys, zs)))
                points = np.array(points)

                ### Run principal axis aligner
                R, points_rotated, principal, explained, centroid = principal_axis_to_z(points)

                self.print_term("Principal axis (unit vector, original coords):", principal, spaces=3, debug=True, debug_keys=["molecule_defs_checker"])
                self.print_term("Variance explained by principal axis:", explained,          spaces=3, debug=True, debug_keys=["molecule_defs_checker"])
                self.print_term("Centroid (subtracted):", centroid,                          spaces=3, debug=True, debug_keys=["molecule_defs_checker"])
                self.print_term("Rotation matrix R:\n", R,                                   spaces=3, debug=True, debug_keys=["molecule_defs_checker"])

                ### Check if both of up/down beads have been given
                assert not ("upbeads" in cur_dict.keys() and "downbeads" in cur_dict.keys()), "Both 'upbeads' and 'downbeads' have been specified for molecule '{cur_name}' while the 'principal' alignment method is being used. Please only use one of the two beads with 'principal' alignment.".format(cur_name=cur_name)
                
                ### Check if either of up/down beads has been given
                assert "upbeads" in cur_dict.keys() or "downbeads" in cur_dict.keys(), "Neither 'upbeads' nor 'downbeads' have been specified for molecule '{cur_name}' while the 'principal' alignment method is being used. Please specify exactly one of the two with 'principal' alignment.".format(cur_name=cur_name)
                
                ### Check which of up/down beads has given
                pointer_direction = False
                if "upbeads" in cur_dict.keys():
                    pointer_direction = "upbeads"
                elif "downbeads" in cur_dict.keys():
                    pointer_direction = "downbeads"

                ### Assign new coordinates to general dict and get up/down bead coords
                bead_i = 0
                pointer_coords = []
                for res_nr, res_dict in enumerate(residues_list, 0):
                    new_xs, new_ys, new_zs = [], [], []
                    bead_i_in_res = 0
                    for x, y, z in zip(xs, ys, zs):
                        new_xs.append(points_rotated[bead_i][0])
                        new_ys.append(points_rotated[bead_i][1])
                        new_zs.append(points_rotated[bead_i][2])
                        if (res_nr, bead_i_in_res) in cur_dict[pointer_direction]:
                            pointer_coords.append((points_rotated[bead_i][0], points_rotated[bead_i][1], points_rotated[bead_i][2]))
                        bead_i += 1
                        bead_i_in_res += 1
                    res_dict["x"], res_dict["y"], res_dict["z"] = new_xs, new_ys, new_zs

                assert len(pointer_coords) > 0, "The beads {beads} given using '{pointer_direction}' for molecule '{cur_name}' were not found.".format(beads=cur_dict[pointer_direction], pointer_direction=pointer_direction, cur_name=cur_name)

                pointer_coords_array = np.array(pointer_coords)
                pointer_coords_mean  = np.mean(pointer_coords_array, axis=0)

                points_rotated_mean = np.mean(points_rotated, axis=0)

                ### Check if they are on the right side of the center and rotate 180 around the x-axis if not
                new_rotated_points = False
                if pointer_direction == "upbeads" and pointer_coords_mean[2] < points_rotated_mean[2]:
                    new_rotated_points = rotate_x_180(points_rotated)
                elif pointer_direction == "downbeads" and pointer_coords_mean[2] > points_rotated_mean[2]:
                    new_rotated_points = rotate_x_180(points_rotated)

                ### Assign new coordinates to general dict if points have been rotated
                if new_rotated_points is not False:
                    bead_i = 0
                    for res_nr, res_dict in enumerate(residues_list, 0):
                        new_xs, new_ys, new_zs = [], [], []
                        for x, y, z in zip(xs, ys, zs):
                            new_xs.append(new_rotated_points[bead_i][0])
                            new_ys.append(new_rotated_points[bead_i][1])
                            new_zs.append(new_rotated_points[bead_i][2])
                            bead_i += 1
                        res_dict["x"], res_dict["y"], res_dict["z"] = new_xs, new_ys, new_zs

        ### Calculates 'minz' used in 'molecule_beads_checker' with lipids as it must be for all beads not only the beads in a specific residue
        minz = False
        if type_of_molecule == "lipid":
            minz = min([min(residue["z"]) for residue in residues_list])
            self.print_term("minz (for lipids only):", minz, spaces=3, debug=True, debug_keys=["molecule_defs_checker"])

        ### Dictionaries without "residues" key have had the beads and positions converted to "residues" key:val pair
        ### Thus "residues" will always be present
        molbeadnrs = []
        total_number_of_beads = sum([len(res_dict["names"]) for res_dict in residues_list])
        for res_nr, res_dict in enumerate(residues_list, 0):
            resname = res_dict["resname"]
            beads, xs, ys, zs = self.molecule_beads_checker(res_dict, cur_name, resname, bd, type_of_molecule, minz)

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

            tags = self.tags_checker(cur_dict)
            mol_class.add_tags(tags)

        
        self.print_term("", debug=True) ### Spacer for terminal and log file

        
        if type_of_molecule == "solvent":
            mol_class.set_coords_to_center(centering = "axis")
        elif type_of_molecule == "lipid":
            mol_class.set_coords_to_center(centering = "axis", AXs = "xy")
    