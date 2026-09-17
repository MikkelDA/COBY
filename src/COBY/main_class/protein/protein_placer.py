import time
import numpy as np
import math

class protein_placer:
    def protein_placer(self):
        '''
        Places all proteins into the systems internal coordinate system
        Checks if all atoms/beads are within the pbc and moves them if they are outside
        '''
        if len(self.PROTEINS) != 0:
            protein_placer_tic = time.time()
            string = " ".join(["", "PROTEIN PLACEMENT", ""])
            self.print_term("{string:-^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), spaces=0, verbose=1)
            for protein_i, (protein_nr, protein) in enumerate(self.PROTEINS.items()):
                if protein_i != 0:
                    self.print_term("", verbose=2)
                self.print_term("Starting protein nr", protein_nr, spaces=0, verbose=2)
                
                ### Printing charge information
                self.print_term(
                    "Total charge of protein is: {}".format(round(protein["protein"].get_mol_charge(), 6)),
                    spaces=1,
                    verbose=2,
                )

                #################
                ### ALIGNMENT ###
                #################
                def rotation_matrix_to_euler_xyz(R):
                    """
                    Extract Euler angles (degrees) from rotation matrix
                    using the SAME convention as your rotate_coords:
                    R = Rx * Ry * Rz
                    """

                    # Clamp to avoid numerical issues
                    r02 = max(min(R[0][2], 1.0), -1.0)

                    y = math.asin(r02)
                    cy = math.cos(y)

                    # Check for gimbal lock
                    if abs(cy) > 1e-6:
                        x = math.atan2(-R[1][2], R[2][2])
                        z = math.atan2(-R[0][1], R[0][0])
                    else:
                        # Gimbal lock case
                        x = math.atan2(R[2][1], R[1][1])
                        z = 0.0

                    return (
                        math.degrees(x),
                        math.degrees(y),
                        math.degrees(z)
                    )

                ### Rotates a protein such that it is vertically aligned based on the designated upwards and downwards residues
                if protein["alignment"] == "manual":
                    residues_list = self.PROTEINS[protein_nr]["protein"].get_res_beads_info(output_type="tuple")
                    
                    ### list of (beadname, beadnr, x, y, z, resname, resnr, charge) tuples
                    up_coords = []
                    down_coords = []
                    for beadname, beadnr, x, y, z, resname, resnr, charge in residues_list:
                        if resnr in protein["upres"]:
                            x = round(x, 4)
                            y = round(y, 4)
                            z = round(z, 4)
                            up_coords.append((x, y, z))
                        if resnr in protein["downres"]:
                            x = round(x, 4)
                            y = round(y, 4)
                            z = round(z, 4)
                            down_coords.append((x, y, z))

                    assert len(up_coords) > 0 and len(down_coords) > 0, "Zero particles found matching the residues for either 'upres' ({up_coords_len} particles) or 'downres' ({down_coords_len} particles) for manual protein alignment.".format(up_coords_len=len(up_coords), down_coords_len=len(down_coords))
                
                    up_coords_array        = np.array(up_coords)
                    down_coords_array      = np.array(down_coords)
                    up_coords_mean         = np.mean(up_coords_array, axis=0)
                    down_coords_mean       = np.mean(down_coords_array, axis=0)
                    original_vector        = up_coords_mean - down_coords_mean
                    original_vector_length = math.sqrt(original_vector[0]**2 + original_vector[1]**2+original_vector[2]**2)
                    alignment_vector       = [0, 0, original_vector_length]
                    rotation_matrix        = self.rotation_matrix_from_vectors(original_vector, alignment_vector)
                    x_deg, y_deg, z_deg    = rotation_matrix_to_euler_xyz(rotation_matrix)
                    self.PROTEINS[protein_nr]["protein"].rotate_coords(rotation = [x_deg, y_deg, z_deg])

                #################
                ### CENTERING ###
                #################

                ### Centers the protein using the given centering algorithm
                if protein["center_protein"]:
                    ### Centered on the mean of largest/smallest x/y/z coordinate
                    if protein["cen_method"][0] in ["cog", "mean_of_beads"]: # Default
                        centering = "mean_of_beads"
                        target = False

                    ### Centered on the mean coordinate of all beads (center of geometry)
                    elif protein["cen_method"][0] in ["axis", "mean_of_extremes"]:
                        centering = "mean_of_extremes"
                        target = False

                    ### Centered on a single bead
                    elif protein["cen_method"][0].startswith("bead"):
                        centering = "beadnr"
                        if len(protein["cen_method"][0]) > 4:
                            centering = centering + protein["cen_method"][0][4:]
                        target = protein["cen_method"][1] # a list of bead numbers

                    ### Centered on the mean position of all beads in a single residue
                    elif protein["cen_method"][0].startswith("res"):
                        centering = "resnr"
                        if len(protein["cen_method"][0]) > 3:
                            centering = centering + protein["cen_method"][0][3:]
                        target = protein["cen_method"][1] # a list of residue numbers

                    ### Centered on the specific x/y/z coordinates
                    elif protein["cen_method"][0] == "point":
                        centering = "vals"
                        target = protein["cen_method"][1:]
                    
                    ### The "get_center_point" method is only run here to get values for printing.
                    ### It is called inside the "set_coords_to_center" method separately.
                    xcen, ycen, zcen = self.PROTEINS[protein_nr]["protein"].get_center_point(centering = centering, target = target)
                    self.print_term(
                        "Centering protein using", "'" + " ".join([str(i) for i in protein["cen_method"]])+"'",
                        "at x/y/z:", round(xcen, 3), round(ycen, 3), round(zcen, 3), "(Input file coordinate system [Å])",
                        spaces=1,
                        verbose=2
                    )
                    self.PROTEINS[protein_nr]["protein"].set_coords_to_center(centering = centering, target = target)
                else:
                    ### Adjusts the coordinates to account for centrosymmetric box used in COBY
                    xlen, ylen, zlen = self.PROTEINS[protein_nr]["box_info"]["x"], self.PROTEINS[protein_nr]["box_info"]["y"], self.PROTEINS[protein_nr]["box_info"]["z"]
                    self.PROTEINS[protein_nr]["protein"].move_coords(translation = [-xlen/2, -ylen/2, -zlen/2])

                #################
                ### ROTATIONS ###
                #################
                for rotation in protein["rotate"]:
                    x_deg, y_deg, z_deg = rotation["rx"], rotation["ry"], rotation["rz"]
                    self.PROTEINS[protein_nr]["protein"].rotate_coords(rotation = [x_deg, y_deg, z_deg])

                ####################
                ### TRANSLATIONS ###
                ####################
                cx, cy, cz = protein["cx"], protein["cy"], protein["cz"]
                if any([ax != 0 for ax in [cx, cy, cz]]):
                    self.PROTEINS[protein_nr]["protein"].move_coords(translation = [cx, cy, cz])

                ###################################################################
                ### CHECKS IF COORDINATES ARE OUTSIDE THE BOX AND MOVES THEM IN ###
                ###################################################################
                if protein["pbc_check"]:
                    errors_count = 0
                    for ri, residue in enumerate(self.PROTEINS[protein_nr]["protein"].residues):
                        for bi, bead in enumerate(residue.beads):
                            bead_coords = [bead.x, bead.y, bead.z]
                            checked_beads, error = self.coord_checker(bead_coords, self.pbc_box, error_count = True)

                            if error > 0:
                                errors_count += 1
                                self.PROTEINS[protein_nr]["protein"].residues[ri].beads[bi].move_atom(*checked_beads)

                    if errors_count > 0:
                        self.print_term(str(errors_count), "Protein beads are outside pbc. Moved to other side. Expect potential problems from this. Please move the protein such that it fits within the pbc.", warn = True, spaces=2)

                xcen_new, ycen_new, zcen_new = 0, 0, 0
                xcen_new, ycen_new, zcen_new = xcen_new + cx, ycen_new + cy, zcen_new + cz
                self.print_term("New protein center at x/y/z:", round(xcen_new, 3), round(ycen_new, 3), round(zcen_new, 3), "(Internal coordinate system [Å])", spaces=1, verbose=2)
                
                self.protein_beads_in_sys += len(self.PROTEINS[protein_nr]["beads"])

                self.print_term("Finished placing protein nr", protein_nr, spaces=1, verbose=2)
                
            protein_placer_toc = time.time()
            protein_placer_time = round(protein_placer_toc - protein_placer_tic, 4)
            string = " ".join(["", "PROTEIN PLACEMENT COMPLETE", ""])
            self.print_term("{string:-^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), spaces=0, verbose=1)
            string = " ".join(["", "(Time spent:", str(protein_placer_time), "[s])", ""])
            self.print_term("{string:^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), "\n", spaces=0, verbose=1)

