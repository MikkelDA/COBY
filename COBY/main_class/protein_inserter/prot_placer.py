import time

class prot_placer:
    def prot_placer(self):
        '''
        Places all proteins into the systems internal coordinate system
        Checks if all atoms/beads are within the pbc and moves them if they are outside
        '''
        if len(self.PROTEINS) != 0:
            prot_placer_tic = time.time()
            string = " ".join(["", "PROTEIN PLACEMENT", ""])
            self.print_term("{string:-^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), spaces=0, verbose=1)
            for protein_i, (protein_nr, protein) in enumerate(self.PROTEINS.items()):
                if protein_i != 0:
                    self.print_term("", verbose=2)
                self.print_term("Starting protein nr", protein_nr, spaces=0, verbose=2)
                
                ### Printing charge information
                self.print_term(
                    "Total charge of protein is: {}".format(protein["protein"].get_mol_charge()),
                    spaces=1,
                    verbose=2,
                )

                #################
                ### CENTERING ###
                #################

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
                
                xcen, ycen, zcen = self.PROTEINS[protein_nr]["protein"].get_center_point(centering = centering, target = target)
                self.print_term(
                    "Centering protein using", "'" + " ".join([str(i) for i in protein["cen_method"]])+"'",
                    "at x/y/z:", round(xcen, 3), round(ycen, 3), round(zcen, 3), "(Input file coordinate system [Å])",
                    spaces=1,
                    verbose=2
                )
                self.PROTEINS[protein_nr]["protein"].set_coords_to_center(centering = centering, target = target)

                #################
                ### ROTATIONS ###
                #################
                x_deg, y_deg, z_deg = protein["rx"], protein["ry"], protein["rz"]
                
                if any([ang != 0 for ang in [x_deg, y_deg, z_deg]]):
                    self.PROTEINS[protein_nr]["protein"].rotate_coords(rotation = [x_deg, y_deg, z_deg])

                ####################
                ### TRANSLATIONS ###
                ####################
                cx, cy, cz = protein["cx"], protein["cy"], protein["cz"]
                if any([ax != 0 for ax in [cx, cy, cz]]):
                    self.PROTEINS[protein_nr]["protein"].move_coords(translation = [cx, cy, cz])

                #################################################
                ### CHECKS IF COORDINATES ARE OUTSIDE THE BOX ###
                #################################################
                if protein["pbc_check"]:
                    beads = self.PROTEINS[protein_nr]["protein"].get_res_beads_info()
                    
                    errors_count = 0
                    for bead in beads:
                        bead_coords = [bead.x, bead.y, bead.z]
                        checked_beads, error = self.coord_checker(bead_coords, self.pbc_box, error_count = True)
                        if bead_coords != checked_beads:
                            bi = bead.bead
                            ri = bead.resnr
                            self.PROTEINS[protein_nr]["protein"].residues[ri].beads[bi].move_atom(checked_beads)
                            
                        if error > 0:
                            errors_count += 1
                    if errors_count > 0:
                        self.print_term("WARNING:", str(errors_count), "Beads are outside pbc. Moved to other side. Expect potential problems from this.", warn = True, spaces=2)
                        self.print_term("Please move the protein such that it fits within the pbc.", warn = True, spaces=2)

                xcen_new, ycen_new, zcen_new = 0, 0, 0
                xcen_new, ycen_new, zcen_new = xcen_new + cx, ycen_new + cy, zcen_new + cz
                self.print_term("New protein center at x/y/z:", round(xcen_new, 3), round(ycen_new, 3), round(zcen_new, 3), "(Internal coordinate system [Å])", spaces=1, verbose=2)
                
                self.print_term("Finished placing protein nr", protein_nr, spaces=1, verbose=2)
                
            prot_placer_toc = time.time()
            prot_placer_time = round(prot_placer_toc - prot_placer_tic, 4)
            string = " ".join(["", "PROTEIN PLACEMENT COMPLETE", ""])
            self.print_term("{string:-^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), spaces=0, verbose=1)
            string = " ".join(["", "(Time spent:", str(prot_placer_time), "[s])", ""])
            self.print_term("{string:^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), "\n", spaces=0, verbose=1)

