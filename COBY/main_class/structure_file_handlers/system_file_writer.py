class system_file_writer:
    def system_file_writer(self):
        out_stringlist = []
        files_string = ""
        if self.output_system_pdb_file_name:
            out_stringlist.append("PDB")
        if self.output_system_gro_file_name:
            out_stringlist.append("GRO")
        if len(out_stringlist) > 1:
            files_string = "files"
        else:
            files_string = "file"

        string = " ".join(["", "Writing structure", files_string, "(" + "/".join(out_stringlist) + ")", ""])
        self.print_term("{string:-^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), spaces=0, verbose=1)
            
        output_system_pdb_file_lines = []
        output_system_gro_file_lines = []
        
        ###############################
        ### Beginning of file lines ###
        ###############################
        if self.output_system_pdb_file_name:
            a, b, c, alpha, beta, gamma = self.pdb_box_dimension
            output_system_pdb_file_lines = [
                "TITLE     " + self.system_name,
                "REMARK    " + "PLACEHOLDER_REMARK",
                '{Rname:<{RnameL}}{a:>{aL}.3f}{b:>{bL}.3f}{c:>{cL}.3f}{alpha:>{alphaL}.2f}{beta:>{betaL}.2f}{gamma:>{gammaL}.2f} {sGroup:<{sGroupL}}{z:>{zL}}'.format(
                    Rname = "CRYST1", RnameL = 6,
                    a = a,            aL = 9,
                    b = b,            bL = 9,
                    c = c,            cL = 9,
                    alpha = alpha,    alphaL = 7,
                    beta = beta,      betaL = 7,
                    gamma = gamma,    gammaL = 7,
                    sGroup = "P 1",   sGroupL = 11,
                    z = 1,            zL = 4,
                ),
                "MODEL        1",
            ]
        if self.output_system_gro_file_name:
            
            output_system_gro_file_lines = [
                self.system_name,
                "PLACEHOLDER_ATOM_COUNT",
            ]

        atom_count = 0
        self.molecules_for_top = []

        atom_nr    = 0
        res_nr     = 0
        
        #####################
        ### Protein lines ###
        #####################
        if len(self.PROTEINS) != 0:
            for protein_nr, protein in self.PROTEINS.items():
                for mol_name in protein["moleculetypes"]:
                    self.molecules_for_top.append((mol_name, str(1)))
                current_prot_res = 0
                
                original_bead_info = protein["beads"].keys()
                bead_vals = [bead.get_tuple() for bead in protein["protein"].get_res_beads_info()]
                for (i, atom, res), (a_name, beadnr, x, y, z, r_name, resnr, charge) in zip(original_bead_info, bead_vals):
#                 for (i, atom, res), bead_vals in protein["beads_centered"].items():
                    if current_prot_res != res:
                        current_prot_res = res
                        res_nr += 1
                        if res_nr >= 10000:
                            res_nr -= 10000 * (res_nr // 10000)
                    atom_nr += 1
                    atom_count += 1
                    if atom_nr >= 100000:
                        atom_nr -= 100000 * (atom_nr // 100000)
                    x += self.pbc_box[0] / 2
                    y += self.pbc_box[1] / 2
                    z += self.pbc_box[2] / 2
                    if self.output_system_pdb_file_name:
                        output_system_pdb_file_lines.append(self.pdb_atom_writer("ATOM", atom_nr, a_name, " ", r_name, "A", res_nr, " ", float(x), float(y), float(z), float(1), float(0), " ", " ", " "))
                    if self.output_system_gro_file_name: ### gro coordinates are in [nm] not [Å]
                        output_system_gro_file_lines.append(self.gro_atom_writer(res_nr, r_name, a_name, atom_nr, x / 10, y / 10, z / 10, " ", " ", " "))

        ######################
        ### Membrane lines ###
        ######################
        if len(self.MEMBRANES) != 0:
            for memb_key, memb_dict in self.MEMBRANES.items():
                for leaflet_key, leaflet in memb_dict["leaflets"].items():
                    self.molecules_for_top.extend(sorted(leaflet["leaf_lipid_count"], key=lambda m: m[0]))
                    
                    ### Ordering lipids per leaflet for concise topology creation
                    ordered_lipids = sorted(
                        [grid_vals for grid_vals in leaflet["grid_lipids"]],
                        key=lambda gp: gp["lipid"]["name"]
                    )
            
                    for i, grid_vals in enumerate(ordered_lipids):
                        res_nr += 1
                        if res_nr >= 10000:
                            res_nr -= 10000 * (res_nr // 10000)
                        r_name = grid_vals["lipid"]["name"]

                        for x, y, z, bead_name in zip(grid_vals["lipid"]["x"], grid_vals["lipid"]["y"], grid_vals["lipid"]["z"], grid_vals["lipid"]["beads"]):
                            atom_nr += 1
                            atom_count += 1
                            if atom_nr >= 100000:
                                atom_nr -= 100000 * (atom_nr // 100000)
                            x += self.pbc_box[0] / 2
                            y += self.pbc_box[1] / 2
                            z += self.pbc_box[2] / 2
                            a_name = bead_name
                            string = ""
                            if self.output_system_pdb_file_name:
                                output_system_pdb_file_lines.append(self.pdb_atom_writer("ATOM", atom_nr, a_name, " ", r_name, "A", res_nr, " ", float(x), float(y), float(z), float(1), float(0), " ", " ", " "))
                            if self.output_system_gro_file_name: ### gro coordinates are in [nm] not [Å]
                                output_system_gro_file_lines.append(self.gro_atom_writer(res_nr, r_name, a_name, atom_nr, x / 10, y / 10, z / 10, " ", " ", " "))
        
        #########################
        ### Solvent/ion lines ###
        #########################
        if len(self.SOLVATIONS) != 0:
            for solvation_nr, solvation in self.SOLVATIONS.items():
                solvent_count = [(key_name, count) for (key_name, key_type, key_charge), count in solvation["solv_count"].items()]
                self.molecules_for_top.extend(solvent_count)

                for i, grid_point_3D in enumerate(solvation["grid"]):
                    res_nr += 1
                    if res_nr >= 10000:
                        res_nr -= 10000 * (res_nr // 10000)

                    for (x, y, z), r_name, bead_name in zip(grid_point_3D["coords"], grid_point_3D["resnames"], grid_point_3D["beads"]):
                        atom_nr += 1
                        atom_count += 1
                        if atom_nr >= 100000:
                            atom_nr -= 100000 * (atom_nr // 100000)
                        x += self.pbc_box[0] / 2
                        y += self.pbc_box[1] / 2
                        z += self.pbc_box[2] / 2
                        a_name = bead_name
                        string = ""
                        if self.output_system_pdb_file_name:
                            output_system_pdb_file_lines.append(self.pdb_atom_writer("ATOM", atom_nr, a_name, " ", r_name, "A", res_nr, " ", float(x), float(y), float(z), float(1), float(0), " ", " ", " "))
                        if self.output_system_gro_file_name: ### gro coordinates are in [nm] not [Å]
                            output_system_gro_file_lines.append(self.gro_atom_writer(res_nr, r_name, a_name, atom_nr, x / 10, y / 10, z / 10, " ", " ", " "))
        
        #########################
        ### End of file lines ###
        #########################
        if self.output_system_pdb_file_name:
            output_system_pdb_file_lines.append("TER")
            output_system_pdb_file_lines.append("END")
        
        if self.output_system_gro_file_name:
            v1x, v2y, v3z, v1y, v1z, v2x, v2z, v3x, v3y = self.gro_box_vectors
            output_system_gro_file_lines[1] = " " + str(atom_count)
            output_system_gro_file_lines.append( ### gro vectors are in [nm] not [Å]
                '{v1x:>{v1xL}.5f}{v2y:>{v2yL}.5f}{v3z:>{v3zL}.5f}{v1y:>{v1yL}.5f}{v1z:>{v1zL}.5f}{v2x:>{v2xL}.5f}{v2z:>{v2zL}.5f}{v3x:>{v3xL}.5f}{v3y:>{v3yL}.5f}'.format(
                    v1x = v1x, v1xL = 10,
                    v2y = v2y, v2yL = 10,
                    v3z = v3z, v3zL = 10,
                    v1y = v1y, v1yL = 10,
                    v1z = v1z, v1zL = 10,
                    v2x = v2x, v2xL = 10,
                    v2z = v2z, v2zL = 10,
                    v3x = v3x, v3xL = 10,
                    v3y = v3y, v3yL = 10,
                )
            )
            # output_system_gro_file_lines.append("") # Apparantly should not be present

        if self.output_system_pdb_file_name:
            if self.backup:
                if self.output_system_pdb_file_name:
                    self.backupper(self.output_system_pdb_file_name)
            new_file = open(self.output_system_pdb_file_name, "w")
            for line in output_system_pdb_file_lines:
                new_file.write(line + "\n")
            new_file.close()
            self.print_term("PDB file written:", self.output_system_pdb_file_name, verbose=1)
        
        if self.output_system_gro_file_name:
            if self.backup:
                if self.output_system_gro_file_name:
                    self.backupper(self.output_system_gro_file_name)
            new_file = open(self.output_system_gro_file_name, "w")
            for line in output_system_gro_file_lines:
                new_file.write(line + "\n")
            new_file.close()
            self.print_term("GRO file written:", self.output_system_gro_file_name, verbose=1)
        self.print_term("", verbose=1)

