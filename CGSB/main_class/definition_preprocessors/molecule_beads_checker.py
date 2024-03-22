class molecule_beads_checker:
    def molecule_beads_checker(self, res_dict, cur_name, resname, bd, type_of_molecule):
        xs = res_dict["x"]
        ys = res_dict["y"]
        zs = res_dict["z"]
        beads = res_dict["beads"]
        if type(xs) != tuple and type(xs) != list:
            xs = (xs,)
        if type(ys) != tuple and type(ys) != list:
            ys = (ys,)
        if type(zs) != tuple and type(zs) != list:
            zs = (zs,)
        if type(beads) != tuple and type(beads) != list:
            beads = (beads,)
        
        assert len(beads) == len(xs) == len(ys) == len(zs), (
            "Number of beads, x-values, y-values and z-values must be the same within a residue." + "\n"
            "molecule name: " + cur_name + "\n"
            "name of problematic residue: " + resname + "\n"
            "len(beads) " + str(len(beads)) + ", values: " + str(beads) + "\n"
            "len(xs)    " + str(len(xs)) + ", values: " + str(xs) + "\n"
            "len(ys)    " + str(len(ys)) + ", values: " + str(ys) + "\n"
            "len(zs)    " + str(len(zs)) + ", values: " + str(zs) + "\n"
        )
        
        bdx, bdy, bdz = bd

        lx = [xi * bdx * 10 for xi in xs]
        ly = [yi * bdy * 10 for yi in ys]
        if type_of_molecule == "solvent":
            lz = [zi * bdz * 10 for zi in zs]
        elif type_of_molecule == "lipid":
            minz = min(zs)
            inter_leaflet_buffer = 1.5
            lz = [(zi - minz) * bdz * 10 + inter_leaflet_buffer for zi in zs]
        
        return beads, lx, ly, lz
