from operator import itemgetter

class get_lipid_volume:
    def get_lipid_volume(self, cxmin, cxmax, cymin, cymax, czmin, czmax):
        lipid_beads_for_cell_checker = []
        solvent_box_lipids_charge = 0
        if len(self.MEMBRANES) > 0:
            '''
            Finds the lipids beads and estimates their volume
            '''
            for memb_key, memb_dict in self.MEMBRANES.items():
                for leaflet_key, leaflet in memb_dict["leaflets"].items():
                    for grid_point in leaflet["grid_lipids"]:
                        xs, ys, zs, charges = itemgetter("x", "y", "z", "charges")(grid_point["lipid"])
                        for x, y, z, charge in zip(xs, ys, zs, charges):
                            if cxmin < x <= cxmax and cymin < y <= cymax and czmin < z <= czmax:
                                lipid_beads_for_cell_checker.append((x, y, z))
                                solvent_box_lipids_charge += charge
        return lipid_beads_for_cell_checker, solvent_box_lipids_charge

