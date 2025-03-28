from operator import itemgetter

class get_solute_volume:
    def get_solute_volume(self, cxmin, cxmax, cymin, cymax, czmin, czmax):
        solvent_beads_for_cell_checker = []
        solvent_box_solute_charge = 0
        if len(self.SOLVATIONS) > 0:
            '''
            Finds the solvents beads and estimates their volume
            '''
            for solv_key, solv_dict in self.SOLVATIONS.items():
                if "grid" in solv_dict:
                    for grid_point in solv_dict["grid"]:
                        beads, charges = itemgetter("coords", "bead_charges")(grid_point)
                        for (x, y, z), charge in zip(beads, charges):
                            if cxmin < x <= cxmax and cymin < y <= cymax and czmin < z <= czmax:
                                solvent_beads_for_cell_checker.append((x, y, z))
                                solvent_box_solute_charge += charge
        return solvent_beads_for_cell_checker, solvent_box_solute_charge

