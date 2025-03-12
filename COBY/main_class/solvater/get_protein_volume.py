class get_protein_volume:
    def get_protein_volume(self, cxmin, cxmax, cymin, cymax, czmin, czmax):
        prot_beads_for_cell_checker = []
        solvent_box_proteins_charge = 0
        if len(self.PROTEINS) > 0:
            '''
            Finds the proteins beads and estimates their volume
            '''
            for protein_i, (protein_nr, protein) in enumerate(self.PROTEINS.items()):
                beads = protein["protein"].get_beads("xyz")
                charges = protein["protein"].get_bead_charges()
                for (x, y, z), charge in zip(beads, charges):
                    if cxmin < x <= cxmax and cymin < y <= cymax and czmin < z <= czmax:
                        prot_beads_for_cell_checker.append((x, y, z))
                        solvent_box_proteins_charge += charge
        return prot_beads_for_cell_checker, solvent_box_proteins_charge
