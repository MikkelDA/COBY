import numpy as np

from COBY.structure_classes.ATOM_class import ATOM

class RESIDUE:
    def __init__(self, resname, resnumber):
        self.resname = resname
        self.resnr = resnumber
        self.beads = []
        
    def add_bead_to_res(self, bead, beadnr, x, y, z, charge=0):
        self.beads.append(ATOM(bead, beadnr, x, y, z, self.resname, self.resnr, charge))
            
    def add_beads_to_res(self, beads = False, beadnrs = False, xs = False, ys = False, zs = False, charges = False):
        assert beads and beadnrs and xs and ys and zs, "Lacking data for either 'beads', 'beadnr', 'xs', 'ys' or 'zs'"
        if not charges:
            charges = [0 for _ in range(len(beads))]
        for bead, beadnr, x, y, z, charge in zip(beads, beadnrs, xs, ys, zs, charges):
            self.beads.append(ATOM(bead, beadnr, x, y, z, self.resname, self.resnr, charge))
    
    def add_bead_data_to_res(self, bead_data):
        if len(bead_data) == 5:
            ### Adding charge if not given
            bead_data.append([0 for _ in range(len(bead_data[0]))])
        for bead, beadnr, x, y, z, charge in zip(bead_data):
            self.beads.append(ATOM(bead, beadnr, x, y, z, self.resname, self.resnr, charge))
    
    def get_coords_res(self, AXs = "xyz"):
        AXs.lower()
        AXsList = []
        if AXs == "all":
            AXs = "xyz"
        for AX in AXs:
            if AX in ["x"]:
                AXsList.append([bead.x for bead in self.beads])
            if AX in ["y"]:
                AXsList.append([bead.y for bead in self.beads])
            if AX in ["z"]:
                AXsList.append([bead.z for bead in self.beads])
        if len(AXsList) == 1:
            AXsList = AXsList[0]
        return AXsList
    
    def get_center_point(self, centering = "mean_of_extremes", AXs = "xyz"):
        coords_AXs_res = self.get_coords_res(AXs)
        if len(AXs) == 1:
            coords_AXs_res = [coords_AXs_res]
        
        centers = []
        for ci, coords in enumerate(coords_AXs_res):
            if centering  in ["axis", "mean_of_extremes"]:
                centers.append((max(coords) + min(coords)) / 2)
                
            elif centering in ["cog", "mean_of_beads"]:
                centers.append(np.mean(coords))
                
            elif centering.startswith("beadnr"):
                bead_nrs = []
                for bead_nr in target:
                    if "-" in bead_nr:
                        bead_nr1, bead_nr2 = bead_nr.split("-")
                        bead_nrs.extend(list(range(int(bead_nr1), int(bead_nr2)+1)))
                    else:
                        bead_nrs.append(int(bead_nr))
                bead_ax_vals = [coords[bead_nr] for bead_nr in bead_nrs]

                if any([centering.endswith(value) for value in ["resnr", "cog", "mean_of_beads"]]):
                    centers.append(np.mean(bead_ax_vals))
                elif any([centering.endswith(value) for value in ["axis", "mean_of_extremes"]]):
                    centers.append((max(bead_ax_vals) + min(bead_ax_vals)) / 2)
                
            elif centering == "vals":
                centers.append(target[ci])
        return tuple(centers)

