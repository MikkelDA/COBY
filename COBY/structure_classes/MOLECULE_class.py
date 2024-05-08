import math
import numpy as np

from COBY.structure_classes.RESIDUE_class import RESIDUE
from COBY.structure_classes.ATOM_class import ATOM
from COBY.general_functions.flatten import flatten

class MOLECULE:
    def __init__(self, molname = False, moleculetype = False):
        self.residues = []
        self.center = False
        self.n_residues = 0
        self.resnames = []
        self.molname = molname
        self.moleculetype = moleculetype
        self.last_res_n = 0
    
    def add_res(self, resname, resnumber = False):
        if not resnumber:
            resnumber = self.n_residues
#         self.residues[resnumber] = RESIDUE(resname)
        self.residues.append(RESIDUE(resname, resnumber))
        self.resnames.append(resname)
        self.last_res_n = self.n_residues
        self.n_residues += 1
    
    def add_bead(self, bead, beadnr, x, y, z, resnumber = False, charge=0):
        if not resnumber:
            resnumber = self.last_res_n
        self.residues[resnumber].add_bead_to_res(bead, beadnr, x, y, z, charge)
    
    def add_res_and_beads(self, resname, beads = False, beadnrs = False, xs = False, ys = False, zs = False, bead_data = False, resnumber = False, charges = False):
        if not resnumber:
            resnumber = self.n_residues
        self.residues.append(RESIDUE(resname, resnumber))
        self.resnames.append(resname)
        self.last_res_n = self.n_residues
        self.n_residues += 1
        if bead_data:
            self.residues[self.last_res_n].add_bead_data_to_res(bead_data)
        else:
            self.residues[self.last_res_n].add_beads_to_res(beads, beadnrs, xs, ys, zs, charges)
    
    def set_bead_charges(self, charges):
        i = 0
        for ri, res in enumerate(self.residues):
            for bi, bead in enumerate(res.beads):
                self.residues[ri].beads[bi].set_charge(charges[i])
                i += 1
    
    def get_coords(self, AXs = "xyz"):
        AXs.lower()
        AXsList = []
        if AXs == "all":
            AXs = "xyz"
        for AX in AXs:
            AXsList.append(flatten([res.get_coords_res(AX) for res in self.residues]))
#         if len(AXsList) == 1:
#             AXsList = AXsList[0]
        return AXsList
    
    def get_beads(self, AXs = "all"):
        AXsList = self.get_coords(AXs)
#         print(AXsList)
#         print(list(zip(*AXsList)))
        ### Not sure why following assert was made
#         assert len(AXsList) > 1, "Length of AXsList must be greater than 1 to create beads"
        return list(zip(*AXsList))
    
    def get_bead_names(self):
        bead_names = []
        for res in self.residues:
            res_bead_names = []
            for bead in res.beads:
                res_bead_names.append(bead.bead)
            bead_names.append(res_bead_names)
        ### If only 1 residue
        if len(bead_names) == 1:
            bead_names = bead_names[0]
        return bead_names

    def get_mol_charge(self):
        charge = sum(self.get_bead_charges())
        return charge

    def get_bead_charges(self):
        charges = [bead.charge for res in self.residues for bead in res.beads]
        return charges

    def get_centered_coords(self, centering = "mean_of_extremes", AXs = "xyz", target = False):
        coords_AXs = self.get_coords(AXs)
        centered_coords = []
        
        centers = self.get_center_point(centering = centering, AXs = AXs, target = target)

        for coords, center in zip(coords_AXs, centers):
            centered_coords.append([coord - center for coord in coords])
        
        return tuple(centered_coords)
    
    def get_center_point(self, centering = "mean_of_extremes", AXs = "xyz", target = False):
        coords_AXs = self.get_coords(AXs)
        
        centers = []
        for ci, coords in enumerate(coords_AXs):
            if centering in ["axis", "mean_of_extremes"]:
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

                if any([centering.endswith(value) for value in ["beadnr", "cog", "mean_of_beads"]]):
                    centers.append(np.mean(bead_ax_vals))
                elif any([centering.endswith(value) for value in ["axis", "mean_of_extremes"]]):
                    centers.append((max(bead_ax_vals) + min(bead_ax_vals)) / 2)
                
            elif centering.startswith("resnr"):
                res_nrs = []
                ### Find all residues
                for res_nr in target:
                    if "-" in res_nr:
                        res_nr1, res_nr2 = res_nr.split("-")
                        res_nrs.extend(list(range(int(res_nr1), int(res_nr2)+1)))
                    else:
                        res_nrs.append(int(res_nr))
                
                ### Find center of each individual residue
                res_centers = [
                    self.residues[res_nr].get_center_point(centering = "mean_of_beads", AXs = AXs[ci])
                    for res_nr in res_nrs
                ]

                ### Find center of residues
                if any([centering.endswith(value) for value in ["resnr", "cog", "mean_of_beads"]]):
                    centers.append(np.mean(res_centers))
                elif any([centering.endswith(value) for value in ["axis", "mean_of_extremes"]]):
                    centers.append((max(res_centers) + min(res_centers)) / 2)

            elif centering == "vals":
                centers.append(target[ci])
        
        return tuple(centers)
    
    def get_radius(self, AXs="xyz"):
        beads  = self.get_beads(AXs)
        center = self.get_center_point("mean_of_extremes", AXs)
        radius = max([math.dist(bead, center) for bead in beads])
        return radius
    
    def set_coords_to_center(self, centering = "mean_of_extremes", AXs = "xyz", target = False):
        if AXs == "xyz":
            xsc, ysc, zsc = self.get_centered_coords(centering, AXs, target)
        else:
            ### [0] because they are returned as tuples of lists
            if "x" in AXs:
                xsc = self.get_centered_coords(centering, "x", target)[0]
            if "y" in AXs:
                ysc = self.get_centered_coords(centering, "y", target)[0]
            if "z" in AXs:
                zsc = self.get_centered_coords(centering, "z", target)[0]

        i = 0
        for ri, res in enumerate(self.residues):
            for bi, bead in enumerate(res.beads):
                if "x" in AXs:
                    self.residues[ri].beads[bi].movex(xsc[i])
                if "y" in AXs:
                    self.residues[ri].beads[bi].movey(ysc[i])
                if "z" in AXs:
                    self.residues[ri].beads[bi].movez(zsc[i])
                i += 1
    
    def move_coords(self, translation = [0, 0, 0]):
        tx, ty, tz = translation
        i = 0
        for ri, res in enumerate(self.residues):
            for bi, bead in enumerate(res.beads):
                self.residues[ri].beads[bi].move_atom(bead.x+tx, bead.y+ty, bead.z+tz)
                i += 1
    
    def move_atom(self, ri, bi, translation = [0, 0, 0]):
        tx, ty, tz = translation
        bead = self.residues[ri].beads[i]
        self.residues[ri].beads[bi].move_atom(bead.x+tx, bead.y+ty, bead.z+tz)
    
    def rotate_coords(self, rotation = [0, 0, 0]):
        x_deg, y_deg, z_deg = rotation
        i = 0
        for ri, res in enumerate(self.residues):
            for bi, bead in enumerate(res.beads):
                
                ### Get positions to limit calls to bead
                x = bead.x
                y = bead.y
                z = bead.z
                
                ### First convert angles to radians
                x_rad = math.radians(x_deg)
                y_rad = math.radians(y_deg)
                z_rad = math.radians(z_deg)

                ### Calculate rotation matrices (https://en.wikipedia.org/wiki/Rotation_matrix)
                rm_x = [
                    [1, 0,                0              ],
                    [0, math.cos(x_rad), -math.sin(x_rad)],
                    [0, math.sin(x_rad),  math.cos(x_rad)],
                ]
                rm_y = [
                    [math.cos(y_rad),  0, math.sin(y_rad)],
                    [0,                1, 0              ],
                    [-math.sin(y_rad), 0, math.cos(y_rad)],
                ]
                rm_z = [
                    [math.cos(z_rad), -math.sin(z_rad), 0],
                    [math.sin(z_rad),  math.cos(z_rad), 0],
                    [0,                0,               1],
                ]

                ### Combine rotation matrices
                rm_xy  = [
                    [sum(r * c for r, c in zip(row, col)) for col in zip(*rm_y)]
                    for row in rm_x
                ]
                rm_xyz = [
                    [sum(r * c for r, c in zip(row, col)) for col in zip(*rm_z)]
                    for row in rm_xy
                ]

                ### New positions
                new_x = rm_xyz[0][0] * x + rm_xyz[0][1] * y + rm_xyz[0][2] * z
                new_y = rm_xyz[1][0] * x + rm_xyz[1][1] * y + rm_xyz[1][2] * z
                new_z = rm_xyz[2][0] * x + rm_xyz[2][1] * y + rm_xyz[2][2] * z
                
                self.residues[ri].beads[bi].move_atom(new_x, new_y, new_z)
                i += 1
    
    def get_max_length(self):
        beads = self.get_beads()
        max_length = max([math.dist(bead1, bead2) for bead1 in beads for bead2 in beads])
        return max_length

    def get_res_beads_info(self, output_type="ATOM"):
        res_beads_info = []
        for res in self.residues:
            for bead in res.beads:
                if output_type == "ATOM":
                    res_beads_info.append(bead)
                elif output_type in ["tuple", "ziptuple"]:
                    res_beads_info.append(bead.get_tuple())
        if output_type == "ziptuple":
            res_beads_info = tuple(zip(*res_beads_info))
        return res_beads_info
