import time
import random

class lipid_inserter:
    def lipid_inserter(self):
        if len(self.MEMBRANES) != 0:
            lipid_inserter_tic = time.time()
            string = " ".join(["", "CREATING LIPIDS", ""])
            self.print_term("{string:-^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), spaces=0, verbose=1)
            for memb_key, memb_dict in self.MEMBRANES.items():
                self.print_term("Starting membrane nr", memb_key, spaces=0, verbose=2)
                for leaflet_i, (leaflet_key, leaflet) in enumerate(memb_dict["leaflets"].items()):
                    self.print_term("Starting", leaflet["leaflet_type"], "leaflet", spaces=1, verbose=2)
                        
                    if leaflet["HG_direction"] == "up":
                        sign = +1
                    if leaflet["HG_direction"] == "down":
                        sign = -1
                    leaflet["grid_lipids"] = []
                    
                    for (slxi, slyi, sli), subleaflet in leaflet["subleaflets"].items():
                        lipids      = subleaflet["lipids"]
                        grid_points = subleaflet["grid_points"]
                        for (l_name, l_radius), (grid_point_x, grid_point_y, grid_point_z) in zip(lipids, grid_points):
                            new_x, new_y, new_z = [], [], []
                            leaflet["grid_lipids"].append({
                                "grid_point_x": grid_point_x,
                                "grid_point_y": grid_point_y,
                                "grid_point_z": grid_point_z,
                            })
                            random_rotaion_angle = random.uniform(0, 360)
                            for x, y, z in leaflet["lipids"][l_name].get_beads("xyz"):
                                nx, ny, nz = self.rotate_point(x, y, z, 0, 0, random_rotaion_angle)
                                new_x.append(nx)
                                new_y.append(ny)
                                new_z.append(nz)
                            
                            lipid_dict = {
                                "name": l_name,
                                "beads": [bead.bead for bead in leaflet["lipids"][l_name].get_res_beads_info()],
                                "x": [x      + grid_point_x + random.uniform(-leaflet["kickxy"], leaflet["kickxy"]) for x in new_x],
                                "y": [y      + grid_point_y + random.uniform(-leaflet["kickxy"], leaflet["kickxy"]) for y in new_y],
                                "z": [z*sign + grid_point_z + random.uniform(-leaflet["kickz"], leaflet["kickz"]) for z in new_z],
                                "charges": [bead.charge for bead in leaflet["lipids"][l_name].get_res_beads_info()],
                            }
                            leaflet["grid_lipids"][-1].update({"lipid": lipid_dict})
                            self.system_charge += leaflet["lipids"][l_name].get_mol_charge()
            lipid_inserter_toc = time.time()
            lipid_inserter_time = round(lipid_inserter_toc - lipid_inserter_tic, 4)
            string = " ".join(["", "LIPID CREATION COMPLETE", ""])
            self.print_term("{string:-^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), spaces=0, verbose=1)
            string = " ".join(["", "(Time spent:", str(lipid_inserter_time), "[s])", ""])
            self.print_term("{string:^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), "\n", spaces=0, verbose=1)
