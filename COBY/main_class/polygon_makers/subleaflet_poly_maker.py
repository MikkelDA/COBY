import time
import numpy as np
import shapely

class subleaflet_poly_maker:
    def subleaflet_poly_maker(self):
        if len(self.MEMBRANES) != 0:
            subleaflet_poly_maker_tic = time.time()
            string = " ".join(["", "CREATING LEAFLET BOUNDARY BOXES", ""])
            self.print_term("{string:-^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), spaces=0, verbose=1)

            ### Creating a base polygon to delimit planar (x/y) pbc boundaries
            pbc_xmin          = -self.pbcx/2
            pbc_xmax          = +self.pbcx/2
            pbc_ymin          = -self.pbcy/2
            pbc_ymax          = +self.pbcy/2
            pbc_point_min_min = (pbc_xmin, pbc_ymin) # lower left
            pbc_point_max_min = (pbc_xmax, pbc_ymin) # lower right
            pbc_point_max_max = (pbc_xmax, pbc_ymax) # upper right
            pbc_point_min_max = (pbc_xmin, pbc_ymax) # upper left
            pbc_points        = [pbc_point_min_min, pbc_point_max_min, pbc_point_max_max, pbc_point_min_max]
            pbc_Points        = [shapely.Point(point) for point in pbc_points]
            pbc_base_polygon  = shapely.Polygon(pbc_Points)

            for memb_i, (memb_key, memb_dict) in enumerate(self.MEMBRANES.items()):
                if memb_i != 0:
                    self.print_term("", verbose=2)
                self.print_term("Starting membrane nr", memb_key, spaces=0, verbose=2)
                
                for leaflet_i, (leaflet_key, leaflet) in enumerate(memb_dict["leaflets"].items()):
                    leaflet["subleaflets"] = {}
                    self.print_term("Starting", leaflet["leaflet_type"], "leaflet", spaces=1, verbose=2)
                    self.print_term("Base parameters:", "x=" + str(leaflet["xlength"] / 10) + "nm", "y=" + str(leaflet["ylength"] / 10) + "nm", "APL=" + str(leaflet["apl"] / 100) + "nm^2", spaces=2, verbose=2)
                    
                    xmin, xmax = leaflet["center"][0] - leaflet["xlength"]/2, leaflet["center"][0] + leaflet["xlength"]/2
                    ymin, ymax = leaflet["center"][1] - leaflet["ylength"]/2, leaflet["center"][1] + leaflet["ylength"]/2
                    
                    point_min_min = (xmin, ymin) # lower left
                    point_max_min = (xmax, ymin) # lower right
                    point_max_max = (xmax, ymax) # upper right
                    point_min_max = (xmin, ymax) # upper left
                    
                    box_points    = [point_min_min, point_max_min, point_max_max, point_min_max]
                    box_Points    = [shapely.Point(point) for point in box_points]
                    box_Poly_full = shapely.Polygon(box_Points)

                    ### Ensuring bbox does not exceed pbc boundaries
                    box_Poly = shapely.intersection(box_Poly_full, pbc_base_polygon)

                    leaflet.update({
                        "xmin":       xmin,
                        "xmax":       xmax,
                        "ymin":       ymin,
                        "ymax":       ymax,
                        "box_points": box_points,
                        "box_poly":   box_Poly,
                    })
                    
                    xmin, ymin, xmax, ymax = box_Poly.bounds

                    if leaflet["gridsplits"][0] == "auto":
                        xlength = abs(xmax) + abs(xmin)
                        ylength = abs(ymax) + abs(ymin)
                        xsplits = xlength // leaflet["gridsplits"][1] + 1
                        ysplits = ylength // leaflet["gridsplits"][1] + 1
                    else:
                        xsplits, ysplits = leaflet["gridsplits"]
                    if xsplits > 1:
                        self.print_term("x-axis split into", xsplits, "subleaflets due to axis length or manual designation", spaces=2, verbose=3)
                    if ysplits > 1:
                        self.print_term("y-axis split into", ysplits, "subleaflets due to axis length or manual designation", spaces=2, verbose=3)
                    if xsplits * ysplits > 1:
                        self.print_term("A total of", xsplits*ysplits, "subleaflets have been made for the current leaflet", spaces=2, verbose=3)

                    xpoint_vals = np.linspace(start=xmin, stop=xmax, num=round(xsplits+1), endpoint=True)
                    ypoint_vals = np.linspace(start=ymin, stop=ymax, num=round(ysplits+1), endpoint=True)
                    
                    sli = 0
                    for xi in range(round(xsplits)):
                        for yi in range(round(ysplits)):
                            xmin = xpoint_vals[xi]
                            xmax = xpoint_vals[xi+1]
                            ymin = ypoint_vals[yi]
                            ymax = ypoint_vals[yi+1]
                            
                            point_min_min = (xmin, ymin) # lower left
                            point_max_min = (xmax, ymin) # lower right
                            point_max_max = (xmax, ymax) # upper right
                            point_min_max = (xmin, ymax) # upper left
                            
                            box_points = [point_min_min, point_max_min, point_max_max, point_min_max]
                            box_Points = [shapely.Point(point) for point in box_points]
                            box_Poly   = shapely.Polygon(box_Points)
                            
                            leaflet["subleaflets"][(xi, yi, sli)] = {
                                "xmin":          xmin, # Cutouts change these values
                                "xmax":          xmax, # Cutouts change these values
                                "ymin":          ymin, # Cutouts change these values
                                "ymax":          ymax, # Cutouts change these values
                                "xmin_original": xmin, # Cutouts DO NOT change these values
                                "xmax_original": xmax, # Cutouts DO NOT change these values
                                "ymin_original": ymin, # Cutouts DO NOT change these values
                                "ymax_original": ymax, # Cutouts DO NOT change these values
                                "box_points":    box_points,
                                "box_poly":      box_Poly,
                            }
                            sli += 1
            
            subleaflet_poly_maker_toc = time.time()
            subleaflet_poly_maker_time = round(subleaflet_poly_maker_toc - subleaflet_poly_maker_tic, 4)
            string = " ".join(["", "LEAFLET BOUNDARY BOXES CREATED", ""])
            self.print_term("{string:-^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), spaces=0, verbose=1)
            string = " ".join(["", "(Time spent:", str(subleaflet_poly_maker_time), "[s])", ""])
            self.print_term("{string:^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), "\n", spaces=0, verbose=1)
