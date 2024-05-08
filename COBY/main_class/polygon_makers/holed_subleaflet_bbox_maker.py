import time
import shapely
import shapely.affinity
import copy
from COBY.general_functions.__init__ import flatten

import alphashape

class holed_subleaflet_bbox_maker:
            
    def holed_subleaflet_bbox_maker(self):
        if len(self.MEMBRANES) != 0:
            holed_subleaflet_bbox_maker_tic = time.time()
            string = " ".join(["", "CREATING HOLED BOUNDARY BOXES", ""])
            self.print_term("{string:-^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), spaces=0, verbose=1)
            for memb_i, (memb_key, memb_dict) in enumerate(self.MEMBRANES.items()):
                if memb_i != 0:
                    self.print_term("", verbose=2)
                self.print_term("Starting membrane nr", memb_key, spaces=0, verbose=2)
                
                ### Defining some default values for union, intersection and holed_bbox
                for leaflet_i, (leaflet_key, leaflet) in enumerate(memb_dict["leaflets"].items()):
                    leaflet["intersection"]                         = False
                    leaflet["intersection_without_protein_removed"] = False
                    leaflet["holed_bbox"]                           = leaflet["box_poly"]
                    leaflet["holed_bbox_without_protein_removed"]   = leaflet["box_poly"]
                    leaflet["remove_union"]                         = False
                    leaflet["remove_union_without_protein"]         = False
                    leaflet["require_union"]                        = False
                    for string1 in ["points", "Points", "Points_buffered", "Points_buffered_union", "polygon_exterior_points", "ALPHASHAPE", "ConcaveHulls_Polygon"]:
                        for string2 in ["inside_membrane", "surrounding_membrane"]:
                            leaflet["protein"+"_"+string1+"_"+string2] = False
                    for (slxi, slyi, sli), subleaflet in leaflet["subleaflets"].items():
                        subleaflet["intersection"] = False
                        subleaflet["holed_bbox"]   = subleaflet["box_poly"]
                
                ######################################
                ### HANDLING PROTEIN RELATED HOLES ###
                ######################################
                ### Getting all the protein bead positions
                if len(self.PROTEINS) != 0:
                    self.print_term("Finding protein beads inside leaflets", spaces=1, verbose=2)

                    for leaflet_i, (leaflet_key, leaflet) in enumerate(memb_dict["leaflets"].items()):
                        self.print_term("Starting", leaflet["leaflet_type"], "leaflet", spaces=2, verbose=2)

                        ####################################################
                        ### Finding protein points contained in leaflets ###
                        ####################################################

                        proteins_beads_inside_membrane = [
                            (beadx, beady)
                            for protein_i, (protein_nr, protein) in enumerate(self.PROTEINS.items())
                            for beadx, beady, beadz in protein["protein"].get_beads("xyz")
                            if (
                                leaflet["HG_direction"] == "up" and (
                                    leaflet["center"][2] - protein["buffer"] < beadz < leaflet["center"][2] + leaflet["lipid_dimensions"]["lipid_height"] + protein["buffer"]
                                )
                            )
                            or (
                                leaflet["HG_direction"] == "down" and (
                                    leaflet["center"][2] + protein["buffer"] > beadz > leaflet["center"][2] - leaflet["lipid_dimensions"]["lipid_height"] - protein["buffer"]
                                )
                            )
                            if not protein["membrane_border"]
                        ]

                        proteins_beads_surrounding_membrane = [
                            (beadx, beady)
                            for protein_i, (protein_nr, protein) in enumerate(self.PROTEINS.items())
                            for beadx, beady, beadz in protein["protein"].get_beads("xyz")
                            if (
                                leaflet["HG_direction"] == "up" and (
                                    leaflet["center"][2] - protein["buffer"] < beadz < leaflet["center"][2] + leaflet["lipid_dimensions"]["lipid_height"] + protein["buffer"]
                                )
                            )
                            or (
                                leaflet["HG_direction"] == "down" and (
                                    leaflet["center"][2] + protein["buffer"] > beadz > leaflet["center"][2] - leaflet["lipid_dimensions"]["lipid_height"] - protein["buffer"]
                                )
                            )
                            if protein["membrane_border"]
                        ]

                        def calc_protein_alphashapes(protein_beads):
                            Points = [shapely.Point(point) for point in protein_beads]
                            
                            Points_buffered = []
                            for Point in Points:
                                Points_buffered.append(
                                    Point.buffer(
                                        distance = leaflet["prot_buffer"], # [Å]
                                        resolution = 3, # "quad_segs" in the docs
                                    )
                                )
                            
                            Points_buffered_union = shapely.ops.unary_union(Points_buffered)
                            
                            Points_buffered_union_polygons = self.get_geoms_list(Points_buffered_union)
                            
                            poly_points = []
                            for poly in list(Points_buffered_union_polygons):
                                xs, ys = poly.exterior.coords.xy
                                poly_points.extend(list(zip(xs, ys)))
                            
                            alpha = 1 / (leaflet["lipid_dimensions"]["lipid_radius"] * 2 * leaflet["alpha_mult"] + leaflet["prot_buffer"])
                            ### alpha = radius^-1
                            ALPHASHAPE = alphashape.alphashape(points = poly_points, alpha = alpha)
                            self.print_term(type(ALPHASHAPE), debug = True)
                            
                            ### Alphashape output can be a bit unpredictable so need to check all possibilities
                            ConcaveHulls_Polygon = self.get_geoms_list(ALPHASHAPE)
                            
                            # leaflet["ConcaveHulls_Polygon_1"] = ConcaveHulls_Polygon

                            out_dict = {
                                "Points": Points,
                                "Points_buffered": Points_buffered,
                                "Points_buffered_union": Points_buffered_union,
                                "polygon_exterior_points": poly_points,
                                "ALPHASHAPE": ALPHASHAPE,
                                "ConcaveHulls_Polygon": ConcaveHulls_Polygon,
                            }

                            return out_dict
                        
                        if leaflet["inside_protein"]:
                            if len(proteins_beads_inside_membrane) == 0:
                                self.print_term("No proteins beads found within leaflet", spaces=3, verbose=2)
                            else:
                                self.print_term(str(len(proteins_beads_inside_membrane)) + " protein beads found within leaflet", spaces=3, verbose=2)
                                out_dict = calc_protein_alphashapes(proteins_beads_inside_membrane)
                                leaflet["protein_points_inside_membrane"] = proteins_beads_inside_membrane
                                for key, val in out_dict.items():
                                    new_key = "protein" + "_" + key + "_" + "inside_membrane"
                                    leaflet[new_key] = val
                            
                            if len(proteins_beads_surrounding_membrane) == 0:
                                self.print_term("No proteins beads found surrounding leaflet", spaces=3, verbose=2)
                            else:
                                self.print_term(str(len(proteins_beads_surrounding_membrane)) + " protein beads found surrounding leaflet", spaces=3, verbose=2)
                                out_dict = calc_protein_alphashapes(proteins_beads_surrounding_membrane)
                                leaflet["protein_points_surrounding_membrane"] = proteins_beads_surrounding_membrane
                                for key, val in out_dict.items():
                                    new_key = "protein" + "_" + key + "_" + "surrounding_membrane"
                                    leaflet[new_key] = val
                        else:
                            proteins_beads_inside_membrane = proteins_beads_inside_membrane + proteins_beads_surrounding_membrane
                            if len(proteins_beads_inside_membrane) == 0:
                                self.print_term("No proteins beads found within leaflet", spaces=3, verbose=2)
                            else:
                                self.print_term(str(len(proteins_beads_inside_membrane)) + " protein beads found within leaflet", spaces=3, verbose=2)
                                out_dict = calc_protein_alphashapes(proteins_beads_inside_membrane)
                                leaflet["protein_points_inside_membrane"] = proteins_beads_inside_membrane
                                for key, val in out_dict.items():
                                    new_key = "protein" + "_" + key + "_" + "inside_membrane"
                                    leaflet[new_key] = val
                        
                        self.print_term("", verbose=2) # verbose=2 to avoid blank line for verbose < 2
                
                ####################################################
                ### COMBINING SHAPELY SHAPES FOR ALL SUBLEAFLETS ###
                ####################################################
                self.print_term("Calculating holed boundary box", spaces=1, verbose=2)
                for leaflet_i, (leaflet_key, leaflet) in enumerate(memb_dict["leaflets"].items()):
                    self.print_term("Starting", leaflet["leaflet_type"], "leaflet", spaces=2, verbose=2)
                    
                    ### All the various polygons to be removed from the bbox
                    ### Unique to each leaflet but the same for each subleaflet
                    poly_remove_list = []
                    poly_require_list = []
                    
                    ################################
                    ### HANDLING REMOVE POLYGONS ###
                    ################################
                    holes_patches = []
                    if leaflet["hole"]:
                        holes_patches.append((leaflet["hole"], "hole"))
                    if leaflet["patch"]:
                        holes_patches.append((leaflet["patch"], "patch"))
                    
                    for hp_values, hp_type in holes_patches:
                        for settings in hp_values:
                            if settings["shape_type"] in ["polygon"]:
                                ### Polygon must contain at least three points
                                polygon = shapely.Polygon(settings["points"])
                            elif settings["shape_type"] in ["circle", "ellipse"]:
                                polygon = shapely.Point(settings["points"])
                                polygon = polygon.buffer(1, cap_style=1) # round
                            elif settings["shape_type"] in ["square", "rectangle"]:
                                polygon = shapely.Point(settings["points"])
                                polygon = polygon.buffer(1, cap_style=3) # square
                            
                            if settings["xscaling"] != 0 or settings["yscaling"] != 0:
                                ### Needs to correct settings because radius=0 would scale the axis into non-existance
                                if settings["xscaling"] == 0:
                                    settings["xscaling"] = 1
                                if settings["yscaling"] == 0:
                                    settings["yscaling"] = 1
                                ### Expands polygon along specific axes
                                polygon = shapely.affinity.scale(polygon, xfact=settings["xscaling"], yfact=settings["yscaling"], origin="center")
                            if settings["rotate"] != 0:
                                polygon = shapely.affinity.rotate(polygon, settings["rotate"])
                            if settings["buffer"] != 0:
                                polygon = polygon.buffer(settings["buffer"], cap_style=settings["buffer_cap"], join_style=settings["buffer_join"])
                            
                            if hp_type == "hole":
                                poly_remove_list.append(polygon)
                            elif hp_type == "patch":
                                poly_require_list.append(polygon)
                    
                    ### ### Combining polygons into single shape
                    ### First make remove union without protein
                    remove_union_without_protein = shapely.unary_union(poly_remove_list)
                    leaflet["remove_union_without_protein"] = shapely.unary_union(remove_union_without_protein)
                    
                    ### Then add protein and make second remove union
                    if leaflet["protein_ConcaveHulls_Polygon_inside_membrane"]:
                        protein_poly = leaflet["protein_ConcaveHulls_Polygon_inside_membrane"]
                        for val in [protein_poly]:
                            if val:
                                poly_remove_list.extend(val)
                    
                    if leaflet["protein_ConcaveHulls_Polygon_surrounding_membrane"] and not leaflet["inside_protein"]:
                        protein_poly = leaflet["protein_ConcaveHulls_Polygon_surrounding_membrane"]
                        for val in [protein_poly]:
                            if val:
                                poly_remove_list.extend(val)
                    
                    remove_union = shapely.unary_union(poly_remove_list)
                    leaflet["remove_union"] = remove_union
                    
                    ### Make require union
                    if poly_require_list:
                        require_union = shapely.unary_union(poly_require_list)
                        leaflet["require_union"] = require_union
                    
                    ### Leaflet box poly
                    leaflet_box_poly = leaflet["box_poly"]
                    if leaflet["require_union"]:
                        leaflet_box_poly = shapely.intersection(leaflet["require_union"], leaflet_box_poly)

                    if leaflet["inside_protein"] and leaflet["protein_ALPHASHAPE_surrounding_membrane"]:
                        ### Makes polygon contained inside of protein
                        ### Small buffer added to prevent odditities near the surface of the polygons
                        inside_protein   = leaflet["protein_ALPHASHAPE_surrounding_membrane"].difference(leaflet["protein_Points_buffered_union_surrounding_membrane"].buffer(0.2))
                        leaflet_box_poly = leaflet_box_poly.intersection(inside_protein)

                    ### Leaflet holed bbox combined with protein
                    leaflet_intersection                            = shapely.intersection(leaflet["remove_union_without_protein"], leaflet_box_poly)
                    leaflet_holed_bbox                              = shapely.difference(leaflet_box_poly, leaflet_intersection)
                    leaflet["intersection_without_protein_removed"] = leaflet_intersection
                    leaflet["holed_bbox_without_protein_removed"]   = leaflet_holed_bbox
                        
                    ### Leaflet holed bbox with protein removed
                    leaflet_intersection    = shapely.intersection(leaflet["remove_union"], leaflet_box_poly)
                    leaflet_holed_bbox      = shapely.difference(leaflet_box_poly, leaflet_intersection)
                    leaflet["intersection"] = leaflet_intersection
                    leaflet["holed_bbox"]   = leaflet_holed_bbox
                    
                    for (slxi, slyi, sli), subleaflet in leaflet["subleaflets"].items():
                        ### The subleaflet box bbox polygon
                        box_poly = subleaflet["box_poly"]

                        if leaflet["inside_protein"] and leaflet["protein_ALPHASHAPE_surrounding_membrane"]:
                            ### Makes polygon contained inside of protein
                            box_poly = box_poly.intersection(inside_protein)

                        if leaflet["require_union"]:
                            box_poly = shapely.intersection(leaflet["require_union"], box_poly)
                        
                        ### With protein
                        intersection               = shapely.intersection(leaflet["remove_union"], box_poly)
                        holed_bbox                 = shapely.difference(box_poly, intersection)
                        subleaflet["intersection"] = intersection
                        subleaflet["holed_bbox"]   = holed_bbox
                                                
                        ### Readjusts the dimensions of the subleaflet
                        if leaflet["readjust_bbox"]:
                            new_xmin, new_ymin, new_xmax, new_ymax = subleaflet["holed_bbox"].bounds
                            subleaflet["xmin"] = new_xmin
                            subleaflet["ymin"] = new_ymin
                            subleaflet["xmax"] = new_xmax
                            subleaflet["ymax"] = new_ymax
                        
                    ### Splits subleaflet into multiple subleaflets if the individual parts are not touching
                    if leaflet["split_bbox"]:
                        try:
                            ### Following code can easily encounter problems
                            ### Using "try" until all potential errors have been handled 
                            new_subleaflets = {}
                            new_sli = 0
                            for (slxi, slyi, sli), subleaflet in leaflet["subleaflets"].items():
                                geoms = self.get_geoms_list(subleaflet["holed_bbox"])
                                if len(geoms) > 1:
                                    for geom in geoms:
                                        new_subleaflet = copy.deepcopy(subleaflet)
                                        new_subleaflet["original_holed_bbox"] = new_subleaflet["holed_bbox"]
                                        new_subleaflet["holed_bbox"] = geom

                                        if leaflet["readjust_bbox"]:
                                            new_xmin, new_ymin, new_xmax, new_ymax = new_subleaflet["holed_bbox"].bounds
                                            new_subleaflet["xmin"] = new_xmin
                                            new_subleaflet["ymin"] = new_ymin
                                            new_subleaflet["xmax"] = new_xmax
                                            new_subleaflet["ymax"] = new_ymax
                                        
                                        new_subleaflets[(slxi, slyi, new_sli)] = new_subleaflet
                                        new_sli += 1
                                else:
                                    new_subleaflet = copy.deepcopy(subleaflet)
                                    new_subleaflets[(slxi, slyi, new_sli)] = new_subleaflet
                                    new_sli += 1
                            leaflet["subleaflets"] = new_subleaflets
                        except:
                            ### If error encounted then don't split subleaflets
                            self.print_term(
                                "COBY tried to split the subleaflet bbox but encountered an error."+"\n",
                                "This bit of code is prone to errors due to many potential inputs."+"\n",
                                "Will continue running without splitting the bbox."+"\n",
                                "Please report this to https://github.com/MikkelDA/COBY with the entire COBY command that you have used.",
                                warn=True
                            )
                    
                    ### Finally checks that at least one lipid can be placed in each subleaflet
                    ### If a subleaflet is too small then it is removed
                    ### Implemented to prevent errors during lipid calculation
                    new_subleaflets = {}
                    new_sli = 0
                    for (slxi, slyi, sli), subleaflet in leaflet["subleaflets"].items():
                        ### Makes sure number of lipids to be inserted is t at least 1 no matter which rounding method
                        if subleaflet["holed_bbox"].area // leaflet["apl"] >= 1:
                            new_subleaflets[(slxi, slyi, new_sli)] = subleaflet
                            new_sli += 1
                        
                    leaflet["subleaflets"] = new_subleaflets
            
            holed_subleaflet_bbox_maker_toc = time.time()
            holed_subleaflet_bbox_maker_time = round(holed_subleaflet_bbox_maker_toc - holed_subleaflet_bbox_maker_tic, 4)
            string = " ".join(["", "HOLED BOUNDARY BOXES CREATED", ""])
            self.print_term("{string:-^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), spaces=0, verbose=1)
            string = " ".join(["", "(Time spent:", str(holed_subleaflet_bbox_maker_time), "[s])", ""])
            self.print_term("{string:^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), "\n", spaces=0, verbose=1)
