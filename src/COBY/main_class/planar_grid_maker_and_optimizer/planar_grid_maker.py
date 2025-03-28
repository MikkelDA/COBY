import time
import math
import numpy as np
from operator import itemgetter
import random
import shapely

from COBY.general_functions.flatten import flatten

class planar_grid_maker:
    def planar_grid_maker(self):
        self.SYS_lipids_dict = {}
        self.GRID_PLOTTING = {}
        if len(self.MEMBRANES) != 0:
            grid_making_tic = time.time()
            string = " ".join(["", "CREATING LIPID GRID", ""])
            self.print_term("{string:-^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), spaces=0, verbose=1)
            for memb_i, (memb_key, memb_dict) in enumerate(self.MEMBRANES.items()):
                if memb_i != 0:
                    self.print_term("", verbose=2)
                self.print_term("Starting membrane nr", memb_key, spaces=0, verbose=2)
                for leaflet_i, (leaflet_key, leaflet) in enumerate(memb_dict["leaflets"].items()):
                    if leaflet_i != 0:
                        self.print_term("", verbose=2)
                    self.print_term("Starting", leaflet["leaflet_type"], "leaflet", spaces=1, verbose=2)
                    
                    ### Only one subleaflet or print for all if multiple
                    ### Specifically not if self.verbose is greater than 2 and multiple subleaflets are present.
                    if (self.verbose == 2 and len(leaflet["subleaflets"]) > 1) or (self.verbose >= 3 and len(leaflet["subleaflets"]) == 1):
                        self.print_term("Starting lipid point creation", spaces=2, verbose=2)
                    
                    ### Whether or not to print time steps during optimization
                    if (self.verbose == 2 and len(leaflet["subleaflets"]) == 1) or self.verbose >= 3:
                        optimizer_print_to_term = True
                    else:
                        optimizer_print_to_term = False
                    
                    for sli, ((slxi, slyi, sli), subleaflet) in enumerate(leaflet["subleaflets"].items()):
                        if self.verbose >= 3 and len(leaflet["subleaflets"]) > 1:
                            if sli != 0:
                                self.print_term("", verbose=3)
                            self.print_term("Starting lipid point creation for subleaflet nr", sli+1, spaces=2, verbose=3)
                        if self.plot_grid:
                            self.GRID_PLOTTING[(memb_key, leaflet_key, slxi, slyi, sli)] = {}
                        radii = [
                            leaflet["lipids"][name].get_radius(AXs="xy") + leaflet["kickxy"] + leaflet["plane_buffer"]
                            for name in subleaflet["lipid_names"]
                        ]
                        lipid_names_nlipids_radii = [
                            (
                                name,
                                ratio,
                                radius,
                            )
                            for name, ratio, radius in zip(subleaflet["lipid_names"], subleaflet["lipid_ratios"], radii)
                        ]
                        
                        initial_grid_maker_tic = time.time()
                        
                        leaflet_area        = subleaflet["holed_bbox"].area
                        lipids_circle_area  = sum([(math.pi*(radius**2))*ratio for name, ratio, radius in lipid_names_nlipids_radii])
                        lipids_square_area  = sum([((radius*2)**2)*ratio for name, ratio, radius in lipid_names_nlipids_radii])
                        lipids_mean_area    = (lipids_circle_area + lipids_square_area)/2
                        occupation_modifier = (leaflet_area-lipids_square_area)/(leaflet_area*2) # *2 to half the modifier
                        ### Following to prevent negative values
                        if occupation_modifier < 0:
                            occupation_modifier = 0
                        self.print_term("leaflet_area       ", leaflet_area,        debug=True, debug_keys=["lipid_grid_creation"])
                        self.print_term("lipids_circle_area ", lipids_circle_area,  debug=True, debug_keys=["lipid_grid_creation"])
                        self.print_term("lipids_square_area ", lipids_square_area,  debug=True, debug_keys=["lipid_grid_creation"])
                        self.print_term("lipids_mean_area   ", lipids_mean_area,    debug=True, debug_keys=["lipid_grid_creation"])
                        self.print_term("occupation_modifier", occupation_modifier, debug=True, debug_keys=["lipid_grid_creation"])
                        
                        if occupation_modifier < 0.2:
                            self.print_term(
                                "WARNING:",
                                "Chosen apl ("+str(leaflet["apl"]/100)+" [nm^2]) and buffer ("+str(leaflet["plane_buffer"]/10)+" [nm]) cause lipids to be very closely packed.",
                                "Optimization may be slow and some lipid overlaps may not be avoidable.",
                                "Consider increasing the apl using the subargument 'apl' or decreasing the plane buffer using the subargument 'plane_buffer'",
                                warn=True
                            )
                        
                        if leaflet["grid_maker_algorithm"] == "3D_matrix":
                            grid_points, lipids, dict_for_plotting = self.make_rect_grid_3D_matrix_based(leaflet, subleaflet, lipid_names_nlipids_radii, occupation_modifier)
                            dict_for_plotting["grid_maker_algorithm"] = "3D_matrix"
                            subleaflet["lipids"]  = lipids
                            grid_points_no_random = grid_points
                            if leaflet["optimize_run"] in [True, "auto"]:
                                optimize_run = True
                            else:
                                optimize_run = False
                        
                        elif leaflet["grid_maker_algorithm"] in ["no_groups", "iterative_groups"]:
                            min_radius, max_radius = min(radii), max(radii)
                            if leaflet["grid_maker_algorithm"] == "iterative_groups":
                                separator = leaflet["grid_maker_separator"]
                            elif leaflet["grid_maker_algorithm"] == "no_groups":
                                separator = 100000 # No lipids should be larger than this anyways

                            ### *0.49 instead of 0.5 because of float problems
                            n_radius_groups         = math.ceil(((max_radius+separator*0.49) - (min_radius-separator*0.49)) / separator)
                            radius_group_separators = np.linspace(start=min_radius-separator*0.5, stop=max_radius+separator*0.5, endpoint=True, num=n_radius_groups+1)
                            radius_groups           = list(reversed([(i, j) for i, j in zip(radius_group_separators[:-1], radius_group_separators[1:])]))
                            lipid_radius_groups     = [[] for i in radius_groups]

                            self.print_term("leaflet['grid_maker_algorithm']:", leaflet["grid_maker_algorithm"], spaces=3, debug=True, debug_keys=["lipid_grid_creation"])
                            self.print_term("separator:                      ", separator,                       spaces=3, debug=True, debug_keys=["lipid_grid_creation"])
                            self.print_term("n_radius_groups:                ", n_radius_groups,                 spaces=3, debug=True, debug_keys=["lipid_grid_creation"])
                            self.print_term("radii:                          ", radii,                           spaces=3, debug=True, debug_keys=["lipid_grid_creation"])
                            self.print_term("radius_group_separators:        ", radius_group_separators,         spaces=3, debug=True, debug_keys=["lipid_grid_creation"])
                            self.print_term("radius_groups:                  ", radius_groups,                   spaces=3, debug=True, debug_keys=["lipid_grid_creation"])

                            ### Iterates over list sorted according to radius
                            total_lipids = 0
                            total_lipid_area = 0
                            for name, nlipids, radius in sorted(lipid_names_nlipids_radii, key=lambda x: (x[2], x[1], x[0]), reverse=True):
                                for rgi, (rmin, rmax) in enumerate(radius_groups):
                                    if rmin < radius <= rmax:
                                        lipid_radius_groups[rgi].append((name, nlipids, radius))
                                        total_lipids += nlipids
                                        total_lipid_area += nlipids * (math.pi * radius**2)

                            lipids = []
                            lipid_groups = []
                            for lgi, lipid_group in enumerate(lipid_radius_groups):
                                ### Checks if any lipid types in group
                                if lipid_group: # if list not empty
                                    nlipids_in_group = sum([nlipids for name, nlipids, radius in lipid_group])
                                    ### Skips group if total number of lipids in group is zero
                                    if nlipids_in_group == 0:
                                        continue
                                    mean_radius = sum([nlipids*radius for name, nlipids, radius in lipid_group]) / nlipids_in_group
                                    lipid_groups.append((mean_radius, nlipids_in_group, nlipids_in_group/total_lipids, (nlipids_in_group * (math.pi * mean_radius**2))/total_lipid_area))

                                    lipids_in_lipid_group = []
                                    for name, nlipids, radius in lipid_group:
                                        lipids_in_lipid_group.append([(name, radius) for _ in range(nlipids)])
                                    if leaflet["grid_maker_lipid_distribution"] == "evenly":
                                        lipids_in_lipid_group = self.n_list_mixer(*lipids_in_lipid_group)
                                    elif leaflet["grid_maker_lipid_distribution"] == "random":
                                        lipids_in_lipid_group = flatten(lipids_in_lipid_group)
                                        random.shuffle(lipids_in_lipid_group)

                                    lipids.extend(lipids_in_lipid_group)

                                    self.print_term("Number of lipids in group nr {lgi}: {nlipids}".format(lgi=lgi, nlipids=len(lipids_in_lipid_group)), spaces=3, debug=True, debug_keys=["lipid_grid_creation"])
                            self.print_term("lipid_groups:", lipid_groups, spaces=3, debug=True, debug_keys=["lipid_grid_creation"])

                            subleaflet["lipids"]           = lipids
                            subleaflet["lipid_groups"]     = lipid_groups
                            subleaflet["total_lipids"]     = total_lipids
                            subleaflet["total_lipid_area"] = total_lipid_area

                            out_dict = self.make_rect_grid_lines_iterative_based(leaflet, subleaflet, occupation_modifier)

                            grid_points                               = out_dict["grid_points"]
                            grid_points_no_random                     = out_dict["grid_points_no_random"]
                            dict_for_plotting                         = out_dict["dict_for_plotting"]
                            dict_for_plotting["grid_maker_algorithm"] = "lines"

                            self.print_term('leaflet["optimize_run"]', leaflet["optimize_run"], spaces=1, debug=True, debug_keys=["lipid_grid_creation"])
                            if leaflet["optimize_run"] == True:
                                optimize_run = True
                            elif leaflet["optimize_run"] == "auto":
                                def coord_to_indices(pos, dim, center, real_gridres):
                                    '''
                                    Converts a coordinate to an index value
                                    '''
                                    posi_dec = (pos+dim/2-center)/real_gridres
                                    posi_int = int(posi_dec)
                                    if posi_int < 0:
                                        posi_int = 0
                                    return posi_int

                                max_space = max(radii)*2.2
                                self.print_term("Maximum space:", max_space,   spaces=3, debug=True, debug_keys=["lipid_grid_creation", "optimizer"])

                                xmin, xmax, ymin, ymax = itemgetter("xmin", "xmax", "ymin", "ymax")(subleaflet)
                                xlen, ylen = xmax - xmin, ymax - ymin

                                xsplits   = round(xlen // max_space)
                                xnbins    = xsplits-1
                                xbins_len = xlen / xnbins

                                ysplits   = round(ylen // max_space)
                                ynbins    = ysplits-1
                                ybins_len = ylen / ynbins

                                checkgrid_no_offset   = np.zeros((xnbins,   ynbins),   dtype=int)
                                checkgrid_half_offset = np.zeros((xnbins-1, ynbins-1), dtype=int)

                                ### Checks if two points are too close to each other or if a point is too close to the bbox boundary.
                                holed_bbox = subleaflet["holed_bbox"]
                                all_contained = True
                                for (px, py), (lipid_name, radius) in zip(grid_points, lipids):
                                    pxi_no_offset = coord_to_indices(px, xlen, (xmax+xmin)/2, xbins_len)
                                    pyi_no_offset = coord_to_indices(py, ylen, (ymax+ymin)/2, ybins_len)
                                    checkgrid_no_offset[pxi_no_offset][pyi_no_offset] += 1

                                    pxi_half_offset = coord_to_indices(px, xlen-max_space, (xmax+xmin)/2, xbins_len)
                                    pyi_half_offset = coord_to_indices(py, ylen-max_space, (ymax+ymin)/2, ybins_len)
                                    if 0 <= pxi_half_offset <= len(checkgrid_half_offset)-1 and 0 <= pyi_half_offset <= len(checkgrid_half_offset[pxi_half_offset])-1:
                                        checkgrid_half_offset[pxi_half_offset][pyi_half_offset] += 1

                                    Point = shapely.Point(px, py)
                                    if not holed_bbox.contains(Point) or holed_bbox.boundary.distance(Point) < radius:
                                        all_contained = False
                                        break

                                self.print_term("Max of grid with no offset:  ", np.max(checkgrid_no_offset),   spaces=3, debug=True, debug_keys=["lipid_grid_creation", "optimizer"])
                                self.print_term("Max of grid with half offset:", np.max(checkgrid_half_offset), spaces=3, debug=True, debug_keys=["lipid_grid_creation", "optimizer"])
                                self.print_term("All lipids fully contained:  ", all_contained,                 spaces=3, debug=True, debug_keys=["lipid_grid_creation", "optimizer"])
                                if np.max(checkgrid_no_offset) > 1 or np.max(checkgrid_half_offset) > 1:
                                    optimize_run = True
                                else:
                                    optimize_run = False

                                if not all_contained:
                                    optimize_run = True

                            else:
                                ### "iterative_groups" with one lipid group is identical to "lines_single" so no need to warn
                                if leaflet["grid_maker_algorithm"] == "iterative_groups" and len(lipid_groups) > 1:
                                    self.print_term(
                                        " ".join([
                                            "Using 'grid_maker_algorithm' algorithm 'iterative_groups' but 'optimize_run' is set to 'False'.",
                                            "Lipids are VERY likely to overlap with 'iterative_groups' without optimization.",
                                            "We strongly suggest you set 'optimize_run' to 'True' or 'auto'.",
                                        ]),
                                        warn=True,
                                    )
                                optimize_run = False
                            self.print_term("optimize_run:", optimize_run, spaces=3, debug=True, debug_keys=["lipid_grid_creation", "optimizer"])
                        
                        initial_grid_maker_toc = time.time()
                        initial_grid_maker_time = round(initial_grid_maker_toc - initial_grid_maker_tic, 4)
                        if optimizer_print_to_term:
                            self.print_term("Initial grid creation time:", initial_grid_maker_time, "[s]", spaces=3, verbose=2)

                        grid_points_arr = np.asarray(grid_points)
                        if self.plot_grid:
                            for key, vals in dict_for_plotting.items():
                                self.GRID_PLOTTING[(memb_key, leaflet_key, slxi, slyi, sli)][key] = vals
                        
                        z_values_arr = np.asarray([[leaflet["center"][2]] for _ in range(len(grid_points_arr))])

                        optimize_time_tic = time.time()
                        if optimize_run:
                            optimized_grid_points_arr, POINT_STEPS, step, steps_time, mean_steps_time, max_push, bin_size = self.plane_grid_point_optimizer(
                                grid_points_arr,
                                lipid_sizes = list(zip(*lipids))[-1],
                                polygon = subleaflet["holed_bbox"],
                                
                                xcenter = np.mean([subleaflet["xmax"], subleaflet["xmin"]]),
                                ycenter = np.mean([subleaflet["ymax"], subleaflet["ymin"]]),
                                
                                xlen = subleaflet["xmax"] - subleaflet["xmin"],
                                ylen = subleaflet["ymax"] - subleaflet["ymin"],
                                
                                max_steps             = leaflet["optimize_max_steps"],
                                push_tolerance        = leaflet["optimize_push_tolerance"],
                                lipid_push_multiplier = leaflet["optimize_lipid_push_multiplier"],
                                edge_push_multiplier  = leaflet["optimize_edge_push_multiplier"],
                                
                                occupation_modifier = occupation_modifier,
                                
                                optimizer_print_to_term = optimizer_print_to_term,
                            )
                            
                            subleaflet["grid_points"] = np.hstack([optimized_grid_points_arr, z_values_arr])
                            
                            if self.plot_grid:
                                generated_keys = []
                                for string1 in ["points", "Points", "Points_buffered", "Points_buffered_union", "polygon_exterior_points", "ALPHASHAPE", "ConcaveHulls_Polygon"]:
                                    for string2 in ["inside_membrane", "surrounding_membrane"]:
                                        generated_keys.append("protein"+"_"+string1+"_"+string2)

                                for key in generated_keys + ["remove_union", "require_union", "holed_bbox_without_protein_removed"]:
                                    self.GRID_PLOTTING[(memb_key, leaflet_key, slxi, slyi, sli)].update({
                                        key: leaflet[key]
                                    })
                                
                                for key in ["original_holed_bbox"]:
                                    if key in subleaflet and subleaflet[key]:
                                        self.GRID_PLOTTING[(memb_key, leaflet_key, slxi, slyi, sli)].update({
                                            key: subleaflet[key]
                                        })
                                    
                                self.GRID_PLOTTING[(memb_key, leaflet_key, slxi, slyi, sli)].update({
                                    ### Inputs
                                    "lipids"    : lipids,
                                    "holed_bbox": subleaflet["holed_bbox"],
                                    
                                    "xdims"                : (subleaflet["xmin"], subleaflet["xmax"]),
                                    "ydims"                : (subleaflet["ymin"], subleaflet["ymax"]),
                                    "xdims_original"       : (subleaflet["xmin_original"], subleaflet["xmax_original"]),
                                    "ydims_original"       : (subleaflet["ymin_original"], subleaflet["ymax_original"]),
                                    "lipid_push_multiplier": leaflet["optimize_lipid_push_multiplier"],
                                    "edge_push_multiplier" : leaflet["optimize_edge_push_multiplier"],
                                    "apl"                  : leaflet["apl"],
                                    "bin_size"             : bin_size,

                                    ### Outputs
                                    "grid_points"          : grid_points_arr,
                                    "grid_points_no_random": grid_points_no_random,
                                    "points_arr"           : optimized_grid_points_arr,
                                    "POINT_STEPS"          : POINT_STEPS,
                                    "step"                 : step,
                                    "steps_time"           : steps_time,
                                    "mean_steps_time"      : mean_steps_time,
                                    "max_push"             : max_push,
                                })
                        else:
                            subleaflet["grid_points"] = np.hstack([grid_points_arr, z_values_arr])
                    optimize_time_toc = time.time()
                    optimize_time_time = round(optimize_time_toc - optimize_time_tic, 4)
                    ### Only print if optimizer does not print anything
                    if not optimizer_print_to_term:
                        self.print_term("Total leaflet optimization time:", optimize_time_time, "[s]", spaces=3, verbose=2)
                        
            grid_making_toc = time.time()
            grid_making_time = round(grid_making_toc - grid_making_tic, 4)
            string = " ".join(["", "LIPID GRID CREATED", ""])
            self.print_term("{string:-^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), spaces=0, verbose=1)
            string = " ".join(["", "(Time spent:", str(grid_making_time), "[s])", ""])
            self.print_term("{string:^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), "\n", spaces=0, verbose=1)
