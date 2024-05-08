import time
import math
import numpy as np
import random

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
                    
                    ### Print for each individual
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
                                # leaflet["lipids"][name].get_radius(AXs="xy") + leaflet["kickxy"] + leaflet["plane_buffer"],
                                radius,
                            )
                            for name, ratio, radius in zip(subleaflet["lipid_names"], subleaflet["lipid_ratios"], radii)
                        ]
                        
                        initial_grid_maker_tic = time.time()
                        
                        leaflet_area = subleaflet["holed_bbox"].area
                        lipids_circle_area = sum([(math.pi*(radius**2))*ratio for name, ratio, radius in lipid_names_nlipids_radii])
                        lipids_square_area = sum([((radius*2)**2)*ratio for name, ratio, radius in lipid_names_nlipids_radii])
                        lipids_mean_area = (lipids_circle_area + lipids_square_area)/2
                        occupation_modifier = (leaflet_area-lipids_square_area)/(leaflet_area*2) # *2 to half the modifier
                        self.print_term("leaflet_area       ", leaflet_area, debug=True)
                        self.print_term("lipids_circle_area ", lipids_circle_area, debug=True)
                        self.print_term("lipids_square_area ", lipids_square_area, debug=True)
                        self.print_term("lipids_mean_area   ", lipids_mean_area, debug=True)
                        self.print_term("occupation_modifier", occupation_modifier, debug=True)
                        
                        if occupation_modifier < 0.2:
                            self.print_term(
                                "WARNING:",
                                "Chosen apl ("+str(leaflet["apl"]/100)+" [nm^2]) and buffer ("+str(leaflet["plane_buffer"]/10)+" [nm]) cause lipids to be very closely packed.",
                                "Optimization may be slow and some lipid overlaps may not be avoidable.",
                                "Consider increasing the apl using the subcommand 'apl' or decreasing the plane buffer using the subcommand 'plane_buffer'",
                                warn=True
                            )
                        
                        if leaflet["grid_maker"] == "3D_matrix":
                            grid_points, lipids, dict_for_plotting = self.make_rect_grid_3D_matrix_based(leaflet, subleaflet, lipid_names_nlipids_radii, occupation_modifier)
                            dict_for_plotting["grid_maker"] = "3D_matrix"
                            subleaflet["lipids"]  = lipids
                            grid_points_no_random = grid_points
                            if leaflet["optimize"] == "auto":
                                optimize = "v5"
                                
                        else:
                            lipids = []
                            for name, nlipids, radius in lipid_names_nlipids_radii:
                                lipids.append([(name, radius) for _ in range(nlipids)])

                            if leaflet["lipid_distribution"] == "evenly":
                                lipids = self.n_list_mixer(*lipids)
                            elif leaflet["lipid_distribution"] == "random":
                                lipids = flatten(lipids)
                                random.shuffle(lipids)
                            
                            subleaflet["lipids"] = lipids

                            grid_points, grid_points_no_random, lines_info, xspace, dict_for_plotting = self.make_rect_grid_lines_based(leaflet, subleaflet, occupation_modifier)
                            dict_for_plotting["grid_maker"] = "lines"

                            if leaflet["optimize"] == "auto":
                                
                                max_space = max(radii)*2.2
                                min_yspace = min([line["real_yspace"] for line in lines_info])
                                self.print_term("Maximum space:    ", max_space, debug=True)
                                self.print_term("Universal x space:", xspace, debug=True)
                                self.print_term("Minimum y space:  ", min_yspace, debug=True)
                                if max_space > xspace:
                                    optimize = "v5"
                                    if (self.verbose == 2 and len(leaflet["subleaflets"]) == 1):
                                        self.print_term("Spacing between x-values too small. Starting lipid optimization.", spaces=3, verbose=2)
                                    elif self.verbose >= 3:
                                        self.print_term("Spacing between x-values too small. Starting lipid optimization.", spaces=3, verbose=3)

                                else:
                                    if max_space > min_yspace:
                                        optimize = "v5"
                                        if (self.verbose == 2 and len(leaflet["subleaflets"]) == 1):
                                            self.print_term("Spacing between y-values too small. Starting lipid optimization.", spaces=3, verbose=2)
                                        elif self.verbose >= 3:
                                            self.print_term("Spacing between y-values too small. Starting lipid optimization.", spaces=3, verbose=3)

                                    else:
                                        optimize = False
                                
                            else:
                                optimize = leaflet["optimize"]


                                # # else:
                                # #     xvals = []
                                # #     ylines = []
                                # #     last_xval = False
                                # #     for line in lines_info:
                                # #         if last_xval is False:
                                # #             ylines.append([])
                                # #             last_xval = line["xval"]
                                # #             xvals.append(line["xval"])
                                # #         else:
                                # #             ### Can't use "==" due to floats. Must approximate two floats equalling each other.
                                # #             ### Two x-values should never be this close to each other anyways.
                                # #             if abs(last_xval - line["xval"]) > 0.01:
                                # #                 ylines.append([])
                                # #                 last_xval = line["xval"]
                                # #                 xvals.append(line["xval"])
                                        
                                # #         ylines[-1].append(line)
                        
                        initial_grid_maker_toc = time.time()
                        initial_grid_maker_time = round(initial_grid_maker_toc - initial_grid_maker_tic, 4)
                        if optimizer_print_to_term:
                            self.print_term("Initial grid creation time:", initial_grid_maker_time, "[s]", spaces=3, verbose=2)

                        grid_points_arr = np.asarray(grid_points)
                        if self.plot_grid:
                            for key, vals in dict_for_plotting.items():
                                self.GRID_PLOTTING[(memb_key, leaflet_key, slxi, slyi, sli)][key] = vals
                        
                        z_values_arr = np.asarray([[leaflet["center"][2]] for _ in range(len(grid_points_arr))])

                        optim_time_tic = time.time()
                        if optimize:
                            optimized_grid_points_arr, POINT_STEPS, step, steps_time, mean_steps_time, max_push, bin_size = self.plane_grid_point_optimizer(
                                grid_points_arr,
                                lipid_sizes = list(zip(*lipids))[-1],
                                polygon = subleaflet["holed_bbox"],
                                
                                xcenter = np.mean([subleaflet["xmax"], subleaflet["xmin"]]),
                                ycenter = np.mean([subleaflet["ymax"], subleaflet["ymin"]]),
                                
                                xlen = subleaflet["xmax"] - subleaflet["xmin"],
                                ylen = subleaflet["ymax"] - subleaflet["ymin"],
                                
                                maxsteps = leaflet["optim_maxsteps"],
                                push_tolerance = leaflet["optim_push_tol"],
                                push_mult = leaflet["optim_push_mult"],
                                buffer = leaflet["plane_buffer"],
                                
                                occupation_modifier = occupation_modifier,
                                optimize = leaflet["optimize"],
                                
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
                                    
                                    "xdims"          : (subleaflet["xmin"], subleaflet["xmax"]),
                                    "ydims"          : (subleaflet["ymin"], subleaflet["ymax"]),
                                    "xdims_original" : (subleaflet["xmin_original"], subleaflet["xmax_original"]),
                                    "ydims_original" : (subleaflet["ymin_original"], subleaflet["ymax_original"]),
                                    "push_mult"      : leaflet["optim_push_mult"],
                                    "apl"            : leaflet["apl"],
                                    "bin_size"       : bin_size,

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
                    optim_time_toc = time.time()
                    optim_time_time = round(optim_time_toc - optim_time_tic, 4)
                    ### Only print if optimizer does not print anything
                    if not optimizer_print_to_term:
                        self.print_term("Total leaflet optimization time:", optim_time_time, "[s]", spaces=3, verbose=2)
                        
            grid_making_toc = time.time()
            grid_making_time = round(grid_making_toc - grid_making_tic, 4)
            string = " ".join(["", "LIPID GRID CREATED", ""])
            self.print_term("{string:-^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), spaces=0, verbose=1)
            string = " ".join(["", "(Time spent:", str(grid_making_time), "[s])", ""])
            self.print_term("{string:^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), "\n", spaces=0, verbose=1)
