import time
import math
import numpy as np
import random
from operator import itemgetter
import shapely

class make_rect_grid_lines_iterative_based:
    def make_rect_grid_lines_iterative_based(self, leaflet, subleaflet, occupation_modifier):
        mean_lipid_radius_buffered_multiplier = leaflet["grid_maker_multiplier"]
        
        bbox_polygon = subleaflet["holed_bbox"]
        dict_for_plotting = {
            "LineStringsInfo": [],
        }
        grid_points = []
        grid_points_no_random = []
        total_portion = 0
        total_area_portion = 0

        for lgi, (mean_lipid_radius, nlipids, nlipids_portion, nlipids_area_portion) in enumerate(subleaflet["lipid_groups"]):
            self.print_term("Lipid group i:", lgi, spaces=3, debug=True, debug_keys=["optimizer"])
            self.print_term("mean_lipid_radius:   ", mean_lipid_radius, spaces=4, debug=True, debug_keys=["optimizer"])
            self.print_term("nlipids:             ", nlipids, spaces=4, debug=True, debug_keys=["optimizer"])
            self.print_term("nlipids_portion:     ", nlipids_portion, spaces=4, debug=True, debug_keys=["optimizer"])
            self.print_term("nlipids_area_portion:", nlipids_area_portion, spaces=4, debug=True, debug_keys=["optimizer"])
            total_portion += nlipids_portion
            total_area_portion += nlipids_area_portion
            xmin, xmax, ymin, ymax = itemgetter("xmin", "xmax", "ymin", "ymax")(subleaflet)
            self.print_term("xmin, xmax, ymin, ymax:    ", xmin, xmax, ymin, ymax, spaces=4, debug=True, debug_keys=["optimizer"])

            ### "edge_buffer" is used to buffer all edges
            edge_buffer = mean_lipid_radius + (leaflet["kickxy"] + leaflet["plane_buffer"] * (1+occupation_modifier*2))
            self.print_term("edge_buffer:               ", edge_buffer, spaces=4, debug=True, debug_keys=["optimizer"])

            ### "mean_lipid_radius_buffered" may be reduced in size later while "mean_lipid_radius" must remain a separate and unchanged value
            mean_lipid_radius_buffered = mean_lipid_radius
            self.print_term("mean_lipid_radius_buffered:", mean_lipid_radius_buffered, spaces=4, debug=True, debug_keys=["optimizer"])

            # total_area_portion_limiter = 0.7
            total_area_portion_limiter = leaflet["grid_maker_area_portion_limiter"]
            self.print_term("total_portion:             ", total_portion, spaces=4, debug=True, debug_keys=["optimizer"])
            self.print_term("total_area_portion:        ", total_area_portion, spaces=4, debug=True, debug_keys=["optimizer"])
            self.print_term("total_area_portion_limiter:", total_area_portion_limiter, spaces=4, debug=True, debug_keys=["optimizer"])
            if total_area_portion < total_area_portion_limiter:
                grid_method = "2D_Grid"
            else:
                grid_method = "LineStrings"
            self.print_term("grid_method:               ", grid_method, spaces=4, debug=True, debug_keys=["optimizer"])

            xmin_edge = xmin + edge_buffer
            xmax_edge = xmax - edge_buffer
            ymin_edge = ymin + edge_buffer
            ymax_edge = ymax - edge_buffer

            ### In case membranes are very narrow and edges end up outside the opposide membrane edge
            if xmin_edge >= xmax or xmax_edge <= xmin:
                xmin_edge = xmin
                xmax_edge = xmax

            if ymin_edge >= ymax or ymax_edge <= ymin:
                ymin_edge = ymin
                ymax_edge = ymax
            
            x_edge_len, y_edge_len = xmax_edge-xmin_edge, ymax_edge-ymin_edge
            
            ### Always make sure one is 1 and the other is smaller to ensure "building-up" rather than "tearing-down" of LineStrings
            edge_ratio = (xmax_edge-xmin_edge) / (ymax_edge-ymin_edge)
            if edge_ratio > 1:
                x_edgelen_ratio = 1
                y_edgelen_ratio = 1 / edge_ratio
            else:
                x_edgelen_ratio = 1 / (edge_ratio**(-1))
                y_edgelen_ratio = 1

            self.print_term("edge_ratio:            ", edge_ratio,      spaces=4, debug=True, debug_keys=["optimizer"])
            self.print_term("x_edgelen_ratio:       ", x_edgelen_ratio, spaces=4, debug=True, debug_keys=["optimizer"])
            self.print_term("y_edgelen_ratio:       ", y_edgelen_ratio, spaces=4, debug=True, debug_keys=["optimizer"])

            if grid_method == "LineStrings":
                nlipids_sqrt = math.sqrt(nlipids)
                self.print_term("nlipids_sqrt:          ", nlipids_sqrt, spaces=4, debug=True, debug_keys=["optimizer"])
                x_grid_points_ideal = nlipids_sqrt * x_edgelen_ratio
                y_grid_points_ideal = nlipids_sqrt * y_edgelen_ratio
            elif grid_method == "2D_Grid":
                x_grid_points_ideal = x_edge_len / (mean_lipid_radius*2)
                y_grid_points_ideal = y_edge_len / (mean_lipid_radius*2)
            self.print_term("x_grid_points_ideal:   ", x_grid_points_ideal, spaces=4, debug=True, debug_keys=["optimizer"])
            self.print_term("y_grid_points_ideal:   ", y_grid_points_ideal, spaces=4, debug=True, debug_keys=["optimizer"])
            
            ### Making pointer for when number of points along x/y-axis should be expanded
            xlines_ideal = x_grid_points_ideal
            ylines_ideal = y_grid_points_ideal
            self.print_term("xlines_ideal:          ", xlines_ideal, spaces=4, debug=True, debug_keys=["optimizer"])
            self.print_term("ylines_ideal:          ", ylines_ideal, spaces=4, debug=True, debug_keys=["optimizer"])

            ratio_tot_ideal = xlines_ideal + ylines_ideal
            xlines_ratio_ideal = xlines_ideal/ratio_tot_ideal
            ylines_ratio_ideal = ylines_ideal/ratio_tot_ideal
            self.print_term("xlines_ratio_ideal:    ", xlines_ratio_ideal, spaces=4, debug=True, debug_keys=["optimizer"])
            self.print_term("ylines_ratio_ideal:    ", ylines_ratio_ideal, spaces=4, debug=True, debug_keys=["optimizer"])

            xlines = int(xlines_ideal+0.5) or 1
            ylines = int(ylines_ideal+0.5) or 1
            self.print_term("xlines pre while loop: ", xlines, spaces=4, debug=True, debug_keys=["optimizer"])
            self.print_term("ylines pre while loop: ", ylines, spaces=4, debug=True, debug_keys=["optimizer"])
            
            while xlines * ylines < nlipids:
                new_ratio_tot = xlines + ylines
                if xlines / new_ratio_tot < xlines_ratio_ideal:
                    xlines += 1
                elif ylines / new_ratio_tot < ylines_ratio_ideal:
                    ylines += 1
                else:
                    xlines += 1
            self.print_term("xlines post while loop:", xlines, spaces=4, debug=True, debug_keys=["optimizer"])
            self.print_term("ylines post while loop:", ylines, spaces=4, debug=True, debug_keys=["optimizer"])
            
            bbox_polygon_Buffered = shapely.buffer(bbox_polygon, -mean_lipid_radius_buffered*mean_lipid_radius_buffered_multiplier)

            enough_points = False
            while enough_points == False:
                if grid_method == "LineStrings":
                    xlinespace, xspace = np.linspace(start=xmin_edge, stop=xmax_edge, num=xlines, endpoint=True, retstep=True)
                    ylinespace, yspace = np.linspace(start=ymin_edge, stop=ymax_edge, num=ylines, endpoint=True, retstep=True)

                    LineStrings            = []
                    LineStringsOverlapping = []
                    LineStringsContained   = []
                    
                    for xval in xlinespace:
                        top_point = (xval, ylinespace[-1])
                        bot_point = (xval, ylinespace[0])
                        LineString = shapely.LineString([top_point, bot_point])
                        LineStrings.append(LineString)
                        
                        ### Finding parts of LineStrings that are contained within legal area and at least a certain distance from BBOX
                        LineStringContained = shapely.intersection(LineString, bbox_polygon_Buffered)
                        
                        ### "self.get_geoms_list" contains a debug "print_term"
                        momentary_LineStrings = self.get_geoms_list(LineStringContained)

                        ### Removes "empty" LineStrings
                        LineStringsOverlapping.extend(momentary_LineStrings)
                        
                        ### Checks if any endpoints on LineStrings are too close to each other
                        ### Cuts a portion of each if they are
                        new_momentary_LineStrings = []
                        for LineString in momentary_LineStrings:
                            curr_xs, curr_ys = LineString.xy
                            curr_xmax, curr_xmin = max(curr_xs), min(curr_xs)
                            curr_ymax, curr_ymin = max(curr_ys), min(curr_ys)
                            curr_length = curr_ymax - curr_ymin
                            if len(new_momentary_LineStrings) > 0:
                                y_diff = last_ymin - curr_ymax
                                if y_diff <= yspace:
                                    tot_LS_lengths  = curr_length + last_length
                                    curr_LS_portion = curr_length / tot_LS_lengths
                                    last_LS_portion = last_length / tot_LS_lengths

                                    yspace_diff = yspace - y_diff
                                    curr_cut = curr_LS_portion * yspace_diff
                                    last_cut = last_LS_portion * yspace_diff
                                    if last_cut >= last_length:
                                        ### Remove last LineString if it is too short
                                        new_momentary_LineStrings = new_momentary_LineStrings[:-1]
                                    else:
                                        ### Else cut a part of it off and add it to list
                                        last_ymin += last_cut
                                        last_top_point, last_bot_point = (last_xmax, last_ymax), (last_xmin, last_ymin)
                                        new_last_LS = shapely.LineString([last_top_point, last_bot_point])
                                        new_momentary_LineStrings[-1] = new_last_LS
                                    if curr_cut >= curr_length:
                                        ### If current LineString too short, then just don't do anything
                                        pass
                                    else:
                                        ### Else cut a part of it off and add it to list
                                        curr_ymax -= curr_cut
                                        curr_top_point, curr_bot_point = (curr_xmax, curr_ymax), (curr_xmin, curr_ymin)
                                        new_curr_LS = shapely.LineString([curr_top_point, curr_bot_point])
                                        new_momentary_LineStrings.append(new_curr_LS)
                                else:
                                    new_momentary_LineStrings.append(LineString)
                            else:
                                new_momentary_LineStrings.append(LineString)
                            
                            ### Could happen that all LineStrings have been removed due to being too short.
                            if len(new_momentary_LineStrings) > 0:
                                last_xs, last_ys = new_momentary_LineStrings[-1].xy
                                last_xmax, last_xmin = max(last_xs), min(last_xs)
                                last_ymax, last_ymin = max(last_ys), min(last_ys)
                                last_length = last_ymax - last_ymin
                        
                        LineStringsContained.extend(new_momentary_LineStrings)
                    
                    ngridpoints = 0
                    lines_info = []
                    for LineString in LineStringsContained:
                        line_length = LineString.length
                        npoints_on_line = round(line_length // yspace) + 1
                        ngridpoints += npoints_on_line

                        if npoints_on_line > 1:
                            real_yspace = (line_length - ((npoints_on_line-1)*(mean_lipid_radius*2))) / (npoints_on_line-1)
                        else:
                            real_yspace = (line_length - ((npoints_on_line-1)*(mean_lipid_radius*2))) / 1
                        
                        xs, ys = LineString.xy
                        xval = xs[0] # both x-values are the same
                        line_ymin, line_ymax = min(ys), max(ys)
                        line_xmin, line_xmax = xmin, xmax
                        
                        lines_info.append({
                            "LineString": LineString,
                            "xval": xval,
                            "ymin": line_ymin,
                            "ymax": line_ymax,
                            "line_length": line_length,
                            "npoints_on_line": npoints_on_line,
                            "real_yspace": real_yspace,
                        })

                    if self.plot_grid:
                        dict_for_plotting["LineStringsInfo"].append({
                            "grid_method": grid_method,
                            "bbox_polygon": bbox_polygon,
                            "bbox_polygon_Buffered": bbox_polygon_Buffered,
                            "xlines": xlines,
                            "ylines": ylines,
                            "xspace": xspace,
                            "yspace": yspace,
                            "xmin": line_xmin,
                            "ymin": line_ymin,
                            "xmax": line_xmax,
                            "ymax": line_ymax,
                            "xmax_edge": xmax_edge,
                            "xmin_edge": xmin_edge,
                            "ymax_edge": ymax_edge,
                            "ymin_edge": ymin_edge,
                            "nlipids": nlipids,
                            "ngridpoints": ngridpoints,
                            "ngridpoints_initial": xlines * ylines,
                            "ngridpoints_final": ngridpoints,
                            "LineStrings": LineStrings,
                            "LineStringsOverlapping": LineStringsOverlapping,
                            "LineStringsContained": LineStringsContained,
                        })
            
                elif grid_method == "2D_Grid":
                    xlinespace, xspace = np.linspace(start=xmin_edge, stop=xmax_edge, num=xlines, endpoint=True, retstep=True)
                    ylinespace, yspace = np.linspace(start=ymin_edge, stop=ymax_edge, num=ylines, endpoint=True, retstep=True)

                    GridPoints = [(x, y) for x in xlinespace for y in ylinespace]
                    GridPoints_MultiPoint = shapely.MultiPoint(GridPoints)
                    GridPoints_MultiPoint_Contained = shapely.intersection(GridPoints_MultiPoint, bbox_polygon_Buffered)
                    GridPoints_Contained = [(p.x, p.y) for p in GridPoints_MultiPoint_Contained.geoms]
                    ngridpoints = len(GridPoints_Contained)

                    if self.plot_grid:
                        dict_for_plotting["LineStringsInfo"].append({
                            "grid_method": grid_method,
                            "bbox_polygon": bbox_polygon,
                            "bbox_polygon_Buffered": bbox_polygon_Buffered,
                            "xlines": xlines,
                            "ylines": ylines,
                            "xspace": xspace,
                            "yspace": yspace,
                            "xmin": xmin_edge,
                            "ymin": ymin_edge,
                            "xmax": xmax_edge,
                            "ymax": ymax_edge,
                            "xmax_edge": xmax_edge,
                            "xmin_edge": xmin_edge,
                            "ymax_edge": ymax_edge,
                            "ymin_edge": ymin_edge,
                            "nlipids": nlipids,
                            "ngridpoints": ngridpoints,
                            "ngridpoints_initial": xlines * ylines,
                            "ngridpoints_final": ngridpoints,
                            "GridPoints": GridPoints,
                            "GridPoints_Contained": GridPoints_Contained,
                        })

                if ngridpoints < nlipids:

                    new_ratio_tot = xlines + ylines
                    if xlines / new_ratio_tot < xlines_ratio_ideal:
                        xlines += 1
                    elif ylines / new_ratio_tot < ylines_ratio_ideal:
                        ylines += 1
                    else:
                        xlines += 1

                    ### Momentary fix to never ending 'CREATING LIPID GRID
                    ideal_limiter = 2
                    if xlines > xlines_ideal*ideal_limiter or ylines > ylines_ideal*ideal_limiter:
                        edge_buffer = edge_buffer * 0.75
                        mean_lipid_radius_buffered = mean_lipid_radius_buffered * 0.75

                        xmin_edge = xmin + edge_buffer
                        xmax_edge = xmax - edge_buffer
                        ymin_edge = ymin + edge_buffer
                        ymax_edge = ymax - edge_buffer

                        ### In case membranes are very narrow and edges end up outside the opposide membrane edge
                        if xmin_edge >= xmax or xmax_edge <= xmin:
                            xmin_edge = xmin
                            xmax_edge = xmax

                        if ymin_edge >= ymax or ymax_edge <= ymin:
                            ymin_edge = ymin
                            ymax_edge = ymax
                        
                        x_edge_len, y_edge_len = xmax_edge-xmin_edge, ymax_edge-ymin_edge
                        
                        x_edgelen_ratio = (xmax_edge-xmin_edge) / (ymax_edge-ymin_edge)
                        y_edgelen_ratio = 1

                        if grid_method == "LineStrings":
                            nlipids_sqrt = math.sqrt(nlipids)
                            x_grid_points_ideal = nlipids_sqrt * x_edgelen_ratio
                            y_grid_points_ideal = nlipids_sqrt * y_edgelen_ratio
                        elif grid_method == "2D_Grid":
                            x_grid_points_ideal = x_edge_len / (mean_lipid_radius*2)
                            y_grid_points_ideal = y_edge_len / (mean_lipid_radius*2)
                        
                        ### Making pointer for when number of points along x/y-axis should be expanded
                        xlines_ideal = x_grid_points_ideal
                        ylines_ideal = y_grid_points_ideal

                        ratio_tot_ideal = xlines_ideal + ylines_ideal
                        xlines_ratio_ideal = xlines_ideal/ratio_tot_ideal
                        ylines_ratio_ideal = ylines_ideal/ratio_tot_ideal

                        bbox_polygon_Buffered = shapely.buffer(bbox_polygon, -mean_lipid_radius_buffered*mean_lipid_radius_buffered_multiplier)
                        if self.plot_grid:
                            dict_for_plotting["LineStringsInfo"][-1]["bbox_polygon_Buffered"] = bbox_polygon_Buffered

                else:
                    enough_points = True
                
                if self.plot_grid:
                    dict_for_plotting["LineStringsInfo"][-1]["enough_points"] = enough_points
            
            ### Finding the points that should be removed due to number of lipids
            ### No duplicate index values as "nlipids" is always equal to or smaller than "ngridpoints"
            if grid_method == "LineStrings":
                lines_info = [line for line in lines_info if line["npoints_on_line"] != 0]
                while ngridpoints > nlipids:
                    line_smallest_yspace = min(lines_info, key=lambda line: line["real_yspace"])

                    line_smallest_yspace["npoints_on_line"] -= 1
                    
                    if line_smallest_yspace["npoints_on_line"] > 1:
                        line_smallest_yspace["real_yspace"] = (line_smallest_yspace["line_length"] - ((line_smallest_yspace["npoints_on_line"]-1)*(mean_lipid_radius*2))) / (line_smallest_yspace["npoints_on_line"]-1)
                    else:
                        line_smallest_yspace["real_yspace"] = 0
                        
                    ### Removes lines not containing any grid points
                    lines_info = [line for line in lines_info if line["npoints_on_line"] != 0]

                    ngridpoints = sum([line["npoints_on_line"] for line in lines_info])

                ### Creating the points
                rand_force = 1
                points_for_lipid_poly = []
                for line in lines_info:
                    xval = line["xval"]
                    if line["npoints_on_line"] == 1:
                        yvals = [np.mean([line["ymin"], line["ymax"]])]
                    else:
                        yvals = np.linspace(start=line["ymin"], stop=line["ymax"], num=line["npoints_on_line"], endpoint=True)
                    
                    for yval in yvals:
                        randx = random.uniform(-leaflet["kickxy"]/rand_force, leaflet["kickxy"]/rand_force)
                        randy = random.uniform(-leaflet["kickxy"]/rand_force, leaflet["kickxy"]/rand_force)
                        grid_points.append((xval+randx, yval+randy))
                        grid_points_no_random.append((xval, yval))
                        points_for_lipid_poly.append((xval+randx, yval+randy))

            elif grid_method == "2D_Grid":
                GridSplitting_delimiter = 50
                xsplit_number_of_splits = round(x_edge_len // GridSplitting_delimiter)
                ysplit_number_of_splits = round(y_edge_len // GridSplitting_delimiter)
                ### +2 to number of splits below to add positive and negative endpoints
                xsplit_delimiters, xsplit_space = np.linspace(start=xmin_edge, stop=xmax_edge, num=xsplit_number_of_splits+2, endpoint=True, retstep=True)
                ysplit_delimiters, ysplit_space = np.linspace(start=ymin_edge, stop=ymax_edge, num=ysplit_number_of_splits+2, endpoint=True, retstep=True)
                split_group_delimiters = [
                    (round(xmin, 4), round(xmax, 4), round(ymin, 4), round(ymax, 4))
                    for ymin, ymax in zip(ysplit_delimiters[:-1], ysplit_delimiters[1:])
                    for xmin, xmax in zip(xsplit_delimiters[:-1], xsplit_delimiters[1:])
                ]
                xsplit_number_of_groups = xsplit_number_of_splits + 1
                ysplit_number_of_groups = ysplit_number_of_splits + 1

                split_group_points = [[] for _ in range(len(split_group_delimiters))]
                for px, py in GridPoints_Contained:
                    px_shifted = round(px-xmin_edge, 4)
                    py_shifted = round(py-ymin_edge, 4)
                    xi = round(px_shifted//xsplit_space)
                    yi = round(py_shifted//xsplit_space)
                    if xi > xsplit_number_of_groups-1:
                        xi = xsplit_number_of_groups-1
                    elif xi < 0:
                        xi = 0
                    if yi > ysplit_number_of_groups-1:
                        yi = ysplit_number_of_groups-1
                    elif yi < 0:
                        yi = 0
                    i = xi + yi*ysplit_number_of_groups

                    split_group_points[i].append((px, py))
                
                ### Removes groups containing no points to prevent errors further down
                split_group_points = [i for i in split_group_points if i]

                split_group_occupancy = [len(i) / ngridpoints for i in split_group_points]

                split_group_lipid_distribution = [int(nlipids * i) for i in split_group_occupancy]
                iters = [split_group_occupancy[:]]
                while sum(split_group_lipid_distribution) < nlipids:
                    ### If-Else statement to avoid division by zero errors
                    if sum(split_group_lipid_distribution) > 0:
                        diffs = [
                            (split_group_occupancy[i] - lip_val / sum(split_group_lipid_distribution), i)
                            for i, lip_val in enumerate(split_group_lipid_distribution)
                        ]
                        biggest_diff, biggest_diff_i = max(diffs, key=lambda d: d[0])

                        split_group_lipid_distribution[biggest_diff_i] += 1
                        iters.append(split_group_lipid_distribution[:])
                    else:
                        split_group_lipid_distribution[0] += 1
                        iters.append(split_group_lipid_distribution[:])

                rand_force = 1
                points_for_lipid_poly = []
                for point_group, nlipids_in_group in zip(split_group_points, split_group_lipid_distribution):
                    chosen_points = random.sample(point_group, k = nlipids_in_group)
                    for xval, yval in chosen_points:
                        randx = random.uniform(-leaflet["kickxy"]/rand_force, leaflet["kickxy"]/rand_force)
                        randy = random.uniform(-leaflet["kickxy"]/rand_force, leaflet["kickxy"]/rand_force)
                        grid_points.append((xval+randx, yval+randy))
                        grid_points_no_random.append((xval, yval))
                        points_for_lipid_poly.append((xval+randx, yval+randy))

            lipid_poly = shapely.MultiPoint(points_for_lipid_poly).buffer(mean_lipid_radius)
            bbox_polygon = shapely.difference(bbox_polygon, lipid_poly)

        out_dict = {
            "grid_points":           grid_points,
            "grid_points_no_random": grid_points_no_random,
            "dict_for_plotting":     dict_for_plotting,
            "bbox_polygon":          bbox_polygon
        }

        return out_dict