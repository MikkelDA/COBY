import time
import math
import numpy as np
import random
from operator import itemgetter
import shapely

class make_rect_grid_lines_based:
    def make_rect_grid_lines_based(self, leaflet, subleaflet, occupation_modifier):
        xmin, xmax, ymin, ymax = itemgetter("xmin", "xmax", "ymin", "ymax")(subleaflet)

        ### "edge_buffer" uses "occupation_modifier*2" while "mean_lipid_radius" uses only "occupation_modifier"
        ### Done to ensure lipids are allowed to be placed right near the border
        edge_buffer       = (leaflet["lipid_dimensions"]["lipid_radius"] + leaflet["kickxy"] + leaflet["plane_buffer"]) * (1+occupation_modifier*2)
        mean_lipid_radius = (leaflet["lipid_dimensions"]["lipid_radius"] + leaflet["kickxy"] + leaflet["plane_buffer"]) * (1+occupation_modifier)
        bbox_polygon      = subleaflet["holed_bbox"]
        lipids            = subleaflet["lipids"]
        
        if leaflet["lipid_optim"] == "abs_val":
            sidelen = mean_lipid_radius*2
        else:
            sidelen = math.sqrt(leaflet["apl"])
        
        xmin_edge = xmin + edge_buffer
        xmax_edge = xmax - edge_buffer
        ymin_edge = ymin + edge_buffer
        ymax_edge = ymax - edge_buffer

        ### In case membranes are very narrow and edges end up outside the opposide membrane edge
        if xmin_edge < xmax or xmax_edge > xmin:
            xmin_edge = xmin
            xmax_edge = xmax

        if ymin_edge < ymax or ymax_edge > ymin:
            ymin_edge = ymin
            ymax_edge = ymax
        
        ### Making pointer for when number of points along x/y-axis should be expanded
        xlines_ideal = (xmax_edge-xmin_edge)/sidelen
        ylines_ideal = (ymax_edge-ymin_edge)/sidelen

        ratio_tot_ideal = xlines_ideal + ylines_ideal
        xlines_ratio_ideal = xlines_ideal/ratio_tot_ideal
        ylines_ratio_ideal = ylines_ideal/ratio_tot_ideal

        xlines = int(xlines_ideal+0.5) or 1
        ylines = int(ylines_ideal+0.5) or 1
        
        while xlines * ylines < len(lipids):
            self.print_term("xlines, ylines:", xlines, ylines, debug=True)
            new_ratio_tot = xlines + ylines
            if xlines / new_ratio_tot < xlines_ratio_ideal:
                xlines += 1
            elif ylines / new_ratio_tot < ylines_ratio_ideal:
                ylines += 1
            else:
                xlines += 1
        
        bbox_polygon_BufferedForLineStrings = shapely.buffer(bbox_polygon, -mean_lipid_radius)
        dict_for_plotting = {
            "bbox_polygon_BufferedForLineStrings": bbox_polygon_BufferedForLineStrings,
            "LineStringsInfo": [],
        }

        enough_points = False
        while enough_points == False:
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
                LineStringContained = shapely.intersection(LineString, bbox_polygon_BufferedForLineStrings)
                
                ### "self.get_geoms_list" contains a debug "print_term"
                momentary_LineStrings = self.get_geoms_list(LineStringContained)

                ### Removes "empty" LineStrings
                LineStringsOverlapping.extend(momentary_LineStrings)
                
                ### Checks if any endpoints on LineStrings are too close to each other
                ### Cuts a portion of each if they are
                new_momentary_LineStrings = []
                for LS, LineString in enumerate(momentary_LineStrings):
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
                    real_yspace = line_length / (npoints_on_line-1)
                else:
                    real_yspace = 0
                
                xs, ys = LineString.xy
                xval = xs[0] # both x-values are the same
                ymin, ymax = min(ys), max(ys)
                
                lines_info.append({
                    "LineString": LineString,
                    "xval": xval,
                    "ymin": ymin,
                    "ymax": ymax,
                    "line_length": line_length,
                    "npoints_on_line": npoints_on_line,
                    "real_yspace": real_yspace,
                })

            dict_for_plotting["LineStringsInfo"].append({
                "xlines": xlines,
                "ylines": ylines,
                "xspace": xspace,
                "yspace": yspace,
                "xmin": xmin,
                "ymin": ymin,
                "xmax": xmax,
                "ymax": ymax,
                "xmax_edge": xmax_edge,
                "xmin_edge": xmin_edge,
                "ymax_edge": ymax_edge,
                "ymin_edge": ymin_edge,
                "ngridpoints": ngridpoints,
                "LineStrings": LineStrings,
                "LineStringsOverlapping": LineStringsOverlapping,
                "LineStringsContained": LineStringsContained,
            })

            if ngridpoints < len(lipids):

                new_ratio_tot = xlines + ylines
                if xlines / new_ratio_tot < xlines_ratio_ideal:
                    xlines += 1
                elif ylines / new_ratio_tot < ylines_ratio_ideal:
                    ylines += 1
                else:
                    xlines += 1

                ### Momentary fix to never ending 'CREATING LIPID GRID
                if xlines > xlines_ideal*1.25 or ylines > ylines_ideal*1.25:
                    edge_buffer = edge_buffer * 0.75
                    mean_lipid_radius = mean_lipid_radius * 0.75

                    xmin_edge = xmin + edge_buffer
                    xmax_edge = xmax - edge_buffer
                    ymin_edge = ymin + edge_buffer
                    ymax_edge = ymax - edge_buffer

                    ### In case membranes are very narrow and edges end up outside the opposide membrane edge
                    if xmin_edge < xmax or xmax_edge > xmin:
                        xmin_edge = xmin
                        xmax_edge = xmax

                    if ymin_edge < ymax or ymax_edge > ymin:
                        ymin_edge = ymin
                        ymax_edge = ymax
                    
                    xlines_ideal = (xmax_edge-xmin_edge)/sidelen
                    ylines_ideal = (ymax_edge-ymin_edge)/sidelen

                    ratio_tot_ideal = xlines_ideal + ylines_ideal
                    xlines_ratio_ideal = xlines_ideal/ratio_tot_ideal
                    ylines_ratio_ideal = ylines_ideal/ratio_tot_ideal

                    bbox_polygon_BufferedForLineStrings = shapely.buffer(bbox_polygon, -mean_lipid_radius)
                    dict_for_plotting["bbox_polygon_BufferedForLineStrings"] = bbox_polygon_BufferedForLineStrings

            else:
                enough_points = True
            
            dict_for_plotting["LineStringsInfo"][-1]["enough_points"] = enough_points
        
        ### Finding the points that should be removed due to number of lipids
        ### No duplicate index values as "len(lipids)" is always equal to or smaller than "ngridpoints"
        
        lines_info = [line for line in lines_info if line["npoints_on_line"] != 0]
        while ngridpoints > len(lipids):
            line_smallest_yspace = min(lines_info, key=lambda line: line["real_yspace"])
            
            line_smallest_yspace["npoints_on_line"] -= 1
            
            if line_smallest_yspace["npoints_on_line"] > 1:
                line_smallest_yspace["real_yspace"] = line_smallest_yspace["line_length"] / (line_smallest_yspace["npoints_on_line"]-1)
            else:
                line_smallest_yspace["real_yspace"] = 0
            
            ### Removes lines not containing any grid points
            lines_info = [line for line in lines_info if line["npoints_on_line"] != 0]

            ngridpoints = sum([line["npoints_on_line"] for line in lines_info])

        ### Creating the points
        grid_points = []
        grid_points_no_random = []
        rand_force = 1
        for line in lines_info:
            xval = line["xval"]
            if line["real_yspace"] == 0:
                yvals = [np.mean([line["ymin"], line["ymax"]])]
            else:
                yvals = np.linspace(start=line["ymin"], stop=line["ymax"], num=line["npoints_on_line"], endpoint=True)
            
            for yval in yvals:
                randx = random.uniform(-leaflet["kickxy"]/rand_force, leaflet["kickxy"]/rand_force)
                randy = random.uniform(-leaflet["kickxy"]/rand_force, leaflet["kickxy"]/rand_force)
                grid_points.append((xval+randx, yval+randy))
                grid_points_no_random.append((xval, yval))

        return grid_points, grid_points_no_random, lines_info, xspace, dict_for_plotting