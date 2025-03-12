import numpy as np
import random
from operator import itemgetter
import shapely

class make_rect_grid_3D_matrix_based:
    def make_rect_grid_3D_matrix_based(self, leaflet, subleaflet, lipid_names_nlipids_radii, occupation_modifier):
        xmin, xmax, ymin, ymax = itemgetter("xmin", "xmax", "ymin", "ymax")(subleaflet)

        dict_for_plotting = {"lipid_insertion": {"lipid_insertion_steps": []}}

        bbox_polygon      = subleaflet["holed_bbox"]

        lipids_sorted_by_size = list(sorted(lipid_names_nlipids_radii, key=lambda x: (x[2], x[1], x[0]), reverse=True))
        lipids = []
        for name, nlipids, radius in lipids_sorted_by_size:
            lipids.extend([(name, radius) for _ in range(nlipids)])
        total_number_of_lipids = len(lipids)

        lipids_sorted_by_size_strs = [("Name", "Number", "Radius")] + [[str(j) for j in i] for i in lipids_sorted_by_size]
        max_str_lengths = [max([len(j) for j in i]) for i in zip(*lipids_sorted_by_size_strs)]
        self.print_term(
            "Lipid information:",
            spaces=3,
            verbose=3,
        )
        for lipid in  lipids_sorted_by_size_strs:
            for_printer = []
            for si, substring in enumerate(lipid):
                for_printer.append('{0: <{L}}'.format(substring, L = max_str_lengths[si]))
            self.print_term(
                *for_printer,
                spaces=4,
                verbose=3,
            )


        radii = [radius for name, nlipid, radius in lipid_names_nlipids_radii]
        radii_sorted_by_size = list(sorted(radii, reverse=True))

        gridres = min(radii_sorted_by_size)/4
        dict_for_plotting["lipid_insertion"]["gridres"] = gridres

        xpoints, ypoints = int((xmax-xmin)/gridres), int((ymax-ymin)/gridres)

        ### Calculates actual coordinate ranges for each axis and calculates the "real" grid resolution
        xcoords, xreal_gridres = np.linspace(xmin-gridres/2, xmax+gridres/2, xpoints+2, retstep=True)
        ycoords, yreal_gridres = np.linspace(ymin-gridres/2, ymax+gridres/2, ypoints+2, retstep=True)
        ### Removes the first and last points as they are the actual edges of the box
        xcoords = xcoords[1:-1]
        ycoords = ycoords[1:-1]

        nlipidtypes = len(lipid_names_nlipids_radii)

        grid_bool_matrix = np.ones((nlipidtypes, xpoints, ypoints), dtype=np.int8)

        dict_for_plotting["lipid_insertion"]["gridres"] = gridres
        dict_for_plotting["lipid_insertion"]["xreal_gridres"] = xreal_gridres
        dict_for_plotting["lipid_insertion"]["yreal_gridres"] = yreal_gridres
        dict_for_plotting["lipid_insertion"]["xcoords"] = xcoords
        dict_for_plotting["lipid_insertion"]["ycoords"] = ycoords

        self.print_term("Number of grid points:", [np.sum(grid_bool_matrix[i,:,:]) for i in range(len(radii_sorted_by_size))], "Initial", spaces=4, verbose=3)

        def coord_to_indices(pos, dim, center, real_gridres, min_buffer, max_buffer):
            ### Buffer limit
            bead_min = pos-min_buffer
            bead_max = pos+max_buffer

            ### Convert solvent buffer limit to decimal index
            beadi_min_dec = (bead_min+dim/2-center)/real_gridres
            beadi_max_dec = (bead_max+dim/2-center)/real_gridres

            ### Round buffer decimal index to integer index
            beadi_min_int = int(beadi_min_dec)
            beadi_max_int = int(beadi_max_dec)+1 # +1 due to last index not being counted
            
            if beadi_min_int < 0:
                beadi_min_int = 0
            if beadi_max_int < 0:
                beadi_max_int = 0
                
            return beadi_min_int, beadi_max_int

        geom = bbox_polygon

        xlen = xmax - xmin
        ylen = ymax - ymin
        cx = (xmax + xmin)/2
        cy = (ymax + ymin)/2

        ### Checks if any parts of grid is outside exterior bounds (bbox)
        poly_xmin, poly_ymin, poly_xmax, poly_ymax = geom.bounds # xmin, xmax, ymin, ymax

        lcx = np.mean([poly_xmax, poly_xmin])
        lcy = np.mean([poly_ymax, poly_ymin])
        llx = poly_xmax - poly_xmin
        lly = poly_ymax - poly_ymin

        xmin_i, xmax_i = coord_to_indices(lcx, xlen, cx, xreal_gridres, llx/2, llx/2)
        ymin_i, ymax_i = coord_to_indices(lcy, ylen, cy, yreal_gridres, lly/2, lly/2)
        xcoords_in_leaflet = xcoords[xmin_i:xmax_i+1]
        ycoords_in_leaflet = ycoords[ymin_i:ymax_i+1]

        points_to_check = [(round(x, 3), round(y, 3)) for x in xcoords_in_leaflet for y in ycoords_in_leaflet]
        Points_to_check = [shapely.Point(point) for point in points_to_check]
        
        geom_buffered = shapely.buffer(geom, gridres)
        points_contained_bools = geom_buffered.contains(Points_to_check)
        points_in_geom = [point for point, point_bool in zip(points_to_check, points_contained_bools) if not point_bool]
        
        for lcx, lcy in points_in_geom:
            xmin, xmax = coord_to_indices(lcx, xlen, cx, xreal_gridres, gridres, gridres)
            ymin, ymax = coord_to_indices(lcy, ylen, cy, yreal_gridres, gridres, gridres)
            grid_bool_matrix[:, xmin:xmax, ymin:ymax] = 0

        self.print_term("Number of grid points:", [np.sum(grid_bool_matrix[i,:,:]) for i in range(len(radii_sorted_by_size))], "Contained grid  points", spaces=4, verbose=3)

        ### Segmentize bbox to ensure points are close to each other
        segmentized_bbox = bbox_polygon.segmentize(gridres)

        ### Exterior points
        exterior = segmentized_bbox.exterior
        xs, ys = exterior.coords.xy
        exterior_points = list(zip(xs, ys))
        for li, (name, nlipid, radius) in enumerate(lipids_sorted_by_size):
            buffer = radius
            for x, y in exterior_points:
                xmin, xmax = coord_to_indices(x, xlen, cx, xreal_gridres, buffer, buffer)
                ymin, ymax = coord_to_indices(y, ylen, cy, yreal_gridres, buffer, buffer)
                grid_bool_matrix[li, xmin:xmax, ymin:ymax] = 0

        self.print_term("Number of grid points:", [np.sum(grid_bool_matrix[i,:,:]) for i in range(len(radii_sorted_by_size))], "Exterior points", spaces=4, verbose=3)
        
        interiors = segmentized_bbox.interiors
        for interior in interiors:
            xs, ys = interior.coords.xy
            interior_points = list(zip(xs, ys))
            for li, (name, nlipid, radius) in enumerate(lipids_sorted_by_size):
                buffer = radius
                for x, y in interior_points:
                    xmin, xmax = coord_to_indices(x, xlen, cx, xreal_gridres, buffer, buffer)
                    ymin, ymax = coord_to_indices(y, ylen, cy, yreal_gridres, buffer, buffer)
                    grid_bool_matrix[li, xmin:xmax, ymin:ymax] = 0

        self.print_term("Number of grid points:", [np.sum(grid_bool_matrix[i,:,:]) for i in range(len(radii_sorted_by_size))], "Interior points", spaces=4, verbose=3)

        multiplier = 4
        nysplits = int(ylen / (max(radii_sorted_by_size)*multiplier))
        nxsplits = int(xlen / (max(radii_sorted_by_size)*multiplier))
        subgrids_x = np.linspace(cx - xlen/2, cx + xlen/2, nxsplits+1)
        subgrids_y = np.linspace(cy - ylen/2, cy + ylen/2, nysplits+1)

        subgrids = [] # (xmin, xmax, ymin, ymax)
        for xi, x in enumerate(subgrids_x):
            if xi == 0:
                continue
            for yi, y in enumerate(subgrids_y):
                if yi == 0:
                    continue
                subgrid_cx = (subgrids_x[xi-1] + subgrids_x[xi])/2
                subgrid_cy = (subgrids_y[yi-1] + subgrids_y[yi])/2
                subgrid_xlen = (subgrids_x[xi] - subgrids_x[xi-1])/2
                subgrid_ylen = (subgrids_y[yi] - subgrids_y[yi-1])/2
                xmin, xmax = coord_to_indices(subgrid_cx, xlen, cx, xreal_gridres, subgrid_xlen, subgrid_xlen)
                ymin, ymax = coord_to_indices(subgrid_cy, ylen, cy, yreal_gridres, subgrid_ylen, subgrid_ylen)
                subgrids.append((
                    xmin,
                    xmax+1,
                    ymin,
                    ymax+1,
                ))
        random.shuffle(subgrids)

        dict_for_plotting["lipid_insertion"]["subgrid_polygons_points"] = []
        dict_for_plotting["lipid_insertion"]["subgrid_polygons"] = []
        for xmin, xmax, ymin, ymax in subgrids:
            subgrid_polygon_points = [(xmax, ymax), (xmax, ymin), (xmin, ymin), (xmin, ymax)]
            subgrid_Polygon = shapely.Polygon(subgrid_polygon_points)
            dict_for_plotting["lipid_insertion"]["subgrid_polygons_points"].append(subgrid_polygon_points)
            dict_for_plotting["lipid_insertion"]["subgrid_polygons"].append(subgrid_Polygon)

        if self.plot_grid:
            grid_points = []
        step = 1
        subgrid_pointer = 0
        for li, (name, nlipid, radius) in enumerate(lipids_sorted_by_size):
            for i in range(nlipid):
                point_found = False
                if point_found == "break":
                    break
                subgrid_pointer_start = subgrid_pointer
                while point_found is False:
                    if subgrid_pointer == len(subgrids):
                        subgrid_pointer = 0
                    xmin, xmax, ymin, ymax = subgrids[subgrid_pointer]
                    indices_where_true = np.nonzero(grid_bool_matrix[li, xmin:xmax, ymin:ymax])
                    ntrues = len(indices_where_true[0])
                    if ntrues > 0:
                        point_found = True
                    else:
                        subgrid_pointer += 1
                        if subgrid_pointer_start == subgrid_pointer:
                            point_found = "break"
                            break
                        else:
                            continue
                    random_int = np.random.randint(len(indices_where_true[0]),size=1)
                    xi, yi = np.take(indices_where_true, random_int, axis=1)
                    x = xcoords[int(xi)+xmin]
                    y = ycoords[int(yi)+ymin]
                    subgrid_pointer += 1
                if point_found == "break":
                    break
                if self.plot_grid:
                    grid_points.append((x, y))

                for li2, radius2 in enumerate(radii_sorted_by_size):

                    buffer = max([radius, radius2])/2
                    geom = shapely.Point((x,y)).buffer(radius+buffer)
                    # geom = shapely.Point((x,y)).buffer(radius)

                    ### Checks if any parts of grid is outside exterior bounds (bbox)
                    poly_xmin, poly_ymin, poly_xmax, poly_ymax = geom.bounds # xmin, xmax, ymin, ymax

                    lcx = np.mean([poly_xmax, poly_xmin])
                    lcy = np.mean([poly_ymax, poly_ymin])
                    llx = poly_xmax - poly_xmin
                    lly = poly_ymax - poly_ymin

                    xmin_i, xmax_i = coord_to_indices(lcx, xlen, cx, xreal_gridres, llx/2, llx/2)
                    ymin_i, ymax_i = coord_to_indices(lcy, ylen, cy, yreal_gridres, lly/2, lly/2)
                    xcoords_in_leaflet = xcoords[xmin_i:xmax_i+1]
                    ycoords_in_leaflet = ycoords[ymin_i:ymax_i+1]

                    points_to_check = [(round(x, 5), round(y, 5)) for x in xcoords_in_leaflet for y in ycoords_in_leaflet]
                    Points_to_check = [shapely.Point(point) for point in points_to_check]
                    
                    points_contained_bools = geom.contains(Points_to_check)
                    points_in_geom = [point for point, point_bool in zip(points_to_check, points_contained_bools) if point_bool]
                    
                    for lcx, lcy in points_in_geom:
                        xmin, xmax = coord_to_indices(lcx, xlen, cx, xreal_gridres, gridres, gridres)
                        ymin, ymax = coord_to_indices(lcy, ylen, cy, yreal_gridres, gridres, gridres)
                        grid_bool_matrix[li2, xmin:xmax, ymin:ymax] = 0
                
                self.print_term("Number of grid points:", [np.sum(grid_bool_matrix[i,:,:]) for i in range(len(radii_sorted_by_size))], "Lipid insertion step nr: " + str(step) + " / " + str(total_number_of_lipids), spaces=4, verbose=3)
                if self.plot_grid:
                    dict_for_plotting["lipid_insertion"]["lipid_insertion_steps"].append({
                        "grid_points": grid_points.copy(),
                        "grid_bool_matrix": grid_bool_matrix.copy(),
                    })
                step += 1
        
        return grid_points, lipids, dict_for_plotting