import time
import math
import numpy as np
import scipy
import shapely

class plane_grid_point_optimizer:
    def plane_grid_point_optimizer(self, grid_points, lipid_sizes, polygon, xcenter, ycenter, xlen, ylen, max_steps, push_tolerance, lipid_push_multiplier, edge_push_multiplier, occupation_modifier, optimizer_print_to_term):

        def get_vector(A, B):
            '''
            Get vector from point A to point B
            '''
            return (round(B[0] - A[0], 3), round(B[1] - A[1], 3))
        
        def get_vector_len_fastsingle(v):
            '''
            Fastest single vector calculator according to the following
            https://stackoverflow.com/questions/37794849/efficient-and-precise-calculation-of-the-euclidean-distance
            https://stackoverflow.com/questions/7370801/how-do-i-measure-elapsed-time-in-python
            '''
            return math.sqrt(v[0]**2 + v[1]**2)

        def update_neighborlist(bins_arr):
            neighborlist = np.zeros((points_arr.shape[0], points_arr.shape[0]), dtype=int)
            nst_tic = time.time()
            for bix, binx in enumerate(bins_arr):
                for biy, biny in enumerate(binx):
                    points_for_nstlist = biny
                    for pj1, pi1 in enumerate(points_for_nstlist):
                        for pj2, pi2 in enumerate(points_for_nstlist[pj1+1:], pj1+1):
                            point1 = points_arr[pi1]
                            point2 = points_arr[pi2]
                            dist = scipy.spatial.distance.cdist([point1], [point2])[0]
                            dist = get_vector_len_fastsingle(get_vector(point1, point2))
                            if dist < bin_size*2: ### Equivalent to largest_lipid*4
                                neighborlist[pi1, pi2] = 1
            nst_toc = time.time()
            if self.debug_prints:
                self.print_term("Time spent calculating neighborlist:", round(nst_toc-nst_tic, 4), spaces=5, debug=True, debug_keys="optimizer")

            return neighborlist

        def coord_to_indices(pos, dim, center, real_gridres):
            '''
            Converts a coordinate to an index value
            '''
            posi_dec = (pos+dim/2-center)/real_gridres
            posi_int = int(posi_dec)
            if posi_int < 0:
                posi_int = 0
            return posi_int

        def binning_func():
            ### Creating bins array
            new_bins_arr = np.empty((xnbins, ynbins), dtype=object)
            ### np.fill([]) makes all array indeces point to the same list so can't be used
            ### For loop is very fast so doesn't matter for time
            bin_tic = time.time()
            for xi in range(xnbins):
                for yi in range(ynbins):
                    new_bins_arr[xi][yi] = []

            indeces = [-1, 0, 1]
            for pi, (px, py) in enumerate(points_arr):
                pix = coord_to_indices(px, xlen, xcenter, xbinlen)
                piy = coord_to_indices(py, ylen, ycenter, ybinlen)
                ### Surrounding bins
                for xi in indeces:
                    for yi in indeces:
                        if xnbins > pix+xi >= 0 and ynbins > piy+yi >= 0:
                            new_bins_arr[pix+xi][piy+yi].append(pi)
            bin_toc = time.time()
            if self.debug_prints:
                self.print_term("Time spent calculating bins:        ", round(bin_toc-bin_tic, 4), spaces=5, debug=True, debug_keys="optimizer")

            return new_bins_arr
        
        points_arr = np.array(grid_points)
        
        largest_lipid  = max(lipid_sizes)
        bin_size       = largest_lipid*2
        
        ### Calculating number of bins along each axis, by rounding down
        xnbins = int(xlen / bin_size)
        ynbins = int(ylen / bin_size)
        ### Finding the actual lengths of the bins
        xbinlen = xlen / xnbins
        ybinlen = ylen / ynbins
        
        if self.debug_prints:
            self.print_term("Optimizer info:",                    debug=True, debug_keys="optimizer", spaces=3)
            self.print_term("len(grid_points)", len(grid_points), debug=True, debug_keys="optimizer", spaces=4)
            self.print_term("xlen            ", xlen,             debug=True, debug_keys="optimizer", spaces=4)
            self.print_term("ylen            ", ylen,             debug=True, debug_keys="optimizer", spaces=4)
            self.print_term("bin_size        ", bin_size,         debug=True, debug_keys="optimizer", spaces=4)
            self.print_term("xnbins          ", xnbins,           debug=True, debug_keys="optimizer", spaces=4)
            self.print_term("ynbins          ", ynbins,           debug=True, debug_keys="optimizer", spaces=4)
            self.print_term("xbinlen         ", xbinlen,          debug=True, debug_keys="optimizer", spaces=4)
            self.print_term("ybinlen         ", ybinlen,          debug=True, debug_keys="optimizer", spaces=4)
        
        if self.plot_grid:
            POINT_STEPS = []
            POINT_STEPS.append(points_arr.copy())
        else:
            POINT_STEPS = False

        steps_tic = time.time()
        
        ### Making initial bins and neighborlists
        if self.debug_prints:
            self.print_term("STEP:", 0,                      debug=True, debug_keys="optimizer", spaces=3)
            self.print_term("Updating bins and neigborlist", debug=True, debug_keys="optimizer", spaces=4)
        
        unmoved_points     = np.array(grid_points)
        bins_arr           = binning_func()
        neighborlist       = update_neighborlist(bins_arr)
        
        step_multiplier_limit = 15
        
        bounce_counter = np.zeros((len(points_arr)))
        max_push = 0

        if optimizer_print_to_term and not self.debug_prints:
            self.print_term("CURRENT STEP:", end=" ", spaces=3, verbose=2)
        
        for si, step in enumerate(np.arange(1, max_steps+1)):
            '''CODE LEGEND

            'dists_traveled':
                Numpy array with shape (nlipids).
                Remembers how much all individual lipids have moved.
                Forces a bin and neighborlist update if largest movement is larger than the bin size.
                Resets to zero-array upon neighborlist updates.

            'bins_arr':
                Numpy array with shape (xbins, ybins).
                Each index is a list of points that are present within the bin or in the neighbouring bins.
                Used for neighborlist updates.

            'bin_pointers': Currently unused
                List with length (nlipids) with each index being a tuple containing (xpointer, ypointer).
                Contains the bin x/y indeces that points are in.

            'neighborlist':
                Numpy array with shape (nlipids, nlipids).
                Contains the neighbor indeces for each grid point.

            'points_arr':
                Numpy array with shape (nlipids, 2).
                Contains the x/y coordinates of all grid points.

            'push_arr':
                Numpy array with shape (nlipids, 2).
                Contains the combined vector push applied during a single time step.
                Modifies 'points_arr' at the end of a time step and then resets to only containing zeros.

            'pushes':
                Numpy array with shape (nlipids).
            '''
            step_tic = time.time()
            
            if step >= step_multiplier_limit:
                step_multiplier = step_multiplier_limit**-1
            else:
                step_multiplier = 1-((step-1)/step_multiplier_limit) # 'step-1' to have first step modifier be '1'
            
            if optimizer_print_to_term and not self.debug_prints:
                self.print_term(step, end=" ", spaces=0, verbose=2)
            
            push_arr = np.empty((len(points_arr)), dtype=object)
            for i in range(len(push_arr)):
                push_arr[i] = []

            pushes = np.zeros((len(points_arr)))

            ### Lipid-lipid pushes
            for pi1, point1 in enumerate(points_arr):
                neighbor_is = np.nonzero(neighborlist[pi1])[0]
                for pi2 in neighbor_is:
                    ### Skipping to avoid double-processing lipids
                    if pi2 <= pi1:
                        continue
                    
                    point2 = points_arr[pi2]

                    vector              = np.array(get_vector(point1, point2))
                    dist                = get_vector_len_fastsingle(vector)
                    combined_lipid_size = lipid_sizes[pi1] + lipid_sizes[pi2]
                    
                    dist_ideal       = combined_lipid_size * (1 + occupation_modifier)
                    dist_ideal_diff  = dist_ideal - combined_lipid_size
                    dist_ideal_lower_limit = combined_lipid_size
                    dist_ideal_upper_limit = dist_ideal + dist_ideal_diff

                    ### Final equation for 'dist_ideal_upper_limit' becomes:
                    ### dist_ideal_upper_limit = combined_lipid_size * (2 * occupation_modifier + 1)

                    ### No point in pushing if they are already far enough away. Also saves some computing power.
                    if dist <= dist_ideal_upper_limit:
                        ### Could cause problems if not included. I forgot to write down what those problems were.
                        ### Update in 0.2.7: The problem is probably just division by zero (happens if two points are on top of each other)
                        ### ### Changed it such that the vector then points in a random direction instead of simply having a length of zero 
                        if not (vector[0] == 0. and vector[1] == 0.):
                            vector /= dist
                        else:
                            rand_vector = np.random.rand(2)
                            vector = rand_vector / get_vector_len_fastsingle(rand_vector)
                        
                        push = 0
                        
                        ### Difference in distance from ideal distance. Forces the two lipids to move out of each others radii.
                        if dist <= dist_ideal_lower_limit: # combined_lipid_size
                            ### Push very hard if way too close
                            push += (combined_lipid_size-dist)/2

                        ### Always add smaller extra push to move lipids closer to dist_ideal
                        ### Step modifier reduces the push for each step of the optimization algorithm to "calm it down"
                        ### Removed "abs" from "abs(dist_ideal_upper_limit-dist)" in v0.2.7 since it seems superfluous since "dist <= dist_ideal_upper" requires that "dist" always be smaller or equal to "dist_ideal_upper"
                        push += (dist_ideal_upper_limit-dist)*step_multiplier
                        
                        vector_push = vector*push*lipid_push_multiplier

                        ### Pseudo-normalizes the force behind the pushes
                        ### The largest lipid is pushed with the full force, but the force on the smallest lipid is decreased
                        max_size = max([lipid_sizes[pi2], lipid_sizes[pi1]])
                        pi1_mult = 1/(max_size/lipid_sizes[pi2]) ### Equivalent to lipid_sizes[pi2] / max_size
                        pi2_mult = 1/(max_size/lipid_sizes[pi1]) ### Equivalent to lipid_sizes[pi1] / max_size
                        push_arr[pi1].append(-vector_push*pi1_mult)
                        push_arr[pi2].append(+vector_push*pi2_mult)

            ### Pushes grid points
            for pi, push in enumerate(push_arr):
                if push:
                    push_vector         = np.mean(np.array(push), axis=0)
                    push_vector_len     = get_vector_len_fastsingle(push_vector)
                    points_arr[pi]     += np.mean(np.array(push), axis=0)
                    pushes[pi]         += push_vector_len

            ### BBOX and protein distance checks
            all_contained = True
            for pi1, point1 in enumerate(points_arr):
                ### Slightly reduce bounce counter at every step
                if bounce_counter[pi1] > 0:
                    bounce_counter[pi1] -= 0.1

                ### Shapely Point
                point_Point = shapely.Point(point1)
                ### polygon.boundary includes both the outer surface and the surface of holes
                point_contained = polygon.contains(point_Point)
                
                dist    = polygon.boundary.distance(point_Point)
                eq_dist = lipid_sizes[pi1]*(1+occupation_modifier*2)

                ### Lipid contained in legal areas but close to edge
                if point_contained and dist < eq_dist:
                    ### Finds nearest surface point (membrane/protein edge)
                    poly_nearest, point_Point = shapely.ops.nearest_points(polygon.boundary, point_Point)
                    vector                    = np.array(get_vector((poly_nearest.x, poly_nearest.y), (point_Point.x, point_Point.y)))
                    
                    diff = eq_dist - dist
                    push = diff/dist
                    
                    ### If the radius overlaps with illegal areas then push it harder
                    ### Push becomes harder and harder every time the radius of a lipid overlaps with illegal areas
                    if dist < lipid_sizes[pi1]:
                        bounce_counter[pi1] += 1
                        push                *= bounce_counter[pi1]

                    vector_push     = vector*push*edge_push_multiplier
                    vector_push_len = get_vector_len_fastsingle(vector_push)
                    
                    points_arr[pi1] += vector_push
                    pushes[pi1]     += vector_push_len
                
                ### Lipid outside legal areas
                ### Push it just enough for it to be contained in legal areas
                elif not point_contained:
                    all_contained = False
                    poly_nearest, point_Point = shapely.ops.nearest_points(polygon.boundary, point_Point)
                    vector                    = np.array(get_vector((poly_nearest.x, poly_nearest.y), (point_Point.x, point_Point.y)))
                    
                    diff            = -(eq_dist + dist)
                    push            = diff/dist
                    vector_push     = vector*push*edge_push_multiplier
                    vector_push_len = get_vector_len_fastsingle(vector_push)
                    
                    points_arr[pi1] += vector_push
                    pushes[pi1]     += vector_push_len
            
            if self.plot_grid:
                POINT_STEPS.append(points_arr.copy())

            ### Checks if anything was pushed
            ### Values are floats so can't use '=='
            max_push = np.max(pushes)
            if max_push < push_tolerance and all_contained:
                break

            step_toc = time.time()

            max_dists_traveled = max([get_vector_len_fastsingle(get_vector(unmoved_point, point)) for unmoved_point, point in zip(unmoved_points, points_arr)])
            
            if self.debug_prints:
                self.print_term("STEP:", step,                                                       debug=True, debug_keys="optimizer", spaces=3)
                self.print_term("Step time:        ",      round(step_toc - step_tic, 4),  "[s]",    debug=True, debug_keys="optimizer", spaces=4)
                self.print_term("Time spent so far:",      round(step_toc - steps_tic, 4), "[s]",    debug=True, debug_keys="optimizer", spaces=4)
                self.print_term("Max push:         ",      round(max_push, 4),             "[nm]",   debug=True, debug_keys="optimizer", spaces=4)
                self.print_term("Max dist traveled:",      round(max_dists_traveled, 4),   "[nm]",   debug=True, debug_keys="optimizer", spaces=4)
                self.print_term("step_multiplier:  ",      round(step_multiplier, 4),      "[mult]", debug=True, debug_keys="optimizer", spaces=4)

            if max_dists_traveled >= bin_size/2:
                if self.debug_prints:
                    self.print_term("Updating bins and neighborlist", debug=True, debug_keys="optimizer", spaces=4)
                unmoved_points = np.array(grid_points)
                bins_arr       = binning_func()
                neighborlist   = update_neighborlist(bins_arr)

        steps_toc = time.time()
        steps_time = steps_toc - steps_tic
        mean_steps_time = steps_time/(step)
        if optimizer_print_to_term:
            self.print_term("", verbose=2)
            self.print_term("Last step:         ", step,                               spaces=3, verbose=2)
            self.print_term("Optimization time: ", round(steps_time, 4),     "[s]",    spaces=3, verbose=2)
            self.print_term("Mean step time:    ", round(steps_time/(step), 4), "[s]", spaces=3, verbose=2)
        
        return points_arr, POINT_STEPS, step, steps_time, mean_steps_time, max_push, bin_size
