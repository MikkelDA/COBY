import time
import math
import numpy as np
import scipy
import shapely

class plane_grid_point_optimizer:
    def plane_grid_point_optimizer(self, grid_points, lipid_sizes, polygon, xcenter, ycenter, xlen, ylen, maxsteps, push_tolerance, push_mult, buffer, occupation_modifier, optimize, optimizer_print_to_term):

        def push_func(dist, power, mult):
            return dist**-power*mult

        def get_vector(A, B):
            ### Get vector from point A to point B
            return (round(B[0] - A[0], 3), round(B[1] - A[1], 3))
        
        def get_vector_len_fastsingle(v):
            '''
            Fastest single vector calculator
            https://stackoverflow.com/questions/37794849/efficient-and-precise-calculation-of-the-euclidean-distance
            https://stackoverflow.com/questions/7370801/how-do-i-measure-elapsed-time-in-python
            '''
            return math.sqrt(v[0]**2 + v[1]**2)

        def rand_sign():
            if random.random() < 0.5:
                return 1
            else:
                return -1

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
#                             if dist < bin_size*3:
                            if dist < largest_lipid*4:
                                neighborlist[pi1, pi2] = 1
            nst_toc = time.time()
            if self.debug_prints:
                self.print_term("Time spent calculating neighborlist:", round(nst_toc-nst_tic, 4), spaces=4, debug=True)

            return neighborlist

        def coord_to_indices(pos, dim, center, real_gridres):
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

#             indeces = [-2, -1, 0, 1, 2]
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
                self.print_term("Time spent calculating bins:        ", round(bin_toc-bin_tic, 4), spaces=4, debug=True)

            return new_bins_arr
        
        points_arr = np.array(grid_points)
        
        largest_lipid  = max(lipid_sizes)
        bin_size       = largest_lipid*2
        
        if self.debug_prints:
            self.print_term("len(grid_points)", len(grid_points), debug=True)
            self.print_term("xlen            ", xlen,             debug=True)
            self.print_term("ylen            ", ylen,             debug=True)
            self.print_term("bin_size        ", bin_size,         debug=True)
        
        ### Calculating number of bins along each axis, by rounding down
        xnbins = int(xlen / bin_size)
        ynbins = int(ylen / bin_size)
        ### Finding the actual lengths of the bins
        xbinlen = xlen / xnbins
        ybinlen = ylen / ynbins
        
        dists_traveled = np.zeros((points_arr.shape[0],))
        
        if self.plot_grid:
            POINT_STEPS = []
            POINT_STEPS.append(points_arr.copy())
        else:
            POINT_STEPS = False

        steps_tic = time.time()
        
        ### Making initial bins and neighborlists
        if self.debug_prints:
            self.print_term("STEP:", 0, debug=True, spaces=2)
            self.print_term("Updating bins and neigborlist", debug=True, spaces=3)
        
        dists_traveled = np.zeros((points_arr.shape[0],))
        bins_arr       = binning_func()
        neighborlist   = update_neighborlist(bins_arr)
        
        step_modifier_limit = 15
        
        bounce_counter = np.zeros((len(points_arr)))
        max_push = 0

        if optimizer_print_to_term and not self.debug_prints:
            self.print_term("CURRENT STEP:", end=" ", spaces=3, verbose=2)
        
        for si, step in enumerate(np.arange(1, maxsteps+1)):
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
            
            if step >= step_modifier_limit:
                step_modifier = step_modifier_limit**-1
            else:
                step_modifier = 1-(step/step_modifier_limit)
            
            if optimizer_print_to_term and not self.debug_prints:
                self.print_term(step, end=" ", spaces=0, verbose=2)
            
            push_arr = np.empty((len(points_arr)), dtype=object)
            for i in range(len(push_arr)):
                push_arr[i] = []

            pushes   = np.zeros((len(points_arr)))
            max_dists_traveled = max(dists_traveled)
            
            if self.debug_prints:
                self.print_term("STEP:", step, "AT TIME:", round(time.time() - steps_tic, 4), debug=True, spaces=2)
                self.print_term("Max push:         ", round(max_push, 4),                     debug=True, spaces=3)
                self.print_term("Max dist traveled:", round(max_dists_traveled, 4),           debug=True, spaces=3)
            
            if max_dists_traveled >= bin_size/4:
                if self.debug_prints:
                    self.print_term("Updating bins and neigborlist", debug=True, spaces=3)
                dists_traveled = np.zeros((points_arr.shape[0],))
                bins_arr       = binning_func()
                neighborlist   = update_neighborlist(bins_arr)
            
            ### Lipid-lipid pushes
            for pi1, point1 in enumerate(points_arr):
                neighbor_is = np.nonzero(neighborlist[pi1])[0]
                for pi2 in neighbor_is:
                    if pi2 <= pi1:
                        continue
                    
                    point2 = points_arr[pi2]

                    vector              = np.array(get_vector(point1, point2))
                    vector_len          = get_vector_len_fastsingle(vector)
                    dist                = vector_len
                    combined_lipid_size = lipid_sizes[pi1] + lipid_sizes[pi2]
                    
                    ### Dev note: Other optimization algorithm versions removed in COBY version 0.0.8. Look in 0.0.7 for the rest.                    
                    ideal_dist = combined_lipid_size*(1+occupation_modifier)
                    ideal_dist_diff = ideal_dist - combined_lipid_size
                    ideal_dist_lower = combined_lipid_size
                    ideal_dist_upper = ideal_dist + ideal_dist_diff
                    if dist <= ideal_dist_upper:
                        if not (vector[0] == 0. and vector[1] == 0.):
                            vector /= vector_len                                
                        push = 0
                        
                        ### Difference in distance from ideal distance 
                        if dist <= ideal_dist_lower: # combined_lipid_size
                            ### Push very hard if way too close
                            push += (combined_lipid_size-dist)/2

                        ### Always add smaller extra push
                        push += abs(ideal_dist_upper-dist)*step_modifier
                        
                        vector_push = vector*push

                        ### Pseudo-normalizes the force behind the pushes
                        max_size = max([lipid_sizes[pi2], lipid_sizes[pi1]])
                        pi1_mult = 1/(max_size/lipid_sizes[pi2])
                        pi2_mult = 1/(max_size/lipid_sizes[pi1])
                        push_arr[pi1].append(-vector_push*pi1_mult)
                        push_arr[pi2].append(+vector_push*pi2_mult)

            ### Pushes grid points
            for pi, push in enumerate(push_arr):
                if push:
                    push_vector         = np.mean(np.array(push), axis=0)
                    push_vector_len     = get_vector_len_fastsingle(push_vector)
                    points_arr[pi]     += np.mean(np.array(push), axis=0)
                    pushes[pi]         += push_vector_len
                    dists_traveled[pi] += push_vector_len

            ### BBOX and protein distance checks
            all_contained = True
            for pi1, point1 in enumerate(points_arr):
                if bounce_counter[pi1] > 0:
                    bounce_counter[pi1] -= 0.1
                ### Shapely Point
                point_Point = shapely.Point(point1)
                ### polygon.boundary includes both the outer surface and the surface of holes
                dist = polygon.boundary.distance(point_Point)
                point_contained = polygon.contains(point_Point)

                eq_dist = lipid_sizes[pi1]*(1+occupation_modifier*2)
                
                ### Lipid contained in legal areas but close to edge
                if point_contained and dist < eq_dist:
                    poly_nearest, point_Point = shapely.ops.nearest_points(polygon.boundary, point_Point)
                    vector = np.array(get_vector((poly_nearest.x, poly_nearest.y), (point_Point.x, point_Point.y)))
                    
                    diff        = eq_dist - dist
                    push        = diff/dist
                    
                    if dist < lipid_sizes[pi1]:
                        bounce_counter[pi1] += 1
                        push *= bounce_counter[pi1]

                    vector_push = vector*push
                    vector_push_len = get_vector_len_fastsingle(vector_push)
                    
                    points_arr[pi1]     += vector_push
                    pushes[pi1]         += vector_push_len
                    dists_traveled[pi1] += vector_push_len
                
                ### Lipid outside legal areas
                ### Push it just enough for it to be contained in legal areas
                elif not point_contained:
                    all_contained = False
                    poly_nearest, point_Point = shapely.ops.nearest_points(polygon.boundary, point_Point)
                    vector = np.array(get_vector((poly_nearest.x, poly_nearest.y), (point_Point.x, point_Point.y)))
                    
                    diff        = -(eq_dist + dist)
                    push        = diff/dist
                    vector_push = vector*push
                    vector_push_len = get_vector_len_fastsingle(vector_push)
                    
                    points_arr[pi1]     += vector_push
                    pushes[pi1]         += vector_push_len
                    dists_traveled[pi1] += vector_push_len
            
            if self.plot_grid:
                POINT_STEPS.append(points_arr.copy())

            ### Checks if anything was pushed
            ### Values are floats so can't use '=='
            max_push = np.max(pushes)
            if max_push < push_tolerance and all_contained:
                break

        steps_toc = time.time()
        steps_time = steps_toc - steps_tic
        mean_steps_time = steps_time/(step)
        if optimizer_print_to_term:
            self.print_term("", verbose=2)
            self.print_term("Last step:         ", step,                            spaces=3, verbose=2)
            self.print_term("Optimization time: ", round(steps_time, 4),     "[s]", spaces=3, verbose=2)
            self.print_term("Mean step time:    ", round(steps_time/(step), 4), "[s]", spaces=3, verbose=2)
        
        return points_arr, POINT_STEPS, step, steps_time, mean_steps_time, max_push, bin_size
