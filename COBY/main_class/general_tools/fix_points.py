class fix_points:
    def fix_points(self, grid, leaf_cmd):
        #################
        ### CENTERING ###
        #################
        for grid_point_nr in grid.keys():
            grid[grid_point_nr]["x"] -= 0.5 * leaf_cmd["x"]
            grid[grid_point_nr]["y"] -= 0.5 * leaf_cmd["y"]
            grid[grid_point_nr]["z"] = 0
            grid[grid_point_nr]["internal_z"] = 0

        #################
        ### ROTATIONS ###
        #################
        ### No rotations around x/y axis as of right now
#         x_ang_deg, y_ang_deg, z_ang_deg = leaf_cmd["rx"], leaf_cmd["ry"], leaf_cmd["rz"]
#         x_ang_deg, y_ang_deg, z_ang_deg = 0, 0, leaf_cmd["rz"]
#         if any([ang != 0 for ang in [x_ang_deg, y_ang_deg, z_ang_deg]]):
#             for grid_point_nr in grid.keys():
#                 grid[grid_point_nr]["x"], grid[grid_point_nr]["y"], grid[grid_point_nr]["z"] = rotate_point(
#                     grid[grid_point_nr]["x"],
#                     grid[grid_point_nr]["y"],
#                     grid[grid_point_nr]["z"],
#                     x_ang_deg,
#                     y_ang_deg,
#                     z_ang_deg,
#                 )
#                 grid[grid_point_nr]["x_ang_deg"], grid[grid_point_nr]["y_ang_deg"], grid[grid_point_nr]["z_ang_deg"] = x_ang_deg, y_ang_deg, z_ang_deg

        ####################
        ### TRANSLATIONS ###
        ####################
        cx, cy, cz = leaf_cmd["center"]
        if any([ang != 0 for ang in [cx, cy, cz]]):
            for grid_point_nr in grid.keys():
                grid[grid_point_nr]["x"] += cx
                grid[grid_point_nr]["y"] += cy
                grid[grid_point_nr]["z"] += cz

        #################################################
        ### CHECKS IF COORDINATES ARE OUTSIDE THE BOX ###
        #################################################
        if leaf_cmd["pbc_check"]:
            errors_count = 0
            for grid_point_nr in grid.keys():
                bead_coords = [grid[grid_point_nr]["x"], grid[grid_point_nr]["y"], grid[grid_point_nr]["z"]]
                checked_beads, error = self.coord_checker(bead_coords, self.pbc_box, error_count = True)
                if error > 0:
                    errors_count += 1
                    grid[grid_point_nr]["x"] = bead_coords[0]
                    grid[grid_point_nr]["y"] = bead_coords[1]
                    grid[grid_point_nr]["z"] = bead_coords[2]

        return grid

