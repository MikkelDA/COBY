class coord_checker:
    def coord_checker(self, positions, axes, error_count = True):
            '''
            Checks if coordinates are within the box and moves them inside if they are not
            '''
            count = 0
            if type(positions) != list:
                positions = [positions]
            if type(axes) != list:
                axes = [axes]
            
            assert len(positions) == len(axes), "\n".join([
                 "Incorrect dimensions for 'positions' and 'axes' in in 'coord_checker' function.",
                 "'positions' length:", len(positions), "'axes' length:", len(axes)
            ])
            new_positions = []
            changes_detected = False
            for pi, (pos, ax) in enumerate(zip(positions, axes)):
                if pos > ax / 2:
                    changes_detected = True
                    new_positions.append(pos - ax)
                    if error_count:
                        count += 1
                    else:
                        self.print_term("    ", "    ", "WARNING: Points are outside pbc. Moved to other side. Expect potential bead overlaps with other system components.", warn = True)
                    
                elif pos < -ax / 2:
                    changes_detected = True
                    new_positions.append(pos + ax)
                    if error_count:
                        count += 1
                    else:
                        self.print_term("    ", "    ", "WARNING: Points are outside pbc. Moved to other side. Expect potential bead overlaps with other system components.", warn = True)
                else:
                    new_positions.append(pos)
            
            ### Recursively checks in case beads are placed so far outside the box that the correction also ends up being outside the box. Should eventually converge towards being inside the box.
            if changes_detected:
                new_positions = self.coord_checker(new_positions, axes, error_count=False)

            if error_count:
                return new_positions, count
            else:
                return new_positions
