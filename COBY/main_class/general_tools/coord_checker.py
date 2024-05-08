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

            if len(positions) == len(axes):
                for pi, (pos, ax) in enumerate(zip(positions, axes)):
                    if pos > ax / 2:
                        positions[pi] -= ax
                        if error_count:
                            count += 1
                        else:
                            self.print_term("    ", "    ", "WARNING: Points are outside pbc. Moved to other side. Expect potential overlap of lipids", warn = True)
                        
                    elif pos < -ax / 2:
                        positions[pi] += ax
                        if error_count:
                            count += 1
                        else:
                            self.print_term("    ", "    ", "WARNING: Points are outside pbc. Moved to other side. Expect potential overlap of lipids", warn = True)
            else:
                self.print_term("    ", "    ", "Incorrect dimensions for 'positions' and 'axes' in in 'coord_checker' function")
                self.print_term("    ", "    ", "'positions' length:", len(positions), "'axes' length:", len(axes))
            
            if error_count:
                return positions, count
            else:
                return positions
