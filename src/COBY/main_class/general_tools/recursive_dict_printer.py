
class recursive_dict_printer():
    #####################################
    ### Specialized printing function ###
    #####################################
    def recursive_dict_printer(self, iterable, current_spacing="", spacing="    "):
        '''
        Prints dictionaries recursively
        '''
        if type(iterable) == dict:
            for key, val in iterable.items():
                if type(val) == dict:
                    self.print_term(current_spacing + str(key) + ": {", verbose=0)
                    self.recursive_dict_printer(iterable=val, current_spacing=current_spacing+spacing, spacing=spacing)
                    self.print_term(current_spacing + "}", verbose=0)
                else:
                    self.print_term(current_spacing + str(key) + ": " + str(val), verbose=0)
        else:
            self.print_term(current_spacing + str(iterable), verbose=0)

