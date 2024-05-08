class print_term:
    #####################################
    ### Specialized printing function ###
    #####################################
    def print_term(self, *string, spaces=0, verbose=0, space="    ", sep=" ", end="\n", debug = False, extra = False, warn = False):
        '''
        Specialized printing function.
        Allows for easy customization of requirements for when a print statement should be printed or added to log file.
        '''
        
        print_true = False
        if type(string) in [tuple, list]:
            string = space*spaces + sep.join([str(i) for i in string])

        ### Check if string should be printed based on settings
        if debug:
            if self.debug_prints:
                print_true = True
                if not string.startswith("DEBUG:"):
                    string = "DEBUG: " + string
        elif extra:
            if self.extra_info:
                print_true = True
        elif warn:
            if self.warnings:
                print_true = True
                if not string.startswith("WARNING:"):
                    string = "WARNING: " + string
        else:
            ### If no special case, then assume it should be printed
            print_true = True

        ### Print to terminal if not quiet
        if print_true and not self.quiet and self.verbose >= verbose:
            print(string, end=end)
        
        ### Appends string to log file, for later writing
        if print_true and self.output_log_file_name:
            self.LOG_FILE.append(string + end)
    
