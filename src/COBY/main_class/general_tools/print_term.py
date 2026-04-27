class print_term:
    #####################################
    ### Specialized printing function ###
    #####################################
    def print_term(self, *string, spaces=0, verbose=0, space="    ", sep=" ", end="\n", debug = False, debug_keys=[], extra = False, warn = False, inp=False, print_string=True, log_string=True):
        '''
        Specialized printing function.
        Allows for easy customization of requirements for when a print statement should be printed or added to the log file.
        '''
        
        valid = False
        if type(string) in [tuple, list]:
            string = space*spaces + sep.join([str(i) for i in string])

        ### Check if string should be printed based on settings
        if debug:
            ### If debug prints defined and either no debug keys given or specific debug key in given debug keys
            if self.debug_prints and ((not self.debug_keys) or (self.debug_keys and debug_keys)):
                if type(debug_keys) == str:
                    debug_keys = [debug_keys]
                if not self.debug_keys or any([key in self.debug_keys for key in debug_keys]):
                    valid = True
                    if not string.startswith("DEBUG:"):
                        string = "DEBUG: " + string
        elif warn:
            if self.warnings:
                valid = True
                if not string.startswith("WARNING:"):
                    string = "WARNING: " + string
        elif extra:
            if self.extra_info:
                valid = True
        elif print_string or log_string:
            valid = True

        if self.debug_prints and ((not self.debug_keys) or ("printer" in self.debug_keys)):
            print("-------------------------------------------------------------")
            print("string:",       "'" + str(string) + "'")
            print("spaces:",       "'" + str(spaces) + "'")
            print("verbose:",      "'" + str(verbose) + "'")
            print("space:",        "'" + str(space) + "'")
            print("sep:",          "'" + str(sep) + "'")
            print("end:",          "'" + str(end) + "'")
            print("debug:",        "'" + str(debug) + "'")
            print("debug_keys:",   "'" + str(debug_keys) + "'")
            print("extra:",        "'" + str(extra) + "'")
            print("warn:",         "'" + str(warn) + "'")
            print("inp:",          "'" + str(inp) + "'")
            print("print_string:", "'" + str(print_string) + "'")
            print("log_string:",   "'" + str(log_string) + "'")
            print("valid:",        "'" + str(valid) + "'")
            print("-------------------------------------------------------------")

        ### Print to terminal if not quiet
        if print_string and valid and self.verbose >= verbose and not self.quiet:
            print(string, end=end, flush=True)
        
        ### Writes to log file
        if log_string and valid and self.LOG_FILE_HANDLE:
            self.LOG_FILE_HANDLE.write(string + end)
            self.LOG_FILE_HANDLE.flush()
        
        if inp:
            inp_answer = input("Your response: ")
            if self.LOG_FILE_HANDLE:
                self.LOG_FILE_HANDLE.write("Your response: " + inp_answer + "\n")
            return inp_answer
