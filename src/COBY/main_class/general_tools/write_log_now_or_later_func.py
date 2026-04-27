class write_log_now_or_later_func:
    def write_log_now_or_later_func(self, *strings, spaces=1, verbose=1, warn=False, debug=False, debug_keys=[], write_log_now_or_later="now"):
        if write_log_now_or_later == "now":
            self.print_term(*strings, spaces=spaces, verbose=verbose, warn=warn, debug=debug, debug_keys=debug_keys)
        elif write_log_now_or_later == "later":
            self.print_term(*strings, spaces=spaces, verbose=verbose, warn=warn, debug=debug, debug_keys=debug_keys, log_string=False)
            ### First print then use LOG_FILE to append to log file laterlater
            strings_joined = " ".join([str(s) for s in strings])
            self.LOG_FILE.append({"string": strings_joined, "for_print": False, "for_log": True, "spaces": spaces, "verbose": verbose, "debug": debug, "debug_keys": debug_keys})
