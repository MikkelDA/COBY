import ast
import os

class itp_read_initiater:
    def itp_read_initiater(self):
        if len(self.ITP_INPUT_cmds) != 0:
            self.print_term("Loading topology file(s)", spaces=0, verbose=2)
            for itp_i, itp_cmd in enumerate(self.ITP_INPUT_cmds, 0):
                write_includes = True
                itp_input_files = []

                ### First gather all includes and files along with settings
                for cmd in itp_cmd.split():
                    sub_cmd = cmd.split(":")

                    if sub_cmd[0].lower() == "write_includes":
                        assert sub_cmd[1] in ["False", "True", "0", "1"], "Value given to 'write_includes' must be True/1 or False/0: " + sub_cmd[1]
                        if type(sub_cmd[1]) == str:
                            sub_cmd[1] = ast.literal_eval(sub_cmd[1])
                        write_includes = sub_cmd[1]
                    
                    elif sub_cmd[0].lower() == "include":
                        assert os.path.exists(sub_cmd[1]) and os.path.isfile(sub_cmd[1]), "'" + sub_cmd[1] + "' does not exist or is not a file for subargument '" + cmd + "'"
                        itp_input_files.append(("include", sub_cmd[1]))

                    elif sub_cmd[0].lower() == "file":
                        assert os.path.exists(sub_cmd[1]) and os.path.isfile(sub_cmd[1]), "'" + sub_cmd[1] + "' does not exist or is not a file for subargument '" + cmd + "'"
                        itp_input_files.append(("file", sub_cmd[1]))

                    elif sub_cmd[0].lower() == "define":
                        self.itp_defs_all_defnames.add(sub_cmd[1])

                    else:
                        assert False, "Unknown subcommand given to 'itp_input': " + str(cmd)
                
                ### Then run the includes and files in the order that they were originally given
                for filetype, file in itp_input_files:
                    if filetype == "include":
                        self.itp_reader(file, recursion_layer = 1, write_includes = write_includes)
                        if write_includes:
                            self.TOP_include_statements.append("#include \"" + file + "\"")

                    elif filetype == "file":
                        self.itp_reader(file, recursion_layer = 0, write_includes = write_includes)

            self.print_term("Finished loading topologies. Number of moleculetypes found:", len(self.itp_moleculetypes), "\n", spaces=1, verbose=2)
    
