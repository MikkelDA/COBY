class topol_file_writer:
    def topol_file_writer(self):
        if self.output_topol_file_name:
            string = " ".join(["", "Writing topology file", ""])
            self.print_term("{string:-^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), spaces=0, verbose=1)
            output_topol_file_lines = []

            for INCLUDE_statement in self.TOP_include_statements:
                if INCLUDE_statement.endswith("\n"):
                    INCLUDE_statement = INCLUDE_statement[:-1]
                output_topol_file_lines.append(INCLUDE_statement)

            output_topol_file_lines.extend([
                "",
                "[ system ]",
                "; name",
                self.system_name or "PLACEHOLDER_TITLE",
                "",
                "[ molecules ]",
                "; name number",
            ])
            molecules_for_top_lengths = [len(str(max(data, key = lambda d: len(str(d))))) for data in zip(*self.molecules_for_top)]
            for n, c in self.molecules_for_top:
                output_topol_file_lines.append(
                    '{NAME:<{Ln}} {COUNT:>{Lc}}'.format(
                        NAME = n, Ln = molecules_for_top_lengths[0],
                        COUNT = c, Lc = molecules_for_top_lengths[1],
                    )
                )
            if self.backup:
                self.backupper(self.output_topol_file_name)

            new_file = open(self.output_topol_file_name, "w")
            for line in output_topol_file_lines:
                new_file.write(line + "\n")
            new_file.close()
            self.print_term("Topology file written:", self.output_topol_file_name, "\n", verbose=1)
    
