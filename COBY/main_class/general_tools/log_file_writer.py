class log_file_writer:
    def log_file_writer(self):
        if self.output_log_file_name:
            string = " ".join(["", "Writing log file", ""])
            self.print_term("{string:-^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), spaces=0, verbose=1)
            if self.backup:
                self.backupper(self.output_log_file_name)
            new_file = open(self.output_log_file_name, "w")
            for line in self.LOG_FILE:
                new_file.write(line)
            new_file.close()
            self.print_term("Log file written:", self.output_log_file_name, "\n", verbose=1)

