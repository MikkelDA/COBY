class log_file_writer:
    def log_file_writer(self):
        if self.output_log_file_name:
            self.print_term("Writing log file (LOG)", spaces=0, verbose=1)

            if self.backup:
                self.backupper(self.output_log_file_name)
            new_file = open(self.output_log_file_name, "w")
            for line in self.LOG_FILE:
                new_file.write(line)
            new_file.close()
            self.print_term("LOG file written:", self.output_log_file_name, spaces=1, verbose=1)

