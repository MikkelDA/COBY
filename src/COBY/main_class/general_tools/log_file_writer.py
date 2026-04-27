class log_file_writer:
    def log_file_writer(self):
        if self.LOG_FILE_HANDLE and not self.LOG_FILE_HANDLE.closed:
            self.print_term("Closing log file", spaces=0, verbose=1)
            self.LOG_FILE_HANDLE.close()
            self.LOG_FILE_HANDLE = None