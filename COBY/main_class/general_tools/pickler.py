import pickle

class pickler:
    def pickler(self):
        if self.PICKLE_cmd:
            self.print_term("-------------------", verbose=1)
            self.print_term("Pickling data into:", self.PICKLE_cmd, verbose=1)
            self.print_term("-------------------", "\n", verbose=1)
            
            ### ### To unpickle use the following and change "pickled_file_path" to your pickled file name
            ### with open(pickled_file_path, 'rb') as pickled_file:
            ###     unpickled_class = pickle.load(pickled_file)
            
            pickled_data = self
            if self.backup:
                self.backupper(self.PICKLE_cmd)
            with open(self.PICKLE_cmd, 'wb') as f:
                pickle.dump(pickled_data, f)

