import os

class backupper:
    def backupper(self, output_file_name):
        """
        Checks if output file already exists and backs it up
        """
        ### Makes sure original path is fine
        output_file_name = os.path.join(output_file_name)
        self.print_term("output_file_name", output_file_name, debug=True, debug_keys=["backupper", "backup"])
        if os.path.exists(output_file_name):
            output_file_split = output_file_name.split("/")
            self.print_term("output_file_split", output_file_split, debug=True, debug_keys=["backupper", "backup"])
            if len(output_file_split) > 1:
                output_path = os.path.join(*output_file_split[:-1])
            else:
                output_path = os.path.join("")

            ### In case full path is given
            if output_file_name.startswith("/"):
                output_path = "/" + output_path
            
            output_name = output_file_split[-1]
            self.print_term("output_name", output_name, debug=True, debug_keys=["backupper", "backup"])
            self.print_term("File " + '"' + output_file_name + '"' + " already exists. Backing it up", spaces=1, verbose=1)
            number = 1
            while True:
                renamed_file = os.path.join(output_path, "#" + output_name + "." + str(number) + "#")
                self.print_term("renamed_file", renamed_file, debug=True, debug_keys=["backup", "backupper"])
                if os.path.exists(renamed_file):
                    number += 1
                else:
                    os.rename(output_file_name, renamed_file)
                    break

