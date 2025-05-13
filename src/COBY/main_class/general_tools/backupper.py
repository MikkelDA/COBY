import os

class backupper:
    def backupper(self, output_file_name):
        """
        Checks if output file already exists and backs it up
        """
        ### Makes sure original path is fine
        output_file_name = os.path.join(output_file_name)
        self.print_term("output_file_name:", output_file_name, debug=True, debug_keys=["backupper", "backup"])
        if os.path.exists(output_file_name):
            output_file_split = os.path.split(output_file_name)
            self.print_term("output_file_split:", output_file_split, debug=True, debug_keys=["backupper", "backup"])
            if len(output_file_split[0]) > 1 and len(output_file_split[1]) > 1:
                output_path = os.path.join(*output_file_split[:-1])
            else:
                output_path = os.path.join("")
            self.print_term("output_path:", output_path, debug=True, debug_keys=["backupper", "backup"])

            output_name = output_file_split[-1]
            self.print_term("output_name:", output_name, debug=True, debug_keys=["backupper", "backup"])
            self.print_term("File " + '"' + output_file_name + '"' + " already exists. Backing it up", spaces=1, verbose=1)
            number = 1
            while True:
                renamed_file = os.path.join(output_path, "#" + output_name + "." + str(number) + "#")
                if os.path.exists(renamed_file):
                    number += 1
                else:
                    self.print_term("renamed_file from:", output_file_name, debug=True, debug_keys=["backup", "backupper"])
                    self.print_term("renamed_file to  :", renamed_file, debug=True, debug_keys=["backup", "backupper"])
                    os.rename(output_file_name, renamed_file)
                    break

