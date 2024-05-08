import os

class backupper:
    def backupper(self, output_file_name):
        """
        Checks if output file already exists and backs it up
        """
        ### Makes sure original path is fine
        output_file_name = os.path.join(output_file_name)
        # print(output_file_name)
        if os.path.exists(output_file_name):
            output_file_split = output_file_name.split("/")
            # print(output_file_split, output_file_split[:-1], *output_file_split[:-1])
            if len(output_file_split) > 1:
                output_path = os.path.join(*output_file_split[:-1])
            else:
                output_path = os.path.join("")
            output_name = output_file_split[-1]
            self.print_term("File " + '"' + output_file_name + '"' + " already exists. Backing it up", spaces=0, verbose=1)
            number = 1
            while True:
                if os.path.exists(os.path.join(output_path, "#" + output_name + "." + str(number) + "#")):
                    number += 1
                else:
                    os.rename(output_file_name, os.path.join(output_path, "#" + output_name + "." + str(number) + "#"))
                    break

