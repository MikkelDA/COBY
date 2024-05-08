import os
from COBY.main_class.topology_handlers.moleculetype_class import MOLECULETYPE

class itp_reader:
    def itp_reader(self, itp_file_dest, recursion_layer, write_includes):
        '''
        Reads itp files
        Understands itp file definitions using "#ifdef" and "#endif"
        Understands references to other files using "#include" (by calling itself recursively on the files)
        Currently only charges are used for anything
        '''
        def line_filter(string):
            return list(filter(None, string.split()))

        interactions_name_list = [
            "bonds",
            "pairs",
            "pairs_nb",
            "angles",
            "dihedrals",
            "exclusions",
            "constraints",
            "settles",
            "virtual_sites1",
            "virtual_sites2",
            "virtual_sites3",
            "virtual_sites4",
            "virtual_sitesn",
            "position_restraints",
            "distance_restraints",
            "dihedral_restraints",
            "orientation_restraints",
            "angle_restraints",
            "angle_restraints_z",
        ]
        types_name_list = [
            "atomtypes",
            "bondtypes",
            "pairtypes",
            "angletypes",
            "dihedraltypes",
            "constrainttypes",
        ]

        with open(itp_file_dest, "r") as input_file:
            moleculetype = False
            topology_type = False
            itp_if = False

            for line_nr, line in enumerate(input_file):
                line = line.split(";")[0] # Removes comments
                line = line.rstrip("\n")  # Removes "\n" from end of string
                line_values = line_filter(line)

                ### space-only lines or empty lines are skipped
                if line.isspace() or line == "":
                    continue
                
                ### Checks for "#ifdef" and "#ifndef" statements
                elif line_values[0] in ["#ifdef", "#ifndef"]:
                    itp_if = tuple(line)
                elif line_values[0] == "#else":
                    itp_if = (itp_if, line[0])
                elif line_values[0] == "#endif":
                    itp_if = False

                ### Recursively calls the function for '#include' statements
                elif line.startswith("#include"):
                    inc_path = os.path.join("/".join(itp_file_dest.split("/")[:-1]), line_values[1].replace('"', ''))
                    
                    ### Adds include statement to written topology file
                    if recursion_layer == 0 and write_includes:
                        self.TOP_include_statements.append(line)
                    
                    self.itp_reader(inc_path, recursion_layer = recursion_layer+1, write_includes = write_includes)
                
                ### Checks if new topology entry type
                elif line.startswith("["):
                    ### Only interactions remain after "moleculetype" and "atoms"
                    topology_type = line[line.find("[")+len("["):line.rfind("]")].replace(" ", "")
                    
                    if topology_type in ["system", "molecules"]:
                        break
                    
                    elif topology_type in ["atoms"] + interactions_name_list:
                        self.itp_moleculetypes[moleculetype].add_topology_type(topology_type)

                ### Continue if unnecessary data types are currently being run through
                elif topology_type in ["defaults", "nonbond_params"]:
                    continue

                ### '#defines'
                elif topology_type in types_name_list and line.startswith("#define"):
                    self.itp_defs[topology_type][line_values[1]] = line_values[2:]

                elif topology_type == "moleculetype":
                    ### Second value is the number of excluded neighbors. We don't need to think about that.
                    moleculetype = line_values[0]
                    self.itp_moleculetypes[moleculetype] = MOLECULETYPE(moleculetype)
                
                elif topology_type == "atoms":
                    entry_id = line_values[0]
                    self.itp_moleculetypes[moleculetype].add_entry(topology_type=topology_type, entry=line_values, entry_id=entry_id, itp_if=itp_if, itp_defs=self.itp_defs)
                
                elif topology_type in interactions_name_list:
                    entry_id = line_values[0]
                    self.itp_moleculetypes[moleculetype].add_entry(topology_type=topology_type, entry=line_values, itp_if=itp_if, itp_defs=self.itp_defs)
