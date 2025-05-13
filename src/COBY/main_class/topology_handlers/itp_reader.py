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

        def MFB_checker(cmd, moleculetype):
            if "params:" not in cmd:
                cmd = " ".join([cmd, "params:TOP"])
            if "name:" not in cmd:
                cmd = " ".join([cmd, "name:{name}".format(name=moleculetype)])
            if "moleculetype:" not in cmd:
                cmd = " ".join([cmd, "moleculetype:{moleculetype}".format(moleculetype=moleculetype)])
            return cmd

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
            itp_if_lineskip = False
            MFB_cmd = False

            for line_nr, line in enumerate(input_file):
                
                if line.startswith(";@COBY") or (line.lstrip().startswith(";") and len(line.split()) > 1 and line.split()[1].startswith("@COBY")):
                    MFB_cmd = " ".join([line.split("@COBY")[1].lstrip(" \"'").rstrip(" \"'"), "FromTopology:True"])
                    ### Is only directly added to molecule builder argument list if defined within the [ moleculetype ] section.
                    ### Else wait for next moleculetype to be defined.
                    if topology_type == "moleculetype":
                        MFB_cmd = MFB_checker(MFB_cmd, moleculetype)
                        self.MOLECULE_FRAGMENT_BUILDER_cmds.append(MFB_cmd)
                        MFB_cmd = False
                
                line = line.split(";")[0] # Removes comments
                line = line.rstrip("\n")  # Removes "\n" from end of string
                line_values = line_filter(line)

                ### space-only lines or empty lines are skipped
                if line.isspace() or line == "":
                    continue
                
                ### Checks for "#ifdef" and "#ifndef" statements
                elif line_values[0] in ["#ifdef", "#ifndef", "#else", "#endif"]:
                    ### If "#ifdef" and not defined then skip lines
                    if line_values[0] == "#ifdef" and line_values[1] not in self.itp_defs_all_defnames:
                        itp_if_lineskip = True

                    ### If "#ifndef" and defined then skip lines
                    elif line_values[0] == "#ifndef" and line_values[1] in self.itp_defs_all_defnames:
                        itp_if_lineskip = True

                    ### Else don't skip lines
                    ### Activates both if any of the above are False or if the line is "#else" or "#endif"
                    else:
                        itp_if_lineskip = False

                ### Skips if "#ifdef" or "#ifndef" should not be processed
                elif itp_if_lineskip:
                    continue

                ### Recursively calls the function for '#include' statements
                elif line.startswith("#include"):
                    inc_path = os.path.join("/".join(itp_file_dest.split("/")[:-1]), line_values[1].replace('"', ''))
                    
                    ### Adds include statement to written topology file
                    if recursion_layer == 0 and write_includes:
                        inc_statement = '#include "{inc_path}"'.format(inc_path=inc_path)
                        self.TOP_include_statements.append(inc_statement)
                    
                    self.itp_reader(inc_path, recursion_layer = recursion_layer+1, write_includes = write_includes)
                
                ### Checks if new topology entry type
                elif line.lstrip(" ").startswith("["):
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
                    self.itp_defs_all_defnames.add(line_values[1])

                elif topology_type == "moleculetype":
                    ### Second value is the number of excluded neighbors. We don't need to think about that.
                    moleculetype = line_values[0]
                    self.itp_moleculetypes[moleculetype] = MOLECULETYPE(moleculetype)

                    ### MFB argument given outside of [ moleculetype ] section. Add argument with next designated [ moleculetype ].
                    if MFB_cmd:
                        MFB_cmd = MFB_checker(MFB_cmd, moleculetype)
                        self.MOLECULE_FRAGMENT_BUILDER_cmds.append(MFB_cmd)
                        MFB_cmd = False
                
                elif topology_type == "atoms":
                    entry_id = line_values[0]
                    self.itp_moleculetypes[moleculetype].add_entry(topology_type=topology_type, entry=line_values, entry_id=entry_id, itp_defs=self.itp_defs)#, itp_if=itp_if)
                
                elif topology_type in interactions_name_list:
                    entry_id = line_values[0]
                    self.itp_moleculetypes[moleculetype].add_entry(topology_type=topology_type, entry=line_values, itp_defs=self.itp_defs)#, itp_if=itp_if)
