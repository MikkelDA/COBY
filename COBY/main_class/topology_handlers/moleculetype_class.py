
class topology_entry_generic:
    def __init__(self, value_dict, topology_type, itp_defs):#, itp_if):
        '''
        Generic topology entry class

        value_dict: [dict] of key:value pairs
            For "atoms" topology type it contains "id", "atomtype", "resnr", "resid", "atom", "cgnr", "charge", (and potentially "mass")
            For all other topology types it contains "atoms" and "parameters"

        topology_type: [str], Must be one of the following topology type entries:
            "atoms", "bonds", "pairs", "pairs_nb", "angles", "dihedrals", "exclusions", "constraints", "settles",
            "virtual_sites1", "virtual_sites2", "virtual_sites3", "virtual_sites4", "virtual_sitesn",
            "position_restraints", "distance_restraints", "dihedral_restraints", "orientation_restraints", "angle_restraints", "angle_restraints_z",
        '''
        ### Pointer list for definitions dictionary
        defs_pointer_dict = {
            "atoms":       "atomtypes",
            "bonds":       "bondtypes",
            "pairs":       "pairtypes",
            "angles":      "angletypes",
            "dihedrals":   "dihedraltypes",
            "constraints": "constrainttypes",
        }
        ### Runs through each value in the value dict to see if values are explicitly given or if they are defined using variables (defs)
        self.values = []
        self.values_dict = {}
        for key, val in value_dict.items():
            ### If single value then it is a string
            if type(val) == str:
                ### If value is defined using variable (defs) then look for it in definitions dict for the given topology type
                if topology_type in defs_pointer_dict.keys() and val in itp_defs[defs_pointer_dict[topology_type]]:
                    val = itp_defs[defs_pointer_dict[topology_type]]
            
            ### If value is a [list or tuple] of values, then iterate over each and add to new value list
            elif type(val) in [list, tuple]:
                new_val = []
                for subval in val:
                    ### If value is defined using variable (defs) then look for it in definitions dict for the given topology type
                    if topology_type in defs_pointer_dict.keys() and subval in itp_defs[defs_pointer_dict[topology_type]]:
                        subval = itp_defs[defs_pointer_dict[topology_type]]
                    
                    ### Adding subval to new_val list
                    if type(subval) in [str, int, float]:
                        new_val.append(subval)
                    ### Definitions can technically be lists so should be extended onto new_val instead of appended
                    elif type(subval) == list:
                        new_val.extend(subval)

                if type(val) == list:
                    val = new_val
                elif type(val) == tuple:
                    val = list(new_val)
            
            ### Sets class attribute
            setattr(self, key, val)

            ### Sets values_dict attribute
            self.values.append(val)
            self.values_dict[key] = val

        # ### Sets whether entry is ifdef or not
        # self.itp_if = itp_if
    
    def get_list(self):
        List = []
        for val in self.values:
            if type(val) in [list, tuple]:
                for subval in val:
                    if subval is not False:
                        List.append(subval)
            else:
                if val is not False:
                    List.append(val)
        return List
    
    def get_string(self, lengths=[]):
        List = self.get_list()
        if len(lengths) == 0:
            string = " ".join(List)
        else:
            assert len(lengths) == len(List), "\n".join([
                "Number of string 'lengths' is not the same as the number of values to be turned into string.",
                "Current topology: " + self.topology_type,
                "value list: " + str(List),
                "lengths list: " + str(lengths),
            ])
            
            string = ""
            for val, length in zip(List, lengths):
                string = string + val.rjust(length, " ")

        return string

class topology_type_generic:
    def __init__(self, topology_type):
        '''
        Generic topology type class

        topology_type: [str], Must be one of the following topology type entries:
            "atoms", "bonds", "pairs", "pairs_nb", "angles", "dihedrals", "exclusions", "constraints", "settles",
            "virtual_sites1", "virtual_sites2", "virtual_sites3", "virtual_sites4", "virtual_sitesn",
            "position_restraints", "distance_restraints", "dihedral_restraints", "orientation_restraints", "angle_restraints", "angle_restraints_z",
        '''
        self.topology_type = topology_type
        self.entries = {}

        if self.topology_type != "atoms":
            ### Using dictionary lookup instead of 19 individual if/elif statements to set number of atoms
            number_of_atoms_dict = {
                "bonds"                  : 2,
                "pairs"                  : 2,
                "pairs_nb"               : 2,
                "angles"                 : 3,
                "dihedrals"              : 4,
                "exclusions"             : 1,
                "constraints"            : 2,
                "settles"                : 1,
                "virtual_sites1"         : 2,
                "virtual_sites2"         : 3,
                "virtual_sites3"         : 4,
                "virtual_sites4"         : 5,
                "virtual_sitesn"         : 1,
                "position_restraints"    : 1,
                "distance_restraints"    : 2,
                "dihedral_restraints"    : 4,
                "orientation_restraints" : 2,
                "angle_restraints"       : 4,
                "angle_restraints_z"     : 2,
            }
            self.number_of_atoms = number_of_atoms_dict[topology_type]

    def add_entry(self, entry, entry_id, itp_defs):#, itp_if):
        '''
        entry: [list] of entry values
        
        entry_id: [int or False], Numbered entry ID. Entry ID will be generated if not given.
        '''
        ### Atoms must be treated differently than all other topology types
        if self.topology_type == "atoms":
            keys = ["id", "atomtype", "resnr", "resid", "atom", "cgnr", "charge", "mass"]
            ### zip ensures that "mass" is exempted if the value is not given in the topology
            entry_dict = {key: value for key, value in zip(keys, entry)}
        ### All other topology types
        else:
            ### Entry dict split into affected atoms and entry parameters
            entry_dict = {
                "atoms":      entry[:self.number_of_atoms],
                "parameters": entry[self.number_of_atoms:], # parameters includes "func_type"
            }

        ### Generate entry id if it is not given
        if entry_id is False:
            ### If first entry, start at 1 (don't use python 0-indexing)
            if len(self.entries.keys()) == 0:
                entry_id = 1
            ### Otherwise just add 1 to max entry id
            else:
                max_id = max(list(self.entries.keys()))
                entry_id = max_id + 1
        
        ### Define the entry using the generic entry class
        self.entries[entry_id] = topology_entry_generic(entry_dict, self.topology_type, itp_defs)#, itp_if)

    def get_value_type(self, value_type):
        List = []
        for entry in self.entries.values():
            if hasattr(entry, value_type):
                List.append(getattr(entry, value_type))
        return List

    def get_entries_dict(self):
        return self.entries
    
    def get_entries_list(self):
        return [val for val in self.entries.values()]
    
    def get_ids(self):
        return [key for key in self.entries.keys()]
    
    # def get_itp_ifs(self, itp_if_type="all"):
    #     '''
    #     Gets all itp ifdefs and ifndefs.

    #     Currently unused. Kept in case it is needed some day.
    #     '''
    #     ### Examples of "itp_if_type"
    #     ### False
    #     ### ("#ifdef", "M3_CBT")
    #     ### (("#ifdef", "M3_CBT"), "#else")
    #     ### ("#ifndef", "M3_CBT")
    #     ### (("#ifndef", "M3_CBT"), "#else")

    #     itp_if_dict = {}
    #     for entry_id, entry in self.entries.items():
    #         if (itp_if_type == "all") or (itp_if_type == entry.itp_if):
    #             if entry.itp_if not in if_dict.keys():
    #                 itp_if_dict[entry.itp_if] = {}
    #             itp_if_dict[entry.itp_if].update({entry_id: entry})
    #     return itp_if_dict

class MOLECULETYPE:
    def __init__(self, moleculetype):
        '''
        Initializes a new MOLECULETYPE

        moleculetype: [str], Name of the moleculetype (molecule name)
        '''
        ### Sets name / moleculetype
        self.moleculetype               = moleculetype
        self.topology_types_in_molecule = {}

        ### Would create all topology types at object creation. Probably not necessary.
        ### Less data clutter if topology types are created as they are needed.
        ### Kept commented out for overview
        # topology_types = [
        #     "atoms",
        #     "bonds",
        #     "pairs",
        #     "pairs_nb",
        #     "angles",
        #     "dihedrals",
        #     "exclusions",
        #     "constraints",
        #     "settles",
        #     "virtual_sites1",
        #     "virtual_sites2",
        #     "virtual_sites3",
        #     "virtual_sites4",
        #     "virtual_sitesn",
        #     "position_restraints",
        #     "distance_restraints",
        #     "dihedral_restraints",
        #     "orientation_restraints",
        #     "angle_restraints",
        #     "angle_restraints_z",
        # ]
        # for topology_type in topology_types:
        #     self.topology_types_in_molecule[topology_type] = topology_type_generic(topology_type)

    def add_topology_type(self, topology_type):
        '''
        Adds new topology type if it does not already exist.

        topology_type: [str], Must be one of the following topology type entries:
            "atoms", "bonds", "pairs", "pairs_nb", "angles", "dihedrals", "exclusions", "constraints", "settles",
            "virtual_sites1", "virtual_sites2", "virtual_sites3", "virtual_sites4", "virtual_sitesn",
            "position_restraints", "distance_restraints", "dihedral_restraints", "orientation_restraints", "angle_restraints", "angle_restraints_z",
        '''
        self.topology_types_in_molecule[topology_type] = topology_type_generic(topology_type)
    
    def add_entry(self, topology_type, entry, entry_id=False, itp_defs=False):#, itp_if=False):
        '''
        Adds a new line entry to the specified topology type under the assumption that the topology type has already been added.

        entry: [str, tuple or list] of values for entry.
            If string or tuple then it is turned into a list of values. Strings are split by whitespaces.
        
        entry_id: [int or False], Numbered entry ID. Entry ID will be generated if False.
        '''
        if type(entry) == str:
            entry = entry.split()
        elif type(entry) == tuple:
            entry = list(entry)
        self.topology_types_in_molecule[topology_type].add_entry(entry=entry, entry_id=entry_id, itp_defs=itp_defs)#, itp_if=itp_if)
        
    def get_total_charge(self):
        assert "atoms" in self.topology_types_in_molecule.keys(), "'atoms' not defined for moleculetype: " + str(self.moleculetype)
        charges = self.topology_types_in_molecule["atoms"].get_value_type("charge")
        total_charge = sum(map(float, charges))
        return total_charge
