
class topology_entry_generic:
    def __init__(self, value_dict, itp_if, itp_defs, topology_type):

        defs_pointer_dict = {
            "atoms":       "atomtypes",
            "bonds":       "bondtypes",
            "pairs":       "pairtypes",
            "angles":      "angletypes",
            "dihedrals":   "dihedraltypes",
            "constraints": "constrainttypes",
        }
        self.values = []
        self.values_dict = {}
        for key, val in value_dict.items():
            if type(val) == str:
                if topology_type in defs_pointer_dict.keys() and val in itp_defs[defs_pointer_dict[topology_type]]:
                    val = itp_defs[defs_pointer_dict[topology_type]]
            
            elif type(val) in [list, tuple]:
                new_val = []
                for subval in val:
                    if topology_type in defs_pointer_dict.keys() and subval in itp_defs[defs_pointer_dict[topology_type]]:
                        subval = itp_defs[defs_pointer_dict[topology_type]]
                    new_val.append(subval)
                if type(val) == list:
                    val = new_val
                elif type(val) == tuple:
                    val = tuple(new_val)
            
            setattr(self, key, val)

            self.values.append(val)
            self.values_dict[key] = val

        self.itp_if = itp_if
    
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

    def add_entry(self, entry, entry_id, itp_if, itp_defs):
        if self.topology_type == "atoms":
            keys = ["id", "atomtype", "resnr", "resid", "atom", "cgnr", "charge", "mass"]
            ### zip ensures that "mass" is exempted if the value is not given in the topology
            entry_dict = {key: value for key, value in zip(keys, entry)}
        else:
            entry_dict = {
                "atoms":      entry[:self.number_of_atoms],
                "parameters": entry[self.number_of_atoms:], # parameters includes "func_type"
            }
        if entry_id is False:
            if len(self.entries.keys()) == 0:
                entry_id = 1
            else:
                max_id = max(list(self.entries.keys()))
                entry_id = max_id + 1
        
        self.entries[entry_id] = topology_entry_generic(entry_dict, itp_if, itp_defs, self.topology_type)

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
    
    def get_itp_ifs(self, itp_if_type="all"):
        ### Examples of "itp_if_type"
        ### False
        ### ("#ifdef", "M3_CBT")
        ### (("#ifdef", "M3_CBT"), "#else")
        ### ("#ifndef", "M3_CBT")
        ### (("#ifndef", "M3_CBT"), "#else")

        itp_if_dict = {}
        for entry_id, entry in self.entries.items():
            if (itp_if_type == "all") or (itp_if_type == entry.itp_if):
                if entry.itp_if not in if_dict.keys():
                    itp_if_dict[entry.itp_if] = {}
                itp_if_dict[entry.itp_if].update({entry_id: entry})
        return itp_if_dict

class MOLECULETYPE:
    def __init__(self, moleculetype):
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
        self.topology_types_in_molecule[topology_type] = topology_type_generic(topology_type)
    
    def add_entry(self, topology_type, entry, entry_id=False, itp_if=False, itp_defs=False):
        if type(entry) == str:
            entry = entry.split()
        elif type(entry) == tuple:
            entry = list(entry)
        self.topology_types_in_molecule[topology_type].add_entry(entry=entry, entry_id=entry_id, itp_if=itp_if, itp_defs=itp_defs)
        
    def get_total_charge(self):
        assert "atoms" in self.topology_types_in_molecule.keys(), "'atoms' not defined for moleculetype: " + str(self.moleculetype)
        charges = self.topology_types_in_molecule["atoms"].get_value_type("charge")
        total_charge = sum(map(float, charges))
        return total_charge
