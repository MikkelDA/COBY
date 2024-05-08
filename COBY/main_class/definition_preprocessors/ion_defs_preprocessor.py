from COBY.structure_classes.SOLVENT_class import SOLVENT
import copy

class ion_defs_preprocessor:
    def ion_defs_preprocessor(self):
        '''
        Preprocesses ion defintions by rearranging the data into a class
        '''
        self.print_term("Preprocessing ion definitions", verbose=2)
        ### Creates parameter libraries in the solvent defs for positive and negative ions
        ### In case people want to flood a system with a specific number of ions
        for structure_origin, defs in [("defs", self.ion_defs), ("import", self.ion_defs_imported)]:
            if "pos_ions" not in self.solvent_defs.keys():
                self.solvent_defs["pos_ions"] = {}
            if "neg_ions" not in self.solvent_defs.keys():
                self.solvent_defs["neg_ions"] = {}
                
            for params, charge_type_dict in defs.items():
                if params not in self.ion_dict.keys():
                    self.ion_dict[params] = {}
                for charge_name, ion_type_dict in charge_type_dict.items():
                    if charge_name not in self.ion_dict[params].keys():
                        self.ion_dict[params][charge_name] = {}
                    ### Loop over indexes
                    for cur_name, cur_dict in ion_type_dict.items():
                        if cur_name not in self.ion_dict[params][charge_name].keys():
                            ion = SOLVENT(mapping_ratio = 1, molname = cur_name)
                            self.ion_dict[params][charge_name][cur_name] = ion
                            if params == "default":
                                if charge_name not in self.solvent_defs["pos_ions"].keys():
                                    self.solvent_defs["pos_ions"][cur_name] = copy.deepcopy(cur_dict)
                                if charge_name not in self.solvent_defs["neg_ions"].keys():
                                    self.solvent_defs["neg_ions"][cur_name] = copy.deepcopy(cur_dict)
                        else:
                            if charge_name == "pos_ions":
                                ion_name = "positive"
                            if charge_name == "neg_ions":
                                ion_name = "negative"
                            self.print_term(
                                "The name \"{name}\" has been detected multiple times for {ion_name} ions within the same parameter library.".format(name=cur_name, ion_name=ion_name),
                                "If multiple {ion_name} ions have the same name then only the first is loaded into the library.".format(ion_name=ion_name),
                                warn=True,
                            )
                            continue
                        
                        self.molecule_defs_checker(self.ion_dict[params][charge_name][cur_name], cur_name, cur_dict, type_of_molecule = "solvent", structure_origin=structure_origin)

        tot_pos_ions = sum([len(vals["positive"]) for vals in self.ion_dict.values()])
        self.print_term("Number of positive ions preprocessed:", tot_pos_ions, spaces=1, verbose=2)
        tot_neg_ions = sum([len(vals["negative"]) for vals in self.ion_dict.values()])
        self.print_term("Number of negative ions preprocessed:", tot_neg_ions, spaces=1, verbose=2)
