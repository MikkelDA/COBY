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
        molecules_for_preprocessing = []
        if hasattr(self, "pos_ion_defs") and self.pos_ion_defs:
            molecules_for_preprocessing.append(("defs", self.pos_ion_defs, "pos_ions"))

        if hasattr(self, "pos_ion_defs_built") and self.pos_ion_defs_built:
            molecules_for_preprocessing.append(("defs", self.pos_ion_defs_built, "pos_ions"))

        if hasattr(self, "pos_ion_defs_imported") and self.pos_ion_defs_imported:
            molecules_for_preprocessing.append(("import", self.pos_ion_defs_imported, "pos_ions"))

        if hasattr(self, "neg_ion_defs") and self.neg_ion_defs:
            molecules_for_preprocessing.append(("defs", self.neg_ion_defs, "neg_ions"))

        if hasattr(self, "neg_ion_defs_built") and self.neg_ion_defs_built:
            molecules_for_preprocessing.append(("defs", self.neg_ion_defs_built, "neg_ions"))

        if hasattr(self, "neg_ion_defs_imported") and self.neg_ion_defs_imported:
            molecules_for_preprocessing.append(("import", self.neg_ion_defs_imported, "neg_ions"))

        for structure_origin, defs, charge_type in molecules_for_preprocessing:
            if charge_type == "pos_ions":
                curr_ion_dict = self.pos_ion_dict
            if charge_type == "neg_ions":
                curr_ion_dict = self.neg_ion_dict

            for params, params_dict in defs.items():
                ions_params = "_".join([charge_type, params])
                if ions_params not in self.solvent_defs.keys():
                    self.solvent_defs[ions_params] = {}

                if params not in curr_ion_dict.keys():
                    curr_ion_dict[params] = {}
            
                ### Loop over indexes
                for cur_name, cur_dict in params_dict.items():
                    if cur_name not in curr_ion_dict[params].keys():
                        ion = SOLVENT(mapping_ratio = 1, molname = cur_name)
                        curr_ion_dict[params][cur_name] = ion
                        ### Adding ion libraries to solvent libraries
                        if charge_type not in self.solvent_defs[ions_params].keys():
                            self.solvent_defs[ions_params][cur_name] = copy.deepcopy(cur_dict)
                        if charge_type not in self.solvent_defs[ions_params].keys():
                            self.solvent_defs[ions_params][cur_name] = copy.deepcopy(cur_dict)
                    else:
                        if charge_type == "pos_ions":
                            ion_name = "positive"
                        if charge_type == "neg_ions":
                            ion_name = "negative"
                        self.print_term(
                            "The name \"{name}\" has been detected multiple times for {ion_name} ions within the same parameter library.".format(name=cur_name, ion_name=ion_name),
                            "If multiple {ion_name} ions have the same name then only the first is loaded into the library.".format(ion_name=ion_name),
                            warn=True,
                        )
                        continue
                    
                    self.molecule_defs_checker(curr_ion_dict[params][cur_name], cur_name, cur_dict, type_of_molecule = "solvent", structure_origin=structure_origin)

        tot_pos_ions = sum([len(params_dict) for params_dict in self.pos_ion_dict.values()])
        self.print_term("Number of positive ions preprocessed:", tot_pos_ions, spaces=1, verbose=2)
        tot_neg_ions = sum([len(params_dict) for params_dict in self.neg_ion_dict.values()])
        self.print_term("Number of negative ions preprocessed:", tot_neg_ions, "\n", spaces=1, verbose=2)
