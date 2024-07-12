from COBY.structure_classes.LIPID_class import LIPID

class lipid_defs_preprocessor:
    def lipid_defs_preprocessor(self):
        '''
        Preprocesses lipid defintions by rearranging the data into a class
        '''
        self.print_term("Preprocessing lipid definitions", verbose=2)
        molecules_for_preprocessing = []
        if hasattr(self, "lipid_defs") and self.lipid_defs:
            molecules_for_preprocessing.append(("defs", self.lipid_defs))

        if hasattr(self, "lipid_defs_built") and self.lipid_defs_built:
            molecules_for_preprocessing.append(("defs", self.lipid_defs_built))

        if hasattr(self, "lipid_defs_imported") and self.lipid_defs_imported:
            molecules_for_preprocessing.append(("import", self.lipid_defs_imported))

        tot_lipids = 0
        for structure_origin, defs in molecules_for_preprocessing:
            for params, lipid_type_dict in defs.items():
                if params not in self.lipid_dict.keys():
                    self.lipid_dict[params] = {}
                
                ### Loop over indexes
                for cur_name, cur_dict in lipid_type_dict.items():
                    tot_lipids += 1
                    if cur_name not in self.lipid_dict[params].keys():
                        self.lipid_dict[params][cur_name] = LIPID(molname = cur_name)
                    else:
                        self.print_term(
                            "The name \"{name}\" has been detected multiple times for lipids within the same parameter library.".format(name=cur_name),
                            "If multiple lipids have the same name then only the first is loaded into the library.",
                            warn=True,
                        )
                        continue
                    
                    self.molecule_defs_checker(self.lipid_dict[params][cur_name], cur_name, cur_dict, type_of_molecule = "lipid", structure_origin=structure_origin)
                
        self.print_term("Number of lipids preprocessed:", tot_lipids, "\n", spaces=1, verbose=2)
