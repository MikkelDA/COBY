from COBY.structure_classes.LIPID_class import LIPID

class lipid_defs_preprocessor:
    def lipid_defs_preprocessor(self):
        '''
        Preprocesses lipid defintions by rearranging the data into a class
        '''
        self.print_term("Preprocessing lipid definitions", verbose=2)
        for structure_origin, defs in [("defs", self.lipid_defs), ("import", self.lipid_defs_imported)]:
            for params, lipid_type_dict in defs.items():
                if params not in self.lipid_dict.keys():
                    self.lipid_dict[params] = {}
                
                ### Loop over indexes
                for cur_name, cur_dict in lipid_type_dict.items():
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
                
        tot_lipids = sum([len(vals) for vals in self.lipid_dict.values()])
        self.print_term("Number of lipids preprocessed:", tot_lipids, "\n", spaces=1, verbose=2)
