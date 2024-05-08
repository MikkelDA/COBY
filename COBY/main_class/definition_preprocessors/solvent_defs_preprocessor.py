from COBY.structure_classes.SOLVENT_class import SOLVENT

class solvent_defs_preprocessor:
    def solvent_defs_preprocessor(self):
        '''
        Preprocesses solvent defintions by rearranging the data into a class
        '''
        self.print_term("Preprocessing solvent definitions", verbose=2)
        for structure_origin, defs in [("defs", self.solvent_defs), ("import", self.solvent_defs_imported)]:
            for params, solvent_type_dict in defs.items():
                if params not in self.solvent_dict.keys():
                    self.solvent_dict[params] = {}
                
                ### Loop over indexes
                for cur_name, cur_dict in solvent_type_dict.items():
                    if cur_name not in self.solvent_dict[params].keys():
                        self.solvent_dict[params][cur_name] = SOLVENT(mapping_ratio = 1, molname = cur_name)
                    else:
                        self.print_term(
                            "The name \"{name}\" has been detected multiple times for solvents within the same parameter library.".format(name=cur_name),
                            "If multiple solvents have the same name then only the first is loaded into the library.",
                            warn=True,
                        )
                        continue
                    
                    self.molecule_defs_checker(self.solvent_dict[params][cur_name], cur_name, cur_dict, type_of_molecule = "solvent", structure_origin=structure_origin)
                    
        tot_solvents = sum([len(vals) for vals in self.solvent_dict.values()])
        self.print_term("Number of solvents preprocessed:", tot_solvents, "\n", spaces=1, verbose=2)
