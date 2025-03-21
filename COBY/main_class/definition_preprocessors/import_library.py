import importlib

class import_library():
    def import_library(self, subcmd):
        spec = importlib.util.spec_from_file_location(
            name="defs_module",
            location=subcmd,
        )
        defs_module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(defs_module)

        def DictDeepMerger(d1, d2, cur_d=""):
            '''
            Inspired by / adapted from https://stackoverflow.com/questions/7204805/deep-merge-dictionaries-of-dictionaries-in-python/7205107#7205107
            
            Takes two dicts, 'd1' and 'd2', and recursively merges them into 'd1'. Does not return any values as 'd1' is modified.
            '''
            for key in d2:
                if key in d1:
                    if type(d1[key]) == dict and type(d2[key]) == dict:
                        DictDeepMerger(d1[key], d2[key], cur_d + str([key]))
                    elif d1[key] != d2[key]:
                        d1[key] = d2[key]
                        self.print_term("Duplicate dictionary key '{key}' with non-identical value found for '{cur_d}'. Will overwrite prior data.".format(key=key, cur_d=cur_d), warn=True)
                else:
                    d1[key] = d2[key]
        
        for defstype in ["lipid_scaffolds", "lipid_defs", "solvent_defs", "ion_defs", "prot_defs", "fragment_defs"]:
            ### "and hasattr(self, defstype)" prevents "prot_defs" from being run for MoleculeBuilder program
            if hasattr(defs_module, defstype) and hasattr(self, defstype):
                DictDeepMerger(getattr(self, defstype), getattr(defs_module, defstype), "lipid_scaffolds")
        