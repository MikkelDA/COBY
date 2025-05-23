from COBY.structure_classes.MOLECULE_class import MOLECULE

class LIPID(MOLECULE):
    def __init__(self, molname = False, moleculetype = False):
        super().__init__(molname = molname, moleculetype = moleculetype)
        self.ratio      = 0
        self.apl        = 0.6
        self.extra_type = "addition"
        self.extra_val  = 0

    
    def ratio_add(self, ratio):
        self.ratio = ratio

    def apl_add(self, apl):
        self.apl = apl

    def extra_type_add(self, extra_type):
        self.extra_type = extra_type

    def extra_val_add(self, extra_val):
        self.extra_val = extra_val
