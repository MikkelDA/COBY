from COBY.structure_classes.MOLECULE_class import MOLECULE

class LIPID(MOLECULE):
    def __init__(self, molname = False, moleculetype = False):
        super().__init__(molname = molname, moleculetype = moleculetype)
        self.ratio = 0
        self.apl = 0.6
    
    def ratio_add(self, ratio):
        self.ratio = ratio

    def apl_add(self, apl):
        self.apl = apl
