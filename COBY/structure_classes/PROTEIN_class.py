from COBY.structure_classes.MOLECULE_class import MOLECULE

class PROTEIN(MOLECULE):
    def __init__(self, molname = False, moleculetype = False):
        super().__init__(molname = molname, moleculetype = False)
