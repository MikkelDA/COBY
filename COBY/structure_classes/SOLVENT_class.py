from COBY.structure_classes.MOLECULE_class import MOLECULE

class SOLVENT(MOLECULE):
    def __init__(self, mapping_ratio = 1, molname = False, moleculetype = False):
        super().__init__(molname = molname, moleculetype = False)
        self.mapping_ratio = mapping_ratio ### AA-to-CG convertion
        self.molarity = False
        self.density = False
        self.molar_mass = False
        
    def mapping_ratio_set(self, mapping_ratio):
        self.mapping_ratio = mapping_ratio

    def molarity_set(self, molarity):
        self.molarity = molarity
    
    def density_set(self, density):
        self.density = density
    
    def ratio_set(self, ratio):
        self.ratio = ratio
    
    def molar_mass_set(self, molar_mass):
        self.molar_mass = molar_mass
    
    def count_set(self, count):
        self.count = count
