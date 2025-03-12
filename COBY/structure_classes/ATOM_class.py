
class ATOM:
    def __init__(self, bead, beadnr, x, y, z, resname, resnumber, charge=0):
        self.bead   = bead
        self.beadnr = beadnr
        self.x = x
        self.y = y
        self.z = z
        self.resname = resname
        self.resnr   = resnumber
        self.charge  = charge
        
    def movex(self, x):
        self.x = x
    def movey(self, y):
        self.y = y
    def movez(self, z):
        self.z = z
        
    def move_atom(self, x, y, z):
        self.movex(x)
        self.movey(y)
        self.movez(z)
    
    def set_charge(self, charge):
        self.charge = charge
    
    def get_tuple(self):
        return (self.bead, self.beadnr, self.x, self.y, self.z, self.resname, self.resnr, self.charge)
    
