class gro_atom_writer:
    def gro_atom_writer(self, r_nr, r_name, a_name, a_nr, x, y, z, vx, vy, vz):
        GCL = [5, 5, 5, 5, 8, 8, 8, 8, 8, 8] # gro atom column lengths
        if len(r_name) > 5:
            r_name = r_name[:5]
        if len(a_name) > 5:
            a_name = a_name[:5]
        string = '{r_nr:>{L0}}{r_name:<{L1}}{a_name:>{L2}}{a_nr:>{L3}}{x:>{L4}.3f}{y:>{L5}.3f}{z:>{L6}.3f}{vx:>{L7}}{vy:>{L8}}{vz:>{L9}}'.format(
            r_nr = r_nr, L0 = GCL[0], # int
            r_name = r_name, L1 = GCL[1], # str
            a_name = a_name, L2 = GCL[2], # str
            a_nr = a_nr, L3 = GCL[3], # int
            x = float(x), L4 = GCL[4], # float
            y = float(y), L5 = GCL[5], # float
            z = float(z), L6 = GCL[6], # float
            vx = vx, L7 = GCL[7], # float # str placeholder
            vy = vy, L8 = GCL[8], # float # str placeholder
            vz = vz, L9 = GCL[9], # float # str placeholder
        )
        return string
    
