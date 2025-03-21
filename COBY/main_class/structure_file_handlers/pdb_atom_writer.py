class pdb_atom_writer:
    def pdb_atom_writer(self, ATOM, a_nr, a_name, aLoc, r_name, chain, r_nr, aChar, x, y, z, oc, T, JUMP, E, C):
        '''
        PCL: pdb ATOM column lengths. "+" means that multiple column groups are combined
        '''
        PCL = [6, 5, 4, 1, 3 + 1, 1, 4, 1 + 3, 8, 8, 8, 6, 6, 10, 2, 2]
        ### Missing "+1" at third position is added manually as "col12spacer"
        if len(r_name) > 4:
            r_name = r_name[:4]
        if len(r_name) < 4:
            r_name = r_name + " "
        if len(a_name) > 5:
            a_name = a_name[:5]
        if len(a_name) < 4:
            a_name = " " + a_name
        string = '{ATOM:<{L0}}{a_nr:>{L1}}{col12spacer:>{LC12S}}{a_name:<{L2}}{aLoc:>{L3}}{r_name:>{L4}}{chain:^{L5}}{r_nr:>{L6}}{aChar:>{L7}}{x:>{L8}.3f}{y:>{L9}.3f}{z:>{L10}.3f}{oc:>{L11}.2f}{T:>{L12}.2f}{JUMP:>{L13}}{E:>{L14}}{C:>{L15}}'.format(
            ATOM        = ATOM,          L0  = PCL[0],
            a_nr        = a_nr,          L1  = PCL[1], # int
            col12spacer = " ",           LC12S = 1,
            a_name      = a_name,        L2  = PCL[2],
            aLoc        = aLoc,          L3  = PCL[3],
#             aLoc        = "",            L3  = 0,
            r_name      = r_name,        L4  = PCL[4],
            chain       = chain,         L5  = PCL[5],
            r_nr        = r_nr,          L6  = PCL[6], # int
            aChar       = aChar,         L7  = PCL[7],
            x           = float(x),      L8  = PCL[8], # float
            y           = float(y),      L9  = PCL[9], # float
            z           = float(z),      L10 = PCL[10], # float
            oc          = float(oc),     L11 = PCL[11], # float
            T           = float(T),      L12 = PCL[12], # float
            JUMP        = JUMP,          L13 = PCL[13],
            E           = E,             L14 = PCL[14],
            C           = C,             L15 = PCL[15],
        )
        return string


