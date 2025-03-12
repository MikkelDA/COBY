class cif_atom_writer:
    def cif_atom_tupler(self, ATOM, a_nr, a_name, r_name, r_nr, x, y, z):
        '''
        Writes atom tuples for pdbx/mmCIF files

        ### PDB-to-mmCIF translation below:
        _atom_site.group_PDB     = ATOM
        _atom_site.id            = a_nr
        _atom_site.auth_atom_id  = a_name
        _atom_site.auth_comp_id  = r_name
        _atom_site.auth_seq_id   = r_nr
        _atom_site.Cartn_x       = x
        _atom_site.Cartn_y       = y
        _atom_site.Cartn_z       = z
        ### Will always assign "A" as the chain.
        _atom_site.label_asym_id = chain

        ### The following are required to load cif files into ChimeraX for some reason, so they are also written.
        ## ChimeraX still refuses to load CG system properly from cif files, but i will leave the values in, so they are at least present.
        _atom_site.label_atom_id = a_name
        _atom_site.label_comp_id = r_name
        _atom_site.label_seq_id  = r_nr
        ### The following does not make the most sense for CG but is still required by ChimeraX.
        ### It represents the atom type (e.g. C for carbon or O for oxygen). Will just have it write "C" for all atoms.
        _atom_site.type_symbol   = a_name[0]
        '''

        atom_tuple = (
            ATOM,
            str(a_nr),
            a_name,
            r_name,
            str(r_nr),
            str(round(float(x), 3)), # Rounded to 3 decimals to mimic accuracy in .pdb and .gro files.
            str(round(float(y), 3)), # Rounded to 3 decimals to mimic accuracy in .pdb and .gro files.
            str(round(float(z), 3)), # Rounded to 3 decimals to mimic accuracy in .pdb and .gro files.
            "A",

            ### ChimeraX stuff below:
            a_name,
            r_name,
            str(r_nr),
            a_name[0],
        )

        return atom_tuple
    
    def cif_atom_writer(self, atom_tuple, lengths):
        '''
        Writes atom lines for pdbx/mmCIF files

        ### PDB-to-mmCIF translation below:
        _atom_site.group_PDB     = ATOM
        _atom_site.id            = a_nr
        _atom_site.auth_atom_id  = a_name
        _atom_site.auth_comp_id  = r_name
        _atom_site.auth_seq_id   = r_nr
        _atom_site.Cartn_x       = x
        _atom_site.Cartn_y       = y
        _atom_site.Cartn_z       = z
        ### Will always assign "A" as the chain.
        _atom_site.label_asym_id = chain

        ### The following are required to load cif files into ChimeraX for some reason, so they are also written.
        ## ChimeraX still refuses to load CG system properly from cif files, but i will leave the values in, so they are at least present.
        _atom_site.label_atom_id = a_name
        _atom_site.label_comp_id = r_name
        _atom_site.label_seq_id  = r_nr
        ### The following does not make the most sense for CG but is still required by ChimeraX.
        ### It represents the atom type (e.g. C for carbon or O for oxygen). Will just have it write "C" for all atoms.
        _atom_site.type_symbol   = "C"
        '''
        #              0    1    2    3    4    5    6    7    8    9    10   11   12
        centerings = ["<", ">", "<", ">", ">", ">", ">", ">", ">", "<", ">", ">", ">"]
        string = " ".join(['{V:{C}{L}}'.format(V=V, C=C, L=L) for V, C, L in list(zip(atom_tuple, centerings, lengths))])

        return string

