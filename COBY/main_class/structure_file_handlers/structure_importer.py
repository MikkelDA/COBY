class structure_importer:
    def structure_importer(self, file_name):
        ### First check if file ends with .pdb or .gro
        if file_name.endswith(".pdb"):
            structure = self.pdb_reader(file_name)
        elif file_name.endswith(".gro"):
            structure = self.gro_reader(file_name)
        elif file_name.endswith(".cif") or file_name.endswith(".mmcif"):
            structure = self.cif_reader(file_name)

        ### If neither then check if it contains .pdb or .gro (e.g. #file.pdb.1#)
        ### Done after file extension checks due to possible false-positives (e.g. pdb_struct_to_gro.gro)
        elif ".pdb" in file_name.lower():
            structure = self.pdb_reader(file_name)
        elif ".gro" in file_name.lower():
            structure = self.gro_reader(file_name)
        elif ".cif" in file_name.lower() or ".mmcif" in file_name.lower():
            structure = self.cif_reader(file_name)
        
        ### Finally assume .pdb format if neither .pdb nor .gro is found
        else:
            structure = self.pdb_reader(file_name)
        
        return structure
    
