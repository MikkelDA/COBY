import ast
from operator import itemgetter
import copy

class molecule_importer:
    def molecule_importer(self):
        if len(self.MOLECULE_IMPORT_cmds) != 0:
            for imp_struc in self.MOLECULE_IMPORT_cmds:
                
                structure      = False
                params         = "default" # str
                moleculetype   = False
                refname        = False
                name           = False
                upbeads_cmds   = []
                downbeads_cmds = []
                charge_cmds    = []
                library_types  = []
                ### Solvent-specific values
                molar_mass     = False
                density        = False
                mapping_ratio  = False
                ### Used to scale imported x/y/z coordinates
                scalings = {
                    "x": 1,
                    "y": 1,
                    "z": 1,
                }

                ### ### settings dictionary:
                ### name: reference name for solvent selection
                ### charge: whether to use name for topology or given value
                
                ### Check for use of V1 or V2 system. User must choose one version for a given command.
                V1_check = [1 for i in imp_struc.split() if i.startswith("moleculetype:")]
                V2_check = [1 for i in imp_struc.split() if any([i.startswith(j) for j in ["names:", "charge:", "charges:"]])]
                assert not (sum(V1_check) > 0 and sum(V2_check) > 0), "Both V1 and V2 import commands used. Stick to one."
                
                check_structure_files = len([1 for i in imp_struc.split() if i.startswith("file:")]) == 1
                assert check_structure_files, "\n".join([
                    "Exactly one structure file (not zero, not multiple) must be supplied per command/string.",
                    "Below are examples of how to import multiple molecules.",
                    "Within python:",
                    "    "+"molecule_import = [",
                    "    "+"    "+"\"file:molecule1.pdb moleculetype:MOL1\",",
                    "    "+"    "+"\"file:molecule2.pdb name:MOL2 charge:2\",",
                    "    "+"    "+"\"file:molecule3.pdb name:MOL3 charge:res:2:bead:3:charge:-1\",",
                    "    "+"]",
                    "Command line:",
                    "    "+"-molecule_import file:mol1.pdb moleculetype:MOL1 -molecule_import file:mol2.pdb name:MOL2 charge:2 -molecule_import file:mol3.pdb name:MOL2 charge:res:2:bead:3:charge:-1",
                ])

                nmoleculetypes = [1 for i in imp_struc.split() if (i.startswith("moleculetype:") or i.startswith("moleculetypes:"))]
                nnames         = [1 for i in imp_struc.split() if i.startswith("name:")]
                assert not (len(nmoleculetypes) > 0 and len(nnames) > 0), "\n".join([
                    "Both 'moleculetype(s)' and 'name' subarguments have been used. You may only use one of the two.",
                    "Below are examples of correct commands. Note that charge(s) must be explicitly defined when using 'name'.",
                    "    "+"file:molecule1.pdb moleculetype:MOL1",
                    "    "+"file:molecule2.pdb name:MOL2 charge:2",
                    "    "+"file:molecule3.pdb name:MOL3 charge:res:2:bead:3:charge:-1",
                ])

                if len(nmoleculetypes) > 0:
                    assert len(nmoleculetypes) == 1, "Only one 'moleculetype' may be given per command/string"
                if len(nnames) > 0:
                    assert len(nnames) == 1, "Only one 'name' may be given per command/string"

                for cmd in imp_struc.split():
                    sub_cmd = cmd.split(":")
                    ### Read pdb/gro file
                    if sub_cmd[0] == "file":

                        ### Find file name
                        file_name = sub_cmd[1]

                        structure = self.structure_importer(file_name)
                    
                    elif sub_cmd[0] == "params":
                        params = sub_cmd[1]
                    
                    elif sub_cmd[0] in ["moleculetype", "moleculetypes"]:
                        moleculetype = sub_cmd[1]
                        
                    elif sub_cmd[0] in ["refname"]:
                        refname = sub_cmd[1]
                        
                    elif sub_cmd[0] == "name":
                        name = sub_cmd[1]
                    
                    elif sub_cmd[0] == "charge":
                        charge_cmds.append(sub_cmd)
                    
                    elif sub_cmd[0] in ["upbead", "upbeads"]:
                        upbeads_cmds.append(sub_cmd[1:])
                    
                    elif sub_cmd[0] in ["downbead", "downbeads"]:
                        downbeads_cmds.append(sub_cmd[1:])
                    
                    elif sub_cmd[0] in ["library_type", "library_types"]:
                        assert all([string in ["solvent", "solute", "ions", "neg_ions", "pos_ions", "lipid"] for string in sub_cmd[1:]]), "\n".join([
                            "Incorrect value(s) given to 'library_types' subargument: " + str(sub_cmd[1:]),
                            "Valid values are:",
                            "    "+"'solvent'/'solute' (same library)",
                            "    "+"'neg_ions'",
                            "    "+"'pos_ions'",
                            "    "+"'ions' (counts as both 'neg_ions' and 'pos_ions')",
                            "    "+"'lipid'",
                        ])
                        library_types.extend(sub_cmd[1:])
                    
                    ### Solvent-specific values
                    elif sub_cmd[0] == "molar_mass":
                        val = self.get_number_from_string(sub_cmd[1])
                        assert val is not False, "Value given to 'molar_mass' must be a number"
                        molar_mass = val
                    
                    elif sub_cmd[0] == "density":
                        val = self.get_number_from_string(sub_cmd[1])
                        assert val is not False, "Value given to 'density' must be a number"
                        density = val
                    
                    elif sub_cmd[0] == "mapping_ratio":
                        val = self.get_number_from_string(sub_cmd[1])
                        assert val is not False, "Value given to 'mapping_ratio' must be a number"
                        mapping_ratio = val
                    
                    ### Used to scale imported x/y/z values
                    elif sub_cmd[0] == "scale":
                        curr_axis = "all"
                        for val in sub_cmd[1:]:
                            subval = self.get_number_from_string(val)
                            if subval is False:
                                assert val in ["all", "x", "y", "z"], "Values given to sub-argument 'scale' must be either the axis ('all', 'x', 'y', 'z') or the value to scale the axis by."
                                curr_axis = val
                            else:
                                if curr_axis == "all":
                                    scalings["x"], scalings["y"], scalings["z"] = subval, subval, subval
                                elif curr_axis in ["x", "y", "z"]:
                                    scalings[curr_axis] = subval
                    
                    else:
                        assert False, "Unknown subargument given to molecule_import: '{sub_cmd}'".format(sub_cmd=sub_cmd)
                    
                assert len(structure) > 0, "No structure file given (.pdb, .gro and .cif are accepted file formats)"
                nres                             = len(set(list(zip(*structure.keys()))[2]))

                charge                           = False
                charges                          = []
                check_no_charge                  = len(charge_cmds) == 0
                check_single_unversal_charge     = len(charge_cmds) == 1 and len(charge_cmds[0]) == 2
                check_one_res_any_charges_method = len(charge_cmds) > 0 and nres == 1 and all([len(sub_cmd) in [4, 6] for sub_cmd in charge_cmds])
                check_mult_res_only_res_charges  = len(charge_cmds) > 0 and nres >= 1 and all([len(sub_cmd) == 6 for sub_cmd in charge_cmds])

                self.print_term("moleculetype and name:",                   moleculetype, name,               debug=True)
                self.print_term("charge commands:      ",                   charge_cmds,                      debug=True)
                self.print_term("    ", "charge                          ", charge,                           debug=True)
                self.print_term("    ", "charges                         ", charges,                          debug=True)
                self.print_term("    ", "check_no_charge                 ", check_no_charge,                  debug=True)
                self.print_term("    ", "check_single_unversal_charge    ", check_single_unversal_charge,     debug=True)
                self.print_term("    ", "check_one_res_any_charges_method", check_one_res_any_charges_method, debug=True)
                self.print_term("    ", "check_mult_res_only_res_charges ", check_mult_res_only_res_charges,  debug=True)

                assert any([check_no_charge, check_single_unversal_charge, check_one_res_any_charges_method, check_mult_res_only_res_charges]), "\n".join([
                    "Improper charge subarguments have been given. The following charge subarguments are accepted.",
                    "Always usable:",
                    "    " + "No charge subargument",
                    "    " + "charge:[float/int] (a maximum of one charge command is permitted if this one is used)",
                    "    " + "charge:res:[resnr]:bead:[beadnr]:charge:[float/int] (any number of these may be used)",
                    "Only usable for single-residue molecules:",
                    "    " + "charge:bead:[beadnr]:charge:[float/int] (any number of these may be used)",
                ])

                for sub_cmd in charge_cmds:
                    if len(sub_cmd) == 2:
                        charge = ast.literal_eval(sub_cmd[1])
                    else:
                        cur_res = False
                        cur_bead = False
                        i = 0
                        while i < len(sub_cmd):
                            if sub_cmd[i] in ["res", "resnr"]:
                                cur_res = ast.literal_eval(sub_cmd[i+1])
                            elif sub_cmd[i] in ["bead", "beadnr"]:
                                cur_bead = ast.literal_eval(sub_cmd[i+1])
                            elif sub_cmd[i] == "charge":
                                cur_charge = ast.literal_eval(sub_cmd[i+1])
                            ### Always +2
                            i += 2
                        if cur_res is False:
                            cur_res = 0
                        charges.append((cur_res, cur_bead, cur_charge))
                        ### Reset residue and bead designation. Charge designation is pointless to redesignate.
                        cur_res    = False
                        cur_bead   = False
                        cur_charge = False
                
                assert (len(upbeads_cmds) == 0 and len(downbeads_cmds) == 0) or (len(upbeads_cmds) > 0 and len(downbeads_cmds) > 0), "\n".join([
                    "If you wish to designate the top-most and bottom-most beads then you need to supply at least one bead for each.",
                    "Top-most beads:    " + str(upbeads_cmds),
                    "Bottom-most beads: " + str(downbeads_cmds),
                ])

                upbeads   = []
                downbeads = []
                if len(upbeads_cmds) > 0 and len(downbeads_cmds) > 0:
                    check_one_res_any_method_upbead = nres == 1 and all([len(sub_cmd) in [2, 4] for sub_cmd in upbeads_cmds])
                    check_mult_res_only_res_upbead  = nres >= 1 and all([len(sub_cmd) == 4 for sub_cmd in upbeads_cmds])
                    assert any([check_one_res_any_method_upbead, check_mult_res_only_res_upbead,]), "\n".join([
                        "Improper upbead subarguments have been given. The following subarguments are accepted.",
                        "Your subarguments: " + str(upbeads_cmds),
                        "Always usable:",
                        "    " + "No upbead subargument",
                        "    " + "upbead:res:[resnr]:bead:[beadnr] (any number of these may be used)",
                        "Only usable for single-residue molecules:",
                        "    " + "upbead:bead:[beadnr] (any number of these may be used)",
                    ])
                    check_one_res_any_method_downbead = nres == 1 and all([len(sub_cmd) in [2, 4] for sub_cmd in downbeads_cmds])
                    check_mult_res_only_res_downbead  = nres >= 1 and all([len(sub_cmd) == 4 for sub_cmd in downbeads_cmds])
                    assert any([check_one_res_any_method_downbead, check_mult_res_only_res_downbead,]), "\n".join([
                        "Improper downbead subarguments have been given. The following subarguments are accepted.",
                        "Your subarguments: " + str(downbeads_cmds),
                        "Always usable:",
                        "    " + "No downbead subargument",
                        "    " + "downbead:res:[resnr]:bead:[beadnr] (any number of these may be used)",
                        "Only usable for single-residue molecules:",
                        "    " + "downbead:bead:[beadnr] (any number of these may be used)",
                    ])

                    for beads_cmds, beads in [(upbeads_cmds, upbeads), (downbeads_cmds, downbeads)]:
                        for sub_cmd in beads_cmds:
                            cur_res = False
                            cur_bead = False
                            i = 0
                            while i < len(sub_cmd):
                                if sub_cmd[i] in ["res", "resnr"]:
                                    cur_res = ast.literal_eval(sub_cmd[i+1])
                                elif sub_cmd[i] in ["bead", "beadnr"]:
                                    cur_bead = ast.literal_eval(sub_cmd[i+1])
                                ### Always +2
                                i += 2
                            
                            if cur_res is False:
                                cur_res = 0
                            assert cur_bead is not False, "If using 'upbead' or 'downbead', then the bead nr must be designated using 'bead:[int]"
                            beads.append((cur_res, cur_bead))
                            
                            ### Reset residue and bead designation. Charge designation is pointless to redesignate.
                            cur_res = False
                            cur_bead = False

                self.print_term("upbeads", upbeads, debug=True)
                self.print_term("downbeads", downbeads, debug=True)

                assert (moleculetype is not False) or (name is not False), "\n".join([
                    "Either 'moleculetype' or 'name' must be given.",
                    "If 'name' is given then 'charge/charges' should also be given, otherwise molecule charge is assumed to be zero."
                ])

                mol_dict = {
                    "residues": [],
                }

                ### Charge information extracted later.
                ### "charge" must be either number or string
                ### "solute_defs_checker" checks if "charge" is string and interprets string as moleculetype
                ### "charges" must be (resnr, beadnr, chargeval)
                if moleculetype is not False:
                    mol_dict["moleculetype"] = moleculetype

                if charge is not False:
                    mol_dict["charge"] = charge
                elif len(charges) > 0:
                    mol_dict["charges"] = charges
                else:
                    mol_dict["charge"] = 0

                if len(upbeads) > 0 and len(downbeads) > 0:
                    mol_dict["upbeads"]   = upbeads
                    mol_dict["downbeads"] = downbeads

                ### Solvent-specific values
                if molar_mass is not False:
                    mol_dict["molar_mass"] = molar_mass
                if density is not False:
                    mol_dict["density"] = density
                if mapping_ratio is not False:
                    mol_dict["mapping_ratio"] = mapping_ratio

                ### Finds all the residues
                residues = []
                last_resnr = False
                for (index, beadnr, resnr), values in structure.items():
                    if (last_resnr is False) or (last_resnr is not False and resnr != last_resnr):
                        last_resnr = resnr
                        residues.append({
                            "resname": values["res_name"],
                            "names":   [],
                            "x":       [],
                            "y":       [],
                            "z":       [],
                        })
                    residues[-1]["names"].append(values["atom_name"])
                    residues[-1]["x"].append(values["x"] / 10 * scalings["x"]) # Divided by 10 because "solute_beads_check" multiplies them by 10 due to library values needing it
                    residues[-1]["y"].append(values["y"] / 10 * scalings["y"]) # Divided by 10 because "solute_beads_check" multiplies them by 10 due to library values needing it
                    residues[-1]["z"].append(values["z"] / 10 * scalings["z"]) # Divided by 10 because "solute_beads_check" multiplies them by 10 due to library values needing it

                for residue in residues:
                    mol_dict["residues"].append({
                            "resname": residue["resname"],
                            "names":   tuple(residue["names"]),
                            "x":       tuple(residue["x"]),
                            "y":       tuple(residue["y"]),
                            "z":       tuple(residue["z"]),
                    })
                
                if moleculetype is not False:
                    if refname is not False:
                        defs_name = refname
                    else:
                        defs_name = moleculetype
                else:
                    defs_name = name

                ### Solvent
                if len(library_types) == 0 or any([string in library_types for string in ["solvent", "solute"]]):
                    if params not in self.solvent_defs_imported.keys():
                        self.solvent_defs_imported[params] = {}
                    self.solvent_defs_imported[params][defs_name] = copy.deepcopy(mol_dict)
                
                ### Positive ions
                if len(library_types) == 0 or any([string in library_types for string in ["ions", "pos_ions"]]):
                    if params not in self.pos_ion_defs_imported.keys():
                        self.pos_ion_defs_imported[params] = {}
                    self.pos_ion_defs_imported[params][defs_name] = copy.deepcopy(mol_dict)
                
                ### Negative ions
                if len(library_types) == 0 or any([string in library_types for string in ["ions", "neg_ions"]]):
                    if params not in self.neg_ion_defs_imported.keys():
                        self.neg_ion_defs_imported[params] = {}
                    self.neg_ion_defs_imported[params][defs_name] = copy.deepcopy(mol_dict)
                
                ### Lipids
                if len(library_types) == 0 or any([string in library_types for string in ["lipid"]]):
                    if params not in self.lipid_defs_imported.keys():
                        self.lipid_defs_imported[params] = {}
                    self.lipid_defs_imported[params][defs_name] = copy.deepcopy(mol_dict)
