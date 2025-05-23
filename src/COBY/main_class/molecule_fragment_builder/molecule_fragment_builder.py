import copy

class molecule_fragment_builder:
    def molecule_fragment_builder(self):
        
        if self.MOLECULE_FRAGMENT_BUILDER_cmds:
            self.print_term("Building molecules from fragments", verbose=2)
            for cmdi, cmd in enumerate(self.MOLECULE_FRAGMENT_BUILDER_cmds, 1):
                cmd = cmd.lstrip().rstrip()
                cmd_dict = {}
                if ":" in cmd:
                    for cmd_part in cmd.split():
                        part, val = cmd_part.split(":", maxsplit=1)
                        assert part not in cmd_dict.keys(), "Molecule part '{part}' given multiple times for molecule builder argument '{cmd}'".format(part=part, cmd=cmd)
                        cmd_dict[part] = val
                
                self.print_term("MOLECULE_FRAGMENT_BUILDER cmd:     ", cmd,      debug=True, debug_keys=["MFB"])
                self.print_term("MOLECULE_FRAGMENT_BUILDER cmd_dict:", cmd_dict, debug=True, debug_keys=["MFB"])

                assert "moltype" in cmd_dict.keys(), "'moltype' key:value pair must be designated for molecule fragment builder"
                assert cmd_dict["moltype"] in self.fragment_defs.keys(), "\n".join([
                    "moltype '{moltype}' not found in fragment library".format(moltype=cmd_dict["moltype"]),
                    "Currently allowed molecule types are:",
                    "    "+" ".join(list(self.fragment_defs.keys()))
                ])

                ### Sets moltype
                moltype = cmd_dict.pop("moltype")

                ### Sets reference name and resname to be the same
                assert "name" in cmd_dict.keys(), "'name' key:value pair must be designated for molecule fragment builder"
                name = cmd_dict.pop("name")
                
                ### If resname explicitly given, then overwrite
                resname = name
                if resname in cmd_dict.keys():
                    resname = cmd_dict.pop("resname")
                
                ### Sets moleculetype (for .itps) if given. Does not add if not explicitly given.
                moleculetype = False
                if "moleculetype" in cmd_dict.keys():
                    moleculetype = cmd_dict.pop("moleculetype")
                
                ### Sets default params to be "MFB"

                ### Disables adding the molecule to the list of MoleculeBuilder arguments if it is obtained from topology.
                ### Params subargument "params:TOP" is automatically added to argument in 'itp_reader.py' method if 'params' is missing so no need to check for it here.
                FromTopology = False
                if "FromTopology" in cmd_dict.keys():
                    FromTopology = True
                    del cmd_dict["FromTopology"]

                ### Overwrites params if given
                params = "MFB" # Molecule Fragment Builder
                if "params" in cmd_dict.keys():
                    params = cmd_dict.pop("params")

                if "default_parts" in self.fragment_defs[moltype].keys():
                    for part, val in self.fragment_defs[moltype]["default_parts"].items():
                        if part not in cmd_dict.keys():
                            cmd_dict[part] = val

                ### Adds built lipids to molecule builder argument list if currently running "MoleculeBuilder" class and lipid not obtained from topology (checked with FromTopology)
                if self.PROGRAM == "Crafter" and not FromTopology:
                    self.MOLECULE_cmds.append("moltype:lipid params:{params} name:{name}".format(params=params, name=name))
                
                ### Constructs lipid based on given part:fragment pairs
                part_dicts = {}
                for part, val in cmd_dict.items():
                    part_dicts[part] = {}
                    part_dicts[part]["connected_parts"] = []

                    ### Checks if part is among accepted parts
                    assert part in self.fragment_defs[moltype]["accepted_parts"], "\n".join([
                        "Part '{part}' specified for lipid builder argument '{cmd}' not found in list of accepted parts".format(part=part, cmd=cmd),
                        "Currently accepted parts for moltype '{moltype}' are:".format(moltype=moltype),
                        "    "+" ".join(self.fragment_defs[moltype]["accepted_parts"]),
                    ])

                    ### Checks if part is in the part library for the lipid type
                    assert part in self.fragment_defs[moltype]["parts"].keys(), "\n".join([
                        "Part '{part}' specified for lipid builder argument '{cmd}' not found in list of parts".format(part=part, cmd=cmd),
                        "Currently available parts for moltype '{moltype}' are:".format(moltype=moltype),
                        "    "+" ".join(list(self.fragment_defs[moltype]["parts"].keys())),
                    ])

                    ### Checks if "function" is defined for part in which case the part should be dynamically built from an argument string
                    if "function" in self.fragment_defs[moltype]["parts"][part].keys():
                        if "kwargs" in self.fragment_defs[moltype]["parts"][part].keys():
                            kwargs = self.fragment_defs[moltype]["parts"][part]["kwargs"]
                            beads = self.fragment_defs[moltype]["parts"][part]["function"](val, **kwargs)
                        else:
                            beads = self.fragment_defs[moltype]["parts"][part]["function"](val)
                        
                        # assert all([len(bead) in [2, 3] for bead in beads]), "\n".join([
                        assert all(key in bead.keys() for key in ["x", "y", "z"] for bead in beads) and any(key in bead.keys() for key in ["beadname", "name"] for bead in beads), "\n".join([
                            "Beads returned from a molecule part builder function must be written as a list of dictionaries.",
                            "    "+"The beads were returned as:",
                            "    "+"    "+str(beads),
                            "    "+"The dictionaries MUST contain the following key:value pairs:",
                            "    "+"    "+"'beadname' or 'name'",
                            "    "+"    "+"'x'",
                            "    "+"    "+"'y'",
                            "    "+"    "+"'z'",
                            "    "+"The dictionaries MAY contain the following key:value pairs:",
                            "    "+"    "+"'charge' (assumed zero if not given)",
                            "    "+"    "+"'resname' (taken from molecule name/resname subargument)",
                            "    "+"    "+"'resnr' (assumed to be 0 (python indexing starts at 0))",
                        ])

                        part_dicts[part]["beaddata"] = beads.copy()

                        if "join_to" in self.fragment_defs[moltype]["parts"][part].keys():
                            part_dicts[part]["join_to"] = copy.deepcopy(self.fragment_defs[moltype]["parts"][part]["join_to"])

                        if "join_from" in self.fragment_defs[moltype]["parts"][part].keys():
                            part_dicts[part]["join_from"] = copy.deepcopy(self.fragment_defs[moltype]["parts"][part]["join_from"])

                        
                    else:
                        assert val in self.fragment_defs[moltype]["parts"][part].keys(), "\n".join([
                            "Fragment '{fragment}' specified for part '{part}' for molecule builder argument '{cmd}' not found in list of fragments".format(fragment=val, part=part, cmd=cmd),
                            "Currently available fragments for part '{part}' for moltype '{moltype}' are:".format(part=part, moltype=moltype),
                            "    "+" ".join(list(self.fragment_defs[moltype]["parts"][part].keys())),
                        ])

                        assert "beads" in self.fragment_defs[moltype]["parts"][part][val].keys(), "\n".join([
                            "'beads' key not found in dictionary for fragment '{fragment}' for part '{part}' for molecule builder argument '{cmd}'".format(fragment=val, part=part, cmd=cmd),
                            "Current keys for fragment '{fragment}' for part '{part}' for molecule builder argument '{cmd}' are:".format(fragment=val, part=part, cmd=cmd),
                            "    "+" ".join(list(self.fragment_defs[moltype]["parts"][part][val].keys())),
                        ])
                        
                        part_dicts[part]["beaddata"] = copy.deepcopy(self.fragment_defs[moltype]["parts"][part][val]["beads"])

                        if "join_to" in self.fragment_defs[moltype]["parts"][part][val].keys():
                            part_dicts[part]["join_to"] = copy.deepcopy(self.fragment_defs[moltype]["parts"][part][val]["join_to"])

                        if "join_from" in self.fragment_defs[moltype]["parts"][part][val].keys():
                            part_dicts[part]["join_from"] = copy.deepcopy(self.fragment_defs[moltype]["parts"][part][val]["join_from"])

                for part_key in part_dicts.keys():
                    part_dicts[part_key].update({
                        "beadnames": [],
                        "charges":   [],
                        "coords":    [],
                        "resnames":  [],
                        "resnrs":    [],
                    })
                    for bi, bead in enumerate(part_dicts[part_key]["beaddata"]):
                        if "name" in bead:
                            part_dicts[part_key]["beadnames"].append(bead["name"])
                        elif "beadname" in bead:
                            part_dicts[part_key]["beadnames"].append(bead["beadname"])

                        part_dicts[part_key]["coords"].append((bead["x"], bead["y"], bead["z"]))

                        charge = (bi, 0)
                        if "charge" in bead:
                            charge = (bi, bead["charge"])
                        part_dicts[part_key]["charges"].append(charge)

                        bead_resname = resname # resname defined further up
                        if "resname" in bead and bead["resname"] != False:
                            bead_resname = bead["resname"]
                        part_dicts[part_key]["resnames"].append(bead_resname)

                        resnr = 0
                        if "resnr" in bead:
                            resnr= bead["resnr"]
                        part_dicts[part_key]["resnrs"].append(resnr)



                ### Sorting dict based on order needed for correct bead order
                part_dicts = {key:part_dicts[key] for key in self.fragment_defs[moltype]["order"] if key in part_dicts.keys()}

                ### Recursively adjusts coordinates and join_to/join_from values for connected parts
                def connected_coord_adjuster(part, connector_difference):
                    for bi, coords in enumerate(part_dicts[part]["coords"]):
                        coords_new = tuple([round(coord - joiner, 3) for coord, joiner in zip(coords, connector_difference)])
                        part_dicts[part]["coords"][bi] = coords_new
                    
                    part_dicts[part]["join_to"] = (join_to_part, tuple([round(coord - joiner, 3) for coord, joiner in zip(part_dicts[part]["join_to"][1], connector_difference)]))

                    if "join_from" in part_dicts[part].keys():
                        part_dicts[part]["join_from"] = {
                            join_from_part: tuple([round(coord - joiner, 3) for coord, joiner in zip(join_from_position, connector_difference)])
                            for join_from_part, join_from_position in part_dicts[part]["join_from"].items()
                        }
                    
                    for connected_part in part_dicts[part]["connected_parts"]:
                        connected_coord_adjuster(connected_part, connector_difference)


                ### Connects the different parts of the molecule
                for part, part_dict in part_dicts.items():
                    if "join_to" in part_dict:
                        join_to_part, join_to_position = part_dict["join_to"]
                        
                        connector_difference = tuple([
                            round(coord1 - coord2, 3) 
                            for coord1, coord2 
                            in zip(
                                join_to_position,
                                part_dicts[join_to_part]["join_from"][part],
                            )
                        ])

                        part_dicts[join_to_part]["connected_parts"].append(part)

                        connected_coord_adjuster(part, connector_difference)

                joined_beadnames    = []
                joined_charges  = []
                joined_coords   = []
                joined_resnames = []
                joined_resnrs   = []

                for part, part_dict in part_dicts.items():
                    ### Calculate charges first as it uses the length of "joined_coords" from previous parts but should not include current part
                    if "charges" in part_dict.keys() and part_dict["charges"]:
                        charges = []
                        if type(part_dict["charges"]) in (tuple, list) and type(part_dict["charges"][0]) in (tuple, list):
                            charges.extend(list([tuple(charge) for charge in part_dict["charges"]]))
                        else:
                            charges.append(tuple(part_dict["charges"]))
                        
                        for charge_i, charge_value in charges:
                            charge_i_new = charge_i + len(joined_coords)
                            joined_charges.append((charge_i_new, charge_value))
                    
                    joined_beadnames.extend(part_dict["beadnames"])
                    joined_coords.extend(part_dict["coords"])
                    joined_resnames.extend(part_dict["resnames"])
                    joined_resnrs.extend(part_dict["resnrs"])
                
                joined_charges = tuple(sorted(joined_charges, key=lambda x: x[0]))
                joined_resnrs = [i-min(joined_resnrs) for i in joined_resnrs] # ensures that residue numbers start at 0

                xs, ys, zs = zip(*joined_coords)

                residues = []
                cur_resnr = -1
                for bi, (beadname, charge, x, y, z, resname, resnr) in enumerate(zip(joined_beadnames, joined_charges, xs, ys, zs, joined_resnames, joined_resnrs)):
                    if resnr != cur_resnr:
                        residues.append({
                            "resname": resname,
                        })
                        if bi != 0:
                            residues[cur_resnr].update({
                                "names":   tuple(res_beads),
                                "x":       tuple(res_x),
                                "y":       tuple(res_y),
                                "z":       tuple(res_z),
                                "charges": tuple(res_charges),
                            })
                        res_beads   = []
                        res_x       = []
                        res_y       = []
                        res_z       = []
                        res_charges = []
                        cur_resnr += 1

                    res_beads.append(beadname)
                    res_x.append(x)
                    res_y.append(y)
                    res_z.append(z)
                    res_charges.append(charge)
                
                residues[cur_resnr].update({
                    "names":   tuple(res_beads),
                    "x":       tuple(res_x),
                    "y":       tuple(res_y),
                    "z":       tuple(res_z),
                    "charges": tuple(res_charges),
                })
                
                molecule = {
                    name: {
                        "residues": residues.copy(),
                    },
                }
                
                if moleculetype:
                    molecule[name]["moleculetype"] = moleculetype

                if params not in self.lipid_defs_built.keys():
                    self.lipid_defs_built[params] = {}
                self.lipid_defs_built[params].update(copy.deepcopy(molecule))
                
                if params not in self.solvent_defs_built.keys():
                    self.solvent_defs_built[params] = {}
                self.solvent_defs_built[params].update(copy.deepcopy(molecule))
                
                if params not in self.pos_ion_defs_built.keys():
                    self.pos_ion_defs_built[params] = {}
                self.pos_ion_defs_built[params].update(copy.deepcopy(molecule))
                
                if params not in self.neg_ion_defs_built.keys():
                    self.neg_ion_defs_built[params] = {}
                self.neg_ion_defs_built[params].update(copy.deepcopy(molecule))

            self.print_term("Number of molecules built:", cmdi, "\n", spaces=1, verbose=2)