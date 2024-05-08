import ast
import copy
import numpy as np

class memb_preprocessor:
    def memb_preprocessor(self, return_self = "self", cmds_given = False, extra_verbose=0):
        '''
        Preprocesses membrane commands for later ease of use
        '''
        ### ### Preprocessing membranes
        
        if return_self == "return":
            MEMBRANES = {}
        
        if cmds_given:
            MEMBRANES_cmds = cmds_given
        else:
            MEMBRANES_cmds = self.MEMBRANES_cmds
            
        if len(MEMBRANES_cmds) != 0:
            self.print_term("Preprocessing membrane requests", verbose=2+extra_verbose)
            for cmd_nr, memb_cmd in enumerate(MEMBRANES_cmds, 1):
                self.print_term("Starting membrane command:", cmd_nr, spaces=1, verbose=2+extra_verbose)
                ### "membrane" settings always apply to the whole membrane instead of individual leaflets
                ### "default" settings can be overwritten by individual leaflets
                
                settings_dict = {
                    "upper_leaf": {},
                    "lower_leaf": {},
                    "membrane"  : {
                        "xlength": self.pbcx / 10, # [nm] converted to [Å]
                        "ylength": self.pbcy / 10, # [nm] converted to [Å]
                        "center": [0, 0, 0], # [nm] converted to [Å]
                        
                        ### Two-value choices: ("auto", [float/int]) or ([float/int/False], [float/int/False])
                        ### One-value choices [float/int/False]
                        "gridsplits": ("auto", 500), # Å 
                        
                        "pbc_check": True, # [bool]
                        
                        "grid_maker": "lines", # "lines" or "3D_matrix"
                        "lipid_distribution": "random",
                        
                        ### Different optimize version removed in COBY version 0.0.8.
                        ### "optimize" value kept as True/False settings
                        # "optimize": "v5", # False/"no"/0, "limited", True/"full"/1
                        "optimize": "auto", # False/"no"/0, "limited", True/"full"/1
                        "optim_maxsteps": 100,
                        "optim_push_tol": 0.5,
                        "optim_push_mult": 1.0,
                        
                        ### Readjusts the max/min x/y-values if cutouts/holes are made
                        "readjust_bbox": True,
                        ### Readjusts the max/min x/y-values if cutouts/holes are made
                        "split_bbox": True,
                        ### Allows empty space (holes/cutouts) to be solvated
                        "solvate_hole": True,
                        ### Used for nanodiscs. Sets whether the resulting membrane should be inside or outside of a protein.
                        "inside_protein": False,
                    },
                    "default": {
                        "lipids": {}, # empty
                        "lipids_preprocessing": [], # empty

                        "kickxy": 0.025, # [nm] converted to [Å]
                        "kickz": 0.025, # [nm] converted to [Å]

                        "apl": 0.6, # [nm^2] converted to [Å^2]
                        
                        ### 0.264 = vdw of regular beads. 0.264/16*6 = 0.099, 0.264/4=0.066, # 0.264/2=0.132
                        "plane_buffer": 0.264/16*6, # [nm] converted to [Å]
                        "height_buffer": 0.264/16*6, # [nm] converted to [Å]

                        "prot_buffer": 0.132, # [nm] converted to [Å], default = (vdw of regular beads) / 2
                        "alpha_mult": 1.0,
                        
                        "lip_round_func": ("int", int), # int, round, math.floor, math.ceil

                        "lipid_optim": 'force_fill', # str: see beloww
                        ### 'force_fill':  Fills up to the allowed number number of lipids
                        ### 'fill:         Same as 'force_fill' but stops if perfect lipid ratio has been achieved
                        ### 'avg_optimal': Ratios are optimised based on average lipid distance from optimal composition
                        ### 'no':          No optimization.
                        ### 'abs_val':     Ratio-values are treated as absolute number of lipids
                        ### 'insane':      'no' and initial lipid calculation is intentionally wrong like in insane
                        "params": False, # False or str
                        "charge": "top", # "lib" or "top"
                        
                        "hole": [],
                        "patch": [],
                    }
                }
                hole_types = ["circle", "ellipse", "square", "rectangle", "polygon"]
                holes_patches_abbrevs_general_key_vals = [
                    (("rotate"), "rotate"),
                    (("buffer"), "buffer"),
                    (("buffer_cap"), "buffer_cap"),
                    (("buffer_join"), "buffer_join"),
                ]
                holes_patches_abbrevs_polygon_key_vals = [
                    (("p", "point"), "point"),

                    (("scaling"), "scaling"),
                    (("xscaling"), "xscaling"),
                    (("yscaling"), "yscaling"),
                ]
                holes_patches_abbrevs_nonpolygon_key_vals = [
                    (("cx"), "cx"),
                    (("cy"), "cy"),
                ]                
                holes_patches_abbrevs_circle_key_vals = [
                    (("radius"), "radius"),
                ]
                holes_patches_abbrevs_ellipse_key_vals = [
                    (("xradius"), "xradius"),
                    (("yradius"), "yradius"),
                ]
                holes_patches_abbrevs_square_key_vals = [
                    (("length"), "length"),
                ]
                holes_patches_abbrevs_rectangle_key_vals = [
                    (("xlength"), "xlength"),
                    (("ylength"), "ylength"),
                ]

                def abbrevs_keys_vals_converter(keys_vals):
                    outdict = {}
                    for keys, val in keys_vals:
                        if type(keys) == str:
                            keys = (keys,)
                        outdict.update(dict.fromkeys(keys, val))
                    return outdict

                holes_patches_abbrevs_general    = abbrevs_keys_vals_converter(holes_patches_abbrevs_general_key_vals)
                holes_patches_abbrevs_polygon    = abbrevs_keys_vals_converter(holes_patches_abbrevs_polygon_key_vals)
                holes_patches_abbrevs_nonpolygon = abbrevs_keys_vals_converter(holes_patches_abbrevs_nonpolygon_key_vals)
                holes_patches_abbrevs_circle     = abbrevs_keys_vals_converter(holes_patches_abbrevs_circle_key_vals)
                holes_patches_abbrevs_ellipse    = abbrevs_keys_vals_converter(holes_patches_abbrevs_ellipse_key_vals)
                holes_patches_abbrevs_square     = abbrevs_keys_vals_converter(holes_patches_abbrevs_square_key_vals)
                holes_patches_abbrevs_rectangle  = abbrevs_keys_vals_converter(holes_patches_abbrevs_rectangle_key_vals)

                def hp_general_abbrevs_func(self, key, sub_cmd, i):
                    if key in ["rotate", "buffer"]:
                        nvals = 2
                        keyvals = [(key, self.get_number_from_string(sub_cmd[i+1]))]
                    elif key == "buffer_cap":
                        nvals = 2
                        buffer_cap = sub_cmd[i+1]
                        if buffer_cap in [1, "1", "round"]:
                            keyvals = [(key, 1)]
                        if buffer_cap in [2, "2", "flat"]:
                            keyvals = [(key, 2)]
                        if buffer_cap in [3, "3", "square"]:
                            keyvals = [(key, 3)]
                    elif key == "buffer_join":
                        nvals = 2
                        buffer_join = sub_cmd[i+1]
                        if buffer_join in [1, "1", "round"]:
                            keyvals = [(key, 1)]
                        if buffer_join in [2, "2", "mitre"]:
                            keyvals = [(key, 2)]
                        if buffer_join in [3, "3", "bevel"]:
                            keyvals = [(key, 3)]
                        keyvals = [(key, buffer_join)]
                    return nvals, keyvals

                def hp_polygon_abbrevs_func(self, key, sub_cmd, i):
                    ### Point designations
                    if key in ["point"]:
                        nvals = 3
                        point = (self.get_number_from_string(sub_cmd[i+1]), self.get_number_from_string(sub_cmd[i+2]))
                        keyvals = [(key, point)]
                    
                    ### Scalings
                    elif key in ["scaling"]:
                        nvals = 2
                        keyvals = [
                            ("xscaling", self.get_number_from_string(sub_cmd[i+1])),
                            ("yscaling", self.get_number_from_string(sub_cmd[i+1])),
                        ]
                    elif key in ["xscaling", "yscaling"]:
                        nvals = 2
                        keyvals = [(key[0] + "scaling", self.get_number_from_string(sub_cmd[i+1]))]
                    return nvals, keyvals

                def hp_nonpolygon_abbrevs_func(self, key, sub_cmd, i):
                    ### Center designation
                    if key in ["cx"]:
                        nvals = 2
                        keyvals = [("cx", self.get_number_from_string(sub_cmd[i+1]))]
                    elif key in ["cy"]:
                        nvals = 2
                        keyvals = [("cy", self.get_number_from_string(sub_cmd[i+1]))]
                    return nvals, keyvals

                def hp_circle_abbrevs_func(self, key, sub_cmd, i):
                    if key in ["radius"]:
                        nvals = 2
                        keyvals = [
                            ("xscaling", self.get_number_from_string(sub_cmd[i+1])),
                            ("yscaling", self.get_number_from_string(sub_cmd[i+1])),
                        ]
                    return nvals, keyvals
                
                def hp_ellipse_abbrevs_func(self, key, sub_cmd, i):
                    if key in ["xradius", "yradius"]:
                        nvals = 2
                        keyvals = [(key[0] + "scaling", self.get_number_from_string(sub_cmd[i+1]))]
                    return nvals, keyvals
                    
                def hp_square_abbrevs_func(self, key, sub_cmd, i):
                    if key == "length":
                        nvals = 2
                        keyvals = [
                            ("xscaling", self.get_number_from_string(sub_cmd[i+1])/2),
                            ("yscaling", self.get_number_from_string(sub_cmd[i+1])/2),
                        ]
                    return nvals, keyvals
                
                def hp_rectangle_abbrevs_func(self, key, sub_cmd, i):
                    if key in ["xlength", "ylength"]:
                        nvals = 2
                        keyvals = [(key[0] + "scaling", self.get_number_from_string(sub_cmd[i+1])/2)]
                    return nvals, keyvals
                
                ### ### Membrane mono/bilayer names
                ### Mono as explicitly upper added by request
                monolayer_upper_designation = ["u", "up", "upper", "m", "mo", "mono", "monolayer", "mono_upper"]
                monolayer_lower_designation = ["d", "do", "down", "l", "lo", "lower", "mono_lower", "mono_down"]
                bilayer_designation = ["b", "bi", "bilayer", "memb", "membrane", "both"]
                
                layer_definition = "bilayer"
                ### Some values must always be for the membrane, as they make little sense to do per-leaflet
                dict_target      = "membrane"
                
                ### ### Check leaf command
                ### Liberal use of ast.literal_eval() below to convert strings to int/float
                ### ast.literal_eval checks if code is a valid python datatype before interpreting it
                for cmd in memb_cmd.split():
                    sub_cmd = cmd.split(":")
                    
                    ### Bilayer/monolayer defition
                    if sub_cmd[0].lower() in ["type"]:
                        if sub_cmd[1].lower() in bilayer_designation:
                            layer_definition = "bilayer"
                            dict_target      = "membrane"
                        elif sub_cmd[1].lower() in monolayer_upper_designation:
                            layer_definition = "upper"
                            dict_target      = "upper_leaf"
                        elif sub_cmd[1].lower() in monolayer_lower_designation:
                            layer_definition = "lower"
                            dict_target      = "lower_leaf"
                    
                    ### Sets whether following sub commands are for a specific leaflet or the whole membrane
                    elif sub_cmd[0].lower() in ["leaflet", "leaf", "side"]:
                        if sub_cmd[1].lower() in bilayer_designation:
                            dict_target = "membrane"
                        elif sub_cmd[1].lower() in monolayer_upper_designation:
                            dict_target = "upper_leaf"
                        elif sub_cmd[1].lower() in monolayer_lower_designation:
                            dict_target = "lower_leaf"

                    ### Leaflet shape
                    elif sub_cmd[0].lower() == "gridsplits":
                        ELSE = False
                        if len(sub_cmd[1:]) == 1:
                            val = sub_cmd[1]
                            isnumber, isint = self.is_number(val)
                            if isnumber:
                                if val == "0":
                                    val = "1"
                                settings_dict["membrane"]["gridsplits"] = (int(ast.literal_eval(val)), int(ast.literal_eval(val)))
                            elif val == "False":
                                settings_dict["membrane"]["gridsplits"] = (1, 1)
                            elif val.lower() == "auto":
                                settings_dict["membrane"]["gridsplits"] = ("auto", 500)
                            else:
                                ELSE = True
                        
                        elif len(sub_cmd[1:]) == 2:
                            val1, val2 = sub_cmd[1:]
                            isnumber1, isint1 = self.is_number(val1)
                            isnumber2, isint2 = self.is_number(val2)
                            if val1 == "auto" and isnumber2:
                                settings_dict["membrane"]["gridsplits"] = ("auto", ast.literal_eval(val2))
                            elif isnumber1 or isnumber2:
                                ### If values are numbers or "False"
                                if val1 in ["False", "0"]:
                                    val1 = "1"
                                if val2 in ["False", "0"]:
                                    val2 = "1"
                                settings_dict["membrane"]["gridsplits"] = (ast.literal_eval(val1), ast.literal_eval(val2))
                            else:
                                ELSE = True
                        else:
                            ELSE = True
                        if ELSE:
                            self.print_term(
                                "WARNING:",
                                "Subcommand:", "'" + sub_cmd[0] + "'", "was given with invalid value(s): ", "'" + sub_cmd[1:] + "'", "\n",
                                "For subcommand gridsplits:values; values must be one of: 'number/False', 'number1/False:number2/False', 'auto', 'auto:number'",
                                warn=True
                            )
                    
                    elif sub_cmd[0].lower() in ["hole", "patch"]:
                        assert sub_cmd[1].lower() in hole_types, "\n".join([
                            "You have used the subcommand (" + sub_cmd[0].lower() + ") but have specified the shape (" + sub_cmd[1].lower() + ")",
                            "which is not a valid shape for holes and patches",
                            "The valid shapes are: circle, ellipse, square, rectangle, polygon"
                        ])

                        hole_or_patch, shape_type, sub_cmd = sub_cmd[0].lower(), sub_cmd[1].lower(), sub_cmd[2:]
                        
                        if hole_or_patch not in settings_dict[dict_target]:
                            settings_dict[dict_target][hole_or_patch] = []
                        
                        settings_dict[dict_target][hole_or_patch].append({
                            "shape_type":  shape_type,

                            "points":      [],

                            "xscaling":    0, # xscaling=xradius for circle/ellipse and xlength for square/rectangle
                            "yscaling":    0, # yscaling=yradius for circle/ellipse and ylength for square/rectangle

                            "rotate":      0,

                            "buffer":      0,
                            "buffer_cap":  1, # 1=round, 2=flat,  3=square
                            "buffer_join": 1, # 1=round, 2=mitre, 3=bevel
                        })

                        cmd_vals = {}
                        if shape_type in ["circle", "ellipse"]:
                            if shape_type == "circle":
                                assert all([val in sub_cmd for val in ["cx", "cy", "radius"]]), "\n".join([
                                    "If you use the (circle) shape in a (" + hole_or_patch + ") subcommand,",
                                    "then all of the following settings must be specified: cx, cy and radius"
                                ])
                            elif shape_type == "ellipse":
                                assert all([val in sub_cmd for val in ["cx", "cy", "xradius", "yradius"]]), "\n".join([
                                    "If you use the (ellipse) shape in a (" + hole_or_patch + ") subcommand,",
                                    "then all of the following settings must be specified: cx, cy, xradius and yradius"
                                ])
                            i = 0
                            while i < len(sub_cmd):
                                key = sub_cmd[i]
                                ### General
                                if key in holes_patches_abbrevs_general:
                                    abbrev_key = holes_patches_abbrevs_general[sub_cmd[i]]
                                    nvals, keyvals = hp_general_abbrevs_func(self, abbrev_key, sub_cmd, i)

                                ### Center designation
                                elif key in holes_patches_abbrevs_nonpolygon:
                                    abbrev_key = holes_patches_abbrevs_nonpolygon[sub_cmd[i]]
                                    nvals, keyvals = hp_nonpolygon_abbrevs_func(self, abbrev_key, sub_cmd, i)

                                ### Circle specific
                                elif key in holes_patches_abbrevs_circle and shape_type == "circle":
                                    abbrev_key = holes_patches_abbrevs_circle[sub_cmd[i]]
                                    nvals, keyvals = hp_circle_abbrevs_func(self, abbrev_key, sub_cmd, i)

                                ### Ellipse specific
                                elif key in holes_patches_abbrevs_ellipse and shape_type == "ellipse":
                                    abbrev_key = holes_patches_abbrevs_ellipse[sub_cmd[i]]
                                    nvals, keyvals = hp_ellipse_abbrevs_func(self, abbrev_key, sub_cmd, i)
                                
                                ### Wrong subcommand given
                                else:
                                    assert False, "subcommand '" + key + "' not usable with circle/ellipse"
                                for subkey, subval in keyvals:
                                    cmd_vals[subkey] = subval
                                i += nvals
                            
                            settings_dict[dict_target][hole_or_patch][-1].update({
                                "points": [(cmd_vals["cx"], cmd_vals["cy"])],
                            })
                            del cmd_vals["cx"]
                            del cmd_vals["cy"]
                        
                        elif shape_type in ["square", "rectangle"]:
                            if shape_type == "square":
                                assert all([val in sub_cmd for val in ["cx", "cy", "length"]]), "\n".join([
                                    "If you use the (square) shape in a (" + hole_or_patch + ") subcommand,",
                                    "then all of the following settings must be specified: cx, cy and length"
                                ])
                            if shape_type == "rectangle":
                                assert all([val in sub_cmd for val in ["cx", "cy", "xlength", "ylength"]]), "\n".join([
                                    "If you use the (rectangle) shape in a (" + hole_or_patch + ") subcommand,",
                                    "then all of the following settings must be specified: cx, cy, xlength and ylength"
                                ])
                            i = 0
                            while i < len(sub_cmd):
                                key = sub_cmd[i]
                                ### General
                                if key in holes_patches_abbrevs_general:
                                    abbrev_key = holes_patches_abbrevs_general[sub_cmd[i]]
                                    nvals, keyvals = hp_general_abbrevs_func(self, abbrev_key, sub_cmd, i)

                                ### Center designation
                                elif key in holes_patches_abbrevs_nonpolygon:
                                    abbrev_key = holes_patches_abbrevs_nonpolygon[sub_cmd[i]]
                                    nvals, keyvals = hp_nonpolygon_abbrevs_func(self, abbrev_key, sub_cmd, i)

                                ### Square specific
                                elif key in holes_patches_abbrevs_square and shape_type == "square":
                                    abbrev_key = holes_patches_abbrevs_square[sub_cmd[i]]
                                    nvals, keyvals = hp_square_abbrevs_func(self, abbrev_key, sub_cmd, i)

                                ### Rectangle specific
                                elif key in holes_patches_abbrevs_rectangle and shape_type == "rectangle":
                                    abbrev_key = holes_patches_abbrevs_rectangle[sub_cmd[i]]
                                    nvals, keyvals = hp_rectangle_abbrevs_func(self, abbrev_key, sub_cmd, i)
                                
                                ### Wrong subcommand given
                                else:
                                    assert False, "subcommand '" + key + "' not usable with square/rectangle"
                                for subkey, subval in keyvals:
                                    cmd_vals[subkey] = subval
                                i += nvals
                            
                            settings_dict[dict_target][hole_or_patch][-1].update({
                                "points": [(cmd_vals["cx"], cmd_vals["cy"])],
                            })
                            del cmd_vals["cx"]
                            del cmd_vals["cy"]
                        
                        elif shape_type in ["polygon"]:
                            assert len([val for val in sub_cmd if val in ["point", "p"]]) >= 3, "\n".join([
                                "If you use the (polygon) shape in a (" + hole_or_patch + ") subcommand,",
                                "then you must supply a minimum of three points."
                            ])
                            cmd_vals["points"] = []
                            i = 0
                            while i < len(sub_cmd):
                                key = sub_cmd[i]
                                ### General
                                if key in holes_patches_abbrevs_general:
                                    abbrev_key = holes_patches_abbrevs_general[sub_cmd[i]]
                                    nvals, keyvals = hp_general_abbrevs_func(self, abbrev_key, sub_cmd, i)

                                ### Center designation
                                elif key in holes_patches_abbrevs_nonpolygon:
                                    abbrev_key = holes_patches_abbrevs_nonpolygon[sub_cmd[i]]
                                    nvals, keyvals = hp_nonpolygon_abbrevs_func(self, abbrev_key, sub_cmd, i)

                                ### Polygon specific
                                elif key in holes_patches_abbrevs_polygon:
                                    abbrev_key = holes_patches_abbrevs_polygon[sub_cmd[i]]
                                    nvals, keyvals = hp_polygon_abbrevs_func(self, abbrev_key, sub_cmd, i)
                                
                                ### Wrong subcommand given
                                else:
                                    assert False, "subcommand '" + key + "' not usable with poly/polygon"
                                
                                for subkey, subval in keyvals:
                                    if subkey == "point":
                                        cmd_vals["points"].append(subval)
                                    else:
                                        cmd_vals[subkey] = subval
                                
                                i += nvals

                        for key, val in cmd_vals.items():
                            settings_dict[dict_target][hole_or_patch][-1][key] = val
                        
                    ### Center [nm]
                    elif any(sub_cmd[0].lower() == cen for cen in ["center"]):
                        if len(sub_cmd[1:]) == 3:
                            settings_dict["membrane"]["center"] = [ast.literal_eval(i) for i in sub_cmd[1:]]
                        elif len(sub_cmd[1:]) == 2:
                            settings_dict["membrane"]["center"] = [ast.literal_eval(i) for i in sub_cmd[1:]] + [0]
                    
                    elif sub_cmd[0].lower() == "cx":
                        settings_dict["membrane"]["center"][0] = ast.literal_eval(sub_cmd[1])
                    
                    elif sub_cmd[0].lower() == "cy":
                        settings_dict["membrane"]["center"][1] = ast.literal_eval(sub_cmd[1])
                    
                    elif sub_cmd[0].lower() == "cz":
                        settings_dict["membrane"]["center"][2] = ast.literal_eval(sub_cmd[1])

                    ### Area per lipid [nm^2] converted to [Å^2]
                    elif sub_cmd[0].lower() == "apl":
                        settings_dict[dict_target]["apl"] = ast.literal_eval(sub_cmd[1])

                    ### Wiggle room. Minimum distance from a bead to the edge of a bin [Å]
                    elif sub_cmd[0].lower() in ["plane_wr", "plane_buffer"]:
                        settings_dict[dict_target]["plane_buffer"] = ast.literal_eval(sub_cmd[1])

                    elif sub_cmd[0].lower() in ["height_wr", "height_buffer"]:
                        settings_dict[dict_target]["height_buffer"] = ast.literal_eval(sub_cmd[1])

                    ### Size of protein points when creating initial 2D-protein polygons
                    elif sub_cmd[0].lower() == "prot_buffer":
                        prot_dict[dict_target] = ast.literal_eval(sub_cmd[1])

                    ### Integer multiplier for radius used in alphashape function [multiplier]
                    elif sub_cmd[0].lower() == "alpha_mult":
                        prot_dict[dict_target] = ast.literal_eval(sub_cmd[1])

                    ### Rounding method for rounding number of lipids in leaflet [function]
                    elif sub_cmd[0].lower() == "round": # "int", "round", "min" or "max"
                        rounding_types = {"int": ("int", int), "round": ("round", round), "floor": ("math.floor", math.floor), "ceil": ("math.ceil", math.ceil)}
                        settings_dict[dict_target]["lip_round_func"] = rounding_types[sub_cmd[1]] 

                    elif sub_cmd[0].lower() in ["kick", "kickxy", "kickz"]:
                        ### Random kick to beads x/y/z positions [nm] #[Å]
                        if sub_cmd[0].lower() == "kick":
                            settings_dict[dict_target]["kickxy"] = ast.literal_eval(sub_cmd[1])
                            settings_dict[dict_target]["kickz"] = ast.literal_eval(sub_cmd[1])
                        else:
                            settings_dict[dict_target][sub_cmd[0].lower()] = ast.literal_eval(sub_cmd[1])

                    ### Bead distance scaling for z is 0.3 by default and 0.25 by default for x and y [multiplier]
                    elif sub_cmd[0].lower() in ["bdx", "bdy", "bdz"]:
                        settings_dict[dict_target][sub_cmd[0].lower()] = ast.literal_eval(sub_cmd[1])

                    ### ### Rectangle specific subcommands:
                    ### xy dimensions of leaflet [nm]. Converted to [Å]
                    elif sub_cmd[0].lower() in ["xlength", "ylength"]:
                        settings_dict["membrane"][sub_cmd[0].lower()] = ast.literal_eval(sub_cmd[1])

                    ### True/False whether to check for pbc crossing
                    elif sub_cmd[0].lower() == "pbc_check":
                        settings_dict["membrane"]["pbc_check"] = ast.literal_eval(sub_cmd[1])

                    ### Lipid optimization method
                    elif sub_cmd[0].lower() == "lipid_optim":
                        valid_lipid_optim = sub_cmd[1] in ["avg_optimal", "abs_val", "force_fill", "fill", "no", "insane"]
                        assert valid_lipid_optim, "Invalid lipid selection method: '" + str(sub_cmd[1]) + "'"
                        settings_dict[dict_target]["lipid_optim"] = sub_cmd[1]

                    ### Designates which force field to collect lipids from
                    elif sub_cmd[0].lower() == "params":
                        settings_dict[dict_target]["params"] = sub_cmd[1]

                    ### Designates how lipid charges should be collected
                    elif sub_cmd[0].lower() == "charge": # "lib" or "top"
                        assert sub_cmd[1] in ["topology", "top", "library", "lib"], "The 'charge' subcommand must be given either 'topology'/'top' or 'library'/'lib' as value. You have given: " + str(sub_cmd[1])
                        if sub_cmd[1] in ["topology", "top"]:
                            memb_dict["charge"] = "top"
                        elif sub_cmd[1] in ["library", "lib"]:
                            memb_dict["charge"] = "lib"
                    
                    elif sub_cmd[0].lower() in ["grid_maker"]:
                        assert sub_cmd[1] in ["lines", "3D_matrix"], "'grid_maker' setting only accepts 'lines' and '3D_matrix'"
                        settings_dict[dict_target]["grid_maker"] = sub_cmd[1]

                    elif sub_cmd[0].lower() in ["ld", "lipid_dist", "lipid_distribution"]:
                        val = sub_cmd[1]
                        if val in ["e", "even", "evenly"]:
                            settings_dict[dict_target]["lipid_distribution"] = "evenly"
                        elif val in ["r", "ran", "rand", "random"]:
                            settings_dict[dict_target]["lipid_distribution"] = "random"
                    
                    #############################
                    ### OPTIMIZATION SETTINGS ###
                    #############################
                    elif sub_cmd[0].lower() in ["optim", "optimize", "optimization", "minim", "minimize", "minimization"]:
                        if sub_cmd[1] in ["False", "0", "no"]:
                            settings_dict[dict_target]["optimize"] = False
                        else:
                            settings_dict[dict_target]["optimize"] = sub_cmd[1].lower()
                    
                    elif sub_cmd[0].lower() in ["optim_maxsteps", "minim_maxsteps"]:
                        isnumber, isint = self.is_number(sub_cmd[1])
                        if isnumber:
                            settings_dict[dict_target]["optim_maxsteps"] = int(ast.literal_eval(sub_cmd[1]))
                    
                    elif sub_cmd[0].lower() in ["optim_push_tol", "minim_push_tol"]:
                        isnumber, isint = self.is_number(sub_cmd[1])
                        if isnumber:
                            settings_dict[dict_target]["optim_push_tol"] = ast.literal_eval(sub_cmd[1])
                    
                    elif sub_cmd[0].lower() in ["optim_push_mult", "minim_push_mult"]:
                        isnumber, isint = self.is_number(sub_cmd[1])
                        if isnumber:
                            settings_dict[dict_target]["optim_push_mult"] = ast.literal_eval(sub_cmd[1])
                    
                    ### Readjusts the max/min x/y-values if cutouts/holes are made
                    elif sub_cmd[0].lower() == "readjust_bbox":
                        settings_dict["membrane"]["readjust_bbox"] = ast.literal_eval(sub_cmd[1])
                    
                    ### Splits subleaflet bbox into multiple subleaflets if they are not all connected
                    elif sub_cmd[0].lower() == "split_bbox":
                        settings_dict["membrane"]["split_bbox"] = ast.literal_eval(sub_cmd[1])
                    
                    ### Splits subleaflet bbox into multiple subleaflets if they are not all connected
                    elif sub_cmd[0].lower() == "solvate_hole":
                        settings_dict["membrane"]["solvate_hole"] = ast.literal_eval(sub_cmd[1])
                    
                    ### Used for nanodiscs. Sets whether the resulting membrane should be inside or outside of a protein
                    elif sub_cmd[0].lower() == "inside_protein":
                        settings_dict["membrane"]["inside_protein"] = ast.literal_eval(sub_cmd[1])
                    
                    elif sub_cmd[0].lower() == "lipid":
                        if "lipids_preprocessing" not in settings_dict[dict_target].keys():
                            settings_dict[dict_target]["lipids_preprocessing"] = []
                        settings_dict[dict_target]["lipids_preprocessing"].append(sub_cmd[1:])
                    
                    else:
                        assert False, "Unknown subcommand given to '-membrane'. The subcommand is: '" + str(sub_cmd) + "'"
                
                if layer_definition == "bilayer" and len(settings_dict["default"]["lipids_preprocessing"]) == 0:
                    if len(settings_dict["upper_leaf"]) > 0 and len(settings_dict["lower_leaf"]) == 0:
                        layer_definition = "upper"
                    if len(settings_dict["lower_leaf"]) > 0 and len(settings_dict["upper_leaf"]) == 0:
                        layer_definition = "lower"
                
                memb_dict = {"leaflets": {}, "membrane_type": layer_definition}
                
                if layer_definition == "upper" or layer_definition == "bilayer":
                    memb_dict["leaflets"]["upper_leaf"] = settings_dict["upper_leaf"]
                    memb_dict["leaflets"]["upper_leaf"].update({
                        "HG_direction": "up",
                        "HG_sign":      +1,
                        "leaflet_type": "upper",
                    })
                
                if layer_definition == "lower" or layer_definition == "bilayer":
                    memb_dict["leaflets"]["lower_leaf"] = settings_dict["lower_leaf"]
                    memb_dict["leaflets"]["lower_leaf"].update({
                        "HG_direction": "down",
                        "HG_sign":      -1,
                        "leaflet_type": "lower",
                    })
                
                ### Adding membrane-wide settings to specific leaflets if they are not already given for the specific leaflet
                for key, vals in settings_dict["membrane"].items():
                    for leaflet in memb_dict["leaflets"].values():
                        if key not in leaflet:
                            try:
                                leaflet[key] = copy.deepcopy(vals)
                            except:
                                leaflet[key] = vals
                
                ### Adding default settings for leaflets if none were given
                for key, vals in settings_dict["default"].items():
                    for leaflet in memb_dict["leaflets"].values():
                        if key not in leaflet:
                            try:
                                leaflet[key] = copy.deepcopy(vals)
                            except:
                                leaflet[key] = vals
                
                ### fixing a couple of values such that they are in angstrom
                for leaflet_type, leaflet in memb_dict["leaflets"].items():
                    leaflet["center"] = [val*10 for val in leaflet["center"]]
                    leaflet["xlength"] *= 10
                    leaflet["ylength"] *= 10
                    ### Ensuring pbc boundaries are not exceeded
                    cx, cy, cz = leaflet["center"]
                    ### Checking if x and y centers are inside the box. Exits if not.
                    if cx > self.pbcx/2+0.0001: # small addition to pbc to prevent activation by float errors
                        self.print_term("Leaflet center exceeds x-bounds of box in the positive x-direction. Exiting.", warn=True)
                        self.print_term("Leaflet center:", leaflet["center"], warn=True)
                        sys.exit()
                    if cx < -(self.pbcx/2+0.0001): # small addition to pbc to prevent activation by float errors
                        self.print_term("Leaflet center exceeds x-bounds of box in the negative x-direction. Exiting.", warn=True)
                        self.print_term("Leaflet center:", leaflet["center"], warn=True)
                        sys.exit()
                    if cy > self.pbcy/2+0.0001: # small addition to pbc to prevent activation by float errors
                        self.print_term("Leaflet center exceeds y-bounds of box in the positive y-direction. Exiting.", warn=True)
                        self.print_term("Leaflet center:", leaflet["center"], warn=True)
                        sys.exit()
                    if cy < -(self.pbcy/2+0.0001): # small addition to pbc to prevent activation by float errors
                        self.print_term("Leaflet center exceeds y-bounds of box in the negative y-direction. Exiting.", warn=True)
                        self.print_term("Leaflet center:", leaflet["center"], warn=True)
                        sys.exit()
                    
                    ### Checking if x and y lengths stretch outside the box. Warns and continues.
                    ### The problem is corrected in "subleaflet_poly_maker" where a polygon stretching the pbc is used to delimit membranes
                    ### small addition to pbc included to prevent activation by float rounding errors
                    xupper, xlower = cx + leaflet["xlength"]/2, cx - leaflet["xlength"]/2
                    if (xupper > self.pbcx/2+0.0001) or (xlower < -(self.pbcx/2+0.0001)):
                        if xupper > self.pbcx/2+0.0001: 
                            self.print_term("Leaflet exceeds x-bounds of box in the positive x-direction. Will be corrected to fit pbc.", warn=True)
                        if xlower < -(self.pbcx/2+0.0001): # small addition to pbc to prevent activation by float errors
                            self.print_term("Leaflet exceeds x-bounds of box in the negative x-direction. Will be corrected to fit pbc.", warn=True)

                    yupper, ylower = cy + leaflet["ylength"]/2, cy - leaflet["ylength"]/2
                    if (yupper > self.pbcy/2+0.0001) or (ylower < -(self.pbcy/2+0.0001)):
                        if yupper > self.pbcy/2+0.0001: # small addition to pbc to prevent activation by float errors
                            self.print_term("Leaflet exceeds y-bounds of box in the positive y-direction. Will be corrected to fit pbc.", warn=True)
                        if ylower < -(self.pbcy/2+0.0001): # small addition to pbc to prevent activation by float errors
                            self.print_term("Leaflet exceeds y-bounds of box in the negative y-direction. Will be corrected to fit pbc.", warn=True)

                    leaflet["apl"] *= 100
                    leaflet["prot_buffer"] *= 10
                    leaflet["plane_buffer"] *= 10
                    leaflet["height_buffer"] *= 10
                    leaflet["kickxy"] *= 10
                    leaflet["kickz"] *= 10
                    
                    for holes_patches in [leaflet["hole"], leaflet["patch"]]:
                        for polygon in holes_patches:
                            polygon["buffer"] *= 10
                            ### Added if statement because scaling factors are treated as radii and sidelengths for non-polygon
                            ### But should just be treated as actual scaling value specifically for polygons
                            if polygon["shape_type"] != "polygon":
                                polygon["xscaling"] *= 10
                                polygon["yscaling"] *= 10
                            for pi, (xval, yval) in enumerate(polygon["points"]):
                                polygon["points"][pi] = (xval*10 + leaflet["center"][0], yval*10 + leaflet["center"][1])

                ################################
                ### Lipid data incorporation ###
                ################################
                
                ### Reconfigures lipid-specific data for leaflet according to subcommands
                ### Adding zero to "memb_directional_heights" to ensure "min" or "max" defaults to zero
                memb_directional_heights = [0]
                memb_mean_directional_heights = [0]

                for leaflet_key, leaflet in memb_dict["leaflets"].items():
                    tot_ratio = 0
                    assert len(leaflet["lipids_preprocessing"]) > 0, "No lipids given to '-membrane' flag. Please specify at least one lipid if you use it."
                    for sub_cmd in leaflet["lipids_preprocessing"]:
                        """
                        Checks if parameters specified for lipid.
                        Order: lipid specific parameters, leaflet specific parameters, general system parameters (defaults to "default")
                        """

                        name = False
                        ratio = 1
                        params = leaflet["params"] or self.lipid_params or self.sys_params
                        charge = leaflet["charge"]
                        
                        all_subcmd_names = ["params", "charge", "ratio", "name"]

                        i = 0
                        while i < len(sub_cmd):
                            if sub_cmd[i] == "params":
                                params = sub_cmd[i+1]
                                i += 2
                            elif sub_cmd[i] == "charge":
                                assert sub_cmd[i+1] in ["top", "topology", "lib", "library"], "Only 'top' and 'lib' are allowed values for lipid specific 'charge'" + "\n" + "Command: " + str(sub_cmd)
                                if sub_cmd[i+1] in ["top", "topology"]:
                                    charge = "top"
                                elif sub_cmd[i+1] in ["lib", "library"]:
                                    charge = "lib"
                                i += 2
                            elif sub_cmd[i] == "ratio":
                                ratio = ast.literal_eval(sub_cmd[i+1])
                                i += 2
                            elif sub_cmd[i] == "name":
                                name = sub_cmd[i+1]
                                i += 2
                            else: ### Assumes that it is the lipid name and potentially ratio
                                name = sub_cmd[i]
                                if len(sub_cmd[i:]) > 1 and sub_cmd[i+1] not in all_subcmd_names and self.get_number_from_string(sub_cmd[i+1]) is not False:
                                    ratio = ast.literal_eval(sub_cmd[i+1])
                                    i += 1
                                i += 1
                        
                        ### Checks if lipid name has been given and if it exists in the specified parameter library
                        assert name is not False, "A name has not been given for a lipid. Full lipid command: " + str(sub_cmd)
                        assert name in self.lipid_dict[params].keys(), "Lipid name '{name}' was not found in the parameter library '{params}'".format(name=name, params=params)
                        
                        leaflet["lipids"][name] = copy.deepcopy(self.lipid_dict[params][name])

                        leaflet["lipids"][name].ratio_add(ratio)
                        tot_ratio += leaflet["lipids"][name].ratio
                        
                        if leaflet["lipids"][name].moleculetype is not False:
                            moleculetype = leaflet["lipids"][name].moleculetype
                        else:
                            moleculetype = False

                        ### Finds charge data from topology files
                        if len(self.ITP_INPUT_cmds) > 0 and charge == "top":
                            assert moleculetype is not False, "\n".join([
                                "Charges are set to be obtained from topology, but this molecule does not have a moleculetype: {name}".format(name=name),
                                "If you have imported the molecule using 'molecule_import' without specifying a moleculetype, then please change membrane command to the following:",
                                "'lipid:" + ":".join(sub_cmd) + ":charge:lib'",
                            ])
                            assert moleculetype in self.itp_moleculetypes.keys(), "\n".join([
                                "Lipid name '{moleculetype}' could not be found in the topology.".format(moleculetype=moleculetype),
                                "You can exempt all lipids (in this membrane command) from being searched for in the topology by adding the following to the membrane command:",
                                "    "+"charge:lib",
                                "Alternatively, you can exempt a lipid from being searched for in the topology by adding the following to the lipid command:",
                                "    "+"charge:lib",
                                "Example:",
                                "    "+"lipid:POPC:5:charge:lib",
                            ])
                            bead_charges = list(map(
                                self.get_number_from_string,
                                self.itp_moleculetypes[moleculetype].topology_types_in_molecule["atoms"].get_value_type("charge")
                            ))
                            beads_in_structure = leaflet["lipids"][name].get_bead_names()
                            try:
                                beads_in_structure_len = len(beads_in_structure)
                            except:
                                beads_in_structure_len = "NaN"
                            beads_in_topology = [entry.atom for entry in self.itp_moleculetypes[moleculetype].topology_types_in_molecule["atoms"].entries.values()]
                            try:
                                beads_in_topology_len = len(beads_in_topology)
                            except:
                                beads_in_topology_len = "NaN"
                            assert len(bead_charges) == len(leaflet["lipids"][name].get_beads()), "\n".join([
                                "Mismatch found between number of beads in structure and number of beads in topology for \"{name}\" during preprocessing of membrane command".format(name=name),
                                "Beads in structure:",
                                "    " + "Number of beads: " + str(beads_in_structure_len),
                                "    " + "Beads: " + str(beads_in_structure),
                                "Beads in topology:",
                                "    " + "Number of beads: " + str(beads_in_topology_len),
                                "    " + "Beads: " + str(beads_in_topology),
                            ])
                            leaflet["lipids"][name].set_bead_charges(bead_charges)
                    
                    lipid_extremes_dict = {"max_xys": [], "max_zs": [], "min_xys": [], "min_zs": []}
                    for lipid in leaflet["lipids"].keys():
                        for func, func_str, sign in [(min, "min", +1), (max, "max", -1)]:
                            for ax, wr in [("xy", "plane"), ("z", "height")]:
                                lipid_val = round(func([
                                    val + (leaflet["kick" + ax] + leaflet[wr + "_buffer"]) * sign
                                    for val in leaflet["lipids"][lipid].get_coords(ax)[0]
                                ]), 2)
                                lipid_extremes_dict[func_str+"_"+ax+"s"].append(lipid_val)

                    max_max_xys = max(lipid_extremes_dict["max_xys"])
                    max_max_zs  = max(lipid_extremes_dict["max_zs"])
                    max_min_xys = max(lipid_extremes_dict["min_xys"])
                    max_min_zs  = max(lipid_extremes_dict["min_zs"])
                    min_max_xys = min(lipid_extremes_dict["max_xys"])
                    min_max_zs  = min(lipid_extremes_dict["max_zs"])
                    min_min_xys = min(lipid_extremes_dict["min_xys"])
                    min_min_zs  = min(lipid_extremes_dict["min_zs"])
                    
                    lipid_radii = [
                            lipid_class.get_radius(AXs="xy")
                            for lipid_class in leaflet["lipids"].values()
                    ]
                    
                    lipid_heights = [
                            lipid_class.get_radius(AXs="z") * 2
                            for lipid_class in leaflet["lipids"].values()
                    ]
                    
                    memb_directional_heights.append(max(lipid_heights) * leaflet["HG_sign"])
                    
                    lipid_heights_ratios = [
                            lipid_class.get_radius(AXs="z") * 2 * lipid_class.ratio / tot_ratio
                            for lipid_class in leaflet["lipids"].values()
                    ]
                    
                    memb_mean_directional_heights.append(np.sum(lipid_heights_ratios) * leaflet["HG_sign"])
                    
                    ### Update leaflet dictionary. Some are duplicated here, but potentially overwritten later
                    leaflet["lipid_dimensions"] = {
                        "max_max_xys": max_max_xys,
                        "max_max_zs":  max_max_zs,
                        "max_min_xys": max_min_xys,
                        "max_min_zs":  max_min_zs,
                        "min_max_xys": min_max_xys,
                        "min_max_zs":  min_max_zs,
                        "min_min_xys": min_min_xys,
                        "min_min_zs":  min_min_zs,
                        "lipid_height": max(lipid_heights),
                        "lipid_radius": max(lipid_radii),
                    }
                    
                memb_dict.update({
                    "bead_maxz"      : max(memb_directional_heights),
                    "bead_minz"      : min(memb_directional_heights),
                    "bead_mean_maxz" : max(memb_mean_directional_heights),
                    "bead_mean_minz" : min(memb_mean_directional_heights),
                    "membrane_height": sum([leaflet["lipid_dimensions"]["lipid_height"] for leaflet in memb_dict["leaflets"].values()]),
                })
                
                if return_self == "self":
                    self.MEMBRANES[cmd_nr] = memb_dict.copy()
                elif return_self == "return":
                    MEMBRANES[cmd_nr] = memb_dict.copy()
        
            if return_self == "return":
                return MEMBRANES[cmd_nr]
            
            self.print_term("Number of membranes preprocessed:", len(self.MEMBRANES), spaces=1, verbose=2+extra_verbose)

