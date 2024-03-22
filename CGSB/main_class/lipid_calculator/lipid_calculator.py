import time
import math
import numpy as np
import copy

class lipid_calculator:
    def lipid_calculator(self):
        self.SYS_lipids_dict = {}
        if len(self.MEMBRANES) != 0:
            lipid_calculator_tic = time.time()
            string = " ".join(["", "CALCULATING LIPID RATIOS", ""])
            self.print_term("{string:-^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), spaces=0, verbose=1)
            printer_spacing = 1
            for memb_i, (memb_key, memb_dict) in enumerate(self.MEMBRANES.items()):
                membrane_lipid_count_dict = {}
                if memb_i != 0:
                    self.print_term("", verbose=2)
                self.print_term("Starting membrane nr", memb_key, spaces=0, verbose=2)
                
                ### Checks if any leaflets contain multiple subleaflets for printer spacing
                for leaflet_i, (leaflet_key, leaflet) in enumerate(memb_dict["leaflets"].items()):
                    if len(leaflet["subleaflets"]) > 1 and self.verbose >= 3:
                        printer_spacing = 2
                
                for leaflet_i, (leaflet_key, leaflet) in enumerate(memb_dict["leaflets"].items()):
                    if leaflet_i != 0:
                        self.print_term("", verbose=2)
                    self.print_term("Starting", leaflet["leaflet_type"], "leaflet", spaces=1, verbose=2)
                    
                    ### Lipid optimization printing
                    ### Printing here so it isn't spammed if multiple subleafs
                    self.print_term("Lipid optimization 'lipid_optim' setting:", leaflet["lipid_optim"], spaces=2, verbose=2)
                    if leaflet["lipid_optim"] in ["force_fill", "fill"]:
                        if leaflet["lipid_optim"] == "fill":
                            self.print_term("Filling leaflet until perfect ratio between lipids is achieved or leaflet is full", spaces=3, verbose=2)
                        elif leaflet["lipid_optim"] == "force_fill":
                            self.print_term("Forcefully filling leaflet until all grid points have a lipid", spaces=3, verbose=2)
                    elif leaflet["lipid_optim"] == "avg_optimal":
                        self.print_term("Optimizing based on the average deviation from expected ratios", spaces=3, verbose=2)
                    self.print_term("", verbose=3)
                    
                    main_settings_dict = {
                        "leaf_lipid_count_dict": {},
                        "lipid_optim":           leaflet["lipid_optim"],
                        "apl":                   leaflet["apl"],
                        "apl_sqrt":              math.sqrt(leaflet["apl"]),
                        "lipids":                leaflet["lipids"],
                        "lipid_names":           [lipid_name for lipid_name in leaflet["lipids"].keys()],
                        "lipid_ratios":          [0 for lipid_name in leaflet["lipids"].keys()],
                        "lip_round_func":        leaflet["lip_round_func"],
                    }
                    
                    settings_dict = copy.deepcopy(main_settings_dict)
                    settings_dict.update({
                        "xmin":         leaflet["xmin"],
                        "xmax":         leaflet["xmax"],
                        "ymin":         leaflet["ymin"],
                        "ymax":         leaflet["ymax"],
                        "holed_bbox":   leaflet["holed_bbox"],
                        "intersection": leaflet["intersection"],
                        "box_poly":     leaflet["box_poly"],
                    })
                    
                    leaflet_lipid_data = self.lipid_optimizer(settings_dict)
                    
                    leaflet["lipid_names"]           = leaflet_lipid_data["lipid_names"]
                    leaflet["lipid_ratios"]          = leaflet_lipid_data["lipid_ratios"]
                    leaflet["leaf_lipid_count_dict"] = {}
                        
                    for subleaflet_i, ((slxi, slyi, sli), subleaflet) in enumerate(leaflet["subleaflets"].items()):
                        if subleaflet_i != 0:
                            self.print_term("", verbose=3)
                        if len(leaflet["subleaflets"]) > 1:
                            self.print_term("Starting subleaflet", str(subleaflet_i+1)+"/"+str(len(leaflet["subleaflets"])), spaces=2, verbose=3)
                        
                        settings_dict = copy.deepcopy(main_settings_dict)
                        settings_dict.update({
                            "xmin":         subleaflet["xmin"],
                            "xmax":         subleaflet["xmax"],
                            "ymin":         subleaflet["ymin"],
                            "ymax":         subleaflet["ymax"],
                            "holed_bbox":   subleaflet["holed_bbox"],
                            "intersection": subleaflet["intersection"],
                            "box_poly":     subleaflet["box_poly"],
                        })
                        subleaflet_lipid_data = self.lipid_optimizer(settings_dict)

                        self.print_term("Maximum number of lipids allowed in subleaflet:", subleaflet_lipid_data["max_lipids_possible"], spaces=3, verbose=3)
                        
                        subleaflet["lipid_names"]  = subleaflet_lipid_data["lipid_names"]
                        subleaflet["lipid_ratios"] = subleaflet_lipid_data["lipid_ratios"]
                        
                        for name, count in zip(subleaflet["lipid_names"], subleaflet["lipid_ratios"]):
                            if name not in leaflet["leaf_lipid_count_dict"]:
                                leaflet["leaf_lipid_count_dict"][name] = count
                            else:
                                leaflet["leaf_lipid_count_dict"][name] += count
                        
                        if self.extra_info:
                            if len(leaflet["subleaflets"]) > 1 and self.verbose >= 3:
                                ### Printing mean, min and max deviation from wanted ratios
                                ratios_final_decimal  = [round(subleaflet_lipid_data["lipid_ratios"][i] / sum(subleaflet_lipid_data["lipid_ratios"]) * 100, 3) for i in range(len(subleaflet_lipid_data["lipid_ratios"]))]
                                ratios_final_abs_diff = [abs(ratios_final_decimal[i] - subleaflet_lipid_data["expected_ratios_decimal"][i]) for i in range(len(ratios_final_decimal))]
                                ratios_final_avg      = round(np.mean(ratios_final_abs_diff), 5) # Don't multiply by 100 as they are already percentages
                                ratios_final_min      = round(min(ratios_final_abs_diff), 5) # Don't multiply by 100 as they are already percentages
                                ratios_final_max      = round(max(ratios_final_abs_diff), 5) # Don't multiply by 100 as they are already percentages

                                self.print_term("Final deviations from expected ratios:", "Difference in %-values", spaces=1+printer_spacing, verbose=3)
                                self.print_term("Maximum:", ratios_final_max, spaces=2+printer_spacing, verbose=3)
                                self.print_term("Average:", ratios_final_avg, spaces=2+printer_spacing, verbose=3)
                                self.print_term("Minimum:", ratios_final_min, spaces=2+printer_spacing, verbose=3)
                        
                        ##################################################
                        ### PRINTING SUBLEAFLET SPECIFIC LIPID DETAILS ###
                        ##################################################
                        SUBLEAF_final_ratios_decimal = [round(ratio / sum(subleaflet_lipid_data["lipid_ratios"]) * 100, 3) for ratio in subleaflet_lipid_data["lipid_ratios"]]
                        SUBLEAF_headers = ["Name", "Final lipids", "Final %", "Expected %", "Starting %", "Starting lipids"]
                        SUBLEAF_lipids_for_printer = list(zip(
                            subleaflet_lipid_data["lipid_names"],
                            subleaflet_lipid_data["lipid_ratios"],
                            SUBLEAF_final_ratios_decimal,
                            subleaflet_lipid_data["expected_ratios_decimal"],
                            subleaflet_lipid_data["original_ratios_decimal"],
                            subleaflet_lipid_data["original_ratios"],
                        ))
                        ### Sorts based on Final % --> Final lipids --> Name
                        SUBLEAF_lipids_for_printer = sorted(SUBLEAF_lipids_for_printer, key=lambda x: (float(x[2]),float(x[1]),x[0]), reverse=True)
                        SUBLEAF_printer = [tuple(SUBLEAF_headers)] + SUBLEAF_lipids_for_printer
                        SUBLEAF_max_lengths = [len(str(max(data, key = lambda d: len(str(d))))) for data in zip(*SUBLEAF_printer)]

                        if self.extra_info and len(leaflet["subleaflets"]) > 1 and self.verbose >= 3:
                            self.print_term(
                                "Subleaflet specific lipid data",
                                "(Max lipids: "+str(subleaflet_lipid_data["max_lipids_possible"])+", Final lipids: "+str(sum(subleaflet_lipid_data["lipid_ratios"]))+")",
                                spaces=1+printer_spacing,
                                verbose=3,
                            )
                            
                            for string in SUBLEAF_printer:
                                for_printer = []
                                for si, substring in enumerate(string):
                                    for_printer.append('{0: <{L}}'.format(substring, L = SUBLEAF_max_lengths[si]))
                                    if si+1 < len(string):
                                        for_printer.append(":")
                                self.print_term(
                                    *for_printer,
                                    spaces=2+printer_spacing,
                                    verbose=3,
                                )
                    
                    mismatch_found = False
                    for li, (l_ratio, (sb_name, sb_ratio)) in enumerate(zip(leaflet["lipid_ratios"], leaflet["leaf_lipid_count_dict"].items())):
                        lipid_deviation = l_ratio - sb_ratio
                        if lipid_deviation != 0:
                            mismatch_found = True
                    
                    self.print_term("", spaces=3, verbose=3)
                    self.print_term("Maximum number of lipids allowed in leaflet:", leaflet_lipid_data["max_lipids_possible"], spaces=2, verbose=2)
                    if mismatch_found:
                        self.print_term("Number of allowed lipids in leaflet does not match the number of lipids in subleaflets for all lipid types", spaces=3, verbose=3)
                        self.print_term("Adjusting number of lipids in subleaflets", spaces=3, verbose=3)
                        for li, (l_ratio, (sb_name, sb_ratio)) in enumerate(zip(leaflet["lipid_ratios"], leaflet["leaf_lipid_count_dict"].items())):
                            lipid_deviation = l_ratio - sb_ratio
                            if lipid_deviation != 0 and len(leaflet["subleaflets"]) > 1:
                                if lipid_deviation > 0:
                                    val = 1
                                elif lipid_deviation < 0:
                                    val = -1
                                for _ in range(val*lipid_deviation):
                                    subleaflet_occupancy = []
                                    for (slxi, slyi, sli), subleaflet in leaflet["subleaflets"].items():
                                        area    = subleaflet["holed_bbox"].area
                                        apl     = leaflet["apl"]
                                        nlipids = np.sum(subleaflet["lipid_ratios"])
                                        subleaflet_occupancy.append(((slxi, slyi, sli), (apl * (nlipids+val)) / area))
                                    if lipid_deviation > 0:
                                        subleaflet_occupancy_sorted = sorted(subleaflet_occupancy, key=lambda x: x[1], reverse=False)[0]
                                    elif lipid_deviation < 0:
                                        subleaflet_occupancy_sorted = sorted(subleaflet_occupancy, key=lambda x: x[1], reverse=True)[0]
                                    subleaflet_key = subleaflet_occupancy_sorted[0]
                                    leaflet["subleaflets"][subleaflet_key]["lipid_ratios"][li] += val
                                    leaflet["leaf_lipid_count_dict"][sb_name] += val
                        
                    leaflet["leaf_lipid_count"] = []
                    for n, c in leaflet["leaf_lipid_count_dict"].items():
                        leaflet["leaf_lipid_count"].append((n, c))
                        if n not in membrane_lipid_count_dict:
                            membrane_lipid_count_dict[n] = c
                            self.SYS_lipids_dict[n]      = c
                        else:
                            membrane_lipid_count_dict[n] += c
                            self.SYS_lipids_dict[n]      += c
                    
                    ###############################################
                    ### PRINTING LEAFLET SPECIFIC LIPID DETAILS ###
                    ###############################################
                    if self.extra_info:
                        ### Printing mean, min and max deviation from wanted ratios
                        ratios_final_decimal = [round(leaflet_lipid_data["lipid_ratios"][i] / sum(leaflet_lipid_data["lipid_ratios"]) * 100, 3) for i in range(len(leaflet_lipid_data["lipid_ratios"]))]
                        ratios_final_abs_diff = [abs(ratios_final_decimal[i] - leaflet_lipid_data["expected_ratios_decimal"][i]) for i in range(len(ratios_final_decimal))]
                        ratios_final_avg = round(np.mean(ratios_final_abs_diff), 5) # Don't multiply by 100 as they are already percentages
                        ratios_final_min = round(min(ratios_final_abs_diff), 5) # Don't multiply by 100 as they are already percentages
                        ratios_final_max = round(max(ratios_final_abs_diff), 5) # Don't multiply by 100 as they are already percentages

                        self.print_term("Final deviations from expected ratios:", "Difference in %-values", spaces=1+printer_spacing, verbose=2)
                        self.print_term("Maximum:", ratios_final_max, spaces=2+printer_spacing, verbose=2)
                        self.print_term("Average:", ratios_final_avg, spaces=2+printer_spacing, verbose=2)
                        self.print_term("Minimum:", ratios_final_min, spaces=2+printer_spacing, verbose=2)
                        
                        LEAF_lipid_names, LEAF_lipid_vals = zip(*leaflet["leaf_lipid_count"])
                        LEAF_tot_lipids = sum(LEAF_lipid_vals)
                        LEAF_lipid_percentages = [round(val / LEAF_tot_lipids * 100, 3) for val in LEAF_lipid_vals]
                        
                        if len(leaflet["subleaflets"]) > 1:
                            LEAF_headers = ["Name", "Total lipids", "Final %", "Expected %"]
                            LEAF_lipids_for_printer = list(zip(
                                LEAF_lipid_names,
                                LEAF_lipid_vals,
                                LEAF_lipid_percentages,
                                leaflet_lipid_data["expected_ratios_decimal"],
                            ))
                        
                        else:
                            LEAF_headers = ["Name", "Final lipids", "Final %", "Expected %", "Starting %", "Starting lipids"]
                            LEAF_lipids_for_printer = list(zip(
                                LEAF_lipid_names,
                                LEAF_lipid_vals,
                                LEAF_lipid_percentages,
                                leaflet_lipid_data["expected_ratios_decimal"],
                                leaflet_lipid_data["original_ratios_decimal"],
                                leaflet_lipid_data["original_ratios"],
                            ))
                        ### Sorts based on Final % --> Final lipids --> Name
                        LEAF_lipids_for_printer = sorted(LEAF_lipids_for_printer, key=lambda x: (float(x[2]),float(x[1]),x[0]), reverse=True)
                        LEAF_printer = [tuple(LEAF_headers)] + LEAF_lipids_for_printer
                        LEAF_max_lengths = [len(str(max(data, key = lambda d: len(str(d))))) for data in zip(*LEAF_printer)]
                        
                        self.print_term("Leaflet specific lipid data (Combined subleafs)", spaces=2, verbose=2)
                        for string in LEAF_printer:
                            for_printer = []
                            for si, substring in enumerate(string):
                                for_printer.append('{0: <{L}}'.format(substring, L = LEAF_max_lengths[si]))
                                if si+1 < len(string):
                                    for_printer.append(":")
                            self.print_term(
                                *for_printer,
                                spaces=2+printer_spacing,
                                verbose=2,
                            )

                ################################################
                ### PRINTING MEMBRANE SPECIFIC LIPID DETAILS ###
                ################################################
                if self.extra_info and len(memb_dict["leaflets"]) > 1:
                    self.print_term("", verbose=2)
                    membrane_lipid_count_tuples = [(n, c) for n, c in membrane_lipid_count_dict.items()]
                    MEMB_lipid_names, MEMB_lipid_vals = zip(*membrane_lipid_count_tuples)
                    MEMB_tot_lipids = sum(MEMB_lipid_vals)
                    MEMB_lipid_percentages = [round(val / MEMB_tot_lipids * 100, 3) for val in MEMB_lipid_vals]

                    MEMB_headers = ["Name", "Total lipids", "Total %"]
                    MEMB_lipids_for_printer = list(zip(MEMB_lipid_names, MEMB_lipid_vals, MEMB_lipid_percentages))
                    ### Sorts based on Final % --> Final lipids --> Name
                    MEMB_lipids_for_printer = sorted(MEMB_lipids_for_printer, key=lambda x: (float(x[2]),float(x[1]),x[0]), reverse=True)
                    MEMB_printer = [tuple(MEMB_headers)] + MEMB_lipids_for_printer
                    MEMB_max_lengths = [len(str(max(data, key = lambda d: len(str(d))))) for data in zip(*MEMB_printer)]

                    self.print_term("Membrane-wide lipid data", spaces=1, verbose=2)
                    for string in MEMB_printer:
                        for_printer = []
                        for si, substring in enumerate(string):
                            for_printer.append('{0: <{L}}'.format(substring, L = MEMB_max_lengths[si]))
                            if si+1 < len(string):
                                for_printer.append(":")
                        self.print_term(
                            *for_printer,
                            spaces=2+printer_spacing,
                            verbose=2,
                        )
            
            #########################################
            ### PRINTING SYSTEMWIDE LIPID DETAILS ###
            #########################################
            if self.extra_info and len(self.MEMBRANES) > 1:
                self.print_term("", verbose=2)
                SYS_lipid_names = [name for name in self.SYS_lipids_dict.keys()]
                SYS_lipid_vals = [val for val in self.SYS_lipids_dict.values()]
                SYS_tot_lipids = sum(SYS_lipid_vals)
                SYS_lipid_percentages = [round(val / SYS_tot_lipids * 100, 3) for val in SYS_lipid_vals]

                SYS_headers = ["Name", "Total lipids", "Total %"]
                SYS_lipids_for_printer = list(zip(SYS_lipid_names, SYS_lipid_vals, SYS_lipid_percentages))
                SYS_printer = [tuple(SYS_headers)] + SYS_lipids_for_printer
                SYS_max_lengths = [len(str(max(data, key = lambda d: len(str(d))))) for data in zip(*SYS_printer)]
                
                self.print_term("\nLipid data for whole system", spaces=0, verbose=2)
                for string in SYS_printer:
                    for_printer = []
                    for si, substring in enumerate(string):
                        for_printer.append('{0: <{L}}'.format(substring, L = SYS_max_lengths[si]))
                        if si+1 < len(string):
                            for_printer.append(":")
                    self.print_term(
                        *for_printer,
                        spaces=2+printer_spacing,
                        verbose=2,
                    )
            
            lipid_calculator_toc = time.time()
            lipid_calculator_time = round(lipid_calculator_toc - lipid_calculator_tic, 4)
            string = " ".join(["", "LIPID RATIO CALCULATIONS COMPLETE", ""])
            self.print_term("{string:-^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), spaces=0, verbose=1)
            string = " ".join(["", "(Time spent:", str(lipid_calculator_time), "[s])", ""])
            self.print_term("{string:^{string_length}}".format(string=string, string_length=self.terminalupdate_string_length), "\n", spaces=0, verbose=1)
