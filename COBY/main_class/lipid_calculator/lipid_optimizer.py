import numpy as np
import math

class lipid_optimizer:
    def lipid_optimizer(self, settings):
        ### Rounding max number of possible lipids according to requested rounding method
        if settings["lipid_optim"] == "insane":
            nlipidsx = abs(settings["xmax"] - settings["xmin"]) / math.sqrt(settings["apl"])
            nlipidsy = abs(settings["ymax"] - settings["ymin"]) / math.sqrt(settings["apl"])
            if "intersection" in settings and settings["intersection"]:
                area_ratio = math.sqrt(settings["holed_bbox"].area / settings["box_poly"].area)
                ### Using 'round' here instead of 'int' because insane-protein compatability part is very iffy
                ### Definitely not the same as insane, but too much has to be changed for it to be possible to make systems identical to insane systems
                nlipidsx = round(nlipidsx * area_ratio)
                nlipidsy = round(nlipidsy * area_ratio)
            max_lipids_possible_decimal = int(nlipidsx) * int(nlipidsy)
            max_lipids_possible         = int(nlipidsx) * int(nlipidsy)
        else:
            subleaflet_area = settings["holed_bbox"].area
            max_lipids_possible_decimal = subleaflet_area / settings["apl"]
            if settings["lip_round_func"][1] == int: # int rounding
                max_lipids_possible = int(max_lipids_possible_decimal + 0.5)
            else: # round(), math.floor() or math.ceil() rounding
                max_lipids_possible = settings["lip_round_func"][1](max_lipids_possible_decimal)

        ### Initial rounded estimations
        lipids = [(lipid_name, lipid_vals.ratio) for lipid_name, lipid_vals in settings["lipids"].items()]
        lipid_names = [name for name, ratio in lipids]
        lipid_ratios_decimal = [ratio for name, ratio in lipids]
        lipids_tot = sum(lipid_ratios_decimal)
        lipid_ratios = [round(int(i / lipids_tot * max_lipids_possible), 3) for i in lipid_ratios_decimal]

        original_ratios = lipid_ratios[:]
        original_ratios_decimal = [round(ratio / sum(original_ratios) * 100, 3) for ratio in original_ratios]

        expected_ratios_decimal = [round(ratio / sum(lipid_ratios_decimal) * 100, 3) for ratio in lipid_ratios_decimal]

        ### Just take integer lipid values and remove excess grid points
        if settings["lipid_optim"] == "abs_val":
            lipid_ratios = [ratio for name, ratio in lipids]
            ### Error out if values are not integers
            if not all([type(val) == int for val in lipid_ratios]):
                self.print_term(
                    "Not all supplied lipid values are integers.",
                    "\n" + "rounding all values to integers.",
                    warn=True,
                )
            lipid_ratios = [int(ratio + 0.5) for ratio in lipid_ratios]
            ### Warn if more lipids than available grid points based on current apl
            if sum(lipid_ratios) > max_lipids_possible:
                self.print_term(
                    "You have specified more lipids than would be typically allowed based on the leaflets area and area-per-lipid.",
                    "\n" + "Maximum number of lipids: ", max_lipids_possible,
                    "\n" + "Supplied number of lipids:", sum(lipid_ratios),
                    warn=True,
                )

        elif settings["lipid_optim"] in ["force_fill", "fill", "avg_optimal"]:
            ###############################
            ### OPTIMIZING LIPID RATIOS ###
            ###############################
            iters = [lipid_ratios[:]]

            while sum(lipid_ratios) < max_lipids_possible:
                diffs = [(round(lipids[i][1], 3) - round(lip_val / (sum(lipid_ratios) / lipids_tot), 3), i) for i, lip_val in enumerate(lipid_ratios)]
                biggest_diff, biggest_diff_i = max(diffs, key=lambda d: d[0])

                ### In case ideal ratio is found, only continue if force_fill
                if biggest_diff == float(0) and settings["lipid_optim"] == "force_fill":
                    lipid_ratios[lipid_ratios.index(min(lipid_ratios))] += 1
                    iters.append(lipid_ratios)
                elif biggest_diff == float(0):
                    break
                else:
                    lipid_ratios[biggest_diff_i] += 1
                    iters.append(lipid_ratios[:])

            ### Filling in lipids
            if settings["lipid_optim"] in ["force_fill", "fill"]:
                lipid_ratios = iters[-1]

            ### Optimizing lipids
            elif settings["lipid_optim"] == "avg_optimal": # lipid_ratios_decimal
                iters_decimal = [[round(iteration[i] / sum(iteration) * 100, 3) for i in range(len(iteration))] for iteration in iters]
                iters_abs_diff = [
                    [abs(iteration[i] - expected_ratios_decimal[i]) for i in range(len(iteration))]
                    for iteration in iters_decimal
                ]

                iters_avg = [(np.mean(iters), i) for i, iters in enumerate(iters_abs_diff)]
                iters_avg_optimal = min(iters_avg, key=lambda i: i[0])
                lipid_ratios = iters[iters_avg_optimal[1]]

        elif settings["lipid_optim"] in ["no", "insane"]:
            ### E.g. do nothing. Just here to show that the options have not been overlooked.
            pass
        
        settings["lipids"]                      = lipids
        settings["lipid_names"]                 = lipid_names
        settings["lipid_ratios_decimal"]        = lipid_ratios_decimal
        settings["lipids_tot"]                  = lipids_tot
        settings["lipid_ratios"]                = lipid_ratios
        settings["original_ratios"]             = original_ratios
        settings["original_ratios_decimal"]     = original_ratios_decimal
        settings["expected_ratios_decimal"]     = expected_ratios_decimal
        settings["max_lipids_possible"]         = max_lipids_possible
        settings["max_lipids_possible_decimal"] = max_lipids_possible_decimal
        
        return settings
    
