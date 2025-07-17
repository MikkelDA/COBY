import time
import sys
import os
import shapely
import shapely.affinity
import shapely.plotting
import matplotlib.pyplot as plt
from matplotlib import patches
import numpy as np
from operator import itemgetter

class grid_plotter:
    def grid_plotter(self):
        ### Imports inside the method to prevent most people from having to install pyrecorder
        try:
            import pyrecorder
            from pyrecorder.recorder import Recorder
            from pyrecorder.writers.video import Video
            from pyrecorder.writers.gif import GIF
            from pyrecorder.converters.matplotlib import Matplotlib
        except:
            self.print_term("You are trying to use the 'grid plotting' functionality without having the package 'pyrecorder' installed. Please install the package.", warn=True)
            sys.exit()

        l, b, r, t = 80, 80, 40, 40
        dpi = 140
        # l, b, r, t = 110, 60, 30, 20
        # dpi = 280
        pbc_mod = 35
        nrows = 1
        ncols = 1
        row_nr = 0
        col_nr = 0
        size = 5

        do_zfill                      = False
        plot_AllLineStringsIterations = True
        plot_AllLipidInsertions       = True
        article_figures               = False

        ### First make dict with all possible lipids present across all membranes/leaflets
        color_lipid_circle_dict = {}
        for membrane in self.MEMBRANES.values():
            for leaflet in membrane["leaflets"].values():
                for lipid in leaflet["lipids"]:
                    color_lipid_circle_dict[lipid] = False
        
        ### Then update dict with colors
        if self.PLOT_cmd["lipid_colors"] == "type":
            nlipidtypes = len(color_lipid_circle_dict)
            lipid_colors = iter(plt.cm.prism(np.linspace(0, 1, nlipidtypes)))
            for lipid in color_lipid_circle_dict.keys():
                if lipid in self.PLOT_cmd["given_lipid_colors"].keys():
                    color_lipid_circle_dict[lipid] = self.PLOT_cmd["given_lipid_colors"][lipid]
                else:
                    color_lipid_circle_dict[lipid] = next(lipid_colors)
        elif self.PLOT_cmd["lipid_colors"] == "same":
            for lipid in color_lipid_circle_dict.keys():
                color_lipid_circle_dict[lipid] = "blue"

        if article_figures:
            l, b, r, t = 20, 20, 20, 20
            pbc_mod = 30

            color_bbox = "#ff5200ff" # Yellow-red-orange compromise
            alpha_bbox = 0.5

            color_protein_points            = "green"
            color_protein_points_buffered   = "green"
            color_protein_points_exterior   = "green"
            color_protein_points_alphashape = "green"

            color_linestrings  = "black"
            # color_lipid_circle = "black"
            color_lipid_center = "black"
        else:
            l, b, r, t = 60, 60, 40, 40
            pbc_mod = 35

            color_bbox = "red"
            alpha_bbox = 0.3

            color_protein_points            = "orange"
            color_protein_points_buffered   = "blue"
            color_protein_points_exterior   = "blue"
            color_protein_points_alphashape = "green"

            color_linestrings  = "green"
            # color_lipid_circle = "black"
            color_lipid_center = "black"
            # color_lipid_center = "blue"

        def adjuster(figure = None, ax = None, l = None, b = None, r = None, t = None):
            ### Ensures figure size is defined
            try:
                figsize = figure.get_size_inches()*figure.dpi
            except:
                figsize = figs.get_size_inches()*figs.dpi
            assert figsize[0] and figsize[1]

            ### Esnures axis position is defined
            try:
                ax_pos = ax._position
            except:
                try:
                    ax_pos = axs[row_nr, col_nr]._position
                except:
                    ax_pos = [0.07,0.07,0.97,0.97]

            out = {}
            if l != None:
                out["l"] = 100/figsize[0]*l/100 # /100 for fraction
            else:
                out["l"] = ax_pos[0]

            if b != None:
                out["b"] = 100/figsize[1]*b/100 # /100 for fraction
            else:
                out["b"] = ax_pos[1]

            if r != None:
                out["r"] = 1-100/figsize[0]*r/100 # /100 for fraction
            else:
                out["r"] = ax_pos[2]

            if t != None:
                out["t"] = 1-100/figsize[1]*t/100 # /100 for fraction
            else:
                out["t"] = ax_pos[3]

            return out

        def make_figure(given_nrows = False, given_ncols = False):
            if given_nrows is not False:
                nrows = given_nrows
            else:
                nrows = 1

            if given_ncols is not False:
                ncols = given_ncols
            else:
                ncols = 1

            figs, axs = plt.subplots(
                nrows = nrows,
                ncols = ncols,
                figsize = (ncols * (pbcx/pbc_mod + (l+r)/dpi), nrows * (pbcy/pbc_mod + (b+t)/dpi)),
                dpi = dpi,
                squeeze = False
            )
            return figs, axs

        def plot_bbox(i=False):
            if i is not False and vals["LineStringsInfo"][i]["bbox_polygon"]:
                shapely.plotting.plot_polygon(
                    vals["LineStringsInfo"][i]["bbox_polygon"],
                    ax=axs[row_nr, col_nr],
                    add_points=False,
                    color=color_bbox,
                    alpha=alpha_bbox,
                )

            else:
                for ax_row in axs:
                    for ax in ax_row:
                        if vals["holed_bbox"]:
                            shapely.plotting.plot_polygon(
                                vals["holed_bbox"],
                                ax=ax,
                                add_points=False,
                                color=color_bbox,
                                alpha=alpha_bbox,
                            )
        def plot_bbox_original():
            if vals["original_holed_bbox"]:
                shapely.plotting.plot_polygon(
                    vals["original_holed_bbox"],
                    ax=axs[row_nr, col_nr],
                    add_points=False,
                    color=color_bbox,
                    alpha=alpha_bbox,
                )
        def plot_bbox_buffered_for_LineStrings(i=False):
            if i is not False and vals["LineStringsInfo"][i]["bbox_polygon_Buffered"]:
                shapely.plotting.plot_polygon(
                    vals["LineStringsInfo"][i]["bbox_polygon_Buffered"],
                    ax=axs[row_nr, col_nr],
                    add_points=False,
                    color=color_bbox,
                    alpha=alpha_bbox,
                )
    #         elif vals["bbox_polygon_BufferedForLineStrings"]:
    #             shapely.plotting.plot_polygon(
    #                 vals["bbox_polygon_BufferedForLineStrings"],
    #                 ax=axs[row_nr, col_nr],
    #                 add_points=False,
    #                 color=color_bbox,
    #                 alpha=alpha_bbox,
    #             )

        def plot_LineStrings(i):
            if vals["LineStringsInfo"][i]["LineStrings"]:
                for LineString in vals["LineStringsInfo"][i]["LineStrings"]:
                    axs[row_nr, col_nr].plot(
                        *LineString.xy,
                        color=color_linestrings,
                    )
        def plot_GridPoints(i):
            if vals["LineStringsInfo"][i]["GridPoints"]:
                axs[row_nr, col_nr].scatter(
                    *list(zip(*vals["LineStringsInfo"][i]["GridPoints"])),
                    color=color_linestrings,
                )

        def plot_LineStringsXLines(i):
            if vals["LineStringsInfo"][i]["grid_method"] == "LineStrings":
                if vals["LineStringsInfo"][i]["LineStrings"]:
                    info = vals["LineStringsInfo"][i]
                    xlinespace, xspace = np.linspace(start=info["xmin_edge"], stop=info["xmax_edge"], num=info["xlines"], endpoint=True, retstep=True)
                    ylinespace, yspace = np.linspace(start=info["ymin_edge"], stop=info["ymax_edge"], num=info["ylines"], endpoint=True, retstep=True)
                    LineStrings = []

                    for yval in ylinespace:
                        left_point = (xlinespace[-1], yval)
                        right_point = (xlinespace[0], yval)
                        LineString = shapely.LineString([left_point, right_point])
                        LineStrings.append(LineString)
                    for LineString in LineStrings:
                        axs[row_nr, col_nr].plot(
                            *LineString.xy,
                            color=color_linestrings,
                        )
        def plot_LineStringsOverlapping(i):
            if vals["LineStringsInfo"][i]["grid_method"] == "LineStrings":
                if vals["LineStringsInfo"][i]["LineStringsOverlapping"]:
                    for LineStringOverlapping in vals["LineStringsInfo"][i]["LineStringsOverlapping"]:
                        axs[row_nr, col_nr].plot(
                            *LineStringOverlapping.xy,
                            color=color_linestrings,
                        )
        def plot_LineStringsContained(i):
            if vals["LineStringsInfo"][i]["grid_method"] == "LineStrings":
                if vals["LineStringsInfo"][i]["LineStringsContained"]:
                    for LineStringOverlapping in vals["LineStringsInfo"][i]["LineStringsContained"]:
                        axs[row_nr, col_nr].plot(
                            *LineStringOverlapping.xy,
                            color=color_linestrings,
                        )
        def plot_GridPoints_Contained(i):
            if vals["LineStringsInfo"][i]["GridPoints_Contained"]:
                axs[row_nr, col_nr].scatter(
                    *list(zip(*vals["LineStringsInfo"][i]["GridPoints_Contained"])),
                    color=color_linestrings,
                )

        def plot_LipidInsertionStep(i):
            step = vals["lipid_insertion"]["lipid_insertion_steps"][i]
            xcoords = vals["lipid_insertion"]["xcoords"]
            ycoords = vals["lipid_insertion"]["ycoords"]

            ### Grid points
            for lipidtype_i, lipidname in enumerate(list(dict.fromkeys((list(zip(*vals["lipids"]))[0])))):
                grid_bool_matrix    = step["grid_bool_matrix"][lipidtype_i]
                free_indexes        = np.nonzero(grid_bool_matrix)
                free_xs, free_ys    = free_indexes
                free_xs_coords      = xcoords[free_xs]
                free_ys_coords      = ycoords[free_ys]
                nonfree_indexes     = np.nonzero(grid_bool_matrix == 0)
                nonfree_xs, nonfree_ys = nonfree_indexes
                nonfree_xs_coords   = xcoords[nonfree_xs]
                nonfree_ys_coords   = ycoords[nonfree_ys]

                axs[row_nr, lipidtype_i].set_title(lipidname)

                ### Free points
                axs[row_nr, lipidtype_i].scatter(
                    free_xs_coords,
                    free_ys_coords,
                    color="black",
                    s=vals["lipid_insertion"]["gridres"]*8,
                )

                ### Non-free points
                axs[row_nr, lipidtype_i].scatter(
                    nonfree_xs_coords,
                    nonfree_ys_coords,
                    color="green",
                    s=vals["lipid_insertion"]["gridres"]*8,
                )

                for i, xy in enumerate(step["grid_points"]):
                    cur_lipidname = vals["lipids"][i][0]
                    if cur_lipidname in ["POPC", "POPE"]:
                        center_color = "blue"
                    elif cur_lipidname in ["CHOL"]:
                        center_color = "cyan"
                    else:
                        center_color = "red"

                    ### Lipid centers
                    xpoints, ypoints = zip(*step["grid_points"])
                    axs[row_nr, lipidtype_i].scatter(
                        xy[0],
                        xy[1],
                        color=center_color,
                        s=size,
                    )

                    ### Lipid name and radius
                    name, radius = vals["lipids"][i]
                    circle = plt.Circle(
                        xy,
                        radius = radius,
                        fill   = False,
                        color  = color_lipid_circle_dict[name],
                    )
                    axs[row_nr, lipidtype_i].add_artist(circle)

        def plot_PreRandomLipidPlacements():
            xpoints, ypoints = list(zip(*vals["grid_points_no_random"]))
            axs[row_nr, col_nr].scatter(
                xpoints,
                ypoints,
                color = color_lipid_center,
                s     = size,
            )
            for i, xy in enumerate(vals["grid_points_no_random"]):
                name, radius = vals["lipids"][i]
                circle = plt.Circle(
                    xy,
                    radius = radius,
                    fill   = False,
                    color  = color_lipid_circle_dict[name],
                )
                axs[row_nr, col_nr].add_artist(circle)
        def plot_PostRandomLipidPlacements():
            frame0_post_random = vals["POINT_STEPS"][0]
            xpoints, ypoints = frame0_post_random[:,0], frame0_post_random[:,1]
            axs[row_nr, col_nr].scatter(
                xpoints,
                ypoints,
                color = color_lipid_center,
                s     = size,
            )
            for i, xy in enumerate(frame0_post_random):
                name, radius = vals["lipids"][i]
                circle = plt.Circle(
                    xy,
                    radius = radius,
                    fill   = False,
                    color  = color_lipid_circle_dict[name],
                )
                axs[row_nr, col_nr].add_artist(circle)
        def plot_FinalLipidPlacements():
            final_frame = vals["POINT_STEPS"][-1]
            xpoints, ypoints = final_frame[:,0], final_frame[:,1]
            axs[row_nr, col_nr].scatter(
                xpoints,
                ypoints,
                color = color_lipid_center,
                s     = size,
            )
            for i, xy in enumerate(final_frame):
                name, radius = vals["lipids"][i]
                circle = plt.Circle(
                    xy,
                    radius = radius,
                    fill   = False,
                    color  = color_lipid_circle_dict[name],
                )
                axs[row_nr, col_nr].add_artist(circle)

        def plot_protein_points_inside():
            if vals["protein_points_inside_membrane"]:
                xpoints, ypoints = list(zip(*vals["protein_points_inside_membrane"]))
                axs[row_nr, col_nr].scatter(
                    xpoints,
                    ypoints,
                    color=color_protein_points,
                    s=size/2,
                )
        def plot_protein_points_surrounding():
            if vals["protein_points_surrounding_membrane"]:
                xpoints, ypoints = list(zip(*vals["protein_points_surrounding_membrane"]))
                axs[row_nr, col_nr].scatter(
                    xpoints,
                    ypoints,
                    color=color_protein_points,
                    s=size/2,
                )
        def plot_protein_points():
            plot_protein_points_inside()
            plot_protein_points_surrounding()

        def plot_protein_Points_buffered_union_inside():
            if vals["protein_Points_buffered_union_inside_membrane"]:
                shapely.plotting.plot_polygon(
                    vals["protein_Points_buffered_union_inside_membrane"],
                    ax=axs[row_nr, col_nr],
                    add_points=False,
                    color=color_protein_points_buffered,
                )
        def plot_protein_Points_buffered_union_surrounding():
            if vals["protein_Points_buffered_union_surrounding_membrane"]:
                shapely.plotting.plot_polygon(
                    vals["protein_Points_buffered_union_surrounding_membrane"],
                    ax=axs[row_nr, col_nr],
                    add_points=False,
                    color=color_protein_points_buffered,
                )
        def plot_protein_Points_buffered_union():
            plot_protein_Points_buffered_union_inside()
            plot_protein_Points_buffered_union_surrounding()

        def plot_protein_polygon_exterior_points_inside():
            if vals["protein_polygon_exterior_points_inside_membrane"]:
                xpoints, ypoints = list(zip(*vals["protein_polygon_exterior_points_inside_membrane"]))
                axs[row_nr, col_nr].scatter(
                    xpoints,
                    ypoints,
                    color=color_protein_points_exterior,
                    s=size/4,
                )
        def plot_protein_polygon_exterior_points_surrounding():
            if vals["protein_polygon_exterior_points_surrounding_membrane"]:
                xpoints, ypoints = list(zip(*vals["protein_polygon_exterior_points_surrounding_membrane"]))
                axs[row_nr, col_nr].scatter(
                    xpoints,
                    ypoints,
                    color=color_protein_points_exterior,
                    s=size/4,
                )
        def plot_protein_polygon_exterior_points():
            plot_protein_polygon_exterior_points_inside()
            plot_protein_polygon_exterior_points_surrounding()

        def plot_protein_alphashape_inside():
            if vals["protein_ALPHASHAPE_inside_membrane"]:
                shapely.plotting.plot_polygon(
                    vals["protein_ALPHASHAPE_inside_membrane"],
                    ax=axs[row_nr, col_nr],
                    add_points=False,
                    color=color_protein_points_alphashape,
                )
        def plot_protein_alphashape_surrounding():
            if vals["protein_ALPHASHAPE_surrounding_membrane"]:
                shapely.plotting.plot_polygon(
                    vals["protein_ALPHASHAPE_surrounding_membrane"],
                    ax=axs[row_nr, col_nr],
                    add_points=False,
                    color=color_protein_points_alphashape,
                )
        def plot_protein_alphashapes():
            plot_protein_alphashape_inside()
            plot_protein_alphashape_surrounding()

        def plot_protein_ConcaveHulls_Polygon_inside():
            if vals["protein_ConcaveHulls_Polygon_inside_membrane"]:
                shapely.plotting.plot_polygon(
                    shapely.ops.unary_union(vals["protein_ConcaveHulls_Polygon_inside_membrane"]),
                    ax=axs[row_nr, col_nr],
                    add_points=False,
                    color="green",
                )
        def plot_protein_ConcaveHulls_Polygon_surrounding():
            if vals["protein_ConcaveHulls_Polygon_surrounding_membrane"]:
                shapely.plotting.plot_polygon(
                    shapely.ops.unary_union(vals["protein_ConcaveHulls_Polygon_surrounding_membrane"]),
                    ax=axs[row_nr, col_nr],
                    add_points=False,
                    color="green",
                )
        def plot_protein_ConcaveHulls_Polygon():
            plot_protein_ConcaveHulls_Polygon_inside()
            plot_protein_ConcaveHulls_Polygon_surrounding()

        def plot_fixer(plot_minor = True, given_row_nr = False, given_col_nr = False):
            if article_figures:
                plot_fixer_article(plot_minor, given_row_nr, given_col_nr)
            else:
                for ax_row in axs:
                    for ax in ax_row:
                        ### Only plots bins if optimization was run ('bin_size' key:val pair added by optimization)
                        if "bin_size" in vals:
                            xnbins = int(pbcx / vals["bin_size"])
                            ynbins = int(pbcy / vals["bin_size"])
                            ### Finding the actual lengths of the bins
                            xlen, ylen = pbcx / xnbins, pbcy / ynbins

                            ### Minor axes ### Bins
                            if plot_minor:
                                xminors = [round((i*xlen)-(pbcx/2), 3)+np.mean([xmin, xmax]) for i in range(1, xnbins)]
                                yminors = [round((i*ylen)-(pbcy/2), 3)+np.mean([ymin, ymax]) for i in range(1, ynbins)]
                                ax.set_xticks(xminors, labels=["" for _ in xminors], minor=True)
                                ax.set_yticks(yminors, labels=["" for _ in yminors], minor=True)
                                ax.grid(True, which="minor", color="grey")

                        ### Major axes
                    #     ax.set_xticks(range(xmin, xmax+1, pbcx//10), minor=False)
                    #     ax.set_yticks(range(ymin, ymax+1, pbcy//10), minor=False)
                        ax.grid(True, which="major", color="black")

                        ### Axis limits
                    #     ax.set_xlim(xmin-pbcx/100, xmax+pbcx/100)
                    #     ax.set_ylim(ymin-pbcy/100, ymax+pbcy/100)
                        ax.set_xlim(xmin, xmax)
                        ax.set_ylim(ymin, ymax)


                        adj_vals = adjuster(figure=figs, ax=ax, l=l, b=b, r=r, t=t)
                        # print(ax._position)
                        plt.subplots_adjust(
                            left   = adj_vals["l"],
                            bottom = adj_vals["b"],
                            right  = adj_vals["r"],
                            top    = adj_vals["t"],
                    #                 wspace=0.1,
                    #                 hspace=0.0,
                        )
        def plot_fixer_article(plot_minor = True, given_row_nr = False, given_col_nr = False):
            for ax_row in axs:
                for ax in ax_row:
                    xnbins = int(pbcx / vals["bin_size"])
                    ynbins = int(pbcy / vals["bin_size"])
                    ### Finding the actual lengths of the bins
                    xlen, ylen = pbcx / xnbins, pbcy / ynbins

                    ### Minor axes ### Bins
                    if plot_minor:
                        xminors = [round((i*xlen)-(pbcx/2), 3)+np.mean([xmin, xmax]) for i in range(1, xnbins)]
                        yminors = [round((i*ylen)-(pbcy/2), 3)+np.mean([ymin, ymax]) for i in range(1, ynbins)]
                        ax.set_xticks(xminors, labels=["" for _ in xminors], minor=True)
                        ax.set_yticks(yminors, labels=["" for _ in yminors], minor=True)
                        ax.grid(True, which="minor", color="grey")

                    ### Major axes
        #             ax.set_xticks(range(xmin, xmax+1, pbcx//10), minor=False)
        #             ax.set_yticks(range(ymin, ymax+1, pbcy//10), minor=False)
                    ax.set_xticks([], minor=False)
                    ax.set_yticks([], minor=False)
        #             ax.grid(True, which="major", color="black")

                    ### Axis limits
                #     ax.set_xlim(xmin-pbcx/100, xmax+pbcx/100)
                #     ax.set_ylim(ymin-pbcy/100, ymax+pbcy/100)
                    ax.set_xlim(xmin, xmax)
                    ax.set_ylim(ymin, ymax)


        #             l, b, r, t = 20, 20, 20, 20
                    adj_vals = adjuster(figure=figs, ax=ax, l=l, b=b, r=r, t=t)
                    # print(ax._position)
                    plt.subplots_adjust(
                        left   = adj_vals["l"],
                        bottom = adj_vals["b"],
                        right  = adj_vals["r"],
                        top    = adj_vals["t"],
                #                 wspace=0.1,
                #                 hspace=0.0,
                    )
        def plot_fixer_old(plot_minor = True, given_row_nr = False, given_col_nr = False):
            if given_row_nr is not False:
                row_nr = given_row_nr
            else:
                row_nr = 0
            if given_col_nr is not False:
                col_nr = given_col_nr
            else:
                col_nr = 0

            xnbins = int(pbcx / vals["bin_size"])
            ynbins = int(pbcy / vals["bin_size"])
            ### Finding the actual lengths of the bins
            xlen, ylen = pbcx / xnbins, pbcy / ynbins

            ### Minor axes ### Bins
            if plot_minor:
                xminors = [round((i*xlen)-(pbcx/2), 3)+np.mean([xmin, xmax]) for i in range(1, xnbins)]
                yminors = [round((i*ylen)-(pbcy/2), 3)+np.mean([ymin, ymax]) for i in range(1, ynbins)]
                axs[row_nr, col_nr].set_xticks(xminors, labels=["" for _ in xminors], minor=True)
                axs[row_nr, col_nr].set_yticks(yminors, labels=["" for _ in yminors], minor=True)
                axs[row_nr, col_nr].grid(True, which="minor", color="grey")

            ### Major axes
        #     axs[row_nr, col_nr].set_xticks(range(xmin, xmax+1, pbcx//10), minor=False)
        #     axs[row_nr, col_nr].set_yticks(range(ymin, ymax+1, pbcy//10), minor=False)
            axs[row_nr, col_nr].grid(True, which="major", color="black")

            ### Axis limits
        #     axs[row_nr, col_nr].set_xlim(xmin-pbcx/100, xmax+pbcx/100)
        #     axs[row_nr, col_nr].set_ylim(ymin-pbcy/100, ymax+pbcy/100)
            axs[row_nr, col_nr].set_xlim(xmin, xmax)
            axs[row_nr, col_nr].set_ylim(ymin, ymax)


            adj_vals = adjuster(figure=figs, ax=axs[row_nr, col_nr], l=l, b=b, r=r, t=t)
            # print(axs[row_nr, col_nr]._position)
            plt.subplots_adjust(
                left   = adj_vals["l"],
                bottom = adj_vals["b"],
                right  = adj_vals["r"],
                top    = adj_vals["t"],
        #                 wspace=0.1,
        #                 hspace=0.0,
            )

        main_dir = os.path.join(self.PLOT_cmd["path"], "grid_plots")
        os.makedirs(main_dir, exist_ok=True)

        for pi, (keys, vals) in enumerate(list(self.plot_data.items())[:]):
            if pi != 0:
                self.print_term("", verbose=2)
            # subleaflet_dict = self.MEMBRANES[keys[0]]["leaflets"][keys[1]]["subleaflets"][(keys[2], keys[3], keys[4])]

            ### The 'key' is a tuple with 5 values (membrane number, leaflet designation, subleaflet x-index, subleaflet y-index, subleaflet general index)
            self.print_term("Making plot series", pi,                                         spaces=0, verbose=2)

            self.print_term("Plot details:      ",                                            spaces=1, verbose=2)
            self.print_term("Membrane nr:       ", keys[0],                                   spaces=2, verbose=2)
            self.print_term("Leaflet:           ", " ".join(keys[1].capitalize().split("_")), spaces=2, verbose=2)
            self.print_term("Subleaflet x-index:", keys[2],                                   spaces=2, verbose=2)
            self.print_term("Subleaflet y-index:", keys[3],                                   spaces=2, verbose=2)
            self.print_term("Subleaflet nr:     ", keys[4],                                   spaces=2, verbose=2)
            
            plot_dir_name = "_".join([str(pi)] + [str(key) for key in keys])
            self.print_term("Plot directory name:", plot_dir_name, spaces=2, verbose=2)

            key_dir       = os.path.join(main_dir, plot_dir_name)
            key_pngs_dir  = os.path.join(key_dir, "optimization")
            gif_file      = os.path.join(key_dir, "animation_of_optimization.gif")
            mp4_file      = os.path.join(key_dir, "animation_of_optimization.mp4")
            png_file_pre  = os.path.join(key_dir, "grid")
            png_main_name = os.path.join(key_pngs_dir, "grid")
            os.makedirs(key_dir, exist_ok=True)

            if vals["grid_maker_grouping_algorithm"] == "3D_matrix":
                key_lipid_steps_dir  = os.path.join(key_dir, "lipid_steps")
                png_lipid_steps_name = os.path.join(key_lipid_steps_dir, "grid")
                os.makedirs(key_lipid_steps_dir, exist_ok=True)

            xmin, xmax = vals["xdims_original"]
            ymin, ymax = vals["ydims_original"]

            pbcx = xmax - xmin
            pbcy = ymax - ymin

            self.print_term("Leaflet details:",                                                                       spaces=1, verbose=2)
            self.print_term("X- and Y-axial bounds:", round(xmin, 3), round(xmax, 3), round(ymin, 3), round(ymax, 3), spaces=2, verbose=2)
            self.print_term("X- and Y- lengths:    ", round(pbcx, 3), round(pbcy, 3),                                 spaces=2, verbose=2)
            # print("row_length", "col_height", round(ncols * (pbcx/pbc_mod), 3), round(nrows * (pbcy/pbc_mod), 3))

            self.print_term("Optimization details:", spaces=1, verbose=2)
            if "step" in vals:
                self.print_term("Number of steps:      ", vals["step"],                  spaces=2, verbose=2)
                self.print_term("Optimization run time:", round(vals["steps_time"], 4),  spaces=2, verbose=2)
                self.print_term("Last push value:      ", round(vals["max_push"], 4),    spaces=2, verbose=2)
                self.print_term("Lipid push multiplier:", vals["lipid_push_multiplier"], spaces=2, verbose=2)
                self.print_term("Edge push multiplier: ", vals["edge_push_multiplier"],  spaces=2, verbose=2)
                os.makedirs(key_pngs_dir, exist_ok=True)
            else:
                self.print_term("No optimization was run",                               spaces=2, verbose=2)

            converter  = Matplotlib(dpi=dpi)
            gif_writer = GIF(gif_file, duration = 0.1)
            mp4_writer = Video(mp4_file, fps = 8)

            def protein_vals_checker(key_prefix):
                keys = [key_prefix + "_" + string for string in ["inside_membrane", "surrounding_membrane"]]
                return any([vals[key] for key in keys if key in vals])

            self.print_term("General steps:", spaces=1, verbose=2, end=" ")
            frame_counter = 0
            ################################### First empty image
            if True:
                plot_name = "Empty"

                frame_nr = str(frame_counter)
                self.print_term(frame_nr, verbose=2, end=" ")
                figs, axs = make_figure()

                plot_fixer(plot_minor = False)
                frame = frame_nr
                plt.savefig(png_file_pre + "_" + str(frame) + "_" + plot_name + ".png")
                plt.close()
                frame_counter += 1

            ################################### All protein points
            if protein_vals_checker("protein"+"_"+"points"):
                plot_name = "PP"

                frame_nr = str(frame_counter)
                self.print_term(frame_nr, verbose=2, end=" ")
                figs, axs = make_figure()

                plot_protein_points()

                plot_fixer(plot_minor = False)
                frame = frame_nr
                plt.savefig(png_file_pre + "_" + str(frame) + "_" + plot_name + ".png")
                plt.close()
                frame_counter += 1

            ################################### Protein points buffered
            if protein_vals_checker("protein"+"_"+"Points_buffered_union"):
                plot_name = "PPB"

                frame_nr = str(frame_counter)
                self.print_term(frame_nr, verbose=2, end=" ")
                figs, axs = make_figure()

                plot_protein_Points_buffered_union()

                plot_fixer(plot_minor = False)
                frame = frame_nr
                plt.savefig(png_file_pre + "_" + str(frame) + "_" + plot_name + ".png")
                plt.close()
                frame_counter += 1

            ################################### Exterior protein points
            if protein_vals_checker("protein"+"_"+"polygon_exterior_points"):
                plot_name = "PXP"

                frame_nr = str(frame_counter)
                self.print_term(frame_nr, verbose=2, end=" ")
                figs, axs = make_figure()

                plot_protein_polygon_exterior_points()

                plot_fixer(plot_minor = False)
                frame = frame_nr
                plt.savefig(png_file_pre + "_" + str(frame) + "_" + plot_name + ".png")
                plt.close()
                frame_counter += 1

            ################################### ConcaveHulls
            if protein_vals_checker("protein"+"_"+"ConcaveHulls_Polygon"):
                plot_name = "PCH_PXP"

                frame_nr = str(frame_counter)
                self.print_term(frame_nr, verbose=2, end=" ")
                figs, axs = make_figure()

                plot_protein_polygon_exterior_points()
                plot_protein_ConcaveHulls_Polygon()

                plot_fixer(plot_minor = False)
                frame = frame_nr
                plt.savefig(png_file_pre + "_" + str(frame) + "_" + plot_name + ".png")
                plt.close()
                frame_counter += 1

            ################################### Original BBOX before being split and outer protein points
            if "original_holed_bbox" in vals and vals["original_holed_bbox"]:
                plot_name = "OriginalBBOX_PXP"

                frame_nr = str(frame_counter)
                self.print_term(frame_nr, verbose=2, end=" ")
                figs, axs = make_figure()

                plot_bbox_original()
                plot_protein_polygon_exterior_points()

                plot_fixer(plot_minor = False)
                frame = frame_nr
                plt.savefig(png_file_pre + "_" + str(frame) + "_" + plot_name + ".png")
                plt.close()
                frame_counter += 1

            ################################### Outer protein points, concave hulls and BBOX
            if vals["holed_bbox"] and protein_vals_checker("protein"+"_"+"polygon_exterior_points"):
                plot_name = "BBOX_PCH_PXP"

                frame_nr = str(frame_counter)
                self.print_term(frame_nr, verbose=2, end=" ")
                figs, axs = make_figure()

                plot_bbox()
                plot_protein_polygon_exterior_points()
                plot_protein_ConcaveHulls_Polygon()

                plot_fixer(plot_minor = False)
                frame = frame_nr
                plt.savefig(png_file_pre + "_" + str(frame) + "_" + plot_name + ".png")
                plt.close()
                frame_counter += 1

            ################################### BBOX and protein points
            if vals["holed_bbox"] and protein_vals_checker("protein"+"_"+"points"):
                plot_name = "BBOX_PP"

                frame_nr = str(frame_counter)
                self.print_term(frame_nr, verbose=2, end=" ")
                figs, axs = make_figure()

                plot_bbox()
                plot_protein_points()

                plot_fixer(plot_minor = False)
                frame = frame_nr
                plt.savefig(png_file_pre + "_" + str(frame) + "_" + plot_name + ".png")
                plt.close()
                frame_counter += 1

            ################################### Just the BBOX
            if vals["holed_bbox"]:
                plot_name = "BBOX"

                frame_nr = str(frame_counter)
                self.print_term(frame_nr, verbose=2, end=" ")
                figs, axs = make_figure()

                plot_bbox()

                plot_fixer(plot_minor = False)
                frame = frame_nr
                plt.savefig(png_file_pre + "_" + str(frame) + "_" + plot_name + ".png")
                plt.close()
                frame_counter += 1

            ################################### All LineStrings iterations and BBOX
            if self.PLOT_cmd["plot_initial_placement"] and vals["grid_maker_grouping_algorithm"] == "lines" and vals["LineStringsInfo"] and plot_AllLineStringsIterations:
                for i in range(len(vals["LineStringsInfo"])):
                    plot_name = "_".join([
                        "BBOX",
                        "IterationLS" + str(i),
                        "PreOverlap",
                        "xlines" + str(vals["LineStringsInfo"][i]["xlines"]),
                        "ylines" + str(vals["LineStringsInfo"][i]["ylines"]),
                    ])

                    frame_nr = str(frame_counter)
                    self.print_term(frame_nr, verbose=2, end=" ")
                    figs, axs = make_figure()

                    plot_bbox(i)
                    if vals["LineStringsInfo"][i]["grid_method"] == "LineStrings":
                        plot_LineStrings(i)
                        plot_LineStringsXLines(i)
                    elif vals["LineStringsInfo"][i]["grid_method"] == "2D_Grid":
                        plot_GridPoints(i)

                    if "nlipids" in vals["LineStringsInfo"][i]:
                        nlipids = vals["LineStringsInfo"][i]["nlipids"]
                    else:
                        nlipids = len(vals["lipids"])

                    axs[row_nr, col_nr].set_title(
                        ", ".join([
                            "Lipids: " + str(nlipids),
                            "xlines: " + str(vals["LineStringsInfo"][i]["xlines"]),
                            "ylines: " + str(vals["LineStringsInfo"][i]["ylines"]),
                            "Grid points: " + str(vals["LineStringsInfo"][i]["ngridpoints_initial"]),
                        ]),
                        fontsize=8,
                    )

                    plot_fixer(plot_minor = False)
                    frame = frame_nr
                    plt.savefig(png_file_pre + "_" + str(frame) + "_" + plot_name + ".png")
                    plt.close()
                    frame_counter += 1

                    plot_name = "_".join([
                        "BBOX",
                        "IterationLS" + str(i),
                        "PostOverlap",
                        "xlines" + str(vals["LineStringsInfo"][i]["xlines"]),
                        "ylines" + str(vals["LineStringsInfo"][i]["ylines"]),
                    ])

                    frame_nr = str(frame_counter)
                    self.print_term(frame_nr, verbose=2, end=" ")
                    figs, axs = make_figure()

                    plot_bbox(i)
                    if vals["LineStringsInfo"][i]["grid_method"] == "LineStrings":
                        plot_LineStringsContained(i)
                    elif vals["LineStringsInfo"][i]["grid_method"] == "2D_Grid":
                        plot_GridPoints_Contained(i)

                    if "nlipids" in vals["LineStringsInfo"][i]:
                        nlipids = vals["LineStringsInfo"][i]["nlipids"]
                    else:
                        nlipids = len(vals["lipids"])

                    axs[row_nr, col_nr].set_title(
                        ", ".join([
                            "Lipids: " + str(nlipids),
                            "xlines: " + str(vals["LineStringsInfo"][i]["xlines"]),
                            "ylines: " + str(vals["LineStringsInfo"][i]["ylines"]),
                            "Grid points: " + str(vals["LineStringsInfo"][i]["ngridpoints_final"]),
                        ]),
                        fontsize=8,
                    )

                    plot_fixer(plot_minor = False)
                    frame = frame_nr
                    plt.savefig(png_file_pre + "_" + str(frame) + "_" + plot_name + ".png")
                    plt.close()
                    frame_counter += 1

            ################################### LineStrings and BBOX
            if self.PLOT_cmd["plot_initial_placement"] and vals["grid_maker_grouping_algorithm"] == "lines" and vals["LineStringsInfo"]:
                plot_name = "BBOX_FinalGridLS"

                frame_nr = str(frame_counter)
                self.print_term(frame_nr, verbose=2, end=" ")
                figs, axs = make_figure()

                plot_bbox(i)
                if vals["LineStringsInfo"][-1]["grid_method"] == "LineStrings":
                    plot_LineStringsXLines(-1)
                    plot_LineStrings(-1)
                elif vals["LineStringsInfo"][-1]["grid_method"] == "2D_Grid":
                    plot_GridPoints(-1)

                plot_fixer(plot_minor = False)
                frame = frame_nr
                plt.savefig(png_file_pre + "_" + str(frame) + "_" + plot_name + ".png")
                plt.close()
                frame_counter += 1

            ################################### LineStrings and BBOX
            if self.PLOT_cmd["plot_initial_placement"] and vals["grid_maker_grouping_algorithm"] == "lines" and vals["LineStringsInfo"]:
                plot_name = "BBOX_FinalLS"

                frame_nr = str(frame_counter)
                self.print_term(frame_nr, verbose=2, end=" ")
                figs, axs = make_figure()

                plot_bbox(i)
                if vals["LineStringsInfo"][-1]["grid_method"] == "LineStrings":
                    plot_LineStrings(-1)
                elif vals["LineStringsInfo"][-1]["grid_method"] == "2D_Grid":
                    plot_GridPoints(-1)

                plot_fixer(plot_minor = False)
                frame = frame_nr
                plt.savefig(png_file_pre + "_" + str(frame) + "_" + plot_name + ".png")
                plt.close()
                frame_counter += 1

            ################################### Overlapping LineStrings, normal BBOX and reduced BBOX
            if self.PLOT_cmd["plot_initial_placement"] and vals["grid_maker_grouping_algorithm"] == "lines" and vals["LineStringsInfo"]:
#                     c1 = (not all([o==c for o, c in zip(vals["LineStringsInfo"][-1]["LineStringsOverlapping"], vals["LineStringsInfo"][-1]["LineStringsContained"])]))
#                     c2 = (("bbox_polygon_BufferedForLineStrings" in vals and vals["bbox_polygon_BufferedForLineStrings"]) or ("bbox_polygon_BufferedForLineStrings" in vals["LineStringsInfo"][-1] and vals["LineStringsInfo"][-1]["bbox_polygon_BufferedForLineStrings"]))
#                     print("c1, c2", c1, c2)
                c1, c2 = True, True
                if all([c1, c2]):
                    plot_name = "BBOX_LSBufferedBBOX_OverlappingLS"

                    frame_nr = str(frame_counter)
                    self.print_term(frame_nr, verbose=2, end=" ")
                    figs, axs = make_figure()

                    plot_bbox(-1)
                    plot_bbox_buffered_for_LineStrings(-1)
                    plot_LineStringsOverlapping(-1)

                    plot_fixer(plot_minor = False)
                    frame = frame_nr
                    plt.savefig(png_file_pre + "_" + str(frame) + "_" + plot_name + ".png")
                    plt.close()
                    frame_counter += 1

            ################################### Contained LineStrings, normal BBOX and reduced BBOX
            if self.PLOT_cmd["plot_initial_placement"] and vals["grid_maker_grouping_algorithm"] == "lines" and vals["LineStringsInfo"] and (("bbox_polygon_BufferedForLineStrings" in vals and vals["bbox_polygon_BufferedForLineStrings"]) or ("bbox_polygon_BufferedForLineStrings" in vals["LineStringsInfo"][-1] and vals["LineStringsInfo"][-1]["bbox_polygon_BufferedForLineStrings"])):
                plot_name = "BBOX_LSBufferedBBOX_ContainedLS"

                frame_nr = str(frame_counter)
                self.print_term(frame_nr, verbose=2, end=" ")
                figs, axs = make_figure()

                plot_bbox(i)
                plot_bbox_buffered_for_LineStrings(-1)
                plot_LineStringsContained(-1)

                plot_fixer(plot_minor = False)
                frame = frame_nr
                plt.savefig(png_file_pre + "_" + str(frame) + "_" + plot_name + ".png")
                plt.close()
                frame_counter += 1

            ################################### Contained LineStrings and normal BBOX
            if self.PLOT_cmd["plot_initial_placement"] and vals["grid_maker_grouping_algorithm"] == "lines" and vals["LineStringsInfo"] and (("bbox_polygon_BufferedForLineStrings" in vals and vals["bbox_polygon_BufferedForLineStrings"]) or ("bbox_polygon_BufferedForLineStrings" in vals["LineStringsInfo"][-1] and vals["LineStringsInfo"][-1]["bbox_polygon_BufferedForLineStrings"])):
                plot_name = "BBOX_ContainedLS"

                frame_nr = str(frame_counter)
                self.print_term(frame_nr, verbose=2, end=" ")
                figs, axs = make_figure()

                plot_bbox(i)
                plot_LineStringsContained(-1)

                plot_fixer(plot_minor = False)
                frame = frame_nr
                plt.savefig(png_file_pre + "_" + str(frame) + "_" + plot_name + ".png")
                plt.close()
                frame_counter += 1

            ################################### BBOX and intial lipid placements pre random
            if self.PLOT_cmd["plot_initial_placement"] and "grid_points_no_random" in vals and vals["grid_points_no_random"]:
                plot_name = "BBOX_PreRandomLP"

                frame_nr = str(frame_counter)
                self.print_term(frame_nr, verbose=2, end=" ")
                figs, axs = make_figure()

                plot_bbox()
                plot_PreRandomLipidPlacements()

                plot_fixer(plot_minor = False)
                frame = frame_nr
                plt.savefig(png_file_pre + "_" + str(frame) + "_" + plot_name + ".png")
                plt.close()
                frame_counter += 1

            ################################### BBOX and intial lipid placements post random
            if self.PLOT_cmd["plot_initial_placement"] and "POINT_STEPS" in vals and vals["POINT_STEPS"]:
                plot_name = "BBOX_PostRandomLP"

                frame_nr = str(frame_counter)
                self.print_term(frame_nr, verbose=2, end=" ")
                figs, axs = make_figure()

                plot_bbox()
                plot_PostRandomLipidPlacements()

                plot_fixer(plot_minor = False)
                frame = frame_nr
                plt.savefig(png_file_pre + "_" + str(frame) + "_" + plot_name + ".png")
                plt.close()
                frame_counter += 1

            ################################### BBOX and final lipid placements
            if "POINT_STEPS" in vals and vals["POINT_STEPS"]:
                plot_name = "BBOX_FinalLP_PP"

                frame_nr = str(frame_counter)
                self.print_term(frame_nr, verbose=2, end=" ")
                figs, axs = make_figure()

                plot_bbox()
                plot_FinalLipidPlacements()

                plot_fixer(plot_minor = False)
                frame = frame_nr
                plt.savefig(png_file_pre + "_" + str(frame) + "_" + plot_name + ".png")
                plt.close()
                frame_counter += 1            

            ################################### BBOX, final lipid placements and protein
            if "POINT_STEPS" in vals and vals["POINT_STEPS"] and protein_vals_checker("protein"+"_"+"points"):
                plot_name = "BBOX_FinalLP"

                frame_nr = str(frame_counter)
                self.print_term(frame_nr, verbose=2, end=" ")
                figs, axs = make_figure()

                plot_bbox()
                plot_protein_points()
                plot_FinalLipidPlacements()

                plot_fixer(plot_minor = False)
                frame = frame_nr
                plt.savefig(png_file_pre + "_" + str(frame) + "_" + plot_name + ".png")
                plt.close()
                frame_counter += 1            

            self.print_term("", verbose=2)
            ################################### Lipid optimization frames
            if self.PLOT_cmd["plot_optimization_steps"] and "POINT_STEPS" in vals and vals["POINT_STEPS"]:
                with Recorder(gif_writer, converter = converter) as gif_rec:
                    with Recorder(mp4_writer, converter = converter) as mp4_rec:
                        self.print_term("Optimization steps:", spaces=1, verbose=2, end=" ")
                        ### First plot the layout before any steps
                        frame_nr = 0
                        self.print_term(frame_nr, verbose=2, end=" ")
                        figs, axs = plt.subplots(
                            nrows = nrows,
                            ncols = ncols,
                            figsize = (ncols * (pbcx/pbc_mod + (l+r)/dpi), nrows * (pbcy/pbc_mod + (b+t)/dpi)),
                            dpi = dpi,
                            squeeze = False
                        )
                        plot_bbox()

                        xpoints, ypoints = list(zip(*vals["grid_points_no_random"]))
                        axs[row_nr, col_nr].scatter(
                            xpoints,
                            ypoints,
                        #     label=label,
                            color=color_lipid_center,
                        #     alpha=alpha,
                            s=size,
                        )

                        for i, xy in enumerate(vals["grid_points_no_random"]):
                            name, radius = vals["lipids"][i]
                            circle = plt.Circle(
                                xy,
                                radius = radius,
                                fill   = False,
                                color  = color_lipid_circle_dict[name],
                            )
                            axs[row_nr, col_nr].add_artist(circle)

                        plot_fixer(plot_minor = False)

                        if do_zfill:
                            frame = str(frame_nr).zfill(len(str(len(vals["POINT_STEPS"]))))
                        else:
                            frame = str(frame_nr)

                        plt.savefig(png_main_name + "_" + str(frame) + ".png")
                        
                        ### The recorder requires at least 1 frame otherwise it will crash, so the first frame is always recorded. No way to exit the recorder without attempts to make gif/mp4.
                        ### The gif/mp4 files are deleted afterwords if they have been turned off.
                        gif_rec.record(fig = figs)
                        mp4_rec.record(fig = figs)
                        plt.close()

                        ### Then plot the steps
                        for frame_nr, arr in enumerate(vals["POINT_STEPS"], 1):
                            self.print_term(frame_nr, verbose=2, end=" ")
                            figs, axs = plt.subplots(
                                nrows = nrows,
                                ncols = ncols,
                #                 figsize = (ncols * 6 + (l+r)/dpi, nrows * 6 + (b+t)/dpi),
                                figsize = (ncols * (pbcx/pbc_mod + (l+r)/dpi), nrows * (pbcy/pbc_mod + (b+t)/dpi)),
                                dpi = dpi,
                                squeeze = False
                            )
                            plot_bbox()

                            xpoints, ypoints = arr[:,0], arr[:,1]
                            axs[row_nr, col_nr].scatter(
                                xpoints,
                                ypoints,
                            #     label=label,
                                color=color_lipid_center,
                            #     alpha=alpha,
                                s=size,
                            )

                            for i, xy in enumerate(arr):
                                name, radius = vals["lipids"][i]
                                circle = plt.Circle(
                                    xy,
                                    radius = radius,
                                    fill   = False,
                                    color  = color_lipid_circle_dict[name],
                                )
                                axs[row_nr, col_nr].add_artist(circle)

                            plot_fixer(plot_minor = False)

                            if do_zfill:
                                frame = str(frame_nr).zfill(len(str(len(vals["POINT_STEPS"]))))
                            else:
                                frame = str(frame_nr)

                            plt.savefig(png_main_name + "_" + str(frame) + ".png")

                            if frame_nr > 1:
                                if self.PLOT_cmd["make_gif"]:
                                    gif_rec.record(fig = figs)
                                if self.PLOT_cmd["make_mp4"]:
                                    mp4_rec.record(fig = figs)
                            plt.close()
                        
                    if self.PLOT_cmd["make_mp4"]:
                        self.print_term("MP4", end=" ", verbose=2)
                    else:
                        mp4_rec.close()
                        os.remove(mp4_file)

                if self.PLOT_cmd["make_gif"]:
                    self.print_term("GIF", end=" ", verbose=2)
                else:
                    gif_rec.close()
                    os.remove(gif_file)
                
                self.print_term("", verbose=2)
            
            ################################### BBOX and intial lipid placements post random
            if vals["grid_maker_grouping_algorithm"] == "3D_matrix" and "lipid_insertion" in vals and plot_AllLipidInsertions:
                self.print_term("Insertion steps:", spaces=1, verbose=2, end=" ")
                for i in range(len(vals["lipid_insertion"]["lipid_insertion_steps"])):
                    plot_name = "BBOX_LipidInsertionStep" + str(i)

                    frame_nr = str(i)
                    self.print_term(frame_nr, verbose=2, end=" ")
                    figs, axs = make_figure(given_ncols = len(list(dict.fromkeys((list(zip(*vals["lipids"]))[0])))))

                    plot_LipidInsertionStep(i)
                    plot_bbox() # plot bbox after otherwise it is hidden under other stuff

                    plot_fixer(plot_minor = False, given_col_nr = i)
                    frame = frame_nr
                    plt.savefig(png_lipid_steps_name + "_" + str(frame) + "_" + plot_name + ".png")
                    plt.close()
                self.print_term("", verbose=2)
                        
