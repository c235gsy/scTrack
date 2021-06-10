import time
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.covariance import OAS
from sklearn.decomposition import PCA
import sys
from mpl_toolkits.mplot3d import Axes3D
import get_parameters
from utils import *

start = time.time ()

parameters = get_parameters.parse_args_ananlyze_paths ()

expression_file = parameters.expMatrix
cell_type_file = parameters.clusterArray
cell_id_file = parameters.idArray
output = parameters.outputName
path_file = parameters.pathFile
smoothness = parameters.smoothness

expression_matrix = np.loadtxt (expression_file, dtype=float)
cell_type_array = np.loadtxt (cell_type_file, dtype=str)
all_cell_number = len (cell_type_array)
all_cell_types = sorted (list (set (cell_type_array)), key=sort_key)

cell_id_array = None
if cell_id_file == False:
    cell_id_array = np.array (["Cell{}".format (i) for i in range (0, all_cell_number)])
else:
    cell_id_array = np.loadtxt (cell_id_file, dtype=str)

with open (path_file, 'r') as pathFile:
    paths = {l.strip ().split (",")[0]: l.strip ().split (",")[2:] for l in pathFile.readlines ()}

used_cell_types = set ()
for pv in paths.values ():
    used_cell_types.update (set (pv))
unused_cell_types = sorted (list (set (all_cell_types) - used_cell_types), key=sort_key)
used_cell_types = sorted (list (used_cell_types), key=sort_key)

used_cell_array = np.where ([ct in used_cell_types for ct in cell_type_array])[0]
expression_matrix = expression_matrix[used_cell_array]
cell_type_array = cell_type_array[used_cell_array]

colors = get_color_list (n=len (used_cell_types))
color_dic_celltype = {used_cell_types[i]: colors[i] for i in range (len (used_cell_types))}

gene_number = len (expression_matrix[0])
used_cell_number = len (expression_matrix)
used_cell_cluster_number = len (used_cell_types)
print ("Data information :")
print ("gene : ", gene_number)
print ("all cell: ", all_cell_number)
print ("used cell : ", used_cell_number)
print ("used clusters : ", ", ".join (used_cell_types))
print ("unused clusters : ", ", ".join (unused_cell_types))
print ("\nPath information : ")
for k, v in paths.items ():
    print ("{} :  {}".format (k, "-->".join (v)))

n_components_pca = parameters.nComponentsPCA
if n_components_pca == False:
    n_components_pca = min (5 * (used_cell_cluster_number - 1), int(0.01 * gene_number))

n_components_lda = parameters.nComponentsLDA
if n_components_lda == False:
    n_components_lda = min (used_cell_cluster_number - 1, n_components_pca)

if n_components_pca < n_components_lda:
    print ("--nComponentsPCA(-np) should not be less than --nComponentsLDA(-nl)")
    sys.exit (0)

expression_matrix = PCA (n_components=n_components_pca, svd_solver="full").fit_transform (expression_matrix)
print ("Matrix shape after PCA: ", expression_matrix.shape)
oa = OAS(store_precision=False, assume_centered=False)
expression_matrix = LDA (n_components=n_components_lda,
                         covariance_estimator=OAS(store_precision=False, assume_centered=False),
                         solver='eigen').fit_transform (expression_matrix, cell_type_array)
print ("Matrix shape after LDA: ", expression_matrix.shape)


lab = {path_name: {"cells": np.where ([i in paths[path_name] for i in cell_type_array])[0], "curve": None} for path_name
       in paths.keys ()}

for path_name in paths.keys ():
    print ()
    print ("#" * 40)
    print ("#" * 40)
    print (path_name)

    cells_in_path = expression_matrix[lab[path_name]["cells"]]

    cluster_centers = np.array ([expression_matrix[cell_type_array == clu].mean (axis=0) for clu in paths[path_name]])

    print ("\nFitting the principal curve:")

    pcurve_obj, D2_list = pypcurve.get_princurve_curve (points=cells_in_path, start=cluster_centers, cutoff_iter=0.001,
                                                        n_curve_seq=10 * len (cells_in_path), inter_dimension_init=1,
                                                        smoothness=smoothness)

    lab[path_name]["curve"] = pcurve_obj["projection"][pcurve_obj["order"]]
    lab[path_name]["cells"] = lab[path_name]["cells"][pcurve_obj["order"]]

    distance_array = np.sqrt (np.sum (np.square (cells_in_path - pcurve_obj["projection"]), axis=1))
    distance_array = distance_array[pcurve_obj["order"]]
    celltype_in_path = cell_type_array[lab[path_name]["cells"]]
    filter_array = np.array ([True for i in range (len (lab[path_name]["cells"]))])

    # print (sum (filter_array))
    print ("\nFirst filtration: ")
    print ("{} cells in this path before filtering.".format (len (filter_array)))
    for clu in paths[path_name]:
        clu_cells = np.where (celltype_in_path == clu)[0]
        ddd = distance_array[clu_cells]
        Q1, Q3 = np.percentile (ddd, [25, 75])
        IQR = Q3 - Q1
        dis_cut_off = Q3 + 1.5 * IQR
        filter_array[clu_cells[ddd > dis_cut_off]] = False
        # print (sum (filter_array), (1 - sum (filter_array) / len (filter_array)) * 100)
    print ("{} cells in this path after filtering.".format (len (filter_array)))
    print ("{:.4f} % of cells are filtered.".format ((1 - sum (filter_array) / len (filter_array)) * 100))

    filter_cells = lab[path_name]["cells"][filter_array]
    filter_cells_in_path = expression_matrix[filter_cells]

    print ("\nFitting the principal curve:")
    filter_cluster_centers = np.array (
        [expression_matrix[cell_type_array == clu].mean (axis=0) for clu in paths[path_name]])
    filter_pcurve_obj, D2_list = pypcurve.get_princurve_curve (points=filter_cells_in_path,
                                                               start=filter_cluster_centers, cutoff_iter=0.001,
                                                               n_curve_seq=10 * len (cells_in_path),
                                                               inter_dimension_init=1, smoothness=smoothness)

    filter_curve = filter_pcurve_obj["projection"][filter_pcurve_obj["order"]]

    filter_pcurve_obj_all = pypcurve.project_points (points=expression_matrix[lab[path_name]["cells"]],
                                                     curve=filter_curve, n_curve_seq=10 * len (cells_in_path),
                                                     smoothness=smoothness)

    lab[path_name]["cells"] = lab[path_name]["cells"][filter_pcurve_obj_all["order"]]
    lab[path_name]["curve"] = filter_pcurve_obj_all["projection"][filter_pcurve_obj_all["order"]]

    distance_array_2 = np.sqrt (
        np.sum (np.square (expression_matrix[lab[path_name]["cells"]] - lab[path_name]["curve"]), axis=1))
    celltype_in_path = cell_type_array[lab[path_name]["cells"]]
    filter_array = np.array ([True for i in range (len (lab[path_name]["cells"]))])

    print ("\nSecond filtration: ")
    print ("{} cells in this path before filtering.".format (len (filter_array)))
    for clu in paths[path_name]:
        clu_cells = np.where (celltype_in_path == clu)[0]
        ddd = distance_array_2[clu_cells]
        Q1, Q3 = np.percentile (ddd, [25, 75])
        IQR = Q3 - Q1
        dis_cut_off = Q3 + 1.5 * IQR
        filter_array[clu_cells[ddd > dis_cut_off]] = False
        # print (sum (filter_array), (1 - sum (filter_array) / len (filter_array)) * 100)
    print ("{} cells in this path after filtering.".format (len (filter_array)))
    print ("{:.4f} % of cells are filtered.".format ((1 - sum (filter_array) / len (filter_array)) * 100))

    lab[path_name]["filter"] = filter_array

if parameters.combineBeginnings:
    print ("\nCombining beginnings")
    lab = combine_paths (the_lab=lab, the_paths=paths, cell_type=cell_type_array, smoothness=smoothness)

if parameters.combineDestinations:
    print ("\nCombining destinations")
    for k in lab.keys ():
        paths[k] = paths[k][::-1]
        lab[k]["cells"] = lab[k]["cells"][::-1]
        lab[k]["curve"] = lab[k]["curve"][::-1]
    lab = combine_paths (the_lab=lab, the_paths=paths, cell_type=cell_type_array, smoothness=smoothness)
    for k in lab.keys ():
        paths[k] = paths[k][::-1]
        lab[k]["cells"] = lab[k]["cells"][::-1]
        lab[k]["curve"] = lab[k]["curve"][::-1]

path_names = list (paths.keys ())
pseudotime_table_filter = pd.DataFrame (np.full ([len (cell_id_array), len (paths) + 1], np.nan),
                                        columns=["Cell_id"] + path_names)
pseudotime_table_filter["Cell_id"] = cell_id_array

max_pseudotime_filter = 0
for p in range (len (path_names)):
    path_name = path_names[p]
    lab[path_name]["pseudotime_filter"] = curve_2_pseudotime (lab[path_name]["curve"][lab[path_name]["filter"]])
    max_pseudotime_filter = max (max_pseudotime_filter, lab[path_name]["pseudotime_filter"].max ())

for p in range (len (path_names)):
    path_name = path_names[p]
    lab[path_name]["pseudotime_filter"] = (lab[path_name]["pseudotime_filter"] / max_pseudotime_filter) * 100
    pseudotime_table_filter.iloc[used_cell_array[lab[path_name]["cells"][lab[path_name]["filter"]]], p + 1] = \
        lab[path_name]["pseudotime_filter"]

pseudotime_table_filter.to_csv ("{}_filter_pseudotime.csv".format (output), na_rep="NAN", index=False)
del pseudotime_table_filter["Cell_id"]

pseudotime_array_for_plot = np.array (pseudotime_table_filter.iloc[used_cell_array].mean (axis=1))

un_filtered_cell = np.array ([])
for b in lab.values ():
    un_filtered_cell = np.hstack ((un_filtered_cell, b["cells"][b["filter"]]))

un_filtered_cell = np.array (list (set (un_filtered_cell)), dtype=int)
un_filtered_cell_type_array = cell_type_array[un_filtered_cell]
un_filtered_expression_matrix = expression_matrix[un_filtered_cell]
un_filtered_pseudotime_array_for_plot = pseudotime_array_for_plot[un_filtered_cell]

#######################################################################################################################
# plot
#######################################################################################################################
points_and_lines = np.concatenate ([un_filtered_expression_matrix] + [v["curve"][v["filter"]] for v in lab.values ()],
                                   axis=0)
PCA_obj_for_plot = PCA (svd_solver="full")
PCA_obj_for_plot.fit (points_and_lines)
points = PCA_obj_for_plot.transform (un_filtered_expression_matrix)

#######################################################################################################################
# output files
#######################################################################################################################

print ("\nWriting data for plotting ......")

output_points_table = pd.DataFrame (np.full ([len (cell_id_array), n_components_lda + 1], np.nan),
                                    columns=["Cell_id"] + [
                                        "PCA{}_LDA{}_PC{}".format (n_components_pca, n_components_lda, i) for i in
                                        range (1, n_components_lda + 1)])
output_points_table["Cell_id"] = cell_id_array
output_points_table.iloc[used_cell_array[un_filtered_cell], 1:] = points
output_points_table.to_csv ("{}_point_for_plot.csv".format (output), na_rep="NAN", index=False)

for p in range (len (path_names)):
    path_name = path_names[p]
    output_curve_table = pd.DataFrame (np.full ([len (cell_id_array), n_components_lda + 1], np.nan),
                                       columns=["Cell_id"] + [
                                           "PCA{}_LDA{}_PC{}".format (n_components_pca, n_components_lda, i) for i in
                                           range (1, n_components_lda + 1)])
    lab_this_path = lab[path_name]
    curve = PCA_obj_for_plot.transform (lab_this_path["curve"][lab_this_path["filter"]])
    output_curve_table["Cell_id"] = cell_id_array
    output_curve_table.iloc[used_cell_array[lab_this_path["cells"][lab_this_path["filter"]]], 1:] = curve
    output_curve_table.to_csv ("{}_curve{}_for_plot.csv".format (output, p + 1), na_rep="NAN", index=False)

if parameters.plot2D:
    print ("Plotting 2D figures ......")
    plt.figure (figsize=(8, 6))
    for ide in used_cell_types:
        cells_in_type = points[un_filtered_cell_type_array == ide].T
        color = color_dic_celltype[ide]
        plt.scatter (cells_in_type[0], cells_in_type[1], color=color, s=2, alpha=0.5)

    leg = plt.legend (used_cell_types, bbox_to_anchor=(1.0, 0.5), loc=6, borderaxespad=0, ncol=1, frameon=False,
                      fontsize=15, markerscale=5, labelspacing=0.7)

    for lh in leg.legendHandles:
        lh.set_alpha (1)

    for p in lab.values ():
        line = PCA_obj_for_plot.transform (p["curve"][p["filter"]]).T
        plt.plot (line[0], line[1], color="black", linewidth=3)

    plt.xlabel ("PCA{}_LDA{}_PC1".format (n_components_pca, n_components_lda))
    plt.ylabel ("PCA{}_LDA{}_PC2".format (n_components_pca, n_components_lda))

    plt.savefig ("{}_plot.pdf".format (output), bbox_inches='tight')
    plt.clf ()

    plt.figure (figsize=(10, 6))
    cm = plt.cm.get_cmap ('jet')
    plt.scatter (points[:, 0], points[:, 1], c=un_filtered_pseudotime_array_for_plot, s=2, alpha=0.5, vmin=0, vmax=100,
                 cmap=cm)
    for p in lab.values ():
        line = PCA_obj_for_plot.transform (p["curve"][p["filter"]]).T
        plt.plot (line[0], line[1], color="black", linewidth=3)

    color_bar = plt.colorbar (label="pseudotime")
    color_bar.set_alpha (1)
    color_bar.draw_all ()

    plt.xlabel ("PCA{}_LDA{}_PC1".format (n_components_pca, n_components_lda))
    plt.ylabel ("PCA{}_LDA{}_PC2".format (n_components_pca, n_components_lda))
    plt.savefig ("{}_plot_pseudotime.pdf".format (output), bbox_inches='tight')
    plt.clf ()

if parameters.plot3D and n_components_lda < 3:
    print ("There is only 2 dimensions, program cannnot plot 3D figures.")

if parameters.plot3D and n_components_lda >= 3:
    print ("Plotting 3D figures ......")
    fig = plt.figure (figsize=(8, 6))
    for azi in [15, 30, 45, 60, 75, -15, -30, -45, -60, -75, 105, 120, 135, 150, 165, -105, -120, -135, -150, -165]:
        for ele in [30, 45, 60]:
            ax = Axes3D (fig)
            ax.view_init (azim=azi, elev=ele)
            for ide in used_cell_types:
                cells_in_type = points[un_filtered_cell_type_array == ide].T
                color = color_dic_celltype[ide]
                ax.scatter (cells_in_type[0], cells_in_type[1], cells_in_type[2], color=color, s=2, alpha=0.5)

            leg = ax.legend (used_cell_types, bbox_to_anchor=(1, 0.5), loc=6, borderaxespad=0, ncol=1, frameon=False,
                             fontsize=15, markerscale=5, labelspacing=0.7)

            for lh in leg.legendHandles:
                lh.set_alpha (1)

            for p in lab.values ():
                line = PCA_obj_for_plot.transform (p["curve"][p["filter"]]).T
                ax.plot (line[0], line[1], line[2], color="black", linewidth=3, alpha=1)

            ax.set_xlabel ("PCA{}_LDA{}_PC1".format (n_components_pca, n_components_lda))
            ax.set_ylabel ("PCA{}_LDA{}_PC2".format (n_components_pca, n_components_lda))
            ax.set_zlabel ("PCA{}_LDA{}_PC3".format (n_components_pca, n_components_lda))

            plt.savefig ("{}_plot3D_e{}_a{}.pdf".format (output, ele, azi), bbox_inches='tight')
            fig.clf ()

    fig = plt.figure (figsize=(10, 6))
    cm = plt.cm.get_cmap ('jet')
    for azi in [15, 30, 45, 60, 75, -15, -30, -45, -60, -75, 105, 120, 135, 150, 165, -105, -120, -135, -150, -165]:
        for ele in [30, 45, 60]:
            ax = Axes3D (fig)
            ax.view_init (azim=azi, elev=ele)
            im = ax.scatter (points[:, 0], points[:, 1], points[:, 2], c=un_filtered_pseudotime_array_for_plot, s=2,
                             alpha=0.5, vmin=0, vmax=100, cmap=cm)

            color_bar = plt.colorbar (im, ax=ax, label="pseudotime")
            color_bar.set_alpha (1)
            color_bar.draw_all ()

            for p in lab.values ():
                line = PCA_obj_for_plot.transform (p["curve"][p["filter"]]).T
                ax.plot (line[0], line[1], line[2], color="black", linewidth=3, alpha=1)

            ax.set_xlabel ("PCA{}_LDA{}_PC1".format (n_components_pca, n_components_lda))
            ax.set_ylabel ("PCA{}_LDA{}_PC2".format (n_components_pca, n_components_lda))
            ax.set_zlabel ("PCA{}_LDA{}_PC3".format (n_components_pca, n_components_lda))

            plt.savefig ("{}_plot3D_pseudotime_e{}_a{}.pdf".format (output, ele, azi), bbox_inches='tight')
            fig.clf ()
