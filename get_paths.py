import time
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.covariance import OAS
from sklearn.decomposition import PCA
from utils import *
import sys
import get_parameters


start = time.time ()
parameters = get_parameters.parse_args_get_paths()

expression_file = parameters.expMatrix
cell_type_file = parameters.clusterArray
output_file = "{}_paths.txt".format(parameters.output)
all_output_file = "{}_all_paths.txt".format(parameters.output)

expression_matrix = np.loadtxt (expression_file, dtype=float)
cell_type_array = np.loadtxt (cell_type_file, dtype=str)
all_cell_types = sorted (list (set (cell_type_array)), key=sort_key)
cell_subtype_array = np.array ([None] * len (cell_type_array))

gene_number = len(expression_matrix[0])
cell_number = len(expression_matrix)
cell_cluster_number = len(all_cell_types)
print ("Data information :")
print ("gene : ", gene_number)
print ("cell : ", cell_number)
print ("clusters : ", all_cell_types)

n_components_pca = parameters.nComponentsPCA
if n_components_pca == False:
    n_components_pca = min(5*(cell_cluster_number-1), int(0.01*gene_number))

n_components_lda = parameters.nComponentsLDA
if n_components_lda == False:
    n_components_lda = min(cell_cluster_number - 1, n_components_pca)

if n_components_pca < n_components_lda:
    print("--nComponentsPCA(-np) should not be less than --nComponentsLDA(-nl)")
    sys.exit (0)

expression_matrix = PCA (n_components=n_components_pca, svd_solver="full").fit_transform (expression_matrix)
print("Matrix shape after PCA: ", expression_matrix.shape)
oa = OAS(store_precision=False, assume_centered=False)
expression_matrix = LDA (n_components=n_components_lda,
                         covariance_estimator=OAS(store_precision=False, assume_centered=False),
                         solver='eigen').fit_transform (expression_matrix, cell_type_array)
print("Matrix shape after LDA: ", expression_matrix.shape)

average_size_subclusters = parameters.sizeSubcluster
celltype2subtype = {}
for celltype in all_cell_types:
    idx = np.where (cell_type_array == celltype)[0]
    n_clu = int (len (idx) / average_size_subclusters) + 1
    cells_of_this_celltype = expression_matrix[idx]
    predict_of_subtype = kmeans (cells_of_this_celltype, n_clusters=n_clu)
    subcelltype = np.array (["{}*SustechJinLab*{}".format (celltype, p) for p in predict_of_subtype.labels_])
    celltype2subtype[celltype] = sorted (list (set (subcelltype)), key=sort_key)
    cell_subtype_array[idx] = subcelltype

all_cell_subtypes = sorted (list (set (cell_subtype_array)), key=sort_key)

subcelltype_centers = np.array (
    [expression_matrix[cell_subtype_array == subtype].mean (axis=0) for subtype in all_cell_subtypes])
subcelltype_centers_dis_matrix = get_distance_matricx (subcelltype_centers, all_cell_subtypes)

celltype_centers = np.array (
    [expression_matrix[cell_type_array == cell_type].mean (axis=0) for cell_type in all_cell_types])
celltype_centers_dis_matrix = get_distance_matricx (celltype_centers, all_cell_types)

celltype_centers_bf_matrix = get_bayes_factor_matrix(expression_matrix, cell_type_array, all_cell_types,
                                                     num_calculation=parameters.nCalculationBayesFactor)
print("Bayes Factor Matirx: ")
print(celltype_centers_bf_matrix)
celltype_centers_bf_matrix.rename (index=lambda x: "to_" + x).to_csv("{}_BF.txt".format(output_file), sep="\t",
                                                                     header=True, index=True, float_format="%f")

if parameters.plotBayesFactorMatrix:
    plot_bf_matrix(celltype_centers_bf_matrix.rename (index=lambda x: "to_" + x), path="{}_BF.pdf".format(output_file))

search_dic_sub = matrix2dic (subcelltype_centers_dis_matrix, celltype2subtype)

search_dic = {}
for k1, v1 in search_dic_sub.items ():
    k1 = k1.split ("*SustechJinLab*")[0]
    # print(k1)
    if k1 not in search_dic.keys ():
        search_dic[k1] = []
    for k2, v2 in v1.items ():
        k2 = k2.split ("*SustechJinLab*")[0]
        if k2 != k1 and k2 not in search_dic[k1]:
            search_dic[k1].append(k2)

k_neighbor_clu = parameters.kNeighbor
for k1, v1 in search_dic.items ():
    if len (v1) > k_neighbor_clu:
        cutoff = sorted(list(map(lambda clu: celltype_centers_bf_matrix[k1][clu], v1)))[-k_neighbor_clu]
        search_dic[k1] = [clu for clu in v1 if celltype_centers_bf_matrix[k1][clu] >= cutoff]

mandatory_link = parameters.mandatoryLink
only_mandatory_link = parameters.onlyMandatoryLink

cancel_link = parameters.cancelLink

print("The original search graph: ")
print_dict(search_dic)

undirected_link = parameters.undirectedLink
if undirected_link == True:
    for kk in search_dic.keys():
        for vv in search_dic[kk]:
            if kk not in search_dic[vv]:
                search_dic[vv].append(kk)
    print ("\nTransform links into undirected links: ")
    print_dict (search_dic)

for k, v in mandatory_link.items():
    if only_mandatory_link:
        search_dic[k] = list(v)
    else:
        search_dic[k] = list(set(search_dic[k]+v))

if len(mandatory_link) != 0:
    print ("\nMandatory links: ")
    print_dict (mandatory_link)
    if only_mandatory_link:
        print("Only mandatory links !!!!")
    print("The search graph: ")
    print_dict (search_dic)

for k, v in cancel_link.items ():
    for p in v:
        if p in search_dic[k]:
            search_dic[k].remove(p)

if len(cancel_link) != 0:
    print ("\nCancel links: ")
    print_dict (cancel_link)
    print("The search graph: ")
    print_dict (search_dic)

start_points = parameters.beginning
print("Beginnings: {}".format(start_points))
end_points = parameters.destination
print("Destinations: {}".format(end_points))
find_more_end_point = parameters.newDestination
new_find_end_points = None

paths = []
if len (end_points) == 0:

    Queue = [[start_point] for start_point in start_points]
    while len (Queue) != 0:
        path = Queue.pop (0)
        last_point = path[-1]
        next_points = [point for point in search_dic[last_point] if point not in path]
        if len (next_points) == 0:
            paths.append (path)
        else:
            for point in next_points:
                newpath = path + [point]
                Queue.append (newpath)
else:

    Queue = [[start_point] for start_point in start_points]
    while len (Queue) != 0:
        path = Queue.pop (0)
        last_point = path[-1]
        if last_point not in end_points:
            next_points = [point for point in search_dic[last_point] if point not in path]
        else:
            next_points = []
        if len (next_points) == 0:
            paths.append (path)
        else:
            for point in next_points:
                newpath = path + [point]
                Queue.append (newpath)

    if find_more_end_point:
        print ("Finding new ends.")
        new_find_end_points = [path[-1] for path in paths if path[-1] not in end_points]
        new_find_end_points = list (set (new_find_end_points))

        if len (new_find_end_points) != 0:
            print ("Find new ends: {}".format (" ".join (new_find_end_points)))
        else:
            print ("No new end is found.")

    else:
        print ("Ignoring unlisted end points: {} paths -->>".format(len(paths)), end=" ")
        paths = [path for path in paths if path[-1] in end_points]
        print("{} paths".format(len(paths)))

if len(paths) == 0:
    print("Sorry, no path was left with the listed end points, "
          "maybe you can add more points in the list "
          "or use --newDestination(-nd)")
    sys.exit (0)

paths = [(get_score_of_a_path(celltype_centers_bf_matrix, path), path) for path in paths]
paths.sort (key=lambda x: len(x[1]))

final_paths = []
final_paths_sets = [set(paths[0][1])]
final_paths_indexs = [{0}]

for i in range(1, len (paths)):
    the_path = paths[i][1]
    the_set = set (the_path)
    flag_found = False
    for j in range (len (final_paths_sets)):
        old_set = final_paths_sets[j]
        if the_set <= old_set:
            final_paths_indexs[j].add(i)
            flag_found = True
        elif the_set > old_set:
            final_paths_sets[j] = the_set
            final_paths_indexs[j].add(i)
            flag_found = True
    if not flag_found:
        final_paths_sets.append (the_set)
        final_paths_indexs.append ({i})

temp_final_paths_sets = [final_paths_sets[0]]
temp_final_paths_indexs = [final_paths_indexs[0]]

for temp1 in range(1, len(final_paths_sets)):
    new_set = final_paths_sets[temp1]
    flag_found = False
    for temp2 in range(len(temp_final_paths_sets)):
        old_set = temp_final_paths_sets[temp2]
        if new_set == old_set:
            temp_final_paths_indexs[temp2] |= final_paths_indexs[temp1]
            flag_found = True
        # elif new_set > old_set:
        #     temp_final_paths_sets[temp2] = new_set
        #     temp_final_paths_indexs[temp2] += final_paths_indexs[temp1]
        #     flag_found = True
    if not flag_found:
        temp_final_paths_sets.append (new_set)
        temp_final_paths_indexs.append (final_paths_indexs[temp1])

print("\nFinal paths: ")

final_paths_indexs = temp_final_paths_indexs
final_paths_sets = temp_final_paths_sets

for indexs_list in final_paths_indexs:
    the_path = [paths[m] for m in indexs_list]
    the_avgbf = [paths[m][0] for m in indexs_list]
    final_paths.append (the_path[the_avgbf.index (max (the_avgbf))])

final_paths.sort(reverse=True)
paths.sort(reverse=True)
for p in final_paths:
    print(p)

write_path_file(path=output_file, infor=final_paths)
write_path_file(path=all_output_file, infor=paths)

