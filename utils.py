import numpy as np
import pandas as pd
import re
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import seaborn as sns
import scipy
import scipy.special
import pypcurve
import matplotlib as mpl

mpl.use('Agg')

plt.style.use ('ggplot')
plt.rcParams['figure.dpi'] = 400


def sort_key(s):
    if s:
        try:
            c = re.findall ('\d+', s)[0]
        except:
            c = -1
        return int (c)


def get_bayes_factor_matrix(expression, celltype, clusters, num_calculation=100000):
    n = len (clusters)
    distance_between_points = np.zeros (shape=(n, n, num_calculation))
    all_indexes = [np.where (celltype == clu)[0] for clu in clusters]
    all_random_expression = [np.random.choice (clu_indexes, num_calculation) for clu_indexes in all_indexes]
    for i in range (n):
        for j in range (i + 1, n):
            distance_between_points[i][j] = distance_between_points[j][i] = np.sqrt (
                np.sum (np.square (expression[all_random_expression[i]] - expression[all_random_expression[j]]),
                        axis=1))

    # print(distance_between_points)

    bayes_factor_matrix = pd.DataFrame (np.zeros (shape=(n, n)), index=clusters, columns=clusters)
    for ii in range (n):
        clu1 = clusters[ii]
        for jj in range (n):
            clu2 = clusters[jj]
            if ii == jj:
                bayes_factor_matrix[clu1][clu2] = 1
            else:
                clu1_to_clu2_median = np.median (distance_between_points[ii][jj])
                probility_obj = 0.5
                if n > 2:
                    probility_refs = np.sum(distance_between_points[ii][
                                            [kk for kk in range (n) if kk != ii and kk != jj]] <= clu1_to_clu2_median)
                    probility_refs = probility_refs / ((n - 2) * num_calculation)
                else:
                    probility_refs = 0
                bayes_factor_matrix[clu1][clu2] = probility_obj / (probility_obj + probility_refs)

    # bayes_factor_matrix = bayes_factor_matrix.rename (index=lambda x: "to_" + x)
    # print (bayes_factor_matrix)
    return bayes_factor_matrix


def kmeans(x, n_clusters, init='k-means++', algorithm='auto'):
    """
    :return: clustering

    clustering.cluster_centers_ : ndarray of shape (n_clusters, n_features)
    Coordinates of cluster centers. If the algorithm stops before fully converging (see tol and max_iter),
    these will not be consistent with labels_.

    clustering.labels_ : ndarray of shape (n_samples,)
    Labels of each point

    clustering.inertia_ : float
    Sum of squared distances of samples to their closest cluster center.

    clustering.n_iter_ : int
    Number of iterations run.
    """

    clustering = KMeans (n_clusters, init=init, algorithm=algorithm).fit (x)
    return clustering


def get_distance_matricx(x, names):
    """
    :param x:
    :param names:
    :return:
    """
    n_cells = len (x)
    G = np.dot (x, x.T)
    H = np.tile (np.diag (G), (n_cells, 1))
    D = H + H.T - G * 2
    D = np.sqrt (D)
    return pd.DataFrame (D, columns=names, index=names)


def matrix2dic(matrix, celltype2subtype):
    """

    :param matrix:
    :param celltype2subtype:
    :return:
    """
    # print (matrix)
    dictionary = {}
    clusters = matrix.columns
    # print (clusters)
    for clu in clusters:
        dictionary[clu] = {}
        result = matrix[clu][clusters[clusters != clu]].sort_values (ascending=True)
        # print(result)
        relative = celltype2subtype[clu.split ("*SustechJinLab*")[0]]
        # print(relative)
        max_index_1 = np.max (np.where ([sub in relative for sub in result.index]))
        # print(sorted(np.where ([sub not in relative for sub in result.index])[0]))
        max_index_2 = np.min (np.where ([sub not in relative for sub in result.index]))
        # max_index = min (max_nexts - 1, max_index_1)
        max_index = max(max_index_1, max_index_2)
        # max_index = max_index_1
        # print (max_index, len (result))
        cut_off = result[max_index]
        # print (cut_off)
        result = result[result <= cut_off]
        #print (result)
        for i in range (len (result)):
            dictionary[clu][result.index[i]] = result[i]

    return dictionary


def get_score_of_a_path(score_dict, a_path):
    score = 0
    for j in range(len(a_path)-1):
        # score += np.log(1+score_dict[a_path[j]][a_path[j+1]])
        score += np.log (score_dict[a_path[j]][a_path[j + 1]])
        # score += score_dict[a_path[j]][a_path[j + 1]]
    return score/(len(a_path)-1)


def print_dict(a_dict):
    for key, value in a_dict.items():
        print("{}: {}".format(key, value))


def plot_bf_matrix(bf_matirx, title="", path="./", cmap="YlGnBu", fmt='.6f'):
    ax = sns.heatmap (bf_matirx, cmap=cmap, annot=True, fmt=fmt, annot_kws={'size': 4})
    ax.xaxis.tick_top ()
    ax.set_yticklabels (ax.get_yticklabels (), rotation=0)
    plt.xticks (fontsize=10)
    plt.yticks (fontsize=10)
    if title:
        ax.set_title (title, fontsize=10, position=(0.5, 1.10))
    plt.savefig (path, bbox_inches='tight')
    plt.close ()


def write_path_file(path, infor):
    final_path_infor = {"Path{}".format (i): infor[i] for i in range (len (infor))}
    with open (path, 'w') as outFile:
        for p, c in final_path_infor.items ():
            outFile.write ("{},{},{}\n".format (p, c[0], ",".join (c[1])))


def sigmoid(x):
    output = 1 - scipy.special.expit (5 * x)
    return (output - output.min ()) / (output.max () - output.min ())


def get_shared_cells(paths):
    sub_dic = {}
    for path, clus in paths.items ():
        for i in range (1, len (clus)):
            try:
                sub_dic[str (clus[:i])].add (path)
            except:
                sub_dic[str (clus[:i])] = {path}

    sub_paths = list (sub_dic.keys ())
    for sub_path in sub_paths:
        if len (sub_dic[sub_path]) == 1:
            del sub_dic[sub_path]
    sub_paths = list (sub_dic.keys ())

    crossed_paths = []
    for s in sub_dic.values ():
        if s not in crossed_paths:
            crossed_paths.append (s)
    shared_clusters = [[] for i in range (len (crossed_paths))]

    for co in range (len (crossed_paths)):
        for sub_path in sub_paths:
            if sub_dic[sub_path] == crossed_paths[co]:
                if len (eval (sub_path[1:-1])) > len (shared_clusters[co]):
                    shared_clusters[co] = eval (sub_path)

    output = [[list (crossed_paths[i]), shared_clusters[i]] for i in range (len (crossed_paths))]
    output = sorted (output, key=lambda pair: len (pair[1]), reverse=True)
    return output


def get_tmin_tmax(avg_curve, curve):
    distance = np.sqrt (np.sum (np.square (avg_curve - curve), axis=1))

    Q1, Q3 = np.percentile (distance, [25, 75])

    small = np.where (distance <= Q1)[0]
    big = np.where ((distance > Q1) & (distance < Q3))[0]
    tmin = max(small[-1], int(len(distance)*0.2))
    tmax = min(big[-1], int(len(distance)*0.8))
    if tmin >= tmax:
        tmax = len (distance) - 1
    return tmin, tmax


def combine_paths(the_lab, the_paths, cell_type, smoothness):
    shared_clusters = get_shared_cells (the_paths)
    for j in range (len (shared_clusters)):
        these_paths = shared_clusters[j][0]
        share_those = shared_clusters[j][1]
        print (these_paths, "share : ", share_those)

        sub_curves_indexs = []

        for pa in these_paths:
            share = np.where ([ident in share_those for ident in cell_type[the_lab[pa]["cells"]]])[0]
            max_share = share.max()

            unshare = np.where ([ident not in share_those for ident in cell_type[the_lab[pa]["cells"]]])[0]
            min_unshare = unshare.min ()

            if max_share > min_unshare:
                a = np.median(share[(max_share >= share) & (share >= min_unshare)])
                b = np.median(unshare[(max_share >= unshare) & (unshare >= min_unshare)])
                c = int((a+b)/2)
                sub_curves_indexs.append(c)
            else:
                sub_curves_indexs.append(int(max_share))

        sub_curves_indexs_max = max(sub_curves_indexs)
        sub_curves_indexs_max = min([sub_curves_indexs_max]+[len(the_lab[pat]["curve"])-1 for pat in these_paths])
        # print([pat for pat in these_paths])
        sub_curves = np.array ([the_lab[pat]["curve"][:sub_curves_indexs_max + 1] for pat in these_paths])

        # print(sub_curves.shape)
        avg_curve = np.mean (sub_curves, axis=0)
        avg_curve = pypcurve.soomth(curve=avg_curve, smoothness=smoothness)

        for m in range (len (sub_curves)):
            sub_curves_indexs[m] = sub_curves_indexs_max
            # sub_curves_indexs[m] = min(sub_curves_indexs_max, sub_curves_indexs[m])
            sub_curve = sub_curves[m][:sub_curves_indexs[m]+1]
            # sub_point = cell_data[sub_points[m]]
            tmin, tmax = get_tmin_tmax (avg_curve[:sub_curves_indexs[m]+1], sub_curve)
            # print (0, tmin, tmax, len (sub_curve) - 1)

            wm = None
            if tmin < tmax:
                wm = np.array (
                    [1.0] * tmin + list (sigmoid ((np.arange (tmin, tmax + 1) - tmin) / (tmax - tmin) - 0.5)) + [0.0] * (
                            len (sub_curve) - tmax - 1)).reshape (-1, 1)
            elif tmin == tmax:
                wm = np.array ([1.0] * tmin + [0.0] * (len (sub_curve) - tmax)).reshape (-1, 1)

            new_curve = sub_curve * (1 - wm) + avg_curve[:sub_curves_indexs[m]+1] * wm

            the_lab[these_paths[m]]["curve"][:sub_curves_indexs[m]+1] = new_curve
    return the_lab


def curve_2_pseudotime(curve):
    pseudotime = [0]
    a = 0
    for i in range(len(curve)-1):
        b = np.sqrt(np.sum(np.square(curve[i+1] - curve[i])))
        a = a + b
        pseudotime.append(a)
    pseudotime = np.array(pseudotime)
    return pseudotime


def get_color_list(n, hue_start=1/24, hue_end=1, saturation=1, luminance=0.65):
    import colour
    if n <= 0:
        print("n should be bigger than 0")
        return
    elif 0 <= hue_end <= 1 and 0 <= hue_start <= 1:
        out = []
        diss = (hue_end - hue_start) / n
        for i in range(n):
            out.append(colour.Color(hsl=(hue_start + i*diss, saturation, luminance)).hex_l)
        return out



# End
