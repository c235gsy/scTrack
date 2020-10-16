import numpy as np
from scipy.interpolate import UnivariateSpline
import sys


def np_put(p):
    """
    a[order][np_put(order)] = a

    """
    n = p.size
    s = np.zeros (n, dtype=np.int32)
    i = np.arange (n, dtype=np.int32)
    return np.put (s, p, i)  # s[p[i]] = i


def soomth(curve, k=3, bbox=[None, None], smoothness=1.0):
    n_samples = len(curve)
    n_features = len(curve[0])
    interpolation_functions = [
        UnivariateSpline (x=list(range(0, n_samples)),
                          y=curve[:, n_f],
                          k=k,
                          bbox=bbox,
                          s=smoothness*n_samples)
        for n_f in range (n_features)]
    curve_seq = np.array ([function (list(range(0, n_samples))) for function in interpolation_functions]).T
    return curve_seq


def project_points(points, curve, init_lambda=None, inter_dimension=3, extend=2, n_curve_seq=int (1e4),
                   return_curve_seq=False, lambda_model="pseudotime", smoothness=1.0):
    """
    points : the points used to project
    curve : the curve points project to
    inter_dimension : the interpolation dimension of the interpolation functions
    extend : the rate of the extension of the interpolation prediction
    n_curve_seq : the number of point used to general the curve sequences
    return_curve_seq : default False, if True, the dictionary returned contains a key "curve_seq"

    return:
        a dictionary:
            "projection" : the projections of the points
            "order" : the order array of the points in the curve
            "lambda_points" : the lambda of the points
            "extend" : the rate of the extension of the interpolation prediction
            "n_curve_seq" : the number of point used to general the curve sequences
            "curve_seq" : the whole curve sequences
    """

    n_points = len (points)  # the number of the points used to be projected
    n_curves = len (curve)  # the number of the points in the curve
    n_features = len (points[0])  # the number of the features for every point
    n_curve_seq_all = int (n_curve_seq * (100 + extend * 2) / 100)
    # the number of the points in the predicted curve sequences

    if init_lambda is None:
        lambda_curve = np.linspace (0, 100, n_curves)
    else:
        lambda_curve = 100 * (init_lambda - init_lambda.min ()) / (init_lambda.max () - init_lambda.min ())
    # calculate the initial lambda of the points used to be projected
    lambda_seq = np.linspace (0 - extend, 100 + extend, n_curve_seq_all)
    # calculate the initial lambda of the points in the predicted curve

    sorted_lambda_curve, lambda_curve_idx = np.unique (lambda_curve, return_index=True)
    interpolation_functions = [
        UnivariateSpline (x=sorted_lambda_curve,
                          y=curve[lambda_curve_idx, n_f],
                          k=inter_dimension,
                          s=smoothness*len(sorted_lambda_curve))
        for n_f in range (n_features)]
    # a list of interpolation functions used to predict the curve sequences
    curve_seq = np.array ([function (lambda_seq) for function in interpolation_functions]).T
    # print (curve_seq[0])
    # points in the predicted curve sequences

    min_dist_idx_ = np.array([np.argmin(np.sum((p - curve_seq)**2, axis=1), axis=0) for p in points])
    projection = curve_seq[min_dist_idx_]  # get the projection of points
    order = np.argsort (lambda_seq[min_dist_idx_])  # the order array of projections of points in the curve

    if lambda_model == "arc":
        ord_projection = projection[order]
        lambda_points = [0]
        arc = 0
        for i in range (len (ord_projection) - 1):
            arc += np.sqrt (np.sum ((ord_projection[i] - ord_projection[i + 1]) ** 2))
            lambda_points.append (arc)
        lambda_points = np.array (lambda_points)[np_put (order)]
        # print(lambda_points)
        # print(lambda_points[np.argsort (lambda_seq[min_dist_idx])])

    elif lambda_model == "pseudotime":
        lambda_points = lambda_seq[min_dist_idx_]
        lambda_points = 100 * (lambda_points - lambda_points.min ()) / (lambda_points.max () - lambda_points.min ())

    else:
        print ("The lambda_model must be chosen from \"arc\" and \"pseudotime\" ")
        sys.exit ()

    # lambda_points = lambda_seq[min_dist_idx]
    # lambda_points = 100 * (lambda_points - lambda_points.min ()) / (lambda_points.max () - lambda_points.min ())
    # get the relevant lambda of points

    output = {"projection": projection, "order": order, "lambda_points": lambda_points, "extend": extend,
              "n_curve_seq": n_curve_seq}

    if return_curve_seq:
        output["curve_seq"] = curve_seq

    return output


def get_princurve_curve(points, start=None, extend_init=20, extend_iter=2, cutoff_iter=0.001, max_iter=10,
                        n_curve_seq=int (1e4), lambda_model="pseudotime", return_curve_seq=False,
                        inter_dimension_init=3, inter_dimension_iter=3, smoothness=1.0):
    """


    """

    n_points = len (points)  # the number of points
    projection = None
    lambda_curve = None
    # order = None

    # initialization
    if start is None:  # use the first principal component line
        u, s, v = np.linalg.svd (points)  # svd decomposition
        expand_v = np.expand_dims (v[0], axis=1)  # preparing for project
        projection = points @ expand_v @ expand_v.T  # the projection in the first principal component line
        lambda_curve = u[:, 0] * s[0]  # initialization of the lambda
        lambda_curve = 100 * (lambda_curve - lambda_curve.min ()) / (lambda_curve.max () - lambda_curve.min ())
    else:  # use the given points
        first_projection = project_points (points=points, curve=start, init_lambda=None,
                                           lambda_model=lambda_model,
                                           inter_dimension=inter_dimension_init,
                                           extend=extend_init,
                                           n_curve_seq=n_curve_seq,
                                           return_curve_seq=return_curve_seq,
                                           smoothness=smoothness)
        projection = first_projection["projection"]
        lambda_curve = first_projection["lambda_points"]

    D2_list = []
    D2 = np.sum ((points - projection) ** 2) / n_points
    print ("initialization D^2/n:", D2)
    D2_list.append (D2)

    iter_projection = None
    for i in range (max_iter):
        iter_projection = project_points (points=points, curve=points, init_lambda=lambda_curve,
                                          lambda_model=lambda_model, inter_dimension=inter_dimension_iter,
                                          extend=extend_iter, n_curve_seq=n_curve_seq,
                                          return_curve_seq=return_curve_seq, smoothness=smoothness)
        projection = iter_projection["projection"]
        lambda_curve = iter_projection["lambda_points"]

        D2 = np.sum ((points - projection) ** 2) / n_points
        print (i+1, " iteration D^2/n:", D2)
        D2_list.append (D2)
        if abs (D2_list[-2] - D2_list[-1]) / D2_list[-2] < cutoff_iter:
            break

    return iter_projection, D2_list