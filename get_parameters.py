import argparse
import os
import sys

"""parsing and configuration"""


def parse_args_get_paths():
    desc = "The parameters for scTrack !!!!"
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument('-e', '--expMatrix', type=str, help='File path of expression matrix.')
    parser.add_argument('-c', '--clusterArray', type=str,
                        help='File path of cell cluster marker array, '
                             'the number of clusters should not be less than 3.')
    parser.add_argument('-o', '--output', type=str, default='./YOUR_OUTPUT_NAME',
                        help='Paths and names of the output files.')

    parser.add_argument('-np', '--nComponentsPCA', type=int, default=False,
                        help='Integer, the data dimension retained by PCA dimensionality reduction. '
                             'Note that this parameter should be greater than [cell_cluster_number]-1. '
                             'By default it is equal to '
                             'max(5*([cell_cluster_number]-1), int(0.01*[gene_number])).')
    parser.add_argument('-nl', '--nComponentsLDA', type=int, default=False,
                        help='Integer, the data dimension retained by LDA dimensionality reduction. '
                             'Note that this parameter should not be greater than [cell_cluster_number]-1. '
                             'By default it is equal to min([cell_cluster_number]-1, --nComponentsPCA).')

    parser.add_argument('-s', '--sizeSubcluster', type=int, default=50,
                        help='The average size of the subclusters, default=50.')

    parser.add_argument('-k', '--kNeighbor', type=int, default=3,
                        help='Integer, a cell cluster will establish connections with the surrounding k clusters, '
                             'default: 3.')
    parser.add_argument('-nbf', '--nCalculationBayesFactor', type=int, default=100000,
                        help='Integer, the number of cell pairs in calculation to get the distribution of '
                             'distances between cells from 2 cell clusters, '
                             'default: 10e5.')

    parser.add_argument('-b', '--beginning', type=str,
                        help='The format is [clusterA,clusterB,...]. '
                             'Select one or more clusters as the beginning.')
    parser.add_argument('-d', '--destination', type=str, default="[]",
                        help='The format is [clusterA,clusterB,...]. '
                             'Select any number of clusters as the destination. '
                             '(type "[]" if there is no destination). default="[]"')

    parser.add_argument('-nd', '--newDestination', action="store_true", default=False,
                        help='Only work while --destination != []. '
                             'If --destination == [], it will be used automatically'
                             'If it is used, the program will output the paths '
                             'whose end points are not contained in --destination. '
                             'Otherwise, those paths will be ignored.')

    parser.add_argument('-ml', '--mandatoryLink', type=str, default="{}",
                        help='The format is {clusterA:[clusterB,clusterC,...], '
                             '               clusterL:[clusterM,clusterN,...],...}. Default value is {} '
                             'The command in the example means to force the establishment '
                             'of links from clusterA to clusterB and clusterC, '
                             'and links from clusterL to clusterM and clusterN.')

    parser.add_argument('-oml', '--onlyMandatoryLink', action="store_true", default=False,
                        help='If used, the links form clusters included in [--mandatoryLink] to those '
                             ' not included in [--mandatoryLink] will be canceled.')

    parser.add_argument('-cl', '--cancelLink', type=str, default="{}",
                        help='The format is {clusterA:[clusterB,clusterC,...], '
                             '               clusterL:[clusterM,clusterN,...],...}. Default value is {} '
                             'The command in the example means to force the cancel '
                             'of links from clusterA to clusterB and clusterC, '
                             'and links from clusterL to clusterM and clusterN.')

    parser.add_argument('-ul', '--undirectedLink', action="store_true", default=False,
                        help='If used, the links between clusters will be transformed into bidirectional links.')

    parser.add_argument('-pbf', '--plotBayesFactorMatrix', action="store_true", default=False,
                        help='If used, plot the heat map of BF values.')

    return check_args_get_paths(parser.parse_args())


def check_args_get_paths(args):

    # check --expMatrix(-e)
    if not os.path.isfile(args.expMatrix):
        print("File path of expression matrix (--expMatrix(-e)) is wrong, no such file.")
        sys.exit (0)

    # check --clusterArray(-c)
    if not os.path.isfile(args.clusterArray):
        print("File path of cell cluster marker array (--clusterArray(-c)) is wrong, no such file.")
        sys.exit (0)

    # check --nComponentsPCA(-np) & --nComponentsLDA(-nl)
    if type(args.nComponentsPCA) == int:
        if args.nComponentsPCA <= 0:
            print ("--nComponentsPCA(-np) should not be negative.")
            sys.exit (0)

    if type(args.nComponentsLDA) == int:
        if args.nComponentsLDA <= 0:
            print ("--nComponentsLDA(-lp) should not be negative.")
            sys.exit (0)

    if type(args.nComponentsPCA) == int and type(args.nComponentsLDA) == int:
        if args.nComponentsPCA < args.nComponentsLDA:
            print ("--nComponentsPCA(-np) should be greater than --nComponentsLDA(-nl).")
            sys.exit (0)

    # check --sizeSubcluster(-s)
    if args.sizeSubcluster <= 0:
        print ("--sizeSubcluster(-s) should be positive integer.")
        sys.exit (0)

    # check --kNeighbor(-k)
    if args.kNeighbor <= 0:
        print ("--kNeighbor(-k) should be positive integer.")
        sys.exit (0)

    # check --nCalculationBayesFactor(-nbf)
    if args.nCalculationBayesFactor <= 0:
        print ("--nCalculationBayesFactor(-nbf) should be positive integer.")
        sys.exit (0)

    # check --beginning(-b)
    args.beginning = [str(e) for e in eval(str(args.beginning))]

    # check --destination(-d)
    args.destination = [str (e) for e in eval (str (args.destination))]

    # check --newDestination(-nd)
    if args.destination == []:
        args.newDestination = True

    # check --mandatoryLink(-ml)
    args.mandatoryLink = {str(k): [str(v0) for v0 in v] for k, v in eval (str (args.mandatoryLink)).items()}

    # check --cancelLink(-cl)
    args.cancelLink = {str (k): [str (v0) for v0 in v] for k, v in eval (str (args.cancelLink)).items ()}

    return args


def parse_args_ananlyze_paths():
    desc2 = "The parameters for scTrack !!!"
    parser = argparse.ArgumentParser(description=desc2)

    parser.add_argument('-e', '--expMatrix', type=str, default=False, help='File path of expression matrix')
    parser.add_argument('-c', '--clusterArray', type=str, default=False,
                        help='File path of cell cluster marker array, the number of clusters should not be less than 3')
    parser.add_argument('-i', '--idArray', type=str, default=False,
                        help='File path of cell-id array. If missing, cell-id defaults to Cell1, Cell2, Cell3, ...')
    parser.add_argument('-o', '--outputName', type=str, default='YOUR_OUTPUT_NAME',
                        help='Name of the output, default=\'YOUR_OUTPUT_NAME\'.')

    parser.add_argument('-np', '--nComponentsPCA', type=int, default=False,
                        help='Integer, the data dimension retained by PCA dimensionality reduction. '
                             'Note that this parameter should be greater than [cell_cluster_number]-1. '
                             'By default it is equal to '
                             'max(5*([cell_cluster_number]-1), int(0.01*[gene_number])).')
    parser.add_argument('-nl', '--nComponentsLDA', type=int, default=False,
                        help='Integer, the data dimension retained by LDA dimensionality reduction. '
                             'Note that this parameter should not be greater than [cell_cluster_number]-1. '
                             'By default it is equal to [cell_cluster_number]-1.')

    parser.add_argument('-mi', '--maxIter', type=int, default=10,
                        help='Integer, '
                             'the maximum number of iterations in the principal curve fitting process, default=10.')

    parser.add_argument('-ci', '--cutoffIter', type=float, default=0.001,
                        help='Float, the iterations in the principal curve fitting process will stop if '
                             'the rate of change between the point and its corresponding point on the curve '
                             'is less than this value, default=0.001.')

    parser.add_argument('-p', '--pathFile', type=str,
                        help='The path of the path-file inferred in the previous step.')

    parser.add_argument('-cb', '--combineBeginnings', action="store_true",
                        help='If used, combining the shared parts of paths in the beginning.')

    parser.add_argument('-cd', '--combineDestinations', action="store_true",
                        help='If used, combining the shared parts of paths in the destination.')

    parser.add_argument('-p2D', '--plot2D', action="store_true",
                        help='If used, plotting 2D data figures.')

    parser.add_argument('-p3D', '--plot3D', action="store_true",
                        help='If used, plotting a lot 3D data figures.')

    parser.add_argument('-cm', '--colorMap', type=str, default="R", choices=["R", "r", "Python", "py"],
                        help='If R or r, you also need to install Python package rpy2 in order to'
                             'use the function hue_pal() in R to get colormap. '
                             'But the hue_pal() do not use standard hue space. '
                             'If Python or py, program will use Python package colour, '
                             'where the standard hue space is used.')

    parser.add_argument('-s', '--smoothness', type=float, default=1.0,
                        help='Positive smoothing factor used to choose the number of knots. '
                             's = smoothness*[point_number_in_curve]'
                             'More information: https:'
                             '//docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.UnivariateSpline.html')

    return check_args_analyze_paths(parser.parse_args())


def check_args_analyze_paths(args):

    # check --expMatrix(-e)
    if not os.path.isfile(args.expMatrix):
        print("File path of expression matrix (--expMatrix(-e)) is wrong, no such file.")
        sys.exit (0)

    # check --clusterArray(-c)
    if not os.path.isfile(args.clusterArray):
        print("File path of cell cluster marker array (--clusterArray(-c)) is wrong, no such file.")
        sys.exit (0)

    # check --idArray(-c)
    if args.idArray != False:
        if not os.path.isfile (args.idArray):
            print ("File path of cell-id array (--idArray(-i)) is wrong, no such file.")
            sys.exit (0)

    # check --nComponentsPCA(-np) & --nComponentsLDA(-nl)
    if type(args.nComponentsPCA) == int:
        if args.nComponentsPCA <= 0:
            print ("--nComponentsPCA(-np) should not be negative.")
            sys.exit (0)

    if type(args.nComponentsLDA) == int:
        if args.nComponentsLDA <= 0:
            print ("--nComponentsLDA(-lp) should not be negative.")
            sys.exit (0)

    if type(args.nComponentsPCA) == int and type(args.nComponentsLDA) == int:
        if args.nComponentsPCA < args.nComponentsLDA:
            print ("--nComponentsPCA(-np) should be greater than --nComponentsLDA(-nl).")
            sys.exit (0)

    # check --maxIter(-mi)
    if args.maxIter <= 1:
        print ("--maxIter(-mi) should be larger than 1.")
        sys.exit (0)

    # check --cutoffIter(-ci)
    if args.cutoffIter <= 0:
        print ("--cutoffIter(-ci) should not be positive.")
        sys.exit (0)

    # check --pathFile(-p)
    if not os.path.isfile (args.pathFile):
        print ("File path of the path of the path-file inferred in the previous step (--pathFile(-p)) is wrong, "
               "no such file.")
        sys.exit (0)

    # check --smoothness(-s)
    if args.smoothness is not None and args.smoothness <= 0:
        print ("--smoothness(-s) should be positive.")
        sys.exit (0)

    return args


def main():
    print ("Run get_parameters.py")
    paaa = parse_args_get_paths ()
    print (paaa)


if __name__ == '__main__':
    main()
    #parse_args_ananlyze_paths()