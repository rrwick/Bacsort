#!/usr/bin/env python3

import argparse
import sys
import numpy as np


def get_arguments():
    parser = argparse.ArgumentParser(description='Extract rMLST sequences to file')

    parser.add_argument('fastani_matrix', type=str,
                        help='PHYLIP distance matrix from FastANI')
    parser.add_argument('rmlst_matrix', type=str,
                        help='PHYLIP distance matrix from rMLST')

    parser.add_argument('--t1', type=float, required=False, default=0.08,
                        help='threshold 1')
    parser.add_argument('--t2', type=float, required=False, default=0.16,
                        help='threshold 2')
    parser.add_argument('--t3', type=float, required=False, default=0.2,
                        help='threshold 3')

    args = parser.parse_args()
    return args


def main():
    args = get_arguments()

    fastani_distances, fastani_assemblies = load_distance_matrix(args.fastani_matrix)
    rmlst_distances, rmlst_assemblies = load_distance_matrix(args.rmlst_matrix)

    assert fastani_assemblies == rmlst_assemblies
    assemblies = sorted(fastani_assemblies)

    # First we do a regression between FastANI distances and rMLST distance, using the relatively
    # close FastANI distances that we trust.
    slope, intercept = fastani_rmlst_regression(fastani_distances, rmlst_distances, assemblies,
                                                args.t1, args.t2)

    # Then we build a new matrix, using FastANI for low distances and rMLST for large distances
    # (adjusted using our regression slope)
    combined_matrix = build_combined_matrix(fastani_distances, rmlst_distances, assemblies,
                                            args.t2, args.t3, slope, intercept)
    print_matrix(combined_matrix, assemblies)


def load_distance_matrix(matrix_filename):
    print('\nLoading {}'.format(matrix_filename), end=' ', file=sys.stderr, flush=True)
    assemblies = []

    try:
        with open(matrix_filename, 'rt') as matrix:
            assembly_count = int(next(matrix).strip())
            for line in matrix:
                assemblies.append(line.split('\t', 1)[0])
        assert len(assemblies) == assembly_count
        print('({} assemblies)'.format(assembly_count), end='', file=sys.stderr, flush=True)

        matrix = {}
        n = 0
        with open(matrix_filename, 'rt') as matrix_file:
            next(matrix_file)  # skip the count line
            for line in matrix_file:
                parts = line.strip().split('\t')
                assert len(parts) == assembly_count + 1
                assembly = parts[0]
                distances = [float(x) for x in parts[1:]]
                assert len(distances) == assembly_count
                matrix[assembly] = dict(zip(assemblies, distances))
                n += 1
                if n % 100 == 0:
                    print('.', end='', file=sys.stderr, flush=True)

        print(' done', file=sys.stderr, flush=True)

    except AssertionError:
        sys.exit('\nError: failed to load {}\n'
                 'Is this a valid PHYLIP distance matrix?'.format(matrix_filename))

    return matrix, assemblies


def fastani_rmlst_regression(fastani_distances, rmlst_distances, assemblies, t1, t2):
    print('\nPerforming regression:', end='', file=sys.stderr, flush=True)
    x, y = [], []

    for i, assembly_1 in enumerate(assemblies):
        for j in range(i, len(assemblies)):
            assembly_2 = assemblies[j]
            fastani_distance = fastani_distances[assembly_1][assembly_2]
            assert fastani_distance == fastani_distances[assembly_2][assembly_1]
            if t1 < fastani_distance < t2:
                rmlst_distance = rmlst_distances[assembly_1][assembly_2]
                assert rmlst_distance == rmlst_distances[assembly_2][assembly_1]
                x.append(rmlst_distance)
                y.append(fastani_distance)

    x = np.array(x)
    y = np.array(y)
    a = np.vstack([x, np.ones(len(x))]).T
    slope, intercept = np.linalg.lstsq(a, y)[0]

    print('slope:    ', slope, file=sys.stderr, flush=True)
    print('intercept:', intercept, file=sys.stderr, flush=True)
    print('adjusted_rMLST_distance = {} * rMLST_distance + {}'.format(slope, intercept),
          file=sys.stderr, flush=True)
    return slope, intercept


def make_empty_distance_matrix(assemblies):
    matrix = {}
    for i in assemblies:
        matrix[i] = {}
        for j in assemblies:
            matrix[i][j] = None
    return matrix


def build_combined_matrix(fastani_distances, rmlst_distances, assemblies, t2, t3, slope,
                          intercept):
    combined_matrix = make_empty_distance_matrix(assemblies)
    print('\nBuilding combined matrix', end='', file=sys.stderr, flush=True)
    for i in range(len(assemblies)):
        a1 = assemblies[i]
        for j in range(i, len(assemblies)):
            a2 = assemblies[j]
            fastani_distance = fastani_distances[a1][a2]
            rmlst_distance = (rmlst_distances[a1][a2] * slope) + intercept
            if fastani_distance <= t2:
                distance = fastani_distance
            elif fastani_distance >= t3:
                distance = rmlst_distance
            else:  # in between t2 and t3
                distance = blend(fastani_distance, rmlst_distance, t2, t3)
            combined_matrix[a1, a2] = distance
            combined_matrix[a2, a1] = distance
        if i % 100 == 0:
            print('.', end='', file=sys.stderr, flush=True)
    print(' done', file=sys.stderr, flush=True)
    return combined_matrix


def blend(fastani_distance, rmlst_distance, t2, t3):
    rmlst_weight = (fastani_distance - t2) / (t3 - t2)
    fastani_weight = 1.0 - rmlst_weight
    return (fastani_weight * fastani_distance) + (rmlst_weight * rmlst_distance)


def print_matrix(matrix, assemblies):
    print('\nPrinting matrix to stdout', end='', file=sys.stderr, flush=True)
    print(len(assemblies))
    for i in assemblies:
        print(i, end='')
        for j in assemblies:
            print('\t%.6f' % matrix[(i, j)], end='')
        print()
        if i % 100 == 0:
            print('.', end='', file=sys.stderr, flush=True)
    print(' done', file=sys.stderr, flush=True)


if __name__ == '__main__':
    main()
