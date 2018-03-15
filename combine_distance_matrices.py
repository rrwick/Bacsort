#!/usr/bin/env python3

import argparse
import sys
import textwrap
import numpy as np


def get_arguments():
    parser = argparse.ArgumentParser(description='Combine two different distance matrices')

    parser.add_argument('matrix_1', type=str,
                        help='First PHYLIP distance matrix (better at shorter distances)')
    parser.add_argument('matrix_2', type=str,
                        help='Second PHYLIP distance matrix (better at longer distances)')

    parser.add_argument('--regression_min', type=float, required=False, default=0.0,
                        help='Lower end of the regression window')
    parser.add_argument('--regression_max', type=float, required=False, default=0.20,
                        help='Upper end of the regression window')
    parser.add_argument('--blend_min', type=float, required=False, default=0.17,
                        help='Lower end of the blend window')
    parser.add_argument('--blend_max', type=float, required=False, default=0.2,
                        help='Upper end of the blend window')

    args = parser.parse_args()
    return args


def main():
    args = get_arguments()
    print_intro_message(args.matrix_1, args.matrix_2)

    matrix_1_distances, matrix_1_assemblies = load_distance_matrix(args.matrix_1)
    matrix_2_distances, matrix_2_assemblies = load_distance_matrix(args.matrix_2)

    assert matrix_1_assemblies == matrix_2_assemblies
    assemblies = sorted(matrix_1_assemblies)

    # First we do a regression between distances, using an overlapping range where we trust both.
    slope, intercept = distance_regression(matrix_1_distances, matrix_2_distances, assemblies,
                                           args.regression_min, args.regression_max,
                                           args.matrix_1, args.matrix_2)

    # Then we build a new matrix, using matrix 1 for low distances, matrix 2 for large distances
    # (adjusted using our regression values) and a blended region in between.
    combined_matrix = build_combined_matrix(matrix_1_distances, matrix_2_distances, assemblies,
                                            args.blend_min, args.blend_max, slope, intercept)

    print_matrix(combined_matrix, assemblies)


def load_distance_matrix(matrix_filename):
    print('Loading {}'.format(matrix_filename), end=' ', file=sys.stderr, flush=True)
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


def distance_regression(matrix_1_distances, matrix_2_distances, assemblies, regression_min,
                        regression_max, matrix_1_filename, matrix_2_filename):
    print('\nPerforming regression:', file=sys.stderr, flush=True)
    print('gathering distances in the range of {} to {} in '
          '{}'.format(regression_min, regression_max, matrix_1_filename),
          end='', file=sys.stderr, flush=True)
    x, y = [], []
    for i, assembly_1 in enumerate(assemblies):
        for j in range(i, len(assemblies)):
            assembly_2 = assemblies[j]
            m1_distance = matrix_1_distances[assembly_1][assembly_2]
            assert m1_distance == matrix_1_distances[assembly_2][assembly_1]
            if regression_min <= m1_distance < regression_max:
                m2_distance = matrix_2_distances[assembly_1][assembly_2]
                assert m2_distance == matrix_2_distances[assembly_2][assembly_1]
                x.append(m2_distance)
                y.append(m1_distance)
        if i % 100 == 0:
            print('.', end='', file=sys.stderr, flush=True)
    print(' done', file=sys.stderr, flush=True)

    print('linear regression (x = {} distance, y = {} '
          'distance):'.format(matrix_2_filename, matrix_1_filename), file=sys.stderr, flush=True)
    x = np.array(x)
    y = np.array(y)
    a = np.vstack([x, np.ones(len(x))]).T
    slope, intercept = np.linalg.lstsq(a, y, rcond=None)[0]
    print('  slope:     {:.6f}'.format(slope), file=sys.stderr, flush=True)
    print('  intercept: {:.6f}'.format(intercept), file=sys.stderr, flush=True)
    print('  {:.6f} * ({} distance) + {:.6f} = {} distance adjusted to fit '
          '{}'.format(slope, matrix_2_filename, intercept, matrix_2_filename, matrix_1_filename),
          file=sys.stderr, flush=True)
    return slope, intercept


def make_empty_distance_matrix(assemblies):
    matrix = {}
    for i in assemblies:
        matrix[i] = {}
        for j in assemblies:
            matrix[i][j] = None
    return matrix


def build_combined_matrix(matrix_1_distances, matrix_2_distances, assemblies, blend_min, blend_max,
                          slope, intercept):
    combined_matrix = make_empty_distance_matrix(assemblies)
    print('\nBuilding combined matrix', end='', file=sys.stderr, flush=True)
    for i in range(len(assemblies)):
        a1 = assemblies[i]
        for j in range(i, len(assemblies)):
            a2 = assemblies[j]
            m1_distance = matrix_1_distances[a1][a2]
            m2_distance = (matrix_2_distances[a1][a2] * slope) + intercept
            if m1_distance <= blend_min:
                distance = m1_distance
            elif m1_distance >= blend_max:
                distance = m2_distance
            else:  # in between blend_min and blend_max
                distance = blend(m1_distance, m2_distance, blend_min, blend_max)
            combined_matrix[a1, a2] = distance
            combined_matrix[a2, a1] = distance
        if i % 100 == 0:
            print('.', end='', file=sys.stderr, flush=True)
    print(' done', file=sys.stderr, flush=True)
    return combined_matrix


def blend(m1_distance, m2_distance, t2, t3):
    m2_weight = (m1_distance - t2) / (t3 - t2)
    m1_weight = 1.0 - m2_weight
    return (m1_weight * m1_distance) + (m2_weight * m2_distance)


def print_matrix(matrix, assemblies):
    print('Printing matrix to stdout', end='', file=sys.stderr, flush=True)
    print(len(assemblies))
    for i, a1 in enumerate(assemblies):
        print(a1, end='')
        for a2 in assemblies:
            print('\t%.6f' % matrix[(a1, a2)], end='')
        print('')
        if i % 100 == 0:
            print('.', end='', file=sys.stderr, flush=True)
    print(' done\n', file=sys.stderr, flush=True)


def print_intro_message(m1, m2):
    intro = 'This script will create a distance matrix using a combination of distances from ' \
            '{} and {}. Short distances will come from {} and longer distances from {}, with ' \
            'intermediate distances a blend between the two. To ensure a smooth transition ' \
            'between them, a linear regression will be used to adjust {} distances to match ' \
            'those from {}.'.format(m1, m2, m1, m2, m2, m1)
    intro = '\n'.join(textwrap.wrap(intro, 80))
    print('', file=sys.stderr)
    print(intro, file=sys.stderr)
    print('', file=sys.stderr, flush=True)


if __name__ == '__main__':
    main()
