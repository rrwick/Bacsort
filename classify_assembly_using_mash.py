#!/usr/bin/env python3
"""
Copyright 2018 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Bacsort

This script uses Mash to classify an assembly using Bacsorted assemblies as a reference.

Before it is used, you must create a Mash sketch of assemblies organised by either the
copy_assemblies.py or copy_clusters.py script:
    mash sketch -o sketches -s 100000 */*/*.fna.gz

Then you can use this script like so:
    classify_assembly_using_mash.py sketches.msh query.fasta

This file is part of Bacsort. Bacsort is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Bacsort is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Bacsort. If
not, see <http://www.gnu.org/licenses/>.
"""

import argparse
import subprocess
import pathlib


def get_arguments():
    parser = argparse.ArgumentParser(description='Classify assembly using Mash')

    parser.add_argument('mash_sketch', type=str,
                        help='Mash sketch file of Bacsorted assemblies/clusters directory')
    parser.add_argument('assembly', type=str,
                        help='Assembly FASTA file')

    parser.add_argument('--sketch_size', type=int, required=False, default=100000,
                        help='Mash sketch size')
    parser.add_argument('--threshold', type=float, required=False, default=0.05,
                        help='Mash distances at or below this threshold count as a match')
    parser.add_argument('--contamination_threshold', type=float, required=False, default=0.02,
                        help='If two different genera have Mash distances with a difference of '
                             'less than this value, the assembly is considered contaminated')
    args = parser.parse_args()
    return args


def main():
    args = get_arguments()
    mash_out = subprocess.check_output('mash dist -s 100000 -p 4 {} {}'.format(args.mash_sketch,
                                                                               args.assembly),
                                       shell=True).decode()

    assembly_name = pathlib.Path(args.assembly).name
    if assembly_name.endswith('.gz'):
        assembly_name = assembly_name[:-3]
    if assembly_name.endswith('.fna'):
        assembly_name = assembly_name[:-4]
    if assembly_name.endswith('.fasta'):
        assembly_name = assembly_name[:-6]

    best_species, best_distance = 'none', 1.0
    best_distance_per_genus = {}
    for line in mash_out.splitlines():
        parts = line.split('\t')
        distance = float(parts[2])

        genus = parts[0].split('/')[0]
        species = parts[0].split('/')[1]
        binomial = genus + ' ' + species

        if distance <= args.threshold and distance < best_distance:
            best_species, best_distance = binomial, distance

        if genus not in best_distance_per_genus:
            best_distance_per_genus[genus] = distance
        else:
            best_distance_per_genus[genus] = min(best_distance_per_genus[genus], distance)

    best_distance_per_genus = sorted(best_distance_per_genus.values())
    if len(best_distance_per_genus) > 1:
        genus_dist_1, genus_dist_2 = best_distance_per_genus[:2]
        if abs(genus_dist_1 - genus_dist_2) < args.contamination_threshold:
            best_species += ' (contaminated)'

    identity = '%.2f' % (100 * (1.0 - best_distance)) + '%'
    print('\t'.join([assembly_name, best_species, identity]))


if __name__ == '__main__':
    main()
