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
import sys
import gzip
import os


def get_arguments():
    parser = argparse.ArgumentParser(description='Classify assembly/reads using Mash',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('mash_sketch', type=str,
                        help='Mash sketch file of Bacsorted assemblies/clusters directory')
    parser.add_argument('input', type=str,  nargs='+',
                        help='Assembly FASTA file or read FASTQ files')

    parser.add_argument('--input_type', type=str, required=False, default='auto',
                        choices=['assembly', 'reads', 'auto'],
                        help='Type of input file(s)')
    parser.add_argument('--sketch_size', type=int, required=False, default=100000,
                        help='Mash sketch size')
    parser.add_argument('--threshold', type=float, required=False, default=5.0,
                        help='Mash distances at or below this threshold count as a match '
                             '(expressed as a percent)')
    parser.add_argument('--threads', type=int, required=False, default=4,
                        help='Number of threads to use with Mash')
    parser.add_argument('-m', type=int, required=False, default=3,
                        help='-m option for Mash (only applies when using reads as input)')
    args = parser.parse_args()

    if args.input_type == 'auto':
        input_types = set(get_sequence_filetype(x) for x in args.input)
        if len(input_types) != 1:
            sys.exit('Error: could not determine input type - use --input_type to specify')
        input_type = list(input_types)[0]
        if len(args.input) == 1 and input_type == 'FASTA':
            args.input_type = 'assembly'
        elif input_type == 'FASTQ':
            args.input_type = 'reads'
        else:
            sys.exit('Error: could not determine input type - use --input_type to specify')
    return args


def main():
    args = get_arguments()

    cmd = 'cat ' + ' '.join(args.input) + ' | mash dist '
    cmd += '-s {} '.format(args.sketch_size)
    if args.input_type == 'reads':
        cmd += '-m {} '.format(args.m)
    cmd += '-p {} '.format(args.threads)
    cmd += '{} -'.format(args.mash_sketch)

    with open(os.devnull, 'w') as devnull:
        mash_out = subprocess.check_output(cmd, shell=True, stderr=devnull).decode()

    best_species, best_distance = 'none', 1.0
    for line in mash_out.splitlines():
        parts = line.split('\t')
        distance = float(parts[2])

        genus = parts[0].split('/')[0]
        species = parts[0].split('/')[1]
        binomial = genus + ' ' + species

        if distance <= (args.threshold / 100.0) and distance < best_distance:
            best_species, best_distance = binomial, distance

    if best_species == 'none':
        identity = ''
    else:
        identity = '%.2f' % (100 * (1.0 - best_distance)) + '%'
    sample_name = get_sample_name(args)
    print('\t'.join([sample_name, best_species, identity]))


def get_sample_name(args):
    input = args.input[0]
    sample_name = pathlib.Path(input).name
    if sample_name.endswith('.gz'):
        sample_name = sample_name[:-3]
    if sample_name.endswith('.fna'):
        sample_name = sample_name[:-4]
    if sample_name.endswith('.fasta'):
        sample_name = sample_name[:-6]
    if sample_name.endswith('.fastq'):
        sample_name = sample_name[:-6]
    if sample_name.endswith('.fq'):
        sample_name = sample_name[:-3]
    if args.input_type == 'reads' and sample_name.endswith('_1'):
        sample_name = sample_name[:-2]
    if args.input_type == 'reads' and sample_name.endswith('_R1'):
        sample_name = sample_name[:-3]
    return sample_name


def get_compression_type(filename):
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)
    with open(filename, 'rb') as unknown_file:
        file_start = unknown_file.read(max_len)
    compression_type = 'plain'
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == 'bz2':
        sys.exit('Error: cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        sys.exit('Error: cannot use zip format - use gzip instead')
    return compression_type


def get_open_function(filename):
    if get_compression_type(filename) == 'gz':
        return gzip.open
    else:  # plain text
        return open


def get_sequence_filetype(filename):
    with get_open_function(filename)(filename, 'rt') as fh:
        for line in fh:
            if line.startswith('>'):
                return 'FASTA'
            elif line.startswith('@'):
                return 'FASTQ'
            else:
                sys.exit('Could not determine type of {} (must be FASTA or FASTQ'.format(filename))


if __name__ == '__main__':
    main()
