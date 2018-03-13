#!/usr/bin/env python3

import argparse
import edlib
import glob
import gzip
import os
import pathlib
import re
import sys
import multiprocessing as mp


def get_arguments():
    parser = argparse.ArgumentParser(description='FastANI-style identities from rMLST genes')

    parser.add_argument('assembly_dir', type=str,
                        help='Directory containing assembly fasta files and rMLST files')

    parser.add_argument('--threads', type=int, required=False, default=8,
                        help='Number of CPU threads to use')
    parser.add_argument('--min_cov', type=float, required=False, default=95.0,
                        help='Minimum coverage to use in gene search')
    parser.add_argument('--min_id', type=float, required=False, default=90.0,
                        help='Minimum identity to use in gene search')

    args = parser.parse_args()
    return args


def main():
    args = get_arguments()
    assembly_files = get_fastas(args.assembly_dir)
    gene_seqs = {}

    print('\nLoading rMLST genes:', file=sys.stderr)
    for assembly in assembly_files:
        rmlst_file = assembly + '.rmlst'
        assembly_name = os.path.basename(assembly)
        if not pathlib.Path(rmlst_file).is_file():
            sys.exit('Error: {} is missing'.format(rmlst_file))
        print('  {}'.format(rmlst_file), file=sys.stderr)
        gene_seqs[assembly_name] = load_fasta(rmlst_file)

    processes = [mp.Process(target=get_assembly_identity_group,
                            args=(assembly_files, gene_seqs, i, args.threads))
                 for i in range(args.threads)]
    for p in processes:
        p.start()
    for p in processes:
        p.join()


def get_assembly_identity_group(assembly_files, gene_seqs, subset_num, subset_total):
    n = 0
    for i in range(len(assembly_files)):
        assembly_1 = os.path.basename(assembly_files[i])
        for j in range(i, len(assembly_files)):
            assembly_2 = os.path.basename(assembly_files[j])
            if n % subset_total == subset_num:
                id = get_assembly_identity(assembly_1, assembly_2, gene_seqs)
                out_line = '\t'.join([assembly_1, assembly_2, '%.6f' % id]) + '\n'

                # Using sys.stdout.write instead of print helps to avoid some multiprocess
                # printing collisions.
                sys.stdout.write(out_line)
            n += 1


def get_assembly_identity(assembly_1, assembly_2, gene_seqs):
    if assembly_1 == assembly_2:
        return 100.0
    pattern = re.compile(r'\d+[\w=]')
    common_genes = set(gene_seqs[assembly_1]) & set(gene_seqs[assembly_2])
    if not common_genes:
        return 0.0
    alignments = []
    for gene in common_genes:
        gene_1 = gene_seqs[assembly_1][gene]
        gene_2 = gene_seqs[assembly_2][gene]
        result = edlib.align(gene_1, gene_2, 'NW', 'path')
        cigar = [(int(x[:-1]), x[-1]) for x in pattern.findall(result['cigar'])]
        alignment_length = sum(x[0] for x in cigar)
        match_count = sum(x[0] for x in cigar if x[1] == '=')
        identity = match_count / alignment_length
        alignments.append((identity, match_count, alignment_length))
    return 100.0 * sum(x[1] for x in alignments) / sum(x[2] for x in alignments)


def get_fastas(fasta_dir):
    fastas = glob.glob(fasta_dir + '/*.fna.gz')
    fastas += glob.glob(fasta_dir + '/*.fa.gz')
    fastas += glob.glob(fasta_dir + '/*.fas.gz')
    fastas += glob.glob(fasta_dir + '/*.fasta.gz')
    fastas += glob.glob(fasta_dir + '/*.fna')
    fastas += glob.glob(fasta_dir + '/*.fa')
    fastas += glob.glob(fasta_dir + '/*.fas')
    fastas += glob.glob(fasta_dir + '/*.fasta')
    return sorted(fastas)


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
        quit_with_error('cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        quit_with_error('cannot use zip format - use gzip instead')
    return compression_type


def get_open_function(filename):
    if get_compression_type(filename) == 'gz':
        return gzip.open
    else:  # plain text
        return open


def load_fasta(filename):
    fasta_seqs = {}
    open_func = get_open_function(filename)
    with open_func(filename, 'rt') as fasta_file:
        name = ''
        sequence = ''
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    contig_name = name.split()[0]
                    if contig_name in fasta_seqs:
                        sys.exit('Error: duplicate contig names in {}'.format(filename))
                    fasta_seqs[contig_name] = sequence
                    sequence = ''
                name = line[1:]
            else:
                sequence += line
        if name:
            contig_name = name.split()[0]
            if contig_name in fasta_seqs:
                sys.exit('Error: duplicate contig names in {}'.format(filename))
            fasta_seqs[contig_name] = sequence
    return fasta_seqs


if __name__ == '__main__':
    main()
