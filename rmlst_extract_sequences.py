#!/usr/bin/env python3

import argparse
import glob
import gzip
import os
import re
import subprocess
import sys
from multiprocessing.pool import ThreadPool


def get_arguments():
    parser = argparse.ArgumentParser(description='Distance matrix from rMLST gene identity')

    parser.add_argument('assembly_dir', type=str,
                        help='Directory containing assembly fasta files (can be gzipped)')
    parser.add_argument('rmlst_dir', type=str,
                        help='Directory containing rMLST fasta files')

    parser.add_argument('--search_tool', type=str, required=False, default='minimap',
                        help='Which tool to use for finding rMLST genes (must be "blast" or "minimap"')
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
    gene_files = get_fastas(args.rmlst_dir)

    for assembly in assembly_files:
        db_name = build_database(assembly, args.search_tool)
        if args.search_tool == 'minimap':
            assembly_seqs = load_fasta(assembly)
        else:
            assembly_seqs = None

        assembly_name = os.path.basename(assembly)

        print('\nSearching for genes in {}'.format(assembly_name))
        gene_seqs = {}

        if args.threads == 1:
            for gene in gene_files:
                find_gene_seq_for_assembly(assembly_name, db_name, gene, gene_seqs, assembly_seqs, args)
        else:
            pool = ThreadPool(args.threads)
            pool.map(lambda gene: find_gene_seq_for_assembly(assembly_name, db_name, gene, gene_seqs,
                                                             assembly_seqs, args), gene_files)
        clean_up_database(db_name, args.search_tool)

        with open(assembly + '.rmlst', 'wt') as rmlst_genes:
            for name in sorted(gene_seqs.keys()):
                if gene_seqs[name] is not None:
                    rmlst_genes.write('>')
                    rmlst_genes.write(name)
                    rmlst_genes.write('\n')
                    rmlst_genes.write(gene_seqs[name])
                    rmlst_genes.write('\n')




def find_gene_seq_for_assembly(assembly, db_name, gene, gene_seqs, assembly_seqs, args):
    query_name = query_name_from_filename(gene)
    hit = get_best_match(db_name, gene, args.search_tool, assembly_seqs, args.min_cov, args.min_id)
    if hit is None:
        print('    {}: none'.format(query_name))
        gene_seqs[query_name] = None
    else:
        print('    {}: {}, {:.2f}% cov, {:.2f}% id, {} bp '
              '({}...{})'.format(query_name, hit.name, hit.coverage,hit.identity, len(hit.seq),
                                 hit.seq[:6], hit.seq[-6:]))
        gene_seqs[query_name] = hit.seq


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


def build_database(assembly, search_tool):
    print()
    if search_tool == 'blast':
        return build_blast_database(assembly)
    elif search_tool == 'minimap':
        return build_minimap_database(assembly)
    else:
        assert False


def build_blast_database(assembly):
    if assembly.endswith('.gz'):
        unzipped_name = assembly[:-3]
        subprocess.run('gunzip -c {} > {}'.format(assembly, unzipped_name), shell=True)
    else:
        unzipped_name = assembly
    subprocess.run('makeblastdb -dbtype nucl -in {}'.format(unzipped_name), shell=True)
    if unzipped_name != assembly:
        os.remove(unzipped_name)
    return unzipped_name


def build_minimap_database(assembly):
    if assembly.endswith('.gz'):
        unzipped_name = assembly[:-3]
    else:
        unzipped_name = assembly
    subprocess.run('minimap2 -d {}.mmi {}'.format(unzipped_name, assembly), shell=True)
    return unzipped_name


def clean_up_database(db_name, search_tool):
    if search_tool == 'blast':
        clean_up_blast_database(db_name)
    elif search_tool == 'minimap':
        clean_up_minimap_database(db_name)
    else:
        assert False


def clean_up_blast_database(db_name):
    os.remove(db_name + '.nhr')
    os.remove(db_name + '.nin')
    os.remove(db_name + '.nsq')


def clean_up_minimap_database(db_name):
    os.remove(db_name + '.mmi')


def get_best_match(db_name, query, search_tool, assembly_seqs, min_cov, min_id):
    if search_tool == 'blast':
        return get_best_match_using_blast(db_name, query, min_cov, min_id)
    elif search_tool == 'minimap':
        return get_best_match_using_minimap(db_name, query, min_cov, min_id, assembly_seqs)
    else:
        assert False


def get_best_match_using_blast(db_name, query, min_cov, min_id):
    blast_out = \
        subprocess.check_output('blastn -db {} -query {} '
                                '-outfmt "6 bitscore pident qcovs sseq qseqid"'.format(db_name, query),
                                shell=True).decode()
    hits = sorted((BlastHit(x) for x in blast_out.split('\n') if x), key=lambda x: x.bitscore)
    if not hits:
        return None
    best = hits[-1]
    if best.coverage < min_cov or best.identity < min_id:
        return None
    return best


def get_best_match_using_minimap(db_name, query, min_cov, min_id, assembly_seqs):
    with open(os.devnull, 'w') as devnull:
        minimap_out = subprocess.check_output('minimap2 -c {}.mmi {}'.format(db_name, query),
                                              shell=True, stderr=devnull).decode()
    hits = sorted((MinimapHit(x) for x in minimap_out.split('\n') if x), key=lambda x: x.score)
    if not hits:
        return None
    best = hits[-1]
    if best.coverage < min_cov or best.identity < min_id:
        return None
    best.seq = get_target_seq(assembly_seqs, best)
    return best


def get_target_seq(assembly_seqs, hit):
    contig_seq = assembly_seqs[hit.contig_name]
    if hit.start < hit.end:
        target_seq = contig_seq[hit.start:hit.end]
        if hit.strand == '+':
            return target_seq
        else:
            return reverse_complement(target_seq)
    else:
        assert False


def query_name_from_filename(query_filename):
    return os.path.basename(os.path.splitext(os.path.basename(query_filename))[0])


class BlastHit(object):
    def __init__(self, blast_line):
        parts = blast_line.split('\t')
        self.bitscore = float(parts[0])
        self.identity = float(parts[1])
        self.coverage = float(parts[2])
        self.seq = parts[3]
        self.name = parts[4]


class MinimapHit(object):
    def __init__(self, minimap_line):
        parts = minimap_line.split('\t')

        self.name = parts[0]

        q_length = int(parts[1])
        q_start = int(parts[2])
        q_end = int(parts[3])
        self.coverage = 100.0 * (q_end - q_start) / q_length

        self.strand = parts[4]
        self.contig_name = parts[5]
        self.start = int(parts[7])
        self.end = int(parts[8])

        matches = int(parts[9])
        alignment_length = int(parts[10])
        self.identity = 100.0 * matches / alignment_length

        self.score = int([x for x in parts if x.startswith('AS:i:')][0][5:])

        self.seq = ''


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
    contig_names = set(fasta_seqs.keys())
    return fasta_seqs


REV_COMP_DICT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                 'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                 'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
                 'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
                 'D': 'H', 'H': 'D', 'N': 'N',
                 'r': 'y', 'y': 'r', 's': 's', 'w': 'w',
                 'k': 'm', 'm': 'k', 'b': 'v', 'v': 'b',
                 'd': 'h', 'h': 'd', 'n': 'n',
                 '.': '.', '-': '-', '?': '?'}


def reverse_complement(seq):
    return ''.join([complement_base(x) for x in seq][::-1])


def complement_base(base):
    try:
        return REV_COMP_DICT[base]
    except KeyError:
        return 'N'


if __name__ == '__main__':
    main()
