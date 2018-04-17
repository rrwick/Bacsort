#!/usr/bin/env python3
"""
Copyright 2018 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Bacsort

This merges Bacsorted assemblies into a Centrifuge library build.

This file is part of Bacsort. Bacsort is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Bacsort is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Bacsort. If
not, see <http://www.gnu.org/licenses/>.
"""

import argparse
import gzip
import os
import pathlib
import sys


def get_arguments():
    parser = argparse.ArgumentParser(description='Cluster assemblies in each genus')

    parser.add_argument('binned_assembly_dir', type=str,
                        help='Directory of Bacsort-binned assemblies')
    parser.add_argument('centrifuge_db_dir', type=str,
                        help='Directory of Bacsort-binned assemblies')

    parser.add_argument('--min_contig_len', type=int, default=10000,
                        help='Contigs shorter than this will not be included in the Centrifuge '
                             'library')

    args = parser.parse_args()
    return args


def main():
    args = get_arguments()

    bin_dir = pathlib.Path(args.binned_assembly_dir)
    genus_dirs = [x for x in bin_dir.iterdir() if x.is_dir()]
    genera = sorted(str(x).split('/')[-1] for x in genus_dirs)

    centrifuge_assemblies = \
        sorted(pathlib.Path(args.centrifuge_db_dir).glob('library/bacteria/*.fna'))

    original_seqid2taxid = args.centrifuge_db_dir + '/seqid2taxid_original.map'
    new_seqid2taxid = args.centrifuge_db_dir + '/seqid2taxid.map'
    os.rename(str(new_seqid2taxid), str(original_seqid2taxid))
    library_dir = pathlib.Path(args.centrifuge_db_dir) / 'library' / 'bacteria'

    ids_to_names, names_to_ids, ids_to_nodes = get_ncbi_taxonomy(args.centrifuge_db_dir)

    print('\nProcessing Bacsorted assemblies')
    print('-------------------------------------------------------------------------------')
    print('This step looks at all assemblies that were binned by Bacsort and does two')
    print('things: 1) makes a seqid2taxid record for each of the assembly\'s contigs and')
    print('2) unzips the assembly into Centrifuge\'s library.')
    print('-------------------------------------------------------------------------------')
    with open(new_seqid2taxid, 'wt') as seqid2taxid:
        for genus in genera:
            if genus == 'Unknown':
                continue
            if genus not in names_to_ids:
                print('WARNING: {} not in NCBI taxonomy, cannot use'.format(genus))
                continue

            genus_dir = bin_dir / pathlib.Path(genus)
            species_dirs = [x for x in genus_dir.iterdir() if x.is_dir()]
            species_names = sorted(str(x).split('/')[-1] for x in species_dirs)
            for species in species_names:
                species_dir = genus_dir / pathlib.Path(species)
                binomial = genus + ' ' + species
                if species == 'unknown':
                    tax_id = names_to_ids[genus]
                    print('{} -> {} (genus level)'.format(binomial, tax_id))
                elif binomial in names_to_ids:
                    tax_id = names_to_ids[binomial]
                    print('{} -> {} (species level)'.format(binomial, tax_id))
                else:
                    tax_id = names_to_ids[genus]
                    print('{} -> {} (genus level)'.format(binomial, tax_id))
                assemblies = sorted(x for x in species_dir.iterdir()
                                    if x.is_file() and str(x).endswith('.fna.gz'))
                for assembly in assemblies:
                    assembly_name = str(assembly).split('/')[-1].replace('.gz', '')
                    new_location = str(library_dir / assembly_name)
                    print('  {} -> {}'.format(assembly, new_location))
                    with open(new_location, 'wt') as new_fasta:
                        for contig_name, seq in load_fasta(str(assembly)):
                            if len(seq) >= args.min_contig_len:
                                seqid2taxid.write(' '.join([contig_name, str(tax_id)]))
                                seqid2taxid.write('\n')
                                new_fasta.write('>')
                                new_fasta.write(contig_name)
                                new_fasta.write('\n')
                                new_fasta.write(seq)
                                new_fasta.write('\n')

    print('\n\n\nAdding existing Centrifuge seqid2taxid entries')
    print('-------------------------------------------------------------------------------')
    print('This step determines which of Centrifuge\'s existing assemblies to exclude')
    print('(those covered by Bacsorted genera) and which to include (those in different)')
    print('genera. Excluded assemblies are renamed and included assemblies are added to')
    print('the seqid2taxid file.')
    print('-------------------------------------------------------------------------------')
    ids_to_remove = set()
    removed_contigs = set()

    # For the purposes of Bacsort, Shigella is part of E. coli.
    if 'Escherichia' in genera and 'Shigella' not in genera:
        genera.append('Shigella')
        genera = sorted(genera)

    for genus in genera:
        if genus not in names_to_ids:
            continue
        tax_id = names_to_ids[genus]
        ids_to_remove.add(tax_id)
        genus_node = ids_to_nodes[tax_id]
        descendant_ids = genus_node.get_descendant_ids(ids_to_nodes)
        for descendant_id in descendant_ids:
            ids_to_remove.add(descendant_id)
        id_str = ','.join(str(x) for x in sorted([tax_id] + descendant_ids))
        print('Excluding {} tax IDs: {}\n'.format(genus, id_str))
    with open(original_seqid2taxid, 'rt') as original:
        with open(new_seqid2taxid, 'at') as filtered:
            for line in original:
                parts = line.strip().split()
                contig = parts[0]
                tax_id = int(parts[1])
                if tax_id in ids_to_remove:
                    removed_contigs.add(contig)
                else:
                    filtered.write(line)
    for assembly in centrifuge_assemblies:
        contig_names = load_contig_names(str(assembly))
        if all(c in removed_contigs for c in contig_names):
            new_name = str(assembly) + '.excluded'
            os.rename(str(assembly), new_name)
            print('{} -> {}'.format(str(assembly), new_name))


class NcbiTaxonomyNode(object):
    def __init__(self, line, ids_to_names):
        parts = [x.strip() for x in line.split('|')]
        self.id = int(parts[0])
        self.name = ids_to_names[self.id]
        self.parent_id = int(parts[1])
        self.child_ids = set()
        self.rank = parts[2]

    def get_descendant_ids(self, ids_to_nodes):
        visited, queue = set(), [self.id]
        while queue:
            v = queue.pop(0)
            if v not in visited:
                visited.add(v)
                queue += sorted(ids_to_nodes[v].child_ids)
        visited.remove(self.id)
        return sorted(visited)


def get_ncbi_taxonomy(centrifuge_db_dir):
    names_filename = centrifuge_db_dir + '/taxonomy/names.dmp'
    nodes_filename = centrifuge_db_dir + '/taxonomy/nodes.dmp'
    print('\nLoading NCBI taxonomy... ', end='', flush=True)
    ids_to_names, names_to_ids = get_ncbi_name_dicts(names_filename)
    ids_to_nodes = {}
    with open(nodes_filename, 'rt') as nodes_file:
        for line in nodes_file:
            node = NcbiTaxonomyNode(line, ids_to_names)
            ids_to_nodes[node.id] = node
    for node in ids_to_nodes.values():
        parent_node = ids_to_nodes[node.parent_id]
        parent_node.child_ids.add(node.id)
    print('done')
    return ids_to_names, names_to_ids, ids_to_nodes


def get_ncbi_name_dicts(names_filename):
    ids_to_names, names_to_ids = {}, {}
    with open(names_filename, 'rt') as names_file:
        for line in names_file:
            parts = [x.strip() for x in line.split('|')]
            if parts[3] != 'scientific name':
                continue
            tax_id = int(parts[0])
            name = parts[1]
            ids_to_names[tax_id] = name
            names_to_ids[name] = tax_id
    return ids_to_names, names_to_ids


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


def load_contig_names(filename):
    contig_names = set()
    open_func = get_open_function(filename)
    with open_func(filename, 'rt') as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                contig_name = line[1:].split()[0]
                if contig_name in contig_names:
                    sys.exit('Error: duplicate contig names in {}'.format(filename))
                contig_names.add(contig_name)
    return sorted(contig_names)


def load_fasta(filename):
    fasta_seqs = []
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
                    fasta_seqs.append((contig_name, sequence))
                    sequence = ''
                name = line[1:]
            else:
                sequence += line
        if name:
            contig_name = name.split()[0]
            fasta_seqs.append((contig_name, sequence))
    return fasta_seqs


if __name__ == '__main__':
    main()
