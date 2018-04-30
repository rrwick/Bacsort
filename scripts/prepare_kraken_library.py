#!/usr/bin/env python3
"""
Copyright 2018 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Bacsort

This merges Bacsorted assemblies into a Kraken library build.

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
    parser = argparse.ArgumentParser(description='Add Bacsort assemblies to Kraken database')

    parser.add_argument('binned_assembly_dir', type=str,
                        help='Directory of Bacsort-binned assemblies')
    parser.add_argument('kraken_db_dir', type=str,
                        help='Directory of Kraken database')

    parser.add_argument('--min_contig_len', type=int, default=10000,
                        help='Contigs shorter than this will not be included in the Kraken '
                             'library')

    args = parser.parse_args()
    return args


def main():
    args = get_arguments()

    bin_dir = pathlib.Path(args.binned_assembly_dir)
    genus_dirs = [x for x in bin_dir.iterdir() if x.is_dir()]
    genera = sorted(str(x).split('/')[-1] for x in genus_dirs)

    names_to_ids, ids_to_nodes = get_ncbi_taxonomy(args.kraken_db_dir)

    print('\nFinding assembly taxIDs')
    print('-------------------------------------------------------------------------------')
    print('This step looks at all assemblies that were binned by Bacsort and gets the\n'
          'appropriate taxonomy ID for each assembly.')
    print('-------------------------------------------------------------------------------')

    assemblies_to_taxids = {}

    for genus in genera:
        genus_dir = bin_dir / pathlib.Path(genus)

        if genus == 'Unknown':
            continue
        if genus not in names_to_ids:
            print('WARNING: {} not in NCBI taxonomy, cannot use'.format(genus))
            continue
        genus_id = get_taxid(genus, names_to_ids, ids_to_nodes, 'genus')
        if genus_id is None:
            print('WARNING: {} is ambiguous, cannot use'.format(genus))
            continue

        species_dirs = [x for x in genus_dir.iterdir() if x.is_dir()]
        species_names = sorted(str(x).split('/')[-1] for x in species_dirs)

        for species in species_names:
            species_dir = genus_dir / pathlib.Path(species)

            binomial = genus + ' ' + species
            if species == 'unknown':
                tax_id = genus_id
                print('{} -> {} (genus level)'.format(binomial, tax_id))
            elif binomial in names_to_ids:
                tax_id = get_taxid(binomial, names_to_ids, ids_to_nodes, 'species')
                if tax_id is None:
                    print('WARNING: {} is ambiguous, cannot use'.format(binomial))
                    continue
                print('{} -> {} (species level)'.format(binomial, tax_id))
            else:
                tax_id = genus_id
                print('{} -> {} (genus level)'.format(binomial, tax_id))

            assemblies = sorted(x for x in species_dir.iterdir()
                                if x.is_file() and str(x).endswith('.fna.gz'))
            for assembly in assemblies:
                assemblies_to_taxids[str(assembly)] = tax_id

    print('\nCleaning Kraken library')
    print('-------------------------------------------------------------------------------')
    print('This step goes through the existing Kraken library and removes any contigs')
    print('which are covered by Bacsort\'s genera.')
    print('-------------------------------------------------------------------------------')
    ids_to_remove = set()

    # For the purposes of Bacsort, Shigella is part of E. coli.
    if 'Escherichia' in genera and 'Shigella' not in genera:
        genera.append('Shigella')
        genera = sorted(genera)

    for genus in genera:
        if genus not in names_to_ids:
            continue
        tax_id = get_taxid(genus, names_to_ids, ids_to_nodes, 'genus')
        if tax_id is None:
            print('WARNING: {} is ambiguous, cannot use'.format(genus))
            continue
        ids_to_remove.add(tax_id)
        genus_node = ids_to_nodes[tax_id]
        descendant_ids = genus_node.get_descendant_ids(ids_to_nodes)
        for descendant_id in descendant_ids:
            ids_to_remove.add(descendant_id)
        id_str = ','.join(str(x) for x in sorted([tax_id] + descendant_ids))
        print('Excluding {} tax IDs: {}'.format(genus, id_str))

    library_dir = pathlib.Path(args.kraken_db_dir) / 'library' / 'bacteria'
    original_library = str(library_dir / 'library_original.fna')
    new_library = str(library_dir / 'library.fna')
    os.rename(new_library, original_library)
    with open(original_library, 'rt') as original:
        with open(new_library, 'wt') as new:
            include_contig = True
            for line in original:
                if line.startswith('>'):
                    parts = line.split(' ')[0].split('|')
                    assert parts[0] == '>kraken:taxid'
                    tax_id = int(parts[1])
                    contig_id = parts[2]
                    include_contig = tax_id not in ids_to_remove
                    if not include_contig:
                        print('Excluding contig: {}'.format(contig_id))
                if include_contig:
                    new.write(line)
    print()

    print('\nAdding new assemblies to Kraken library')
    print('-------------------------------------------------------------------------------')
    print('This step prepares Bacsorted assemblies for the Kraken library by adding the')
    print('necessary taxonomic information.')
    print('-------------------------------------------------------------------------------')
    if not pathlib.Path('additional_assemblies').is_dir():
        os.makedirs('additional_assemblies')
    for assembly, tax_id in assemblies_to_taxids.items():
        assembly_name = assembly.split('/')[-1].replace('.gz', '')
        print(assembly_name)
        new_assembly_filename = 'additional_assemblies/' + assembly_name
        if not new_assembly_filename.endswith('.fna'):
            new_assembly_filename += '.fna'
        fasta_seqs = load_fasta(assembly)
        with open(new_assembly_filename, 'wt') as new:
            for contig_name, seq in fasta_seqs:
                if len(seq) > args.min_contig_len:
                    new.write('>')
                    new.write('kraken:taxid|')
                    new.write(str(tax_id))
                    new.write('|')
                    new.write(contig_name)
                    new.write('\n')
                    new.write(seq)
                    new.write('\n')


def get_taxid(name, names_to_ids, ids_to_nodes, rank):
    ids = [x for x in names_to_ids[name] if ids_to_nodes[x].rank == rank]
    if len(ids) == 1:
        return ids[0]

    # If there is more than one match for the name and rank (e.g. Buchnera is a genus in
    # Enterobacterales and in plants), then we look to see if only one is in Bacteria/Archaea.
    prokaryote_ids = []
    for i in ids:
        node = ids_to_nodes[i]
        superkingdom = node.get_superkingdom(ids_to_nodes)
        if superkingdom == 'Bacteria' or superkingdom == 'Archaea':
            prokaryote_ids.append(i)
    if len(prokaryote_ids) == 1:
        return prokaryote_ids[0]
    return None


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

    def get_superkingdom(self, ids_to_nodes):
        if self.rank == 'superkingdom':
            return self.name
        else:
            parent = ids_to_nodes[self.parent_id]
            return parent.get_superkingdom(ids_to_nodes)


def get_ncbi_taxonomy(kraken_db_dir):
    names_filename = kraken_db_dir + '/taxonomy/names.dmp'
    nodes_filename = kraken_db_dir + '/taxonomy/nodes.dmp'
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
    return names_to_ids, ids_to_nodes


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
            if name not in names_to_ids:
                names_to_ids[name] = []
            names_to_ids[name].append(tax_id)
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
