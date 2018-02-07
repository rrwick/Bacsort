#!/usr/bin/env python3

import sys
from Bio import Phylo
from pathlib import Path
from collections import Counter


def main():
    cluster_accessions = load_cluster_accessions()
    accession_species = load_accession_species()

    tree = Phylo.read('tree/tree_1.newick', 'newick')
    tree.root_at_midpoint()

    for clade in tree.find_clades():
        if clade.name is None:
            continue
        cluster_name = clade.name[:-7]
        accessions = cluster_accessions[cluster_name]
        species = [accession_species[a] for a in accessions]
        species_counts = Counter(species)
        formatted_species_counts = []
        for species, count in sorted(species_counts.items(), key=lambda x: (1 / x[1], x[0])):
            formatted_species_counts.append(str(count) + ' x ' + species)
        formatted_species_counts = ', '.join(formatted_species_counts)
        clade.name = cluster_name + ' (' + formatted_species_counts + ')'

    Phylo.write(tree, 'tree/tree_2.newick', 'newick')


def load_cluster_accessions():
    cluster_accessions = {}
    with open('clusters/cluster_accessions', 'rt') as accessions_file:
        for line in accessions_file:
            parts = line.split('\t')
            cluster_name = parts[0]
            accessions = [x[:-7] for x in parts[1].split(',')]
            cluster_accessions[cluster_name] = accessions
    return cluster_accessions


def load_accession_species():
    accession_species = {}

    # Load from NCBI metadata first...
    data_files = [str(x) for x in Path.cwd().glob('assemblies/*/data.tsv')]
    for data_file in data_files:
        with open(data_file, 'rt') as data:
            for line in data:
                parts = line.split('\t')
                accession = parts[0]
                if accession == 'assembly_accession':
                    continue
                species = parts[9]
                species_parts = species.split(' ')[0:2]
                if species_parts[1] == 'sp.':
                    species_parts[1] = 'unknown'
                species = ' '.join(species_parts)
                accession_species[accession] = species

    # ...and then load from the user-defined file, so they can overwrite NCBI species.
    if Path('user-defined_accession_species').is_file():
        with open('user-defined_accession_species', 'rt') as user_species:
            for line in user_species:
                parts = line.split('\t')
                if parts[0] == 'Accession':
                    continue
                if len(parts) < 2:
                    continue
                accession, species = parts[0], parts[1]
                accession_species[accession] = species

    # If the user-defined file doesn't exist yet, make an empty one.
    else:
        with open('user-defined_accession_species', 'wt') as user_species:
            user_species.write('Accession\tSpecies\tNotes\n')

    return accession_species

if __name__ == '__main__':
    main()
