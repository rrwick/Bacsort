#!/usr/bin/env python3
"""
Copyright 2018 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Bacsort

This script compares species definitions before and after taking the species_definitions into
account. It prints assemblies where the species differs. It takes no arguments and must be run in
the Bacsort base directory

This file is part of Bacsort. Bacsort is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Bacsort is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Bacsort. If
not, see <http://www.gnu.org/licenses/>.
"""

import pathlib
import os


def main():
    cluster_accessions = load_all_cluster_accessions()
    accession_species_before = load_accession_species(False)
    accession_species_after = load_accession_species(True)
    cluster_files = sorted(str(x) for x in pathlib.Path.cwd().glob('clusters/*.fna.gz'))
    total_count = 0
    for cluster_file in cluster_files:
        cluster_name = os.path.basename(cluster_file)[:-7]
        accessions = cluster_accessions[cluster_name]
        for accession in accessions:
            before_species = accession_species_before[accession]
            after_species = accession_species_after[accession]
            total_count += 1
            category = get_rename_type(before_species, after_species)
            if category != 'match':
                print('\t'.join([accession, before_species, after_species,
                                 category]))
    print('total:', total_count)


def get_rename_type(before, after):
    before_genus, before_species = before.split(' ')[0:2]
    after_genus, after_species = after.split(' ')[0:2]

    # Shigella counts as E. coli.
    if before_genus == 'Shigella':
        before_genus, before_species = 'Escherichia', 'coli'

    unknown_before_genus = (before_genus.endswith('aceae') or before_genus.endswith('ales'))
    unknown_after_genus = (after_genus == 'Unknown')

    unknown_before_species = (before_species == 'sp.' or before_species == 'bacterium' or
                              before_species == 'symbiont' or before_species == 'endosymbiont')
    unknown_after_species = (after_species == 'unknown')

    # If the previous name was something like 'Enterobacter cloacae complex sp. GN04826', then that
    # counts as an unknown species.
    try:
        if before.split(' ')[2] == 'complex':
            unknown_before_species = True
    except IndexError:
        pass

    if unknown_before_genus:
        assert unknown_before_species
    if unknown_after_genus:
        assert unknown_after_species

    # Check for matches.
    if before == after:
        return 'match'
    if unknown_before_genus and unknown_after_genus:
        return 'match'
    if before_genus == after_genus:
        if unknown_before_species and unknown_after_species:
            return 'match'
        if before_species == after_species and (not unknown_before_species) and \
                (not unknown_after_species):
            return 'match'

    if unknown_before_genus:
        if unknown_after_species:
            return 'given genus, no species'
        else:
            return 'given genus and species'
    if unknown_after_genus:
        return 'removed genus'

    assert (not unknown_before_genus) and (not unknown_after_genus)

    if before_genus == after_genus:
        if unknown_before_species:
            return 'same genus, given species'
        elif unknown_after_species:
            return 'same genus, removed species'
        else:
            return 'same genus, changed species'
    else:
        if unknown_before_species and unknown_after_species:
            return 'different genus, no species'
        elif unknown_before_species:
            return 'different genus, given species'
        elif unknown_after_species:
            return 'different genus, removed species'
        else:
            return 'different genus, changed species'


def load_all_cluster_accessions():
    cluster_accessions = {}
    with open('cluster_accessions', 'rt') as accessions_file:
        for line in accessions_file:
            parts = line.strip().split('\t')
            accessions = [x.split('.fna.gz')[0][:13] for x in parts[1].split(',')]
            cluster_accessions[parts[0]] = accessions
    return cluster_accessions


def load_accession_species(use_species_definitions_file):
    accession_species = {}

    # Load from NCBI metadata first...
    data_files = [str(x) for x in pathlib.Path.cwd().glob('assemblies/*/data.tsv')]
    for data_file in data_files:
        with open(data_file, 'rt') as data:
            for line in data:
                parts = line.split('\t')
                if parts[0] == 'assembly_accession':
                    continue
                accession = parts[0][:13]
                species = parts[9]
                accession_species[accession] = species

    # ...and then load from the user-defined file, so they can overwrite NCBI species.
    if use_species_definitions_file and pathlib.Path('species_definitions').is_file():
        with open('species_definitions', 'rt') as user_species:
            for line in user_species:
                if not line.startswith('GCF'):
                    continue
                parts = line.strip().split('\t')
                if parts[0] == 'Accession':
                    continue
                if len(parts) < 2:
                    continue
                accession, species = parts[0][:13], parts[1]
                accession_species[accession] = species

    return accession_species


if __name__ == '__main__':
    main()
