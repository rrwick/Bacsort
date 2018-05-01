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
    for cluster_file in cluster_files:
        cluster_name = os.path.basename(cluster_file)[:-7]
        accessions = cluster_accessions[cluster_name]
        for accession in accessions:
            before_species = accession_species_before[accession]
            after_species = accession_species_after[accession]
            if before_species != after_species:
                print('\t'.join([accession, before_species, after_species]))


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
                species_parts = species.split(' ')[0:2]

                # Some 'species names' aren't really species names.
                if species_parts[1] == 'sp.':
                    species_parts[1] = 'unknown'
                elif species_parts[1] == 'bacterium':
                    species_parts[1] = 'unknown'
                elif species_parts[1] == 'symbiont':
                    species_parts[1] = 'unknown'
                elif species_parts[1] == 'endosymbiont':
                    species_parts[1] = 'unknown'

                species = ' '.join(species_parts)
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

    # Shigella and E. coli are grouped together.
    for accession, species in accession_species.items():
        if species == 'Escherichia coli' or species.startswith('Shigella '):
            accession_species[accession] = 'Escherichia coli / Shigella'

    return accession_species


if __name__ == '__main__':
    main()
