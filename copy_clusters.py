#!/usr/bin/env python3
"""
Copyright 2018 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Bacsort

This script is the fifth step of Bacsort. It uses the curated tree to copy assemblies into
species-specific directories.

This file is part of Bacsort. Bacsort is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Bacsort is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Bacsort. If
not, see <http://www.gnu.org/licenses/>.
"""

import os
import pathlib
import shutil


def main():
    cluster_accessions = load_cluster_accessions()
    accession_species = load_accession_species()

    if not pathlib.Path('clusters_binned').is_dir():
        os.makedirs('clusters_binned')

    cluster_files = sorted(str(x) for x in pathlib.Path.cwd().glob('clusters/*.fna.gz'))

    for cluster_file in cluster_files:
        cluster_name = os.path.basename(cluster_file)[:-7]
        accession = cluster_accessions[cluster_name]
        species = accession_species[accession]

        # For the purposes of this script and its directories, Shigella is in E. coli.
        if species == 'Escherichia coli / Shigella':
            species = 'Escherichia coli'

        genus, species = species.split(' ')[0:2]
        print('{} ({}) -> {} {}'.format(accession, cluster_name, genus, species))

        species_dir = 'clusters_binned/' + genus + '/' + species
        if not pathlib.Path(species_dir).is_dir():
            os.makedirs(species_dir)

        dest_file = species_dir + '/' + accession + '.fna.gz'
        shutil.copyfile(cluster_file, dest_file)


def load_cluster_accessions():
    cluster_accessions = {}
    with open('cluster_accessions', 'rt') as accessions_file:
        for line in accessions_file:
            parts = line.strip().split('\t')
            accessions = [x.split('.fna.gz')[0][:13] for x in parts[1].split(',')
                          if x.endswith('*')]
            assert len(accessions) == 1
            cluster_accessions[parts[0]] = accessions[0]
    return cluster_accessions


def load_accession_species():
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
    if pathlib.Path('species_definitions').is_file():
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

    # If the user-defined file doesn't exist yet, make an empty one with instructions.
    else:
        with open('species_definitions', 'wt') as user_species:
            t = '# This file is where you can define species for particular RefSeq assemblies,\n' \
                '# overriding the RefSeq species labels. Simply add lines to this file which\n' \
                '# have the RefSeq assembly accession (starts with GCF) followed by a tab and\n' \
                '# then the binomial species name. The version number (e.g. final .1) does not\n' \
                'need to be included in the accession.\n' \
                '\n' \
                '# Example:\n' \
                'GCF_000000000	Genus species\n'
            user_species.write(t)

    # Shigella and E. coli are grouped together.
    for accession, species in accession_species.items():
        if species == 'Escherichia coli' or species.startswith('Shigella '):
            accession_species[accession] = 'Escherichia coli / Shigella'

    return accession_species


if __name__ == '__main__':
    main()
