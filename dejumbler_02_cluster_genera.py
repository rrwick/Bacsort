#!/usr/bin/env python3

import sys
import collections
import os
import gzip
import shutil
from pathlib import Path

threshold = 0.0025


def main():
    genera = sys.argv[1].split()
    print()

    for genus in genera:
        print('Clustering ' + genus)
        print('------------------------------------------------')

        distance_filename = 'assemblies/' + genus + '/mash_distances'
        if not Path(distance_filename).is_file():
            sys.exit('Error: could not find pairwise distances file')

        assemblies, graph = create_graph_from_distances(distance_filename)
        clusters = cluster_assemblies(assemblies, graph, threshold)

        if not Path('clusters').is_dir():
            os.makedirs('clusters')

        print()
        cluster_num_digits = len(str(len(clusters)))
        cluster_num_format = '%0' + str(cluster_num_digits) + 'd'
        with open('cluster_accessions', 'at') as cluster_accessions_file:
            for num, assemblies in clusters.items():
                cluster_name = genus + '_' + (cluster_num_format % num)

                # Choose representative for each cluster by N50.
                # TO DO: maybe choose using some other logic?
                if len(assemblies) == 1:
                    representative = assemblies[0]
                else:
                    representative = sorted([(get_assembly_n50('assemblies/' + genus + '/' + a), a) for a in assemblies])[-1][1]

                cluster_info = '\t'.join([cluster_name, ','.join(assemblies), representative])
                cluster_accessions_file.write(cluster_info)
                cluster_accessions_file.write('\n')

                print(cluster_name + '\t', end='')
                print(','.join([(a + '*' if a == representative else a) for a in assemblies]))

                shutil.copyfile('assemblies/' + genus + '/' + representative,
                                'clusters/' + cluster_name + '.fna.gz')

            # Delete the assemblies (to save disk space)
            pass

        print('\n')



def create_graph_from_distances(distance_filename):
    print('Loading distances...', end='', flush=True)

    assemblies = set()
    graph = collections.defaultdict(set)
    all_connections = collections.defaultdict(set)

    with open(distance_filename, 'rt') as distance_file:
        for line in distance_file:
            parts = line.split('\t')
            assembly_1 = parts[0]
            assembly_2 = parts[1]
            distance = float(parts[2])

            if assembly_1 == assembly_2:
                continue

            assemblies.add(assembly_1)
            assemblies.add(assembly_2)

            all_connections[assembly_1].add(assembly_2)
            all_connections[assembly_2].add(assembly_1)

            if distance < threshold:
                graph[assembly_1].add(assembly_2)
                graph[assembly_2].add(assembly_1)

    assemblies = sorted(assemblies)
    assembly_count = len(assemblies)

    print(' found', assembly_count, 'assemblies')

    # Sanity check: make sure we have all the connections.
    for assembly in assemblies:
        assert len(all_connections[assembly]) == assembly_count - 1

    return assemblies, graph


def cluster_assemblies(assemblies, graph, threshold):
    visited = set()
    i = 0
    clusters = {}
    for assembly in assemblies:
        if assembly in visited:
            continue
        i += 1
        connected = dfs(graph, assembly)
        clusters[i] = sorted(connected)
        visited |= connected
    return clusters


def dfs(graph, start):
    visited, stack = set(), [start]
    while stack:
        vertex = stack.pop()
        if vertex not in visited:
            visited.add(vertex)
            stack.extend(graph[vertex] - visited)
    return visited


def get_assembly_n50(filename):
    contig_lengths = sorted(get_contig_lengths(filename), reverse=True)
    total_length = sum(contig_lengths)
    target_length = total_length * 0.5
    length_so_far = 0
    for contig_length in contig_lengths:
        length_so_far += contig_length
        if length_so_far >= target_length:
            return contig_length
    return 0


def get_contig_lengths(filename):
    lengths = []
    with gzip.open(filename, 'rt') as fasta_file:
        name = ''
        sequence = ''
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    lengths.append(len(sequence))
                    sequence = ''
                name = line[1:].split()[0]
            else:
                sequence += line
        if name:
            lengths.append(len(sequence))
    return lengths


if __name__ == '__main__':
    main()
