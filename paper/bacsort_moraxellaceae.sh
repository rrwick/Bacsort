mkdir -p "$BACSORT_GENOMES"/moraxellaceae
cd "$BACSORT_GENOMES"/moraxellaceae

# Download genomes.
"$BACSORT"/download_genomes.sh "Acinetobacter Alkanindiges Cavicella Faucicola Fluviicoccus Moraxella Paraperlucidibaca Perlucidibaca Psychrobacter" > download.out 2>&1

# Copy my previous exclusions and definitions files so they will be used.
cp "$BACSORT"/excluded_assemblies .
cp "$BACSORT"/species_definitions .

# Create genome clusters to reduce the number of assemblies for tree-building.
"$BACSORT"/scripts/cluster_genera.py assemblies

# Distance matrices
"$BACSORT"/scripts/mash_distance_matrix.sh 80
"$BACSORT"/scripts/fastani_in_parallel.sh 80
"$BACSORT"/scripts/pairwise_identities_to_distance_matrix.py tree/fastani_output > tree/fastani.phylip
"$BACSORT"/scripts/combine_distance_matrices.py tree/fastani.phylip tree/mash.phylip > tree/combined.phylip

# Build the tree.
"$BACSORT"/bionj_tree.R tree/combined.phylip tree/tree.newick

# Annotate species in the tree.
# This script assumes an input file of tree/tree.newick
"$BACSORT"/scripts/find_species_clades.py

# Then I refine species definitions and re-run these commands, as many times
# as necessary:
cp "$BACSORT"/species_definitions .
"$BACSORT"/scripts/find_species_clades.py

# Copy curated assemblies to species directories.
"$BACSORT"/scripts/copy_assemblies.py
"$BACSORT"/scripts/copy_clusters.py

# Run genome painter on bins with enough (5 or more) assemblies.
mkdir -p genome_painter
cd genome_painter
mkdir -p counts database results
for d in ../assemblies_binned/*/*; do
    genus=$(echo "$d" | sed 's|../assemblies_binned/||' | sed -r 's|/.+||')
    species=$(echo "$d" | sed 's|../assemblies_binned/||' | sed -r 's|.+/||')
    if [ "$species" != "unknown" ]; then
        count=$(ls ../assemblies_binned/"$genus"/"$species"/*.fna.gz | wc -l)
        if (( $count >= 5)); then
            echo $genus $species": "$count" assemblies"
            count_kmers --threads 16 --genome_fps "$d"/*.fna.gz --output_fp counts/"$genus"_"$species".tsv
        fi
    fi
done
generate_database --alpha 0.0 --threshold 0.95 --count_fps counts/*.tsv --output_fp database/painter_database.bin
for d in ../assemblies_binned/*/*; do
    genus=$(echo "$d" | sed 's|../assemblies_binned/||' | sed -r 's|/.+||')
    species=$(echo "$d" | sed 's|../assemblies_binned/||' | sed -r 's|.+/||')
    if [ "$species" != "unknown" ]; then
        count=$(ls ../assemblies_binned/"$genus"/"$species"/*.fna.gz | wc -l)
        if (( $count >= 5)); then
            paint_genome --genome_fp ../assemblies_binned/"$genus"/"$species"/*.fna.gz --kmer_db_fp database/painter_database.bin --output_dir results
            summarise_species.py --threads 16 results > "$genus"_"$species".tsv.gz
            rm results/*.tsv.gz
        fi
    fi
done
