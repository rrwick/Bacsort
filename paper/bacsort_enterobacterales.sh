mkdir -p "$BACSORT_GENOMES"/enterobacterales
cd "$BACSORT_GENOMES"/enterobacterales

# Download genomes. This list includes all genera in Enterobacterales,
# including ones with no assemblies (currently anyway), but they don't cause
# any issues.
"$BACSORT"/scripts/download_genomes.sh "Annandia Aquamonas Aranicola Arsenophonus Aschnera Atlantibacter Averyella Benitsuchiphilus Biostraticola Brenneria Buchnera Budvicia Buttiauxella Cedecea Chania Citrobacter Coetzeea Cosenzaea Cronobacter Curculioniphilus Cuticobacterium Dickeya Edwardsiella Enterobacillus Enterobacter Erwinia Escherichia Ewingella Franconibacter Fukatsuia Gibbsiella Gillettellia Grimontella Guhaiyinggella Hafnia Hamiltonella Hartigia Ishikawaella Izhakiella Klebsiella Kleidoceria Kluyvera Kosakonia Leclercia Lelliottia Leminorella Levinea Limnobaculum Lonsdalea Macropleicola Mangrovibacter Margalefia Metakosakonia Moellerella Moranella Morganella Nissabacter Obesumbacterium Pantoea Pectobacterium Phaseolibacter Phlomobacter Photorhabdus Phytobacter Plesiomonas Pluralibacter Pragia Profftia Proteus Providencia Pseudescherichia Pseudocitrobacter Purcelliella Rahnella Raoultella Regiella Riesia Rohrkolberia Rosenbergiella Rosenkranzia Rouxiella Saccharobacter Salmonella Samsonia Schneideria Serratia Shigella Shimwellia Siccibacter Sodalis Stammerula Tachikawaea Tatumella Thorsellia Tiedjeia Trabulsiella Wigglesworthia Xenorhabdus Yersinia Yokenella" > download.out 2>&1

# Two Enterobacter genomes (GCF_000164865.1 and GCF_001461805.1) didn't
# download because they were listed as '[Enterobacter] lignolyticus' (the
# brackets confused ncbi-genome-download). I therefore added those two genomes
# manually to the Enterobacter directory, added the data lines (without the
# brackets in the genus) and rebuilt the Mash distances.
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/164/865/GCF_000164865.1_ASM16486v1/GCF_000164865.1_ASM16486v1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/461/805/GCF_001461805.1_ASM146180v1/GCF_001461805.1_ASM146180v1_genomic.fna.gz
mv GCF_000164865.1_ASM16486v1_genomic.fna.gz assemblies/Enterobacter/GCF_000164865.1.fna.gz
mv GCF_001461805.1_ASM146180v1_genomic.fna.gz assemblies/Enterobacter/GCF_001461805.1.fna.gz
cd assemblies/Enterobacter
echo "GCF_000164865.1\tPRJNA224116\tSAMN00116754\t\t\trepresentative genome\tassembly from type material\t701347\t1334193\tEnterobacter lignolyticus SCF1\tstrain=SCF1\t\tlatest\tComplete Genome\tMajor\tFull\t2010/10/15\tASM16486v1\tUS DOE Joint Genome Institute\tGCA_000164865.1\tidentical\tftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/164/865/GCF_000164865.1_ASM16486v1\t./refseq/bacteria/GCF_000164865.1/GCF_000164865.1_ASM16486v1_genomic.fna.gz" >> data.tsv
echo "GCF_001461805.1\tPRJNA224116\tSAMN04074888\t\t\tna\t\t1334193\t1334193\tEnterobacter lignolyticus\tstrain=G5\t\tlatest\tComplete Genome\tMajor\tFull\t2015/12/07\tASM146180v1\tUniversity of Malaya\tGCA_001461805.1\tidentical\tftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/461/805/GCF_001461805.1_ASM146180v1\t./refseq/bacteria/GCF_001461805.1/GCF_001461805.1_ASM146180v1_genomic.fna.gz" >> data.tsv
rm mash_distances mash.msh
mash sketch -p 16 -o mash -s 10000 *.fna.gz
mash dist -p 16 mash.msh mash.msh > mash_distances
cd ../..

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
"$BACSORT"/scripts/bionj_tree.R tree/combined.phylip tree/tree.newick

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

# Run Genome Painter (https://github.com/scwatts/genome_painter) on assemblies
# to identify contaminated assemblies and cross-species hybrids to exclude.
mkdir -p genome_painter
cd genome_painter
mkdir -p counts database results
for d in ../assemblies_binned/*/*; do
    genus=$(echo "$d" | sed 's|../assemblies_binned/||' | sed -r 's|/.+||')
    species=$(echo "$d" | sed 's|../assemblies_binned/||' | sed -r 's|.+/||')
    if [ "$species" != "unknown" ]; then
        count=$(ls ../assemblies_binned/"$genus"/"$species"/*.fna.gz | wc -l)
        if (( $count >= 5)); then  # Only use bins with 5+ assemblies to keep db size down
            echo $genus $species": "$count" assemblies"
            count_kmers --threads 16 --genome_fps "$d"/*.fna.gz --output_fp counts/"$genus"_"$species".tsv
        fi
    fi
done
generate_database --alpha 0.0 --threshold 0.99 --count_fps counts/*.tsv --output_fp database/painter_database.bin
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
