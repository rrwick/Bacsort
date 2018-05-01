mkdir -p "$DATABASES"/centrifuge_databases
cd "$DATABASES"/centrifuge_databases

# This script will build two databases: one using the standard bacterial
# genomes, and one that swaps out Bacsorted genomes where it can.
db_1=bacteria_"$DATABASE_SUFFIX"
db_2=bacsort_"$DATABASE_SUFFIX"

# First, we build the standard bacteral database.
centrifuge-download -o taxonomy taxonomy
module load BLAST+/2.2.31-GCC-4.9.2  # for dustmasker
centrifuge-download -o library -m -d "bacteria" refseq > seqid2taxid.map
cat library/*/*.fna > input-sequences.fna
centrifuge-build -p 16 --conversion-table seqid2taxid.map --taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp input-sequences.fna $db_1 > "$db_1"_build.out 2>&1

# Clean up 
rm input-sequences.fna




















# Now we make Bacsort's modifications and build again.
rm input-sequences.fna
"$BACSORT"/prepare_centrifuge_library.py ../bacsorted_assemblies .
cat library/*/*.fna > input-sequences.fna
centrifuge-build -p 16 --conversion-table seqid2taxid.map --taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp input-sequences.fna $db_2

# Clean up
rm input-sequences.fna seqid2taxid.map
rm -r library taxonomy
