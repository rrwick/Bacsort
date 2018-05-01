mkdir -p "$DATABASES"
cd "$DATABASES"

# Copy the Bacsorted clusters to one directory.
mkdir -p bacsorted_assemblies
cp -r "$BACSORT_GENOMES"/enterobacterales/clusters_binned/* bacsorted_assemblies
cp -r "$BACSORT_GENOMES"/moraxellaceae/clusters_binned/* bacsorted_assemblies
cd ..

mkdir -p kraken_databases
cd kraken_databases

# This script will build two databases: one using the standard bacterial
# genomes, and one that swaps out Bacsorted genomes where it can.
db_1=bacteria_"$DATABASE_SUFFIX"
db_2=bacsort_"$DATABASE_SUFFIX"

# First we prepare genomes for a bacterial database.
kraken-build --download-taxonomy --db $db_1
kraken-build --download-library bacteria --db $db_1

# Now make a copy of the entire thing, so we can modify one with Bacsort.
cp -r $db_1 $db_2

# Build the first database.
module load Jellyfish/1.1.11-GCC-4.9.3
kraken-build --build --max-db-size 64 --threads 16 --db $db_1 > "$db_1"_build.out 2>&1

# Clean up.
cd $db_1
rm accmap_file.tmp database.jdb database.jdb.big lca.complete seqid2taxid.map
rm -r library
cd taxonomy
rm accmap.dlflag citations.dmp delnodes.dmp division.dmp gc.prt gencode.dmp merged.dmp nucl_est.accession2taxid nucl_gb.accession2taxid nucl_gss.accession2taxid nucl_wgs.accession2taxid prelim_map.txt readme.txt taxdump.dlflag taxdump.tar.gz taxdump.untarflag
cd ../..















# Now make the Bacsort modifications and build the second.
"$BACSORT"/prepare_kraken_library.py ../bacsorted_assemblies $db_2
for f in additional_assemblies/*.fna; do
    kraken-build --add-to-library $f --db $DBNAME
done
kraken-build --build --max-db-size 64 --db $db_2
