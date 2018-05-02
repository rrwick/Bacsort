# Build Mash sketch
cd "$DATABASES"/bacsorted_assemblies
mash sketch -o mash_sketches -s 100000 -p 32 */*/*.fna.gz

mkdir -p "$CLASSIFICATION"/mash/assemblies
cd "$CLASSIFICATION"/mash/assemblies

# Grab the sample names from the read files
sample_names=$(for f in "$TEST_READS"/*_1.fastq.gz; do echo ${f##*/} | sed 's/_1.fastq.gz//'; done)

# Run Mash
sketch="$DATABASES"/bacsorted_assemblies/mash_sketches.msh
for sample in $(echo $sample_names); do
    assembly="$TEST_ASSEMBLIES"/"$sample".fasta
    sbatch -p sysgen --nodes=1 --job-name=mash_assembly --ntasks=1 --cpus-per-task=4 --mem=12000 --time=0-0:10:00 --wrap "/usr/bin/time -v "$BACSORT"/scripts/classify_using_mash.py "$sketch" "$assembly
    sleep 0.1
done

# Get top species
for job in slurm-*; do
    head -n 1 $job
done

# Get time and RAM
for job in slurm-*; do
    time=$(grep "Elapsed (wall clock) time" $job | grep -oP "[\d\.:]{4,}")
    memory=$(grep "Maximum resident set size" $job | grep -oP "\d+")
    echo $time"\t"$memory
done
