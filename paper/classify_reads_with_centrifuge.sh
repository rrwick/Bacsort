# Grab the sample names from the read files
sample_names=$(for f in "$TEST_READS"/*_1.fastq.gz; do echo ${f##*/} | sed 's/_1.fastq.gz//'; done)

# Run classification for both Centrifuge databases
for db_name in bacteria bacsort; do
    
    mkdir -p "$CLASSIFICATION"/centrifuge/"$db_name"
    cd "$CLASSIFICATION"/centrifuge/"$db_name"
    db="$DATABASES"/centrifuge_databases/"$db_name"_"$DATABASE_SUFFIX"

    # Run Centrifuge
    for sample in $(echo $sample_names); do
        r1="$TEST_READS"/"$sample"_1.fastq.gz
        r2="$TEST_READS"/"$sample"_2.fastq.gz
        sbatch -p "$SLURM_PARTITION" --nodes=1 --job-name=centrifuge --ntasks=1 --cpus-per-task=4 --mem=40000 --time=0-1:0:00 --wrap "/usr/bin/time -v centrifuge -x "$db" -1 "$r1" -2 "$r2" --report-file "$sample".report -S "$sample".output --threads 4; centrifuge-kreport -x "$db" "$sample".output > "$sample".report; rm "$sample".output"
        sleep 0.1  # slow down because the server can mysteriously lose jobs if you submit them too fast
    done

    # WAIT HERE FOR THE SLURM JOBS TO FINISH

    # Get top species
    for sample in $(echo $sample_names); do
        printf $sample"\t"
        "$BACSORT"/paper/get_top_kraken_species.py "$sample".report
    done

    # Get time and RAM
    for job in slurm-*; do
        time=$(grep "Elapsed (wall clock) time" $job | grep -oP "[\d\.:]{4,}")
        memory=$(grep "Maximum resident set size" $job | grep -oP "\d+")
        echo $time"\t"$memory
    done

done
