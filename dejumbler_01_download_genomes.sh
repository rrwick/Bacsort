#!/usr/bin/env bash

genera=$1
mash_sketch_size=10000

printf "\n"


echo "Checking for Mash"
echo "------------------------------------------------"
if command -v mash; then
    echo "Mash found!"
else
    echo "Error: could not find mash"
    echo "install from https://github.com/marbl/Mash"
    exit 1
fi
printf "\n"


if [ ! -d ncbi-genome-download ]; then
    printf "\n"
    echo "Cloning ncbi-genome-download tool"
    echo "------------------------------------------------"
    git clone https://github.com/kblin/ncbi-genome-download
    printf "\n"
fi
download=$(pwd)/ncbi-genome-download/ncbi-genome-download-runner.py


for genus in $genera; do
    printf "\n"
    echo "Downloading genomes from "$genus
    echo "------------------------------------------------"
    mkdir -p assemblies/$genus
    cd assemblies/$genus
    $download --verbose --genus $genus --metadata-table data.tsv --format fasta bacteria

    if test -n "$(find refseq/bacteria -name '*.fna.gz' -print -quit)"
        then
        for f in refseq/bacteria/*/*.fna.gz; do
            new_name=$(echo $f | grep -oP "\w{3}_\d{9}\.\d" | head -n 1)".fna.gz"
            echo "mv "$f" "$new_name
            mv $f $new_name
        done
        rm -r refseq
        printf "\n\n"

        echo "Finding pairwise distances for "$genus
        echo "------------------------------------------------"
        mash sketch -o mash -s $mash_sketch_size *.fna.gz
        mash dist mash.msh mash.msh > mash_distances
        cd ../..
    else
        echo "No assemblies downloaded for "$genus
        cd ..
        rm -r $genus
        cd ..
    fi
    printf "\n"

done
