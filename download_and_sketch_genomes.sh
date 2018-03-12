#!/usr/bin/env bash

# Copyright 2018 Ryan Wick (rrwick@gmail.com)
# https://github.com/rrwick/Dejumbler

# This script is the first step of Dejumbler. It takes a single argument: a space-delimited list
# of genera. When run, it downloads assemblies from NCBI with their metadata (including the
# species). It then uses Mash to do pairwise distances between all assemblies in each genus.

# This file is part of Dejumbler. Dejumbler is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later version. Dejumbler is
# distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details. You should have received a copy of the GNU General Public
# License along with Dejumbler. If not, see <http://www.gnu.org/licenses/>.


genera=$1
mash_sketch_size=10000
threads=16

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
    $download --verbose --genus $genus --metadata-table data.tsv --format fasta -p $threads bacteria

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
        mash sketch -p $threads -o mash -s $mash_sketch_size *.fna.gz
        mash dist -p $threads mash.msh mash.msh > mash_distances
        cd ../..
    else
        echo "No assemblies downloaded for "$genus
        cd ..
        rm -r $genus
        cd ..
    fi
    printf "\n"

done
