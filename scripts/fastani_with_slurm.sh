#!/usr/bin/env bash

# Copyright 2018 Ryan Wick (rrwick@gmail.com)
# https://github.com/rrwick/Bacsort

# This script is the third step of Bacsort. It takes no arguments. When run, it uses FastANI to
# build a distance matrix between all assemblies.

# This file is part of Bacsort. Bacsort is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free Software Foundation,
# either version 3 of the License, or (at your option) any later version. Bacsort is distributed
# in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# more details. You should have received a copy of the GNU General Public License along with
# Bacsort. If not, see <http://www.gnu.org/licenses/>.

# Check for fastANI
if command -v fastANI >/dev/null 2>&1 ; then
    echo "fastANI found"
    echo "version: $(fastANI -v)"
else
    echo "Error: fastANI not found - please install fastANI and try again"
    exit
fi

group_count=32

mkdir -p tree

printf "\n"
echo "Find pairwise distances with FastANI"
echo "------------------------------------------------"
cd clusters
ls *.fna.gz > cluster_list

total_count=$( wc -l < cluster_list )
count_per_file=$( perl -w -e "use POSIX; print ceil($total_count/$group_count), qq{\n}" )
cat cluster_list | shuf > split.tmp
split -a 4 -dl $count_per_file split.tmp cluster_list_
rm split.tmp

for q in cluster_list_*; do
    q_num=$(echo $q | sed 's/cluster_list_//')
    for r in cluster_list_*; do
        r_num=$(echo $r | sed 's/cluster_list_//')
        sbatch --nodes=1 --job-name=fastANI_"$q_num"_"$r_num" --ntasks=1 --cpus-per-task=1 --mem=4096 --time=0-24:0:00 --wrap "fastANI --rl "$r" --ql "$q" -o ../tree/fastani_output_"$q_num"_"$r_num" &> ../tree/fastani_stdout_"$q_num"_"$r_num
    done
done

cd ..
