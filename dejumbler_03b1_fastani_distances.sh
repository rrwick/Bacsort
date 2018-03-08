#!/usr/bin/env bash

# Copyright 2018 Ryan Wick (rrwick@gmail.com)
# https://github.com/rrwick/Dejumbler

# This script is the third step of Dejumbler. It takes no arguments. When run, it uses FastANI to
# build a distance matrix between all assemblies.

# This file is part of Dejumbler. Dejumbler is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later version. Dejumbler is
# distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details. You should have received a copy of the GNU General Public
# License along with Dejumbler. If not, see <http://www.gnu.org/licenses/>.

threads=32

mkdir -p tree

printf "\n"
echo "Find pairwise distances with FastANI"
echo "------------------------------------------------"
cd clusters
ls *.fna.gz > cluster_list


# Option 1: run FastANI in one thread.
# This is the simplest option, but is slow.
fastANI --rl cluster_list --ql cluster_list -o ../tree/fastani_output
exit 0



# Options 2 and 3 rely on splitting up the cluster list.
total_count=$( wc -l < cluster_list )
count_per_file=$( perl -w -e "use POSIX; print ceil($total_count/$threads), qq{\n}" )
cat cluster_list | shuf > split.tmp
split -a 4 -dl $count_per_file split.tmp cluster_list_
rm split.tmp



# Option 2: run FastANI in parallel.
for q in cluster_list_*; do
    q_num=$(echo $q | sed 's/cluster_list_//')
    for r in cluster_list_*; do
        r_num=$(echo $r | sed 's/cluster_list_//')
        FastANI --rl $r --ql $q -o ../tree/fastani_output_"$q_num"_"$r_num" &> ../tree/fastani_stdout_"$q_num"_"$r_num" &
    done
    wait
done



# Option 3: run FastANI in parallel on a Slurm-managed cluster.
for q in cluster_list_*; do
    q_num=$(echo $q | sed 's/cluster_list_//')
    for r in cluster_list_*; do
        r_num=$(echo $r | sed 's/cluster_list_//')
        sbatch --nodes=1 --job-name=fastANI_"$q_num"_"$r_num" --ntasks=1 --cpus-per-task=1 --mem=4096 --time=0-24:0:00 --wrap "fastANI --rl "$r" --ql "$q" -o ../tree/fastani_output_"$q_num"_"$r_num" &> ../tree/fastani_stdout_"$q_num"_"$r_num
    done
done



# If option 2 or 3 was used, we can now must concatentate all results together
# and can clean up the intermediate files.
cd ..
cat tree/fastani_output_* > tree/fastani_output
rm tree/fastani_output_*
