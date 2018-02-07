#!/usr/bin/env bash

mkdir -p tree

# create a distance matrix
cd clusters
mash sketch -o ../tree/reference -s 100000 *.fna.gz
cd ..
mash dist tree/reference.msh tree/reference.msh -t > tree/distances.tab

# create phylip format file of the distance matrix
tail -n +2 tree/distances.tab > tree/distances.tab.temp  # remove first line
wc -l tree/distances.tab.temp | awk '"[0-9]+ errors" {sum += $1}; END {print sum}' > tree/distances.ndist  # get number for sample/line count
cat tree/distances.ndist tree/distances.tab.temp > tree/distances.phylip  # write matrix

# make tree from phylip distance matrix
quicktree -in m tree/distances.phylip > tree/tree.newick
