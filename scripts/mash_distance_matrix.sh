#!/usr/bin/env bash

# Copyright 2018 Ryan Wick (rrwick@gmail.com)
# https://github.com/rrwick/Bacsort

# This script is the third step of Bacsort. It takes no arguments. When run, it uses Mash to build
# a PHYLIP distance matrix between all assemblies.

# This script runs Mash to build a distance matrix between all assembly clusters. It takes one
# argument: the number of threads.

# This file is part of Bacsort. Bacsort is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free Software Foundation,
# either version 3 of the License, or (at your option) any later version. Bacsort is distributed
# in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# more details. You should have received a copy of the GNU General Public License along with
# Bacsort. If not, see <http://www.gnu.org/licenses/>.

mkdir -p tree

printf "\n"
echo "Mash distance matrix of all clusters"
echo "------------------------------------------------"
cd clusters
mash sketch -p $1 -o ../tree/reference -s 100000 *.fna.gz
cd ..
mash dist -p $1 tree/reference.msh tree/reference.msh -t > tree/distances.tab

# Convert matrix to PHYLIP format
tail -n +2 tree/distances.tab > tree/distances.tab.temp  # remove first line
wc -l tree/distances.tab.temp | awk '"[0-9]+ errors" {sum += $1}; END {print sum}' > tree/distances.ndist  # get number for sample/line count
cat tree/distances.ndist tree/distances.tab.temp > tree/mash.phylip
rm tree/distances.ndist tree/reference.msh tree/distances.tab tree/distances.tab.temp
printf "\n"
