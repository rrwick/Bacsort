#!/usr/bin/env bash

# Copyright 2018 Ryan Wick (rrwick@gmail.com)
# https://github.com/rrwick/Dejumbler

# This script is the third step of Dejumbler. It takes no arguments. When run, it uses Mash to
# build a distance matrix between all assemblies. It then builds a neighbour joining tree with
# Quicktree (github.com/khowe/quicktree).

# This file is part of Dejumbler. Dejumbler is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later version. Dejumbler is
# distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details. You should have received a copy of the GNU General Public
# License along with Dejumbler. If not, see <http://www.gnu.org/licenses/>.

threads=16

mkdir -p tree

# create a distance matrix
cd clusters
mash sketch -p $threads -o ../tree/reference -s 100000 *.fna.gz
cd ..
mash dist -p $threads tree/reference.msh tree/reference.msh -t > tree/distances.tab

# create phylip format file of the distance matrix
tail -n +2 tree/distances.tab > tree/distances.tab.temp  # remove first line
wc -l tree/distances.tab.temp | awk '"[0-9]+ errors" {sum += $1}; END {print sum}' > tree/distances.ndist  # get number for sample/line count
cat tree/distances.ndist tree/distances.tab.temp > tree/distances.phylip  # write matrix

# make tree from phylip distance matrix
quicktree -in m tree/distances.phylip > tree/tree.newick
