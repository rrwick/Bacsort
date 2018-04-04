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

threads=32

mkdir -p tree

printf "\n"
echo "Find pairwise distances with FastANI"
echo "------------------------------------------------"
cd clusters
ls *.fna.gz > cluster_list

fastANI --rl cluster_list --ql cluster_list -o ../tree/fastani_output
