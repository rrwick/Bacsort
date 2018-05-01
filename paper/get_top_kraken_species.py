#!/usr/bin/env python3
"""
Copyright 2018 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Bacsort

This script prints the top species in a Kraken report. Doing so is a bit more complex than you
might expect, because there can be different levels of 'species', and we want the deepest one.

Here's the top of a Kraken report that is straightforward. The correct answer is 'Klebsiella
oxytoca', which is simply the first 'S' line:

  1.84  73825   73825   U   0       unclassified
 58.36  2341944 0       -   1       root
 58.36  2341944 0       -   131567    cellular organisms
 58.36  2341871 0       D   2           Bacteria
 58.33  2340848 0       P   1224          Proteobacteria
 58.30  2339472 0       C   1236            Gammaproteobacteria
 58.23  2336662 35      O   91347             Enterobacterales
 57.78  2318490 3894    F   543                 Enterobacteriaceae
 56.23  2256544 1110    G   570                   Klebsiella
 23.98  962213  962213  S   571                     Klebsiella oxytoca
 18.22  731142  731142  S   1134687                 Klebsiella michiganensis
 12.46  499908  499908  S   1905288                 Klebsiella sp. LTGPAF-6F
  1.38  55319   55319   S   573                     Klebsiella pneumoniae

Here's a more complex example (and the reason this script exists). The first 'S' line is
'Enterobacter cloacae complex', but that's not right answer (for this script's purposes, anyway)
because there is another level of species nested underneath it. The answer we want is instead
'Enterobacter hormaechei':

  1.17  61597   61597   U   0       unclassified
 51.48  2715943 0       -   1       root
 51.48  2715943 0       -   131567    cellular organisms
 51.48  2715942 0       D   2           Bacteria
 51.44  2713965 0       P   1224          Proteobacteria
 51.22  2702022 30      C   1236            Gammaproteobacteria
 50.53  2665534 502     O   91347             Enterobacterales
 50.14  2645044 51586   F   543                 Enterobacteriaceae
 41.42  2185093 33030   G   547                   Enterobacter
 40.76  2150151 0       S   354276                  Enterobacter cloacae complex
 29.08  1533870 1533870 S   158836                    Enterobacter hormaechei
  7.66  404017  404017  S   1296536                   Enterobacter xiangfangensis
  2.44  128461  128461  S   550                       Enterobacter cloacae

This file is part of Bacsort. Bacsort is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Bacsort is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Bacsort. If
not, see <http://www.gnu.org/licenses/>.
"""

import sys

# Gather up all the taxa labels in the first group of lines that are designated 'S' (for species).
top_species = []
with open(sys.argv[1], 'rt') as kraken_results:
    for line in kraken_results:
        parts = line.split('\t')
        if parts[3] != 'S':
            if not top_species:
                continue
            else:
                break
        top_species.append(parts[5])

# Find the depth of the deepest level in this group.
max_indent = 0
for species in top_species:
    indent = len(species) - len(species.lstrip(' '))
    max_indent = max(max_indent, indent)

# Return the first species in the group with the maximum indent.
for species in top_species:
    indent = len(species) - len(species.lstrip(' '))
    if indent == max_indent:
        print(species.strip())
        sys.exit(0)
