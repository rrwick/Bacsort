#!/usr/bin/env Rscript

library(ape)
library(phangorn)

args = commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Error: you must supply two arguments: PHYLIP distance matrix filename (input) and a Newick tree filename (output)", call.=FALSE)
}

distance_matrix_filename = args[1]
newick_tree_filename = args[2]

cat("Loading distance matrix...")
distances <- readDist(distance_matrix_filename)
cat(" done\n")

cat("Building BIONJ tree...")
tree <- bionj(distances)
cat(" done\n")

write.tree(tree, newick_tree_filename)

