#!/usr/bin/env Rscript

# Install ape and phangorn if they aren't already installed.
list_of_packages <- c("ape", "phangorn")
new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

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

tree <- midpoint(tree)
write.tree(tree, newick_tree_filename)
