# Bacsort


## Requirements

Running Bacsort requires that you have [Mash](https://github.com/marbl/Mash) installed and available in your PATH. If you can type `mash -h` into your terminal and not get an error, you should be good! To build trees, you'll need either [Quicktree](https://github.com/khowe/quicktree) or [RapidNJ](http://birc.au.dk/software/rapidnj/).

You'll also need Python3 and [BioPython](http://biopython.org/). If `python3 -c "import Bio"` doesn't give you an error, you should be good! If you need to install BioPython, it's easiest to do with pip: `pip3 install biopython`

Finally, depending on how you want to compute pairwise distances, you may need [FastANI](https://github.com/ParBLiSS/FastANI), [minimap2](https://github.com/lh3/minimap2) and/or the [edlib Python library](https://github.com/Martinsos/edlib/tree/master/bindings/python).




## Installation

You don't need to install Bacsort - it's just a collection of independent scripts which can be run directly. However, copying them to somewhere in your path (e.g. `~/.local/bin`) will make it easier to run:

```
git clone https://github.com/rrwick/Bacsort
install_dir="$HOME"/.local/bin
mkdir -p "$install_dir"
cp Bacsort/*.sh Bacsort/*.py "$install_dir"
export PATH="$install_dir":"$PATH"  # Add this line to your .bashrc file (or equivalent) if needed
```




## How to use

All Bacsort commands need to be run in the same directory. They will create directories and files where it is run, so I would recommended running it from a new directory:

```
mkdir bacsort_results
cd bacsort_results
```

Bacsort will download _all_ NCBI assemblies for your genera of interest, so make sure you have enough free disk space! Each assembly is about 2-3 MB (gzipped), so you'll need many GB for common genera with many assemblies.


### Step 1: download assemblies

```
download_genomes.sh "Citrobacter Klebsiella Salmonella Yersinia"
```

### Step 2: cluster assemblies

```
cluster_genera.py assemblies
```

### Step 3: distance matrix

As input for the neighbour joining tree algorithm, we need a matrix of all pairwise distances between clusters. There a few alternative ways to produce such a matrix: using Mash, FastANI or rMLST.

#### Option 1: Mash

Advantages:
* Simple and fast.

Disadvantages:
* Since Mash uses the entire set of k-mers for each pairwise comparison, both core sequence (shared by both assemblies) and accessory sequence (only in one assembly) affect the result. This means it may produce less accurate (compared to vertical evolution) trees.

This one script will take care of running Mash and converting its output to a PHYLIP distance matrix. Run it with no arguments in your working Bacsort directory.
```
mash_distance_matrix.sh
```

#### Option 2: FastANI

Advantages:
* Produces pairwise ANI measurements using only the sequence shared by two assemblies. This makes it less swayed by the accessory genome and it may produce more accurate trees.

Disadvantages:
* FastANI is designed for 80-100% nucleotide identity, and so this approach is only recommended if you are working within a genus or closely related genera.
* It is slower than Mash.
* FastANI is not intrinsically parallel, so you'll need to parallelise it one way or another to cope with large datasets.

To run FastANI on your Bacsort clusters in serial (only appropriate for small datasets):
```
fastani_one_thread.sh
```

Or use this script to run it in parallel:
```
fastani_in_parallel.sh
```

Or if you have a Slurm-managed cluster, this may be the fastest approach:
```
fastani_with_slurm.sh
# Wait for Slurm jobs to finish
cat tree/fastani_output_* > tree/fastani_output
rm tree/fastani_output_* tree/fastani_stdout_*
```

Once the distances are computed, they must be converted into a PHYLIP distance matrix, which is relatively quick and carried out using this command. We use a maximum distance of 0.2 because FastANI wasn't designed to quantify ANI less than 80%.
```
pairwise_identities_to_distance_matrix.py --max_dist 0.2 tree/fastani_output > tree/fastani_distances.phylip
```

#### Option 3: rMLST

Advantages:
* By using a mostly shared core genome, this approach should result in a tree that's decently compatible with vertical evolution.
* Can deal with much more sequence divergence than FastANI, and is therefore suitable for more diverse sets.

Disadvantages:
* You'll need to get the rMLST alleles yourself ([pubmlst.org/rmlst](https://pubmlst.org/rmlst/) - account required)

This approach has three steps. First, you must extract the rMLST sequences from your assemblies:
```
rmlst_extract_sequences.py clusters rmlst_allele_dir
```

Next, conduct pairwise alignments to get identities. Double check that this script produced a file with (n^2 + n) / 2 lines (all non-redundant pairwise comparisons).
```
rmlst_identities.py clusters > tree/rmlst_identities
```

Finally, convert the identities into a distance matrix (same as for FastANI):
```
pairwise_identities_to_distance_matrix.py tree/rmlst_identities > tree/distances.phylip
```

#### Option 4: combination

Finally, you can create two different distance matrices and combine them together. I wrote this for combining FastANI distances (which are very good up to 20% divergence) with Mash/rMLST distances (which can handle greater divergence).

```
combine_distance_matrices.py fastani.phylip mash.phylip > combined.phylip
```



### Step 4: build tree

There are many tools available for building a Newick-format tree from a PHYLIP distance matrix, and any can be used here.

I prefer the [BIONJ algorithm](https://academic.oup.com/mbe/article/14/7/685/1119804) based on the recommendation of the [this paper](https://wellcomeopenresearch.org/articles/3-33). A simple script is included in this repo to build one from a PHYLIP distance matrix using R's [ape package](http://ape-package.ird.fr/):
```
bionj_tree.R tree/distances.phylip tree/tree.newick
```

Alternatively, there are plenty of other tools for building a Newick-format tree from a PHYLIP distance matrix:
* [RapidNJ](http://birc.au.dk/software/rapidnj/) is a particularly fast tool using the neighbour-joining (NJ) algorithm.
* [Quicktree](https://github.com/khowe/quicktree) is another NJ tool which produces similar results.
* [BIONJ](http://www.atgc-montpellier.fr/bionj/download.php) is a commmand line tool for building BIONJ trees.
* The [ape](http://ape-package.ird.fr/) and [phangorn](https://github.com/KlausVigo/phangorn) R packages have a number of tree-building algorithms, including UPGMA, NJ and BIONJ.
* The [PHYLIP package](http://evolution.genetics.washington.edu/phylip.html) has a number of programs which may be relevant.


### Step 5: curate tree

Run this command to make a tree file which contains the species labels in the node names:
```
find_species_clades.py
```

It produces the tree in two formats: newick and PhyloXML. The PhyloXML tree is suitable for viewing in [Archaeopteryx](https://sites.google.com/site/cmzmasek/home/software/archaeopteryx). Clades which perfectly define a species (i.e. all instances of that species are in a clade which contains no other species) will be coloured in the tree. You can then focus on the uncoloured parts where mislabellings may occur.

As you find species label errors, add them to the `species_definitions` file in this format:
```
GCF_000000000	Genus species
```

You can then run `find_species_clades.py` again to generate a new tree with your updated definitions. This process can be repeated (fix labels, make tree, fix labels, make tree, etc) until you have no more changes to make.


### Step 6: copy assemblies to species directories

When you are happy with the species labels in the tree, you can run this command to copy assemblies into species directories. It will make a directory titled `assemblies_binned` with subdirectories for each genus and then subdirectories for each species.
```
copy_assemblies.py
```





## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
