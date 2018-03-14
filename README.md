# Bacterial species Dejumbler


## Requirements

Running Dejumbler requires that you have [Mash](https://github.com/marbl/Mash) installed and available in your PATH. If you can type `mash -h` into your terminal and not get an error, you should be good! To build trees, you'll need either [Quicktree](https://github.com/khowe/quicktree) or [RapidNJ](http://birc.au.dk/software/rapidnj/).

You'll also need Python3 and [BioPython](http://biopython.org/). If `python3 -c "import Bio"` doesn't give you an error, you should be good! If you need to install BioPython, it's easiest to do with pip: `pip3 install biopython`

Finally, depending on how you want to compute pairwise distances, you may need [FastANI](https://github.com/ParBLiSS/FastANI), [minimap2](https://github.com/lh3/minimap2) and/or the [edlib Python library](https://github.com/Martinsos/edlib/tree/master/bindings/python).




## Installation

You don't need to install Dejumbler - it's just a collection of independent scripts which can be run directly. However, copying them to somewhere in your path (e.g. `~/.local/bin`) will make it easier to run:

```
git clone https://github.com/rrwick/Dejumbler
install_dir="$HOME"/.local/bin
mkdir -p "$install_dir"
cp Dejumbler/dejumbler_* "$install_dir"
export PATH="$install_dir":"$PATH"  # Add this line to your .bashrc file (or equivalent) if needed
```




## How to use

All Dejumbler commands need to be run in the same directory. They will create directories and files where it is run, so I would recommended running it from a new directory:

```
mkdir dejumbler_results
cd dejumbler_results
```

Dejumbler will download _all_ NCBI assemblies for your genera of interest, so make sure you have enough free disk space! Each assembly is about 2-3 MB (gzipped), so you'll need many GB for common genera with many assemblies.


### Step 1: download assemblies

```
download_genomes.sh "Citrobacter Klebsiella Salmonella Yersinia"
```

### Step 2: cluster assemblies

```
cluster_genera.py "Citrobacter Klebsiella Salmonella Yersinia"
```

### Step 3: distance matrix

As input for the neighbour joining tree algorith, we need a matrix of all pairwise distances between clusters. There a few alternative ways to produce such a matrix: using Mash, FastANI or rMLST.

#### Option 1: Mash

Advantages:
* Simple and fast.

Disadvantages:
* Since Mash uses the entire set of k-mers for each pairwise comparison, both core sequence (shared by both assemblies) and accessory sequence (only in one assembly) affect the result. This means it may produce less accurate (compared to vertical evolution) trees.

This one script will take care of running Mash and converting its output to a PHYLIP distance matrix. Run it with no arguments in your working Dejumbler directory.
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

To run FastANI on your Dejumbler clusters, choose one of these scripts. You may need to modify them to suit your system (particularly the script which uses Slurm).
```
run_fastani_one_thread.sh
run_fastani_in_parallel.sh
run_fastani_with_slurm.sh
```

Once the distances are computed, they must be converted into a PHYLIP distance matrix, which is relatively quick and carried out using this command. We use a maximum distance of 0.2 because FastANI wasn't designed to quantify ANI less than 80%.
```
pairwise_identities_to_distance_matrix.py --max_dist 0.2 tree/fastani_output > tree/distances.phylip
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

### Step 4: build tree

Building the tree with [Quicktree](https://github.com/khowe/quicktree) is relatively quick and easy:
```
quicktree -in m tree/distances.phylip > tree/tree.newick
```

Alternatively, you can use [RapidNJ](http://birc.au.dk/software/rapidnj/):
```
rapidnj -i pd tree/distances.phylip > tree/tree.newick
```

### Step 5: curate tree

```
dejumbler_05_find_species_clades.py
```

### Step 6: copy assemblies to species directories

```
dejumbler_06_copy_assemblies.py
```




## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
