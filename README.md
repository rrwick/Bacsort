# Bacsort


## Intro

RefSeq is a wonderful public repository of bacterial genome assemblies, but unfortunately many of its assemblies are mislabelled. This means that you cannot simply use a tool like [ncbi-genome-download](https://github.com/kblin/ncbi-genome-download) to get all assemblies of a particular genus, e.g. _Klebsiella_. If you did so, you would certainly get many _Klebsiella_ assemblies, but you would also get some _E. coli_ and _S. marcescens_ assemblies that were mislabelled as _Klebsiella_. Furthermore, you would miss some _Klebsiella_ assemblies that were mislabelled as something else, e.g. _Enterobacter_.

Bacsort is a collection of script to help you solve this problem. It assists in downloading assemblies and constructing a tree. Using this tree you can curate the species labels, resulting in a consistently named collection.


## Types of problems

RefSeq species mislabellings happen for a number of reasons, and Bacsort can help with each. Sometimes assemblies are just given completely wrong labels. These are usually quite easy to fix in Bacsort:

<p align="center"><img src="images/problem_1.png" alt="Problem 1" width="600"></p>

Sometimes a new species name has been coined for a group but hasn't yet fully caught on, leaving an awkward mix of the old and new names. Bacsort can bring that new name to the whole group:

<p align="center"><img src="images/problem_2.png" alt="Problem 2" width="600"></p>

Sometimes a group is not well studied and many species are not yet named. Bacsort can help to assign names to samples which are closely related to those already named:

<p align="center"><img src="images/problem_3.png" alt="Problem 3" width="600"></p>

And sometimes a group is very inconsistently named. These are the toughest problems to fix! The best solution is to find an authoritative piece of literature that defines the species in question so you can apply that scheme to your assemblies:

<p align="center"><img src="images/problem_4.png" alt="Problem 4" width="600"></p>



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




## Bacsorting assemblies

What follows are instructions for refining the species labels for one or more genera of interest. The end result will be assemblies organised into genus/species directories. For example:
```
Moraxella/
  atlantae/
    GCF_001591265.fna.gz
    GCF_001678995.fna.gz
    GCF_001679065.fna.gz
  boevrei/
    GCF_000379845.fna.gz
Psychrobacter/
  faecalis/
    GCF_001652315.fna.gz
    GCF_002836335.fna.gz
  glacincola/
    GCF_001411745.fna.gz
```

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

This script will make a new directory, `clusters` which contains assemblies with redundancy removed. For example, if there are 10 very similar assemblies, they will form one cluster and have only a single representative in `clusters`. Cluster representatives are chosen based on assembly N50 so more completed assemblies are preferred.

This step also produces a file, `cluster_accessions`, which lists the cluster name, followed by a tab, followed by a comma-delimited list of the assemblies in that cluster, with the representative assembly marked with a `*`:
```
Edwardsiella_01	GCF_000020865.1.fna.gz,GCF_000146305.1.fna.gz,GCF_000711155.1.fna.gz,GCF_000804515.1.fna.gz,GCF_002075835.1.fna.gz*
Edwardsiella_02	GCF_000022885.2.fna.gz*,GCF_000264785.1.fna.gz
Edwardsiella_03	GCF_000163955.1.fna.gz*,GCF_000341505.1.fna.gz
Edwardsiella_04	GCF_000264765.2.fna.gz*,GCF_000711175.1.fna.gz,GCF_000800725.1.fna.gz,GCF_000804575.1.fna.gz,GCF_000804595.1.fna.gz,GCF_000804615.1.fna.gz,GCF_001186215.1.fna.gz
Edwardsiella_05	GCF_000264805.1.fna.gz*
Edwardsiella_06	GCF_000264825.1.fna.gz*
Edwardsiella_07	GCF_000348565.1.fna.gz*
Edwardsiella_08	GCF_000474215.1.fna.gz*
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
* FastANI is designed for 80-100% nucleotide identity and will therefore struggle with greater divergences (see [Combining Mash and FastANI distances](#combining-mash-and-fastani-distances) below).
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
pairwise_identities_to_distance_matrix.py --max_dist 0.2 tree/fastani_output > tree/fastani.phylip
```

#### Combining Mash and FastANI distances

If you would like, you can use this script to combine FastANI distances (which are very good up to 20% divergence) with Mash distances (which can handle greater divergence) to get a best-of-both-worlds distance matrix:

```
combine_distance_matrices.py tree/fastani.phylip tree/mash.phylip > tree/distances.phylip
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

It produces the tree in two formats: newick and PhyloXML. The PhyloXML tree is suitable for viewing in [Archaeopteryx](https://sites.google.com/site/cmzmasek/home/software/archaeopteryx). Clades which perfectly define a species (i.e. all instances of that species are in a clade which contains no other species) will be coloured in the tree. 

To illustrate, here an Archaeopteryx visualisation of the genus _Edwardsiella_, before any curation:

<p align="center"><img src="images/Edwardsiella_before.png" alt="Edwardsiella before" width="100%"></p>

Each tip on the tree is labelled with the cluster name (e.g. 'Edwardsiella_01') and the species contained in that cluster (e.g. '3 x Edwardsiella tarda, 2 x Edwardsiella piscicida'). The species _ictaluri_ and _hoshinae_ are coloured because they are in consistent clades, but the other species are not. _E. tarda_, for example, is mostly contained in the bottom clade, but cluster Edwardsiella_01 also contains three _E. tarda_ assemblies. Since a single clade cannot define that entire species, it is not consistent.

You can then focus on the uncoloured parts where mislabellings may occur. As you find species label errors, add them to the `species_definitions` file in this format:
```
GCF_000020865   Edwardsiella piscicida
GCF_000146305   Edwardsiella piscicida
GCF_000804515   Edwardsiella piscicida
GCF_001186215   Edwardsiella anguillarum
GCF_000711175   Edwardsiella anguillarum
GCF_000800725   Edwardsiella anguillarum
```

You may have to refer back to the `cluster_accessions` file to see exactly which assemblies are in a particular cluster and need to be renamed.

You can then run `find_species_clades.py` again to generate a new tree with your updated definitions. This process can be repeated (fix labels, make tree, fix labels, make tree, etc) until you have no more changes to make.

Here is the _Edwardsiella_ tree after curation, now with all consistent species:
<p align="center"><img src="images/Edwardsiella_after.png" alt="Edwardsiella after" width="80%"></p>


### Step 6: copy assemblies and/or clusters to species directories

When you are happy with the species labels in the tree, you can run this command to copy assemblies into species directories:
```
copy_assemblies.py
```

It will make a directory titled `assemblies_binned` with subdirectories for each genus and then subdirectories for each species.

Alternatively, you can create a `clusters_binned` directory:
```
copy_clusters.py
```

This is the same as above, except that the redundancy-removed clusters are used instead of all assemblies. The `clusters_binned` directory will therefore be smaller than the `assemblies_binned` directory.



## Classify new assemblies using Mash

Mash can use your newly organised assemblies to query unknown assemblies and give the best match.

First, build a sketch of all your organised assemblies. You can do this either with `assemblies_binned` (all assemblies) or `clusters_binned` (redundancy-removed), but I'd recommend the latter for performance:
```
cd clusters_binned
mash sketch -o sketches -s 100000 -p 4 */*/*.fna.gz
```

Now you can use Mash directly to find the distances between your query and all of your reference assemblies:
```
mash dist -s 100000 -p 4 sketches.msh query.fasta | sort -gk3,3
```

Or you can use this script with comes with Bacsort:
```
classify_assembly_using_mash.py sketches.msh query.fasta
```

This script has some additional logic to help with classification:
* The `--threshold` option controls how close the query must be to a reference to count as a match (default: 5%)
* The `--contamination_threshold` option helps to spot contaminated assemblies. If the top two genera have matches closer than this (default: 2%), the assembly is considered contaminated. E.g. if your assembly is a strong match to both _Klebsiella_ and _Citrobacter_, then something is probably not right!



## Using Bacsort with Centrifuge

You can use Bacsort's assemblies to build a [Centrifuge](http://www.ccb.jhu.edu/software/centrifuge/) database which incorporates your revised species labels. For these instructions, I assume that you've only run Bacsort on _some_ bacterial genera, but you'd still like to include _all_ bacterial genera in the Centrifuge database. I also assume that you've already run all of the Bacsort steps above and produced a `clusters_binned` directory.

To begin, use Centrifuge to download genomes and taxonomy information (more detail is available in the [Centrifuge docs](http://www.ccb.jhu.edu/software/centrifuge/manual.shtml#database-download-and-index-building)):
```
mkdir -p centrifuge_bacsort
cd centrifuge_bacsort
centrifuge-download -o taxonomy taxonomy
centrifuge-download -o library -m -d bacteria refseq > seqid2taxid.map
```

Now from the same directory, run this script to incorporate your Bacsorted genomes into the Centrifuge library:
```
prepare_centrifuge_library.py /path/to/clusters_binned .
```

This script does two things:
1. It copies your Bacsorted genomes into the Centrifuge library and makes appropriate entries in the `seqid2taxid.map` file.
2. It _removes_ any genomes already in the Centrifuge library which are in the same genera as your Bacsorted assemblies. This prevent conflicts between the existing species labels and your refined species labels.

Now you can continue building the Centrifuge library (again, head to the [Centrifuge docs](http://www.ccb.jhu.edu/software/centrifuge/manual.shtml#database-download-and-index-building) for more detail):
```
cat library/*/*.fna > input-sequences.fna
centrifuge-build -p 16 --conversion-table seqid2taxid.map --taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp input-sequences.fna bacsort
```

This should create `*.cf` files which comprise your new Centrifuge database - use it like you would any other Centrifuge database! You can now delete any other files made along the way to save disk space.



## Using Bacsort with Kraken

Kraken is a similar tool to Centrifuge, and you can similarly build a database using Bacsorted assemblies. The assumptions stated above in the [Using Bacsort with Centrifuge](#using-bacsort-with-centrifuge) section apply here too.

First create a bacterial database as described in the [Kraken docs](http://ccb.jhu.edu/software/kraken/MANUAL.html#custom-databases):
```
kraken-build --download-taxonomy --db bacsort
kraken-build --download-library bacteria --db bacsort
```

Now run this script which does two things: 1) adjusts the taxonomy IDs in the Kraken database to match Bacsort's classifications, and 2) prepares additional Bacsorted assemblies for inclusion in the database:
```
prepare_kraken_library.py /path/to/clusters_binned bacsort
```

Now we can add the remaining assemblies to the library:
```
for f in additional_assemblies/*.fna; do
    kraken-build --add-to-library $f --db bacsort
done
```

Finally, we build the library (again, please read the [Kraken docs](http://ccb.jhu.edu/software/kraken/MANUAL.html#custom-databases) for more info):
```
kraken-build --build --threads 16 --db bacsort
```



## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
