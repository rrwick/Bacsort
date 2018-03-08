# Bacterial species Dejumbler


## Requirements

Running Dejumbler requires that you have [Mash](https://github.com/marbl/Mash) and [Quicktree](https://github.com/khowe/quicktree) installed and available in your PATH. If you can type `mash -h` and `quicktree -h` into your terminal and not get an error, you should be good!

You'll also need Python3 and [BioPython](http://biopython.org/). If `python3 -c "import Bio"` doesn't give you an error, you should be good! If you need to install BioPython, it's easiest to do with pip: `pip3 install biopython`




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
dejumbler_01_download_genomes.sh "Citrobacter Klebsiella Salmonella Yersinia"
```

### Step 2: cluster assemblies

```
dejumbler_02_cluster_genera.py "Citrobacter Klebsiella Salmonella Yersinia"
```

### Step 3: distance matrix for all clusters

There are two alternative ways to accomplish this step: using Mash and using FastANI.

#### Step 3a: Mash distance matrix

This is the simpler (and probably faster) of the two options. Mash estimates pairwise genomic distance by comparing the set of k-mers in the two assemblies. This means that both core sequence (shared by both assemblies) and accessory sequence (only in one assembly) affect the result.

```
dejumbler_03a_mash_distances.sh
```

#### Step 3a: FastANI distance matrix

FastANI, while fast compared to other methods of getting average nucleotide identity, is not as fast as Mash. Unless you are working with a small dataset, it will probably be necessary to parallelise this process in a method appropriate to your system.

```
dejumbler_03b1_fastani_distances.sh
```

Once the distances are computed, they must be converted into a PHYLIP distance matrix, which is relatively quick and carried out using this command:

```
dejumbler_03b2_fastani_distances.py
```

### Step 4: build tree

```
dejumbler_04_build_tree.sh
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
