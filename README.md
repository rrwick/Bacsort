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

### Step 3: build tree

```
dejumbler_03_build_tree.sh
```

### Step 4: curate tree

```
dejumbler_04_find_species_clades.py
```

### Step 5: copy assemblies to species directories

```
dejumbler_05_copy_assemblies.py
```




## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
