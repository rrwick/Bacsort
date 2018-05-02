# Bacsort paper scripts

This directory contains scripts and information regarding the data in the Bacsort manuscript:
[TBA link to manuscript](URL)

Overall, these scripts are _not_ meant to be one-size-fits-all. They should illustrate the methods used, but if you want to run them yourself, you'll probably need to look inside them and modify them to suit your system. In particular, the classification scripts are written for a Slurm-managed cluster, so if you don't use Slurm, you'll need to change those. That being said, I tried to make the scripts somewhat generalisable by using variables instead of hard-coded paths. So please set these before running them:
* `BACSORT`: path to the Bacsort repo directory
* `BACSORT_GENOMES`: path to where the Bacsorting will happen
* `TEST_ASSEMBLIES`: path to a directory of the test assemblies (`*.fasta`)
* `TEST_READS`: path to a directory of the test reads (two per sample: `*_1.fastq.gz` and `*_2.fastq.gz`)
* `DATABASES`: path to where the Kraken and Centrifuge databases will be built
* `DATABASE_SUFFIX`: a string appended to Kraken/Centrifuge database names (I used 2018-04-30)
* `CLASSIFICATION`: path to where classification results are stored
* `SLURM_PARTITION`: name of your Slurm partition (for `-p`)

Also, you unfortunately can't run these and get the exact same results as in the manuscript, because they involve downloading publicly available genomes. Therefore, if you run the scripts, you'll get more genomes than I did when I ran them on 30 April 2018.



## Bacsorting scripts

These are the commands I used for doing the actual Bacsorting of the genomes:

* [`bacsort_enterobacterales.sh`](bacsort_enterobacterales.sh): Enterobacterales order
* [`bacsort_moraxellaceae.sh`](bacsort_moraxellaceae.sh): Moraxellaceae family



## Classification scripts:

These are the scripts for classifying the assemblies used in the manuscript.

* Building classification databases:
  * [`build_kraken_databases.sh`](build_kraken_databases.sh): Commands to construct two [Kraken](http://ccb.jhu.edu/software/kraken/) databases, one using the RefSeq complete bacterial/archaeal genomes (as instructed in the [Kraken docs](http://ccb.jhu.edu/software/kraken/MANUAL.html#custom-databases)) and another where I modified the database using my Bacsorted assemblies.
  * [`build_centrifuge_databases.sh`](build_centrifuge_databases.sh): Same as above, but for [Centrifuge](http://www.ccb.jhu.edu/software/centrifuge/).
* Classifying genomes:
  * [`classify_reads_with_kraken.sh`](classify_reads_with_kraken.sh): Run Kraken on the test reads, using both the RefSeq and the Bacsort databases.
  * [`classify_reads_with_centrifuge.sh`](classify_reads_with_centrifuge.sh): Same as above, but for Centrifuge.
  * [`classify_assemblies_with_mash.sh`](classify_assemblies_with_mash.sh): Uses Mash find the Bacsorted assembly which best matches each test genome assembly.
  * [`classify_reads_with_mash.sh`](classify_reads_with_mash.sh): Same as above, but uses test genome reads instead of assemblies.



## Other scripts

* [`get_top_kraken_species.py`](get_top_kraken_species.py): Extract the top species match from a Kraken-style report. See the script's header for an example that illustrates why it's needed. Used by the `classify_reads_with_kraken.sh` and `classify_reads_with_centrifuge.sh` scripts described above.
* [`find_name_changes.py`](find_name_changes.py): Used to generate the RefSeq renamings table.
