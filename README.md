Stacks Pipeline
---------------

### Overview

The goal of this workflow is to automate all major steps involved with processing common ddRADseq datasets, using the newest version of [**Stacks**](http://catchenlab.life.illinois.edu/stacks/) (v2.4). In particular, this workflow is designed to process single end (SE) read data generated from ddRADseq libraries prepared with the *SbfI* and *MspI* restriction enzymes. 

The gzipped fastq files from the sequencer are optionally trimmed for UMI sites, demultiplexed using `process_radtags`, and RAD cutsites are trimmed from all resulting fastq files. The partial or full Stacks pipeline can then be run (involving the `ustacks`, `cstacks`, `sstacks`, `tsv2bam`, `gstacks`, and `populations` modules), and the user can specify several key parameters for these steps. A range of missing data values is automatically used for the `populations` module, resulting in several `populations.haplotypes.tsv` files. These `haplotypes.tsv` files resulting from the `populations` module are then subjected to additional filtering. This includes removing samples exceeding a user-specified threshold of missing data, optionally removing singletons, and selecting one SNP site per locus (either the first site or a random site). This step can be run multiple times using different per-sample missing data thresholds, and distinct filtered tsv files are created for each run. Finally, all the filtered tsv files are converted into corresponding phylip, fasta, nexus, structure, ped, and map files. Summaries of all datasets created from the filtered tsv files are provided, allowing the user to choose which settings resulted in the highest quality dataset (in terms of number of samples, loci, and missing data). 

### Dependencies

The Stacks Pipeline relies on **seqtk** (**NEW as of v2.1**) and **Stacks v2.4**. These programs must be installed in path. They can be downloaded from the following sources:
+ [**seqtk**](https://github.com/lh3/seqtk)
+ [**Stacks**](http://catchenlab.life.illinois.edu/stacks/)

The Stacks Pipeline scripts can be run using Mac OSX (10.10+) and Linux, and can also work with Windows using a program like Cygwin. 

### Instructions

Documentation and usage instructions are available on the [**wiki page here**](https://github.com/dportik/Stacks_pipeline/wiki/Stacks-Pipeline-Instructions). 

The general order of the workflow is as follows:

1. **`Demultiplex_Trim.py`**: Demultiplexes fastq.gz files using `process_radtags` and trims RAD cutsites with `seqtk`. Offers an option to remove UMI sites of any length prior to demultiplexing.

2. **`Run_Stacks.py`**: Automates the full Stacks pipeline (`ustacks`, `cstacks`, `sstacks`, `tsv2bam`, `gstacks`, `populations`) or a partial Stacks run (post-ustacks or individual modules), based on a variety of user-selected options and parameter settings.

3. **`Filter_All_tsv.py`**: Applies filtering to loci contained in `populations.haplotypes.tsv` files resulting from independent runs of the `populations` module. Selects on SNP per locus (first or random site) and optionally removes singletons. Calculates per-sample missing data and removes samples above user-selected thresholds. Writes filtered tsv files.

4. **`Convert_All_tsv.py`**: Converts filtered tsv files to phylip, fasta, nexus, structure, ped, and map formats. Summarizes dataset metrics for convenient comparisons.

5. **`Convert_Stacks_Fasta_to_Loci`**: **Optional**. Parses the phased sequences of all samples and loci in a given `populations.samples.fa` fasta file and writes to locus-specific fasta files. Offers option to write both alleles per sample, the first allele, a random allele, or a consensus sequence. An optional filter option is included that only writes loci with at least one variable site. 

### Version

The current release of the Stacks Pipeline is [**v2.1**](https://github.com/dportik/Stacks_pipeline/releases). 

#### Major changes in v2.1:
  - `seqtk` is now used in place of `fastx_trimmer`. It is much faster and easier to install.

Changes in v2.0:
  - Now uses Stacks v2.41 (vs. 1.35).
  - All modules are now compatible with Python 2.7 and Python 3.7.
  - Offers new custom filtering and output file options.
  - Allows specification of key parameters for individual Stacks modules (including `-M`, `-m`, and `-n`). 

### License

GNU Lesser General Public License v3.0

### Contact

The Stacks Pipeline is maintained by Daniel Portik (daniel.portik@gmail.com)
