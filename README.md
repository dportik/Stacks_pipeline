Stacks Pipeline
---------------

## Overview

Coming soon!

## Dependencies

The Stacks Pipeline is relies on fastx_trimmer and Stacks v2.4. These programs must be installed in path. They can be downloaded from the following sources:
+ [**fastx_trimmer**](http://hannonlab.cshl.edu/fastx_toolkit/download.html)
+ [**Stacks**](http://catchenlab.life.illinois.edu/stacks/)

The Stacks Pipeline scripts can be run using Mac OSX (10.10+) and Linux, and can also work with Windows using a program like Cygwin. 

## Instructions

Documentation and usage instructions are available on the wiki page [here](https://github.com/dportik/Stacks_pipeline/wiki). The general order of the workflow is as follows:

1. `Demultiplex_Trim.py`: Demultiplexes fastq.gz files using `process_radtags` and trims RAD cutsites with `fastx_trimmer`. Offers an option to remove UMI sites of any length prior to demultiplexing.
2. `Run_Stacks.py`: Automates the full Stacks pipeline (`ustacks`, `cstacks`, `sstacks`, `tsv2bam`, `gstacks`, `populations`) or a partial Stacks run (post-ustacks or individual modules), based on a variety of user-selected options and parameter settings.
3. `Filter_All_tsv.py`: Applies filtering to loci contained in `populations.haplotypes.tsv` files resulting from independent runs of the `populations` module. Calculates per-sample missing data and removes samples above user-selected thresholds. Writes filtered tsv files.
4. `Convert_All_tsv.py`: Converts filtered tsv files to phylip, fasta, nexus, structure, ped, and map formats. Summarizes dataset metrics.


## Version

The current release of the Stacks Pipeline is [**v2.0**](https://github.com/dportik/Stacks_pipeline/releases). 

### Major changes in v2.0:
  - Now uses Stacks v2.41 (vs. 1.35).
  - All modules are now compatible with Python 2.7 and Python 3.7.
  - Offers new filtering and output options.
  - Allows specification of key parameters for Stacks modules (including `-M`, `-m`, and `-n`). 

## License

GNU Lesser General Public License v3.0

## Contact

The Stacks Pipeline is maintained by Daniel Portik (daniel.portik@gmail.com)
