# 1_stacks_pipeline

Usage: python 1_stacks_pipeline.py [full path to directory containing all de-multiplexed trimmed fq files]

Example:

`python 1_stacks_pipeline.py /Volumes/Project1/3_Trimmed_Output/`

Example file naming scheme (must follow the'SAMPLE.trim.fq' to work!):
```
4_MB312.trim.fq
```

Will result in many intermediate files:
```
4_MB312.trim.snps.tsv
4_MB312.trim.tags.tsv
4_MB312.trim.alleles.tsv
4_MB312.trim.matches.tsv
```

Full pipeline:

Script will run ustacks across all SAMPLE.trim.fq files in directory, then cstacks with all samples,
then sstacks on a per sample basis. Finally, the population script is run taking all samples as one
population, outputs several folders with results from population script with varying -r values 
(0, 50, 60, 70, 80, 90, 100).

User decides whether to run a) full pipeline (above), b) ustacks portion only, c) cstacks sstacks and
populations portion only, or d) populations only. As long as you stick to the same directory, you can 
run ustacks first (decision b), then the remaining pipeline later (decision c), or re-run populations (d).

ustacks call:

	ustacks -t fastq -f SAMPLE.trim.fq -i {iterable} -o . -r -m 5 -M 2 -p {user decision}

cstacks call:

	cstacks -b 1 {-s all samples together} -p {user decision} -o .

sstacks call:

	sstacks -b 1 -c batch_1 -s SAMPLE.trim -p {user decision} -o .

populations call:

	populations -b 1 -P . -M {path to script-generated file} -r {0, 50, 60, 70, 80, 90, 100} -t {user decision} --min_maf 0.05

This script was written for and tested with Stacks version 1.35. Note you can edit the call for any of the above modules in the code below for a specific application, these settings are most appropriate for the data sets we have been processing. Familiarize yourself with the different thresholds and parameters to understand the predefined values and how they may affect your processed data.


Written for Python 2.7.3

External Dependencies: Stacks v1.35

ustacks, cstacks, sstacks, populations (called from command line)

# Dan Portik

daniel.portik@uta.edu

February 2016



If you decide to use these scripts or modify the code for your purposes, please cite:

*Portik, D.M., Leaché, A.D., Rivera, D., Blackburn, D.C., Rödel, M.-O., Barej, M.F., 
Hirschfeld, M., Burger, M., and M.K. Fujita. Evaluating mechanisms of diversification 
in a Guineo-Congolian forest frog using demographic model selection. 
In Review, Molecular Ecology.*
