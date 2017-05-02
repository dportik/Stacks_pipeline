# 5_stacks_seq_depth

Usage: python 5_stacks_seq_depth.py [full path to a fully filtered 'batch_#.haplotypes.tsv' file] [directory with all "SAMPLE.matches.tsv" files]

Example:

`python 5_stacks_seq_depth.py Volume/Stacks_Outputs/population_r70/batch_1.singleRandomSNP_haplotypes.tsv Volume/Stacks_Outputs/3_Trimmed_Output`


Designed to work with output of '2_haplo_summary.py' script, particularly the:

    'batch_1.singleFirstSNP_haplotypes.tsv' 

or 

    'batch_1.singleRandomSNP_haplotypes.tsv'
    
This script will calculate the sequencing depth of each locus across all samples,
as well as the sequencing depth of each sample. These values will come from the
SAMPLE.matches.tsv files created for all samples created during the sstacks step.
Because these files contain ALL matched loci to the catalog, the loci included
will come from a filtered haplotypes.tsv file, which presumably will contain your
final set of loci to be used in analyses. The loci contained in this file will
be used to search the SAMPLE.matches.tsv files, so as to only include the
relevant loci in the sequencing depth calculations. Note that because the 
SAMPLE.matches.tsv files contain depth per allele, rather than locus, this
script will add depth values across alleles of the same locus for the
calculation, and if a locus is missing it will receive a value of 0. This can be
verified in the output, which reports how many loci were inspected per sample.
It should be the exact number of loci in the haplotypes.tsv file.

Unfortunately to process files this way increases the time it takes to run, so a
time stamp is shown at the beginning of each sample to provide an estimate. For
a data set containing ~7,500 loci, this script took ~2-4 minutes per sample. 


The outputs take two forms, and include the average, minimum, maximum, median
and standard deviation of the sequencing depths.


Sequencing depth stats on a per sample basis. Based on 1) only loci present, or 
2) including missing loci as depth value of zero.

1. Output_Sample_NoMissingData_Coverages.txt

2. Output_Sample_MissingData_Coverages.txt

example output:

```
Sample	Avg	Min	Max	Median	SD
1_AM224	42.0	1.0	537.0	16.0	54.7
1_AM225	52.7	1.0	558.0	20.0	70.4
1_AM226	54.2	1.0	492.0	20.0	73.7
1_AM227	46.4	1.0	526.0	17.0	62.3
```


Sequencing depth stats on a per locus basis. Based on 3) only when the locus is
present in the samples, or 4) when the locus is missing in a sample it is included 
as depth value of zero.

3. Output_Loci_NoMissingData_Coverages.txt

4. Output_Loci_MissingData_Coverages.txt

example output:

```
Locus	Avg	Min	Max	Median	SD
1	31.2	3.0	135.0	23.0	26.5
5	17.6	4.0	47.0	15.0	10.3
6	179.8	7.0	673.0	166.0	137.7
9	146.8	6.0	745.0	128.0	109.1
11	24.5	5.0	76.0	20.0	15.4
12	14.0	3.0	38.0	12.0	8.0
```


These ouptut files can be used to make plots in R to visualize the data. 
    
    
Written for Python 2.7

External Dependencies: Numpy (Numerical Python)

# Dan Portik

daniel.portik@uta.edu

April 2017
