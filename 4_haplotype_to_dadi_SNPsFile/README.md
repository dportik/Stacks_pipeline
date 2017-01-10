# 4_haplotype_to_dadi_SNPsFile

Usage: python 4_haplo_to_dadiSNPsV2.py [full path to haplotypes file] [full path to an output directory]

Takes a single SNP per locus haplotype.tsv file and prepares for SNP text file input to dadi.

Will make SNPs input file for one population, two populations, three populations, or 
four populations, based on lists created by user and functions called.
The major and minor alleles will be assigned based on all the samples in populations that
are included for a particular output file (and not based on the entire input file). If 
sites become constant by narrowing down individuals for a set of populations, they are 
ignored. 


Designed to work with output of '1_haplo_summary.py' script, particularly the:

    'batch_1.singleFirstSNP_haplotypes.tsv' 

or 
    'batch_1.singleRandomSNP_haplotypes.tsv'
    
output files. 
		
This script requires hand editing at the #*********** indicators.


Written for Python 2.7.3

External Dependencies: Numpy (Numerical Python)

# Dan Portik

daniel.portik@uta.edu

August 2016
