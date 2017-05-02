# 3_haplotype_converter

Usage: python 3_haplotype_converter.py [full path to filtered haplotypes file] [full path to an output directory]

Example:

`python 3_haplotype_converter.py Volume/Stacks_Outputs/batch_1.singleRandomSNP_haplotypes.tsv Volume/Stacks_Outputs/Input_Files`

Designed to work with output of '2_haplo_summary.py' script, particularly the:

	'batch_1.singleFirstSNP_haplotypes.tsv' 

or
      
	'batch_1.singleRandomSNP_haplotypes.tsv'

output files. 
    
Assumes naming scheme is batch_1.SOMETHING_haplotypes.tsv,
and will rename files with the SOMETHING_haplotypes portion
of the filename. That is, there must be two '.' in the filename
or this won't name the outputs properly.

Opens the filtered haplotype.tsv file, filters out samples above a max missing
data threshold selected by user, and converts it into the following input formats:

Structure

Adegenet (two types)

Phylip

Fasta

Nexus

Plink (.bed and .map) - mainly to use for EEMS


Written for Python 2.7.3

External Dependencies: Numpy (Numerical Python)

# Dan Portik

daniel.portik@uta.edu

February 2016
