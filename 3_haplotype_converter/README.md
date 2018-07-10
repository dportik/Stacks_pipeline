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


**Citation Information:**

***Using the pipeline.***
The scripts involved with this pipeline were originally published as part of the following work:

+ *Portik, D.M., Leache, A.D., Rivera, D., Blackburn, D.C., Rodel, M.-O., Barej, M.F., Hirschfeld, M., Burger, M., and M.K. Fujita. 2017. Evaluating mechanisms of diversification in a Guineo-Congolian forest frog using demographic model selection. Molecular Ecology, 26: 5245-5263. doi: 10.1111/mec.14266*

If you use or modify this script for your own purposes, please cite this publication.


**Contact:**

Daniel Portik, PhD

Postdoctoral Researcher

University of Arizona

daniel.portik@gmail.com
