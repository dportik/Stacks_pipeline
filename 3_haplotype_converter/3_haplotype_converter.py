import sys
import os
import shutil
import datetime
import numpy as np
'''
Usage: python 3_haplotype_converter.py [full path to filtered haplotypes file] [full path to an output directory]

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

############################################
Written for Python 2.7.3
External Dependencies: Numpy (Numerical Python)
############################################

Dan Portik
daniel.portik@uta.edu
February 2016
'''

#get haplotype file path
haplo_file = sys.argv[1]
fh_haplo = open(haplo_file, 'r')

#get output dir
out_dir = sys.argv[2]
os.chdir(out_dir)

#break up path name to get appropriate parts of haplotype file name
file_name = haplo_file.split('/')
file_name = file_name[-1]
file_names = file_name.split('.')
batch = file_names[0]
haplo_name = file_names[1]

#load file lines into memory for downstream functions
haplo_lines = fh_haplo.readlines()

#obtain missing data threshold for individuals
print '\n','\n'
decision_perc = None
while decision_perc is None:
    try:
        decision_perc = int(raw_input("Please enter the maximum percentage (%) of missing data allowed for a sample (ex. 30): "))
    except ValueError:
        print "That wasn't a number."

#=================================================================================================
#simple function to calculate percents
def percent_calc(x,y):
    percentage = float( ( float(x) / float(y)) * float(100))
    percentage = np.around(percentage, decimals = 1)
    return percentage

#=================================================================================================
#function for examining missing data across individuals
def missing_data(lines,x):
    print '\n', "========================================================"
    print "Examining missing data for each sample..."

    passed = []
    not_passed = []
    
    #get number of columns to search through appropriate indices per line list slice
    columns = lines[1].strip()
    columns = columns.split('\t')
    
    for i in range(len(columns)):
        #skip first two indices, or columns (catalogID and counts)
        if i < 2:
            pass
        #for every column onwards
        elif i >= 2:
            ind_list = []
            missing_count = int(0)
            #identify the individual (column header)
            header = lines[0].strip()
            header = header.split('\t')
            #add individual name to temporary individual list
            ind_list.append(header[i])
            locus_count = int(0)

            #go through every line for this column (matching by index on the line split list)
            for line in lines[1:]:
                #count each locus for total
                locus_count+=1
                line = line.strip()
                line_list = line.split('\t')
                #add locus to list of this individual
                ind_list.append(line_list[i])
                #look for the missing data and keep count
                if line_list[i] == '-':
                    missing_count+=1

            #calculate percentage missing data using function above
            percent = percent_calc(missing_count, locus_count)

            #make not of data with >50% missing data on screen, write all info to screen
            if percent > float(x):
                not_passed.append(ind_list)
                print "***{0} is missing data for {1} loci ({2}% missing data)".format(header[i],missing_count,percent)
            else:
                print '\t', "{0} is missing data for {1} loci ({2}% missing data)".format(header[i],missing_count,percent)
                passed.append(ind_list)

                
    print '\n', "Number of individuals passed = {}".format(len(passed))
    print "Number of individuals not passed = {}".format(len(not_passed))
    return passed, not_passed
#passed is a list, element[0] is sample name and elements[1:] are loci (A, A/T, -, etc.)


#=================================================================================================
#function for getting names of loci
def loci(lines):
    loci_list = []
    for line in lines[1:]:
        line = line.strip()
        line = line.split('\t')
        loci_list.append(line[0])
    return loci_list

#=================================================================================================
#create structure file
def structure(passed,prefix,loci):
    print '\n', "========================================================"
    print "Writing STRUCTURE input file..."

    structure_out = prefix+"_structure_classic.str"
    fh_struct = open(structure_out, 'a')
    fh_struct.write('\t')
    for locus in loci:
        fh_struct.write(locus+'\t')
    fh_struct.write('\n')
    
    file_writing_list = []
    
    #for every individual list
    for ind_list in passed:
        #make two lists:
        sample_name = ind_list[0]
        
        sample_list_a = []
        a_name = sample_name+"_a"
        sample_list_a.append(a_name)
        
        sample_list_b = []
        b_name = sample_name+"_b"
        sample_list_b.append(b_name)

        #sort through SNPs now
        for snp in ind_list[1:]:
            #if it's not a heterozygous site
            if '/' not in snp:
                #treat gap as -9 notation for structure
                if snp == "-":
                    sample_list_a.append("-9")
                    sample_list_b.append("-9")

                else:
                    if snp == "A":
                        sample_list_a.append("1")
                        sample_list_b.append("1")
                    elif snp == "T":
                        sample_list_a.append("2")
                        sample_list_b.append("2")
                    elif snp == "C":
                        sample_list_a.append("3")
                        sample_list_b.append("3")
                    elif snp == "G":
                        sample_list_a.append("4")
                        sample_list_b.append("4")

            elif '/' in snp:
                snp = snp.split('/')
                if snp[0] == "A":
                    sample_list_a.append("1")
                elif snp[0] == "T":
                    sample_list_a.append("2")
                elif snp[0] == "C":
                    sample_list_a.append("3")
                elif snp[0] == "G":
                    sample_list_a.append("4")

                if snp[1] == "A":
                    sample_list_b.append("1")
                elif snp[1] == "T":
                    sample_list_b.append("2")
                elif snp[1] == "C":
                    sample_list_b.append("3")
                elif snp[1] == "G":
                    sample_list_b.append("4")

        file_writing_list.append(sample_list_a)
        file_writing_list.append(sample_list_b)

    for ind in file_writing_list:
        print '\t', "Writing sample {}".format(ind[0])
        for info in ind:
            fh_struct.write(info+'\t')
        fh_struct.write('\n')
    fh_struct.close()
    print '\n', "Wrote '{0}' to output directory.".format(structure_out), '\n'

#=================================================================================================
#create adegenet file
def adegenet_1(passed,prefix):
    print '\n', "========================================================"
    print "Writing Adegenet input file 1..."

    ad_out = prefix+"_adegenet_1.str"
    fh_ad = open(ad_out, 'a')
    
    file_writing_list = []
    
    #for every individual list
    for ind_list in passed:
        #make two lists:
        sample_name = ind_list[0]
        
        sample_list_a = []
        sample_list_a.append(sample_name)
        
        sample_list_b = []
        sample_list_b.append(sample_name)

        #sort through SNPs now
        for snp in ind_list[1:]:
            #if it's not a heterozygous site
            if '/' not in snp:
                #treat gap as -9 notation for structure
                if snp == "-":
                    sample_list_a.append("-9")
                    sample_list_b.append("-9")

                else:
                    if snp == "A":
                        sample_list_a.append("1")
                        sample_list_b.append("1")
                    elif snp == "T":
                        sample_list_a.append("2")
                        sample_list_b.append("2")
                    elif snp == "C":
                        sample_list_a.append("3")
                        sample_list_b.append("3")
                    elif snp == "G":
                        sample_list_a.append("4")
                        sample_list_b.append("4")

            elif '/' in snp:
                snp = snp.split('/')
                if snp[0] == "A":
                    sample_list_a.append("1")
                elif snp[0] == "T":
                    sample_list_a.append("2")
                elif snp[0] == "C":
                    sample_list_a.append("3")
                elif snp[0] == "G":
                    sample_list_a.append("4")

                if snp[1] == "A":
                    sample_list_b.append("1")
                elif snp[1] == "T":
                    sample_list_b.append("2")
                elif snp[1] == "C":
                    sample_list_b.append("3")
                elif snp[1] == "G":
                    sample_list_b.append("4")

        file_writing_list.append(sample_list_a)
        file_writing_list.append(sample_list_b)

    for ind in file_writing_list:
        print '\t', "Writing sample {}".format(ind[0])
        for info in ind:
            fh_ad.write(info+'\t')
        fh_ad.write('\n')
    fh_ad.close()
    print '\n', "Wrote '{0}' to output directory.".format(ad_out), '\n'

#=================================================================================================
#create adegenet file
def adegenet_2(passed,prefix,loci):
    print '\n', "========================================================"
    print "Writing Adegenet input file 2..."

    ad_out = prefix+"_adegenet_2.in"
    fh_ad = open(ad_out, 'a')
    fh_ad.write('\t')
    for locus in loci:
        fh_ad.write(locus+'\t')
    fh_ad.write('\n')

    
    file_writing_list = []
    
    #for every individual list
    for ind_list in passed:
        #make two lists:
        sample_name = ind_list[0]
        
        sample_list = []
        sample_list.append(sample_name)

        #sort through SNPs now
        for snp in ind_list[1:]:
            #if it's not a heterozygous site
            if '/' not in snp:
                #treat gap as -9 notation for structure
                if snp == "-":
                    sample_list.append("NA,NA")

                else:
                    if snp == "A":
                        sample_list.append("1,1")
                    elif snp == "T":
                        sample_list.append("2,2")
                    elif snp == "C":
                        sample_list.append("3,3")
                    elif snp == "G":
                        sample_list.append("4,4")

            elif '/' in snp:
                snp = snp.split('/')
                if snp[0] == "A" and snp[1] == "T":
                    sample_list.append("1,2")
                elif snp[0] == "A" and snp[1] == "C":
                    sample_list.append("1,3")
                elif snp[0] == "A" and snp[1] == "G":
                    sample_list.append("1,4")
                    
                elif snp[0] == "T" and snp[1] == "A":
                    sample_list.append("2,1")
                elif snp[0] == "T" and snp[1] == "C":
                    sample_list.append("2,3")
                elif snp[0] == "T" and snp[1] == "G":
                    sample_list.append("2,4")
                    
                elif snp[0] == "C" and snp[1] == "A":
                    sample_list.append("3,1")
                elif snp[0] == "C" and snp[1] == "T":
                    sample_list.append("3,2")
                elif snp[0] == "C" and snp[1] == "G":
                    sample_list.append("3,4")
                    
                elif snp[0] == "G" and snp[1] == "A":
                    sample_list.append("4,1")
                elif snp[0] == "G" and snp[1] == "T":
                    sample_list.append("4,2")
                elif snp[0] == "G" and snp[1] == "C":
                    sample_list.append("4,3")

        file_writing_list.append(sample_list)

    for ind in file_writing_list:
        print '\t', "Writing sample {}".format(ind[0])
        for info in ind:
            fh_ad.write(info+'\t')
        fh_ad.write('\n')
    fh_ad.close()
    print '\n', "Wrote '{0}' to output directory.".format(ad_out), '\n'

#=================================================================================================
#create phylip file
def phylip(passed,prefix):
    print '\n', "========================================================"
    print "Writing PHYLIP file..."

    phy_out = prefix+"_out.phy"
    fh_phy = open(phy_out, 'a')

    #count taxon number and base pairs
    taxon_no = len(passed)
    aln_len = ( len(passed[1]) - int(1))

    fh_phy.write("{0} {1}".format(taxon_no, aln_len)+'\n')

    for ind_list in passed:
        #write sample name followed by two whitespaces
        fh_phy.write(ind_list[0]+"  ")
        #loop through all the SNPs for this sample
        for snp in ind_list[1:]:
            #if it's not a heterozygous site
            if '/' not in snp:
                #treat gap as "N" for this alignment
                if snp == "-":
                    fh_phy.write("N")
                #otherwise will be an A, T, G, or C, so write as is
                else:
                    fh_phy.write(snp)

            #now to insert ambiguity codes for SNPs
            elif '/' in snp:
                snp = snp.split('/')
                if snp[0] == "A" and snp[1] == "T":
                    fh_phy.write("W")
                elif snp[0] == "A" and snp[1] == "C":
                    fh_phy.write("M")
                elif snp[0] == "A" and snp[1] == "G":
                    fh_phy.write("R")
                elif snp[0] == "C" and snp[1] == "A":
                    fh_phy.write("M")
                elif snp[0] == "C" and snp[1] == "G":
                    fh_phy.write("S")
                elif snp[0] == "C" and snp[1] == "T":
                    fh_phy.write("Y")
                elif snp[0] == "G" and snp[1] == "A":
                    fh_phy.write("R")
                elif snp[0] == "G" and snp[1] == "C":
                    fh_phy.write("S")
                elif snp[0] == "G" and snp[1] == "T":
                    fh_phy.write("K")
                elif snp[0] == "T" and snp[1] == "A":
                    fh_phy.write("W")
                elif snp[0] == "T" and snp[1] == "C":
                    fh_phy.write("Y")                    
                elif snp[0] == "T" and snp[1] == "G":
                    fh_phy.write("K")
                    
        #indicate where we are in the list
        print '\t', "Writing sample {}".format(ind_list[0])
        fh_phy.write('\n')
        
    fh_phy.close()
    print '\n', "Wrote '{0}' to output directory.".format(phy_out), '\n'



#=================================================================================================
#create nexus file
def nexus(passed,prefix):
    print '\n', "========================================================"
    print "Writing NEXUS file..."

    nex_out = prefix+"_out.nex"
    fh_nex = open(nex_out, 'a')

    #count taxon number and base pairs
    taxon_no = len(passed)
    aln_len = ( len(passed[1]) - int(1))

    fh_nex.write('''#NEXUS 

BEGIN DATA;
	DIMENSIONS  NTAX={0} NCHAR={1};
	FORMAT DATATYPE=DNA  MISSING=N GAP=- interleave;
MATRIX

'''.format(taxon_no, aln_len))

    for ind_list in passed:
        #write sample name followed by a tab
        fh_nex.write(ind_list[0]+'\t')
        #loop through all the SNPs for this sample
        for snp in ind_list[1:]:
            #if it's not a heterozygous site
            if '/' not in snp:
                #treat gap as "N" for this alignment
                if snp == "-":
                    fh_nex.write("N")
                #otherwise will be an A, T, G, or C, so write as is
                else:
                    fh_nex.write(snp)

            #now to insert ambiguity codes for SNPs
            elif '/' in snp:
                snp = snp.split('/')
                if snp[0] == "A" and snp[1] == "T":
                    fh_nex.write("W")
                elif snp[0] == "A" and snp[1] == "C":
                    fh_nex.write("M")
                elif snp[0] == "A" and snp[1] == "G":
                    fh_nex.write("R")
                elif snp[0] == "C" and snp[1] == "A":
                    fh_nex.write("M")
                elif snp[0] == "C" and snp[1] == "G":
                    fh_nex.write("S")
                elif snp[0] == "C" and snp[1] == "T":
                    fh_nex.write("Y")
                elif snp[0] == "G" and snp[1] == "A":
                    fh_nex.write("R")
                elif snp[0] == "G" and snp[1] == "C":
                    fh_nex.write("S")
                elif snp[0] == "G" and snp[1] == "T":
                    fh_nex.write("K")
                elif snp[0] == "T" and snp[1] == "A":
                    fh_nex.write("W")
                elif snp[0] == "T" and snp[1] == "C":
                    fh_nex.write("Y")                    
                elif snp[0] == "T" and snp[1] == "G":
                    fh_nex.write("K")
                    
        #indicate where we are in the list
        print '\t', "Writing sample {}".format(ind_list[0])
        fh_nex.write('\n')
    fh_nex.write(";"+'\n'+"end;")
    fh_nex.close()
    print '\n', "Wrote '{0}' to output directory.".format(nex_out), '\n'


#=================================================================================================
#create FASTA file
def fasta(passed,prefix):
    print '\n', "========================================================"
    print "Writing FASTA file..."

    fas_out = prefix+"_out.fasta"
    fh_fas = open(fas_out, 'a')

    #count taxon number and base pairs
    taxon_no = len(passed)
    aln_len = ( len(passed[1]) - int(1))

    for ind_list in passed:
        #write sample name followed by a tab
        fh_fas.write(">"+ind_list[0]+'\n')
        #loop through all the SNPs for this sample
        for snp in ind_list[1:]:
            #if it's not a heterozygous site
            if '/' not in snp:
                #treat gap as "N" for this alignment
                if snp == "-":
                    fh_fas.write("N")
                #otherwise will be an A, T, G, or C, so write as is
                else:
                    fh_fas.write(snp)

            #now to insert ambiguity codes for SNPs
            elif '/' in snp:
                snp = snp.split('/')
                if snp[0] == "A" and snp[1] == "T":
                    fh_fas.write("W")
                elif snp[0] == "A" and snp[1] == "C":
                    fh_fas.write("M")
                elif snp[0] == "A" and snp[1] == "G":
                    fh_fas.write("R")
                elif snp[0] == "C" and snp[1] == "A":
                    fh_fas.write("M")
                elif snp[0] == "C" and snp[1] == "G":
                    fh_fas.write("S")
                elif snp[0] == "C" and snp[1] == "T":
                    fh_fas.write("Y")
                elif snp[0] == "G" and snp[1] == "A":
                    fh_fas.write("R")
                elif snp[0] == "G" and snp[1] == "C":
                    fh_fas.write("S")
                elif snp[0] == "G" and snp[1] == "T":
                    fh_fas.write("K")
                elif snp[0] == "T" and snp[1] == "A":
                    fh_fas.write("W")
                elif snp[0] == "T" and snp[1] == "C":
                    fh_fas.write("Y")                    
                elif snp[0] == "T" and snp[1] == "G":
                    fh_fas.write("K")
                    
        #indicate where we are in the list
        print '\t', "Writing sample {}".format(ind_list[0])
        fh_fas.write('\n')
    fh_fas.close()
    print '\n', "Wrote '{0}' to output directory.".format(fas_out), '\n'
#=================================================================================================
#create plink file
def plink(passed,prefix,loci):
    print '\n', "========================================================"
    print "Writing plink input file (.ped)..."

    ped_out = prefix+"_plink.ped"
    fh_ped = open(ped_out, 'a')
    
    file_writing_list = []
    
    #for every individual list
    for ind_list in passed:
        #make two lists:
        sample_name = ind_list[0]

        fh_ped.write("{0} {0} 0 0 0 0  ".format(sample_name))

        #sort through SNPs now
        for snp in ind_list[1:]:
            #if it's not a heterozygous site
            if '/' not in snp:
                #treat gap as -9 notation for structure
                if snp == "-":
                    fh_ped.write("0 0  ")

                else:
                    if snp == "A":
                        fh_ped.write("A A  ")
                    elif snp == "T":
                        fh_ped.write("T T  ")
                    elif snp == "C":
                        fh_ped.write("C C  ")
                    elif snp == "G":
                        fh_ped.write("G G  ")

            elif '/' in snp:
                snp = snp.split('/')
                if snp[0] == "A" and snp[1] == "T":
                    fh_ped.write("A T  ")
                elif snp[0] == "A" and snp[1] == "C":
                    fh_ped.write("A C  ")
                elif snp[0] == "A" and snp[1] == "G":
                    fh_ped.write("A G  ")
                    
                elif snp[0] == "T" and snp[1] == "A":
                    fh_ped.write("T A  ")
                elif snp[0] == "T" and snp[1] == "C":
                    fh_ped.write("T C  ")
                elif snp[0] == "T" and snp[1] == "G":
                    fh_ped.write("T G  ")
                    
                elif snp[0] == "C" and snp[1] == "A":
                    fh_ped.write("C A  ")
                elif snp[0] == "C" and snp[1] == "T":
                    fh_ped.write("C T  ")
                elif snp[0] == "C" and snp[1] == "G":
                    fh_ped.write("C G  ")
                    
                elif snp[0] == "G" and snp[1] == "A":
                    fh_ped.write("G A  ")
                elif snp[0] == "G" and snp[1] == "T":
                    fh_ped.write("G T  ")
                elif snp[0] == "G" and snp[1] == "C":
                    fh_ped.write("G C  ")
        print '\t', "Writing sample {}".format(ind_list[0])
        fh_ped.write('\n')
        
    fh_ped.close()
    print '\n', "Wrote '{0}' to output directory.".format(ped_out), '\n'
    
    print '\n', "========================================================"
    print "Writing plink input file (.map)..."
    map_out = prefix+"_plink.map"
    fh_map = open(map_out, 'a')
    for locus in loci:
        fh_map.write("0"+'\t'+locus+'\t'+"0"+'\t'+"100"+'\n')
    fh_map.close()
    print '\n', "Wrote '{0}' to output directory.".format(map_out), '\n'

#=================================================================================================

#filter samples with missing data function
passed,not_passed = missing_data(haplo_lines,decision_perc)
#get list of loci
loci_list = loci(haplo_lines)
#write structure file
structure(passed,haplo_name,loci_list)
#write adegenet_1 file
adegenet_1(passed,haplo_name)
#write adegenet_2 file
adegenet_2(passed,haplo_name,loci_list)
#write phylip file
phylip(passed,haplo_name)
#write nexus file
nexus(passed,haplo_name)
#write fasta file
fasta(passed,haplo_name)
#write plink files
plink(passed,haplo_name,loci_list)

print '\n',"========================================================"
print "Total number of loci = {}.".format(len(loci_list))
print "Total number of samples included = {}.".format(len(passed))
print "Total number of samples discarded = {}.".format(len(not_passed)), '\n', '\n'
