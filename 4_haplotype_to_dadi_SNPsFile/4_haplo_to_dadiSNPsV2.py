import sys
import os
import shutil
import datetime
import numpy as np
'''
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
		
This script requires hand editing at the *********** indicators.

############################################
Written for Python 2.7.3
External Dependencies: Numpy (Numerical Python)
############################################

Dan Portik
daniel.portik@uta.edu
August 2016
'''

#get haplotype file path
haplo_file = sys.argv[1]
fh_haplo = open(haplo_file, 'r')

#get output dir
out_dir = sys.argv[2]
os.chdir(out_dir)

#load file lines into memory for downstream functions
haplo_lines = fh_haplo.readlines()

###########################################
#Define all configurations of populations for data set, up to four populations to be written
#You can create 2-d style input files from any number of population combinations here
#Name of lists here must match name below when calling dadi functions, but is arbitrary
#Not all individuals need to be included in population sets

#************** You'll want to edit these lists with your sample names and population names

Nigeria = ["8_MVZ253426.trim","8_MVZ253427.trim","8_MVZ253428.trim","8_MVZ253429.trim","8_MVZ253430.trim","8_MVZ253431.trim","8_MVZ253432.trim","8_MVZ253433.trim","9_MVZ253434.trim","9_MVZ253435.trim"]

CVL = ["1_AM224.trim","1_AM225.trim","1_AM226.trim","1_AM227.trim","1_AM228.trim","1_AM229.trim","1_AM230.trim","1_AM232.trim","2_AM233.trim","2_AM234.trim","2_DPL3059.trim","2_DPL3279.trim","2_MB292.trim","2_MB293.trim","3_BMNH2008_625.trim","3_MB295.trim","3_MB297.trim","3_MB304.trim","3_MB310.trim","4_MB311.trim","4_MB312.trim","4_MB313.trim","4_MB316.trim","4_MB318.trim","5_MB320.trim","5_MC11_038.trim","5_MC11_040.trim","5_MCZA_136769.trim","6_MCZA_136770.trim","6_MCZA_136840.trim","6_MCZA_136841.trim","6_MCZA_136842.trim","6_MCZA_136843.trim","6_MH0406.trim","6_MH0419.trim","7_MM036.trim","7_MM050.trim","7_MM051.trim","7_MVZ253421.trim","7_MVZ253422.trim","7_MVZ253423.trim","7_MVZ253424.trim","7_MVZ253425.trim"]

South = ["3_KU113_BLS13571.trim","4_MB314.trim","4_MB315.trim","5_MBUR3400.trim","5_MC11_170.trim","5_MC11_171.trim","6_MCZA_139706.trim","9_NCSM78914.trim","9_NCSM78915.trim","9_TCWC95629.trim","9_TCWC95630.trim","2_YPM6988.trim","2_YPM6989.trim"]

South1 = ["4_MB314.trim","4_MB315.trim","5_MC11_170.trim","5_MC11_171.trim","9_TCWC95629.trim","9_TCWC95630.trim"]

South2 = ["2_YPM6988.trim","2_YPM6989.trim","3_KU113_BLS13571.trim","6_MCZA_139706.trim","9_NCSM78914.trim","9_NCSM78915.trim","5_MBUR3400.trim"]


CVLN = ["1_AM232.trim","2_AM233.trim","2_AM234.trim","2_DPL3059.trim","2_DPL3279.trim","3_BMNH2008_625.trim","4_MB316.trim","4_MB318.trim","5_MB320.trim","5_MCZA_136769.trim","6_MCZA_136770.trim","6_MCZA_136840.trim","6_MCZA_136841.trim","6_MCZA_136842.trim","6_MCZA_136843.trim","7_MM036.trim","7_MM050.trim","7_MM051.trim","7_MVZ253421.trim","7_MVZ253422.trim","7_MVZ253423.trim","7_MVZ253424.trim","7_MVZ253425.trim"]

CVLS = ["2_MB292.trim","2_MB293.trim","3_MB295.trim","3_MB297.trim","3_MB304.trim","3_MB310.trim","4_MB311.trim","4_MB312.trim","4_MB313.trim","5_MC11_038.trim","5_MC11_040.trim","6_MH0406.trim","6_MH0419.trim"]

#=================================================================================================
#function for examining missing data across individuals
def reformat(lines):
    print '\n', "==============================================================================================================="
    print "Reformatting tsv file..."
    print "===============================================================================================================", '\n'

    all_ind=[]
    
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
            #identify the individual (column header)
            header = lines[0].strip()
            header = header.split('\t')
            #add individual name to temporary individual list
            ind_list.append(header[i])

            #go through every line for this column (matching by index on the line split list)
            for line in lines[1:]:
                line = line.strip()
                line_list = line.split('\t')
                #add locus to list of this individual
                ind_list.append(line_list[i])
            all_ind.append(ind_list)

    return all_ind
#passed is a list of lists, element[0] is sample name and elements[1:] are loci (A, A/T, -, etc.)


#=================================================================================================
#function for getting names of loci
def loci(lines):
    loci_list = []
    for line in lines[1:]:
        line = line.strip()
        line = line.split('\t')
        loci_list.append(line[0])
    return loci_list
#passed is a list of loci names, from the tsv file it is a locus number from cstacks

#=================================================================================================
#function for dadi SNPs input file, 1 population

def dadi1(all_ind, loci, pop1, header1):
    outfile = "dadi_1pop_{0}_snps.txt".format(header1)

    print '\n', "==============================================================================================================="
    print "Creating SNPs file {}...".format(outfile)
    print "===============================================================================================================", '\n'
    
    fh_out = open(outfile, 'a')
    fh_out.write("Ingroup\tOutgroup\tAllele1\t{0}\tAllele2\t{0}\tGene\tPosition\n".format(header1))
    
    allele_list = ["alleles"]
    columns = len(all_ind[0])
    
    for i in range(columns):        
        allele_set = set()
        alleles = []
        #skip first index
        if i < 1:
            pass
        #for every column onwards
        elif i >= 1:
            for row in all_ind:
                #only search for major and minor allele within populations chosen
                if (row[0] in pop1):
                    if '/' not in row[i]:
                        if row[i] == "-":
                            pass
                        else:
                            alleles.append(row[i])
                            #add again to account for homozygous site
                            alleles.append(row[i])
                            allele_set.add(row[i])                        
                    elif '/' in row[i]:
                        snp = row[i].split('/')
                        #add each heterozygous allele site
                        alleles.append(snp[0])
                        allele_set.add(snp[0])
                        alleles.append(snp[1])
                        allele_set.add(snp[1])

            #Turn bases from this SNP site into a list (length of 2 if diploid SNP site, 1 if became a constant site) 
            temp_list = list(allele_set)

            #sometimes narrowing down SNP file produces a constant site, deal with this here
            if len(temp_list) == int(2):
                base1 = alleles.count(temp_list[0])
                base2 = alleles.count(temp_list[1])
                
                major_minor = []
                if base1 > base2:
                    major_minor.append(temp_list[0])
                    major_minor.append(temp_list[1])
                else:
                    major_minor.append(temp_list[1])
                    major_minor.append(temp_list[0])
                allele_list.append(major_minor)

            #Put some kind of indicator in 1st index to remind it is a bad site and will exclude
            elif len(temp_list) == int(1):
                major_minor = []
                major_minor.append(temp_list[0])
                major_minor.append("Constant")
                allele_list.append(major_minor)
            elif len(temp_list) == int(0):
                major_minor = []
                major_minor.append("Missing")
                major_minor.append("Constant")
                allele_list.append(major_minor)
   
    for i in range(columns):
        #check if locus has the constant site flag or not
        if (i >= 1) and (allele_list[i][1] == "Constant"):
            print "Excluding constant site or missing locus: {}".format(loci[i])            
        if (i >= 1) and (allele_list[i][1] != "Constant"):            
            pop1_major = int(0)
            pop1_minor = int(0)
            pop2_major = int(0)
            pop2_minor = int(0)
            
            for row in all_ind:
                if row[0] in pop1:
                    if '/' in row[i]:
                        snp = row[i].split('/')
                        if snp[0] == allele_list[i][0]:
                            pop1_major += 1
                        elif snp[0] == allele_list[i][1]:
                            pop1_minor += 1
                        if snp[1] == allele_list[i][0]:
                            pop1_major += 1
                        elif snp[1] == allele_list[i][1]:
                            pop1_minor += 1
                    elif '-' in row[i]:
                        pass
                    elif row[i] == allele_list[i][0]:
                        pop1_major += 2
                    elif row[i] == allele_list[i][1]:
                        pop1_minor += 2
                else:
                	pass

            print "-{0}- --- {0} {1} {2} {3} {4} 15".format(allele_list[i][0],pop1_major,allele_list[i][1],pop1_minor,loci[i])
            fh_out.write("-{0}-\t---\t{0}\t{1}\t{2}\t{3}\t{4}\t15\n".format(allele_list[i][0],pop1_major,allele_list[i][1],pop1_minor,loci[i]))
    fh_out.close()


#=================================================================================================
#function for dadi SNPs input file, 2 pops

def dadi2(all_ind, loci, pop1, pop2, header1, header2):
    outfile = "dadi_2pops_{0}_{1}_snps.txt".format(header1, header2)

    print '\n', "==============================================================================================================="
    print "Creating SNPs file {}...".format(outfile)
    print "===============================================================================================================", '\n'
    
    fh_out = open(outfile, 'a')
    fh_out.write("Ingroup\tOutgroup\tAllele1\t{0}\t{1}\tAllele2\t{0}\t{1}\tGene\tPosition\n".format(header1, header2))
    
    allele_list = ["alleles"]
    columns = len(all_ind[0])
    
    for i in range(columns):        
        allele_set = set()
        alleles = []
        #skip first index
        if i < 1:
            pass
        #for every column onwards
        elif i >= 1:
            for row in all_ind:
                #only search for major and minor allele within populations chosen
                if (row[0] in pop1) or (row[0] in pop2):
                    if '/' not in row[i]:
                        if row[i] == "-":
                            pass
                        else:
                            alleles.append(row[i])
                            #add again to account for homozygous site
                            alleles.append(row[i])
                            allele_set.add(row[i])                        
                    elif '/' in row[i]:
                        snp = row[i].split('/')
                        #add each heterozygous allele site
                        alleles.append(snp[0])
                        allele_set.add(snp[0])
                        alleles.append(snp[1])
                        allele_set.add(snp[1])

            #Turn bases from this SNP site into a list (length of 2 if diploid SNP site, 1 if became a constant site) 
            temp_list = list(allele_set)

            #sometimes narrowing down SNP file produces a constant site, deal with this here
            if len(temp_list) == int(2):
                base1 = alleles.count(temp_list[0])
                base2 = alleles.count(temp_list[1])
                
                major_minor = []
                if base1 > base2:
                    major_minor.append(temp_list[0])
                    major_minor.append(temp_list[1])
                else:
                    major_minor.append(temp_list[1])
                    major_minor.append(temp_list[0])
                allele_list.append(major_minor)

            #Put some kind of indicator in 1st index to remind it is a bad site and will exclude
            elif len(temp_list) == int(1):
                major_minor = []
                major_minor.append(temp_list[0])
                major_minor.append("Constant")
                allele_list.append(major_minor)
                
            elif len(temp_list) == int(0):
                major_minor = []
                major_minor.append("Missing")
                major_minor.append("Constant")
                allele_list.append(major_minor)
                
   
    for i in range(columns):
        #check if locus has the constant site flag or not
        if (i >= 1) and (allele_list[i][1] == "Constant"):
            print "Excluding constant site locus: {}".format(loci[i])            
        if (i >= 1) and (allele_list[i][1] != "Constant"):            
            pop1_major = int(0)
            pop1_minor = int(0)
            pop2_major = int(0)
            pop2_minor = int(0)
            
            for row in all_ind:
                if row[0] in pop1:
                    if '/' in row[i]:
                        snp = row[i].split('/')
                        if snp[0] == allele_list[i][0]:
                            pop1_major += 1
                        elif snp[0] == allele_list[i][1]:
                            pop1_minor += 1
                        if snp[1] == allele_list[i][0]:
                            pop1_major += 1
                        elif snp[1] == allele_list[i][1]:
                            pop1_minor += 1
                    elif '-' in row[i]:
                        pass
                    elif row[i] == allele_list[i][0]:
                        pop1_major += 2
                    elif row[i] == allele_list[i][1]:
                        pop1_minor += 2
                elif row[0] in pop2:
                    if '/' in row[i]:
                        snp = row[i].split('/')
                        if snp[0] == allele_list[i][0]:
                            pop2_major += 1
                        elif snp[0] == allele_list[i][1]:
                            pop2_minor += 1
                        if snp[1] == allele_list[i][0]:
                            pop2_major += 1
                        elif snp[1] == allele_list[i][1]:
                            pop2_minor += 1
                    elif '-' in row[i]:
                        pass
                    elif row[i] == allele_list[i][0]:
                        pop2_major += 2
                    elif row[i] == allele_list[i][1]:
                        pop2_minor += 2
                else:
                	pass

            print "-{0}- --- {0} {1} {2} {3} {4} {5} {6} 15".format(allele_list[i][0],pop1_major,pop2_major,allele_list[i][1],pop1_minor,pop2_minor,loci[i])
            fh_out.write("-{0}-\t---\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t15\n".format(allele_list[i][0],pop1_major,pop2_major,allele_list[i][1],pop1_minor,pop2_minor,loci[i]))
    fh_out.close()


#=================================================================================================
#function for dadi SNPs input file, 3 pops

def dadi3(all_ind, loci, pop1, pop2, pop3, header1, header2, header3):
    outfile = "dadi_3pops_{0}_{1}_{2}_snps.txt".format(header1, header2, header3)
    
    print '\n', "==============================================================================================================="
    print "Creating SNPs file {}...".format(outfile)
    print "===============================================================================================================", '\n'
    
    fh_out = open(outfile, 'a')
    fh_out.write("Ingroup\tOutgroup\tAllele1\t{0}\t{1}\t{2}\tAllele2\t{0}\t{1}\t{2}\tGene\tPosition\n".format(header1, header2, header3))

    allele_list = ["alleles"]
    columns = len(all_ind[0])
    
    for i in range(columns):        
        allele_set = set()
        alleles = []
        if i < 1:
            pass
        elif i >= 1:
            for row in all_ind:
                if (row[0] in pop1) or (row[0] in pop2) or (row[0] in pop3):
                    if '/' not in row[i]:
                        if row[i] == "-":
                            pass
                        else:
                            alleles.append(row[i])
                            alleles.append(row[i])
                            allele_set.add(row[i])                        
                    elif '/' in row[i]:
                        snp = row[i].split('/')
                        alleles.append(snp[0])
                        allele_set.add(snp[0])
                        alleles.append(snp[1])
                        allele_set.add(snp[1])

            temp_list = list(allele_set)

            if len(temp_list) == int(2):
                base1 = alleles.count(temp_list[0])
                base2 = alleles.count(temp_list[1])
                
                major_minor = []
                
                if base1 > base2:
                    major_minor.append(temp_list[0])
                    major_minor.append(temp_list[1])
                else:
                    major_minor.append(temp_list[1])
                    major_minor.append(temp_list[0])
                allele_list.append(major_minor)
                
            elif len(temp_list) < int(2):
                major_minor = []
                major_minor.append(temp_list[0])
                major_minor.append("Constant")
                allele_list.append(major_minor)
    
    for i in range(columns):
        if (i >= 1) and (allele_list[i][1] == "Constant"):
            print "Excluding constant site locus: {}".format(loci[i])            
        if (i >= 1) and (allele_list[i][1] != "Constant"):            
            pop1_major = int(0)
            pop1_minor = int(0)
            pop2_major = int(0)
            pop2_minor = int(0)
            pop3_major = int(0)
            pop3_minor = int(0)
            
            for row in all_ind:
                if row[0] in pop1:
                    if '/' in row[i]:
                        snp = row[i].split('/')
                        if snp[0] == allele_list[i][0]:
                            pop1_major += 1
                        if snp[0] == allele_list[i][1]:
                            pop1_minor += 1
                        if snp[1] == allele_list[i][0]:
                            pop1_major += 1
                        if snp[1] == allele_list[i][1]:
                            pop1_minor += 1
                    elif '-' in row[i]:
                        pass
                    else:
                        if row[i] == allele_list[i][0]:
                            pop1_major += 2
                        elif row[i] == allele_list[i][1]:
                            pop1_minor += 2
                elif row[0] in pop2:
                    if '/' in row[i]:
                        snp = row[i].split('/')
                        if snp[0] == allele_list[i][0]:
                            pop2_major += 1
                        elif snp[0] == allele_list[i][1]:
                            pop2_minor += 1
                        if snp[1] == allele_list[i][0]:
                            pop2_major += 1
                        elif snp[1] == allele_list[i][1]:
                            pop2_minor += 1
                    elif '-' in row[i]:
                        pass
                    else:
                        if row[i] == allele_list[i][0]:
                            pop2_major += 2
                        elif row[i] == allele_list[i][1]:
                            pop2_minor += 2
                elif row[0] in pop3:
                    if '/' in row[i]:
                        snp = row[i].split('/')
                        if snp[0] == allele_list[i][0]:
                            pop3_major += 1
                        elif snp[0] == allele_list[i][1]:
                            pop3_minor += 1
                        if snp[1] == allele_list[i][0]:
                            pop3_major += 1
                        elif snp[1] == allele_list[i][1]:
                            pop3_minor += 1
                    elif '-' in row[i]:
                        pass
                    else:
                        if row[i] == allele_list[i][0]:
                            pop3_major += 2
                        elif row[i] == allele_list[i][1]:
                            pop3_minor += 2

            print "-{0}- --- {0} {1} {2} {7} {3} {4} {5} {8} {6} 15".format(allele_list[i][0],pop1_major,pop2_major,allele_list[i][1],pop1_minor,pop2_minor,loci[i],pop3_major,pop3_minor)
            fh_out.write("-{0}-\t---\t{0}\t{1}\t{2}\t{7}\t{3}\t{4}\t{5}\t{8}\t{6}\t15\n".format(allele_list[i][0],pop1_major,pop2_major,allele_list[i][1],pop1_minor,pop2_minor,loci[i],pop3_major,pop3_minor))
    fh_out.close()

#=================================================================================================
#function for dadi SNPs input file, 4 pops

def dadi4(all_ind, loci, pop1, pop2, pop3, pop4, header1, header2, header3, header4):
    outfile = "dadi_4pops_{0}_{1}_{2}_{3}_snps.txt".format(header1, header2, header3, header4)

    print '\n', "==============================================================================================================="
    print "Creating SNPs file {}...".format(outfile)
    print "===============================================================================================================", '\n'
    
    fh_out = open(outfile, 'a')
    fh_out.write("Ingroup\tOutgroup\tAllele1\t{0}\t{1}\t{2}\t{3}\tAllele2\t{0}\t{1}\t{2}\t{3}\tGene\tPosition\n".format(header1, header2, header3, header4))
    
    allele_list = ["alleles"]
    columns = len(all_ind[0])
    
    for i in range(columns):        
        allele_set = set()
        alleles = []
        if i < 1:
            pass
        elif i >= 1:
            for row in all_ind:
                if (row[0] in pop1) or (row[0] in pop2) or (row[0] in pop3) or (row[0] in pop4):
                    if '/' not in row[i]:
                        if row[i] == "-":
                            pass
                        else:
                            alleles.append(row[i])
                            alleles.append(row[i])
                            allele_set.add(row[i])                        
                    elif '/' in row[i]:
                        snp = row[i].split('/')
                        alleles.append(snp[0])
                        allele_set.add(snp[0])
                        alleles.append(snp[1])
                        allele_set.add(snp[1])

            temp_list = list(allele_set)

            if len(temp_list) == int(2):
                base1 = alleles.count(temp_list[0])
                base2 = alleles.count(temp_list[1])
                
                major_minor = []
                
                if base1 > base2:
                    major_minor.append(temp_list[0])
                    major_minor.append(temp_list[1])
                else:
                    major_minor.append(temp_list[1])
                    major_minor.append(temp_list[0])
                allele_list.append(major_minor)
                
            elif len(temp_list) < int(2):
                major_minor = []
                major_minor.append(temp_list[0])
                major_minor.append("Constant")
                allele_list.append(major_minor)

    for i in range(columns):
        if (i >= 1) and (allele_list[i][1] == "Constant"):
            print "Excluding constant site locus: {}".format(loci[i])            
        if (i >= 1) and (allele_list[i][1] != "Constant"):            
            pop1_major = int(0)
            pop1_minor = int(0)
            pop2_major = int(0)
            pop2_minor = int(0)
            pop3_major = int(0)
            pop3_minor = int(0)
            pop4_major = int(0)
            pop4_minor = int(0)
            
            for row in all_ind:
                if row[0] in pop1:
                    if '/' in row[i]:
                        snp = row[i].split('/')
                        if snp[0] == allele_list[i][0]:
                            pop1_major += 1
                        if snp[0] == allele_list[i][1]:
                            pop1_minor += 1
                        if snp[1] == allele_list[i][0]:
                            pop1_major += 1
                        if snp[1] == allele_list[i][1]:
                            pop1_minor += 1
                    elif '-' in row[i]:
                        pass
                    else:
                        if row[i] == allele_list[i][0]:
                            pop1_major += 2
                        elif row[i] == allele_list[i][1]:
                            pop1_minor += 2
                elif row[0] in pop2:
                    if '/' in row[i]:
                        snp = row[i].split('/')
                        if snp[0] == allele_list[i][0]:
                            pop2_major += 1
                        elif snp[0] == allele_list[i][1]:
                            pop2_minor += 1
                        if snp[1] == allele_list[i][0]:
                            pop2_major += 1
                        elif snp[1] == allele_list[i][1]:
                            pop2_minor += 1
                    elif '-' in row[i]:
                        pass
                    else:
                        if row[i] == allele_list[i][0]:
                            pop2_major += 2
                        elif row[i] == allele_list[i][1]:
                            pop2_minor += 2
                elif row[0] in pop3:
                    if '/' in row[i]:
                        snp = row[i].split('/')
                        if snp[0] == allele_list[i][0]:
                            pop3_major += 1
                        elif snp[0] == allele_list[i][1]:
                            pop3_minor += 1
                        if snp[1] == allele_list[i][0]:
                            pop3_major += 1
                        elif snp[1] == allele_list[i][1]:
                            pop3_minor += 1
                    elif '-' in row[i]:
                        pass
                    else:
                        if row[i] == allele_list[i][0]:
                            pop3_major += 2
                        elif row[i] == allele_list[i][1]:
                            pop3_minor += 2
                elif row[0] in pop4:
                    if '/' in row[i]:
                        snp = row[i].split('/')
                        if snp[0] == allele_list[i][0]:
                            pop4_major += 1
                        elif snp[0] == allele_list[i][1]:
                            pop4_minor += 1
                        if snp[1] == allele_list[i][0]:
                            pop4_major += 1
                        elif snp[1] == allele_list[i][1]:
                            pop4_minor += 1
                    elif '-' in row[i]:
                        pass
                    else:
                        if row[i] == allele_list[i][0]:
                            pop4_major += 2
                        elif row[i] == allele_list[i][1]:
                            pop4_minor += 2


            print "-{0}- --- {0} {1} {2} {7} {9} {3} {4} {5} {8} {10} {6} 15".format(allele_list[i][0],pop1_major,pop2_major,allele_list[i][1],pop1_minor,pop2_minor,loci[i],pop3_major,pop3_minor,pop4_major,pop4_minor)
            fh_out.write("-{0}-\t---\t{0}\t{1}\t{2}\t{7}\t{9}\t{3}\t{4}\t{5}\t{8}\t{10}\t{6}\t15\n".format(allele_list[i][0],pop1_major,pop2_major,allele_list[i][1],pop1_minor,pop2_minor,loci[i],pop3_major,pop3_minor,pop4_major,pop4_minor))
    fh_out.close()


#=================================================================================================
#=================================================================================================
#execute functions with appropriate arguments

all_ind = reformat(haplo_lines)

loci_list = loci(haplo_lines)
loci_list.insert( 0, "locus")

###########################################
#************** You'll want to edit your population names

#keep the arguments all_ind, loci_list as they are 
#the header refers to the column label in the resulting SNPs file for a population
#arguments for: dadi1(all_ind, loci,  pop1, header1)
dadi1(all_ind, loci_list, Nigeria, "Nigeria")
dadi1(all_ind, loci_list, CVL, "CVL")
dadi1(all_ind, loci_list, CVLN, "CVLN")
dadi1(all_ind, loci_list, CVLS, "CVLS")
dadi1(all_ind, loci_list, South, "South")

#the header refers to the column label in the resulting SNPs file for a population
#arguments for: dadi2(all_ind, loci,  pop1, pop2, header1, header2)
dadi2(all_ind, loci_list, Nigeria, CVLN, "Nigeria", "CVLN")
dadi2(all_ind, loci_list, CVLN, CVLS, "CVLN", "CVLS")
dadi2(all_ind, loci_list, CVL, South, "CVL", "South")
dadi2(all_ind, loci_list, South1, South2, "NorthCoast", "SouthCoast")

#arguments for: dadi3(all_ind, loci, pop1, pop2, pop3, header1, header2, header3)
dadi3(all_ind, loci_list, Nigeria, CVLN, CVLS, "Nigeria", "CVLN", "CVLS")

#arguments for: dadi4(all_ind, loci, pop1, pop2, pop3, pop4, header1, header2, header3, header4)


print '\n', '\n', '\n'
#=================================================================================================
