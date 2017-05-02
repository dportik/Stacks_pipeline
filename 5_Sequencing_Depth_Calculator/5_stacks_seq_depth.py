import sys
import os
import shutil
import datetime
import numpy as np
from datetime import datetime
'''
Usage: python stacks_seq_depth.py [full path to a fully filtered 'batch_#.haplotypes.tsv' file] [directory with all "SAMPLE.matches.tsv" files]

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

Sample	Avg	Min	Max	Median	SD
1_AM224	42.0	1.0	537.0	16.0	54.7
1_AM225	52.7	1.0	558.0	20.0	70.4
1_AM226	54.2	1.0	492.0	20.0	73.7
1_AM227	46.4	1.0	526.0	17.0	62.3



Sequencing depth stats on a per locus basis. Based on 3) only when the locus is
present in the samples, or 4) when the locus is missing in a sample it is included 
as depth value of zero.

3. Output_Loci_NoMissingData_Coverages.txt
4. Output_Loci_MissingData_Coverages.txt

example output:

Locus	Avg	Min	Max	Median	SD
1	31.2	3.0	135.0	23.0	26.5
5	17.6	4.0	47.0	15.0	10.3
6	179.8	7.0	673.0	166.0	137.7
9	146.8	6.0	745.0	128.0	109.1
11	24.5	5.0	76.0	20.0	15.4
12	14.0	3.0	38.0	12.0	8.0



These ouptut files can be used to make plots in R to visualize the data. 

############################################
Written for Python 2.7
External Dependencies: Numpy (Numerical Python)
############################################

Dan Portik
daniel.portik@uta.edu
April 2017
'''
t_begin = datetime.now()

#=================================================================================================
#function to convert lists into arrays and do very basic stats with numpy
def quickstats(x):
    x_array = np.asarray(x,dtype=np.float64)
    x_avg = np.average(x_array)
    x_avg = np.around(x_avg, decimals = 1)
    x_min = np.amin(x_array)
    x_max = np.amax(x_array)
    x_med = np.median(x_array)
    x_std = np.std(x_array)
    x_std = np.around(x_std, decimals = 1)
    return x_avg, x_min, x_max, x_med, x_std
#=================================================================================================

#get haplotype file path and get catalog IDs to search
hap_file = sys.argv[1]
loci_list = []
with open(hap_file, 'r') as fh_temp:
    hap_lines = fh_temp.readlines()
    hap_lines = [line.strip() for line in hap_lines]
    for line in hap_lines[1:]:
        line = line.split('\t')
        loci_list.append(line[0])

#get input dir location
in_dir = sys.argv[2]
os.chdir(in_dir)

#get list of files to search
files_list = []
for filetype in os.listdir('.'):
    if filetype.endswith("matches.tsv"):
        files_list.append(filetype)

#====================================================
#start with calculating coverage per individual
print '\n','\n',"====================================================",'\n',"Beginning sequencing depth calculations for {} samples...".format(len(files_list)), '\n'
fh_ind1 = open("Output_Sample_NoMissingData_Coverages.txt", 'a')
fh_ind1.write("Sample\tLociPresent\tLociAbsent\tAvg\tMin\tMax\tMedian\tSD\n")
fh_ind2 = open("Output_Sample_WithMissingData_Coverages.txt", 'a')
fh_ind2.write("Sample\tLociPresent\tLociAbsent\tAvg\tMin\tMax\tMedian\tSD\n")

#make list of lists (locus name and depth), across all individuals, for later 
locdep_missing_list = []
locdep_nomissing_list = []

for files in files_list:
    #split filename to get sample name access
    names = files.split(".matches")
    #list to store only depth values for this individual
    depth_ind_missing = []
    depth_ind_nomissing = []
    #search file, match catalog ID to loci list and record the stack depth/coverage
    print datetime.now().strftime("%b-%d %H:%M")
    print "  Processing {}".format(files)
    with open(files, 'r') as fh_temp:
        lines = fh_temp.readlines()
        lines = [line.strip() for line in lines]
        miss_count = int(0)
        for locus in loci_list:
            #list to store locus name and depth values for this individual
            temp_list = []
            depth_val = int(0)
            for line in lines[1:]:
                line = line.split('\t')
                if line[2] == locus:
                    depth_val += int(line[6])
            if depth_val == int(0):
            	miss_count += 1
            	depth_ind_missing.append(depth_val)
            	temp_list.append(locus)
            	temp_list.append(depth_val)
            	locdep_missing_list.append(temp_list)
            else:
            	depth_ind_nomissing.append(depth_val)
            	depth_ind_missing.append(depth_val)
            	temp_list.append(locus)
            	temp_list.append(depth_val)
            	locdep_missing_list.append(temp_list)
            	locdep_nomissing_list.append(temp_list)

    print '\t',"Recorded {0} locus depth values and {1} missing loci for sample: {2}".format(len(depth_ind_nomissing),miss_count,names[0])
    #stats on depth/coverage
    avg1, min1, max1, med1, std1 = quickstats(depth_ind_nomissing)
    avg2, min2, max2, med2, std2 = quickstats(depth_ind_missing)
    print '\t',"Average sequencing depth for only loci present: {}".format(avg1)
    print '\t',"Average sequencing depth when including missing data: {}".format(avg2), '\n', '\n'
    #add relevant results to output file
    fh_ind1.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(names[0],len(depth_ind_nomissing),miss_count,avg1,min1,max1,med1,std1))
    fh_ind2.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(names[0],len(depth_ind_nomissing),miss_count,avg2,min2,max2,med2,std2))

fh_ind1.close()
fh_ind2.close()
    
#====================================================
#finish with calculating coverage per locus
print '\n','\n',"====================================================",'\n',"Beginning sequencing depth calculations of each locus...", '\n'
fh_loc1 = open("Output_Loci_NoMissingSamples_Coverages.txt", 'a')
fh_loc1.write("Locus\tAvg\tMin\tMax\tMedian\tSD\n")
fh_loc2 = open("Output_Loci_WithMissingSamples_Coverages.txt", 'a')
fh_loc2.write("Locus\tAvg\tMin\tMax\tMedian\tSD\n")

for locus in loci_list:
    depth_loc_nomissing = []
    depth_loc_missing = []
    #search through all file lines from all samples, match locus and record depths
    for item in locdep_nomissing_list:
        if item[0] == locus:
            depth_loc_nomissing.append(item[1])
    for item in locdep_missing_list:
        if item[0] == locus:
            depth_loc_missing.append(item[1])        
                    
    print '\t',"Found {0} present depth values for locus: {1}".format((len(depth_loc_nomissing)),locus)
    #stats on depth/coverage
    if len(depth_loc_nomissing) == int(0):
    	avg1=min1=max1=med1=std1=int(0)
    else:
    	avg1, min1, max1, med1, std1 = quickstats(depth_loc_nomissing)
    	
    if len(depth_loc_missing) == int(0):
    	avg2=min2=max2=med2=std2=int(0)
    else:
    	avg2, min2, max2, med2, std2 = quickstats(depth_loc_missing)

    print '\t',"Average sequencing depth for locus only when present: {}".format(avg1)
    print '\t',"Average sequencing depth for locus when including missing data: {}".format(avg2), '\n', '\n'
    #add relevant results to output file
    fh_loc1.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(locus,avg1,min1,max1,med1,std1))
    fh_loc2.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(locus,avg2,min2,max2,med2,std2))
    
fh_loc1.close()
fh_loc2.close()

#====================================================
#clean up output files, move to new directory, show elapsed time
out_dir = "Coverage_Calculations"
if not os.path.exists(out_dir):
    os.mkdir(out_dir)
for filetype in os.listdir('.'):
    if filetype.endswith("_Coverages.txt"):
        shutil.move(filetype, out_dir)

t_finish = datetime.now()
elapsed = t_finish - t_begin

print '\n', '\n', "-----------------------------------------------------------------------------------------------------"
print "Total time: {0} (H:M:S)".format(elapsed)
print "-----------------------------------------------------------------------------------------------------", '\n', '\n'
