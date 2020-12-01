import argparse
import os
import subprocess as sp
import numpy as np
import shutil
from datetime import datetime
from collections import Counter

def get_args():
    """
    Get arguments from command line.
    """
    parser = argparse.ArgumentParser(
            description="""---------------------------------------------------------------------------
    Filter_Single_tsv - Filter a single populations.haplotypes.tsv file to remove invariant, 'blank', and 
    non-biallelic loci, and then select a single SNP per locus. The SNP selection method is set by -s, and  
    allows selecting the first SNP or a random SNP. Missing data levels are calculated on a per sample basis, 
    and a maximum level of missing data per sample can be enforced using the -m flag. After sample removal,
    the loci are re-examined and all newly invariant loci are removed. The --remove_singletons will also
    cause singletons to be removed at this step.  The output files include a filtered tsv file ("Filtered.tsv"), 
    a general log file ("Filter_Single_tsv.log"), a snp site distribution file ("SNP_distributions.log"), 
    and two missing data files (Initial_Missing_Data_Per_Sample.log, Final_Missing_Data_Per_Sample.log). 
 
    DEPENDENCIES: Python: numpy.
    ---------------------------------------------------------------------------""")
    
    parser.add_argument("-i", "--input",
                            required=True,
                            help="REQUIRED: The full path to an unfiltered populations.haplotypes.tsv file.")
    
    parser.add_argument("-m", "--missingdata",
                            required=True,
                            type=int,
                            help="REQUIRED: The maximum missing data percent to allow in "
                            "a sample. For example, entering 40 would include all samples with "
                            "<40 perc missing data and exclude all samples with >40 perc missing data.")
    
    parser.add_argument("-s", "--snpselection",
                            required=True,
                            choices=["first", "random", "all"],
                            help="REQUIRED: Specify whether to choose the first SNP, a "
                            "random SNP, or all SNPs from each available locus.")
    
    parser.add_argument("--remove_singletons",
                            required=False,
                            action='store_true',
                            help="Optional: After removal of samples failing missing data "
                            "threshold, exclude loci with singleton SNPs.")
        
    return parser.parse_args()

def print_log(LOG, string):
    """
    Function to print the input string and
    write it to the log file. Will return 
    the input string.

    Arguments:
    LOG - Full path to the log file.
    string - The string to write and print.
    """
    print(string)
    with open(LOG, 'a') as fh:
        fh.write("{}\n".format(string))
        
    return string

def tsv_to_dict(f):
    """
    Function to convert a populations.haplotypes.tsv file
    to a dictionary structure. Each key is a locus numbers and 
    the val is a list containing the 'haplotypes' (CCTT/CGTT, CCTT/CCTA, etc.)
    for all samples included. 

    Arguments:
    f - Name (not full path) of the haplotypes.tsv file.
    """
    # initiate empty dict
    tsv_dict = {}
    
    # open file
    with open(f, 'r') as fh:
        # skip header - '# Catalog Locus ID', 'Cnt', 'Sample 1', Sample 2'...
        next(fh)
        # iterate over file lines
        for line in fh:
            # use locus ID as key, set val to list of all haplotype strings
            # the int let's us sort keys in numerical order (vs strings)
            tsv_dict[int(line.split('\t')[0])] = line.strip().split('\t')[2:]
            
    return tsv_dict

def check_poly_haplo(v):
    """
    Function to check if more than 2 alleles occur
    for an individual. This happened with Stacks v1.35. 
    For example, "CCTT/CGTT/CCAT". A simple check is to
    split by '/' and get the resulting list length. So far, 
    it doesn't look this is a problem for Stacks v2+, but
    it is here just in case. Returns True if any sample has
    >2 alleles, False if all have 2 or less.

    Arguments:
    v - A list containing the 'haplotypes' of all samples for a locus.
    """
    # set variable to False
    test = False
    
    # iterate over samples/haplotypes
    for i in v:
        # if three or more haplotypes included,
        # set variable to True
        if len(i.split('/')) >= 3:
            test = True
            
    return test

def check_blank(v):
    """
    Function to check if loci consist of only '/' and/or '-'
    symbols. In Stacks v2+, loci can be filtered out during 
    populations but they will be in the haplotypes.tsv file, 
    but will be 'blank'. The list of haplotypes is converted 
    to a set and compared to sets expected from blank loci.
    Returns True if locus is blank, False if locus is present.

    Arguments:
    v - A list containing the 'haplotypes' of all samples for a locus.
    """
    # set variable to False
    test = False
    
    # sets expected from blank loci
    blank1 = set(["/"])
    blank2 = set(["/", "-"])
    
    # compare locus set to blank sets
    if set(v) == blank1 or set(v) == blank2:
        # if match found (is blank locus), set variable to True
        test = True
        
    return test

def check_biallelic(v):
    """
    Function to check if any of the SNP sites in the alleles
    are non-biallelic (e.g., three or more bases present at a single
    site). First, all unique alleles for the locus are obtained by
    putting the first allele and second allele. The alleles are all the
    same length, so it then iterates over the length of the allele 
    (e.g., iterates over each base). All the bases of that position are
    obtained from the allele set, converted to a set, and checked to see 
    if the set contains more than 2 bases (non-biallelic) or less than or
    equal to two bases (biallelic). Non-biallelic loci happened with Stacks
    v1.35, but so far this doesn't seem to occur in Stacks v2+. Returns 
    True if any site is non-biallelic and False if all sites are biallelic. 

    Arguments:
    v - A list containing the 'haplotypes' of all samples for a locus.
    """
    # set variable to False
    test = False
    
    # get all alleles from all samples
    # each list comprehension is for the first (XXX) or second (YYY) allele (XXX/YYY)
    # the resulting lists are concatenated, then reduced to a set, then turned into a sorted list
    alleles = sorted(set([x.split('/')[0] for x in v if x != "-"] + [x.split('/')[1] for x in v if x != "-"]))
    
    # iterate over bases in the allele length
    for i in range(0, len(alleles[0])):
        # get bases from all alleles in unique allele list for this position
        bases = set([a[i] for a in alleles if a[i] != "N"])
        # if nonbiallelic, set variable to True
        if len(bases) > 2:
            test = True
            
    return test
    
def filter_dict(tsv_dict, LOG):
    """
    Function to filter loci in the tsv_dict {locus_number: [haplo1, haplo2, etc]}
    to remove invariant loci, 'blank' loci, loci with > 2 alleles in a sample, 
    and non-biallelic loci. Uses dictionary comprehension for filtering, and
    returns a new dictionary each time. Returns final filtered dictionary 
    and a list of strings that summarize the filtering of each step. Takes 
    advantage of the print_log() function to print/log strings and optionally
    return the strings for particular statements.

    Arguments:
    tsv_dict - Dictionary structure of the haplotypes.tsv file, produced
               by tsv_to_dict(). Keys = loci names, vals = list of allele combos.
    LOG - Full path to the log file.
    """
    print_log(LOG, ("{}".format("-"*50)))
    print_log(LOG, ("Filtering loci in tsv file."))
    print_log(LOG, ("{}\n".format("-"*50)))
    tb = datetime.now()
    print_log(LOG, ("Beginning filtering..."))
    sum1 = print_log(LOG, ("\nNumber of starting loci: {:,}\n".format(len(tsv_dict))))

    # remove 'blank' loci using dict comprehension and check_blank() function
    rmblanks = {k:v for (k,v) in tsv_dict.items() if check_blank(v) is False}
    sum2 = print_log(LOG, ("Number of 'blank' loci: {:,}\n".format(len(tsv_dict)-len(rmblanks))))
    # remove 'blank' loci using dict comprehension and check_blank() function
    
    # remove invariant loci using dict comprehension and search for "consensus" string
    rmconsensus = {k:v for (k,v) in rmblanks.items() if "consensus" not in v}
    sum3 = print_log(LOG, ("Number of invariant loci: {:,}\n".format(len(rmblanks)-len(rmconsensus))))
    
    # remove loci with >2 alleles at a sample using dict comprehension and check_poly_haplo() function
    rmpoly = {k:v for (k,v) in rmconsensus.items() if check_poly_haplo(v) is False}
    sum4 = print_log(LOG, ("Number of loci with >2 haplotypes per sample: {:,}\n".format(len(rmconsensus)-len(rmpoly))))

    # remove non-biallelic loci using dict comprehension and check_biallelic() function
    filtered_dict = {k:v for (k,v) in rmpoly.items() if check_biallelic(v) is False}
    sum5 = print_log(LOG, ("Number of non-biallelic loci: {:,}\n".format(len(rmpoly)-len(filtered_dict))))
    
    sum6 = print_log(LOG, ("Remaining loci: {:,}\n".format(len(filtered_dict))))
    
    tf = datetime.now()
    print_log(LOG, ("\nDone.\nElapsed time: {} (H:M:S)\n".format(tf - tb)))

    return filtered_dict, [sum1, sum2, sum3, sum4, sum5, sum6]
    
def first_base(v):
    """
    Function to select the first base position of the haplotypes
    to narrow down dataset to one SNP per locus.

    Arguments:
    v - A list containing the 'haplotypes' of all samples for a locus.
    """
    # initiate empty list to store new haplotypes values in 
    new_v = []
    # get length of a single allele - the number of SNP sites
    # reverse sort haplotypes list (longest first), split first item
    # to isolate single allele, take length
    bases = len(sorted(v, reverse=True)[0].split('/')[0])
    # iterate over items in haplotypes list
    for calls in v:
        # append  all '-' characters (missing data) as is
        if "/" not in calls:
            new_v.append(calls)
        # otherwise append first position of alleles (ex., X/X from XATC/XATG)
        else:
            # if this happens to produce 'N/N', write '-' instead
            if "{}/{}".format(calls.split('/')[0][0], calls.split('/')[1][0]) == "N/N":
                new_v.append("-")
            # otherwise write new allele combo
            else:
                new_v.append("{}/{}".format(calls.split('/')[0][0], calls.split('/')[1][0]))
                
    return new_v, bases
    
def rando_base(v):
    """
    Function to select a random base position of the haplotypes
    to narrow down dataset to one SNP per locus. Uses numpy 
    random.choice() function to select random position from those
    available for a given locus. 

    Arguments:
    v - A list containing the 'haplotypes' of all samples for a locus.
    """
    # initiate empty list to store new haplotypes values in 
    new_v = []
    
    # get length of a single allele - the number of SNP sites.
    # reverse sort haplotypes list to get longest entries first (avoids '-' characters),
    # split first item by '/' to isolate a single allele, take its length
    bases = len(sorted(v, reverse=True)[0].split('/')[0])
    
    # first check to see if > 1 SNP site present
    if bases > 1:
        # select random number from range of allele length for this locus
        i = int(np.random.choice(bases, 1))
        #iterate over items in haplotypes list
        for calls in v:
            # append '-' characters (missing data) as is
            if "/" not in calls:
                new_v.append(calls)
            #otherwise take position 'i' from alleles
            else:
                # if this happens to produce 'N/N', write '-' instead
                if "{}/{}".format(calls.split('/')[0][i], calls.split('/')[1][i]) == "N/N":
                    new_v.append("-")
                # otherwise write the new allele combo
                else:
                    new_v.append("{}/{}".format(calls.split('/')[0][i], calls.split('/')[1][i]))
                    
    # if only a single SNP site present, take all haplotypes as is
    else:
        for calls in v:
            new_v.append(calls)
            
    return new_v, bases

def all_base(v):
    """
    Function to select ALL base positions of the haplotypes
    to retain all SNPs per locus.

    Arguments:
    v - A list containing the 'haplotypes' of all samples for a locus.
    """
    # initiate empty list to store new haplotypes values in 
    new_v_list = []
    
    # get length of a single allele - the number of SNP sites.
    # reverse sort haplotypes list to get longest entries first (avoids '-' characters),
    # split first item by '/' to isolate a single allele, take its length
    bases = len(sorted(v, reverse=True)[0].split('/')[0])
    
    # first check to see if > 1 SNP site present
    if bases > 1:
        # select random number from range of allele length for this locus
        for i in range(0, bases):
            new_v = []
            #iterate over items in haplotypes list
            for calls in v:
                # append '-' characters (missing data) as is
                if "/" not in calls:
                    new_v.append(calls)
                #otherwise take position 'i' from alleles
                else:
                    # if this happens to produce 'N/N', write '-' instead
                    if "{}/{}".format(calls.split('/')[0][i], calls.split('/')[1][i]) == "N/N":
                        new_v.append("-")
                    # otherwise write the new allele combo
                    else:
                        new_v.append("{}/{}".format(calls.split('/')[0][i], calls.split('/')[1][i]))
            new_v_list.append(new_v)
                    
    # if only a single SNP site present, take all haplotypes as is
    else:
        new_v = []
        for calls in v:
            new_v.append(calls)
        new_v_list.append(new_v)
            
    return new_v_list, bases

def write_base_dist(base_dist, snpdist):
    """
    Function to write the number of SNP sites per
    locus to a log file. 

    Arguments:
    base_dist - A list containing sublists for each locus: [locus name, number of snps].
    snpdist - The name of the log file to write to.
    """
    with open(snpdist, 'a') as fh:
        fh.write("Locus\tNumber_SNPs\n")
        for b in base_dist:
            fh.write("{}\t{}\n".format(b[0], b[1]))

def get_SNPs(filtered_dict, snpmethod, LOG, snpdist):
    """
    Function automate selection of one SNP per locus. Writes 
    number of SNP sites across loci to separate log file.
    Returns the tsv dictionary with one SNP per locus. 

    Arguments:
    filtered_dict - Dictionary structure of the haplotypes.tsv file, produced
                    by filter_dict().
    snpmethod - Choice from args.snpmethod: "first", or "random"
    LOG - Full path to the log file.
    snpdist - Name of snp summary log file.
    """
    print_log(LOG, ("\n\n{}".format("-"*50)))
    print_log(LOG, ("Selecting final SNP sites."))
    print_log(LOG, ("{}\n".format("-"*50)))
    print_log(LOG, ("Choosing {} SNP site for all loci...".format(snpmethod)))
    tb = datetime.now()
    
    #initiate empty dict
    snp_dict = {}
    #initiate empty list for storing loci/snp sites sublists
    base_dist = []
    
    #iterate over tsv dictionary
    for k, v in filtered_dict.items():
        
        # if first SNP, use first_base() function
        if snpmethod == "first":
            single_snps, bases = first_base(v)
            # add locus name as key and single SNP list as val to new dict
            snp_dict[k] = single_snps
            
        # if random SNP, use rando_base() function
        elif snpmethod == "random":
            single_snps, bases = rando_base(v)
            # add locus name as key and single SNP list as val to new dict
            snp_dict[k] = single_snps

        elif snpmethod == "all":
            snps_list, bases = all_base(v)
            if len(snps_list) > 1:
                for i in range(0, len(snps_list)):
                    snp_dict["{}_{}".format(k, i)] = snps_list[i]
            else:
                snp_dict["{}_0".format(k)] = snps_list[0]
        
        # append sublist info of locus name, SNP sites
        base_dist.append([k, bases])
        
    # sort list of sublists by number of SNPs per locus, low to high
    base_dist.sort(key=lambda x: int(x[1]))
    # write SNP summary log file
    write_base_dist(base_dist, snpdist)
    
    tf = datetime.now()
    print_log(LOG, ("\nDone.\nElapsed time: {} (H:M:S)\n".format(tf - tb)))
       
    return snp_dict

def remove_samples(v, bad_indices):
    """
    Function to remove all samples failing missing data threshold 
    from a given list of 'haplotypes'. The list of haplotypes is 
    always in order of the original sample order, so the bad samples 
    are located and removed based on an index. The indices for the 
    samples to remove are stored in a list, bad_indices. Returns a
    filtered haplotypes list with all bad samples removed. 

    Arguments:
    v - A list containing the 'haplotypes' of all samples for a 
        locus. At this point will be single SNPs.
    bad_indices - A list of integers which represent the indices of
                  samples which failed missing data thresholds.
    """
    # initiate empty list to store new haplotypes values in 
    new_v = []
    # use enumerate for easy indexing of haplotypes list
    for i, j in enumerate(v):
        # check if current haplotype index is in the bad sample index list
        if i not in bad_indices:
            # if not, add haplotype to new list
            new_v.append(j)
            
    return new_v

def screen_new_invariants(samples_removed_dict, remove_singletons, LOG):

    code_dict = {"A/A":"A/A",
                     "C/C":"C/C",
                     "G/G":"G/G",
                     "T/T":"T/T",
                     "A/C":"A/C",
                     "A/G":"A/G",
                     "A/T":"A/T",
                     "C/A":"A/C",
                     "C/G":"C/G",
                     "C/T":"C/T",
                     "G/A":"A/G",
                     "G/C":"C/G",
                     "G/T":"G/T",
                     "T/A":"A/T",
                     "T/C":"C/T",
                     "T/G":"G/T"}
    
    final_dict = {}
    #iterate over tsv dictionary
    for k, v in sorted(samples_removed_dict.items()):
        # get all SNP combos here, transform so all in same hetero order (so A/T == T/A)
        # this will ensure homo/hetero counts are accurate; exclude Ns and -s
        alleles = [code_dict[x] for x in v if "N" not in x and "-" not in x]

        # remove invariant loci and loci with singletons if they occur
        if remove_singletons:
            # at this point, all loci are biallelic
            # therefore, any SNP site will have at most three possible combinations
            # homozygous1, homozygous2, heterozygous: ex. A/A, T/T, A/T
            # use Counter to count number of occurrences for all unique SNP combos in the alleles list
            c = Counter(alleles)

            # if there are three combos present, then there can't be a singleton SNP, 
            # even with a minimum number for combos ('C/C', 90), ('C/G', 1), ('G/G', 1)
            # so, add key val pair to final_dict
            if len(c.most_common()) == 3:
                final_dict[k] = v
                #print("Locus {} passes:\n\t\t{}".format(k, c.most_common()))

            # if there are two combos present and the less frequent combo = 1, then
            # this is a singleton SNP
            # ex.  [('C/C', 90), ('G/G', 1)] or  [('C/C', 90), ('C/G', 1)]
            elif len(c.most_common()) == 2:
                # check if less frequent SNP combo occurs more than once, if so add key val pair
                if c.most_common(2)[1][1] > 1:
                    final_dict[k] = v
                    #print("Locus {} passes:\n\t\t{}".format(k, c.most_common()))
                else:
                    #print("Locus {} fails:\n\tsingleton: {}".format(k, c.most_common()))
                    pass
                    
            # if there is only one combo present, then it is invariant and also
            # needs to be removed
            elif len(c.most_common()) == 1:
                #print("Locus {} fails:\n\tinvariant: {}".format(k, c.most_common()))
                pass

        # just remove invariant loci, leaving loci with singletons if they occur
        else:
            # use Counter to count number of occurrences for all unique SNP combos in the alleles list
            c = Counter(alleles)

            # if there is only one SNP combo present, then it is invariant
            if len(c.most_common()) == 1:
                #print("Locus {} fails:\n\tinvariant: {}".format(k, c.most_common()))
                pass
            
            # all other possibilities involve variation
            else:
                #print("Locus {} passes:\n\t\t{}".format(k, c.most_common()))
                final_dict[k] = v
                
    print_log(LOG, ("\n\nRe-filtering loci after sample removal...\n"))
    if remove_singletons:    
        str1 = print_log(LOG, ("Number of newly invariant or singleton loci removed: {:,}\n".format(len(samples_removed_dict)-len(final_dict))))
    else:
        str1 = print_log(LOG, ("Number of newly invariant loci removed: {:,}\n".format(len(samples_removed_dict)-len(final_dict))))
    str2 = print_log(LOG, ("Final number of loci: {:,}\n".format(len(final_dict))))

    return final_dict, [str1, str2]

def missing_data(f, filtered_dict, thresh, remove_singletons, outname, LOG):
    """
    Function to calculate missing data levels on a per
    sample basis, compare values to threshold, and if 
    necessary remove samples from the tsv dictionary.
    Takes advantage of the print_log() function to 
    print/log strings and optionally return the strings 
    for particular statements. Returns a tsv dictionary 
    with bad samples removed, a list of the final sample 
    names (in order), and list of strings to write to main 
    log file. 

    Arguments:
    f - Name of haplotypes.tsv file. Used to get a list of the 
        sample names in correct order. 
    filtered_dict - Dictionary structure of the haplotypes.tsv file, produced
                    by get_SNPS().
    thresh - Missing data threshold value (int) from args.missingdata. 
    LOG - Full path to the main log file.
    """
    print_log(LOG, ("\n\n{}".format("-"*50)))
    print_log(LOG, ("Examining per-sample missing data."))
    print_log(LOG, ("{}\n".format("-"*50)))
    print_log(LOG, ("Calculating missing data..."))
    tb = datetime.now()
    
    # obtain list of sample names from haplotypes.tsv file
    with open(f, 'r') as fh:
        samples = fh.readline().strip().split('\t')[2:]
    # initiate empty list to store missing data sublists: [sample index, sample name, missing data value]
    sample_info = []
    # iterate over indices in length of sample number
    for i in range(0, len(samples)):
        # empty list to store all haplotypes for this sample
        loci = []
        # iterate over tsv dictionary
        for k, v in filtered_dict.items():
            # add haplotype strings of this sample for each locus
            # (using index value i) to the loci list
            loci.append(v[i])
        # count missing data (-) and divide by total number of loci to get % missing data
        missing = round(((float((loci.count('-') + loci.count('N/N'))) / float(len(loci))) *100), 1)
        # add [sample index, sample name, missing data value] to sample_info list
        sample_info.append([i, samples[i], missing])

    # sort larger list by sample name in sublists (element 1)
    sample_info.sort(key = lambda x: x[1])
    # write initial missing data info to output file
    with open(outname, 'a') as fh:
        fh.write("{}\t{}\n".format("Sample", "Perc_Missing_Data"))
        for s in sample_info:
            fh.write("{}\t{}\n".format(s[1], s[2]))
        
    #list comprehension to find all samples passing missing data filter
    passed = [x for x in sample_info if x[2] < thresh]
    #list comprehension to find all samples failing missing data filter
    failed = [x for x in sample_info if x[2] >= thresh]
    #both lists above also follow sublist structure: [sample index, sample name, missing data value]

    #print and log info
    sum1 = print_log(LOG, ("\n\n{0} samples PASSED missing data threshold (<{1}% missing):".format(len(passed), thresh)))
    for i in passed:
        print_log(LOG, ("\t{0}:\t{1}% missing data".format(i[1], i[2])))
    sum2 = print_log(LOG, ("\n\n{0} samples FAILED missing data threshold (>={1}% missing):".format(len(failed), thresh)))
    for i in failed:
        print_log(LOG, ("\t{0}:\t{1}% missing data".format(i[1], i[2])))
        
    print_log(LOG, ("\n\nRemoving {} failed samples from dataset...\n".format(len(failed))))
    # get bad indices from failed list sublists: [sample index, sample name, missing data value]
    bad_indices = [x[0] for x in failed]
    # initiate empty dict
    samples_removed_dict = {}
    #iterate over tsv dictionary
    for k, v in filtered_dict.items():
        # filter haplotypes list using remove_samples() function
        new_v = remove_samples(v, bad_indices)
        # add locus name as key and filtered haplotypes list as val to new dict 
        samples_removed_dict[k] = new_v

    # re-filter subset of loci to remove newly invariant loci (mandatory) and singletons (if selected)
    final_dict, summary_strs = screen_new_invariants(samples_removed_dict, remove_singletons, LOG)
        
    # get list of final samples included
    final_samples = []
    # enumerate sample_info for easy indexing
    for s in sample_info:
        # check if current index is in the bad sample index list
        if s[0] not in bad_indices:
            # if not, get sample name from sublist: [sample index, sample name, missing data value]
            final_samples.append(s[1])
            
    tf = datetime.now()
    print_log(LOG, ("\nDone.\nElapsed time: {} (H:M:S)\n".format(tf - tb)))
    
    return final_dict, final_samples, [sum1, sum2] + summary_strs

def write_final_missing_data(final_dict, final_samples, outname):
    
    with open(outname, 'a') as fh:
        fh.write("{0}\t{1}\t{2}\t{3}\n".format("Sample", "Loci_Present", "Loci_Absent", "Perc_Missing_Data"))
        # iterate over indices in length of sample number
        for i in range(0, len(final_samples)):
            # initiate empty list to store allele combos for sample
            loci = []
            # iterate over tsv dictionary
            for k, v in final_dict.items():
                # add haplotype strings of this sample for each locus
                # (using index value i) to the loci list
                loci.append(v[i])
            # count missing data
            missing = loci.count('-') + loci.count('N/N')
            # count missing data (-) and divide by total number of loci to get % missing data
            perc_missing = round(((float(missing)/float(len(loci)))*100), 1)
            # write to output
            fh.write("{0}\t{1}\t{2}\t{3}\n".format(final_samples[i], len(loci)-missing, missing, perc_missing))


def write_tsv(f, final_dict, final_samples, snpselection, missingdata, LOG):
    """
    Function to write the final filtered tsv file.

    Arguments:
    f - Name of haplotypes.tsv file.
    final_dict - The fully filtered tsv dictionary.
    final_samples - A list of the final samples included, in correct order.
    snpselection - Value of args.snpselection: "first" or "random".
    missingdata - Missing data threshold value (int) from args.missingdata. 
    LOG - Full path to the log file.
    """
    print_log(LOG, ("\n\n{}".format("-"*50)))
    print_log(LOG, ("Writing filtered loci to new file."))
    print_log(LOG, ("{}\n".format("-"*50)))
    tb = datetime.now()
    
    # create output file name
    outname = "{}.haplotypes.filtered_m{}_{}SNP.tsv".format(f.split('.')[0], missingdata, snpselection)
    print_log(LOG, ("\nWriting data to file: {}".format(outname)))
    # open output file
    with open(outname, 'a') as fh:
        # write header line
        fh.write("# Catalog Locus ID\tCnt\t{}\n".format("\t".join(final_samples)))
        # iterate over SORTED tsv dictionary
        # because keys are ints, will sort loci in numerical order 
        for k, v in sorted(final_dict.items()):
            # write locus name, updated sample count, and all haplotypes (single SNP site per sample)
            fh.write("{0}\t{1}\t{2}\n".format(k, (len(v)-v.count('-')), "\t".join(v)))
            
    print_log(LOG, ("\nFinal number of loci written: {:,}".format(len(final_dict))))
    tf = datetime.now()
    print_log(LOG, ("\nDone.\nElapsed time: {} (H:M:S)\n".format(tf - tb)))
    
    
def main():
    args = get_args()
    tb1 = datetime.now()
    argd = vars(args)

    
    #create main log file, write analysis settings
    LOG = "Filter_Single_tsv.m{}_{}SNP.log".format(args.missingdata, args.snpselection)
    if os.path.isfile(LOG):
        raise ValueError("\n\n\nERROR: A log file called {} exists already!".format("Filter_Single_tsv.m{}_{}SNP.log".format(args.missingdata, args.snpselection)))
    else:
        with open(LOG, 'a') as fh:
            fh.write("Filter_Single_tsv settings:\n\n-i (input): {0}\n-m (missingdata): {1}\n-s "
                         "(snpselection): {2}\n\n{3}\n\n".format(argd["input"], argd["missingdata"], argd["snpselection"], "="*80))
        

    print_log(LOG, ("\n\n{}".format("="*80)))
    print_log(LOG, ("\nProcessing haplotypes.tsv file: {}\n".format(args.input)))
    print_log(LOG, ("{}\n".format("="*80)))

    # convert tsv file to dictionary structure
    tsv_dict = tsv_to_dict(args.input)

    # filter the loci to remove invariants, blanks, polyhaplos, and nonbiallelics
    filtered_dict, summary = filter_dict(tsv_dict, LOG)

    # write locus filtering summary to main log file
    with open(LOG, 'a') as fh:
        fh.write("\n\nFiltering summary for: {}\n\n".format(args.input))
        for s in summary:
            fh.write("{}\n".format(s.strip()))

    if filtered_dict:
        # select one SNP per locus - first SNP or random SNP
        snp_dict = get_SNPs(filtered_dict, args.snpselection, LOG, "SNP_distributions.m{}.log".format(args.missingdata))

        # calculate per sample missing data and remove samples below threshold
        imdout = "Initial_Missing_Data_Per_Sample.m{}_{}SNP.log".format(args.missingdata, args.snpselection)
        final_dict, final_samples, sum_missing = missing_data(args.input, snp_dict, args.missingdata, args.remove_singletons, imdout, LOG)

        # write missing data filtering summary to main log file
        with open(LOG, 'a') as fh:
            for s in sum_missing:
                fh.write("\n{}\n".format(s.strip().strip(':')))
            fh.write("\n{}\n".format("-"*90))

        if len(final_dict) > 1:
            # write final missing data summary
            fmdout = "Final_Missing_Data_Per_Sample.m{}_{}SNP.log".format(args.missingdata, args.snpselection)
            write_final_missing_data(final_dict, final_samples, fmdout)

        # write final filtered tsv file
        write_tsv(args.input, final_dict, final_samples, args.snpselection, args.missingdata, LOG)
         
        
    tf1 = datetime.now()
    print("\n\n{}".format("="*80))
    print("\nTotal elapsed time: {} (H:M:S)\n".format(tf1 - tb1))
    print("{}\n\n".format("="*80))
    
if __name__ == '__main__':
    main()
