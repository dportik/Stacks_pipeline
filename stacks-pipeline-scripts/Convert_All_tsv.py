import argparse
import os
import subprocess as sp
import shutil
from random import shuffle
from datetime import datetime
from collections import Counter

def get_args():
    """
    Get arguments from command line.
    """
    parser = argparse.ArgumentParser(
            description="""---------------------------------------------------------------------------
    Convert_All_tsv - Converts all available FILTERED haplotypes.tsv files into several types of output 
    files, including phylip, fasta, nexus, structure, ped, map, and occupancy files. The structure file produced is
    compatible with the program Structure and the R package Adegenet. Two versions of the nexus file are produced,
    one containing nucleotides and the other containing integer coding (0, 1, 2) allowing use with SNAPP (via Beauti). 
    Two versions of the ped and map files 
    are also created: a 'simple' ped file in which alleles are represented by nucleotides, and a 'recoded' ped 
    file in which alleles are represented by 2 (major allele) and 1 (minor allele). The 'recoded' ped file 
    is required for the program Admixture. This script is meant to run on a directory (-i) containing the 
    outputs of the Run_Stacks.py and Filter_All_tsv.py scripts. It expects separate subdirectories to be present 
    for each population analysis (one for each -r or -R value). Each of these subdirectories must contain a FILTERED 
    version of the original haplotypes.tsv file. More than one filtered tsv file can be present, but the 
    filtered haplotypes file(s) within each subdirectory must follow this naming scheme: 
    'populations_[#].haplotypes.filtered_[#].tsv'. Here, the [#] indicates additional text (automatically 
    added by Filter_All_tsv.py). If using your own files, this naming structure must be present in order 
    for the files to be read by this script. The output files for each filtered tsv file in a given subdirectory 
    are written to a new directory ('Output-Files'). In addition, a summary file ('Convert_All_tsv.summary.txt') is 
    written to the input directory (-i), which contains the number of loci and samples contained in each tsv file.

    DEPENDENCIES: None.
    ---------------------------------------------------------------------------""")
    
    parser.add_argument("-i", "--indir",
                            required=True,
                            help="REQUIRED: The full path to the directory which contains "
                            "all of the populations output subdirectories produced by Run_Stacks.py "
                            "These subdirectories must contain at least one filtered haplotypes.tsv "
                            "file that was produced using Filter_All_tsv.py.")
            
    return parser.parse_args()

    return string

def find_dirs(indir):
    """
    Function to get full paths to all directories 
    produced by the Run_Stacks populations module. 
    These directories will be labeled as:
    Populations_R100
    Populations_R90
    Populations_R80
    ...

    Arguments:
    indir - Full path to the directory that contains the
            populations subdirectories.
    """
    # move to indir
    os.chdir(indir)
    
    # list comprehension to get all subdirs here (should only be populations-related)
    paths = sorted([os.path.abspath(f) for f in os.listdir('.') if os.path.isdir(f) and f.startswith('Populations')])
    
    print("\n\nFound {} populations directories.".format(len(paths)))
    
    return paths

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

    total_md = int(0)
    
    # open file
    with open(f, 'r') as fh:
        # skip header - '# Catalog Locus ID', 'Cnt', 'Sample 1', Sample 2'...
        next(fh)
        # iterate over file lines
        for line in fh:
            # use locus ID as key, set val to list of all haplotype strings
            # the int let's us sort keys in numerical order (vs strings)
            # to make things easier downstream, any "-" characters will be
            # transformed to "-/-"
            if line.strip().split('\t')[2:]:
                tsv_dict[int(line.split('\t')[0])] = [i.replace("-", "-/-") for i in line.strip().split('\t')[2:]]
                
                # get missing data for locus, add to total
                total_md += ([i.replace("-", "-/-") for i in line.strip().split('\t')[2:]].count("-/-") +
                                 [i.replace("-", "-/-") for i in line.strip().split('\t')[2:]].count("N/N"))
                
    # obtain list of sample names from haplotypes.tsv file
    with open(f, 'r') as fh:
        # access first line with readline, strip, create list from
        # split by tabs, grab only element 2 and onwards
        samples = fh.readline().strip().split('\t')[2:]

    print("\n\tFound {:,} loci.".format(len(tsv_dict)))
    print("\tFound {:,} samples.".format(len(samples)))
    
    if len(tsv_dict) == 0:
        print("\tTotal possible SNPs: {}".format("NA"))
        print("\tMissing SNPs: {}".format("NA"))
        print("\tTotal percent missing data: {}".format("NA"))
        
        return tsv_dict, samples, [f, len(tsv_dict), len(samples), "NA", "NA", "NA"]
        
    else:
        print("\tTotal possible SNPs: {:,}".format(len(tsv_dict) * len(samples)))
        print("\tMissing SNPs: {:,}".format(total_md))
        print("\tTotal percent missing data: {}".format( round(((float(total_md) / float(len(tsv_dict) * len(samples))) * 100 ), 1) ))
        
        return tsv_dict, samples, [f, len(tsv_dict), len(samples), (len(tsv_dict) * len(samples)),
                                       total_md, round(((float(total_md) / float(len(tsv_dict) * len(samples))) * 100 ), 1)]

def dict_to_samples(tsv_dict, samples):
    """
    Function to convert the tsv_dict (keys = loci names; 
    vals = lists of the allele combos for all samples) to a 
    dictionary in which the keys are sample names and the vals 
    are lists of the allele combos across all loci (in numerical 
    order). Returns the sample dictionary

    Arguments:
    tsv_dict - Dictionary structure of the haplotypes.tsv file, produced
               by tsv_to_dict(). Keys = loci names, vals = list of allele combos.
    samples - A list of the samples included in the tsv file, in correct order.
    """
    # initiate empty dict for sample:loci pairs
    sample_dict = {}
    # iterate over indices in length of sample number
    for i in range(0, len(samples)):
        # empty list to store all haplotypes for this sample
        loci = []
        # iterate over tsv dictionary sorted by locus number
        for k, v in sorted(tsv_dict.items()):
            # add haplotype strings of this sample for each locus
            # (using index value i) to the loci list
            loci.append(v[i])
        # add sample name as key, and new loci list as val
        sample_dict[samples[i]] = loci
        
    return sample_dict
        
def write_fasta_phy_nex(sample_dict, samples, label, outdir):
    """
    Function to write phylip, fasta, and nexus files to the 
    specified output directory using the sample_dict derived 
    from a particular haplotypes.tsv file. 

    Arguments:
    sample_dict - Dictionary structure of the haplotypes.tsv file in which 
                  key = sample name, val = list of allele combos for sample 
                  across all loci (loci in list are in numerical order).
    samples - A list of the samples included in the tsv file, in correct order.
    label - A string for output file naming, constructed from parsed tsv name.
    outdir - Directory to write output files in.
    """
    # dictionary of conversion values
    code_dict = {"A/A":"A",
                     "C/C":"C",
                     "G/G":"G",
                     "T/T":"T",
                     "A/C":"M",
                     "A/G":"R",
                     "A/T":"W",
                     "C/A":"M",
                     "C/G":"S",
                     "C/T":"Y",
                     "G/A":"R",
                     "G/C":"S",
                     "G/T":"K",
                     "T/A":"W",
                     "T/C":"Y",
                     "T/G":"K",
                     "N/N":"N",
                     "-/-":"N"}

    # initiate empty dictionary for converted codes
    converted_dict = {}
    # iterate over sample dict
    for k, v in sample_dict.items():
        # assign sample as k, assign list of code_dict vals as v
        converted_dict[k] = [code_dict[i] for i in v]
        
    # get number of loci from a sample in dict
    num_loci = len(sample_dict[samples[0]])

    # write phylip
    phy = "{}.phy".format(label)
    with open(phy, 'a') as fh:
        fh.write("{0} {1}\n".format(len(samples), num_loci))
        for k, v in sorted(converted_dict.items()):
            fh.write("{0} {1}\n".format(k, "".join(v)))
    print("\n\tWrote phylip file:\t{}".format(phy))

    # write fasta
    fas = "{}.fasta".format(label)
    with open(fas, 'a') as fh:
        for k, v in sorted(converted_dict.items()):
            fh.write(">{0}\n{1}\n".format(k, "".join(v)))
    print("\tWrote fasta file:\t{}".format(fas))

    # write nexus
    nex = "{}.nex".format(label)
    with open(nex, 'a') as fh:
        fh.write('''
#NEXUS 
BEGIN DATA;
	DIMENSIONS  NTAX={0} NCHAR={1};
	FORMAT DATATYPE=DNA  MISSING=N GAP=-;
MATRIX
'''.format(len(samples), num_loci))
        for k, v in sorted(converted_dict.items()):
            fh.write("{0} {1}\n".format(k, "".join(v)))
        fh.write("\n;\nEnd;")
    print("\tWrote nexus file:\t{}".format(nex))

    # move output(s) to correct directory
    for f in [phy, fas, nex]:
        shutil.move(f, outdir)
        
def write_struct(sample_dict, samples, label, outdir):
    """
    Function to write a structure output file to the 
    specified output directory using the sample_dict derived 
    from a particular haplotypes.tsv file. 

    Arguments:
    sample_dict - Dictionary structure of the haplotypes.tsv file in which 
                  key = sample name, val = list of allele combos for sample 
                  across all loci (loci in list are in numerical order).
    samples - A list of the samples included in the tsv file, in correct order.
    label - A string for output file naming, constructed from parsed tsv name.
    outdir - Directory to write output files in.
    """
    
    # dictionary of conversion values
    code_dict = {"A":"1",
                       "T":"2",
                       "C":"3",
                       "G":"4",
                       "N":"-9",
                       "-":"-9"}
        
    # initiate empty dictionary for converted codes
    converted_dict = {}
    
    # iterate over sample dict
    for k, v in sample_dict.items():
        # get SNPs from first allele of haplotypes for this sample, labeled A
        converted_dict["{}_A".format(k)] = [code_dict[i.split('/')[0]] for i in v]
        # get SNPs from second allele of haplotypes for this sample, labeled A
        converted_dict["{}_B".format(k)] = [code_dict[i.split('/')[1]] for i in v]

    # write structure
    struct = "{}.str".format(label)
    with open(struct, 'a') as fh:
        for k, v in sorted(converted_dict.items()):
            fh.write("{0}\t{1}\n".format(k.replace("_A", "").replace("_B",""), "\t".join(v)))
            
    print("\tWrote structure file:\t{}".format(struct))
    
    # move output(s) to correct directory
    shutil.move(struct, outdir)

def write_simple_ped(sample_dict, tsv_dict, samples, label, outdir):
    """
    Function to write a simple ped and map file to the 
    specified output directory using the sample_dict derived 
    from a particular haplotypes.tsv file. The ped file will 
    have actual nucleotide codes in it. 

    Arguments:
    sample_dict - Dictionary structure of the haplotypes.tsv file in which 
                  key = sample name, val = list of allele combos for sample 
                  across all loci (loci in list are in numerical order).
    tsv_dict - Dictionary structure of the haplotypes.tsv file, produced
               by tsv_to_dict(). Keys = loci names, vals = list of allele combos.
    samples - A list of the samples included in the tsv file, in correct order.
    label - A string for output file naming, constructed from parsed tsv name.
    outdir - Directory to write output files in.
    """
    
    # dictionary of conversion values
    code_dict = {"A/A":"A A",
                     "C/C":"C C",
                     "G/G":"G G",
                     "T/T":"T T",
                     "A/C":"A C",
                     "A/G":"A G",
                     "A/T":"A T",
                     "C/A":"C A",
                     "C/G":"C G",
                     "C/T":"C T",
                     "G/A":"G A",
                     "G/C":"G C",
                     "G/T":"G T",
                     "T/A":"T A",
                     "T/C":"T C",
                     "T/G":"T G",
                     "N/N":"0 0",
                     "-/-":"0 0"}

    # initiate empty dictionary for converted codes
    converted_dict = {}
    
    # iterate over sample dict
    for k, v in sample_dict.items():
        # assign sample as k, assign list of code_dict vals as v
        converted_dict[k] = [code_dict[i] for i in v]
    
    # write ped file
    ped = "{}.ped".format(label)
    with open(ped, 'a') as fh:
        for k, v in sorted(converted_dict.items()):
            fh.write("{0} {0} {1} {1} {1} {1} {2}\n".format(k, "0", " ".join(v)))

    # write useless map file
    dumbmap = "{}.map".format(label)
    with open(dumbmap, 'a') as fh:
        for k, v in sorted(tsv_dict.items()):
            fh.write("{0} Locus_{1} {0} {1}\n".format("0", k))
    
    print("\tWrote ped file:\t\t{}".format(ped))
    print("\tWrote ped map file:\t{}".format(dumbmap))

    # move output(s) to correct directory
    for f in [ped, dumbmap]:
        shutil.move(f, outdir)
    
def write_recoded_ped(tsv_dict, samples, label, outdir):
    """
    Function to write a recoded (012) ped and map file to the 
    specified output directory using the sample_dict derived 
    from a particular haplotypes.tsv file. This ped file is 
    coded, with 2 representing major allele, 1 representing 
    minor allele, 0 representing missing data. Major and minor 
    alleles are determined for each locus based on allele counts.

    Arguments:
    tsv_dict - Dictionary structure of the haplotypes.tsv file, produced
               by tsv_to_dict(). Keys = loci names, vals = list of allele combos.
    samples - A list of the samples included in the tsv file, in correct order.
    label - A string for output file naming, constructed from parsed tsv name.
    outdir - Directory to write output files in.
    """

    # initiate empty dictionary for converted codes
    converted_dict = {}
    
    # iterate over tsv dictionary sorted by locus number
    for k, v in sorted(tsv_dict.items()):
        
        # get all real SNP bases for this locus (exclude N's and -'s)
        bases = ([x.split('/')[0] for x in v if x != "N/N" and x != "-/-"] +
                     [x.split('/')[1] for x in v if x != "N/N" and x != "-/-"])
        
        # use Counter to count number of occurrences for each unique element in list
        c = Counter(bases)
        
        # if no bases found
        # if tsv's were filtered with Filter_All_tsv.py, this shouldn't occur
        # because these loci are all missing data; include just in case
        if len(c.most_common(2)) == 0:
            # list comprehension to replace all bases with numerical codes
            newv = [x.replace("N", "0").replace("-", "0") for x in v]
            
        # if only one base found
        # if tsv's were filtered with Filter_All_tsv.py, this shouldn't occur
        # because these loci are invariant; include just in case
        elif len(c.most_common(2)) == 1:
            # assign most common base to major ("2")
            major = c.most_common(1)[0][0]
            # list comprehension to replace all bases with numerical codes
            newv = [x.replace(major, "2").replace("N", "0").replace("-", "0") for x in v]
            
        # if two bases found
        elif len(c.most_common(2)) == 2:
            # assign most common base to major ("2")
            major = c.most_common(2)[0][0]
            # assign less common base to minor ("1")
            minor = c.most_common(2)[1][0]
            # list comprehension to replace all bases with numerical codes
            newv = [x.replace(major, "2").replace(minor, "1").replace("N", "0").replace("-", "0") for x in v]

        # add new key val pair to converted dict
        converted_dict[k] = newv

    # convert locus dict to sample dict
    sample_dict = dict_to_samples(converted_dict, samples)

    # dictionary of conversion values
    code_dict = {"1/1":"1 1",
                     "1/2":"1 2",
                     "2/1":"1 2",
                     "2/2":"2 2",
                     "0/0":"0 0"}
        
    # initiate empty dictionary for converted codes
    converted_sample_dict = {}
    
    # iterate over sample dict
    for k, v in sample_dict.items():
        # assign sample as k, assign list of code_dict vals as v
        converted_sample_dict[k] = [code_dict[i] for i in v]
    
    # write ped file
    ped = "{}_recoded.ped".format(label)
    with open(ped, 'a') as fh:
        for k, v in sorted(converted_sample_dict.items()):
            fh.write("{0} {0} {1} {1} {1} {1} {2}\n".format(k, "0", " ".join(v)))
            
    # write useless map file
    dumbmap = "{}_recoded.map".format(label)
    with open(dumbmap, 'a') as fh:
        for k, v in sorted(tsv_dict.items()):
            fh.write("{0} Locus_{1} {0} {1}\n".format("0", k))
    
            
    print("\tWrote recoded ped file:\t{}".format(ped))
    print("\tWrote recoded map file:\t{}".format(dumbmap))
    
    # move output(s) to correct directory
    for f in [ped, dumbmap]:
        shutil.move(f, outdir)
        
def write_snapp_nexus(tsv_dict, samples, label, outdir):
    """
    Function to write SNAPP version of nexus file, in which
    homozygous SNPs are written as 0 or 2, and heterozygous
    SNPs are written as 1. The bases are assigned randomly 
    to prevent any biases (e.g., always turning A into a 0). 

    Arguments:
    tsv_dict - Dictionary structure of the haplotypes.tsv file, produced
               by tsv_to_dict(). Keys = loci names, vals = list of allele combos.
    samples - A list of the samples included in the tsv file, in correct order.
    label - A string for output file naming, constructed from parsed tsv name.
    outdir - Directory to write output files in.
    """
    
    # initiate empty dictionary for converted codes
    converted_dict = {}
    
    # iterate over tsv dictionary sorted by locus number
    for k, v in sorted(tsv_dict.items()):
        
        # get all real SNP bases for this locus (exclude N's and -'s)
        # convert to set and sort, should be length of 2 now
        bases = sorted(set([x.split('/')[0] for x in v if x != "N/N" and x != "-/-"] +
                     [x.split('/')[1] for x in v if x != "N/N" and x != "-/-"]))

        # shuffle the bases
        shuffle(bases)
        # define new codes for homo and hetero SNPs, missing data
        temp_dict = {"{0}/{0}".format(bases[0]):"0",
                         "{0}/{0}".format(bases[1]):"2",
                         "{0}/{1}".format(bases[0], bases[1]):"1",
                         "{1}/{0}".format(bases[0], bases[1]):"1",
                         "N/N":"-",
                         "-/-":"-"}
        # replace old SNPs with new codes, add to dict
        converted_dict[k] = [temp_dict[i] for i in v]
        
    # convert locus dict to sample dict
    sample_dict = dict_to_samples(converted_dict, samples)
    
    # get number of loci from a sample in dict
    num_loci = len(sample_dict[samples[0]])

    # write nexus
    nex = "{}_SNAPP.nex".format(label)
    with open(nex, 'a') as fh:
        fh.write('''
#NEXUS 
BEGIN DATA;
	DIMENSIONS  NTAX={0} NCHAR={1};
	FORMAT DATATYPE=INTEGERDATA SYMBOLS="012" GAP=- ;
MATRIX
'''.format(len(samples), num_loci))
        for k, v in sorted(sample_dict.items()):
            fh.write("{0} {1}\n".format(k, "".join(v)))
        fh.write("\n;\nEnd;")
    print("\tWrote SNAPP nexus file:\t{}".format(nex))

    # move output to correct directory
    shutil.move(nex, outdir)
        
def write_occupancy(sample_dict, samples, label, outdir):
    """
    Function to write an "occupancy" file, for use with specific 
    data visualization tool: 
    https://bmedeiros.shinyapps.io/matrix_condenser/

    Arguments:
    sample_dict - Dictionary structure of the haplotypes.tsv file in which 
                  key = sample name, val = list of allele combos for sample 
                  across all loci (loci in list are in numerical order).
    samples - A list of the samples included in the tsv file, in correct order.
    label - A string for output file naming, constructed from parsed tsv name.
    outdir - Directory to write output files in.
    """
    # dictionary of conversion values
    code_dict = {"A/A":"1",
                     "C/C":"1",
                     "G/G":"1",
                     "T/T":"1",
                     "A/C":"1",
                     "A/G":"1",
                     "A/T":"1",
                     "C/A":"1",
                     "C/G":"1",
                     "C/T":"1",
                     "G/A":"1",
                     "G/C":"1",
                     "G/T":"1",
                     "T/A":"1",
                     "T/C":"1",
                     "T/G":"1",
                     "N/N":"0",
                     "-/-":"0"}

    # initiate empty dictionary for converted codes
    converted_dict = {}
    # iterate over sample dict
    for k, v in sample_dict.items():
        # assign sample as k, assign list of code_dict vals as v
        converted_dict[k] = [code_dict[i] for i in v]
        
    # get number of loci from a sample in dict
    num_loci = len(sample_dict[samples[0]])
    header = ",".join([str(i) for i in list(range(0, num_loci+1))])

    # write file
    occ = "{}.occupancy.csv".format(label)
    with open(occ, 'a') as fh:
        fh.write("{}\n".format(header))
        for k, v in sorted(converted_dict.items()):
            fh.write("{0},{1}\n".format(k, ",".join(v)))
    print("\tWrote occupancy file:\t{}".format(occ))
    shutil.move(occ, outdir)
    
def check_outputs(maindir, outdir, label):
    """
    Function to check if output files already exist 
    in the output directory (to prevent over-writing).
    Throws error if putative output files are found.

    Arguments:
    maindir - The populations subdirectory.
    outdir - Directory to write output files in.
    label - A string for output file naming, constructed from parsed tsv name.
    """
    # move to output directory
    os.chdir(outdir)
    # find all files that begin with the output file prefix
    check = [f for f in os.listdir('.') if f.startswith(label)]
    # move back to main directory
    os.chdir(maindir)
    # if output files found, raise error
    if check:
        raise ValueError("\n\n\nERROR: Output files already exist for {} in "
                             "directory:\n\n\t'{}'.\n\nPlease remove before running!\n".format(f, outdir))        
    
def main():
    args = get_args()
    tb1 = datetime.now()

    #find all populations subdirectories
    paths = find_dirs(args.indir)

    summary_info = []
    
    # iterate over paths to populations subdirectories
    for p in paths:
        print("\n\n{}".format("-"*60))
        print("Converting filtered tsv files in directory: {}".format(p.split('/')[-1]))
        print("{}\n".format("-"*60))
        
        # move to subdir
        os.chdir(p)
        # locate all the filtered haplotypes.tsv files
        # second list comprehension (to find filtered tsv) performed on initial list comprehension (to find all tsv files)
        hapfiles = [f for f in [f for f in os.listdir('.') if f.endswith('.tsv')] if f.split('.')[-2].startswith("filtered_")]

        # create output directory
        outdir = os.path.join(p, "Output-Files")
        if not os.path.exists(outdir):
            os.mkdir(outdir)            

        #iterate over haplotypes files
        for f in hapfiles:
            tb2 = datetime.now()
            print("\nProcessing {}...".format(f))
            
            # get label for output file naming
            label = "_".join([f.split('.')[0], f.split('.')[-2].replace("filtered_", "")])
            
            # ensure output files don't yet exist
            check_outputs(p, outdir, label)
            
            # convert tsv file to dictionary structure
            tsv_dict, samples, info = tsv_to_dict(f)
            summary_info.append(info)

            if tsv_dict:
                # convert tsv_dict into one with samples as keys
                sample_dict = dict_to_samples(tsv_dict, samples)

                # write fasta, phylip, and nexus files
                write_fasta_phy_nex(sample_dict, samples, label, outdir)
                
                # write SNAPP nexus file
                write_snapp_nexus(tsv_dict, samples, label, outdir)

                # write structure file
                write_struct(sample_dict, samples, label, outdir)

                # write simple ped file
                write_simple_ped(sample_dict, tsv_dict, samples, label, outdir)

                # write recoded ped file
                write_recoded_ped(tsv_dict, samples, label, outdir)
                
                # write occupancy file
                write_occupancy(sample_dict, samples, label, outdir)
                
                tf2 = datetime.now()
                print("\n\tDone. Elapsed time:\t{}\n".format(tf2 - tb2))
               
            else:
                print("\n\n\tThis tsv file appears to be empty!\n\t\tSkipping.\n\n")
            
    os.chdir(args.indir)
    with open("Convert_All_tsv.summary.txt", 'a') as fh:
        fh.write("{}\t{}\t{}\t{}\t{}\t{}\n".format("tsv_file", "Loci",
                                                       "Samples", "Total_SNP_sites",
                                                       "Missing_SNP_sites", "Perc_Missing_Data"))
        for i in summary_info:
            fh.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(i[0], i[1], i[2], i[3], i[4], i[5]))
    
    tf1 = datetime.now()
    print("\n\n{}".format("="*80))
    print("\nTotal elapsed time: {} (H:M:S)\n".format(tf1 - tb1))
    print("{}\n\n".format("="*80))
    
if __name__ == '__main__':
    main()
