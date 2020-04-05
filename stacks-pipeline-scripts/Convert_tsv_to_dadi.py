import argparse
import os
from datetime import datetime
from collections import Counter

def get_args():
    """
    Get arguments from command line.
    """
    parser = argparse.ArgumentParser(
            description="""---------------------------------------------------------------------------
    Convert_tsv_to_dadi: Converts a FILTERED haplotypes.tsv file into a SNPs input file 
    that can be used with dadi and moments, two demographic modeling programs. This should be 
    used after running the Filter_All_tsv.py script. In addition to the FILTERED tsv file (-i), you 
    will need to provide a population assignment file (-p). This simple text file contains one 
    line per population. On each line, the name of the population should be given, followed by a colon, 
    followed by all the sample names included in the population, which should be separated by commas. 
    Lines with a hash (#) symbol will be ignored. The script will automatically determine the number of 
    populations based on the assignment file, and write a SNPs input file from those populations. To 
    quickly write a file containing a list of the samples in a given tsv file, use the --getnames flag. 
    ---------------------------------------------------------------------------""")
    parser.add_argument("-i", "--input",
                            required=True,
                            help="Required: The full path to a filtered tsv file to convert "
                            "into a SNPs input file for dadi or moments.")
    
    parser.add_argument("-o", "--outdir",
                            required=True,
                            help="REQUIRED: The full path to an existing directory "
                            "to write the output files. Must be used with -p flag.")
    
    parser.add_argument("-p", "--popfile",
                            required=False,
                            default=None,
                            help="OPTIONAL: Provide the full path to a population assignment file to "
                            "use to generate the SNPS input file. Must be used with -o flag.")
    
    parser.add_argument("--getnames",
                            required=False,
                            action='store_true',
                            help="OPTIONAL: Write an output file containing all the sample names "
                             "in the tsv file, then quit.")
    
    parser.add_argument("--checkpopfile",
                            required=False,
                            action='store_true',
                            help="OPTIONAL: Print population assignments from file to screen and quit. "
                            "Must be used with -p.")
    
    parser.add_argument("--quiet",
                            required=False,
                            action='store_true',
                            help="OPTIONAL: Make output less verbose (do not per-locus summary to screen).")

           
    return parser.parse_args()


def find_pops(popfile):
    """
    Function to read the popfile to identify the population 
    labels and the samples they include. Returns a dictionary 
    where key is popname and val is list of sample names. Also 
    returns a list containing the population names.

    Popfile content is one line per population, beginning with the name 
    of the population ending with colon, followed by all sample names separated 
    by commas. For example:
    Popname1: Sample1, Sample2, Sample4
    Popname2: Sample3, Sample5, Sample6, Sample12

    Arguments: 
    popfile - Full path to the pop assignment file.
    """
    # initialize empty dicts
    popdict = {}
    sampledict = {}
    # initialize empty list
    popnames = []

    # open pop file
    with open(popfile, 'r') as fh:
        # iterate through lines
        for line in fh:
            if line.strip() and not line.startswith('#'):
                # get list of samples from line
                samples = [x.strip() for x in line.split(':')[1].split(',') ]
                pop = line.split(':')[0].strip()
                # get popname as key and list of sample names as val
                popdict[pop] = samples
                # add sample name to sampledict, with popname as val
                for s in samples:
                    sampledict[s] = pop
                # add population name to popname list
                popnames.append(pop)
            
    print("\n\n{}".format("-"*80))
    print("\n\nFound {} populations in the assignment file.\n".format(len(popnames)))
    for p in popnames:
        print("\n{} contains {} samples:".format(p, len(popdict[p])))
        for s in sorted(popdict[p]):
            print("\t{}".format(s))
    print("\n{}\n\n".format("-"*80))
    
    return popnames, popdict, sampledict

def get_sample_names(tsv):
    """
    Function to read header line of tsv and return list of sample names.

    tsv header format:
    # Catalog Locus ID	Cnt	SAMPLE1	SAMPLE2	SAMPLE3
    where split on tabs = ["# Catalog Locus ID", "Cnt", "SAMPLE1", "SAMPLE2", "SAMPLE3"]

    Arguments:
    tsv - Full path to the filtered haplotypes tsv file.
    """
    # open tsv
    with open(tsv, 'r') as fh:
        # get first line with readline, get cleaned list of sample names
        samples = fh.readline().strip().split('\t')[2:]
        
    return samples

def write_sample_names(tsv, outdir):
    """
    Function to write sample names to file.

    Arguments:
    tsv - Full path to the filtered haplotypes tsv file.
    outdir - Directory to write output file to.
    """
    # use function to return list of sample names
    samples = get_sample_names(tsv)
    # move to output directory
    os.chdir(outdir)
    # open output file
    with open("Sample_Names.txt", 'a') as fh_out:
        # iterate over list of sample names
        for s in samples:
            # write sample name and line break to output file
            fh_out.write("{}\n".format(s))

def compare_sampling(popnames, popdict, samples):
    """
    Function to ensure that each sample name from the pop assignment
    file has a corresponding match to the tsv file. If any names do 
    not match for any of the populations, raises an error and prints the
    unmatched names. Uses set comparisons. Does not return anything.

    Arguments:
    popnames - List of population names.
    popdict - Dictionary with (key = population name) and (val = list of sample names)
    samples - List of samples in tsv file.
    """
    print("\n\n{}".format("-"*80))
    print("\nChecking sample name compatibility between assignment file and tsv file:\n")

    # create set from sample list obtained from tsv file
    allsamples = set(samples)

    # iterate over population names, in order to access dict keys
    for pop in popnames:
        # create set from sample list for this population name
        popset = set(popdict[pop])
        # check if this population set is completely contained within the tsv sample set
        if popset.issubset(allsamples):
            # if so, print something friendly
            print("\t{1}: All {0} sample names match a tsv sample name.\n".format(len(popset), pop))
        else:
            # if not, identify the offending names through set subtraction
            unmatched = list(popset - allsamples)
            # raise error message and print population name and unmatched sample names
            raise ValueError("\n\n\nERROR in population {1}: {0} names do not match a tsv name.\n"
                                 "\tThe names are: {2}.\n\n\n".format(len(unmatched), pop, unmatched))
    print("{}\n\n".format("-"*80))

def prep_outfile(popnames):
    """
    Create output file with header line. Performs a check 
    to ensure file doesn't already exist.
    """
    # create file name based on population names included
    outname = "SNPs_file_{}.txt".format("_".join(popnames))
    # check to ensure file doesn't already exist
    if os.path.exists(outname):
        raise ValueError("\n\n\nWARNING: Output file already exists in the output directory. \nPlease "
                                 "specify a different output directory or remove this file: \n\t'{0}'.\n\n\n".format(outname))
    # create a string of popnames separated by tabs
    cols = "\t".join(popnames)
    # open output file
    with open(outname, 'a') as fh:
        # write SNPS input file header
        fh.write("Ingroup\tOutgroup\tAllele1\t{0}\tAllele2\t{0}\tGene\tPosition\n".format(cols))
        
    return outname

def find_major_minor(snps):
    """
    snps is a list of nucleotide characters, in which there are only
    two distinct bases present. This will determine which is the major
    and minor allele by a count of the bases. It will return the nucleotide
    strings assigned to major and minor. Uses Counter function from collections.
    """
    # create counter object
    c = Counter(snps)
    # assign most common base to major allele
    major = c.most_common(2)[0][0]
    # assign less common base to minor allele
    minor = c.most_common(2)[1][0]
    
    return major, minor

def parse_alleles(locusname, snplist, tsvsamples, popnames, sampledict, outfile, quiet):
    """
    Creates lists of nucleotides for each population and the combined populations.
    The combined populations bases are used to assess if locus is variable, if so 
    the major and minor alleles are determined using the find_major_minor() function. 
    From there the counts of the major and minor alleles are found for each population. 
    The SNPS input file line for this locus is created and written using that info. 
    Returns counts for variable and invariant, where one will be = 1 and the other = 0.

    Arguments:
    locusname - string name of locus
    tsvsamples - list (in order of tsv) = [Sample1, Sample2, Sample3]
    sampledict - dictionary = {Sample1:popname, Sample2:popname}
    popnames - list = [popname1, popname2]
    snplist - list = [A/A, T/A, N/N, etc.]
    outfile - name of output file
    quiet - True or False, affects printing to screen
    """
    # create counters for tracking variant and invariant sites
    invarcount, varcount = int(0), int(0)

    # create a dictionary to populate with bases
    # need a key "!Full" for all samples
    alleledict = {"!Full":[]}
    # also create a key for each popname included
    for p in popnames:
        alleledict[p] = []
    
    # iterate over range covering length of snplist
    for i in range(0, len(snplist)):
        # check if sample name of that column is in the samples desired
        if tsvsamples[i] in sampledict:
            # if snps are not missing data
            if snplist[i] != "N/N":
                # deals with stacks v1 cleaned outputs
                # stacks1 loci format = G/T G/T T T T T
                if "/" not in snplist[i]:
                    # basically add the homozygous base twice (T = T, T) to population specific bases list
                    alleledict[sampledict[tsvsamples[i]]].append(snplist[i])
                    alleledict[sampledict[tsvsamples[i]]].append(snplist[i])
                    # add the homozygous base twice (T = T, T) to the full sampling bases list
                    alleledict["!Full"].append(snplist[i])
                    alleledict["!Full"].append(snplist[i])
                # usage with stacks v2 and heterozygous from stacks1
                # stacks2 loci format = C/C C/G C/G G/G -
                else:
                    # add both bases to the population specific bases list
                    alleledict[sampledict[tsvsamples[i]]].extend(snplist[i].split('/'))
                    # add both bases to the full sampling bases list
                    alleledict["!Full"].extend(snplist[i].split('/'))
                    
    # check if only one base type present
    if len(set(alleledict["!Full"])) == 1:
        if not quiet:
            print("Locus {}: invariant!".format(locusname))
        # add 1 to invariant count
        invarcount += 1
    # else if more than one base type present (should only ever be two)
    else:
        # identify major and minor alleles from combined population sampling
        major, minor = find_major_minor(alleledict["!Full"])
        #print major, minor
        if not quiet:
            print("Locus {}: major = {}, minor = {}".format(locusname, major, minor))

        # initiate empty lists for allele counts for all pops
        # list will contain counts in order of the populations in popnames
        majorcounts, minorcounts = [], []
        # iterate over popnames
        for p in popnames:
            #TEST: print(Counter(alleledict[p]))
            # add major allele count of this population to list
            majorcounts.append(str(alleledict[p].count(major)))
            # add minor allele count of this population to list
            minorcounts.append(str(alleledict[p].count(minor)))
            
        # create output line entry for this locus
        write_row = "-{0}-\t---\t{0}\t{1}\t{2}\t{3}\t{4}\t15\n".format(major, "\t".join(majorcounts), minor,
                                                                            "\t".join(minorcounts), locusname)
        # open output file and write this line
        with open(outfile, 'a') as fh:
            fh.write(write_row)
            
        # add 1 to variant count
        varcount += 1

    return invarcount, varcount
    
def tsv_to_dict(f, popnames, popdict, sampledict, tsvsamples, outfile, quiet):
    """
    Parses tsv file line by line and processes.
    Output is written by parse_alleles function along the way.

    Arguments:
    f - Name (not full path) of the haplotypes.tsv file.
    """
    print("\n\n{}".format("-"*80))
    print("\n\nParsing tsv file to obtain allele counts.\n".format(len(popdict)))

    # create counters to track invariant loci and variable loci
    snps, blanks = int(0), int(0)
    
    # open file
    with open(f, 'r') as fh:
        # skip header - '# Catalog Locus ID', 'Cnt', 'Sample 1', Sample 2'...
        next(fh)
        # iterate over file lines
        for line in fh:
            # Get all snp sites from this line
            if line.strip().split('\t')[2:]:
                # To make things easier downstream, all "-" characters will be transformed to "N/N"
                snplist = [i.replace("-", "N/N") for i in line.strip().split('\t')[2:]]
                # get name of locus from column 1
                locusname = line.split('\t')[0]
                # run parse_alleles for this locus, return counts
                invarcount, varcount = parse_alleles(locusname, snplist, tsvsamples, popnames, sampledict, outfile, quiet)
                # update counters 
                snps += varcount
                blanks += invarcount                
    print("\tDone.\n\n{}\n\n".format("-"*80))
    
    print("\n\n{}".format("-"*80))
    print("\n\nFound {:,} invariant loci with sampling scheme.\n".format(blanks))
    print("Wrote {:,} loci to file: {}.\n\n".format(snps, outfile))
    print("{}\n\n".format("-"*80))
                
def main():
    args = get_args()
    tb = datetime.now()

    # if getnames flag is used, just write sample names to file
    if args.getnames:
         write_sample_names(args.input, args.outdir)

    elif args.checkpopfile:
        popnames, popdict, sampledict = find_pops(args.popfile)
         
    elif args.popfile:
        os.chdir(args.outdir)
        popnames, popdict, sampledict = find_pops(args.popfile)
        tsvsamples = get_sample_names(args.input)
        compare_sampling(popnames, popdict, tsvsamples)
        outfile = prep_outfile(popnames)
        tsv_to_dict(args.input, popnames, popdict, sampledict, tsvsamples, outfile, args.quiet)

    tf = datetime.now()
    print("\n\n{}".format("="*80))
    print("\nTotal elapsed time: {} (H:M:S)\n".format(tf - tb))
    print("{}\n\n".format("="*80))
    
if __name__ == '__main__':
    main()
