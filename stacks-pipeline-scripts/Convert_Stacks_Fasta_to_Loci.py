import os
import argparse
import numpy as np
from datetime import datetime

def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""---------------------------------------------------------------------------
    Convert_Stacks_Fasta_to_Loci: A tool for converting a single populations.samples.fa
    file, which contains the full phased sequences for each sample for each locus. This 
    script will write a fasta file for each locus found. For each sample, the -a flag 
    chooses whether both alleles are written per sample, the first allele only, a 
    randomly selected allele, or a consensus sequence from the two alleles. The -f 
    flag specifies the input fasta file, and the -o flag specifies a directory to write 
    all the output files. 
	---------------------------------------------------------------------------""")
    parser.add_argument("-f", "--fasta",
                            required=True,
                            help="REQUIRED: The full path to a Stacks samples fasta file.")
    
    parser.add_argument("-o", "--outdir",
                            required=True,
                            help="REQUIRED: The full path to an existing directory "
                            "to write output files.")
    
    parser.add_argument("-a", "--alleles",
                            required=True,
                            choices=["both","first", "random", "consensus"],
                            help="REQUIRED: Specify whether to write both alleles of a "
                            "sample to the output file, the first allele, or select one randomly.")
    
    parser.add_argument("--variable",
                            required=False,
                            action='store_true',
                            help="OPTIONAL: If used, retains loci that contain "
                            "at least one variable site.")
    
    return parser.parse_args()

def fasta_dict(f):
    """
    Function to convert a fasta file into
    dictionary structure with custom key name
    and sequence as value.
    
    For example, the following fasta description lines:
    
    >CLocus_6926_Sample_1_Locus_6926_Allele_0 [oc_780LG]
    >CLocus_6926_Sample_1_Locus_6926_Allele_1 [oc_780LG]
    >CLocus_24360_Sample_35_Locus_24360_Allele_0 [oc_GFMJ982]
    >CLocus_24360_Sample_35_Locus_24360_Allele_1 [oc_GFMJ982]

    are converted into:
    
    >Locus-6926|oc_780LG
    >Locus-6926|oc_780LG
    >Locus-24360|oc_GFMJ982
    >Locus-24360|oc_GFMJ982
    """
    
    print("\n\nProcessing {}...".format(f))
    #initiate empty dictionary to store fasta contents
    fdict = {}
    #initiate empy set to store loci names
    loci = set()
    
    #put cleaned lines from fasta file into list
    with open(f, 'r') as fh:
        next(fh)
        lines = [l.strip() for l in fh if l.strip()]
    #iterate over lines and create fasta dictionary
    for line in lines:
        if line.startswith(">"):
            #get locus name from description line
            locus = "Locus-{}".format(line.split('_')[1])
            #add locus name to set
            loci.add(locus)
            if locus not in fdict:
                fdict[locus] = {}
    print("\n\nFound {:,} loci.".format(len(loci)))


    #iterate over lines and create fasta dictionary
    for line in lines:
        if line.startswith(">"):
            #get locus name from description line
            locus = "Locus-{}".format(line.split('_')[1])
            #get sample name from description line
            sample = line.split('[')[-1].replace("]","")
            #create new description line that fits the general structure:
            
            #sample
            #check if key already in dict
            if sample not in fdict[locus]:
                #if not, add key with empty list as val
                fdict[locus][sample] = []
            #if key present, skip it
            else:
                pass
            
        #add alleles to val list for key
        else:
            fdict[locus][sample].append(line.upper())

    return fdict, loci

def make_consensus(v):
    """
    Function to create a consensus sequence from a 
    list of two sequences of equal length. Returns
    consensus sequence. 
    """
    con_dict = {"A-A":"A",
                    "A-C":"M",
                    "A-G":"R",
                    "A-T":"W",
                    "A-N":"A",
                    "C-A":"M",
                    "C-C":"C",
                    "C-G":"S",
                    "C-T":"Y",
                    "C-N":"C",
                    "G-A":"R",
                    "G-C":"S",
                    "G-G":"G",
                    "G-T":"K",
                    "G-N":"G",
                    "T-A":"W",
                    "T-C":"Y",
                    "T-G":"K",
                    "T-T":"T",
                    "T-N":"T",
                    "N-A":"A",
                    "N-C":"C",
                    "N-G":"G",
                    "N-T":"T",
                    "N-N":"N"}
    # intiate empty string
    consensus = ""
    # iterate over length of alleles
    for i in range(0, len(v[0])):
        # construct dict key from allele positions
        # concatenate val to con string
        consensus += con_dict["{0}-{1}".format(list(v[0])[i], list(v[1])[i])]

    return consensus

def write_loci_fastas(fdict, loci, alleles, outdir):
    """
    Function to write output fastas, using alleles option
    selected by the user. 
    """
    print("\n\nWriting output fasta files:")
    # move to output directory
    os.chdir(outdir)
    # iterate over loci names
    for l in sorted(loci):
        # iterate over key, val pairs in dict for locus
        for k, v in sorted(fdict[l].items()):
            # create output file name
            outname = "{}.fasta".format(l)
            # open output file
            with open(outname, 'a') as fh:
                # write information, based on allele option
                if alleles == "both":
                    fh.write(">{0}_Allele0\n{1}\n>{0}_Allele1\n{2}\n".format(k, v[0], v[1]))
                elif alleles == "first":
                    fh.write(">{}\n{}\n".format(k, v[0]))
                elif alleles == "random":
                    i = int(np.random.choice(2, 1))
                    fh.write(">{}\n{}\n".format(k, v[i]))
                elif alleles == "consensus":
                    # make consensus seq first
                    consensus = make_consensus(v)
                    fh.write(">{}\n{}\n".format(k, consensus))
        print("\tWrote {}...".format(outname))

def filter_invariant(outdir):
    # move to output directory
    os.chdir(outdir)
    loci = [f for f in os.listdir('.') if f.endswith('.fasta')]
    print("\n\nInspecting loci for presence of variable sites:\n\n\tFound {:,} loci.".format(len(loci)))
    removed = int(0)
    for fasta in loci:
        with open(fasta, 'r') as fh:
            seqs = set([line.strip() for line in fh if not line.startswith('>') and line.strip()])
            if len(seqs) == 1:
                removed += 1
                os.remove(fasta)
    print("\tRemoved {:,} invariant loci.".format(removed))
    print("\tKept {:,} variable loci.\n\n".format(len(loci)-removed))
    
def main():
    args = get_args()
    tb = datetime.now()
    
    fdict, loci = fasta_dict(args.fasta)
    write_loci_fastas(fdict, loci, args.alleles, args.outdir)
    if args.variable:
        filter_invariant(args.outdir)
    
    tf = datetime.now()
    te = tf - tb
    print("\nFinished.\nElapsed time: {} (H:M:S)\n".format(te))

if __name__ == '__main__':
    main()
