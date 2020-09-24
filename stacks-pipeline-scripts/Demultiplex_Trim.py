
import argparse
import os
import subprocess as sp
import shutil
from datetime import datetime

def get_args():
    """
    Get arguments from command line.
    """
    parser = argparse.ArgumentParser(
            description="""---------------------------------------------------------------------------
    Demultiplex_Trim - A script to demultiplex and trim ddRADseq sequencing files. This is designed to 
    run on ddRADseq libraries prepared with 'sbfI', and single-end (SE) reads of any length. 
    The input directory specified (-i) should contain several subdirectories. Each subdirectory should 
    contain a sequencing file (gzipped fastq, ending with 'fastq.gz') and a barcode file ('barcode.txt'). 
    The gzipped fastq is demultiplexed using the process_radtags module in Stacks (v2.4) with the 
    associated barcode file. The RAD cutsites of the reads in the demultiplexed fastq files are 
    then trimmed using fastx_trimmer, based on the number supplied (-b). After, all trimmed fastq files 
    are moved to the output directory (-o). If the libraries were prepared using a unique molecular 
    identifier (UMI), this should be trimmed before demultiplexing by using the -u flag. The gzipped 
    fastq files are trimmed based on the number supplied (-u), then demultiplexed and trimmed as described 
    above. 
    DEPENDENCIES: fastx_trimmer, process_radtags (from Stacks v2.4). Both must be in path.
    ---------------------------------------------------------------------------""")
    
    parser.add_argument("-i", "--indir",
                            required=True,
                            help="REQUIRED: The full path to a directory which contains "
                            "the subdirectories containing the fastq.gz sequencing files.")
    
    parser.add_argument("-o", "--outdir",
                            required=True,
                            help="REQUIRED: The full path to an existing directory "
                            "to write the demultiplexed and trimmed fastq files.")
    
    parser.add_argument("-b", "--bases",
                            required=True,
                            type=int,
                            help="REQUIRED: Supply an integer representing the number of "
                            "bases to remove from the demultiplexed fq files, in order "
                            "to remove the RAD cutsites from the processed reads. If using "
                            "'sbfI', it is 6 bp long, so enter 6.")

    parser.add_argument("-u", "--umi",
                            required=False,
                            type=int,
                            default=None,
                            help="Optional: Supply an integer representing the number of "
                            "the first base to keep in the fastq.gz files, in order to remove "
                            "the UMI (unique molecular identifier) from the reads. For "
                            "example, if the UMI is 8 bp, enter 9.")
    
    return parser.parse_args()

def get_directories(indir):
    """
    Function to find all directories contained in the
    input directory.

    Arguments
    indir: Full path to an input directory containing the subdirectories.
    """
    print("\n\n======================================================================================")
    print("Locating subdirectories in main directory.")
    print("======================================================================================\n")
    # change working directory to input directory
    os.chdir(indir)
    # list comprehension to find all directories, get full paths
    paths = sorted([os.path.abspath(f) for f in os.listdir('.') if os.path.isdir(f)])
    
    # check if paths actually found
    if paths:
        # print number of subdirs found
        print("\nFound {} subdirectories to examine:".format(len(paths)))
        for p in paths:
            print("\t{}".format(p))
    else:
        raise ValueError("\n\n\nERROR: No subdirectories were found in the directory: {}".format(indir))
            
    return paths

def trim_umi(indir, paths, umi):
    """
    Function to trim the UMI from all gzipped
    fastq files present in the subdirectories
    contained in the input directory.

    Arguments
    indir: Full path to an input directory containing the subdirectories.
    paths: List of full paths to all subdirectories in 'indir'.
    umi: An integer representing the length of the UMI, and the 
          number of front bases to trim.
    """
    print("\n\n======================================================================================")
    print("Trimming UMI from all sequencing files.")
    print("======================================================================================\n")

    # iterate over subdir paths
    for p in paths:
        # set working dir to subdir
        os.chdir(p)
        # get input gzip fastq via list comprehension
        seqf = sorted([f for f in os.listdir('.') if f.endswith(".fastq.gz")])
        
        # if gzip fastq file(s) found in subdir
        if seqf:
            for f in seqf:
                b = datetime.now()
                print("\nProcessing {}:\n".format(f))
                # split file by period and take first element (prefix)
                prefix = f.split('.')[0]
                
                # system call for fastx_trimmer
                call_str = ("seqtk trimfq -b {0} {1} > {2}_UMI_trimmed.fastq.gz"
                                   .format(umi, f, prefix))
                print(call_str)
                # use subprocess to execute system call using shell
                proc = sp.call(call_str, shell=True)

                f = datetime.now()
                print("\nElapsed time: {0} (H:M:S)\n".format(f-b))
                print("---------------------------------------------------")
        else:
            print("\n\nERROR: No sequencing file ending with '.fastq.gz' was found in directory: {}\n\n"
                      .format(p.split('/')[-1]))
            
    # set working directory to indir
    os.chdir(indir)

def trim(fq_list, bases, outdir):
    """
    Function to trim the demultiplexed fastq files
    using fastx_trimmer.

    Arguments:
    fq_list: List of fq filenames.
    bases: An integer representing the first base position of the reads to keep.
    """        
    print("\n\nTrimming first {0} bases in {1} demultiplexed fastq files:\n".format(bases-1, len(fq_list)))

    # iterate over fastq files
    for fq in fq_list:
        tb = datetime.now()
        print("\tTrimming {}...".format(fq))
        
        # create output file name
        outname = "{}.trim.fq".format(fq.split('.')[0])
        
        # system call for fastx_trimmer
        call_str = ("seqtk trimfq -b {0} {1} > {2}"
                               .format(bases, fq, outname))
            
        # use subprocess to execute system call using shell
        print(call_str)
        proc = sp.call(call_str, shell=True)
        
        # move output file to output dir
        shutil.move(outname, outdir)
        
        # remove untrimmed fq
        os.remove(fq)
        
        tf = datetime.now()
        te = tf-tb
        print("\t\tElapsed time: {0} (H:M:S)\n".format(te))    
    
def demultiplex(indir, paths, umi, bases, outdir):
    """
    Function to demultiplex all gzipped fastq files
    using the 'process_radtags' module of Stacks v2.0+.

    Arguments:
    indir: Full path to an input directory containing the subdirectories.
    paths: List of full paths to all subdirectories in 'indir'.
    umi: An integer representing the length of the UMI. 
          If not provided, default = None.
    bases: An integer representing the first base position of the reads to keep.
    outdir: Full path to an output directory to move the 
            demultiplexed and trimmed fastq files.
    """
    print("\n\n======================================================================================")
    print("Demultiplexing and trimming all sequencing files.")
    print("======================================================================================\n")
    
    # name of barcode file that should be present in each subdir
    BCODE = "barcode.txt"
    # iterate over subdir paths
    for p in paths:
        # set working dir to subdir
        os.chdir(p)
        # get input gzip fastq depending on whether UMI was trimmed or not
        if umi:
            seqf = sorted([f for f in os.listdir('.') if f.endswith("UMI_trimmed.fastq.gz")])
        else:
            seqf = sorted([f for f in os.listdir('.') if f.endswith("fastq.gz")])
            
        # if gzip fastq file(s) found in subdir
        if seqf:
            for gz in seqf:
                tb = datetime.now()
                print("Demultiplexing: {}\n\n".format(gz))

                # system call for process_radtags
                call_str = ("process_radtags -f {0} -i {1} -y {2} -o {3} -b {4} -c -q -r -e {5} --inline_null"
                                   .format(gz, "gzfastq", "fastq", ".", BCODE, 'sbfI'))
                    
                print("{}\n".format(call_str))

                # use subprocess to execute system call using shell
                proc = sp.call(call_str, shell=True)

                # get fastq files
                fq_list = sorted([f for f in os.listdir('.') if f.endswith('.fq')])
                
                # if fastq files were found
                if fq_list:
                    # call function to trim all demultiplexed fastq files found here
                    trim(fq_list, bases, outdir)
                
                else:
                    print("\n\nERROR: No demultiplexed fastq files found after process_radtags!\n\n")
                
                tf = datetime.now()
                te = tf-tb
                print("---------------------------------------------------------------------")
                print("\nElapsed time for {0}: {1} (H:M:S)\n".format(gz, te))
                print("---------------------------------------------------------------------\n\n\n")  
                
        else:
            if umi:
                print("\n\nERROR: No sequencing file ending with 'UMI_trimmed.fastq.gz' was found in directory: {}\n\n"
                          .format(p.split('/')[-1]))
            else:
                print("\n\nERROR: No sequencing file ending with '.fastq.gz' was found in directory: {}\n\n"
                          .format(p.split('/')[-1]))
                
    # set working directory to indir
    os.chdir(indir)

def base_checker(text):
    """
    Simple function to check if all elements in
    a string are contained in a list of characters.
    In this case, whether the string is a valid bardcode.
    Returns True or False.
    
    text: the string to test
    """
    outcome = True
    bases = ["A", "T", "C", "G"]
    for t in text:
        if t not in bases:
            outcome = False
    return outcome
    
def log_combiner(indir, paths):
    """
    Function to combine all log files produced for fastq files
    using the 'process_radtags' module of Stacks v2.0+.

    Arguments:
    indir: Full path to an input directory containing the subdirectories.
    paths: List of full paths to all subdirectories in 'indir'.
    """
    print("\n\n======================================================================================")
    print("Combining log files.")
    print("======================================================================================\n")

    # name of log file that should be present in each subdir
    LOG = "process_radtags.log"
    # final list of contents to write
    loglist = []
    # iterate over subdir paths
    for p in paths:
        print("Getting contents from: /{}/{}".format(p.split('/')[-1], LOG))
        # set working dir to subdir
        os.chdir(p)
        with open(LOG, 'r') as fh:
            # get sublists from lines if they have more than 4 elements after tab split
            contents = [line.strip().split('\t') for line in fh if len(line.split('\t')) > 4]
            
            # keep sublists if element one is a barcode (contains only A, T, C, or G), join by tabs
            filtered = ["\t".join(i) for i in contents if base_checker(i[0]) is True]
            
            # add filtered, tab-joined strings to other list
            for i in filtered:
                loglist.append(i)

    # set working directory to indir
    os.chdir(indir)
    # open output file
    with open("Stacks_ProcessRadtags.log", 'a') as fh:
        # write headers
        fh.write("Barcode\tFilename\tTotal\tNo_RadTag\tLow_Quality\tRetained\n")
        # write main contents
        for i in loglist:
            fh.write("{}\n".format(i))
            
    # set working directory to indir
    os.chdir(indir)
            
def main():
    args = get_args()
    tb = datetime.now()
    
    paths = get_directories(args.indir)
    
    if args.umi:
        trim_umi(args.indir, paths, args.umi)
        
    demultiplex(args.indir, paths, args.umi, args.bases, args.outdir)
    log_combiner(args.indir, paths)
    
    tf = datetime.now()
    te = tf - tb
    print("\n\n--------------------------------------------------------------------------------------")
    print("\nTotal elapsed time: {0} (H:M:S)\n".format(te))
    print("--------------------------------------------------------------------------------------\n\n")    

if __name__ == '__main__':
    main()
