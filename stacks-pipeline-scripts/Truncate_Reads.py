
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
    Truncate - A script to trim reads to a desired length (by removing tail bases). The input 
    directory specified (-i) should contain all uncompressed fastq files. Reads are trimmed 
    using fastx_trimmer, based on the number supplied (-b). For example, if the final desired length
    is 38 bases, enter -b 38. After, all truncated fastq files are moved to the specified 
    output directory (-o). Fastq files must end with '.fq' to be recognized.
    DEPENDENCIES: fastx_trimmer. Must be in path.
    ---------------------------------------------------------------------------""")
    
    parser.add_argument("-i", "--indir",
                            required=True,
                            help="REQUIRED: The full path to a directory which contains "
                            "the uncompressed fastq files.")
    
    parser.add_argument("-o", "--outdir",
                            required=True,
                            help="REQUIRED: The full path to an existing directory "
                            "to write the trimmed fastq files.")
    
    parser.add_argument("-b", "--bases",
                            required=True,
                            type=int,
                            help="REQUIRED: Supply an integer representing the number of "
                            "the LAST base to keep in the fq files. This value represents the "
                            "final length of the trimmed reads.")
    
    return parser.parse_args()


def trim(indir, bases, outdir):
    """
    Function to trim the demultiplexed fastq files
    using fastx_trimmer.

    Arguments:
    indir: Full path to a directory containing the fq files.
    bases: An integer representing the position of the last base 
           to keep.
    """    
    # set working dir to subdir
    os.chdir(indir)
    
    # get fastq files
    fq_list = sorted([f for f in os.listdir('.') if f.endswith('.fq')])

    # if fastq files were found
    if fq_list:
        print("\n\nTrimming reads to {0} bases total in {1} fastq files:\n".format(bases, len(fq_list)))
        
        # iterate over fastq files
        for f in fq_list:
            b = datetime.now()
            print("\tTrimming {}...".format(f))
            
            # split file by period and take first element (prefix)
            outname = "{}.truncated.fq".format(f.split('.')[0])
            # system call for fastx_trimmer
            call_str = ("fastx_trimmer -Q33 -l {0} -i {1} -o {2}"
                                   .format(bases, f, outname))
            # use subprocess to execute system call using shell
            proc = sp.call(call_str, shell=True)
            # move trimmed fastq to output dir
            shutil.move(outname, outdir)
            
            f = datetime.now()
            e = f-b
            print("\t\tElapsed time: {0} (H:M:S)\n".format(e))
        print("---------------------------------------------------\n\n")

    else:
        print("\n\nERROR: No demultiplexed fastq files found after process_radtags!\n\n")
            
def main():
    args = get_args()
    tb = datetime.now()
    
    trim(args.indir, args.bases, args.outdir)
    
    tf = datetime.now()
    te = tf - tb
    print("\n\n--------------------------------------------------------------------------------------")
    print("\nTotal elapsed time: {0} (H:M:S)\n".format(te))
    print("--------------------------------------------------------------------------------------\n\n")    

if __name__ == '__main__':
    main()
