import argparse
import os
import subprocess as sp
import shutil
import numpy as np
from datetime import datetime

def get_args():
    """
    Get arguments from command line.
    """
    parser = argparse.ArgumentParser(
            description="""---------------------------------------------------------------------------
    Run_Stacks - Process RADseq data using Stacks v2.4+. Automates running the full Stacks
    pipeline (ustacks, cstacks, sstacks, tsv2bam, gstacks, and populations) for a set of 
    fastq files. The input fastq files should be located in the output directory where all output 
    files are generated (-o), or they can optionally be located in a separate directory (-i). If located in 
    a different directory (-i), the fastq files are first copied (and renamed) to the output directory 
    (-o) before any steps are run in Stacks. The pipeline can be run with all steps (-a full), only post-ustacks 
    steps (-a post-ustacks), or each of the Stacks modules can be accessed separately. Several key parameters 
    can be set for ustacks (-M, -m), cstacks (-n, --catnum), and populations (--maf, --mac, --mintype, --whitelist, 
    --blacklist). The populations module is run using a set of values for the --mintype selection 
    (implementing -r or -R in populations). These values include 10 - 100, with increments of 10. This 
    script can be used to run ustacks independently for separate sets of fastq files for the same 
    project (different -i, same -o), which requires using the --intid and -l flags. The combined 
    output files can then be analyzed together for all subsequent steps. This strategy can be used 
    to speed up ustacks processing (typically the longest step when run sequentially). See documentation 
    for additional details.
 
    DEPENDENCIES: Stacks v2.4. All modules must be in path.
    ---------------------------------------------------------------------------""")
    
    parser.add_argument("-o", "--outdir",
                            required=True,
                            help="REQUIRED: The full path to an existing directory "
                            "to run the Stacks pipeline steps. If fastq files are not "
                            "present here, their location must be specified with -i.")
    
    parser.add_argument("-a", "--analysis",
                            required=True,
                            choices=["full", "post-ustacks", "ustacks","cstacks", "sstacks",
                                         "tsv2bam", "gstacks","populations"],
                            help="REQUIRED: Specify running the complete Stacks pipeline "
                            "(full), all steps beyond ustacks (post-ustacks), or "
                            "one of the specific modules.")

    parser.add_argument("-i", "--indir",
                            required=False,
                            default=None,
                            help="Optional: The full path to a directory which contains "
                            "the starting fastq files, which are copied to the location "
                            "specified by -o to run Stacks.")
    
    parser.add_argument("-t", "--threads",
                            required=False,
                            type=int,
                            default=1,
                            help="OPTIONAL: Specifies number of threads to use.")
    
    parser.add_argument("-p", "--popmap",
                            required=False,
                            default=None,
                            help="OPTIONAL: Provide the full path to a custom population map file to "
                            "use instead of the default (all auto-assigned to pop 1).")
    
    parser.add_argument("-M", "--maxdist",
                            required=False,
                            type=int,
                            default=2,
                            help="ustacks OPTIONAL: Maximum distance (in nucleotides) allowed "
                            "between stacks. Default = 2.")

    parser.add_argument("-m", "--mindepth",
                            required=False,
                            type=int,
                            default=3,
                            help="ustacks OPTIONAL: Minimum depth of coverage required to "
                            "create a stack. Default = 3.")
    
    parser.add_argument("--intid",
                            required=False,
                            type=int,
                            default=None,
                            help="ustacks OPTIONAL: Specify an integer ID to start at when labeling "
                            "samples (overrides starting at 1). Useful for running ustacks on multiple "
                            "fastq file batches for a single project (different -i, same -o).")
    
    parser.add_argument("-n", "--nmismatch",
                            required=False,
                            type=int,
                            default=2,
                            help="cstacks OPTIONAL: Number of mismatches allowed between sample "
                            "loci when building the catalog. Default = 2.")
    
    parser.add_argument("--catnum",
                            required=False,
                            type=int,
                            default=None,
                            help="cstacks OPTIONAL: Specify the number of samples to build the catalog "
                            "from, rather than use all samples. This many samples will be drawn randomly "
                            "from the population map file.")
    
    parser.add_argument("--maf",
                            required=False,
                            type=float,
                            default=None,
                            help="populations OPTIONAL: Specify a minimum minor allele FREQUENCY "
                            "(0 - 0.5) required to process a nucleotide site at a locus. "
                            "Default = None.")
    
    parser.add_argument("--mac",
                            required=False,
                            type=int,
                            default=None,
                            help="populations OPTIONAL: Specify a minimum minor allele COUNT "
                            "(an integer) required to process a nucleotide site at a locus. "
                            "Setting to 2 eliminates singletons. Default = None.")
    
    parser.add_argument("--mintype",
                            required=False,
                            default="r",
                            choices=["r", "R"],
                            help="populations OPTIONAL: Specify whether the minimum number of samples "
                            "within a population (r) or across populations (R) should be used to "
                            "process loci in populations. The R option should only be used when multiple "
                            "populations are defined in a custom population map using -p. Default = r.")
    
    parser.add_argument("--whitelist",
                            required=False,
                            default=None,
                            help="populations OPTIONAL: Full path to a file containing whitelisted "
                            "markers to be included in the export.")
    
    parser.add_argument("--blacklist",
                            required=False,
                            default=None,
                            help="populations OPTIONAL: Full path to a file containing blacklisted "
                            "markers to be excluded from the export.")
        
    parser.add_argument("-l", "--logstring",
                            required=False,
                            default=None,
                            help="OPTIONAL: Add unique label to end of Stacks pipeline log file"
                            "name. Must use if running ustacks on multiple fastq file batches for a "
                            "single project (different -i, same -o). to prevent log file overwriting.")

    return parser.parse_args()

def pipeline_logger(LOG, argd):
    """
    Function to record argument settings in
    the log file. 

    Arguments:
    LOG: Full path to log file.
    argd: Dictionary of argparse entries.
    """
    with open(LOG, 'a') as fh:
        fh.write("Run executed: {}\n\nStacks Pipeline settings:\n-o (outdir): {}\n-a (analysis): {}\n"
                     "-i (indir): {}\n-t (threads): {}\n-p (popmap): {}\n-M (maxdist): {}\n"
                     "-m (mindepth): {}\n--intid (id int start): {}\n-n (nmismatch): {}\n--catnum: {}\n"
                     "--maf (minor allele freq): {}\n--mac (minor allele count): {}\n"
                     "--mintype (r or R): {}\n--whitelist: {}\n--blacklist: {}\n-l (logstring): {}\n\n"
                     .format(datetime.now(), argd["outdir"], argd["analysis"],
                                 argd["indir"], argd["threads"], argd["popmap"], argd["maxdist"],
                                 argd["mindepth"], argd["intid"], argd["nmismatch"], argd["catnum"],
                                 argd["maf"],  argd["mac"], argd["mintype"], argd["whitelist"],
                                 argd["blacklist"], argd["logstring"]))
    
def write_id_log(dst_fastq, intid):
    """
    Function to record the integer identity assigned
    to all fastq files processed in Stacks. If intid 
    flag was included, assigns starting with this number. 
    If not included, it is auto-assigned using enumerate.

    Arguments:
    dst_fastq: List of fastq filenames.
    intid: An integer to start sample labeling. 
    """
    if intid:
        i = intid
        if not os.path.exists("Identifier.log"):
            with open("Identifier.log", 'a') as fh:
                fh.write("Sample\tIdentifier_code\n")
                for f in dst_fastq:
                    fh.write("{0}\t{1}\n".format(i, f))
                    i=i+1
        else:
            with open("Identifier.log", 'a') as fh:
                for f in dst_fastq:
                    fh.write("{0}\t{1}\n".format(i, f))
                    i=i+1
    else:
        with open("Identifier.log", 'a') as fh:
            fh.write("Sample\tIdentifier_code\n")
            for i, f in enumerate(dst_fastq, start=1):
                fh.write("{0}\t{1}\n".format(i, f))

def copy_fastq(indir, outdir):
    """
    Function to copy (and rename) all fastq files
    from input directory to main directory. Returns
    list (dst_fastq) of new fastq files names found 
    in main directory.

    Arguments:
    indir: Full path to an input directory containing the 
           fastq files.
    outdir: Full path to an output directory to rename and 
            move the fastq files, and write outputs.
    """
    print("\nCopying and renaming all fastq files to output directory:")
    os.chdir(indir)
    src_fastq = sorted([os.path.abspath(f) for f in os.listdir('.') if f.endswith(".fq")])
    dst_fastq = []
    if src_fastq:
        for f in src_fastq:
            print("\t{}...".format(f.split('/')[-1]))
            outname = os.path.join(outdir, "{}.fq".format(f.split('/')[-1].split('.')[0]))
            dst_fastq.append("{}.fq".format(f.split('/')[-1].split('.')[0]))
            if not os.path.exists(outname):
                shutil.copyfile(f, outname)
            else:
                raise ValueError("\n\n\nERROR: File {} already exists! Please remove before running again.\n".format(outname))
        print("\nFinished actions for {} fastq files.".format(len(src_fastq)))   
    else:
        raise ValueError("\n\n\nERROR: No fastq files found!\n\n")

    return dst_fastq

def show_elapsed(tb, tf, text, LOG, outline="thick"):
    """
    Function to calculate elapsed time, print
    message, then record to log file.

    Arguments:
    tb: start time (from datetime.now())
    tf: end time (from datetime.now())
    text: string to insert into message (method)
    LOG: Full path to log file.
    """
    te = tf - tb
    tprint = "Elapsed time for {0}: {1} (H:M:S)".format(text, te)
    
    if outline == "thick":
        print("\n\n{}".format("="*80))
        print("\n{}\n".format(tprint))
        print("{}\n\n".format("="*80))
        
    elif outline == "thin":
        print("\n\n{}".format("="*70))
        print("\n{}\n".format(tprint))
        print("{}\n\n\n".format("="*70))
        
    with open(LOG, 'a') as fh:
        fh.write("{0}: Completed {1}: Elapsed time: {2}  (H:M:S)\n\n".format( tf, text, te))
    
            
def run_ustacks(indir, outdir, threads, maxdist, mindepth, intid, LOG):
    """
    Function to process all trimmed fastq files
    using the 'ustacks' module of Stacks v2.0+.

    Arguments:
    indir: Full path to an input directory containing the 
           fastq files.
    outdir: Full path to an output directory to rename and 
            move the fastq files, and write outputs.
    threads: Integer. Specifies number of threads to use.
    maxdist: Maximum distance (in nucleotides) allowed 
             between stacks. Default = 2.
    mindepth: Minimum depth of coverage required to create 
              a stack. Default = 3.
    intid: An integer to start sample labeling. 
    LOG: Full path to log file.
    """
    print("\n\n{}".format("="*80))
    print("\nRunning ustacks.\n")
    print("{}\n".format("="*80))
    tb1 = datetime.now()

    if indir:
        dst_fastq = copy_fastq(indir, outdir)
        
    else:
        os.chdir(outdir)
        dst_fastq = sorted([f for f in os.listdir('.') if f.endswith(".fq")])
        
    if dst_fastq:
        os.chdir(outdir)
        write_id_log(dst_fastq, intid)
        print("\n\n\nRunning ustacks for {} fastq files.\n\n".format(len(dst_fastq)))
        if intid:
            i = intid
            for f in dst_fastq:
                tb2 = datetime.now()
                print("Processing {}:\n\n".format(f))
                # system call for ustacks
                call_str = ("ustacks -f {0} -i {1} -o {2} -M {3} -m {4} -p {5}"
                                .format(f, i, outdir, maxdist, mindepth, threads))
                print("{}\n".format(call_str))
                with open(LOG, 'a') as fh:
                    fh.write("{0}: Executed ustacks: {1}\n".format(datetime.now(), call_str))
                # use subprocess to execute system call using shell
                proc = sp.call(call_str, shell=True)
                tf2 = datetime.now()
                show_elapsed(tb2, tf2, "ustacks {}".format(f), LOG, outline="thin")
                i=i+1
        else:
            for i, f in enumerate(dst_fastq, start=1):
                tb2 = datetime.now()
                print("Processing {}:\n\n".format(f))
                # system call for ustacks
                call_str = ("ustacks -f {0} -i {1} -o {2} -M {3} -m {4} -p {5}"
                                .format(f, i, outdir, maxdist, mindepth, threads))
                print("{}\n".format(call_str))
                with open(LOG, 'a') as fh:
                    fh.write("{0}: Executed ustacks: {1}\n".format(datetime.now(), call_str))
                # use subprocess to execute system call using shell
                proc = sp.call(call_str, shell=True)
                tf2 = datetime.now()
                show_elapsed(tb2, tf2, "ustacks {}".format(f), LOG, outline="thin")
            
    else:
        raise ValueError("\n\n\nERROR: No fastq files found to run ustacks!\n\n")

    tf1 = datetime.now()
    show_elapsed(tb1, tf1, "all ustacks processing", LOG)
    
def generate_pop_map(outdir):
    """
    Function to create a population map to execute steps.
    Auto assigns all samples to population 1.

    Arguments:
    outdir: Full path to an output directory with the renamed 
            fastq files.
    """
    print("\n\n{}".format("="*80))
    print("\nGenerating population map.\n")
    print("{}\n".format("="*80))
    print("Adding samples to file:")
    # change to output directory
    os.chdir(outdir)
    # name output file
    outname = "population_map.txt"
    # remove if present, before writing
    if os.path.exists(outname):
        os.remove(outname)
    # open output file
    with open(outname, 'a') as fh:
        # iterate over file prefixes and write to file
        for i in [f.split('.')[0] for f in os.listdir('.') if f.endswith('.fq')]:
            print("\t{}".format(i))
            fh.write("{0}\t1\n".format(i))
    # send back full path to this file
    return os.path.abspath(outname)
    
def run_cstacks(outdir, threads, pop_map, nmismatch, catnum, LOG):
    """
    Function to process all fastq & output files 
    using the 'cstacks' module of Stacks v2.0+.

    Arguments:
    outdir: Full path to an output directory with the renamed fastq 
            files, containing the outputs (#.alleles.tsv, #.snps.tsv, #.tags.tsv).
    threads: Integer. Specifies number of threads to use.
    pop_map: Full path to the population map file.
    nmismatch: Integer. Number of mismatches allowed between sample
               loci when building the catalog. Default = 2.
    catnum: Integer. Number of samples to use for building the catalog. This
            number will be drawn randomly from the population map file.
    LOG: Full path to log file.
    """
    print("\n\n{}".format("="*80))
    print("\nRunning cstacks.\n")
    print("{}\n".format("="*80))
    tb = datetime.now()
    # change to output directory
    os.chdir(outdir)
    
    # check if outputs already present, raise error
    if os.path.exists("catalog.tags.tsv"):
        raise ValueError("\n\n\nERROR: Loci catalog files already exists! Please remove before running.\n")

    # if a subset of samples are to be used to build the catalog
    if catnum:
        # open the population map file
        with open(pop_map, 'r') as fh:
            # get sample names via list comprehension
            samples = [l.split()[0] for l in fh]
        # create a list of X randomly generated numbers from length of sample names,
        # sampling without replacement
        indices = np.random.choice(len(samples), catnum, replace=False)
        # initiate empty list to store name flags
        cat_samples = []
        # enumerate sample name list
        for i, j in enumerate(samples):
            # see if index is in the random numbers 
            if i in indices:
                # if so, append the cstacks flag with this sample name
                cat_samples.append("-s {}".format(j))
        # create single string for use in cstacks call
        samples_str = " ".join(cat_samples)
        
        # system call for cstacks
        call_str = ("cstacks {0} -p {1} -n {2}".format(samples_str, threads, nmismatch))

    else:
        # system call for cstacks
        call_str = ("cstacks -P {0} -M {1} -p {2} -n {3}".format(outdir, pop_map, threads, nmismatch))
    
    print("{}\n".format(call_str))
    with open(LOG, 'a') as fh:
        fh.write("{0}: Executed cstacks: {1}\n".format(datetime.now(), call_str))
    
    # use subprocess to execute system call using shell
    proc = sp.call(call_str, shell=True)
    
    tf = datetime.now()
    show_elapsed(tb, tf, "cstacks", LOG)
    
def run_sstacks(outdir, threads, pop_map, LOG):
    """
    Function to process all fastq & output files 
    using the 'sstacks' module of Stacks v2.0+.

    Arguments:
    outdir: Full path to an output directory with the renamed fastq 
            files, containing the outputs (#.alleles.tsv, #.snps.tsv, #.tags.tsv).
    threads: Integer. Specifies number of threads to use.
    pop_map: Full path to the population map file.
    LOG: Full path to log file.
    """
    print("\n\n{}".format("="*80))
    print("\nRunning sstacks.\n")
    print("{}\n".format("="*80))
    tb = datetime.now()
    # change to output directory
    os.chdir(outdir)
    
    # check if catalog present. if not, raise error
    if not [f for f in os.listdir('.') if f.endswith(("atalog.snps.tsv", "atalog.tags.tsv", "atalog.alleles.tsv"))]:
        raise ValueError("\n\n\nERROR: No catalog files present. Quitting.\n")
    
    # check if outputs already present, raise error
    if [f for f in os.listdir('.') if f.endswith((".matches.tsv"))]:
        raise ValueError("\n\n\nERROR: Loci match files (.matches.tsv) already exist! Please remove before running.\n")
    
    # system call for sstacks
    call_str = ("sstacks -P {0} -M {1} -p {2}".format(outdir, pop_map, threads))
    print("{}\n".format(call_str))
    with open(LOG, 'a') as fh:
        fh.write("{0}: Executed sstacks: {1}\n".format(datetime.now(), call_str))
    # use subprocess to execute system call using shell
    proc = sp.call(call_str, shell=True)
    
    tf = datetime.now()
    show_elapsed(tb, tf, "sstacks", LOG)
    
def run_tsv2bam(outdir, threads, pop_map, LOG):
    """
    Function to process all fastq & output files 
    using the 'tsv2bam' module of Stacks v2.0+.

    Arguments:
    outdir: Full path to an output directory with the renamed fastq 
            files, containing the outputs (#.alleles.tsv, #.snps.tsv, #.tags.tsv).
    threads: Integer. Specifies number of threads to use.
    pop_map: Full path to the population map file.
    LOG: Full path to log file.
    """
    print("\n\n{}".format("="*80))
    print("\nRunning tsv2bam.\n")
    print("{}\n".format("="*80))
    tb = datetime.now()
    # change to output directory
    os.chdir(outdir)
    
    # check if sstacks output present, if not raise error
    if not [f for f in os.listdir('.') if f.endswith((".matches.tsv"))]:
        raise ValueError("\n\n\nERROR: No loci match files (.matches.tsv) found. Quitting.\n")

    # check if outputs already present, raise error
    if [f for f in os.listdir('.') if f.endswith((".matches.bam"))]:
        raise ValueError("\n\n\nERROR: tsv2bam match files (.matches.bam) already exist! Please remove before running.\n")
    
    # system call for tsv2bam
    call_str = ("tsv2bam -P {0} -M {1} -t {2}".format(outdir, pop_map, threads))
    print("{}\n".format(call_str))
    with open(LOG, 'a') as fh:
        fh.write("{0}: Executed tsv2bam: {1}\n".format(datetime.now(), call_str))
    # use subprocess to execute system call using shell
    proc = sp.call(call_str, shell=True)
    
    tf = datetime.now()
    show_elapsed(tb, tf, "tsv2bam", LOG)

def run_gstacks(outdir, threads, pop_map, LOG):
    """
    Function to process all fastq & output files 
    using the 'gstacks' module of Stacks v2.0+.

    Arguments:
    outdir: Full path to an output directory with the renamed fastq 
            files, containing the outputs (#.alleles.tsv, #.snps.tsv, #.tags.tsv).
    threads: Integer. Specifies number of threads to use.
    pop_map: Full path to the population map file.
    LOG: Full path to log file.
    """
    print("\n\n{}".format("="*80))
    print("\nRunning gstacks.\n")
    print("{}\n".format("="*80))
    tb = datetime.now()
    # change to output directory
    os.chdir(outdir)
    
    # check if outputs already present, raise error
    if not [f for f in os.listdir('.') if f.endswith((".matches.bam"))]:
        raise ValueError("\n\n\nERROR: No tsv2bam match files (.matches.bam) found. Quitting.\n")
    
    # check if outputs already present, raise error
    if [f for f in os.listdir('.') if f.endswith((".fa.gz", ".calls", ".log.distribs"))]:
        raise ValueError("\n\n\nERROR: gstacks output files (catalog.fa.gz and others) already exist! Please remove before running.\n")
    
    # system call for gstacks
    call_str = ("gstacks -P {0} -M {1} -t {2}".format(outdir, pop_map, threads))
    print("{}\n".format(call_str))
    with open(LOG, 'a') as fh:
        fh.write("{0}: Executed gstacks: {1}\n".format(datetime.now(), call_str))
    # use subprocess to execute system call using shell
    proc = sp.call(call_str, shell=True)
    
    tf = datetime.now()
    show_elapsed(tb, tf, "gstacks", LOG)

def run_populations(outdir, threads, pop_map, maf, mac, mintype, whitelist, blacklist, LOG):
    """
    Function to process all fastq & output files 
    using the 'populations' module of Stacks v2.0+.

    Arguments:
    outdir: Full path to an output directory with the renamed fastq 
            files, containing the outputs (#.alleles.tsv, #.snps.tsv, #.tags.tsv).
    threads: Integer. Specifies number of threads to use.
    pop_map: Full path to the population map file.
    LOG: Full path to log file.
    """
    print("\n\n{}".format("="*80))
    print("\nRunning populations.\n")
    print("{}\n".format("="*80))
    tb = datetime.now()
    
    # check if gstacks outputs present, if not raise error
    if not [f for f in os.listdir('.') if f.endswith((".fa.gz", ".calls", ".log.distribs"))]:
        raise ValueError("\n\n\nERROR: gstacks output files (catalog.fa.gz and others) not found. Quitting.\n")

    # set r (or R) values 
    rvals = ["10", "20", "30", "40", "50", "60", "70", "80", "90", "100"]
    # iterate over r vals to run populations
    for r in rvals:
        tb1 = datetime.now()
        print("\n\n{}".format("="*70))
        print("Running population module with {0} value {1}".format(mintype, r))
        print("{}\n\n".format("="*70))
        
        rdir = "Populations_{0}{1}".format(mintype, r)
        os.chdir(outdir)
        if not os.path.exists(rdir):
            os.mkdir(rdir)
        else:
            raise ValueError("\n\n\nERROR: populations output directories already exist! Please remove before running.\n")
        
        # system call for populations
        if maf:
            mafstr = "--min-maf {} ".format(maf)
        else:
            mafstr = ""
        if mac:
            macstr = "--min-mac {} ".format(mac)
        else:
            macstr = ""
        if whitelist:
            whitestr = "-W {} ".format(whitelist)
        else:
            whitestr = ""
        if blacklist:
            blackstr = "-B {} ".format(blacklist)
        else:
            blackstr = ""
        call_str = ("populations -P {0} -O {1} -M {2} -t {3} -{4} {5} {6}{7}{8}{9}--fasta-samples"
                        .format(outdir, os.path.abspath(rdir), pop_map, threads, mintype, r, mafstr, macstr, whitestr, blackstr))
        print("{}\n".format(call_str))
        with open(LOG, 'a') as fh:
            fh.write("{0}: Executed populations -{3} {2}: {1}\n".format(datetime.now(), call_str, r, mintype))
        # use subprocess to execute system call using shell
        proc = sp.call(call_str, shell=True)
        
        # relabel output files using args.mintype and r value
        os.chdir(rdir)
        print("\n\nRenaming output files...\n\n")
        outputs = [f for f in os.listdir('.') if f.startswith("populations.")]
        for o in outputs:
            shutil.copyfile(o, os.path.join(outdir, rdir, o.replace("populations.", "populations_{0}{1}.".format(mintype, r))))
            os.remove(o)
            
        tf1 = datetime.now()
        show_elapsed(tb1, tf1, "populations -{0} {1}".format(mintype, r), LOG, outline="thin")
    
    tf = datetime.now()
    show_elapsed(tb, tf, "populations", LOG)

def main():
    args = get_args()
    argd = vars(args)
    tb = datetime.now()
    
    if args.logstring:
        LOG = os.path.join(args.outdir, "Stacks_pipeline_{}_{}.log".format(args.analysis, args.logstring))
    else:
        LOG = os.path.join(args.outdir, "Stacks_pipeline.log")
    pipeline_logger(LOG, argd)
    print("\n\nStacks Pipeline settings:\n\n-o (outdir): {}\n-a (analysis): {}\n"
              "-i (indir): {}\n-t (threads): {}\n-p (popmap): {}\n-M (maxdist): {}\n"
              "-m (mindepth): {}\n--intid (id int start): {}\n-n (nmismatch): {}\n--catnum: {}\n"
              "--maf (minor allele freq): {}\n--mac (minor allele count): {}\n"
              "--mintype (r or R): {}\n--whitelist: {}\n--blacklist: {}\n-l (logstring): {}\n\n"
              .format(argd["outdir"], argd["analysis"],
                          argd["indir"], argd["threads"], argd["popmap"], argd["maxdist"],
                          argd["mindepth"], argd["intid"], argd["nmismatch"], argd["catnum"],
                          argd["maf"],  argd["mac"],
                          argd["mintype"], argd["whitelist"], argd["blacklist"], argd["logstring"]))
        
    if args.analysis == "full":
        run_ustacks(args.indir, args.outdir, args.threads, args.maxdist, args.mindepth, args.intid, LOG)
        os.chdir(args.outdir)
        if [f for f in os.listdir('.') if f.endswith((".snps.tsv", ".tags.tsv", ".alleles.tsv"))]:
            if args.popmap:
                pop_map = args.popmap
            else:
                pop_map = generate_pop_map(args.outdir)
            run_cstacks(args.outdir, args.threads, pop_map, args.nmismatch, args.catnum, LOG)
            run_sstacks(args.outdir, args.threads, pop_map, LOG)
            run_tsv2bam(args.outdir, args.threads, pop_map, LOG)
            run_gstacks(args.outdir, args.threads, pop_map, LOG)
            run_populations(args.outdir, args.threads, pop_map, args.maf, args.mac, args.mintype,
                                args.whitelist, args.blacklist, LOG)

    elif args.analysis == "post-ustacks":
        os.chdir(args.outdir)
        if [f for f in os.listdir('.') if f.endswith((".snps.tsv", ".tags.tsv", ".alleles.tsv"))]:
            if args.popmap:
                pop_map = args.popmap
            else:
                pop_map = generate_pop_map(args.outdir)
            run_cstacks(args.outdir, args.threads, pop_map, args.nmismatch, args.catnum, LOG)
            run_sstacks(args.outdir, args.threads, pop_map, LOG)
            run_tsv2bam(args.outdir, args.threads, pop_map, LOG)
            run_gstacks(args.outdir, args.threads, pop_map, LOG)
            run_populations(args.outdir, args.threads, pop_map, args.maf, args.mac, args.mintype,
                                args.whitelist, args.blacklist, LOG)
        else:
            raise ValueError("\n\n\nERROR: No ustacks output files present (#.snps.tsv, #.tags.tsv, #.alleles.tsv). Quitting.\n")
        
    elif args.analysis == "ustacks":
        run_ustacks(args.indir, args.outdir, args.threads, args.maxdist, args.mindepth, args.intid, LOG)
        
    elif args.analysis == "cstacks":
        os.chdir(args.outdir)
        if [f for f in os.listdir('.') if f.endswith((".snps.tsv", ".tags.tsv", ".alleles.tsv"))]:
            if args.popmap:
                pop_map = args.popmap
            else:
                pop_map = generate_pop_map(args.outdir)
            run_cstacks(args.outdir, args.threads, pop_map, args.nmismatch, args.catnum, LOG)
        else:
            raise ValueError("\n\n\nERROR:  No ustacks output files present. Quitting.\n")
       
    elif args.analysis == "sstacks":
        os.chdir(args.outdir)
        if [f for f in os.listdir('.') if f.endswith(("atalog.snps.tsv", "atalog.tags.tsv", "atalog.alleles.tsv"))]:
            if args.popmap:
                pop_map = args.popmap
            else:
                pop_map = generate_pop_map(args.outdir)
            run_sstacks(args.outdir, args.threads, pop_map, LOG)
        else:
            raise ValueError("\n\n\nERROR:  No catalog from cstacks present. Quitting.\n")
        
    elif args.analysis == "tsv2bam":
        os.chdir(args.outdir)
        if [f for f in os.listdir('.') if f.endswith((".matches.tsv"))]:
            if args.popmap:
                pop_map = args.popmap
            else:
                pop_map = generate_pop_map(args.outdir)
            run_tsv2bam(args.outdir, args.threads, pop_map, LOG)
        else:
            raise ValueError("\n\n\nERROR:  No sstacks output files present (#.matches.tsv). Quitting.\n")
        
    elif args.analysis == "gstacks":
        os.chdir(args.outdir)
        if [f for f in os.listdir('.') if f.endswith(("matches.bam"))]:
            if args.popmap:
                pop_map = args.popmap
            else:
                pop_map = generate_pop_map(args.outdir)
            run_gstacks(args.outdir, args.threads, pop_map, LOG)
        else:
            raise ValueError("\n\n\nERROR:  No tsv2bam output files present (#.matches.bam). Quitting.\n")
        
    elif args.analysis == "populations":
        os.chdir(args.outdir)
        if [f for f in os.listdir('.') if f.endswith((".fa.gz", ".calls", ".log.distribs"))]:
            if args.popmap:
                pop_map = args.popmap
            else:
                pop_map = generate_pop_map(args.outdir)
            run_populations(args.outdir, args.threads, pop_map, args.maf, args.mac, args.mintype,
                                args.whitelist, args.blacklist, LOG)
        else:
            raise ValueError("\n\n\nERROR:  No gstacks output files present (catalog.fa.gz). Quitting.\n")
        
    tf = datetime.now()
    show_elapsed(tb, tf, "analysis (-a {})".format(args.analysis), LOG)
    
if __name__ == '__main__':
    main()
