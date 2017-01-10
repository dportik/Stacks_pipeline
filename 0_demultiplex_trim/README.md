Usage: python 0_demultiplex_trim.py [full path to main directory containing several subdirectories, each subdir with an INDEX.gz file from the sequencer]

ex. python 0_demultiplex_trim.py /Volumes/Project1/Main_Dir/

Requires the following folder and file structure to run properly:

Main Directory 
	Subdirectory
		FILE.gz
		barcode.txt
	Subdirectory
		FILE.gz
		barcode.txt
	etc.

I name my main directory folder "1_Raw_data", as additional directories are created on the
same level that are called '2_Post_Demultiplex' and '3_Trimmed_Output' during this script.
In each subdirectory you will need to create a barcode file called "barcode.txt" with the relevant 
information for your samples.

Mine are structured like this so the samples get renamed based on their barcode:
CGATC	1_CUMV_15058
GGTTG	1_CUMV_15060
TGCAT	1_CUMV_15175
AAGGA	1_CUMV_15209
TCGAT	1_CAS_258004
AACCA	1_CAS_258161
GCATG	1_DMP_647
CAACC	1_CUMV_14900
T

he 'barcode.txt' file will be used for the particular INDEX.gz (zipped!) file located in the same subdirectory as the file,
which currently has several barcoded individuals within it. In my implementation of the 'process_barcodes' script of stacks v1.35

I assume these are single (not paired-end) reads and the enzyme used is 'sbfI'. The overall structure of the 'process_barcodes'
command is below:
         process_radtags -f INDEX.gz -i gzfastq -y fastq -o . -b barcode.txt -c -q -r -e sbfI
After demultiplexing, the script will run 'fastqc' on each 'SAMPLE.fq' file before moving to next subdirectory. The output of 
each sample is contained in a separate folder in the subdirectory, which can be visualized as an html page.
The script will then collect relevant information from the log files from each subdirectory and concatenate them 
in a tab-delimited text file in the main directory, called 'Stacks_ProcessRadtags.log'. This places all of
the stacks filtering information in one place, which can be opened in excel or a text editor. This is 
the best place to look for the number of reads per sample, and those that pass initial filtration. 

All the de-multiplexed sample files are then moved to a new directory at the same level as the main directory,
called '2_Post_Demultiplex'.

At the beginning of the script the user is prompted to select how many bases to trim off the beginning of the reads. This
corresponds to the restriction site overhang, which is not trimmed by process_radtags (though barcodes are removed). For 'sbfI'
there are 6 bases that need to be trimmed. The selected number of sites will be cut from the beginning of all reads present
using fastx_trimmer.

The trimmed reads are renamed and then moved to a final folder, '3_Trimmed_Output'.  If the demultiplexed sample is called
"SAMPLE_1.fq", it will be renamed to SAMPLE_1.trim.fq during this process.  This directory and contents can be used as the basis for
the next script, '0_stacks_pipeline.py'.

Written for Python 2.7.3
External Dependencies:
fastqc (call from command line)
stacks (tested for v 1.35):
    process_radtags (called from command line, put stacks in path)
fastx (call from command line):
    fastx_trimmer


# Dan Portik
daniel.portik@uta.edu
February 2016
