import sys
import os
import subprocess as sp
import shutil
import datetime
'''
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

The 'barcode.txt' file will be used for the particular INDEX.gz (zipped!) file located in the same subdirectory as the file,
which currently has several barcoded individuals within it. In my implementation of the 'process_radtags' script of stacks v1.35
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


############################################
Written for Python 2.7.3

External Dependencies:
fastqc (call from command line)
stacks (tested for v 1.35):
    process_radtags (called from command line, put stacks in path)
fastx (call from command line):
    fastx_trimmer
############################################

Dan Portik
daniel.portik@uta.edu
February 2016
'''

#Decide on fastqc analyses
print '\n','\n'
decision_fastqc = (raw_input("Perform fastqc quality analysis on all demultiplexed samples (y,n)? "))
print '\n','\n'
#Decide how many bases to cut
decision_bp = None
while decision_bp is None:
    try:
        decision_bp = 1 + int(raw_input("Enter number of base pairs to trim from beginning of reads after removing barcode (for restriction site overhang, ex. 6): "))
    except ValueError:
        print "That wasn't a number."
print '\n','\n'



#record start time
t0 = datetime.datetime.now()

#get path to directory we want to be in
bin_directory = sys.argv[1]
os.chdir(bin_directory)

#Create list of directories here and begin looping
file_list = os.listdir(bin_directory)
for folder in file_list:
    if os.path.isdir(folder) is True:
        os.chdir(folder)
        print "-----------------------------------------------"
        print '\n', "{0}: Moving to {1} to begin demultiplexing.".format(datetime.datetime.now().strftime("%b-%d %H:%M"),folder), '\n'
        print "-----------------------------------------------"

#*******************************************
#Begin demultiplexing routine
        bcode = "barcode.txt"
        for filetype in os.listdir('.'):
            if filetype.endswith('.gz'):
                print "Beginning demultiplexing of {}...".format(filetype)
                
                call_string_recheck = "process_radtags -f {0} -i {1} -y {2} -o {3} -b {4} -c -q -r -e {5}".format(filetype,"gzfastq","fastq",".",bcode,'sbfI')
                #call_string_reoff = "process_radtags -f {0} -i {1} -y {2} -o {3} -b {4} -c -q -r -e {5} --disable_rad_check".format(filetype,"gzfastq","fastq",".",bcode,'sbfI')
                #call_string_window = "process_radtags -f {0} -i {1} -y {2} -o {3} -b {4} -c -q -r -e {5} -w {6} --disable_rad_check".format(filetype,"gzfastq","fastq",".",bcode,'sbfI','0.3')
                
                proc_stacks = sp.call(call_string_recheck, shell=True)
                print '\n', '\n', "-----------------------------------------------",'\n', "Finished demultiplexing {}.".format(filetype)
                print "-----------------------------------------------", '\n', '\n'

#*******************************************
#Begin fastqc routine
        if decision_fastqc == "y":
            for filetype in os.listdir('.'):
                if filetype.endswith('.fq'):
                    print "{0}: Beginning fastqc analysis of {1}...".format(datetime.datetime.now().strftime("%b-%d %H:%M"),filetype), '\n'
                    names = filetype.split('.')
                    dir_name = names[0]
                    if not os.path.exists(dir_name):
                        os.mkdir(dir_name)
                    call_fastqc = "fastqc -t 2 {0} -o {1}".format(filetype,dir_name)
                    proc_fastqc = sp.call(call_fastqc, shell=True)
                    print '\n', '\n'
            print '\n', '\n'
            os.chdir(bin_directory)
        elif decision_fastqc == "n":
            os.chdir(bin_directory)

#*******************************************
#Begin log file routine
print "Combining stacks log files..."

fh_out = open("Stacks_ProcessRadtags.log", 'a')
header_list = []
contents_list = []

for folder in file_list:
    if os.path.isdir(folder) is True:
        os.chdir(folder)
        for filetype in os.listdir('.'):
            if filetype.endswith('.log'):
                log_name = folder+"/"+filetype
                print "Getting contents from {}...".format(log_name)
                temp_fh = open(filetype, 'r')
                lines = temp_fh.readlines()
                for line in lines:
                    line = line.strip()
                    if line.startswith('Barcode'):
                        temp_header_list = line.split('\t')
                        if len(temp_header_list) > int(4):
                            header_list.append(temp_header_list)
                            indices = len(temp_header_list)
                for line in lines:
                    line = line.strip()
                    temp_cont_list = line.split('\t')
                    if len(temp_cont_list) == indices:
                        if temp_cont_list[0] != 'Barcode' and len(temp_cont_list[0]) == int(5):
                            contents_list.append(temp_cont_list) 
        os.chdir(bin_directory)

for item in header_list[0]:
    fh_out.write(item+'\t')
fh_out.write('\n')
for sublist in contents_list:
    for item in sublist:
        fh_out.write(item+'\t')
    fh_out.write('\n')

fh_out.close()

#*******************************************
#Move newly created fastq files to new directory at level of main directory

os.chdir('..')
up_one = os.getcwd()
out_dir = up_one+'/2_Demultiplex_Output'
if not os.path.exists(out_dir):
    os.mkdir(out_dir)
    
os.chdir(bin_directory)
for folder in file_list:
    if os.path.isdir(folder) is True:
        os.chdir(folder)
        for filetype in os.listdir('.'):
            if filetype.endswith('.fq'):
                shutil.move(filetype, out_dir)
        os.chdir(bin_directory)

print '\n', '\n'

#*******************************************
#Begin trimming subroutine

#move to output directory
os.chdir(out_dir)

file_count = int(0)

#call fast_trimmer for all fastq files
for filetype in os.listdir('.'):
    if filetype.endswith('.fq'):
        file_count+=1
        print "Beginning trimming of {}...".format(filetype)
        names = filetype.split('.')
        newfilename = names[0]+'.trim.fq'
        call_fastx = "fastx_trimmer -Q33 -f {0} -i {1} -o {2}".format(decision_bp, filetype, newfilename)
        proc_fastx = sp.call(call_fastx, shell=True)
        print "Output written as {}.".format(newfilename), '\n', '\n'

#Build output directory
out_dir2 = "3_Trimmed_Output"
if not os.path.exists(out_dir2):
    os.mkdir(out_dir2)

#Move files to output directory
filecount = int(0)    
for filetype in os.listdir('.'):
    if filetype.endswith('.trim.fq'):
        shutil.move(filetype, out_dir2)
        filecount += 1
        
#Move output directory up one level
shutil.move(out_dir2, up_one)

print "Finished trimming {} fastq files.".format(filecount)
print "Trimmed files are located in {}.".format(out_dir), '\n', '\n'


#*******************************************
#Time it!
print '\n', '\n'
print "All filtering and trimming is finished!", '\n'
t1 = datetime.datetime.now()
total_t = t1 - t0
hours =  ( float(total_t.seconds) / float(60)) / float(60)
print "Total time: {} hours".format(hours)
print '\n', '\n'
