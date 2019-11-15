import sys
import os
import subprocess as sp
import shutil
from datetime import datetime
'''
Usage: python 1_stacks_pipeline.py [full path to directory containing all de-multiplexed trimmed fq files]
ex. python 1_stacks_pipeline.py /Volumes/Project1/3_Trimmed_Output/

Example file naming scheme (must follow the'SAMPLE.trim.fq' to work!):
4_MB312.trim.fq

Will result in many intermediate files:
4_MB312.trim.snps.tsv
4_MB312.trim.tags.tsv
4_MB312.trim.alleles.tsv

4_MB312.trim.matches.tsv


Full pipeline:
Script will run ustacks across all SAMPLE.trim.fq files in directory, then cstacks with all samples,
then sstacks on a per sample basis. Finally, the population script is run taking all samples as one
population, outputs several folders with results from population script with varying -r values 
(0, 50, 60, 70, 80, 90, 100).

User decides whether to run full pipeline (above), ustacks portion only, or cstacks sstacks and
populations portion only. As long as you stick to the same directory, you can run ustacks
first (decision b), then the remaining pipeline later (decision c). 

ustacks call:
	ustacks -t fastq -f SAMPLE.trim.fq -i {iterable} -o . -r -m 5 -M 2 -p {user decision}

cstacks call:
	cstacks -b 1 {-s all samples together} -p {user decision} -o .

sstacks call:
	sstacks -b 1 -c batch_1 -s SAMPLE.trim -p {user decision} -o .

populations call:
	populations -b 1 -P . -M {path to script-generated file} -r {0, 50, 60, 70, 80, 90, 100} -t {user decision} --min_maf 0.05

This script was written for and tested with Stacks version 1.35.

############################################
Written for Python 2.7.3

External Dependencies: Stacks v1.35
ustacks, cstacks, sstacks, populations (called from command line)
############################################

Dan Portik
daniel.portik@uta.edu
February 2016
'''

print '\n',"#########################################################################"
print "a. Full pipeline: ustacks, cstacks, sstacks, and populations."
print "b. Run ustacks only."
print "c. Run cstacks, sstacks, and populations (ustacks previously completed)."
print "d. Run populations only."
print "#########################################################################",'\n'

decision_analysis = (raw_input("Please select one of the above options (a, b, c), or 'q' to quit: "))
print '\n'

#get to fastq directory we want to be in
main_dir = sys.argv[1]
os.chdir(main_dir)

#****************************************************************
#Begin ustacks function
def ustacks():
    print '\n', '\n', "****************************************************************",'\n', "Beginning ustacks processing...", '\n', "****************************************************************", '\n', '\n'

    #count fastq files
    filecount = int(0)
    for filetype in os.listdir('.'):
        if filetype.endswith('.fq'):
            filecount += 1
    print '\n', "There are {} fastq files in this directory to analyze.".format(filecount), '\n'

    #create output file matching file and identifier code
    fh_out = open("Identifier.log", 'a')
    fh_out.write("Sample"+'\t'+"Identifier_code"+'\n')

    #call fast_trimmer for all fastq files
    identifier = int(1)

    for filetype in os.listdir('.'):
        if filetype.endswith('.fq'):

            #put identifier info in output log
            names = filetype.split('.')
            fh_out.write(names[0]+'\t'+'{}'.format(identifier)+'\n')

            #create string for ustacks system call and execute
            print "{0}: Beginning ustacks processing of {1}...".format(datetime.now().strftime("%b-%d %H:%M"), filetype), '\n'
            call_ustacks = "ustacks -t fastq -f {0} -i {1} -o {2} -r -m {3} -M {4} -p {5}".format(filetype, identifier, '.', '5', '2', processors)
            print call_ustacks, '\n'
            proc_ustacks = sp.call(call_ustacks, shell=True)
            print "Finished analysis of {}.".format(filetype), '\n', '\n'

            identifier += 1

    fh_out.close()

#****************************************************************
#Begin cstacks, sstacks, and populations function
def cstacks_sstacks_populations():

    #---------------------
    #Begin cstacks routine
    print '\n', '\n', "****************************************************************",'\n', "Beginning cstacks processing...", '\n', "****************************************************************", '\n', '\n'

    #create set of file names
    name_set = set()
    for filetype in os.listdir('.'):
        if '.trim' in filetype:
            names = filetype.split('.trim')
            name = names[0]+'.trim'
            name_set.add(name)
    name_list = list(name_set)
    name_list.sort()
    print '\n', "There are {} unique files to include for cstacks: ".format( len(name_list)), '\n'
    for name in name_list:
        print name
    print '\n', '\n'

    #concatenate file names into -s format for system call of cstacks
    name_cat = " "
    for name in name_list:
        name_cat = name_cat+'-s {} '.format(name)

    #call cstacks
    print "{0}: Beginning cstacks processing of files...".format(datetime.now().strftime("%b-%d %H:%M")), '\n'
    call_cstacks = "cstacks -b 1{0}-p {1} -o {2}".format(name_cat, processors, '.')
    print call_cstacks, '\n', '\n'
    proc_cstacks = sp.call(call_cstacks, shell=True)
    print "{0}: Finished cstacks analysis.".format(datetime.now().strftime("%b-%d %H:%M")), '\n', '\n'


    #---------------------
    #Begin sstacks routine
    print '\n', '\n', "****************************************************************",'\n', "Beginning sstacks processing...", '\n', "****************************************************************", '\n', '\n'

    #call sstacks
    for name in name_list:
        #call cstacks
        print "{0}: Beginning sstacks processing of file {1}...".format(datetime.now().strftime("%b-%d %H:%M"), name), '\n'
        call_sstacks = "sstacks -b 1 -c {0} -s {1} -p {2} -o {3}".format('batch_1',name,processors,'.')
        print call_sstacks, '\n'
        proc_sstacks = sp.call(call_sstacks, shell=True)
        print '\n', '\n'

    
    #---------------------
    #Begin population routine
    print '\n', '\n', "****************************************************************",'\n', "Beginning population processing...", '\n', "****************************************************************", '\n', '\n'

    #create single population map file
    pop_dir = 'population_directory'
    if not os.path.exists(pop_dir):
        os.mkdir(pop_dir)

    os.chdir(pop_dir)
    full_path = os.getcwd()
    full_name = full_path+'/population_map.txt'
    fh_pop = open('population_map.txt', 'a')
    for name in name_list:
        fh_pop.write(name+'\t'+'1'+'\n')
    fh_pop.close()
    os.chdir(main_dir)

    #---------------------
    #Begin first call of population
    basic_dir = 'population_no_r'
    if not os.path.exists(basic_dir):
        os.mkdir(basic_dir)
    pop_string_basic = "populations -b {0} -P {1} -M {2} -t {3} --min_maf {4}".format('1', '.', full_name, processors,'0.05')
    print '\n', '\n', "{0}: Beginning basic population assessment...".format(datetime.now().strftime("%b-%d %H:%M")), '\n'
    print pop_string_basic, '\n'
    proc_pop_basic = sp.call(pop_string_basic, shell=True)

    #move files to prevent overwriting
    for filetype in os.listdir('.'):
        if filetype.startswith('batch_1'):
            names = filetype.split('.')
            if names[1] != 'catalog':
                shutil.move(filetype, basic_dir)
                
    #---------------------
    #begin calls of population involving different r values
    r_list = ["50", "60", "70", "80", "90", "100"]
    for r in r_list:
        iter_dir = 'population_r{}'.format(r)
        if not os.path.exists(iter_dir):
            os.mkdir(iter_dir)
        pop_string_iter = "populations -b {0} -P {1} -M {2} -r {3} -t {4} --min_maf {5}".format('1', '.', full_name, r, processors,'0.05')
        print '\n', '\n', "{0}: Beginning population assessment of -r {1}...".format(datetime.now().strftime("%b-%d %H:%M"), r), '\n'
        print pop_string_iter, '\n'
        proc_pop_iter = sp.call(pop_string_iter, shell=True)

        #move files to prevent overwriting
        for filetype in os.listdir('.'):
            if filetype.startswith('batch_1'):
                names = filetype.split('.')
                if names[1] != 'catalog':
                    shutil.move(filetype, iter_dir)

#****************************************************************
def populations():
    #---------------------
    #Begin population routine
    print '\n', '\n', "****************************************************************",'\n', "Beginning population processing...", '\n', "****************************************************************", '\n', '\n'

    #create set of file names
    name_set = set()
    for filetype in os.listdir('.'):
        if '.trim' in filetype:
            names = filetype.split('.trim')
            name = names[0]+'.trim'
            name_set.add(name)
    name_list = list(name_set)
    name_list.sort()
    print '\n', "There are {} unique files to include for populations module: ".format( len(name_list)), '\n'
    for name in name_list:
        print name
    print '\n', '\n'

    #create single population map file, unless already present
    pop_dir = 'population_directory'
    if not os.path.exists(pop_dir):
        os.mkdir(pop_dir)
    os.chdir(pop_dir)
    full_path = os.getcwd()
    full_name = full_path+'/population_map.txt'
    if not os.path.exists(full_name):
    	fh_pop = open('population_map.txt', 'a')
    	for name in name_list:
        	fh_pop.write(name+'\t'+'1'+'\n')
    	fh_pop.close()
    os.chdir(main_dir)
                
    #---------------------
    #begin calls of population involving different r values
    r_list = ["50", "60", "70", "80", "90", "100"]
    for r in r_list:
        iter_dir = 'population_r{}'.format(r)
        if not os.path.exists(iter_dir):
            os.mkdir(iter_dir)
        pop_string_iter = "populations -b {0} -P {1} -M {2} -r {3} -t {4} --min_maf {5}".format('1', '.', full_name, r, processors,'0.05')
        print '\n', '\n', "{0}: Beginning population assessment of -r {1}...".format(datetime.now().strftime("%b-%d %H:%M"), r), '\n'
        print pop_string_iter, '\n'
        proc_pop_iter = sp.call(pop_string_iter, shell=True)

        #move files to prevent overwriting
        for filetype in os.listdir('.'):
            if filetype.startswith('batch_1'):
                names = filetype.split('.')
                if names[1] != 'catalog':
                    shutil.move(filetype, iter_dir)


#****************************************************************
#Run routines based on user decisions

#Full pipeline:
if decision_analysis == "a":
    #Get processor info
    processors = None
    while processors is None:
        try:
            processors = int(raw_input("Number of processors to use (ex. 4): "))
        except ValueError:
            print "That wasn't a number."

    #Record start time
    t0 = datetime.now()

    #run routines
    ustacks()
    cstacks_sstacks_populations()

    #Calculate total time for processing
    print '\n', '\n'
    print "Full stacks pipeline finished!", '\n'
    t1 = datetime.now()
    total_t = t1 - t0
    print "Total time: {}".format(total_t)
    print '\n', '\n'

#ustacks only    
elif decision_analysis == "b":
    #Get processor info
    processors = None
    while processors is None:
        try:
            processors = int(raw_input("Number of processors to use (ex. 4): "))
        except ValueError:
            print "That wasn't a number."

    #Record start time
    t0 = datetime.now()

    #Run routine
    ustacks()
    
    #Calculate total time for processing
    print '\n', '\n'
    print "ustacks processing finished!", '\n'
    t1 = datetime.now()
    total_t = t1 - t0
    print "Total time: {}".format(total_t)
    print '\n', '\n'

#cstacks, sstacks, populations only
elif decision_analysis == "c":
    #Get processor info
    processors = None
    while processors is None:
        try:
            processors = int(raw_input("Number of processors to use (ex. 4): "))
        except ValueError:
            print "That wasn't a number."

    #Record start time
    t0 = datetime.now()

    #Run routines
    cstacks_sstacks_populations()
    
    #Calculate total time for processing
    print '\n', '\n'
    print "All cstacks, sstacks, and populations processing finished!", '\n'
    t1 = datetime.now()
    total_t = t1 - t0
    print "Total time: {}".format(total_t)
    print '\n', '\n'

#populations only
elif decision_analysis == "d":
    #Get processor info
    processors = None
    while processors is None:
        try:
            processors = int(raw_input("Number of processors to use (ex. 4): "))
        except ValueError:
            print "That wasn't a number."

    #Record start time
    t0 = datetime.now()

    #Run routines
    populations()
    
    #Calculate total time for processing
    print '\n', '\n'
    print "All populations processing finished!", '\n'
    t1 = datetime.now()
    total_t = t1 - t0
    print "Total time: {}".format(total_t)
    print '\n', '\n'


elif decision_analysis == "q":
    print "Quitting program now", '\n'

    
else:
    print "Not a valid decision, now quitting", '\n'
