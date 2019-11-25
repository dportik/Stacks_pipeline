# Stacks Pipeline <a name=TOP></a>

---------------

Welcome to the Stacks Pipeline wiki! Instructions for using the pipeline can be found here. Major steps of the pipeline are outlined below, and the links jump the relevant section for easy navigation.

+ ### [**Demultiplexing and trimming**](#DT)
+ ### [**Running Stacks**](#RS)
+ ### [**Filtering tsv files**](#FTF)
+ ### [**Converting filtered tsv files**](#CTF)

These sections contain a summary of the tasks performed, usage instructions for the relevant script, and examples of outputs.

**Please note that helpful information can always be printed to the screen for any script by running it using the `-h` option.**

---------------

# **Demultiplexing and trimming** <a name=DT></a>

## Contents

[**Overview**](#DTOV)

[**Script arguments**](#DTSA)

[**Organizing the sequencing data**](#DTOSD)

[**Running the script: Demultiplex and trim cutsites**](#DTDT)
	
[**Running the script: Remove UMIs, demultiplex, and trim cutsites**](#DTUDT)
	
[**Outputs**](#DTO)

## Overview <a name=DTOV></a>

For this step the `Demultiplex_Trim.py` script will be used. The script will perform the following tasks in order:

+ Trim bases from the front of reads in fastq.gz files to remove unique molecular identifiers (UMIs). The length of the UMI is specified by the user. The `fastx_trimmer` program is used to trim bases. **This step is optional.**
+ Demultiplex all fastq.gz files using user-supplied barcode files and the `process_radtags` module of Stacks v2.4.
+ Trim bases from the front of reads in the demultiplexed fastq files to remove RAD cutsites. The length of the cutsites is specified by the user. The `fastx_trimmer` program is used to trim bases. 
+ All demultiplexed and trimmed fastq files are moved to the specified output directory.

## Script arguments <a name=DTSA></a>

The `Demultiplex_Trim.py` script has three mandatory flags and one optional flag:

**Mandatory Arguments:**

+ `-i <path-to-directory>`: The full path to a directory which contains the subdirectories containing the fastq.gz sequencing files (see section below).

+ `-o <path-to-directory>`: The full path to an *existing* directory to write the demultiplexed and trimmed fastq files.

+ `-b <integer>`: Supply an integer representing the first base position to keep in the demultiplexed fq files. All sites before this position will be removed. This will remove all RAD cutsites from the demultiplexed reads. If you are using *sbfI*, it is 6 bp long, so you would enter 7.

**Optional Arguments**:

+ `-u <integer>`: Supply an integer representing the first base position to keep in the fastq.gz files prior to demultiplexing. This is intended to remove all UMI (unique molecular identifier) sites from the raw reads. If the length of the UMI is 8 bp, you would enter 9.

Before running the script, you will need to organize your sequence data. 

## Organizing the sequencing data <a name=DTOSD></a>

This `Demultiplex_Trim.py` script is designed to run on an input directory that contains subdirectories. Each subdirectory should contain a `fastq.gz` sequencer file, and a corresponding barcode file called `barcode.txt`. The organizational structure should look like this:

```
Raw-Data
│
├── S127
│	├── FP5-ATCACG_S127.fastq.gz
│	└── barcode.txt
│
├── S128
│	├── FP5-CGATGT_S128.fastq.gz
│	└── barcode.txt
│
├── S129
│	├── FP5-TGACCA_S129.fastq.gz
│	└── barcode.txt
│
├── S130
│	├── FP5-TTAGGC_S130.fastq.gz
│	└── barcode.txt
```

Here, the input directory `Raw-Data` contains four subdirectories, each of which contains a gzipped fastq file and a barcode text file. In this example, the `-i ` flag should contain the full path to `Raw-Data`. 

The barcode file is a tab-delimited text file that contains the barcode sequence in column one, and the sample name in column two:

```
AGTAAG	din_9250
ACTAGG	din_9251
AGCATT	din_9252
GTCTAT	hut_CD15010
GACCAA	hut_CD15012
CATCTC	hut_CFS1504g
TCTTAG	hut_CG15146
CAGATG	hut_CG15148
GGGATA	hut_CG15149
```

The barcode file must have the same name in each subdirectory (`barcode.txt`). The barcode file in a given subdirectory will only be used for the fastq.gz file also located in that directory.

## Running the script

### 1. Demultiplex and trim cutsites <a name=DTDT></a>

The basic routine involves demultiplexing all fastq.gz files and removing the RAD cutsites from the resulting demultiplexed fastq files. To perform this requires specifying the input directory (`-i`), output directory (`-o`), and number of bases to trim in the demultiplexed reads (`-b`). For example:

```
python Demultiplex_Trim.py -i /user/Proj/1-Raw-Data -o /user/Proj/2-Demultiplexed-Data -b 7
```
The above command would demultiplex all fastq.gz files found in the subdirectories in the `/1-Raw-Data` folder, trim 6 bases from the front of all reads in the demultiplexed fastq files, then move those trimmed files to the `/2-Demultiplexed-Data` folder. 

When the script completes, you should see the demultiplexed fastq files in the output directory specified. They will be labeled as `sample.trim.fastq`, where the `sample` component is the name provided in the `barcode.txt` file. 

Under the hood, to demultiplex the gzipped fastq files the equivalent command is used internally to call the `process_radtags` module:

```
process_radtags -f <input fastq.gz> -i gzfastq -y fastq -o <output directory> -b barcode.txt -c -q -r -e sbfI --inline_null
```

After demultiplexing, the equivalent command is used internally to call `fastx_trimmer` to trim the RAD cutsites in all resulting fastq files:

```
fastx_trimmer -Q33 -f <base position> -i <input fastq> -o <output fastq>
```

### 2. Remove UMIs, demultiplex, and trim cutsites <a name=DTUDT></a>

In addition to the normal routine, there is an option to remove UMI sites from the raw reads before beginning the demultiplexing. In order to perform this, you will need to supply the `-u` flag with an appropriate integer. For example:

```
python Demultiplex_Trim.py -i /user/Proj/1-Raw-Data -o /user/Proj/2-Demultiplexed-Data -b 7 -u 8
```

The above command would remove 7 bases from all the raw reads in the fastq.gz files found in the subdirectories in the `/1-Raw-Data` folder, demultiplex the trimmed fastq.gz files, trim 6 bases from the front of all reads in the demultiplexed fastq files, then move those trimmed files to the `/2-Demultiplexed-Data` folder. 

When the script completes, you should see the demultiplexed fastq files in the output directory specified. They will be labeled as `sample.trim.fastq`, where the `sample` component is the name provided in the `barcode.txt` file. 

To trim UMI's from the fastq.gz files, the following equivalent command is used internally:

```
gunzip -c <input fastq.gz> | fastx_trimmer -Q33 -f <base position> -z -o <output fastq.gz>
```

This essentially streams the unzipped fastq to `fastx_trimmer`, which trims the appropriate number of bases, and then writes the output to a second gzipped fastq file.
The output file will be labeled with a `_UMI_trimmed` component. For example, trimming the UMIs in:

`FP5-ATCACG_S127.fastq.gz` 

would result in a second file called:

`FP5-ATCACG_S127_UMI_trimmed.fastq.gz`. 

## Outputs <a name=DTO></a>

The main outputs from the script will be the demultiplexed, trimmed fastq files. They will be written to output directory specified by the `-o` flag. These files should be used to for the next major step ([Running Stacks](#RS)). The output fastq files will be labeled as `sample.trim.fastq`, where the `sample` component is the name provided in the `barcode.txt` file. For example:


```
Output-Fastq
│
├── din_9250.trim.fq
├── din_9251.trim.fq
├── din_9252.trim.fq
├── hut_CD15010.trim.fq
├── hut_CD15012.trim.fq
├── hut_CFS1504g.trim.fq
├── hut_CG15146.trim.fq
```


In addition to the fastq files, in each subdirectory will be a `process_radtags.log` file that summarizes information about each sample. In the main input directory (`-i`), a file called `Stacks_ProcessRadtags.log` will be present. It is a combination of all the log files, and will have information about every sample that was demultiplexed. Here is an example of the contents:

```
Barcode	Filename	Total	No_RadTag	Low_Quality	Retained
AGTAAG	din_9250	228057	15443	34	212580
ACTAGG	din_9251	2780362	83405	553	2696404
AGCATT	din_9252	1316292	29878	241	1286173
GTCTAT	hut_CD15010	712783	23747	154	688882
GACCAA	hut_CD15012	1509790	60509	278	1449003
CATCTC	hut_CFS1504g	40737	17595	8	23134
TCTTAG	hut_CG15146	465987	16957	75	448955
CAGATG	hut_CG15148	47307	46627	0	680
GGGATA	hut_CG15149	1280957	54748	231	1225978
CGAATG	hut_CRSN71	1209829	26465	232	1183132
TAAGAC	hut_EBG1990	2326188	15370	450	2310368
ATAACC	hut_MSK013	784996	42653	154	742189
AACGGT	hut_MSK033	4119648	26256	829	4092563
GACGTT	hut_MSK154	1384951	12539	286	1372126
ATGTCC	hut_UMA1267	197453	13355	32	184066
GCAGAA	hut_YOKO0382	1698448	35329	321	1662798
TGGGAT	oc_GFMJ2364	2690018	44224	509	2645285
TGTTGG	oc_GFMJ2427	114137	23849	22	90266
TCAATC	oc_GFMJ2428	1369575	15558	276	1353741
CGAAAC	oc_GFMJ2430	2395667	24170	451	2371046
TCTGCT	oc_GFMJ2076	2693117	38798	513	2653806
ACCAAA	oc_GFMJ2433	1971237	60042	408	1910787
CCCATA	oc_GFMJ2459	2770265	18368	543	2751354
```

This file should be useful for initial troubleshooting. For example, making sure you trimmed the correct number of base pairs for the UMI and/or restriction cutsite!

[Back to top](#TOP)

---------------

# **Running Stacks** <a name=RS></a>

## Contents

[**Overview**](#RSOV)

[**Script arguments**](#RSSA)

[**Considerations for running Stacks**](#RSCRS)

[**Running the full Stacks workflow, the simple way**](#RSRSWS)
	
[**Running the full Stacks workflow for big datasets**](#RSRSWB)
	 
[**Running the partial Stacks workflow**](#RSRSWP)

[**Stacks workflow details...**](#RSSWD)

[**Outputs**](#RSO) <a name=RSO></a>

## Overview <a name=RSOV></a>

For this step the `Run_Stacks.py` script will be used. The intended purpose of the script is to automate the complete Stacks v2.4 workflow for a set of input fastq files. The Stacks workflow involves performing the following steps:

1. Run **ustacks**: Takes a fastq file and aligns the reads into stacks, forms putative loci, and identifies SNPs. This is performed independently for each fastq file. 
2. Run **cstacks**: Builds a catalog of loci from all the ustacks output files. The loci found across all the samples are used to create consensus loci, merging alleles together. 
3. Run **sstacks**: The stacks for a given sample (found using ustacks) are compared to the catalog to identify loci. 
4. Run **tsv2bam**: Converts outputs to a format that is used in the next step.
5. Run **gstacks**: On a per locus basis, aligns reads from all samples to identify SNPs genotype individuals, and phase SNPs into set of haplotypes.
6. Run **populations**: Performs filtering tasks and ouputs phased haplotypes for all samples across all loci. The outputs of the populations module (specifically the `populations.haplotypes.tsv` files) are the targets for the subsequent steps.

The `Run_Stacks.py` script was designed to take as input the directory containing all the demultiplexed and trimmed fastq files (the output from `Demultiplex_Trim.py`), copy those files to an analysis directory, and run the full Stacks workflow in the analysis directory. A population map file is automatically generated, and it assigns all the fastq files to a single population. This map file is used to call all the samples for the post-ustacks steps. However, the script is extremely flexible and allows any Stacks module to be accessed directly, as well as certain chunks of the workflow (for example, run all steps beyond ustacks). A custom population map file can also be used instead. 

Running the script requires only two flags, the path to the directory to perform the Stacks analyses (`-o`), and the type of Stacks analysis to run (`-a`). 

One major requirement is that there are fastq files are in the `-o` directory. There are two ways to accomplish this:

+ The directory containing the trimmed fastq files can be specified  with `-i` flag. If this option is used, the fastq files are copied to the `-o` directory, and then the Stacks workflow is run in the `-o` directory. This is the intended workflow, which involves picking up where `Demultiplex_Trim.py` left off.  During the copy, the fastq files are renamed (removing the `.trim.` label if present). This is also the required usage for parallelizing the ustacks step (see sections below).

+ If copying the fastq files to a new directory is not desired, their location can be specified directly using the `-o` flag. If this option is used, the Stacks workflow will be run in the `-o` directory using the fastq files present. If there are no fastq files in the `-o` directory, it will crash. Also, if you are picking up where `Demultiplex_Trim.py` left off, the `.trim.` label will not be removed from the fastq file labels. This means it will become part of the sample names, including in the output files (e.g. `din_MVZ23535.trim` instead of `din_MVZ23535`). 

There are many optional flags that control the settings in one or more Stacks modules. Before providing examples of how to use the script, the complete set of arguments is explained below.

## Script Arguments <a name=RSSA></a>

**Mandatory Arguments:**

+ `-o <path-to-directory>`: The full path to an existing directory to run the Stacks pipeline steps. If fastq files are not already present here, their location must be specified with the optional `-i` flag.

+ `-a <choice>`: Choices = *full*, *post-ustacks*, *ustacks*, *cstacks*, *sstacks*, *tsv2bam*, *gstacks*, *populations*. Specify running the complete Stacks workflow (full), all steps beyond ustacks (post-ustacks), or one of the specific modules.

**Optional Arguments**:

+ `-i <path-to-directory>`: The full path to a directory which contains the starting fastq files, which are copied to the location specified by `-o` to run Stacks.
    
+ `-t <integer>`: Specifies the number of threads to use.
    
+ `-p <path-to-file>`: **Affects all modules except ustacks.** Provide the full path to a custom population map file to use instead of the default. In the default population map, all samples are auto-assigned to a single population that is labeled 1.
    
+ `-M <integer>`: **For ustacks**. Maximum distance (in nucleotides) allowed between stacks. Default value = 2.

+ `-m <integer>`: **For ustacks**. Minimum depth of coverage required to create a stack. Default value = 3.
    
+ `--intid <integer>`: **For ustacks**. Specify an integer ID to start at when labeling samples (overrides starting at 1). Useful for running ustacks on multiple fastq file batches for a single project (different `-i`, same `-o`).
    
+ `-n <integer>`: **For cstacks**. Number of mismatches allowed between sample loci when building the catalog. Default value = 2.
    
+ `--catnum <integer>`: **For cstacks**. Specify the number of samples to build the catalog from, rather than use all samples. This many samples will be drawn randomly from the population map file.
    
+ `--maf <float>`: **For populations**. Specify a minimum minor allele FREQUENCY (0 - 0.5) required to process a nucleotide site at a locus. 
    
+ `--mac <integer>`: **For populations**. Specify a minimum minor allele COUNT (an integer) required to process a nucleotide site at a locus. Setting to 2 should eliminate singletons. 
    
+ `--mintype <choice>`: **For populations**. Choices = *r*, *R*. Specify whether the minimum percentage of samples **within** a population (r) or **across** populations (R) should be used to process loci. The R option should really only be used when multiple populations are defined in a custom population map using `-p`. Default = r.
    
+ `--whitelist <path-to-file>`: **For populations**. Full path to a file containing whitelisted markers to be included in the export.
    
+ `--blacklist <path-to-file>`: **For populations**. Full path to a file containing blacklisted markers to be excluded from the export.
        
+ `-l <string>`: Add a unique label to end of Stacks pipeline log file name. Must use if running ustacks on multiple fastq file batches for a single project (different `-i`, same `-o`). to prevent log file overwriting.

## Considerations for running Stacks <a name=RSCRS></a>

Building a high quality SNP dataset with Stacks requires some careful thought. As highlighted by [Paris et al. (2017)](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12775), the main parameters affecting the resulting ddRADseq datasets are the `-m` and `-M` parameters used in **ustacks** and the `-n` parameter used in **cstacks**. The `-m` parameter sets the minimum number of raw reads required to form a stack (a putative allele). The `-M` parameter
sets the number of mismatches allowed between stacks (putative alleles) to merge them into a putative locus. The `-n` parameter sets the number of mismatches allowed between stacks (putative loci) during construction of the catalog. In this pipeline, all default values follow the recommendations of [Paris et al. (2017)](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12775) (`-m 3`, `-M 2`, `-n 2`), but you may wish to change these values based on other considerations, including the read length, phylogenetic scope, and/or stringency required. Each of these parameters can be easily changed using the optional flags in this script with a new value (also labeled `-m`, `-M`, and `-n` for consistency). 

In addition to these main parameters, a variety of filtering options are available in the **populations** module. The `-r` parameter controls the percentage of individuals **within** a population that must possess a particular locus for it to be included. In this pipeline, several `-r` values are automatically set, including 10, 20, 30, 40, 50, 60, 70, 80, 90, and 100. Each `-r` value produces a distinct output directory where all the outputs from that particular run can be found, and all the output files are labeled using the specific r value. However, the `-R` option can be used instead, in which `-R` controls the minimum percentage of individuals **across** populations required to process a locus. This option really only makes sense if you have supplied a custom population map (using the `-p` flag) and the custom map contains two or more populations. To switch between `-r` and `-R` settings, the `--mintype` flag should be used. The default is to use `--mintype r`, but supplying `--mintype R` would invoke the `-R` setting in **populations** instead. 

Other filtering options are available in **populations**, including specifying a minor allele frequency (`--maf`) or a minor allele count (`--mac`), and supplying a whitelist of loci (`--whitelist`) or a blacklist of loci (`--blacklist`). There are no defaults set for these options, meaning they will not be invoked unless you supply the flag and the associated value or file path. 


## Running the script

### Running the full Stacks workflow, the simple way <a name=RSRSWS></a>

The original intention of this script was to process all the output fastq files produced by the `Demultiplex_Trim.py` script. Let's say the some output was produced by `Demultiplex_Trim.py` in the output directory called `Output-Fastq`, and we want to run the full Stacks workflow in the `Stacks-Run` directory:

```
Analysis
│
├── Output-Fastq
│	├── din_9250.trim.fq
│	├── din_9251.trim.fq
│	├── din_9252.trim.fq
│	├── hut_CD15010.trim.fq
│	├── hut_CD15012.trim.fq
│	├── hut_CFS1504g.trim.fq
│	└── hut_CG15146.trim.fq
│
├── Stacks-Run
│
```

We would use the following command to accomplish this task:

```
python Run_Stacks.py -i /Analysis/Output-Fastq -o /Analysis/Stacks-Run -a full -t 4
```

In the above command, the `-t` flag simply specifies the number of threads to use. This simple command may also include optional flags that affect **ustacks** (`-M`, `-m`), **cstacks** (`-n`, `--catnum`), and **populations** (`--maf`, `--mac`). Using the `-i` and `-o` flags means that the fastq files in `-i` are first copied to `-o` before Stacks is run in the `-o` directory. 

And that's it! You should end up with a series of populations output directories that can then be used for the next step:

```
Analysis
│
├── Output-Fastq
│	├── din_9250.trim.fq
│	├── din_9251.trim.fq
│	├── din_9252.trim.fq
│	├── hut_CD15010.trim.fq
│	├── hut_CD15012.trim.fq
│	├── hut_CFS1504g.trim.fq
│	└── hut_CG15146.trim.fq
│
├── Stacks-Run
│	├── din_9250.fq
│	├── din_9251.fq
│	├── din_9252.fq
│	├── hut_CD15010.fq
│	├── hut_CD15012.fq
│	├── hut_CFS1504g.fq
│	├── hut_CG15146.fq
│	├──    and..
│	├──    lots...
│	├──    of..
│	├──    other...
│	├──    files...
│	│	
│	├── Populations_r10
│	│	└── outputs...
│	│	
│	├── Populations_r20
│	│	└── outputs...
│	│	
│	├── Populations_r30
│	│	└── outputs...
│	│	
│	├── Populations_r40
│	│	└── outputs...
│	│	
│	├── Populations_r50
│	│	└── outputs...
```

In each of the populations output directories, there will be a `populations.haplotypes.tsv` file that is the target of the next step in the pipeline.



With the usage above, the fastq files in the `Output-Fastq` folder (`-i`) were copied to the `Stacks-Run` folder (`-o`), and the `.trim` component was removed from each fastq label. 

```
Analysis
│
├── Output-Fastq
│	├── din_9250.trim.fq
│	├── din_9251.trim.fq
│	├── din_9252.trim.fq
│	├── hut_CD15010.trim.fq
│	├── hut_CD15012.trim.fq
│	├── hut_CFS1504g.trim.fq
│	└── hut_CG15146.trim.fq
│
├── Stacks-Run
│	├── din_9250.fq 	<- copied and renamed!
│	├── din_9251.fq 	<- copied and renamed!
│	├── din_9252.fq 	<- copied and renamed!
│	├── hut_CD15010.fq 	<- copied and renamed!
│	├── hut_CD15012.fq 	<- copied and renamed!
│	├── hut_CFS1504g.fq <- copied and renamed!
│	├── hut_CG15146.fq 	<- copied and renamed!
```

After the copying was completed, the Stacks workflow was then run on the `Stacks-Run` directory (`-o`) using the newly copied fastq files. This is the intended way to run the script, but it does copy the files (taking up more space!).

If copying the fastq files is not desired, the `-i` flag should be excluded and their location can be specified using the `-o` flag instead. For example:

```
python Run_Stacks.py -o /Analysis/Output-Fastq -a full -t 4
```

This would have resulted in the following outcome:

```
Analysis
│
├── Output-Fastq
│	├── din_9250.trim.fq
│	├── din_9251.trim.fq
│	├── din_9252.trim.fq
│	├── hut_CD15010.trim.fq
│	├── hut_CD15012.trim.fq
│	├── hut_CFS1504g.trim.fq
│	│	
│	├── Populations_r10
│	│	└── outputs...
│	│	
│	├── Populations_r20
│	│	└── outputs...
│	│	
│	├── Populations_r30
│	│	└── outputs...
│	│	
│	├── Populations_r40
│	│	└── outputs...
│	│	
│	├── Populations_r50
│	│	└── outputs...
```

The Stacks workflow was run on the `-o` directory using the fastq files present (no copying took place). However, in this case all the sample names in the output files would contain `.trim` in their label. 


### Running the full Stacks workflow for big datasets <a name=RSRSWB></a>

Running the **ustacks** module sequentially on all the input fastq files can sometimes take a very long time. In general, it is the longest step to run in the Stacks workflow. For large datasets it is much more efficient to split the fastq files into smaller batches, and process them simultaneously. The `Run_Stacks.py` script makes it easy to accomplish this. To run an analysis like this you will first need to create a directory for each batch of fq files. Here is a relatively simple example:

```
Analysis
│
├── Fastq-1
│	 ├── din_9250.trim.fq
│	 └── din_9251.trim.fq
│
├── Fastq-2
│	 ├── din_9252.trim.fq
│	 └── hut_CD15010.trim.fq
│
├── Fastq-3
│	 ├── hut_CFS1504g.trim.fq
│	 └── hut_CG15146.trim.fq
│
├── Stacks-Run
│
```

Then, for each of these directories you can execute the `Run_Stacks.py` script independently to run the **ustacks** component. Here is the command for the first batch:

```
python Run_Stacks.py -i /Analysis/Fastq-1 -o /Analysis/Stacks-Run -a ustacks --intid 1 -t 2 -l run1
```

Although I've left them out here, you can also include the other optional flags that affect ustacks (`-M`, `-m`) for this step. However, two other flags are necessary here: `-l` and `--intid`.

The `-l` flag assigns a unique label to the `Stacks_pipeline.log` file. In this case it would create a log file called `Stacks_pipeline_run1.log`. Since we are going to be writing all outputs to the same directory, this prevents the log files from overwriting or clashing. 

The `--intid` flag is important because every sample must receive a unique integer ID. In a typical run using all the fastq files, the script automatically assigns this number to the samples sequentially, starting with 1. Given that we are running these samples independently, this would mean there would be duplicates for ID 1 and 2. Instead of the automatic assignment of numbers, we'll need to keep track of where to start the next batch. Because there are two samples in the first directory (`Fastq-1`), they will be assigned ID 1 and 2, so we will have to start numbering with 3 for the next directory (`Fastq-2`):

```
python Run_Stacks.py -i /Analysis/Fastq-2 -o /Analysis/Stacks-Run -a ustacks --intid 3 -t 2 -l run2
```

This will assign 3 and 4 to the samples here. And for the last directory we will have to start with and ID integer of 5:

```
python Run_Stacks.py -i /Analysis/Fastq-3 -o /Analysis/Stacks-Run -a ustacks --intid 5 -t 2 -l run3
```

This guarantees that the final ustacks files will all have unique IDs, ranging from 1-6, and that they will be compatible when running the remaining Stacks modules. The ustacks processing would finish and result in the following output files in the `Stacks-Run` directory:

```
Analysis
│
├── Fastq-1
│	 ├── din_9250.trim.fq
│	 └── din_9251.trim.fq
│
├── Fastq-2
│	 ├── din_9252.trim.fq
│	 └── hut_CD15010.trim.fq
│
├── Fastq-3
│	 ├── hut_CFS1504g.trim.fq
│	 └── hut_CG15146.trim.fq
│
├── Stacks-Run
│	 ├── Stacks_pipeline_run1.log
│	 ├── Stacks_pipeline_run2.log
│	 ├── Stacks_pipeline_run3.log
│	 ├── din_9250.fq
│	 ├── din_9250.alleles.tsv
│	 ├── din_9250.snps.tsv
│	 ├── din_9250.tags.tsv
│	 ├── din_9251.fq
│	 ├── din_9251.alleles.tsv
│	 ├── din_9251.snps.tsv
│	 ├── din_9251.tags.tsv
│	 ├── din_9252.fq
│	 ├── din_9252.alleles.tsv
│	 ├── din_9252.snps.tsv
│	 ├── din_9252.fq
│	 ├── hut_CD15010.fq
│	 ├── hut_CD15010.alleles.tsv
│	 ├── hut_CD15010.snps.tsv
│	 ├── hut_CD15010.tags.tsv
│	 ├── hut_CFS1504g.fq
│	 ├── hut_CFS1504g.alleles.tsv
│	 ├── hut_CFS1504g.snps.tsv
│	 ├── hut_CFS1504g.tags.tsv
│	 ├── hut_CG15146.fq
│	 ├── hut_CG15146.alleles.tsv
│	 ├── hut_CG15146.snps.tsv
│	 └── hut_CG15146.tags.tsv
```

Now all the ustacks output files are created, they are compatible, and they are in the same directory (along with copied fastq files). The remaining portion of the Stacks workflow (all steps after ustacks) can then be run using the `-a post-ustacks` option. In this example, we would use the following command:

```
python Run_Stacks.py -o /Analysis/Stacks-Run -a post-ustacks -t 6
```

**Note that in this case the `-i` flag should not be used because the fastq files have already been copied to the `-o` directory in the previous steps.**

And that's it! If all remaining modules run successfully you should end up with a series of populations output directories. In each of the populations output directories, there will be a `populations.haplotypes.tsv` file. These `haplotypes.tsv` files are the main targets of the next major step in the pipeline. You can also include any number of the optional flags to change particular settings in **cstacks** or **populations** when running `-a post-ustacks`.

Note that in this example only a handful of samples were being used. In reality, it may be hundreds of samples. Although this strategy will speed up processing with **ustacks**, there can also be a bottleneck in **cstacks**. Loading all the samples to create the catalog can take a considerable amount of time (and memory). To reduce the time used for **cstacks**, the `--catnum` flag can be supplied with an integer value. The integer specifies the number of samples to use to build the catalog. If you have 300 samples, you can use `--catnum 50` to select 50 samples to build the catalog. The samples are drawn randomly from the `population_map.txt` file (or a user-supplied map file). This will definitely save time, but it could potentially result in the loss of some loci (which are never added to the catalog). The data loss is typically minor, but it depends on the specific characteristics of each dataset and which samples happened to be selected. Regardless, for really big datasets this is strongly recommended.

**Known bug with big datasets:** On many systems there is an upper limit to the number of files that can be open simultaneously. The **gstacks** module has the potential to exceed this limit and crash. A description of the problem, and its solution, are posted on the Stacks user group [here](https://groups.google.com/forum/#!msg/stacks-users/GZqJM_WkMhI/m9Hreg4oBgAJ). In short, to prevent this issue you can type:
```
ulimit -n 4096
```
into the terminal *before* you execute the `Run_Stacks.py` script. If **gstacks** crashes with the error described above, it is quite easy to re-run **gstacks** and subsequently **populations**. They can be run one after another using the `-a` flag with the appropriate module option to finish the analysis. 


### Running the partial Stacks workflow <a name=RSRSWP></a>

For various reasons, you may wish to run a single module in the Stacks workflow. This can be accomplished by specifying the module using the `-a` flag and one of the options (**ustacks**, **cstacks**, **sstacks**, **tsv2bam**, **gstacks**, **populations**). In general, the Stacks workflow generates outputs that are used for the next module. Therefore, any execution of the full or partial Stacks workflow is accompanied by a quick check to make sure that: 1) the requisite files are present, and 2) the output files of the current module do not yet exist. However, the check is not exhaustive. Use discretion when running a particular module.

Unless specified otherwise, the population map is generated automatically and used for any of the steps invoked with the `-a` flag. If the `population_map.txt` is already present, it is removed and created again each time the `Run_Stacks.py` script is used. 

The details for how each of the specific Stacks modules is called is provided below. This can be helpful for understanding how each module is used, especially if you choose to run the workflow in small pieces.


### The details... <a name=RSSWD></a>

There are a lot of internal functions and system calls happening when the `Run_Stacks.py` is run. Some of these details may be of particular interest, or you may be interested in what's happening behind the scenes. Here I provide an explicit walkthrough of the main use-case for the script (e.g., using the `-i` `-o` and `-a full` flags). 

First, the fastq files in the `-i` directory will get copied to the `-o` directory before any other action is taken. They are also relabeled, such that `.trim` is removed from the filename, and just the sample name remains. If you are using your own files (skipping the `Demultiplex_Trim.py` script), their names can be affected if they contain multiple `.` characters. Essentially, the file names are split using the `.` character, and the first component is written with the `.fq` extension. So, `din_9250.trim.fq` would become `din_9250.fq`, `Sample1.Lane2.Proj4.fq` would become `Sample1.fq`, and `MVZ2344XY.fq` would be unaffected and written as `MVZ2344XY.fq`. 

We can expect the `-o` directory to be populated with these relabeled files. In the example below, we set `-i /Analysis/Output-Fastq` and `-o /Analysis/Stacks-Run`, and the following is the expected outcome:

Starting directory structure:
```
Analysis
│
├── Output-Fastq
│	├── din_9250.trim.fq
│	├── din_9251.trim.fq
│	├── din_9252.trim.fq
│	├── hut_CD15010.trim.fq
│	├── hut_CD15012.trim.fq
│	├── hut_CFS1504g.trim.fq
│	└── hut_CG15146.trim.fq
│
├── Stacks-Run
│
```

Directory structure after copying is completed:
```
Analysis
│
├── Output-Fastq
│	├── din_9250.trim.fq
│	├── din_9251.trim.fq
│	├── din_9252.trim.fq
│	├── hut_CD15010.trim.fq
│	├── hut_CD15012.trim.fq
│	├── hut_CFS1504g.trim.fq
│	└── hut_CG15146.trim.fq
│
├── Stacks-Run
│	├── din_9250.fq
│	├── din_9251.fq
│	├── din_9252.fq
│	├── hut_CD15010.fq
│	├── hut_CD15012.fq
│	├── hut_CFS1504g.fq
│	└── hut_CG15146.fq
│
```

After the copying is completed, the following steps will be completed in order using the `-o` directory (`/Analysis/Stacks-Run` in the example above):

**1. ustacks** is run sequentially for every fastq file present. The default command used internally to call ustacks is:
```
ustacks -f <fastq file> -i <automatic> -o <output directory specified> -M 2 -m 3 -p 1
```
The `-i` flag is set automatically, as files will assigned an integer ID in sequential order starting with 1. The ustacks -o flag is the same as the `-o` flag used in the script. The ustacks -M and -m flags have defaults that can be changed using the identical flag names (`-M` and `-m`) in the script. Finally, the ustacks -p flag sets the number of threads, which defaults to 1 unless the `-t` flag is supplied in the script.

After **ustacks** is complete, there will be `sample.alleles.tsv`, `sample.snps.tsv`, and `sample.tags.tsv` files for each fastq file in the starting set. These must be present for the **cstacks** step to run.

**2.** A population map file (`population_map.txt`) is automatically created based on the fastq files present. For example, if the directory contains:
```
├── Stacks-Run
│	├── din_9250.fq
│	├── din_9251.fq
│	├── din_9252.fq
```
then the population map file will contain:

```
din_9250	1
din_9251	1
din_9252	1
```
The `population_map.txt` is used as a convenient way to reference all the samples for all downstream steps, rather than having to specify the samples individually in each module. If the `-p` flag is used with a user-supplied population map file, this step is skipped and that file will be used instead.

**3. cstacks** is run to build the catalog. The default command used internally to call cstacks is: 
```
cstacks -P <output directory specified> -M population_map.txt -p 1 -n 2
```
By default, all samples are added to the catalog using the automatically generated `population_map.txt` file. The `--catnum` flag can be used to specify a particular number of samples to create the catalog instead. They will be randomly sampled from the population map. The cstacks -o flag is the same as the `-o` flag used in the script. The cstacks -n flag can be changed using the identical flag names (`-n`) in the script. The cstacks -p flag sets the number of threads, which defaults to 1 unless the `-t` flag is supplied in the script.

After **cstacks** is complete, there will be a `catalog.alleles.tsv`, `catalog.snps.tsv`, and `catalog.tags.tsv` file present. These must be present for the following step to run.

**4. sstacks** is run to match loci for all samples. The default command used internally to call sstacks is: 
```
sstacks -P <output directory specified> -M population_map.txt -p 1
```
The sstacks -o flag is the same as the `-o` flag used in the script. All samples in the `population_map.txt` file are screened against the catalog using the sstacks -M flag. The sstacks -p flag sets the number of threads, which defaults to 1 unless the `-t` flag is supplied in the script.

After **sstacks** is complete, there will be a `sample.matches.tsv` file for every sample in the `population_map.txt` file. These must be present for the following step to run.

**5. tsv2bam** is run to create properly formatted files for the subsequent step. The default command used internally to call tsv2bam is: 
```
tsv2bam -P <output directory specified> -M population_map.txt -p 1
```
The tsv2bam -o flag is the same as the `-o` flag used in the script. All samples in the `population_map.txt` file are included using the tsv2bam -M flag. The tsv2bam -p flag sets the number of threads, which defaults to 1 unless the `-t` flag is supplied in the script.

After **tsv2bam** is complete, there will be a `sample.matches.bam` file for every sample in the `population_map.txt` file. These must be present for the following step to run.

**6. gstacks** is run to aligns reads and genotype indivuals. The default command used internally to call gstacks is: 
```
gstacks -P <output directory specified> -M population_map.txt -p 1
```
The gstacks -o flag is the same as the `-o` flag used in the script. All samples in the `population_map.txt` file are included using the gstacks -M flag. The gstacks -p flag sets the number of threads, which defaults to 1 unless the `-t` flag is supplied in the script.

After **gstacks** is complete, there will be a `catalog.fa.gz`, `catalog.calls`, and `catalog.log.distribs` file present. These must be present for the following step to run.


**7. populations** is run using several automatically set values for r or R, including 10, 20, 30, 40, 50, 60, 70, 80, 90, and 100. For each of these values, populations is executed and an output directory is created for that run within the `-o` directory. Sticking with the above example directory structure, this would result in:
```
├── Stacks-Run
│	├── din_9250.fq
│	├── din_9251.fq
│	├── din_9252.fq
│	│	
│	├── Populations_r10
│	│	└── outputs...
│	│	
│	├── Populations_r20
│	│	└── outputs...
│	│	
│	├── Populations_r30
│	│	└── outputs...
│	│	
│	├── Populations_r40
│	│	└── outputs...
│	│	
│	├── Populations_r50
│	│	└── outputs...
```
And so on.

The default command will change based on the current r or R value. Here is an example based on the value of 20:
```
populations -P <output directory specified> -O Populations_r20 -M population_map.txt -t 1 -r 20 
```

The output files associated with each run include: `populations.log`, `populations.log.distribs`, `populations.sumstats.tsv`, `populations.hapstats.tsv`, `populations.sumstats_summary.tsv`, and `populations.haplotypes.tsv`. They are relabeled to include the particular r or R value. For example, `populations.haplotypes.tsv` becomes `populations_r20.haplotypes.tsv`. The `populations.haplotypes.tsv` contained in each of the output populations directories are the main target of interest for the downstream steps.

## Outputs <a name=RSO></a>

The main outputs for the next major step are in the directories created from each run of **populations**, as described above. These include the `populations.haplotypes.tsv` files. The outputs from the other individual Stacks modules are described above, but are not used for the next step.

In addition to files generated by Stacks, there is a log file written when the `Run_Stacks.py` script is executed. By default, this file is called `Stacks_pipeline.log`, but this name can be amended by including the `-l` flag with a string. This log file contains information on the settings used (all the flags) to run the script, along with the EXACT commands used to execute each module in Stacks. It also provides the times steps were executed and completed. This information can be valuable for tracking which settings were used, and the relative timing of steps. An example of the contents is provided below:

```
Run executed: 2019-11-20 13:37:22.123817

Stacks Pipeline settings:
-o (outdir): /Volumes/West_Africa/West_Africa/Hyperolius-Merged/4-post-ustacks-L5andOld
-a (analysis): post-ustacks
-i (indir): None
-t (threads): 4
-p (popmap): None
-M (maxdist): 2
-m (mindepth): 3
--intid (id int start): None
-n (nmismatch): 3
--catnum: 50
--maf (minor allele freq): None
--mac (minor allele count): 2
--mintype (r or R): r
--whitelist: None
--blacklist: None
-l (logstring): None

2019-11-20 13:37:24.182604: Executed cstacks: cstacks -s din_CAS254137 -s hut_MSK154 -s hut_UTEP21452 -s hut_UTEP21453 -s hut_UTEP21456 -s hut_UTEP21460 -s hut_UTEP21461 -s hut_UTEP21466 -s hut_YOKO0382 -s hut_ZTNA415 -s oc_GFMJ1074 -s oc_GFMJ1317 -s oc_GFMJ1525 -s oc_GFMJ2185 -s oc_GFMJ2186 -s oc_GFMJ2470 -s oc_GFMJ2649 -s oc_GFMJ983 -s oc_MH0420 -s oc_TKM549 -s occ_CAR023 -s occ_CAS254057 -s occ_CAS254058 -s occ_CD14030 -s occ_CU14893 -s occ_CU15082 -s occ_CU15087 -s occ_CU15205 -s occ_CU15246 -s occ_IRSEN -s occ_JK717 -s occ_MVZ234779 -s occ_MVZ234782 -s occ_NMP6V74700 -s occ_NMP6V775152 -s tub_CAS207703 -s tub_CAS207714 -s tub_CAS207738 -s tub_CAS253398 -s tub_CAS253399 -s tub_CU14912 -s tub_CU14914 -s tub_CU14915 -s tub_CU15144 -s tub_CU15513 -s tub_CU15969 -s tub_GA100 -s tub_GFMJ1076 -s tub_GFMJ1077 -s tub_GFMJ1522 -p 4 -n 3
2019-11-20 14:20:37.581282: Completed cstacks: Elapsed time: 0:43:13.564533  (H:M:S)

2019-11-20 14:20:38.600050: Executed sstacks: sstacks -P /Volumes/West_Africa/West_Africa/Hyperolius-Merged/4-post-ustacks-L5andOld -M /Volumes/West_Africa/West_Africa/Hyperolius-Merged/4-post-ustacks-L5andOld/population_map.txt -p 4
2019-11-20 14:57:16.869404: Completed sstacks: Elapsed time: 0:36:39.249737  (H:M:S)

2019-11-20 14:57:17.716058: Executed tsv2bam: tsv2bam -P /Volumes/West_Africa/West_Africa/Hyperolius-Merged/4-post-ustacks-L5andOld -M /Volumes/West_Africa/West_Africa/Hyperolius-Merged/4-post-ustacks-L5andOld/population_map.txt -t 4
2019-11-20 15:42:06.861269: Completed tsv2bam: Elapsed time: 0:44:49.844840  (H:M:S)

2019-11-20 15:42:11.733182: Executed gstacks: gstacks -P /Volumes/West_Africa/West_Africa/Hyperolius-Merged/4-post-ustacks-L5andOld -M /Volumes/West_Africa/West_Africa/Hyperolius-Merged/4-post-ustacks-L5andOld/population_map.txt -t 4
2019-11-20 15:43:56.254800: Completed gstacks: Elapsed time: 0:01:49.348374  (H:M:S)
```


[Back to top](#TOP)

---------------

# **Filtering tsv files** <a name=FTF></a>

## Contents

[**Overview**](#FTOV)

[**Script arguments**](#FTSA)
	
[**Outputs**](#FTO)

## Overview <a name=FTOV></a>

For this step the `Filter_All_tsv.py` script will be used. After the **populations** module has been run using a variety of r or R values, there are corresponding `populations.haplotypes.tsv` files in each of the resulting populations subdirectories. The purpose of `Filter_All_tsv.py` is to filter all the `populations.haplotypes.tsv` files across these subdirectories to remove invariant, 'blank', and non-biallelic loci, select a single SNP per locus, and filter samples based on a missing data threshold set by the user. 

The SNP selection method is set by the `-s` flag, and allows selection of the first SNP or a random SNP. Missing data levels are calculated on a per sample basis, and a maximum level of percent missing data per sample can be enforced using the `-m` flag. Samples above the missing data threshold are removed. After sample removal, the loci are re-examined to detect if there are new invariant loci. If detected, these invariant loci are discarded. The `--remove_singletons` flag can also be included, and if so it will also cause singletons to be removed at this step. A singleton is when only one sample possesses a variant base for the SNP site (homozygous or heterozygous). Removing singletons appears to reduce errors in model-based clustering methods (Structure, Admixture, etc.), as outlined by [Linck & Battey (2019)](https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.12995).

The `Filter_All_tsv.py` script is intended to run on a directory (`-i`) containing the output of the `Run_Stacks.py` script, in which separate subdirectories were created for each population analysis (one for each r or R value). In other words, the `Filter_All_tsv.py` `-i` directory is the same as the `Run_Stacks.py` `-o` directory. Given the following example:

```
Analysis1
│
├── Stacks-Run
│	├── din_9250.fq
│	├── din_9251.fq
│	├── din_9252.fq
│	│	
│	├── Populations_r10
│	│	└── outputs...
│	│	
│	├── Populations_r20
│	│	└── outputs...
│	│	
│	├── Populations_r30
│	│	└── outputs...
│	│	
│	├── Populations_r40
│	│	└── outputs...
│	│	
│	├── Populations_r50
│	│	└── outputs...
```

The `-i` directory should be specified as: `-i /Anlaysis1/Stacks-Run`. All of the filtering actions described above will be performed for each of the `populations.haplotypes.tsv` files found across the populations subdirectories.

The `Filter_All_tsv.py` script can be run with different missing data values (`-m`) to see the effect on the number of samples and loci retained. The output files from each `-m` value are labeled using the value. This serves to prevent overwriting of output files. Running the next major step ([**Converting filtered tsv files**](#CTF)) will summarize results for all of the filtered datasets. In this regard, it is actually better to run multiple `-m` values for the same input directory (`-i`).

## Script Arguments <a name=FTSA></a>

**Mandatory Arguments:**

+ `-i <path-to-directory>`: The full path to the directory which contains all of the populations output directories produced by `Run_Stacks.py`.

+ `-m <integer>`: The maximum missing data percent to allow in a sample. For example, entering 40 would **include** all samples with **<40%** missing data and **remove** all samples with **>40%** missing data. Setting to 100 will include all samples.

+ `-s <choice>`: Choices = *first*, *random*. Specify whether to choose the first SNP or a random SNP from each available locus.

**Optional Arguments**:

+ `--remove_singletons`: After removal of samples failing missing data threshold, exclude loci with singleton SNPs.


## Outputs <a name=FTO></a>

Portions of the output file names (indicated by `#` below) are filled by the argument values used to run the script, and the other portion (indicated by `XX` below) is filled by the r or R value of the populations directory. This means different values of per sample missing data (`-m`) can be run on the same input directory. Each `-m` value will produce unique output files (rather than causing overwriting).

In each of the populations subdirectory, the output files include:

+ `populations_XX.haplotypes.filtered_#.tsv`: A filtered tsv file. This is the target for the next major step.

+ `Filter_All_tsv.#.log`: A general log file containing the information printed to the screen.

+ `SNP_distributions.#.log`: A SNP site distribution file. This is a record of the number of available SNP sites per locus, across all loci.

+ `Initial_Missing_Data_Per_Sample.#.log`: The calculated missing data level for all samples prior to any filtering steps.

+ `Final_Missing_Data_Per_Sample.#.log`: The calculated missing data level for all samples retained after all filtering steps are completed.

A main summary file (`Filter_All_tsv.summary.#.log`) is written to the input directory (`-i`). This file contains the settings used to run `Filter_All_tsv.py` along with summaries of the filtering steps for each populations tsv file found. Here is an example of the contents:

```
Filter_All_tsv settings:

-i (indir): /Volumes/West_Africa/West_Africa/Hyperolius-Merged/4-post-ustacks-L5andOld
-m (missingdata): 70
-s (snpselection): random

================================================================================


Filtering summary for: populations_r10.haplotypes.tsv

Number of starting loci: 70,676
Number of 'blank' loci: 3,315
Number of invariant loci: 7,544
Number of loci with >2 haplotypes per sample: 0
Number of non-biallelic loci: 0
Remaining loci: 59,817

24 samples PASSED missing data threshold (<70% missing)

232 samples FAILED missing data threshold (>=70% missing)

Number of newly invariant or singleton loci removed: 47,069

Final number of loci: 12,748

------------------------------------------------------------------------------------------


Filtering summary for: populations_r20.haplotypes.tsv

Number of starting loci: 26,279
Number of 'blank' loci: 1,000
Number of invariant loci: 1,210
Number of loci with >2 haplotypes per sample: 0
Number of non-biallelic loci: 0
Remaining loci: 24,069

163 samples PASSED missing data threshold (<70% missing)

93 samples FAILED missing data threshold (>=70% missing)

Number of newly invariant or singleton loci removed: 5,317

Final number of loci: 18,752

------------------------------------------------------------------------------------------

....
```



[Back to top](#TOP)

---------------

# **Converting filtered tsv files** <a name=CTF></a>


## Contents

[**Overview**](#CTOV)

[**Script arguments**](#CTSA)
	
[**Outputs**](#CTO)


## Overview <a name=CTOV></a>

For this step the `Convert_All_tsv.py` script will be used. After the `Filter_All_tsv.py` script has been run once (or multiple times), there will be distinct filtered haplotypes.tsv files present in the populations subdirectories. The purpose of `Convert_All_tsv.py` is to convert all available filtered haplotypes.tsv files into several types of output files. These include the following formats:

+ phylip
+ fasta
+ nexus (two versions: nucleotide and SNAPP)
+ structure
+ ped (two versions: nucleotide and '12' recoded)
+ matrix occupancy 

These formats are described in more detail in the [**Outputs**](#CTO) section. 

The `Convert_All_tsv.py` script is meant to run on a directory (`-i`) containing the outputs of the `Run_Stacks.py` and `Filter_All_tsv.py` scripts. It expects separate subdirectories to be present for each of the populations analyses (one for each r or R value). Each of these subdirectories must contain a **filtered** version of the original `populations.haplotypes.tsv` file. More than one filtered tsv file can be present (resulting from different per sample missing data thresholds `-m`). However, the filtered haplotypes file(s) within each subdirectory must follow this naming scheme: 

`populations_#.haplotypes.filtered_#.tsv`

Here, the `#` indicates additional text which is automatically added by `Filter_All_tsv.py` (for example, `populations_r20.haplotypes.filtered_m50_randomSNP.tsv`). If using your own tsv files, this general naming structure must be intact in order for the files to be read. 

At the end of the script, a main summary file is written to the `-i` directory, called `Convert_All_tsv.summary.txt`. This has a complete summary of every dataset (resulting from every filtered tsv file). This file can be used to compare the datasets and make critical decisions about which one to select for analyses. 


## Script Arguments <a name=CTSA></a>

Running the script is simple. There is only one argument.

**Mandatory Arguments:**

+ `-i <path-to-directory>`: The full path to the directory which contains all of the populations output subdirectories produced by `Run_Stacks.py`. These subdirectories must contain at least one filtered haplotypes.tsv file that was produced using `Filter_All_tsv.py`.


## Outputs <a name=CTO></a>

The output files for each filtered tsv file in a given subdirectory are written to a new directory (`Output-Files`). 

For example, running `Convert_All_tsv.py -i /Analysis1/Stacks-Run` on the following directory:

```
Analysis1
│
├── Stacks-Run
│	├── din_9250.fq
│	├── din_9251.fq
│	├── din_9252.fq
│	│	
│	├── Populations_r10
│	│	└── files...
│	│	
│	├── Populations_r20
│	│	└── files...
│	│	
│	├── Populations_r30
│	│	└── files...
│	│	
│	├── Populations_r40
│	│	└── files...
│	│	
│	├── Populations_r50
│	│	└── files...
```

Would produce output directories like this:

```
Analysis1
│
├── Stacks-Run
│	├── din_9250.fq
│	├── din_9251.fq
│	├── din_9252.fq
│	│	
│	├── Populations_r10
│	│	├── files...
│	│	└── Output-Files
│	│		└── new files...
│	│	
│	├── Populations_r20
│	│	├── files...
│	│	└── Output-Files
│	│		└── new files...
│	│	
│	├── Populations_r30
│	│	├── files...
│	│	└── Output-Files
│	│		└── new files...
│	│	
│	├── Populations_r40
│	│	├── files...
│	│	└── Output-Files
│	│		└── new files...
│	│	
│	├── Populations_r50
│	│	├── files...
│	│	└── Output-Files
│	│		└── new files...
```

For each filtered haplotypes.tsv file found in a given populations subdirectory, the following output files will be written to `Output-Files`:

### phylip

Written as `NAME.phy`. Common file format used in phylogenetics, including [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/index.html). Here, the single SNPs from all loci are concatenated. Consensus SNPs are used for the samples (data are not phased). Example contents:

```
243 371
din_9250 NNGGANTCNTCNNGTANGGGGTNNGNCNARNNCMGNNGAATNNCNNCNCTTGAACANNNAACATCG...
din_9251 TTGGANTCCTCATGTAGNGCGTATGGCNAGNCNCGNAGAATCGCNTCNCTTNANCACAARACATCG...
din_9252 NTGNANTCCTCATGTANNGCGTATGGCNAGCCCCGTAGAANCGCNTNNCTTGANCANAARACATCG...
din_CAS253991 TTGGAATCCYCATGTAGGGCGTATGGCAAGCCCCGTAGAATCGNTTNKCTTGAACACAAAA...
din_CAS254135 TTGGAATCCTCATGTAGGGCGTATGGCAAGNCNMGTNGAATCGCTTCNCTNGAACACAARA...
din_CAS254136 TTGGAATCCTCATGTAGGGSGTATGGCANGNCCCGTAGNATCGCTTCTCTNGAACACAARA...
din_CAS254137 TNGGAATCCYCATGTAGGGCGTATGGCANGCCNCGTAGAATCGCTTCKCTTGAACACAARA...
din_CAS254156 NTGGAATCCTCATGTAGGGSGTATGGCANGCCCCGTAGAATCGCTTCTCTTGAACACAARA...
din_CAS254157 TTGGAATCCNCATGTAGGGSGTATGGCAAGCCCCGTNGAATCGNTTCKCTTGAACACAAAA...
din_CAS254158 TNGGAATCCYCNTGTAGGGCGTATGGCANGCCCCGTNGAATCGCTTNTCTTGANCACAAAA...
din_CAS256693 TNGGAATCCTCATGTAGGGCGTATGGCANGCCCMGTNGANTCGNTTNTCTTGAACACAANA...

```

### fasta

Written as `NAME.fasta`. Common file format for storing sequence data. Here, the single SNPs from all loci are concatenated. Consensus SNPs are used for the samples (data are not phased). Example contents:

```
>din_9250
NNGGANTCNTCNNGTANGGGGTNNGNCNARNNCMGNNGAATNNCNNCNCTTGAACANNNAACATCG...
>din_9251
TTGGANTCCTCATGTAGNGCGTATGGCNAGNCNCGNAGAATCGCNTCNCTTNANCACAARACATCG...
>din_9252
NTGNANTCCTCATGTANNGCGTATGGCNAGCCCCGTAGAANCGCNTNNCTTGANCANAARACATCG...
>din_CAS253991
TTGGAATCCYCATGTAGGGCGTATGGCAAGCCCCGTAGAATCGNTTNKCTTGAACACAAAACANCG...
>din_CAS254135
TTGGAATCCTCATGTAGGGCGTATGGCAAGNCNMGTNGAATCGCTTCNCTNGAACACAARACNTCG...
>din_CAS254136
TTGGAATCCTCATGTAGGGSGTATGGCANGNCCCGTAGNATCGCTTCTCTNGAACACAARACATCG...
```

### nucleotide nexus

Written as `NAME.nex`. Common file format for phylogenetics. Here, the single SNPs from all loci are concatenated. Consensus SNPs are used for the samples (data are not phased). Example contents:

```
#NEXUS 
BEGIN DATA;
	DIMENSIONS  NTAX=243 NCHAR=371;
	FORMAT DATATYPE=DNA  MISSING=N GAP=-;
MATRIX
din_9250 NNGGANTCNTCNNGTANGGGGTNNGNCNARNNCMGNNGAATNNCNNCNCTTGAACANNNAACATCG...
din_9251 TTGGANTCCTCATGTAGNGCGTATGGCNAGNCNCGNAGAATCGCNTCNCTTNANCACAARACATCG...
din_9252 NTGNANTCCTCATGTANNGCGTATGGCNAGCCCCGTAGAANCGCNTNNCTTGANCANAARACATCG...
din_CAS253991 TTGGAATCCYCATGTAGGGCGTATGGCAAGCCCCGTAGAATCGNTTNKCTTGAACACAAAA...
din_CAS254135 TTGGAATCCTCATGTAGGGCGTATGGCAAGNCNMGTNGAATCGCTTCNCTNGAACACAARA...
din_CAS254136 TTGGAATCCTCATGTAGGGSGTATGGCANGNCCCGTAGNATCGCTTCTCTNGAACACAARA...
din_CAS254137 TNGGAATCCYCATGTAGGGCGTATGGCANGCCNCGTAGAATCGCTTCKCTTGAACACAARA...
din_CAS254156 NTGGAATCCTCATGTAGGGSGTATGGCANGCCCCGTAGAATCGCTTCTCTTGAACACAARA...
din_CAS254157 TTGGAATCCNCATGTAGGGSGTATGGCAAGCCCCGTNGAATCGNTTCKCTTGAACACAAAA...
```

### integer nexus

Written as `NAME_SNAPP.nex`. A nexus file that can be imported and used for [SNAPP](https://www.beast2.org/snapp/). Here, instead of DNA nucleotides the characters are integers (0, 1, 2). The characters 0 and 2 represent homozygous SNP sites, and 1 represents heterozygous SNP sites. For this file, the assignment of 0 and 2 is random. That is, A/A is not always coded as 0 across the loci, it could be 0 or 2. This should prevent a known issue with biased assignment in SNAPP. Example contents:

```
#NEXUS 
BEGIN DATA;
	DIMENSIONS  NTAX=105 NCHAR=989;
	FORMAT DATATYPE=INTEGERDATA SYMBOLS="012" GAP=- ;
MATRIX
occ_CAR022 002022020220201022222220211020220-00000220-01002000220-0210000000-0-210...
occ_CAR023 002022020-202000222-2020200020-200-0-00220-020020022200022000001-202210...
occ_CAS207784 -12002000220000022222220--20202--0000-0220-02002002-2-102200-2102202...
occ_CAS207785 01100200-2200000222222202-202020-0000--220002002002020012200-2102202...
occ_CAS207794 -1000200-2201000-22222202-202-1-00000--220002002002-2-00-20002002202...
occ_CAS207795 01-01200-1201000-2222220202020-0-0000-0220-02002002-2-002200-2102202...
occ_CAS207829 012012000120100-222-22202-20-020010000-22000200200202010220002000202...
occ_CAS249970 00202202-2202000222222202000202210000-0220012002-0202000220000000202...
occ_CAS249971 00202-0202202000222--22020002-2200-0-002200020020-2020-0220000000202...
occ_CAS253588 00-122010220200022222220-00020100000000220002002002020-022-0000-0202...
occ_CAS253592 0020220102202000222-2220200-20-0000000022000200200202000220000000202...
occ_CAS254057 012002000220-000-221222020102--000000-0220-02002002-2002220002000202...
occ_CAS254058 011-020002202000222222102020-010-0000-0220-020020020200122-002000202...
```

### structure

Written as `NAME.str`. A file containing phased data for each locus. The structure file produced is compatible with the program [Structure](https://web.stanford.edu/group/pritchardlab/structure.html) and the R package [adegenet](https://github.com/thibautjombart/adegenet/wiki). Example contents:

```
din_9250	-9	-9	4	4	1	-9	2	3	-9	2	3	-9	-9	4	2	1...
din_9250	-9	-9	4	4	1	-9	2	3	-9	2	3	-9	-9	4	2	1...
din_9251	2	2	4	4	1	-9	2	3	3	2	3	1	2	4	2	1...
din_9251	2	2	4	4	1	-9	2	3	3	2	3	1	2	4	2	1...
din_9252	-9	2	4	-9	1	-9	2	3	3	2	3	1	2	4	2	1...
din_9252	-9	2	4	-9	1	-9	2	3	3	2	3	1	2	4	2	1...
din_CAS253991	2	2	4	4	1	1	2	3	3	3	3	1	2	4	2...
din_CAS253991	2	2	4	4	1	1	2	3	3	2	3	1	2	4	2...
din_CAS254135	2	2	4	4	1	1	2	3	3	2	3	1	2	4	2...
din_CAS254135	2	2	4	4	1	1	2	3	3	2	3	1	2	4	2...
din_CAS254136	2	2	4	4	1	1	2	3	3	2	3	1	2	4	2...
```

### standard ped

Written as `NAME.ped`. A ped file in which alleles are represented by nucleotides. Example contents:
```
din_9250 din_9250 0 0 0 0 0 0 0 0 G G G G A A 0 0 T T C C 0 0 T T C C 0 0 0 0 G G...
din_9251 din_9251 0 0 0 0 T T T T G G G G A A 0 0 T T C C C C T T C C A A T T G G...
din_9252 din_9252 0 0 0 0 0 0 T T G G 0 0 A A 0 0 T T C C C C T T C C A A T T G G...
din_CAS253991 din_CAS253991 0 0 0 0 T T T T G G G G A A A A T T C C C C C T C C A...
din_CAS254135 din_CAS254135 0 0 0 0 T T T T G G G G A A A A T T C C C C T T C C A...
din_CAS254136 din_CAS254136 0 0 0 0 T T T T G G G G A A A A T T C C C C T T C C A...
din_CAS254137 din_CAS254137 0 0 0 0 T T 0 0 G G G G A A A A T T C C C C C T C C A...
din_CAS254156 din_CAS254156 0 0 0 0 0 0 T T G G G G A A A A T T C C C C T T C C A...
din_CAS254157 din_CAS254157 0 0 0 0 T T T T G G G G A A A A T T C C C C 0 0 C C A...
din_CAS254158 din_CAS254158 0 0 0 0 T T 0 0 G G G G A A A A T T C C C C C T C C 0...
```

Included is also a **standard map**: Written as `NAME.map`. The map file associated with the standard ped file. It is largely devoid of information, but required by some programs.


### recoded ped

Written as `NAME.ped`. A ped file in which alleles are represented by 2 (major allele) and 1 (minor allele). The 'recoded' ped file is required for the program [Admixture](http://software.genetics.ucla.edu/admixture/). Example contents:
```
din_9250 din_9250 0 0 0 0 0 0 0 0 2 2 2 2 2 2 0 0 2 2 2 2 0 0 2 2 2 2 0 0 0 0 2 2...
din_9251 din_9251 0 0 0 0 2 2 2 2 2 2 2 2 2 2 0 0 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2...
din_9252 din_9252 0 0 0 0 0 0 2 2 2 2 0 0 2 2 0 0 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2...
din_CAS253991 din_CAS253991 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2...
din_CAS254135 din_CAS254135 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2...
din_CAS254136 din_CAS254136 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2...
din_CAS254137 din_CAS254137 0 0 0 0 2 2 0 0 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2...
din_CAS254156 din_CAS254156 0 0 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2...
din_CAS254157 din_CAS254157 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 0 0 2 2 2...
din_CAS254158 din_CAS254158 0 0 0 0 2 2 0 0 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 0...
```

Included is also a **recoded map**: Written as `NAME.map`. The map file associated with the recoded ped file. It is largely devoid of information, but required by some programs.


### matrix occupancy

Written as `NAME.occupancy.csv`. A file that can be used to quickly visualize the data matrix using the tool [here](https://bmedeiros.shinyapps.io/matrix_condenser/). It is just a presence/absence matrix for all the loci per sample. Example contents:
```
0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24...
din_9250,0,0,1,1,1,0,1,1,0,1,1,0,0,1,1,1,0,1,1,1,1,1,0,0,1,0,1,0...
din_9251,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,0...
din_9252,0,1,1,0,1,0,1,1,1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,0...
din_CAS253991,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1...
din_CAS254135,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1...
din_CAS254136,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1...
din_CAS254137,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1...
din_CAS254156,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1...
din_CAS254157,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1...
din_CAS254158,1,0,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1...
```


### Dataset summaries

In addition, a summary file called `Convert_All_tsv.summary.txt` is written to the input directory (`-i`), which contains information on the number of loci and samples contained in every filtered tsv file. This is a complete summary of the final datasets. Here is an example of the contents:

```
tsv_file	Loci	Samples	Total_SNP_sites	Missing_SNP_sites	Perc_Missing_Data
populations_r10.haplotypes.filtered_m30_randomSNP.tsv	0	0	NA	NA	NA
populations_r10.haplotypes.filtered_m50_randomSNP.tsv	0	0	NA	NA	NA
populations_r10.haplotypes.filtered_m70_randomSNP.tsv	31707	48	1521936	956101	62.8
populations_r10.haplotypes.filtered_m90_randomSNP.tsv	45155	119	5373445	3769542	70.2
populations_r20.haplotypes.filtered_m30_randomSNP.tsv	0	0	NA	NA	NA
populations_r20.haplotypes.filtered_m50_randomSNP.tsv	5050	10	50500	18301	36.2
populations_r20.haplotypes.filtered_m70_randomSNP.tsv	26390	104	2744560	1595972	58.2
populations_r20.haplotypes.filtered_m90_randomSNP.tsv	27499	119	3272381	1960144	59.9
populations_r30.haplotypes.filtered_m30_randomSNP.tsv	0	0	NA	NA	NA
populations_r30.haplotypes.filtered_m50_randomSNP.tsv	11881	56	665336	284320	42.7
populations_r30.haplotypes.filtered_m70_randomSNP.tsv	17279	116	2004364	999893	49.9
populations_r30.haplotypes.filtered_m90_randomSNP.tsv	17383	120	2085960	1057607	50.7
populations_r40.haplotypes.filtered_m30_randomSNP.tsv	2346	11	25806	6470	25.1
populations_r40.haplotypes.filtered_m50_randomSNP.tsv	10190	97	988430	379820	38.4
populations_r40.haplotypes.filtered_m70_randomSNP.tsv	10976	119	1306144	538649	41.2
populations_r40.haplotypes.filtered_m90_randomSNP.tsv	10915	120	1309800	544101	41.5
populations_r50.haplotypes.filtered_m30_randomSNP.tsv	3991	42	167622	40737	24.3
populations_r50.haplotypes.filtered_m50_randomSNP.tsv	7155	114	815670	265103	32.5
populations_r50.haplotypes.filtered_m70_randomSNP.tsv	7196	119	856324	284758	33.3
populations_r50.haplotypes.filtered_m90_randomSNP.tsv	7261	120	871320	293377	33.7
```

This file in particular allows a direct comparison of the datasets produced by different r/R values and per-sample missing data thresholds (m). It should be useful in helping to decide which datasets to carry forward for empirical analyses.



[Back to top](#TOP)

-----------

Last updated: November, 2019

For Stacks_pipeline v2.0
