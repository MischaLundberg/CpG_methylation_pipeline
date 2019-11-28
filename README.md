# CpG_methylation_pipeline

The pipeline can be started at different points to allow to use a beforehand created Dataset for the further evaluation.
The different steps of the pipeline are
- Aligning Bisulfite-Sequencing Data using bwameth
- Gathering the sequence of the selected region
- Analyzing the CpG methylation for each read
- Writing data to an Excel sheet containing chromosom, start, stop, UUID, methylated (yes/no)

---

# Getting Started
In order to download `CpG methylation pipeline`, you should clone this repository via the commands
```  
git clone https://github.com/MischaLundberg/CpG_methylation_pipeline
cd CpG_methylation_pipeline
```

In order to install the Python dependencies, you will need the [Anaconda](https://store.continuum.io/cshop/anaconda/) Python distribution and package manager. After installing Anaconda, run the following commands to create an environment with CpG_Meth's dependencies:

```
conda env create --file environment.yml
conda activate CpGmeth

## to deactivate the environment, type
#conda deactivate
```

In case you are updating your current version of the CpG Methylation pipeline, it would be best practice to also update your environment to the updated prerequisetes.
If you are using a Anaconda environment, you can do so by typing
```
conda env update --name CpGmeth --file environment.yml
```


If you receive any errors while running CpG_Meth, please ensure your versioning for the prerequisites is according to the tested versions.

---

# Input
The needed input depends on the stage you want to run with the pipeline.
Full run:
- Bisulfite sequencing data (fastq format, comma separated)
- Reference FASTA (for the alignment)
- Region to calculate the methylation
- (optional) outputfilename
- [not implemented yet](optional) filter the length of insertions (default is 6000)
.....

```
usage: CpG_Meth.py [-h] -i I [-r R] [-o O] [-f F] [--N_color N_COLOR]
                   [--other_color OTHER_COLOR] [--bam BAM] [--region REGION]
                   [--portrait] [--strict_cpg]
```

Determines the Methylation of given position relative to a reference. 

:information_source: Please remember, if you are not natively using Python 2 to load the environment ```conda activate CpGmeth```

You can start a run like: 
```
/DIRECTORY/CpG_Meth.py -i /FILE_DIRECTORY/simplebs_480.sorted_CpG.bedGraph -r /FILE_DIRECTORY/L1HS.rmsk.txt -o /FILE_DIRECTORY/simplebs_480.bedGraph
```

If you want to start your run with *FASTQ* files, your arguments should be set as following:
```
/DIRECTORY/CpG_Meth.py -i /FILE_DIRECTORY/simplebs_480.1.fastq,/FILE_DIRECTORY/simplebs_480.2.fastq -r /FILE_DIRECTORY/L1HS.rmsk.txt -o /FILE_DIRECTORY/simplebs_480.bedGraph
```


optional arguments are:
```
  -h, --help            show this help message and exit
  -i I                  input might be [.bam file, Bedgraph ,comma separated
                        list of fastq, .txt file (with 0 and 1)]
  -r R                  input referece file, depending on input file, might be
                        bed/fasta/empty e.g. L1HS.bed
  -o O                  output file, e.g. ./meth.xlsx (output file is in Excel
                        format).
  -f F                  for filtering the insertions. e.g. 6000 --> minimum
                        6kb length. ONLY APPLIES WHEN ALIGNING A NEW .BAM
                        FILE!
  --N_color N_COLOR     select the color you want to use for illustrating if
                        the read has at the current position a N,
                        default='green'
  --other_color OTHER_COLOR
                        select the color you want to use for illustrating if
                        the read has at the current position something else
                        than CG, TG or N, default='grey'
  --bam BAM             location of .bam file can be set either here or in
                        args.i
  --region REGION       should be e.g. "chr1:10000-20000"
  --portrait            create figure in portrait mode instead of landscape,
                        default=False
  --strict_cpg          Strict CpG site check of bisulfite sequence for
                        repetitive sequence analysis
```

---

# Prerequisites

Additionally to Python (Anaconda2 or Anaconda3 [if environment is used, as the script is written in python 2!]), the following software is needed:
- BWAMeth can be downloaded [here](https://github.com/brentp/bwa-meth). Nevertheless, a tested version is delivered with this git repository
- BWA MEM is needed by BWAMeth and can be downloaded [here](http://bio-bwa.sourceforge.net/)
- ~~SamTools (> 0.1.19) can be downloaded [here](https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2)~~

If not stated otherwise Install the newest version.
Python 2.7 (if running anaconda, everything except pysam, Matplotlib, XLSXwirter and Toolshed is already preinstalled)
- Pandas (install with pip install pandas)
- Numpy (install with pip install numpy)
- Matplotlib (install with pip install matplotlib)
- Pysam (install with pip install pysam)
- Biopython (install with pip instal biopython)
- XLSXwritert (install with pip install xlsxwriter)
- Toolshed (install with pip install toolshed)

Script was tested with the following versions installed:
 - Anaconda2 
 - Pandas: 0.22.0
 - Numpy: 1.14.3
 - Pysam: 0.14.1
 - Biopython: 1.71
 - Xlsxwriter: 1.0.2
 - toolshed: 0.3.9

---

# Further development
There is going to be a future update containing the options
- to derive only sequences for a given BED file

---

# Contact
For questions or errors, either create a issue/pullrequest or contact me via email (mischa(dot)lundberg(at)mater(dot)uq(dot)edu(dot)au

---

# Dont forget to cite
