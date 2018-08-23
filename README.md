# CpG_methylation_pipeline

The pipeline can be started at different points to allow to use a beforehand created Dataset for the further evaluation.
The different steps of the pipeline are
- Aligning Bisulfite-Sequencing Data using bwameth
- Gathering the sequence of the selected region
- Analyzing the CpG methylation for each read
- Writing data to an Excel sheet containing chromosom, start, stop, UUID, methylated (yes/no)

# Input
The needed input depends on the stage you want to run with the pipeline.
Full run:
- Bisulfite sequencing data (fastq format, comma separated)
- Reference FASTA (for the alignment)
- Region to calculate the methylation
- (optional) outputfilename
- [not implemented yet](optional) filter the length of insertions (default is 6000)
.....

# Prerequisites
If not stated otherwise Install the newest version.
Python 2.7 (if running anaconda, everything except pysam and biopython is already preinstalled)
- pandas (install with pip install pandas)
- numpy (install with pip install pnumpy)
- matplotlib (install with pip install matplotlib)
- pysam (install with pip install pysam)
- Biopython (install with pip instal biopython)

# Further development
There is going to be a future update containing the options
- to derive only sequences for a given BED file

# Contact
For questions or errors, either create a issue/pullrequest or contact me via email (mischa(dot)lundberg(at)mater(dot)uq(dot)edu(dot)au

# Dont forget to cite
