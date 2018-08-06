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
- (optional) filter the length of insertions (default is 6000)
.....

# Outlook
There is going to be a future update containing the options
- to derive only sequences for a given BED file
- generating a Excel sheet showing the overall methylation of each sequence/insertion (if used for analyzing e.g. SVA/ALU/L1/LTR.... sequences)

# Contact
For Questions or errors, either create a issue/pullrequest or contact me via email (mischa(dot)lundberg(at)mater(dot)uq(dot)edu(dot)au

# Dont forget to cite
