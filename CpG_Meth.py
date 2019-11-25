#!/usr/bin/python
#-*- coding:utf-8 -*-
################################################
#File Name:CpG_Meth.py                         #
#Python version: 2.7                           #
#Author: Mischa Lundberg                       #
#Mail: mischa.lundberg@mater.uq.edu.au         #
################################################


import argparse
import subprocess
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import pysam
from Bio import SeqIO



def main(args):

    cite = ['Please dont foget to cite: \nhttps://github.com/MischaLundberg/CpG_methylation_pipeline\n','Samtools: https://www.ncbi.nlm.nih.gov/pubmed/19505943','BWAMeth: http://bio-bwa.sourceforge.net/']
    sam = bwameth = 0
    global start
    global end
    global global_read_c
    global global_ref_c
    global nr_of_reads
    start = -1
    end = -1
    global_read_c = 0
    global_ref_c = 0
    nr_of_reads = 0.0
    
    outputFile = args.o
    if outputFile == None:
        outputFile = os.path.splitext(args.i)[0]

    exceloutput = outputFile
    csvoutput = outputFile+".csv"
    if args.region != '':
        regionStart = int(args.region.split(":")[1].split("-")[0])
    else:
        regionStart = 0

    #################################################################
    #   Try to automatically decide where to start in the pipeline  #
    #   by checking the input for a file to start with              #
    #################################################################

    ## if you claim to use preexisting aligned data
    if (".bam" in args.i) or (args.bam is not None):
    
        print "BAM input file supplied"    
        ## if you didnt supply the needed reference fasta
        if ".fasta" not in args.r and ".fa" not in args.r:
            ## TODO: add system.stderr message
            print "!"*47
            print "! Please supply Fasta file (reference Genome) !"
            print "!"*47

        ## if you are using a preexisting .bam file to run you methylation analysis on in
        else:
            if args.bam == None:
                args.bam = args.i
            CpGs = methyl(args)
            CpG_methylation = CpGs[CpGs['chr'].str.endswith("read_c_count")]## get all lines with "read_c_count"
            CpG_methylation.reset_index()
            for index, row in CpG_methylation.iterrows():
                CpG_methylation.loc[index, 'chr'] = str(CpG_methylation.loc[index]['chr'].split('|')[0])
            CpG_methylation.columns = ['chr','start','stop','UUID','[methylated C\'s, C\'s in reference]', 'methylation of this read']
            last_line = pd.DataFrame([['methylation','over','all','reads','together:',(((global_read_c/nr_of_reads)/global_ref_c)*100)]], columns=CpG_methylation.columns)
            CpG_methylation = CpG_methylation.append(last_line, ignore_index=True)
            CpGs = CpGs[~CpGs['chr'].str.endswith("read_c_count")]## get all lines but lines with "read_c_count"
            writer = pd.ExcelWriter(exceloutput+".xlsx")
            CpGs.to_excel(writer,'CpG')
            writer.save()
            writer = pd.ExcelWriter(exceloutput+"_methylation_per_read.xlsx")
            CpG_methylation.to_excel(writer,'CpG methylation')
            writer.save()
            print "*** Excelsheet saved here: %s" %exceloutput+".xlsx"
            CpGs.to_csv(csvoutput, sep=',',header=False,index=False) 
            plot(csvoutput, outputFile+".svg", regionStart, not(args.portrait), args.N_color, args.other_color)

    ## if you use a output file from a prevalent step
    elif ".txt" in args.i or ".csv" in args.i:
        print "Outputfile from previous run supplied. Starting plotting."
        plot(args.i, outputFile+".svg", regionStart, not(args.portrait), args.N_color, args.other_color)

    ## if you want to run the full pipeline
    elif (".fa" in args.r or ".txt" in args.r) and "," in args.i:
        print "Unaligned input supplied, generating BAM using BWAMeth"
        args.bam = make_meth(args)
        sam = 1
        bwameth = 1
        CpGs = methyl(args)
        #print CpGs
        if len(CpGs)==0:
            print "!"*48
            print "! no CpGs found. Please check your input files !"
            print "!"*48
        CpG_methylation = CpGs[CpGs['chr'].str.endswith("read_c_count")]## get all lines with "read_c_count"
        for index, row in CpG_methylation.iterrows():
            CpG_methylation.loc[index, 'chr'] = str(CpG_methylation.loc[index]['chr'].split('|')[0])
        CpG_methylation.columns = ['chr','start','stop','UUID','[methylated C\'s, C\'s in reference]', 'methylation of this read']
        last_line = pd.DataFrame([['methylation','over','all','reads','together:',(((global_read_c/nr_of_reads)/global_ref_c)*100)]], columns=CpG_methylation.columns)
        CpG_methylation = CpG_methylation.append(last_line, ignore_index=True)
        CpG_methylation.reset_index()
        CpGs = CpGs[~CpGs['chr'].str.endswith("read_c_count")]## get all lines but lines with "read_c_count"
        writer = pd.ExcelWriter(exceloutput+".xlsx")
        CpGs.to_excel(writer,'CpG')
        writer.save()
        writer = pd.ExcelWriter(exceloutput+"_methylation.xlsx")
        CpG_methylation.to_excel(writer,'CpG methylation')
        writer.save()
        #print CpG_methylation
        print "*** Excelsheet saved here: %s" %exceloutput+".xlsx"
        CpGs.to_csv(path_or_buf=exceloutput+"_plot.csv", sep=',', header=False,index=False)
        plot(exceloutput+"_plot.csv", outputFile+".svg", regionStart, not(args.portrait), args.N_color, args.other_color)
        os.remove(exceloutput+"_plot.csv")

    elif ".bed" in args.r and ".bed" in args.i:
        print "!"*48
        print "! Run it from a bed is not implemented so far. !"
        print "!"*48
        #TODO: add system.stderr message
        exit(1)
        
    else:
        print "!"*37
        print "! Could not resolve input file type !" 
        print "!"*37
        print ""
        print "your input file was: %s" %(args.i)
        print "Your Input contained several files: %s" %("," in args.i)
        exit(1)

    citation = cite[0]
    if sam == 1:
        citation += cite[1]
    if bwameth == 1:
        citation += cite[2]
    print citation
        
        
##  Creates a Pandas DataFrame containing the parsed
##  input BAM file
def methyl(args):
    

    global global_ref_c
    global nr_of_reads
    #meth_dict = {'Methylation_List','Family','UUID','position'}
    #CG_dict = {'Chr','Start','Stop'}
    print "Loading the bam file"
    samfile = pysam.AlignmentFile(args.bam, "rb")
    ref = SeqIO.to_dict(SeqIO.parse(open(args.r), 'fasta'))
    faChrCheck =  "chr" in SeqIO.parse(open(args.r), 'fasta').next().id 
    header = ["chr","start","stop","UUID","methylated CpGs","position"]
    CpG = pd.DataFrame(columns=header)
    if ":" in args.region and "-" in args.region:
        chrom = args.region.split(":")[0]
        start = int(args.region.split(":")[1].split("-")[0])
        end = int(args.region.split("-")[1])
        if not faChrCheck:
            chrom = args.region.split(":")[0][3:]
        
    else:
        #TODO: add system.stderr message
        print "Either no region given or the region is not given in the needed way chr:start-stop. If it was not your intention, please refer to the argument --region"
        print "We will infer that your reference only contains one element, i.e. an L1 element you want to align agains"
        ## infering chrom and start stop. There should be only one element in the input fasta
        chrom = ref.keys()[0]
        start = 0
        end = len(ref.get(ref.keys()[0]))
        #print "chrom: %s" %chrom
        #print "!"*20
        #print "sequence: %s " %ref.get(chrom)[start:end].seq
        #print "!"*20
        
    ## fill the dataframe with all CpGs from the reference to then check the methylation status of each read for each CpG
    referenceCpGs = initializeCpG(CpG,ref.get(chrom)[start:end].seq.upper(),chrom,start,end)

    ## checking if bam file is indexed by samtools
    if not samfile.check_index():
        print "Input bam file needs to be indexed!"
        print "Please index the bam file using the following command: samtools index "+args.bam
        exit(1)


    global_ref_c = len(referenceCpGs)
    if args.region == "":
        for read in samfile.fetch():
            nr_of_reads += 1
            CpG = updateCpG(read, CpG, faChrCheck, referenceCpGs, start, args.strict_cpg)

    else:
        for read in samfile.fetch(args.region.split(":")[0],int(args.region.split(":")[1].split("-")[0]),int(args.region.split("-")[1])):
            nr_of_reads += 1
            CpG = updateCpG(read, CpG, faChrCheck, referenceCpGs, start, args.strict_cpg)

    CpG = equalizeCpG(CpG)
    return CpG.dropna(thresh=4)

##  Initializes the CpG Dict with all positions of a CG
##  to later on check these positions for their methylation 
##  in the BAM file
def initializeCpG(CpG, referenceSeq, chrom, start, end):


    CpGs = []
    #print "referenceSeq: %s" %referenceSeq
    for element in range(len(referenceSeq)-1):
        if referenceSeq[element] == "C" and referenceSeq[element+1] == "G":
            CpGs.append({"chr":chrom,"start":start+element,"stop":start+element+1,"UUID":np.nan,\
                        "methylated CpGs":-1,"position":start+element})
    #print "CpGs: %s" %CpGs
    return CpGs

##  Updates the given CpG DataFrame for each CG position
##  in the reference file
def updateCpG(read, CpG, faChr, referenceCpGs, regionStart, strict_cpg):

    global start
    global end
    data = []
    current_read_c = 0.0
    current_ref_c = 0.0
    global global_read_c

    ### iterate over all Cs in the referenceCpGs
    position = read.reference_start - regionStart
    chrom = read.reference_name

    #print "referenceCpGs: %s" %referenceCpGs
    for entry in referenceCpGs:#element in range(len(read.seq)-1):
        #### each entry represents one C within the reference
        element = int(entry.get("position"))-regionStart
        methylated = -1
        if (entry.get("position")>=read.reference_start) and (entry.get("position")<=read.reference_end):
            current_ref_c += 1
            if start == -1:
                start = position
            if end < position:
                end = position
            element = int(entry.get("position")-read.reference_start)
            if not read.is_reverse:
                # methylated
                if (strict_cpg and read.seq[element] == "C" and read.seq[element+1] == "G") or (not strict_cpg and read.seq[element] == "C"):
                    methylated = 1 
                    global_read_c += 1
                    current_read_c += 1
                # not methylated
                elif (strict_cpg and read.seq[element] == "T" and read.seq[element+1] == "G") or (not strict_cpg and read.seq[element] == "T"):
                    methylated = 0 
                elif read.seq[element] == "N": 
                    methylated = 2#3
                else:
                    methylated = 3
            else:
                # methylated
                if (strict_cpg and read.seq[element+1] == "C" and read.seq[element] == "G") or (not strict_cpg and read.seq[element+1] == "C"):
                    methylated = 1 
                    global_read_c += 1
                    current_read_c += 1
                # not methylated
                elif (strict_cpg and read.seq[element+1] == "T" and read.seq[element] == "G") or (not strict_cpg and read.seq[element+1] == "T"):
                    methylated = 0#2
                elif read.seq[element+1] == "N": 
                    methylated = 2#3
                else:
                    methylated = 3
            chrom = read.reference_name
            #print "chrom: "+chrom
            if not faChr:
                if len(chrom) <= 2:
                    chrom = "chr"+read.reference_name
            if CpG['UUID'][CpG['chr'] == chrom][CpG['start'] == read.pos].values == []:
                row = CpG.iloc[['chr' == read.reference_name]['start' == read.pos]]
                CpG.at[row,'UUID'] = read.query_name
                CpG.at[row,'methylated CpGs'] = methylated
                CpG.at[row,'position'] = position-start
            else:
                data.append({"chr":chrom,"start":read.reference_start,"stop":read.reference_end,\
                            "UUID":read.query_name,"methylated CpGs":methylated,"position":position-start})
        position += 1
    ##append data with the current methylation rate of the read
    data.append({"chr":str(chrom)+"|read_c_count","start":read.reference_start,\
                "stop":read.reference_end,"UUID":read.query_name,"methylated CpGs":[current_read_c,current_ref_c],\
                "position":((current_read_c/current_ref_c)*100)})
    out = CpG.append(data, ignore_index=True)
    return out

    
def equalizeCpG(CpG):
    
    minPos = CpG['position'].min()
    maxPos = CpG['position'].max()
    for index, item in CpG.iterrows():
        if "read_c_count" not in CpG.loc[index]['chr']:
            CpG.loc[index, 'position'] -= minPos
        
    return CpG

## needs reference (.fa) and fastq as input
### reference has to be indexed by bwameth...
### returns name of BAM file
def make_meth(args):


    if len(args.i.split(",")) == 2 :
       fq1 = args.i.split(",")[0]
       fq2 = args.i.split(",")[1]
    if "fq" in fq1:
        delim = "fq"
    else:
        delim = fq1.split(".")[-1]

    if not os.path.isfile(fq1) or not os.path.isfile(fq2):
    
        print "!"*(37+len(args.i))
        print "! your input files \"%s\" do not exist. !" %(args.i)
        print "!"*(37+len(args.i))
        exit(1) 

    outputSorted = fq1.split(delim)[0]+"sorted.bam"
    outputBam = fq1.split(delim)[0]+".bam"
    outputSam = fq1.split(delim)[0]+".sam"
    ##check if reference is indexed
    if not os.path.isfile(args.r+".bwameth.c2t"):
        command = "bwameth.py index "+args.r
        print "*** Indexing your reference fasta %s" %args.r
        subprocess.call(command, shell=True)
    command = ""

    if ".fa" in args.r and len(args.i.split(",")) == 2 : 
        command += "python bwameth.py --reference "+args.r+" -t "+str(args.t)+" "+fq1+" "+fq2+" > "+outputSam
        command += "; awk 'length($10)>"+args.f+" || $1 ~ /^@/ {print $0}' "+outputSam+"| samtools view -bS > "+outputBam+"; "
        #command += "; awk \'length($10) "+outputSam+" > "+args.f+" | $1 ~ /^@/\' | samtools view -bS > "+outputBam+"; "
        command += "samtools sort "+outputBam+" -o "+outputSorted+"; samtools index "+outputSorted+";"
    if args.region:        
        regionBam = outputSorted.split(".bam")[0]+args.region+".bam"
        print "the full BAM file is located at: "+outputSorted+"\n\
                Creating now a BAM file containing only the calls of the given region, which is located at: "+regionBam
        command += "samtools view -b "+outputSorted+" \""+args.region+"\" > "+regionBam+"; samtools index "+regionBam+"; "
        outputSorted = regionBam   
    
    print "Aligning fastq files"
    print "using the follwoing command for alignment: %s" %(command)
    subprocess.call(command, shell=True)
    print "Alignment of fastq files finished"
    meth = os.path.basename(outputBam)
    meth += "_CpG.bedGraph"
 
    return outputSorted

    
def add_bed(row):


    bed = str(row).split('\n')[1].split(';')[1].split(',')[0].split(':')
    chrom = bed[0]
    bed = bed[1].split('=')
    start = bed[0]
    end = bed[1]
    strand = bed[2]
    row['Chromosome'] = chrom
    row['Start'] = start
    row['End'] = end
    row['Strand'] = strand
    return row


def get_ref_for_pos(ref, pos, inputformat):


    if inputformat == 1:
        location = ref.loc[(ref['genoName'].str.contains(pos[0])) & (ref['genoStart'] <= (int(pos[1]))) & (ref['genoEnd'] >= (int(pos[2]))), 'genoName':'strand']
    else:
        location = ref.loc[(ref['genoName'].str.contains(pos[0])) & (ref['genoStart'] <= (int(pos[1]))) & (ref['genoEnd'] >= (int(pos[2]))), 'bin':'id']
    return location



def plot(plotData, outputFile, regionStart, landscape, N_color, other_color):

    
    
    df = pd.read_csv(plotData, sep=',', names = ["chr","start","stop","UUID","methylated","position"])
    read_grouped = df.groupby('UUID')
    start_grouped = df.groupby('position')
    reads = read_grouped.groups.keys()
    if len(df.index) > 0:

        if landscape:
            f, ax1 = plt.subplots(figsize=(len(start_grouped),len(read_grouped)+0.75))
        else:
            f, ax1 = plt.subplots(figsize=(len(read_grouped),len(start_grouped)))
        bar_width = 0.75
        circles = []
        if landscape:
            ax1.set_xlim((0,len(start_grouped)))
            ax1.set_ylim((0,len(read_grouped)))
        else:
            ax1.set_xlim((0,len(read_grouped)-1))
            ax1.set_ylim((0,len(start_grouped)))
        
        ax1.set_title('CpG methylation')
        ax1.set_xlabel('Bases')
        plt.xticks([])
        plt.yticks([])
        plt.tight_layout()
        offset = 0.492
        emptyCircle = []
        circles = []
        
        if landscape:
            x = 0
        else:
            y = len(start_grouped)

        for position, group in start_grouped:

            if landscape:
                y = 1
            else:
                x = 0
            
            for read in reads:
                circle = plt.Circle((x+offset, y-offset), offset, color='white', linestyle='None', fill=False)
                if group['UUID'].str.contains(read).any():
                    
                    ind = group[group['UUID']==read].index.values.astype(int)[0]
                    methylated = int(group.loc[ind]['methylated'])
                    if methylated == 0: ##not methylated
                        circle = plt.Circle((x+offset, y-offset), offset, color='black', fill=False)
                        label = 'not methylated'
                    elif methylated == 1: ##methylated
                        circle = plt.Circle((x+offset, y-offset), offset, color='black', fill=True)
                        label = 'methylated'
                    elif methylated == 2: ##detected an N
                        circle = plt.Circle((x+offset, y-offset), offset, color=N_color, fill=True)
                        label = 'N'
                    elif methylated == 3: ##something other than CG, TG or N
                        circle = plt.Circle((x+offset, y-offset), offset, color=other_color, fill=False, linestyle='--', hatch='+')
                        label = 'other than CG, TG or N'
                ax1.add_artist(circle)
                ax1.legend([circle], [label],loc='upper right', bbox_to_anchor=(0.7, 1.35), ncol=4)
                if landscape:
                    y += 1
                else:
                    x += 1    
            if landscape:
                x += 1
            else:
                y -= 1
        #ax1.legend(loc='upper right', bbox_to_anchor=(0.75, 1.35), ncol=4, fancybox=True, shadow=True)
        plt.savefig(outputFile, dpi=300)
        print "Output is saved under: %s" %outputFile
        
    
    else:
        #TODO: add stderr Message
        print "#"*80
        print "# no Data to plot. Please check the region of interest for CpGs in the Excel sheet."
        print "#"*80


if __name__ == '__main__':


    parser = argparse.ArgumentParser(description="""Determines the Methylation of given position relative to a 
                                    reference. You can start a run like: 
                                    python /DIRECTORY/methCalc.py 
                                    -i /FILE_DIRECTORY/simplebs_480.sorted_CpG.bedGraph -r /FILE_DIRECTORY/L1HS.rmsk.txt \
                                    -o /FILE_DIRECTORY/simplebs_480.bedGraph\n\n\
\
                                    If you want to start your run with FASTQ files, your arguments should be set as following:\n\
                                    python /DIRECTORY/methCalc.py \
                                    -i /FILE_DIRECTORY/simplebs_480.1.fastq,/FILE_DIRECTORY/simplebs_480.2.fastq \
                                    -r /FILE_DIRECTORY/L1HS.rmsk.txt -o /FILE_DIRECTORY/simplebs_480.bedGraph""")
    parser.add_argument('-i', required=True, help='input might be [.bam file, Bedgraph ,comma separated list of fastq, .txt file (with 0 and 1)]')
    parser.add_argument('-r', required=False, help='input referece file, depending on input file, might be bed/fasta/empty e.g. L1HS.bed')
    parser.add_argument('-o', required=False, help='output file, e.g. ./meth to get the file ./meth.xlsx (output file is in Excel format).')
    parser.add_argument('-t', required=False, default=2, help='Threads to be used for alignment step in BWAmeth. Default = 2 Threads')
    parser.add_argument('-f', required=False, default=6000, help='for filtering the insertions. e.g. 6000 --> \
                        minimum 6kb length. ONLY APPLIES WHEN ALIGNING A NEW .BAM FILE!')
    parser.add_argument('--N_color', required=False, default='green', help='select the color you want to use for \
                        illustrating if the read has at the current position a N, default=\'green\'')
    parser.add_argument('--other_color', required=False, default='grey', help='select the color you want to use \
                        for illustrating if the read has at the current position something \
                        else than CG, TG or N, default=\'grey\'')
    parser.add_argument('--bam', required=False, help='location of .bam file can be set either here or in args.i')
    parser.add_argument('--region', default="", required=False, help='should be e.g. "chr1:10000-20000"')
    parser.add_argument('--portrait', default=False, required=False, action='store_true', help='create figure in \
                        portrait mode instead of landscape, default=False')
    parser.add_argument('--strict_cpg', default=False, required=False, action='store_true', help='Strict CpG site check of bisulfite sequence for repetitive sequence analysis')
#    parser.add_argument('--lookup', required=False, help='needed if -r is a fasta file and you want to create a graph and bedgraph file. e.g. L1HS.bed')

    args = parser.parse_args()

    main(args)
