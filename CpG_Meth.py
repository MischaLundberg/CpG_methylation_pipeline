#!/usr/bin/python
#-*- coding:utf-8 -*-
################################################
#File Name: methCalc.py
#Author: Mischa Lundberg
#Mail: mischa.lundberg@mater.uq.edu.au
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
from collections import defaultdict
from Bio.Seq import Seq
from Bio import SeqIO



def main(args):

    cite = ['Please dont foget to cite: \nhttps://github.com/MischaLundberg/CpG_methylation_pipeline\n','Samtools: https://www.ncbi.nlm.nih.gov/pubmed/19505943','BWAMeth: http://bio-bwa.sourceforge.net/']
    sam = bwameth = 0
    
    outputFile = args.o
    if outputFile == None:
        outputFile = os.path.splitext(args.i)[0]

    exceloutput = outputFile+".xlsx"
    csvoutput = outputFile+".csv"

    #################################################################
    #   Try to automatically decide where to start in the pipeline  #
    #   by checking the input for a file to start with              #
    #################################################################

    ## if you claim to use preexisting aligned data
    if (".bam" in args.i) or (args.bam is not None):
        
        ## if you didnt supply the needed reference fasta
        if ".fasta" not in args.r and ".fa" not in args.r:
            ## TODO: add system.stderr message
            print "Please supply Fasta file (reference Genome)"

        ## if you are using a preexisting .bam file to run you methylation analysis on in
        else:
            if args.bam == None:
                args.bam = args.i
            CpGs = methyl(args)
            writer = pd.ExcelWriter(exceloutput)
            CpGs.to_excel(writer,'CpG')
            writer.save()
            print "excelsheet saved here: "+exceloutput
            CpGs.to_csv(csvoutput, sep=',',header=False,index=False) 
            plot(csvoutput, outputFile+".svg")

    ## if you use a output file from a prevalent step
    elif ".txt" in args.i or ".csv" in args.i:
        plot(args.i, outputFile+".svg")

    ## if you want to run the full pipeline
    elif ".fa" in args.r and "," in args.i:
        args.bam = make_meth(args)
        sam = 1
        bwameth = 1
        CpGs = methyl(args)
        writer = pd.ExcelWriter(exceloutput)
        CpGs.to_excel(writer,'CpG')
        writer.save()
        print "excelsheet saved here: "+exceloutput
        CpGs.to_csv(file_name=csvoutput, sep=',',header=False,index=False)
        

    elif ".bed" in args.r and ".bed" in args.i:
        print "Run it from bed is not implemented so far."
        #TODO: add system.stderr message
        
    citation = cite[0]
    if sam == 1:
        citation += cite[1]
    if bwameth == 1:
        citation += cite[2]
    print citation
        
##  Creates a Pandas DataFrame containing the parsed
##  input BAM file
def methyl(args):
    

    meth_dict = {'Methylation_List','Family','UUID'}
    CG_dict = {'Chr','Start','Stop'}
    samfile = pysam.AlignmentFile(args.bam, "rb")
    ref = SeqIO.to_dict(SeqIO.parse(open(args.r), 'fasta'))
    faChrCheck =  "chr" in SeqIO.parse(open(args.r), 'fasta').next().id 
    header = ["chr","start","stop","UUID","methylated"]
    CpG = pd.DataFrame(columns=header)
    if ":" in args.region and "-" in args.region:
        chrom = args.region.split(":")[0]
        start = int(args.region.split(":")[1].split("-")[0])
        end = int(args.region.split("-")[1])
        if not faChrCheck:
            chrom = args.region.split(":")[0][3:]
    else:
        #TODO: add system.stderr message
        print "the given region is not given in the needed way chr:start-stop"


#    CpG = initializeCpG(CpG,ref.get(chrom)[start:end],chrom,start,end)

    if args.region == "":
        for read in samfile.fetch():
            CpG = updateCpG(read, CpG, faChrCheck)

    else:
        for read in samfile.fetch(args.region.split(":")[0],int(args.region.split(":")[1].split("-")[0]),int(args.region.split("-")[1])):
            CpG = updateCpG(read, CpG, faChrCheck)


    return CpG.dropna(thresh=4)

##  Initializes the CpG DataFrame with all positions of a CG
##  to later on check these positions for their methylation 
##  in the BAM file
def initializeCpG(CpG, referenceSeq, chrom, start, end):


    data = []
    for element in range(len(referenceSeq)-1):
        if referenceSeq[element] == "C" and referenceSeq[element+1] == "G":
            data.append({"chr":chrom,"start":start+element,"stop":start+element+1,"UUID":np.nan,"methylated":-1})
    CpGs = CpG.append(data, ignore_index=True)
    
    return CpGs

##  Appends the given CpG DataFrame with each CG position
##  in the 
def updateCpG(read, CpG, faChr):

    
    data = []
    for element in range(len(read.seq)-1):
        methylated = 0
        if read.seq[element] == "C" and read.seq[element+1] == "G":
            methylated = 1 
        elif read.is_reverse and read.seq[element] == "G" and read.seq[element+1] == "C":
            methylated = 1
        elif read.seq[element] == "T" and read.seq[element+1] == "G":
            methylated = 1
        chrom = read.reference_name
        if not faChr:
            if len(chrom) <= 2:
                chrom = "chr"+read.reference_name
#        print CpG['UUID'][CpG['chr'] == chrom][CpG['start'] == read.pos].values
#        print CpG.loc[(CpG['chr'] == chrom) & (CpG['start'] == read.pos)]
        if CpG['UUID'][CpG['chr'] == chrom][CpG['start'] == read.pos].values == []:
            row = CpG.iloc[['chr' == read.reference_name]['start' == read.pos]]
            CpG.at[row,'UUID'] = read.query_name
            CpG.at[row,'methylated'] = methylated
        else:
            data.append({"chr":chrom,"start":read.pos,"stop":read.pos+1,"UUID":read.query_name,"methylated":methylated})

    tmp = CpG.append(data, ignore_index=True)
    return tmp

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
    outputSorted = fq1.split(delim)[0]+"sorted.bam"
    outputBam = fq1.split(delim)[0]+".bam"
    ##check if reference is indexed
    if not os.path.isfile(args.r+".bwameth.c2t"):
        command = "bwameth.py index "+args.r
    
    command = ""

    if ".fa" in args.r and len(args.i.split(",")) == 2 : 
        command += "bwameth.py --reference "+args.r+" "+fq1+" "+fq2+" -t "+args.t+" | awk \'length(\$10) > "+args.f+" || \$1 ~ /^@/\' | samtools view -bS > "+outputBam+"; "
        command += "samtools sort "+outputBam+" -o "+outputSorted+"; "
    if args.region:        
        regionBam = outputSorted.split(".bam")[0]+args.region+".bam"
        print "the full BAM file is located at: "+outputSorted+"\nCreating now a BAM file containing only the calls of the given region, which is located at: "+regionBam
        command += "samtools view -b "+outputSorted+" \""+args.region+"\" > "+regionBam+"; "
        outputSorted = regionBam   

    meth = os.path.basename(outputBam)
    meth += "_CpG.bedGraph"
 
    return outputSorted


def build_dict(args):


    a = open(args.i, 'r')
    bed = a.readlines()
    a.close()
    a = open(args.r, 'r')
    ref = a.readlines()
    a.close()

    inputformat = 0

    if int(len(ref[0].split('\n')[0].split('\t'))) == 17:
        if ref[0].split('\n')[0][:3] != "bin":
            refpd = pd.read_csv(args.r, sep='\t', header=None, names=["bin","swScore","milliDiv","milliDel","milliIns","genoName","genoStart","genoEnd","genoLeft",\
            "strand","repName","repClass","repFamily","repStart","repEnd","repLeft","id"])
        else:
            refpd = pd.read_csv(args.r, sep='\t', header=0)

    elif args.r[-3:] == "bed":
        refpd = pd.read_csv(args.r, sep='\t', header=None, index_col=False, names=["genoName","genoStart","genoEnd","repName","swScore","strand"])#,"bin"])
        inputformat = 1
    else:
        print "please use either bed format or the native repeatmasker format."
        exit(0)

    outdict = pd.DataFrame(columns=['ID','Name','Chrom','Start','End','Methylation','Count','Primer_pos','Length_of_finding'])

    for lines in bed:
        if lines[:1] != "t":
            refLine = get_ref_for_pos(refpd, lines.split('\n')[0].split('\t'), inputformat)
            if not refLine.empty:
                ID = refLine.iloc[0]['repName']+";"+str(refLine.iloc[0]['genoName'])+":"+str(refLine.iloc[0]['genoStart'])+"="+str(refLine.iloc[0]['genoEnd'])+"="+str(refLine.iloc[0]['strand'])
                outdict = outdict.append({'ID':ID,'Name':refLine.iloc[0]['repName'],'Chrom':str(refLine.iloc[0]['genoName']),'Start':int(refLine.iloc[0]['genoStart']),\
                'End':int(refLine.iloc[0]['genoEnd']),'Methylation':int(lines.split('\n')[0].split('\t')[3]),'Count':1,\
                'Length_of_finding':(int(lines.split('\n')[0].split('\t')[2])-int(lines.split('\n')[0].split('\t')[1]))}, ignore_index=True)

    tmp = pd.DataFrame(columns=['ID','Methylation'])
    count = outdict.groupby('ID')[['Count']].sum()
    meth = outdict.groupby('ID')[['Methylation']].sum()
    meth2 = outdict.groupby('ID').agg({'Count':'sum','Methylation':'sum'})
    tmp = meth2.apply(lambda row: row['Methylation']/row['Count'], axis=1)
    tmp = tmp.to_frame()
    tmp.rename(columns={0:'Methylation'}, inplace=True)
    tmp = tmp.apply(lambda row: add_bed(row), axis=1)
    tmp = tmp.sort_values(by=['Chromosome','Methylation'], ascending=[True, False])

    if "xlsx" in args.o:
        writer = pd.ExcelWriter(args.o, engine='xlsxwriter')
        tmp.to_excel(writer, sheet_name='Methylation')
        writer.save()
    elif ".bedGraph" in args.o or ".bedgraph" in args.o:
        header = "#\tinput bedGraph: "+args.i+"\n#\tinput reference: "+args.r+"\n#\tChromosome\tStart\tEnd\tMethylation\t"
        header += "\ntrack type=bedGraph name=\"methCalc_to_bedGraph\" description=\"Methylation for all overlapping inputs with referece for the input %s\" visibility=full color=200,100,0 altColor=0,100,200 priority=20" %(args.i)
        body = ""
        for index, row in tmp.iterrows():
            body += "\n"+row['Chromosome']+"\t"+row['Start']+"\t"+row['End']+"\t"+str(row['Methylation'])+"\t"+index+"\t"+str(row['Strand'])
        a = open(args.o, 'w')
        a.write(header)
        a.write(body)
        a.close()
    else:
        #TODO: add system.stderr message
        print ""

    return tmp


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



def plot(plotData, outputFile):

    a = open(plotData, 'r')
    data = a.readlines()
    a.close()
    tmp = []
    for ent in data:
        tmp.append(ent.split("\n")[0].split("\r")[0].split(","))
    data = tmp
    # Create the general block and the "subplots" i.e. the bars
    f, ax1 = plt.subplots(1, figsize=(len(data[0]),len(data)))
    # Set the bar width
    bar_width = 0.75
    x = 0
    y = len(data)
    c = 0
    circles = []
    rows = len(data)
    ax1.set_xlim((0,len(data[0])))
    ax1.set_ylim((0,y))
    plt.xticks([])
    plt.yticks([])
    offset = 0.5
    for row in data:
        x = 0 
        for col in row:
            print "row: "+str(row)+", col: "+col+", x: "+str(x)+", y: "+str(y)
            pos = rows - c
            if col == "0":
                circle = plt.Circle((x+offset, y-offset), offset, color='black', fill=False) #color=colors[c], fill=False)

            else:
                circle = plt.Circle((x+offset, y-offset), offset, color='black', fill=True) #colors[c], fill=True)

            circles.append(circle)
            c += 1
            x += 1

        y -= 1
    
    for ent in circles:
        ax1.add_artist(ent)

    plt.savefig(outputFile, dpi=1000)
    print "Output is saved under: %s" %outputFile



if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Determines the Methylation of given position relative to a reference. You can start a run like: ~/newTools/methCalc.py -i ~/Seth/simplebs_480.sorted_CpG.bedGraph -r ~/Seth/L1HS.rmsk.txt -o ~/01_Team_Research_Data/Seth/simplebs_480.bedGraph')
    parser.add_argument('-i', required=True, help='input might be [.bam file, Bedgraph ,comma separated list of fastq, .txt file (with 0 and 1)]')
    parser.add_argument('-r', required=False, help='input referece file, depending on input file, might be bed/fasta/empty e.g. L1HS.bed')
    parser.add_argument('-o', required=False, help='output file, e.g. ./meth.xlsx (output file is in Excel format).')
    parser.add_argument('-f', required=False, default=6000, help='for filtering the insertions. e.g. 6000 --> minimum 6kb length. ONLY APPLIES WHEN ALIGNING A NEW .BAM FILE!')
    parser.add_argument('--bam', required=False, help='location of .bam file can be set either here or in args.i')
    parser.add_argument('--region', default="", required=False, help='should be e.g. "chr1:10000-20000"')
#    parser.add_argument('--lookup', required=False, help='needed if -r is a fasta file and you want to create a graph and bedgraph file. e.g. L1HS.bed')

    args = parser.parse_args()

    main(args)
