import pysam
import sys
import argparse
# """
# Script takes a STAR 2p BAM file and tab-delimited file with splice junctions in the first 3 columns,
# and outputs spliced-only alignments
# @Params:
# @Inputs:
# inbam: str (required)
# 	path to input BAM file
# tab: str (required)
# 	path to tab file with splice junctions in the first 3 columns ... this is typically output from STAR 1p after applying filters.
# @Outputs:
# outbam: str (required)
# 	path to output BAM file
# """

parser = argparse.ArgumentParser(description='extract spliced reads from bam file')
parser.add_argument('--inbam',dest='inbam',required=True,help='STAR bam file with index')
parser.add_argument('--tab',dest='tab',required=True,help='tab file with splice junctions in the first 3 columns')
parser.add_argument('--outbam',dest='outbam',required=True,help='Output bam filename')
args=parser.parse_args()

inbam = pysam.AlignmentFile(args.inbam, "rb" )
outbam = pysam.AlignmentFile(args.outbam, "wb", template=inbam )
tab = open(args.tab)
junctions = tab.readlines()
tab.close()
count=0
threshold=0
incr=5
for l in junctions:
    count+=1
    if count*100/len(junctions)>threshold:
        print("%d %% complete!"% (threshold))
        threshold+=incr
    l=l.strip().split("\t")
    c=l[0]
    s=int(l[1])
    e=int(l[2])
# get chromosome name, start and end positions for the junction
# and fetch reads aligning to this region using "fetch"
# ref: https://pysam.readthedocs.io/en/latest/api.html#pysam.FastaFile.fetch
    for read in inbam.fetch(c,s,e):
# get cigarstring to replace softclips
        cigar=read.cigarstring
# replace softclips with hardclip
        cigar=cigar.replace("S","H")
        cigart=read.cigartuples

# if cigartuple contains
# N	BAM_CREF_SKIP	3
# then it is a spliced read!

# ref: https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples
# cigartuples operation list is
# M	BAM_CMATCH	0
# I	BAM_CINS	1
# D	BAM_CDEL	2
# N	BAM_CREF_SKIP	3
# S	BAM_CSOFT_CLIP	4
# H	BAM_CHARD_CLIP	5
# P	BAM_CPAD	6
# =	BAM_CEQUAL	7
# X	BAM_CDIFF	8
# B	BAM_CBACK	9
# cigartuples returns a list of tuples of (operation, length)
# eg. 30M is returned as [(0, 30)]
# N in CIGAR score is index 3 in tuple represents BAM_CREF_SKIP indicative of spliced read
        if 3 in list(map(lambda z:z[0],cigart)):
# cigart[list(map(lambda z:z[0],cigart)).index(0):]
# get the first item of each tuple in the list of tuples
# first item will be the operation from the above table
# if 3 is in the new list ... means that there was a BAM_CREF_SKIP
# BAM_CREF_SKIP in CIGAR score right after a match (BAM_CMATCH)
# suggests spliced alignment aka spliced read
            cigart=cigart[list(map(lambda z:z[0],cigart)).index(0):]
            if cigart[0][0]==0 and cigart[1][0]==3:
# CIGAR has match ... followed by skip ... aka spliced read
# so gather start and end coordinates
                start=read.reference_start+cigart[0][1]+1
                end=start+cigart[1][1]-1
                if start==s and end==e:
# check if start and end are in the junctions file
# if yes then write to output file
                    outbam.write(read)
                    #print(read.query_name,c,s,e,start,end)
                    #print(read)
inbam.close()
outbam.close()
exit()