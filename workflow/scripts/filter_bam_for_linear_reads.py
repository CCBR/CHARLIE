import pysam
import sys
import argparse
import os
from collections import defaultdict


def read_pair_generator(bam, rids):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(until_eof=True):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue
        qname = read.query_name
        if qname in rids:
            continue
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            elif read.is_read2:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]


parser = argparse.ArgumentParser(description='Filter BAM to exclude BSJs and other chimeric alignments')
parser.add_argument('--inputBAM', dest='inputBAM', type=str, required=True,
                    help='input BAM file')
parser.add_argument('--outputBAM', dest='outputBAM', type=str, required=True,
                    help='filtered output BAM file')
parser.add_argument('-j',dest='junctions',required=True,help='chimeric junctions file')
parser.add_argument('-p',dest='paired', help='bam is paired' action='store_true')
args = parser.parse_args()
rids=list()
inBAM = pysam.AlignmentFile(args.inputBAM, "rb")
outBAM = pysam.AlignmentFile(args.outputBAM, "wb", template=inBAM)

# get a list of the chimeric readids
with open(args.junctions, 'r') as junc_f:
    for line in junc_f:
        if "junction_type" in line:
            continue
        readid=line.split()[9]
        rids.append(readid)

if args.paired:
    for read1, read2 in read_pair_generator(inBAM,rids):
        outBAM.write(read1)
        outBAM.write(read2)
else:
    count=0
    for read in inBAM.fetch():
        if read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue
        qn=read.query_name
        if qn in rids:
            continue
        count+=1
        if count==1 and read.is_paired:
            exit("BAM file is paired. Use -p option!!")
        outBAM.write(read)
inBAM.close()
outBAM.close()

