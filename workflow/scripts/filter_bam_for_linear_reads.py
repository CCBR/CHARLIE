import pysam
import sys
import argparse
import os
from collections import defaultdict

            
class Read:
    def __init__(self):
        self.alignments=list()
        self.read1exists=False
        self.read2exists=False
        
    def append_alignment(self,alignment):
        self.alignments.append(alignment)
        if alignment.is_read1:
            self.read1exists=True
        if alignment.is_read2:
            self.read2exists=True
    
    def is_valid_read(self):
        return(self.read1exists and self.read2exists)
        
        


parser = argparse.ArgumentParser(description='Filter BAM to exclude BSJs and other chimeric alignments')
parser.add_argument('--inputBAM', dest='inputBAM', type=str, required=True,
                    help='input BAM file')
parser.add_argument('--outputBAM', dest='outputBAM', type=str, required=True,
                    help='filtered output BAM file')
parser.add_argument('-j',dest='junctions',required=True,help='chimeric junctions file')
parser.add_argument('-p',dest='paired', help='bam is paired', action='store_true')
args = parser.parse_args()
rids=list()
inBAM = pysam.AlignmentFile(args.inputBAM, "rb")
outBAM = pysam.AlignmentFile(args.outputBAM, "wb", template=inBAM)
reads = defaultdict(lambda: Read())
# reads = dict()

# Chimeric junction file columns:
# 1	chr_donorA
# 2	brkpt_donorA
# 3	strand_donorA
# 4	chr_acceptorB
# 5	brkpt_acceptorB
# 6	strand_acceptorB
# 7	junction_type
# 8	repeat_left_lenA
# 9	repeat_right_lenB
# 10	read_name
# 11	start_alnA
# 12	cigar_alnA
# 13	start_alnB
# 14	cigar_alnB
# 15	num_chim_aln
# 16	max_poss_aln_score
# 17	non_chim_aln_score
# 18	this_chim_aln_score
# 19	bestall_chim_aln_score
# 20	PEmerged_bool
# 21	readgrp
# from STAR manual
# The first 9 columns give information about the chimeric junction:
# column 1: chromosome of the donor
# column 2: first base of the intron of the donor (1-based)
# column 3: strand of the donor
# column 4: chromosome of the acceptor
# column 5: first base of the intron of the acceptor (1-based)
# column 6: strand of the acceptor
# column 7: junction type: -1=encompassing junction (between the mates), 1=GT/AG, 2=CT/AC
# column 8: repeat length to the left of the junction
# column 9: repeat length to the right of the junction
# Columns 10-14 describe the alignments of the two chimeric segments, it is SAM like. Alignments
# are given with respect to the (+) strand
# column 10: read name
# column 11: first base of the first segment (on the + strand)
# column 12: CIGAR of the first segment
# column 13: first base of the second segment
# column 14: CIGAR of the second segment


# get a list of the chimeric readids
with open(args.junctions, 'r') as junc_f:
    for line in junc_f:
        if "junction_type" in line:
            continue
        readid=line.split()[9] # 10th column is read-name
        rids.append(readid)


# rids=list(set(rids))
# convert to dict to speed things up
rids_dict=dict()
for rid in rids:
    rids_dict[rid]=1

print(f"Total chimeric readids:{len(rids)}")

if args.paired: # paired-end
    for read in inBAM.fetch(until_eof=True):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue
        qname = read.query_name
        if qname in rids_dict: # "in" dict is much faster than "in" list
            continue # if readid is in dict then it is a junction read ... so ignore it!
        else:
            outBAM.write(read)
else: # single-end
    incount=0
    outcount=0
    for read in inBAM.fetch(until_eof=True):
        incount+=1
        if incount%1000==0:
            print(f"{incount/1000000:.4f}m reads read in")
            print(f"{outcount/1000000:.4f}m reads written out")
        if read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue
        qname = read.query_name
        if qname in rids_dict:
            continue
        else:
            outcount+=1
            outBAM.write(read)



inBAM.close()

outBAM.close()
