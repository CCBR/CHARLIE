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

# get a list of the chimeric readids
with open(args.junctions, 'r') as junc_f:
    for line in junc_f:
        if "junction_type" in line:
            continue
        readid=line.split()[9]
        rids.append(readid)

rids=list(set(rids))
rids_dict=dict()
for rid in rids:
    rids_dict[rid]=1

print(f"Total chimeric readids:{len(rids)}")

if args.paired:
    for read in inBAM.fetch(until_eof=True):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue
        qname = read.query_name
        if qname in rids:
            continue
        else:
            outBAM.write(read)

#        if not qname in reads:
#            reads[qname]=Read()
#        reads[qname].append_alignment(read)
#        print(qname,len(reads[qname].alignments))
#    output_rids=set(reads.keys())-set(rids)
#    for rid in output_rids:
#        if reads[rid].is_valid_read():
#            for r in reads[rid].alignments:
#                outBAM.write(r)

else:
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
