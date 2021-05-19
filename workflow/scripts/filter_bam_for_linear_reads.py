import pysam
import sys
import argparse
import os

parser = argparse.ArgumentParser(description='Filter BAM to exclude BSJs and other chimeric alignments')
parser.add_argument('--inputBAM', dest='inputBAM', type=str, required=True,
                    help='input BAM file')
parser.add_argument('--outputBAM', dest='outputBAM', type=str, required=True,
                    help='filtered output BAM file')
parser.add_argument('-j',dest='junctions',required=True,help='chimeric junctions file')
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

keep_rids=dict()
count=0
is_paired_bam=0
for read in inBAM.fetch():
    if read.is_unmapped:
        continue
    if read.is_secondary:
        continue
    count+=1
    if count==1:
        if read.is_paired:
            is_paired_bam=1
    if is_paired_bam==1:
        if not read.is_proper_pair:
            continue
    qn=read.query_name
    if qn in rids:
        continue
    if not qn in keep_rids:
        keep_rids[qn]=dict()
        keep_rids[qn]['read1']=""
        keep_rids[qn]['read2']=""
    if read.is_read1:
        keep_rids[qn]['read1']=read
    if read.is_read2:
        keep_rids[qn]['read2']=read
inBAM.close()

#print(rids["SRR1731877.10077876"].keys())
#exit()

for rid in keep_rids.keys():
    if is_paired_bam==1 and keep_rids[rid]['read2']=="":
        continue
    outBAM.write(keep_rids[rid]['read1'])
    if is_paired_bam==1:
        outBAM.write(keep_rids[rid]['read2'])

outBAM.close()

