import pysam
import sys
import argparse
import os
import gzip
# """
# Script takes a BAM file with a list of readids, then 
# filters the input BAM for those readids and outputs
# only those readid alignments into a new BAM file
# @Params:
# @Inputs:
# inputBAM: str (required)
# 	path to input BAM file
# readids: str (required)
# 	path to file with list of readids, one readid per line ... can be gzip-ed
# @Outputs:
# outputBAM: str (required)
# 	path to output BAM file
# """

parser = argparse.ArgumentParser(description='Filter BAM by readids')
parser.add_argument('--inputBAM', dest='inputBAM', type=str, required=True,
                    help='input BAM file')
parser.add_argument('--outputBAM', dest='outputBAM', type=str, required=True,
                    help='filtered output BAM file')
parser.add_argument('--readids', dest='readids', type=str, required=True,
                    help='file with readids to keep (one readid per line)')
args = parser.parse_args()

split_tup = os.path.splitext(args.readids)  
# extract the file name and extension
file_name = split_tup[0]
file_extension = split_tup[1]

rids_dict=dict()
if file_extension==".gz":
	rids=list()
	with gzip.open(args.readids,'rt') as readids:
		for l in readids:
			l = l.strip()
			rids_dict[l]=1
else:
	rids = list(map(lambda x:x.strip(),open(args.readids,'r').readlines()))
	rids = list(set(rids))
	for rid in rids:
		rids_dict[rid]=1
inBAM = pysam.AlignmentFile(args.inputBAM, "rb")
outBAM = pysam.AlignmentFile(args.outputBAM, "wb", template=inBAM)

count=0
for read in inBAM.fetch():
	count+=1
	if count%1000000 == 0:
		print("%d reads read!"%(count))
	qn=read.query_name
	if qn in rids_dict:
		outBAM.write(read)
inBAM.close()
outBAM.close()
