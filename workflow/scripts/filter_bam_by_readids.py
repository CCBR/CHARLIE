import pysam
import sys
import argparse
import os
# """
# Script takes a BAM file with a list of readids, then 
# filters the input BAM for those readids and outputs
# only those readid alignments into a new BAM file
# @Params:
# @Inputs:
# inputBAM: str (required)
# 	path to input BAM file
# readids: str (required)
# 	path to file with list of readids, one readid per line
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
rids = list(map(lambda x:x.strip(),open(args.readids,'r').readlines()))
rids=list(set(rids))
rids_dict=dict()
for rid in rids:
    rids_dict[rid]=1
inBAM = pysam.AlignmentFile(args.inputBAM, "rb")
outBAM = pysam.AlignmentFile(args.outputBAM, "wb", template=inBAM)
# bigdict = dict()

count=0
for read in inBAM.fetch():
	count+=1
	if count%1000000 == 0:
		print("%d reads read!"%(count))
	qn=read.query_name
	if qn in rids_dict:
		outBAM.write(read)
	# if not qn in bigdict:
	# 	bigdict[qn]=list()
	# bigdict[qn].append(read)
inBAM.close()

# for r in rids_dict:
# 	for read in bigdict[r]:
# 		outBAM.write(read)
outBAM.close()
