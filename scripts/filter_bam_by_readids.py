import pysam
import sys
import argparse
import os
parser = argparse.ArgumentParser(description='Filter BAM by readids')
parser.add_argument('--inputBAM', dest='inputBAM', type=str, required=True,
                    help='input BAM file')
parser.add_argument('--outputBAM', dest='outputBAM', type=str, required=True,
                    help='filtered output BAM file')
parser.add_argument('--readids', dest='readids', type=str, required=True,
                    help='file with readids to keep (one readid per line)')
args = parser.parse_args()
rids=list(map(lambda x:x.strip(),open(args.readids,'r').readlines()))
inBAM = pysam.AlignmentFile(args.inputBAM, "rb")
outBAM = pysam.AlignmentFile(args.outputBAM, "wb", template=inBAM)
for read in inBAM.fetch():
	if read.query_name in rids:
		outBAM.write(read)
inBAM.close()
outBAM.close()
