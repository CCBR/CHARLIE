import pysam
import pprint
import sys
import argparse
import os
from itertools import groupby

def split_text(s):
    for k, g in groupby(s, str.isalpha):
        yield ''.join(g)

def get_alt_cigars(c):
	alt_cigars=[]
	x=list(split_text(c))
	if x[1]=="H":
		alt_cigars.append("".join(x[2:]))
	if x[-1]=="H":
		alt_cigars.append("".join(x[:-2]))
	if x[1]=="H" and x[-1]=="H":
		alt_cigars.append("".join(x[2:-2]))
	return alt_cigars

pp = pprint.PrettyPrinter(indent=4)

parser = argparse.ArgumentParser(description='Filter readid filtered BAM file for BSJ alignments')
parser.add_argument('--inputBAM', dest='inputBAM', type=str, required=True,
                    help='input BAM file')
parser.add_argument('--outputBAM', dest='outputBAM', type=str, required=True,
                    help='filtered output BAM file')
parser.add_argument('--readids', dest='readids', type=str, required=True,
                    help='file with readids to keep (tab-delimited with columns:readid,chrom,strand,site1,site2,cigarlist)')
args = parser.parse_args()
rids=dict()
inBAM = pysam.AlignmentFile(args.inputBAM, "rb")
outBAM = pysam.AlignmentFile(args.outputBAM, "wb", template=inBAM)
count=0
for read in inBAM.fetch():
	count+=1
	qn=read.query_name
	if not qn in rids:
		rids[qn]=dict()
	hi=read.get_tag("HI")
	if not hi in rids[qn]:
		rids[qn][hi]=dict()
		rids[qn][hi]['alignments']=list()
		rids[qn][hi]['sites']=list()
		rids[qn][hi]['cigars']=list()
	rids[qn][hi]['alignments'].append(read)
	site=read.get_reference_positions()[1]
	cigar=read.cigarstring
	cigar=cigar.replace("S","H")
	rids[qn][hi]['sites'].append(site)
	rids[qn][hi]['cigars'].append(cigar)
	rids[qn][hi]['cigars'].sort()	
	#print(site)
	#print(cigar)
	#print(len(rids[qn]))
	#if count==25:
	#	pp.pprint(rids)
	#	for rid in rids:
	#		print(rid,len(rids[rid]))
	#	exit()
inBAM.close()

#print(rids["SRR1731877.10077876"].keys())
#exit()

readidfile = open(args.readids,'r')
readids = readidfile.readlines()
readidfile.close()

#print(rids.keys())

for line in readids:
	line=line.strip().split("\t")
	print(line)
## SRR1731877.10077876	chr16	-	16699504	16700478	30H53M,45M,30M53H
	readid=line[0]
	chrom=line[1]
	strand=line[2]
	site1=line[3]
	site2=line[4]
	cigars=line[5].split(",")
	cigars.sort()
	if strand=="-":
		site=int(site1)+1
	else:
		site=int(site2)+1
	if not readid in rids:
		continue
	for hi in rids[readid].keys():
		print(readid,hi,site)
		print(site,"===>>",rids[readid][hi]['sites'])
		print(site in rids[readid][hi]['sites'])
		print(cigars,"====>>",rids[readid][hi]['cigars'])
		print(rids[readid][hi]['cigars'] == cigars)
		if site in rids[readid][hi]['sites']:
			references=[]
			for read in rids[readid][hi]['alignments']:
				references.append(read.reference_name)
			if len(list(set(references)))!=1: # same HI but different aligning to different chromosomes
				continue
			rids[readid][hi]['alignments']=list(set(rids[readid][hi]['alignments']))
			if rids[readid][hi]['cigars'] == cigars:
				for read in rids[readid][hi]['alignments']:
					outBAM.write(read)	
			else:
				aminusb=list(set(rids[readid][hi]['cigars'])-set(cigars))
				if len(aminusb)==1:
					restcigars=list(set(rids[readid][hi]['cigars'])-set(aminusb))
					altcigars=get_alt_cigars(aminusb[0])
					for ac in altcigars:
						newcigars=[]
						newcigars.extend(restcigars)
						newcigars.append(ac)
						newcigars.sort()
						if newcigars == cigars:
							for read in rids[readid][hi]['alignments']:
								outBAM.write(read)
							break
				if len(aminusb)==2:
					commoncigar=list(set(rids[readid][hi]['cigars'])-set(aminusb))
					altcigars1=get_alt_cigars(aminusb[0])
					altcigars2=get_alt_cigars(aminusb[1])
					found=0
					for ac1 in altcigars1:
						if found!=0:
							break
						tmpcigars=[]
						tmpcigars.extend(commoncigar)
						tmpcigars.append(ac1)
						for ac2 in altcigars2:
							newcigars=[]
							newcigars.extend(tmpcigars)
							newcigars.append(ac2)	
							newcigars.sort()
							if newcigars == cigars:
								for read in rids[readid][hi]['alignments']:
									outBAM.write(read)
								found=1
								break
outBAM.close()
