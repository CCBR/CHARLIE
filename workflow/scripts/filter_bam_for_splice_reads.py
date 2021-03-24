import pysam
import sys
import argparse
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
    for read in inbam.fetch(c,s,e):
        cigar=read.cigarstring
        cigar=cigar.replace("S","H")
        cigart=read.cigartuples
        if 3 in list(map(lambda z:z[0],cigart)):
            cigart=cigart[list(map(lambda z:z[0],cigart)).index(0):]
            if cigart[0][0]==0 and cigart[1][0]==3:
                start=read.reference_start+cigart[0][1]+1
                end=start+cigart[1][1]-1
                if start==s and end==e:
                    outbam.write(read)
                    #print(read.query_name,c,s,e,start,end)
                    #print(read)
inbam.close()
outbam.close()
exit()