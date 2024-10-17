#!/usr/bin/env python3

import argparse
import gzip

def main():
    # debug = True
    debug = False
    parser = argparse.ArgumentParser(
    )
    # INPUTs
    parser.add_argument("-i","--inbed",dest="inbed",required=True,type=str,
        help="Input bamtobed bed file")
    # OUTPUTs
    parser.add_argument('-o',"--outbed",dest="outbed",required=True,type=str,
        help="Output bed file")
    parser.add_argument('-l',"--linear",dest="linear",required=True,type=str,
        help="gzip-ed list of linear readids")
    parser.add_argument('-s',"--spliced",dest="spliced",required=True,type=str,
        help="gzip-ed list of spliced readids")
    args = parser.parse_args()
    outbed = open(args.outbed,'w')
    pairtest = 0
    paired = 0
    readname_counts = dict()
    with open(args.inbed,'r') as inbed:
        for l in inbed:
            l=l.strip().split("\t")
            l1=[]
            l2=[]
            l1.append(l[0])
            l2.append(l[0])
            l1.append(l[1])
            l1.append(l[1])
            l2.append(l[2])
            l2.append(l[2])
            if "/" in l[3]: # paired end
                if pairtest == 0:
                    pairtest=1
                    paired=1
                x=l[3].split("/")
                readname=x[0]
                if x[1]=="1":   # pick the strand of mate1 as the read strand
                    strand=l[5]
                else:           # if it is mate2 then reverse the strand
                    if l[5]=="-": 
                        strand="+"
                    elif l[5]=="+":
                        strand="-"
                    else:       # if neither + or - is provided the use whatever is provided
                        strand=l[5]
            else: # single end
                readname=l[3]
                strand=l[5]
            if readname in readname_counts:
                readname_counts[readname]+=1
            else:
                readname_counts[readname]=1
            readname+="##"+strand
            l1.append(readname)
            l2.append(readname)
            l1.append(".")
            l2.append(".")
            l1.append(strand)
            l2.append(strand)
            outbed.write("\t".join(l1)+"\n")
            outbed.write("\t".join(l2)+"\n")
    inbed.close()
    outbed.close()
    # linear = open(args.linear,'w')
    # spliced = open(args.spliced,'w')
    limit = 1
    if paired==1: limit=2
    with gzip.open(args.spliced,'wt') as spliced:
        with gzip.open(args.linear,'wt') as linear:
            for rid,count in readname_counts.items():
                if count>limit:
                    spliced.write("%s\n"%rid)
                else:
                    linear.write("%s\n"%rid)
    spliced.close()
    linear.close()

if __name__ == "__main__":
    main()
