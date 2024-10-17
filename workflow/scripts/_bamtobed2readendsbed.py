#!/usr/bin/env python3

import argparse

def main():
    # debug = True
    debug = False
    parser = argparse.ArgumentParser(
    )
    # INPUTs
    parser.add_argument("-i","--inbed",dest="inbed",required=True,type=str,
        help="Input bamtobed bed file")
    parser.add_argument('-o',"--outbed",dest="outbed",required=True,type=str,
        help="Output bed file")
    args = parser.parse_args()
    outbed = open(args.outbed,'w')
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
            if "/" in l[3]:
                x=l[3].split("/")
                readname=x[0]
                if x[1]=="1":
                    strand=l[5]
                else:
                    if l[5]=="-": 
                        strand="+"
                    elif l[5]=="+":
                        strand="-"
                    else:
                        strand=l[5]
            else:
                strand=l[5]
                readname=l[3]
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

if __name__ == "__main__":
    main()
