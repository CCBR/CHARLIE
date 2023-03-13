import pysam
import sys
import argparse
import gzip
import pprint

def main():
    # debug = True
    debug = False
    parser = argparse.ArgumentParser(
    )
    parser.add_argument("-i","--bedpe",dest="bedpe",required=True,type=str,
        help="Input BEDPE file")
    parser.add_argument("-o","--bed",dest="bed",required=True,type=str,
        help="Output BED file")
    args = parser.parse_args()
    infile = open(args.bedpe,'r')
    outfile = open(args.bed,'w')
    for x in infile.readlines():
        x=x.strip().split("\t")
        chrom=x[0]
        if int(x[1]) < int(x[4]):
            left = x[1]
        else:
            left = x[4]
        if int(x[2]) > int(x[5]):
            right = x[2]
        else:
            right = x[5]
        rid = x[6]
        outfile.write("%s\t%s\t%s\t%s\t0\t.\n"%(chrom,left,right,rid))
        		
    infile.close()
    outfile.close()


if __name__ == "__main__":
    main()