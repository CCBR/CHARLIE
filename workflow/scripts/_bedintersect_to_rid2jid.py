import argparse
import sys
import gzip

def main():
    # debug = True
    debug = False
    parser = argparse.ArgumentParser(
    )
    parser.add_argument("-i","--bedinteresection",dest="bedint",required=True,type=argparse.FileType('r'),default=sys.stdin,
        help="Input BED intersection file")
    parser.add_argument("-o","--rid2jid",dest="outtsv",required=True,type=str,
        help="Output tsv... gziped")
    parser.add_argument("-m","--maxdist",dest="maxdist",required=True,type=int,
        help="Max dist from BSJ coordinate")
    args = parser.parse_args()
    # outfile = open(args.outtsv,'w')
    # for l in args.bedint.readlines():
    with gzip.open(args.outtsv,'wt') as outfile:
        for l in args.bedint:
            l=l.strip().split("\t")
            # print(l)
            if ( abs(int(l[2])-int(l[10])) <= args.maxdist ) or ( abs(int(l[1])-int(l[9])) <= args.maxdist ) or ( abs(int(l[2])-int(l[9])) <= args.maxdist ) or ( abs(int(l[1])-int(l[10])) <= args.maxdist ):
                jid=l[0]+"##"+l[1]+"##"+str(int(l[2])-1)
                # outl=l[3:]
                rid=l[11]
                outl=[rid]
                outl.append(jid)
                outstr="\t".join(outl)
                outfile.write("%s\n"%(outstr))
        args.bedint.close()
    outfile.close()

if __name__ == "__main__":
    main()