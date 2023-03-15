import argparse
import sys
import pandas

def main():
    # debug = True
    debug = False
    parser = argparse.ArgumentParser(
    )
    parser.add_argument("-b","--bsjcounts",dest="bsjcounts",required=True,type=str,
        help="BSJ counts file")
    parser.add_argument("-l","--linearcounts",dest="linearcounts",required=True,type=str,
        help="Linear counts file")
    parser.add_argument("-o","--mergedcounts",dest="mergedcounts",required=True,type=str,
        help="merged counts file")
    args = parser.parse_args()

    bcounts = pandas.read_csv(args.bsjcounts,header=0,sep="\t")
    lcounts = pandas.read_csv(args.linearcounts,header=0,sep="\t")
    print(bcounts.head())
    print(lcounts.head())
    mcounts = bcounts.merge(lcounts,how='outer',on=["#chrom","start","end","strand"])

    mcounts.to_csv(args.mergedcounts,index=False,doublequote=False,sep="\t")

if __name__ == "__main__":
    main()