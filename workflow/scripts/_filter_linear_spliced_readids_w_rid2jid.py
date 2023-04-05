#!/usr/bin/env python3

import argparse
import gzip

def main():
    # debug = True
    debug = False
    parser = argparse.ArgumentParser(
    description=""" Filter read list to only include those that are part of the rid2jid lookup!
    """    
    )
    # INPUTs
    parser.add_argument("--linearin",dest="linearin",required=True,type=str,
        help="gzip-ed input linear readid list")
    parser.add_argument("--splicedin",dest="splicedin",required=True,type=str,
        help="gzip-ed input splicedin readid list")
    parser.add_argument('-r',"--rid2jid",dest="rid2jid",required=True,type=str,
        help="gzip-ed rid2jid lookup")
    # OUTPUTs

    parser.add_argument("--linearout",dest="linearout",required=True,type=str,
        help="gzip-ed output linear readid list")
    parser.add_argument("--splicedout",dest="splicedout",required=True,type=str,
        help="gzip-ed output linear readid list")
    parser.add_argument("--jidcounts",dest="jidcounts",required=True,type=str,
        help="gzip-ed output linear readid list")

    args = parser.parse_args()
# SRR5762377.10004802##-	NC_001806.2##88486##88645##+
# SRR5762377.10008194##+	chrM##1031##1445##+
# SRR5762377.10010198##+	chr45S##8599##9010##+
    linridlist = dict()
    sinridlist = dict()
    with gzip.open(args.linearin,'rt') as inrl:
        for r in inrl:
            r=r.strip()
            linridlist[r]=1
    inrl.close()
    with gzip.open(args.splicedin,'rt') as inrl:
        for r in inrl:
            r=r.strip()
            sinridlist[r]=1
    inrl.close()
    scount=dict()
    lcount=dict()
    with gzip.open(args.rid2jid,'rt') as rid2jid:
        for l in rid2jid:
            l=l.strip().split("\t")
            jid=l[1]
            if jid==".":
                print(">>>>>>>>jid is dot:",l)
            if not jid in scount:
                scount[jid]=dict()
                lcount[jid]=dict()
                scount[jid]["SS"]=0
                scount[jid]["OS"]=0
                scount[jid]["Unknown"]=0
                lcount[jid]["SS"]=0
                lcount[jid]["OS"]=0
                lcount[jid]["Unknown"]=0
            if "##" in l[0]:
                rid,rstrand=l[0].split("##")
            else:
                rid=l[0]
                rstrand="."
            jstrand=jid.split("##")[-1]
            if rstrand=="+" or rstrand=="-":
                if rstrand==jstrand:
                    strandinfo="OS"
                else:
                    strandinfo="SS"
            else:
                strandinfo="Unknown"
            if rid in linridlist:
                linridlist[rid]+=1
                lcount[jid][strandinfo]+=1
            if rid in sinridlist:
                sinridlist[rid]+=1
                scount[jid][strandinfo]+=1
    rid2jid.close()
    with gzip.open(args.linearout,'wt') as outrl:
        for k,v in linridlist.items():
            if v!=1:
                outrl.write("%s\n"%k)
    outrl.close()
    with gzip.open(args.splicedout,'wt') as outrl:
        for k,v in sinridlist.items():
            if v!=1:
                outrl.write("%s\n"%k)
    outrl.close()
    countout=open(args.jidcounts,'w')
    countout.write("#chrom\tstart\tend\tstrand\tlinear_same_strand\tspliced_same_strand\tlinear_opposite_strand\tspliced_opposite_strand\tlinear_unknown_strand\tspliced_unknown_strand\n")
    for k in lcount.keys():
        v1=lcount[k]
        v2=scount[k]
        kstr=k.split("##")
        k="\t".join(kstr)
        countout.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\n"%(k,v1["SS"],v2["SS"],v1["OS"],v2["OS"],v1["Unknown"],v2["Unknown"]))
    countout.close()

if __name__ == "__main__":
    main()
