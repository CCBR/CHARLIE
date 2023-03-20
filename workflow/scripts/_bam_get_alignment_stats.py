#!/usr/bin/env python3
import argparse
import pysam
    
def read_regions(regionsfile):
    infile=open(regionsfile,'r')
    regions=dict()
    for l in infile.readlines():
        l = l.strip().split("\t")
        region_name=l[0]
        regions[region_name]=dict()
        regions[region_name]['sequences']=dict()
        sequence_names=l[1].split()
        for s in sequence_names:
            regions[region_name]['sequences'][s]=dict()
    return regions        


def main():
    parser = argparse.ArgumentParser(description='Find BAM alignment stats for each region.')
    parser.add_argument('--inbam', dest='inbam', type=str, required=True,
        help='Input BAM file')
    parser.add_argument('--regions', dest='regions', type=str, required=True,
        help='regions file eg. ref.fa.regions')
    # parser.add_argument("--out",dest="outjson",required=True,type=str,
    #     help="Output stats in JSON format")
    parser.add_argument('-p',"--pe",dest="pe",required=False,action='store_true', default=False,
        help="set this if BAM is paired end")
    args = parser.parse_args()		
    samfile = pysam.AlignmentFile(args.inbam, "rb")
    regions = read_regions(regionsfile=args.regions)
    region_names = regions.keys()
    for read in samfile.fetch():
        if args.pe and ( read.reference_id != read.next_reference_id ): continue    # only works for PE ... for SE read.next_reference_id is -1
        if args.pe and ( not read.is_proper_pair ): continue
        if read.is_secondary or read.is_supplementary or read.is_unmapped : continue
        rid = read.query_name
        refname = samfile.get_reference_name(read.reference_id)
        for region in region_names:
            if refname in regions[region]['sequences']:
                regions[region]['sequences'][refname][rid]=1
                break
    samfile.close()
    for region in regions:
        counts=0
        for refname in regions[region]['sequences'].keys():
            counts += len(regions[region]['sequences'][refname])
        print("%d\t%s"%(counts,region))


if __name__ == "__main__":
    main()