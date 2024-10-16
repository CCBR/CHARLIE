import pysam
import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Remove all non-proper-pair, chimeric, secondary, supplementary, unmapped alignments from input BAM file"
    )
    parser.add_argument("-i","--inbam",dest="inbam",required=True,type=str,
        help="Input  BAM file")
    parser.add_argument("-o","--outbam",dest="outbam",required=True,type=str,
        help="Output primary alignment only BAM file")
    parser.add_argument('-p',"--pe",dest="pe",required=False,action='store_true', default=False,
        help="set this if BAM is paired end")
    args = parser.parse_args()		
    samfile = pysam.AlignmentFile(args.inbam, "rb")
    outfile = pysam.AlignmentFile(args.outbam, "wb", template=samfile)
    for read in samfile.fetch():
        if args.pe and ( read.reference_id != read.next_reference_id ): continue    # only works for PE ... for SE read.next_reference_id is -1
        if args.pe and ( not read.is_proper_pair ): continue
        if read.is_secondary or read.is_supplementary or read.is_unmapped : continue
        outfile.write(read)
    samfile.close()
    outfile.close()


if __name__ == "__main__":
    main()