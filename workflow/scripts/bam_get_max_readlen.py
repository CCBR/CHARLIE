import pysam
import argparse


def main():
    # debug = True
    debug = False
    parser = argparse.ArgumentParser(
        description="Print out the maximum aligned read length in the input BAM"
    )
    parser.add_argument("-i","--bam",dest="inbam",required=True,type=str,
        help="Input BAM file")
    args = parser.parse_args()
    samfile = pysam.AlignmentFile(args.inbam, "rb")
    maxrl=0
    for read in samfile.fetch():
        rl = int(read.query_length)
        if rl > maxrl: maxrl=rl
    samfile.close()
    print(maxrl)


if __name__ == "__main__":
    main()