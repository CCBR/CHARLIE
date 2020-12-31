import argparse

parser = argparse.ArgumentParser(description="""
Extract readids of BSJ reads from chimeric junctions file generated using STAR.
This script outputs readids to the screen which may be redundant, hence piping 
the output through `sort` and then `uniq` is highly recommended
""")
parser.add_argument('-j',dest='junctions',required=True,help='chimeric junctions file')
# parser.add_argument('-r',dest='readids',required=True,help='Output txt file with a readid per line')
args = parser.parse_args()
# ofile=open(args.readids,'w')
with open(args.junctions, 'r') as junc_f:
    for line in junc_f:
        if "junction_type" in line:
            continue
        flag = int(line.split()[6])
        if flag < 0:
            continue
        chr1, site1, strand1, chr2, site2, strand2 = line.split()[:6]
        if chr1 != chr2 or strand1 != strand2:
            continue
        if strand1 == '+':
            start = int(site2)
            end = int(site1) - 1
        else:
            start = int(site1)
            end = int(site2) - 1
        if start > end:
            continue
        readid=line.split()[9]
        print(readid)
#         ofile.write("%s\n"%(readid))
# ofile.close()