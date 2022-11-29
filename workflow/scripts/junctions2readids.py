import argparse
from itertools import groupby

# chimeric junction file columns
# column 1: chr donorA : chromosome of the donor
# column 2: brkpt donorA : first base of the intron of the donor (1-based)
# column 3: strand donorA : strand of the donor
# column 4: chr acceptorB : chromosome of the acceptor
# column 5: brkpt acceptorB : first base of the intron of the acceptor (1-based)
# column 6: strand acceptorB : strand of the acceptor
# column 7: junction type : -1=encompassing junction (between the mates), 1=GT/AG,2=CT/AC
# column 8: repeat left lenA : repeat length to the left of the junction
# column 9: repeat right lenB : repeat length to the right of the junction
######
# Columns 10-14 describe the alignments of the two chimeric segments, it is SAM like. Alignments are given with respect to the (+) strand
######
# column 10: read name : name of the RNA-seq fragment
# column 11: start alnA : first base of the first segment (on the + strand)
# column 12: cigar alnA : CIGAR of the first segment
# column 13: start alnB : first base of the second segment
# column 14: cigar alnB : CIGAR of the second segment
######
# Columns 15-20 provide alignment score information and relevant metadata. These columns are only output for multimapping chimeriuc algorithm --chimMultimapNmax >0.
######
# column 15: num chim aln : number of sufficiently scoring chimeric alignments reported for this
# RNA-seq fragment.
# column 16: max poss aln score : maximum possible alignment score for this fragment’s read(s).
# column 17: non chim aln score : best non-chimeric alignment score
# column 18: this chim aln score : score for this individual chimeric alignment
# column 19: bestall chim aln score : the highest chimeric alignment score encountered for this RNA-seq fragment among the num chim aln reported chimeric alignments.
# column 20: PEmerged bool : boolean indicating that overlapping PE reads were first merge into a single contiguous sequence before alignment.
# column 21: readgrp : read group assignment for the read as indicated in the BAM file

# output is tab-delimited file with the following columns:
# a. readids
# b. chromosome
# c. strand
# d. site1
# e. site2
# f. list of cigars comma-separated (soft-clips are converted to hard-clips)

def split_text(s):
    for k, g in groupby(s, str.isalpha):
        yield ''.join(g)

def split_cigar(c):
	cigars=[]
	if 'p' in c:
		x=list(split_text(c))
		cigars.append(''.join(x[:x.index('p')-1]).replace('S','H'))
		cigars.append(''.join(x[x.index('p')+1:]).replace('S','H'))
	else:
		cigars.append(c.replace('S','H'))
	return cigars

def get_cigars(l):
	cigars=[]
	cigars.extend(split_cigar(l.split()[11]))
	cigars.extend(split_cigar(l.split()[13]))
	cigars=list(filter(lambda x:x!='',cigars))
	return cigars
	
parser = argparse.ArgumentParser(description="""
Extract readids,strand,site,cigar etc. of BSJ reads from chimeric junctions file generated using STAR.
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
        print("\t".join([readid,chr1,strand1,site1,site2,",".join(get_cigars(line))]))
