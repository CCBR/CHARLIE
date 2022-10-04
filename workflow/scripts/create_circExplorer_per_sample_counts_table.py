import argparse

class BSJ:
    def __init__(self,chrom="",start=-1,end=-1,strand=".",known_novel="novel",read_count=-1,counted=-1):
        self.chrom=chrom
        self.start=start
        self.end=end
        self.strand=strand
        self.known_novel=known_novel
        self.read_count=read_count
        self.counted=counted
    def __str__(self):
        # id="##".join([self.chrom,str(self.start),str(self.end),self.strand])
        return "%s\t%d\t%d\t%s\t%d\t%s\n"%(self.chrom,self.start,self.end,self.strand,self.read_count,self.known_novel)

def read_BSJs(filename,known_novel="novel",counted=-1,threshold=0):
    infile=open(filename,'r')
    BSJdict=dict()
    for l in infile.readlines():
        l=l.strip().split("\t")
        chrom=l[0]
        start=int(l[1])
        end=int(l[2])
        strand=l[5]
        id="##".join([chrom,str(start),str(end)])
        count=int(l[3].split("/")[1])
        if count<=threshold:
            continue
        BSJdict[id]=BSJ(chrom=chrom,start=start,end=end,strand=strand,known_novel=known_novel,read_count=count,counted=counted)
    return(BSJdict)

parser = argparse.ArgumentParser(description='Create CircExplorer2 Per Sample Counts Table')
parser.add_argument('--back_spliced_bed', dest='bsb', type=str, required=True,
                    help='back_spliced.bed')
parser.add_argument('--back_spliced_min_reads', dest='back_spliced_min_reads', type=int, required=True,
                    help='back_spliced minimum read threshold') # in addition to "known" and "low-conf" circRNAs identified by circexplorer, we also include those found in back_spliced.bed file but not classified as known/low-conf only if the number of reads supporting the BSJ call is greater than this number
parser.add_argument('--circularRNA_known', dest='ck', type=str, required=True,
                    help='circularRNA_known.txt')
parser.add_argument('--low_conf', dest='lc', type=str, required=False,
                    help='low_conf.circularRNA_known.txt')
parser.add_argument('-o',dest='outfile',required=True,help='counts TSV table')
args = parser.parse_args()

o=open(args.outfile,'w')
o.write("chrom\tstart\tend\tstrand\tread_count\tknown_novel\n")
all_BSJs=read_BSJs(args.bsb,counted=0,threshold=args.back_spliced_min_reads)
known_BSJs=read_BSJs(args.ck,known_novel="known",counted=0)
if args.lc:
    low_conf_BSJs=read_BSJs(args.lc,known_novel="known",counted=0)
    for k,v in all_BSJs.items():
        if k in low_conf_BSJs:
            all_BSJs[k].known_novel="low_conf"
            all_BSJs[k].strand=v.strand
            all_BSJs[k].counted=1
            low_conf_BSJs[k].counted=1

for k,v in all_BSJs.items():
    if k in known_BSJs:
        all_BSJs[k].known_novel="known"
        all_BSJs[k].strand=v.strand
        all_BSJs[k].counted=1
        known_BSJs[k].counted=1
    o.write(str(all_BSJs[k]))

lst=[known_BSJs]
if args.lc:
    lst.append(low_conf_BSJs)
for l in lst:
    for k,v in l.items():
        if l[k].counted!=1:
            o.write(str(v))
o.close()