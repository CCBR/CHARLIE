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
    
def read_regions(regionsfile,host,additives,viruses):
    host=host.split(",")
    additives=additives.split(",")
    viruses=viruses.split(",")
    infile=open(regionsfile,'r')
    regions=dict()
    for l in infile.readlines():
        l = l.strip().split("\t")
        region_name=l[0]
        regions[region_name]=dict()
        regions[region_name]['sequences']=dict()
        if region_name in host:
            regions[region_name]['host_additive_virus']="host"
        elif region_name in additives:
            regions[region_name]['host_additive_virus']="additive"
        elif region_name in viruses:
            regions[region_name]['host_additive_virus']="virus"
        else:
            exit("%s has unknown region. Its not a host or a additive or a virus!!")
        sequence_names=l[1].split()
        for s in sequence_names:
            regions[region_name]['sequences'][s]=1
    return regions        

def _get_host_additive_virus(regions,seqname):
    for k,v in regions.items():
        if seqname in v['sequences']:
            return v['host_additive_virus']
    else:
        exit("Sequence: %s does not have a region."%(seqname))

def read_BSJs(filename,regions,host_min,host_max,virus_min,virus_max,known_novel="novel",counted=-1,threshold=0):
    infile=open(filename,'r')
    BSJdict=dict()
    for l in infile.readlines():
        l=l.strip().split("\t")
        chrom=l[0]
        start=int(l[1])
        end=int(l[2])
        strand=l[5]
        circid="##".join([chrom,str(start),str(end)])
        count=int(l[3].split("/")[1])
        if count < threshold:
            continue
        host_additive_virus=_get_host_additive_virus(regions=regions,seqname=chrom)
        # if host_additive_virus == "additive": continue
        size = end-start
        if host_additive_virus == "host" or host_additive_virus == "additive":
            if size < host_min: continue
            if size > host_max: continue
        if host_additive_virus == "virus":
            if size < virus_min : continue
            if size > virus_max : continue
        BSJdict[circid]=BSJ(chrom=chrom,start=start,end=end,strand=strand,known_novel=known_novel,read_count=count,counted=counted)
    return(BSJdict)

parser = argparse.ArgumentParser(description='Create CircExplorer2 Per Sample Counts Table')
# INPUTS
parser.add_argument('--back_spliced_bed', dest='bsb', type=str, required=True,
                    help='back_spliced.bed')
parser.add_argument('--back_spliced_min_reads', dest='back_spliced_min_reads', type=int, required=True,
                    help='back_spliced minimum read threshold') # in addition to "known" and "low-conf" circRNAs identified by circexplorer, we also include those found in back_spliced.bed file but not classified as known/low-conf only if the number of reads supporting the BSJ call is greater than this number
parser.add_argument('--circularRNA_known', dest='ck', type=str, required=True,
                    help='circularRNA_known.txt')
parser.add_argument('--low_conf', dest='lc', type=str, required=False,
                    help='low_conf.circularRNA_known.txt')
parser.add_argument('--host', dest='host', type=str, required=True,
                    help='host name eg.hg38... single value...host_filter_min/host_filter_max filters are applied to this region only')
parser.add_argument('--additives', dest='additives', type=str, required=True,
                    help='additive name(s) eg.ERCC... comma-separated list... all BSJs in this region are filtered out')
parser.add_argument('--viruses', dest='viruses', type=str, required=True,
                    help='virus name(s) eg.NC_009333.1... comma-separated list...virus_filter_min/virus_filter_max filters are applied to this region only')
parser.add_argument('--host_filter_min', dest='host_filter_min', type=int, required=False, default=150,
                    help='min BSJ size filter for host')
parser.add_argument('--virus_filter_min', dest='virus_filter_min', type=int, required=False, default=150,
                    help='min BSJ size filter for virus')
parser.add_argument('--host_filter_max', dest='host_filter_max', type=int, required=False, default=5000,
                    help='max BSJ size filter for host')
parser.add_argument('--virus_filter_max', dest='virus_filter_max', type=int, required=False, default=5000,
                    help='max BSJ size filter for virus')
parser.add_argument('--regions', dest='regions', type=str, required=True,
                    help='regions file eg. ref.fa.regions')
# OUTPUTS
parser.add_argument('-o',dest='outfile',required=True,help='counts TSV table')
args = parser.parse_args()

regions=read_regions(regionsfile=args.regions,host=args.host,additives=args.additives,viruses=args.viruses)
o=open(args.outfile,'w')
o.write("#chrom\tstart\tend\tstrand\tread_count\tknown_novel\n")
all_BSJs=read_BSJs(args.bsb,counted=0,threshold=args.back_spliced_min_reads,regions=regions,host_min=args.host_filter_min,host_max=args.host_filter_max,virus_min=args.virus_filter_min,virus_max=args.virus_filter_max)


known_BSJs=read_BSJs(args.ck,known_novel="known",counted=0,threshold=args.back_spliced_min_reads,regions=regions,host_min=args.host_filter_min,host_max=args.host_filter_max,virus_min=args.virus_filter_min,virus_max=args.virus_filter_max)
if args.lc:
    low_conf_BSJs=read_BSJs(args.lc,known_novel="known",counted=0,threshold=args.back_spliced_min_reads,regions=regions,host_min=args.host_filter_min,host_max=args.host_filter_max,virus_min=args.virus_filter_min,virus_max=args.virus_filter_max)
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