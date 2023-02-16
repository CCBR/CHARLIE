#!/usr/bin/env python3
import argparse
import inspect

# CIRI2 output file has following columns:
# | #  | colName              | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
# |----|----------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
# | 1  | circRNA_ID           | ID of a predicted circRNA in the pattern of "chr:start|end"                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
# | 2  | chr                  | chromosome of a predicted circRNA                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
# | 3  | circRNA_start        | start loci of a predicted circRNA on the chromosome                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
# | 4  | circRNA_end          | end loci of a predicted circRNA on the chromosome                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
# | 5  | #junction_reads      | circular junction read (also called as back-spliced junction read) count of a predicted circRNA                                                                                                                                                                                                                                                                                                                                                                                                         |
# | 6  | SM_MS_SMS            | unique CIGAR types of a predicted circRNA. For example, a circRNAs have three junction reads: read A (80M20S, 80S20M), read B (80M20S, 80S20M), read C (40M60S, 40S30M30S, 70S30M), then its has two SM types (80S20M, 70S30M), two MS types (80M20S, 70M30S) and one SMS type (40S30M30S). Thus its SM_MS_SMS should be 2_2_1.                                                                                                                                                                         |
# | 7  | #non_junction_reads  | non-junction read count of a predicted circRNA that mapped across the circular junction but consistent with linear RNA instead of being back-spliced                                                                                                                                                                                                                                                                                                                                                    |
# | 8  | junction_reads_ratio | ratio of circular junction reads calculated by 2#junction_reads/(2#junction_reads+#non_junction_reads). #junction_reads is multiplied by two because a junction read is generated from two ends of circular junction but only counted once while a non-junction read is from one end. It has to be mentioned that the non-junction reads are still possibly from another larger circRNA, so the junction_reads_ratio based on it may be an inaccurate estimation of relative expression of the circRNA. |
# | 9  | circRNA_type         | type of a circRNA according to positions of its two ends on chromosome (exon, intron or intergenic_region; only available when annotation file is provided)                                                                                                                                                                                                                                                                                                                                             |
# | 10 | gene_id              | ID of the gene(s) where an exonic or intronic circRNA locates                                                                                                                                                                                                                                                                                                                                                                                                                                           |
# | 11 | strand               | strand info of a predicted circRNAs (new in CIRI2)                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
# | 12 | junction_reads_ID    | all of the circular junction read IDs (split by ",")                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
# ref: https://ciri-cookbook.readthedocs.io/en/latest/CIRI2.html#an-example-of-running-ciri2
class CIRIOUT:
    def __init__(self,entry,chrom="",start=0,end=0,nreads=0,size=0,host_additive_virus="additive",filter_out=False):
        self.entry=entry
        l=entry.strip().split('\t')
        self.chrom=l[1]
        self.start=int(l[2])-1
        self.end=int(l[3])
        self.nreads=int(l[4])
        self.size=self.start-self.end
        self.filter_out=False
    
    # @classmethod
    def set_host_additive_virus(self,regions):
        self.host_additive_virus=_get_host_additive_virus(regions=regions,seqname=self.chrom)
    
    # @classmethod
    def filter_by_nreads(self,minreads):
        if self.nreads < minreads: self.filter_out=True
    
    # @classmethod
    def filter_by_size(self,host_min,host_max,virus_min,virus_max):
        if self.host_additive_virus=="host":
            if self.size < host_min : self.filter_out=True
            if self.size > host_max : self.filter_out=True
        elif self.host_additive_virus=="virus":
            if self.size < virus_min : self.filter_out=True
            if self.size > virus_max : self.filter_out=True
        else:
            self.filter_out=True

    
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


parser = argparse.ArgumentParser(description='Filter CIRI2 Per Sample Counts Table')
parser.add_argument('--ciriout', dest='ciriout', type=str, required=True,
                    help='ciri out file')
parser.add_argument('--back_spliced_min_reads', dest='back_spliced_min_reads', type=int, required=True,
                    help='back_spliced minimum read threshold') 
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
parser.add_argument('-o',dest='outfile',required=True,help='filtered ciriout file')
args = parser.parse_args()

regions = read_regions(regionsfile=args.regions,host=args.host,additives=args.additives,viruses=args.viruses)
outfile = open(args.outfile,'w')
infile  = open(args.ciriout,'r')
alllines = infile.readlines()
header = alllines.pop(0)
outfile.write("%s"%(header))
infile.close()
for l in alllines:
    out = CIRIOUT(entry=l)
    out.set_host_additive_virus(regions=regions)
    out.filter_by_nreads(args.back_spliced_min_reads)
    if out.filter_out == False:
        out.filter_by_size(host_min=args.host_filter_min,host_max=args.host_filter_max,virus_min=args.virus_filter_min,virus_max=args.virus_filter_max)
        if out.filter_out == True:
            outfile.write(l)
outfile.close()
