#!/usr/bin/env python3
import argparse
import inspect

# DCC counts table input/output file has following columns:
# | #  | colName              | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
# |----|----------------------|-----------------------------------------------------------------|
# | 1  | chr                  |                                                                 |
# | 2  | start                |                                                                 |
# | 3  | end                  |                                                                 |
# | 4  | strand               |                                                                 |
# | 5  | read_count           |                                                                 |
# | 6  | dcc_annotation       | this is JunctionType##Start-End Region from CircCoordinates file|
class DCC:
    def __init__(self,entry,chrom="",start=0,end=0,nreads=0,size=0,host_additive_virus="additive",filter_out=False):
        self.entry=entry
        l=entry.strip().split('\t')
        self.chrom=l[0]
        self.start=int(l[1])
        self.end=int(l[2])
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


parser = argparse.ArgumentParser(description='Filter DCC Per Sample Counts Table')
parser.add_argument('--in_dcc_counts_table', dest='intable', type=str, required=True,
                    help='DCC in file')
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
parser.add_argument('--out_dcc_filtered_counts_table',dest='outfile',required=True,help='filtered DCC out file')
args = parser.parse_args()

regions = read_regions(regionsfile=args.regions,host=args.host,additives=args.additives,viruses=args.viruses)
outfile = open(args.outfile,'w')
infile  = open(args.intable,'r')
alllines = infile.readlines()
header = alllines.pop(0)
outfile.write("%s"%(header))
infile.close()
for l in alllines:
    out = DCC(entry=l)
    out.set_host_additive_virus(regions=regions)
    out.filter_by_nreads(args.back_spliced_min_reads)
    if out.filter_out == False:
        out.filter_by_size(host_min=args.host_filter_min,host_max=args.host_filter_max,virus_min=args.virus_filter_min,virus_max=args.virus_filter_max)
        if out.filter_out == True:
            outfile.write(l)
outfile.close()
