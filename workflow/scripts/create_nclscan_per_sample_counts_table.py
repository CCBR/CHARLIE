import argparse
import pandas

def _annotation_int2str(i):
    if i==0: 
        return "Intergenic"
    elif i==1:
        return "Intragenic"
    else:
        return "Unknown"

# pandas.options.mode.chained_assignment = None

parser = argparse.ArgumentParser(description='Create per sample Counts Table from NCLscan Outputs')
parser.add_argument('--result', dest='resultsfile', type=str, required=True,
                    help='.result file from NCLscan')
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
parser.add_argument('-o',dest='outfile',required=True,help='output table')
parser.add_argument('-fo',dest='filteredoutfile',required=True,help='filtered output table')
args = parser.parse_args()

# sn=args.samplename
# read in regions
hosts = args.host
hosts = hosts.strip().split(",")
viruses = args.viruses
viruses = viruses.strip().split(",")
regions = dict()
regions["host"] = list()
regions["additive"] = list()
regions["virus"] = list()
r = open(args.regions,'r')
rlines = r.readlines()
r.close()
allseqs = list()
for l in rlines:
    l = l.strip().split("\t")
    reg = l[0]
    seq = l[1]
    seq = seq.split()
    host_additive_virus = "additive"
    if reg in hosts:
        host_additive_virus = "host"
    elif reg in viruses:
        host_additive_virus = "virus"
    regions[host_additive_virus].extend(seq)
    allseqs.extend(seq)
regions["additive"] = list((set(allseqs)-set(regions["host"]))-set(regions["virus"]))



# load files
resultsfile=pandas.read_csv(args.resultsfile,sep="\t",header=None)

# file has no column lables ... add them ... the file format is:
# ref: https://github.com/TreesLab/NCLscan
# | #  | Description                                | ColName
# |----|--------------------------------------------|------------ 
# | 1  | Chromosome name of the donor side (5'ss)   | chrd
# | 2  | Junction coordinate of the donor side      | coordd
# | 3  | Strand of the donor side                   | strandd
# | 4  | Chromosome name of the acceptor side (3'ss)| chra
# | 5  | Junction coordinate of the acceptor side   | coorda
# | 6  | Strand of the acceptor side                | stranda
# | 7  | Gene name of the donor side                | gened
# | 8  | Gene name of the acceptor side             | genea
# | 9  | Intragenic (1) or intergenic (0) case      | case --> changed to 2 and 1 as zero means no annotation in merge_per_sample_counts_table.py
# | 10 | Total number of all supporting reads       | reads
# | 11 | Total number of junc-reads                 | jreads
# | 12 | Total number of span-reads                 | sreads

resultsfile.columns=["chrd", "coordd", "strandd", "chra", "coorda", "stranda", "gened", "genea", "case", "reads", "jreads", "sreads"]
resultsfile = resultsfile[resultsfile["chrd"] == resultsfile["chra"]]
resultsfile = resultsfile[resultsfile["strandd"] == resultsfile["stranda"]]

plus_strand = resultsfile[resultsfile['strandd']=='+']
plus_strand = plus_strand[["chrd", "coorda", "coordd", "strandd", "reads", "case"]] # start and end need to be switched!
plus_strand.columns = ['chrom','end','start','strand','read_count', 'nclscan_annotation']

minus_strand = resultsfile[resultsfile['strandd']=='+']
minus_strand = minus_strand[["chrd", "coordd", "coorda", "strandd", "reads", "case"]] 
minus_strand.columns = ['chrom','end','start','strand','read_count', 'nclscan_annotation']

outdf = pandas.concat([plus_strand,minus_strand],ignore_index=True,sort=False)
outdf["nclscan_annotation"] = outdf["nclscan_annotation"] + 1 #change 1 to 2 and 0 to 1 ... as 0 is for no annotation
outdf["nclscan_annotation"] = outdf["nclscan_annotation"].apply(_annotation_int2str)

outdf = outdf.astype({"chrom": str, "start": int, "end": int, "strand": str, "read_count": int, "nclscan_annotation": str})

# create index
outdf['circRNA_id']=outdf['chrom'].astype(str)+"##"+outdf['start'].astype(str)+"##"+outdf['end'].astype(str)+"##"+outdf['strand'].astype(str)
outdf.set_index(['circRNA_id'],inplace=True)

# sort and write out
outdf.sort_values(by=['chrom','start'],inplace=True)
outdf.to_csv(args.outfile,sep="\t",header=True,index=False)

# filter 
# nreads filter
outdf = outdf[~outdf["chrom"].isin(regions["additive"])]
outdf = outdf[outdf["read_count"] >= args.back_spliced_min_reads]

# host distance/size filter
outdf_host = outdf[outdf["chrom"].isin(regions["host"])]
outdf_host["dist"] = abs(outdf_host["start"] - outdf_host["end"])
outdf_host = outdf_host[outdf_host["dist"] > args.host_filter_min]
outdf_host = outdf_host[outdf_host["dist"] < args.host_filter_max]
outdf_host.drop(["dist"],axis=1,inplace=True)

# virus distance/size filter
outdf_virus = outdf[outdf["chrom"].isin(regions["virus"])]
outdf_virus["dist"] = abs(outdf_virus["start"] - outdf_virus["end"])
outdf_virus = outdf_virus[outdf_virus["dist"] > args.virus_filter_min]
outdf_virus = outdf_virus[outdf_virus["dist"] < args.virus_filter_max]
outdf_virus.drop(["dist"],axis=1,inplace=True)

outdf = pandas.concat([outdf_host,outdf_virus])
# sort and write out
outdf.sort_values(by=['chrom','start'],inplace=True)
outdf.to_csv(args.filteredoutfile,sep="\t",header=True,index=False)
