import argparse
import pandas

# pandas.options.mode.chained_assignment = None

parser = argparse.ArgumentParser(description='Create per sample Counts Table from MapSplice Outputs')
parser.add_argument('--circularRNAstxt', dest='circularRNAstxt', type=str, required=True,
                    help='circular_RNAs.txt file from MapSplice')
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
circularRNAstxt=pandas.read_csv(args.circularRNAstxt,sep="\t",header=None)

# file has no column lables ... add them
# ref: https://github.com/Aufiero/circRNAprofiler/blob/master/R/importFilesPredictionTool.R
circularRNAstxt.columns=["chrom", "donor_end", "acceptor_start", "id", "coverage", "strand", "rgb", "block_count", "block_size", "block_distance", "entropy", "flank_case", "flank_string", "min_mismatch", "max_mismatch", "ave_mismatch", "max_min_suffix", "max_min_prefix", "min_anchor_difference", "unique_read_count", "multi_read_count", "paired_read_count", "left_paired_read_count", "right_paired_read_count", "multiple_paired_read_count", "unique_paired_read_count", "single_read_count", "encompassing_read", "doner_start", "acceptor_end", "doner_iosforms", "acceptor_isoforms", "obsolete1", "obsolete2", "obsolete3", "obsolete4", "minimal_doner_isoform_length", "maximal_doner_isoform_length", "minimal_acceptor_isoform_length", "maximal_acceptor_isoform_length", "paired_reads_entropy", "mismatch_per_bp", "anchor_score", "max_doner_fragment", "max_acceptor_fragment", "max_cur_fragment", "min_cur_fragment", "ave_cur_fragment", "doner_encompass_unique", "doner_encompass_multiple", "acceptor_encompass_unique", "acceptor_encompass_multiple", "doner_match_to_normal", "acceptor_match_to_normal", "doner_seq", "acceptor_seq", "match_gene_strand", "annotated_type", "fusion_type", "gene_strand", "annotated_gene_donor", "annotated_gene_acceptor", "dummy"]

# 'chrom' is in the format 'donor_chr~acceptor_chr' ... hence needs to be split
circularRNAstxt[['Donor', 'Acceptor']] = circularRNAstxt['chrom'].str.split('~', expand=True)

# only select rows with ++ or -- strand
circularRNAstxtnew = pandas.concat([circularRNAstxt[circularRNAstxt['strand'] == '++' ],circularRNAstxt[circularRNAstxt['strand'] == '--' ]],ignore_index=True,sort=False)
# strand is either ++ or -- .. needs to be fixed to + or -
circularRNAstxtnew.replace('++','+',inplace=True)
circularRNAstxtnew.replace('--','-',inplace=True)

# subset columns and rename them and fix start/end order
circularRNAstxtnew=circularRNAstxtnew[['Acceptor', 'donor_end', 'acceptor_start', 'strand', 'coverage', 'fusion_type', 'entropy']]

plus_strand = circularRNAstxtnew[circularRNAstxtnew['strand']=='+']
plus_strand.columns = ['chrom','end','start','strand','read_count','fusion_type', 'entropy'] # start and end need to be switched!

minus_strand = circularRNAstxtnew[circularRNAstxtnew['strand']=='-']
minus_strand.columns = ['chrom','start','end','strand','read_count','fusion_type', 'entropy']

circularRNAstxtnew = pandas.concat([plus_strand,minus_strand],ignore_index=True,sort=False)
# circularRNAstxtnew.columns=['chrom','start','end','strand','read_count','fusion_type', 'entropy']

# create mapsplice_annotation column to include "fusion_type" along with "entropy"
circularRNAstxtnew['mapsplice_annotation']=circularRNAstxtnew['fusion_type'].astype(str)+"##"+circularRNAstxtnew['entropy'].astype(str)
circularRNAstxtnew.drop(['fusion_type', 'entropy'],axis=1,inplace=True)
circularRNAstxtnew.fillna(value="-11",inplace=True)
circularRNAstxtnew = circularRNAstxtnew.astype({"chrom": str, "start": int, "end": int, "strand": str, "read_count": int, "mapsplice_annotation": str})

# create index
circularRNAstxtnew['circRNA_id']=circularRNAstxtnew['chrom'].astype(str)+"##"+circularRNAstxtnew['start'].astype(str)+"##"+circularRNAstxtnew['end'].astype(str)+"##"+circularRNAstxtnew['strand'].astype(str)
circularRNAstxtnew.set_index(['circRNA_id'],inplace=True)

# sort and write out
circularRNAstxtnew.sort_values(by=['chrom','start'],inplace=True)
circularRNAstxtnew.to_csv(args.outfile,sep="\t",header=True,index=False)

# filter 
# nreads filter
circularRNAstxtnew = circularRNAstxtnew[~circularRNAstxtnew["chrom"].isin(regions["additive"])]
circularRNAstxtnew = circularRNAstxtnew[circularRNAstxtnew["read_count"] >= args.back_spliced_min_reads]

# host distance/size filter
circularRNAstxtnew_host = circularRNAstxtnew[circularRNAstxtnew["chrom"].isin(regions["host"])]
circularRNAstxtnew_host["dist"] = abs(circularRNAstxtnew_host["start"] - circularRNAstxtnew_host["end"])
circularRNAstxtnew_host = circularRNAstxtnew_host[circularRNAstxtnew_host["dist"] > args.host_filter_min]
circularRNAstxtnew_host = circularRNAstxtnew_host[circularRNAstxtnew_host["dist"] < args.host_filter_max]
circularRNAstxtnew_host.drop(["dist"],axis=1,inplace=True)

# virus distance/size filter
circularRNAstxtnew_virus = circularRNAstxtnew[circularRNAstxtnew["chrom"].isin(regions["virus"])]
circularRNAstxtnew_virus["dist"] = abs(circularRNAstxtnew_virus["start"] - circularRNAstxtnew_virus["end"])
circularRNAstxtnew_virus = circularRNAstxtnew_virus[circularRNAstxtnew_virus["dist"] > args.virus_filter_min]
circularRNAstxtnew_virus = circularRNAstxtnew_virus[circularRNAstxtnew_virus["dist"] < args.virus_filter_max]
circularRNAstxtnew_virus.drop(["dist"],axis=1,inplace=True)

circularRNAstxtnew = pandas.concat([circularRNAstxtnew_host,circularRNAstxtnew_virus])
# sort and write out
circularRNAstxtnew.sort_values(by=['chrom','start'],inplace=True)
circularRNAstxtnew.to_csv(args.filteredoutfile,sep="\t",header=True,index=False)