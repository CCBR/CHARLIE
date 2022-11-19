import argparse
import pandas

# pandas.options.mode.chained_assignment = None

parser = argparse.ArgumentParser(description='Create per sample Counts Table from NCLscan Outputs')
parser.add_argument('--result', dest='resultsfile', type=str, required=True,
                    help='.result file from NCLscan')
parser.add_argument('-o',dest='outfile',required=True,help='output table')
args = parser.parse_args()

# sn=args.samplename

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

outdf = outdf.astype({"chrom": str, "start": int, "end": int, "strand": str, "read_count": int, "nclscan_annotation": str})

# create index
outdf['circRNA_id']=outdf['chrom'].astype(str)+"##"+outdf['start'].astype(str)+"##"+outdf['end'].astype(str)+"##"+outdf['strand'].astype(str)
outdf.set_index(['circRNA_id'],inplace=True)

# sort and write out
outdf.sort_values(by=['chrom','start'],inplace=True)
outdf.to_csv(args.outfile,sep="\t",header=True,index=False)