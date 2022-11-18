import argparse
import pandas

parser = argparse.ArgumentParser(description='Merge per sample Counts from different circRNA detection tools')
parser.add_argument('--circExplorer', dest='circE', type=str, required=True,
                    help='circExplorer2 per-sample counts table')
parser.add_argument('--ciri', dest='ciri', type=str, required=True,
                    help='ciri2 per-sample output')
parser.add_argument('--dcc', dest='dcc', type=str, required=False,
                    help='dcc per-sample counts table')
parser.add_argument('--mapsplice', dest='mapsplice', type=str, required=False,
                    help='mapsplice per-sample counts table')
parser.add_argument('--nclscan', dest='nclscan', type=str, required=False,
                    help='nclscan per-sample counts table')
parser.add_argument('--samplename', dest='samplename', type=str, required=True,
                    help='Sample Name')
parser.add_argument('--min_read_count_reqd', dest='minreads', type=int, required=False, default=2,
                    help='Read count threshold..circRNA with lower than this number of read support are excluded! (default=2)')
parser.add_argument('-o',dest='outfile',required=True,help='merged table')
args = parser.parse_args()

sn=args.samplename

dfs=[]

# load circExplorer
circE=pandas.read_csv(args.circE,sep="\t",header=0)
# columns are:
# chrom	start	end	strand	read_count	known_novel
# | # | ColName     |
# |---|-------------|
# | 1 | chrom       |
# | 2 | start       |
# | 3 | end         |
# | 4 | strand      |
# | 5 | read_count  |
# | 6 | known_novel |
circE['circRNA_id']=circE['chrom'].astype(str)+"##"+circE['start'].astype(str)+"##"+circE['end'].astype(str)+"##"+circE['strand'].astype(str)
circE.rename({'read_count': sn+'_circExplorer_read_count', 'known_novel':'circExplorer_annotation'}, axis=1, inplace=True)
circE.drop(['chrom','start', 'end','strand'], axis = 1,inplace=True)
circE.set_index(['circRNA_id'],inplace=True)
dfs.append(circE)

# load ciri
ciri=pandas.read_csv(args.ciri,sep="\t",header=0,usecols=['chr', 'circRNA_start', 'circRNA_end', '#junction_reads', 'circRNA_type', 'strand'])
# columns are:
# circRNA_ID	chr	circRNA_start	circRNA_end	#junction_reads	SM_MS_SMS	#non_junction_reads	junction_reads_ratio	circRNA_type	gene_id	strand	junction_reads_ID
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
# | 9  | x         | type of a circRNA according to positions of its two ends on chromosome (exon, intron or intergenic_region; only available when annotation file is provided)                                                                                                                                                                                                                                                                                                                                             |
# | 10 | gene_id              | ID of the gene(s) where an exonic or intronic circRNA locates                                                                                                                                                                                                                                                                                                                                                                                                                                           |
# | 11 | strand               | strand info of a predicted circRNAs (new in CIRI2)                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
# | 12 | junction_reads_ID    | all of the circular junction read IDs (split by ",")   
ciri["circRNA_start"]=ciri["circRNA_start"].astype(int)-1
ciri['circRNA_id']=ciri['chr'].astype(str)+"##"+ciri['circRNA_start'].astype(str)+"##"+ciri['circRNA_end'].astype(str)+"##"+ciri['strand'].astype(str)
ciri.rename({'#junction_reads': sn+'_ciri_read_count', 'circRNA_type': 'ciri_annotation'}, axis=1, inplace=True)
ciri.drop(['chr','circRNA_start', 'circRNA_end','strand'], axis = 1,inplace=True)
ciri.set_index(['circRNA_id'],inplace=True)
dfs.append(ciri)

# load dcc
if [args.dcc]:
    dcc=pandas.read_csv(args.dcc,sep="\t",header=0)
    # output dcc.counts_table.tsv has the following columns:
    # | # | ColName        |
    # |---|----------------|
    # | 1 | chr            |
    # | 2 | start          |
    # | 3 | end            |
    # | 4 | strand         |
    # | 5 | read_count     |
    # | 6 | dcc_annotation | --> this is JunctionType##Start-End Region from CircCoordinates file
    dcc["start"]=dcc["start"].astype(int)-1
    dcc['circRNA_id']=dcc['chr'].astype(str)+"##"+dcc['start'].astype(str)+"##"+dcc['end'].astype(str)+"##"+dcc['strand'].astype(str)
    dcc.rename({'read_count': sn+'_dcc_read_count'}, axis=1, inplace=True)
    dcc.drop(['chr','start', 'end','strand'], axis = 1,inplace=True)
    dcc.set_index(['circRNA_id'],inplace=True)
    dfs.append(dcc)


# load mapsplice
if [args.mapsplice]:
    mapsplice=pandas.read_csv(args.mapsplice,sep="\t",header=0)
    # output .mapslice.counts_table.tsv has the following columns:
    # | # | ColName              | Eg.              |
    # |---|----------------------|------------------|
    # | 1 | chrom                | chr1             |
    # | 2 | start                | 1223244          |
    # | 3 | end                  | 1223968          |
    # | 4 | strand               | -                |
    # | 5 | read_count           | 26               |
    # | 6 | mapsplice_annotation | normal##2.811419 | <--"fusion_type"##"entropy" 
    # "fusion_type" is either "normal" or "overlapping" ... higher "entropy" values are better!
    mapsplice["start"]=mapsplice["start"].astype(int)-1
    mapsplice['circRNA_id']=mapsplice['chrom'].astype(str)+"##"+mapsplice['start'].astype(str)+"##"+mapsplice['end'].astype(str)+"##"+mapsplice['strand'].astype(str)
    mapsplice.rename({'read_count': sn+'_mapsplice_read_count'}, axis=1, inplace=True)
    mapsplice.drop(['chrom','start', 'end','strand'], axis = 1,inplace=True)
    mapsplice.set_index(['circRNA_id'],inplace=True)
    dfs.append(mapsplice)

# load nclscan
if [args.nclscan]:
    nclscan=pandas.read_csv(args.nclscan,sep="\t",header=0)
    # output .mapslice.counts_table.tsv has the following columns:
    # | # | ColName              | Eg.              |
    # |---|----------------------|------------------|
    # | 1 | chrom                | chr1             |
    # | 2 | start                | 1223244          |
    # | 3 | end                  | 1223968          |
    # | 4 | strand               | -                |
    # | 5 | read_count           | 26               |
    # | 6 | nclscan_annotation   | 1                | <--1 for intragenic 0 for intergenic
    nclscan["start"]=nclscan["start"].astype(int)-1
    nclscan['circRNA_id']=nclscan['chrom'].astype(str)+"##"+nclscan['start'].astype(str)+"##"+nclscan['end'].astype(str)+"##"+nclscan['strand'].astype(str)
    nclscan.rename({'read_count': sn+'_nclscan_read_count'}, axis=1, inplace=True)
    nclscan.drop(['chrom','start', 'end','strand'], axis = 1,inplace=True)
    nclscan.set_index(['circRNA_id'],inplace=True)
    dfs.append(nclscan)



merged_counts=pandas.concat(dfs,axis=1,join="outer",sort=False)
merged_counts['circRNA_id']=merged_counts.index
merged_counts.fillna(0,inplace=True)
merged_counts[ sn+'_ntools'] = 0
if [args.dcc] and not [args.mapsplice] and not [args.nclscan]:
    merged_counts = merged_counts.astype({'circRNA_id': 'str', sn+'_ntools' : int, sn+'_circExplorer_read_count': int, sn+'_ciri_read_count': int, sn+'_dcc_read_count': int,'circExplorer_annotation' : str,'ciri_annotation' : str, 'dcc_annotation' : str})
elif [args.dcc] and [args.mapsplice] and not [args.nclscan]:
    merged_counts = merged_counts.astype({'circRNA_id': 'str', sn+'_ntools' : int, sn+'_circExplorer_read_count': int, sn+'_ciri_read_count': int, sn+'_dcc_read_count': int, sn+'_mapsplice_read_count': int, 'circExplorer_annotation' : str,'ciri_annotation' : str, 'dcc_annotation' : str, 'mapsplice_annotation' : str})
elif [args.dcc] and [args.mapsplice] and [args.nclscan]:
    merged_counts = merged_counts.astype({'circRNA_id': 'str', sn+'_ntools' : int, sn+'_circExplorer_read_count': int, sn+'_ciri_read_count': int, sn+'_dcc_read_count': int, sn+'_mapsplice_read_count': int, sn+'_nclscan_read_count': int, 'circExplorer_annotation' : str,'ciri_annotation' : str, 'dcc_annotation' : str, 'mapsplice_annotation' : str, 'nclscan_annotation' : int})
elif [args.dcc] and not [args.mapsplice] and [args.nclscan]:
    merged_counts = merged_counts.astype({'circRNA_id': 'str', sn+'_ntools' : int, sn+'_circExplorer_read_count': int, sn+'_ciri_read_count': int, sn+'_dcc_read_count': int, sn+'_nclscan_read_count': int, 'circExplorer_annotation' : str,'ciri_annotation' : str, 'dcc_annotation' : str, 'nclscan_annotation' : int})
elif not [args.dcc] and [args.mapsplice] and [args.nclscan]:
    merged_counts = merged_counts.astype({'circRNA_id': 'str', sn+'_ntools' : int, sn+'_circExplorer_read_count': int, sn+'_ciri_read_count': int, sn+'_mapsplice_read_count': int, sn+'_nclscan_read_count': int, 'circExplorer_annotation' : str,'ciri_annotation' : str, 'mapsplice_annotation' : str, 'nclscan_annotation' : int})
else:
    merged_counts = merged_counts.astype({'circRNA_id': 'str', sn+'_ntools' : int, sn+'_circExplorer_read_count': int, sn+'_ciri_read_count': int, 'circExplorer_annotation' : str,'ciri_annotation' : str})

merged_counts.loc[merged_counts[sn+'_circExplorer_read_count'] > args.minreads, sn+'_ntools'] += 1
merged_counts.loc[merged_counts[sn+'_ciri_read_count'] > args.minreads, sn+'_ntools'] += 1
if [args.dcc]: merged_counts.loc[merged_counts[sn+'_dcc_read_count'] > args.minreads, sn+'_ntools'] += 1
if [args.mapsplice]: merged_counts.loc[merged_counts[sn+'_mapsplice_read_count'] > args.minreads, sn+'_ntools'] += 1
if [args.nclscan]: merged_counts.loc[merged_counts[sn+'_nclscan_read_count'] > args.minreads, sn+'_ntools'] += 1
merged_counts[['chrom', 'start', 'end', 'strand']] = merged_counts['circRNA_id'].str.split('##', expand=True)
merged_counts.drop(['circRNA_id'],axis=1,inplace=True)
merged_counts['circRNA_id']=merged_counts['chrom'].astype(str)+":"+merged_counts['start'].astype(str)+"-"+merged_counts['end'].astype(str)
outcols=['circRNA_id', 'strand', sn+'_ntools', sn+'_circExplorer_read_count', sn+'_ciri_read_count']
if [args.dcc]: outcols.append(sn+'_dcc_read_count')
if [args.mapsplice]: outcols.append(sn+'_mapsplice_read_count')
if [args.nclscan]: outcols.append(sn+'_nclscan_read_count')
outcols.extend(['circExplorer_annotation', 'ciri_annotation'])
if [args.dcc]: outcols.append('dcc_annotation')
if [args.mapsplice]: outcols.append('mapsplice_annotation')
if [args.nclscan]: outcols.append('nclscan_annotation')
merged_counts = merged_counts[outcols]
merged_counts.to_csv(args.outfile,sep="\t",header=True,index=False)