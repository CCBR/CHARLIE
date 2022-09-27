import argparse
import pandas

parser = argparse.ArgumentParser(description='Merge per sample Counts from different circRNA detection tools')
parser.add_argument('--circExplorer', dest='circE', type=str, required=True,
                    help='circExplorer2 table')
parser.add_argument('--ciri', dest='ciri', type=str, required=True,
                    help='ciri2 output')
parser.add_argument('--samplename', dest='samplename', type=str, required=True,
                    help='Sample Name')
parser.add_argument('-o',dest='outfile',required=True,help='merged table')
args = parser.parse_args()

sn=args.samplename

# load circExplorer
circE=pandas.read_csv(args.circE,sep="\t",header=0)
# columns are:
# chrom	start	end	strand	read_count	known_novel
# circE['circRNA_id']="##".join([circE['chrom'],str(circE['start']),str(circE['end']),circE['strand']])
circE['circRNA_id']=circE['chrom'].astype(str)+"##"+circE['start'].astype(str)+"##"+circE['end'].astype(str)+"##"+circE['strand'].astype(str)
circE.rename({'read_count': sn+'_circExplorer_read_count', 'known_novel': sn+'_circExplorer_known_novel', 'strand': sn+'_circExplorer_strand'}, axis=1, inplace=True)
circE.set_index(['circRNA_id'],inplace=True)

# load ciri
ciri=pandas.read_csv(args.ciri,sep="\t",header=0,usecols=['chr', 'circRNA_start', 'circRNA_end', '#junction_reads', 'circRNA_type', 'strand'])
# columns are:
# circRNA_ID	chr	circRNA_start	circRNA_end	#junction_reads	SM_MS_SMS	#non_junction_reads	junction_reads_ratio	circRNA_type	gene_id	strand	junction_reads_ID
# ciri['circRNA_id']="##".join([ciri['chr'],str(ciri['circRNA_start']),str(ciri['circRNA_end']),ciri['strand']])
ciri["circRNA_start"]=ciri["circRNA_start"].astype(int)-1
ciri['circRNA_id']=ciri['chr'].astype(str)+"##"+ciri['circRNA_start'].astype(str)+"##"+ciri['circRNA_end'].astype(str)+"##"+ciri['strand'].astype(str)
ciri.rename({'#junction_reads': sn+'_ciri_read_count', 'circRNA_type': sn+'_circRNA_type','strand': sn+'_ciri_strand'}, axis=1, inplace=True)
ciri.set_index(['circRNA_id'],inplace=True)

merged_counts=pandas.concat([circE,ciri],axis=1,join="outer",sort=False)
merged_counts['circRNA_id']=merged_counts.index
merged_counts.drop(['chrom', 'start', 'end', sn+'_circExplorer_strand', 'chr', 'circRNA_start', 'circRNA_end', sn+'_ciri_strand'], axis = 1,inplace=True)
merged_counts.fillna(0,inplace=True)
merged_counts = merged_counts.astype({'circRNA_id': 'str', sn+'_circExplorer_read_count': int, sn+'_ciri_read_count': int, sn+'_circExplorer_known_novel' : str,sn+'_circRNA_type' : str})
merged_counts[ sn+'_ntools'] = 0
merged_counts.loc[merged_counts[sn+'_circExplorer_read_count'] > 0, sn+'_ntools'] += 1
merged_counts.loc[merged_counts[sn+'_ciri_read_count'] > 0, sn+'_ntools'] += 1
merged_counts[['chrom', 'start', 'end', 'strand']] = merged_counts['circRNA_id'].str.split('##', expand=True)
merged_counts.drop(['circRNA_id'],axis=1,inplace=True)
merged_counts['circRNA_id']=merged_counts['chrom'].astype(str)+":"+merged_counts['start'].astype(str)+"-"+merged_counts['end'].astype(str)
merged_counts = merged_counts[['circRNA_id', 'strand', sn+'_circExplorer_read_count', sn+'_ciri_read_count', sn+'_circExplorer_known_novel', sn+'_circRNA_type', sn+'_ntools']]
merged_counts.to_csv(args.outfile,sep="\t",header=True,index=False)