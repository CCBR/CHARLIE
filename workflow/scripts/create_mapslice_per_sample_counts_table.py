import argparse
import pandas

# pandas.options.mode.chained_assignment = None

parser = argparse.ArgumentParser(description='Create per sample Counts Table from MapSplice Outputs')
parser.add_argument('--circularRNAstxt', dest='circularRNAstxt', type=str, required=True,
                    help='circular_RNAs.txt file from MapSplice')
parser.add_argument('-o',dest='outfile',required=True,help='output table')
args = parser.parse_args()

# sn=args.samplename

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