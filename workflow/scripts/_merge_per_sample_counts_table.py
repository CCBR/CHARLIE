import argparse
import pandas
import HTSeq
import sys
import gzip

def _df_setcol_as_int(df,collist):
    for c in collist:
        df[[c]]=df[[c]].astype(int)
    return df

def _df_setcol_as_str(df,collist):
    for c in collist:
        df[[c]]=df[[c]].astype(str)
    return df

def _df_setcol_as_float(df,collist):
    for c in collist:
        df[[c]]=df[[c]].astype(float)
    return df

class BSJ:
    def __init__(self,chrom,start,end,strand):
        self.chrom=chrom
        self.start=int(start)
        self.end=int(end)
        self.strand=strand
        self.splice_site_flank_5="" #donor
        self.splice_site_flank_3="" #acceptor

    def add_flanks(self,sequences):
        if self.strand == '+':
            coord = int(self.end)
            self.splice_site_flank_5 = sequences[self.chrom][coord:coord+2]
            coord = int(self.start)
            self.splice_site_flank_3 = sequences[self.chrom][coord-2:coord]
        elif self.strand == '-':
            coord = int(self.end)
            seq =  sequences[self.chrom][coord:coord+2]
            seq = seq.replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c")
            seq = seq.upper()
            self.splice_site_flank_3 = seq
            coord = int(self.start)
            seq = sequences[self.chrom][coord-2:coord]
            seq = seq.replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c")
            seq = seq.upper()
            self.splice_site_flank_5 = seq
    
    def get_flanks(self):
        return self.splice_site_flank_5+"##"+self.splice_site_flank_3

def main() :
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
    parser.add_argument('--circrnafinder', dest='circrnafinder', type=str, required=False,
        help='circrnafinder per-sample counts table')
    parser.add_argument('--samplename', dest='samplename', type=str, required=True,
        help='Sample Name')
    parser.add_argument('--min_read_count_reqd', dest='minreads', type=int, required=False, default=2,
        help='Read count threshold..circRNA with lower than this number of read support are excluded! (default=2)')
    parser.add_argument("--reffa",dest="reffa",required=True,type=argparse.FileType('r'),default=sys.stdin,
        help="reference fasta file")
    parser.add_argument('-o',dest='outfile',required=True,help='merged table')
    args = parser.parse_args()

    sn=args.samplename

    dfs=[]

    # load circExplorer
    circE=pandas.read_csv(args.circE,sep="\t",header=0)
    # columns are:
    # | #  | Column                                   |
    # | -- | ---------------------------------------- |
    # | 1  | #chrom                                   |
    # | 2  | start                                    |
    # | 3  | end                                      |
    # | 4  | strand                                   |
    # | 5  | known_novel                              |
    # | 6  | expected_BSJ_reads                       |
    # | 7  | found_BSJ_reads                          |
    # | 8  | linear_same_strand                       |
    # | 9  | spliced_same_strand                      |
    # | 10 | linear_opposite_strand                   |
    # | 11 | spliced_opposite_strand                  |
    # | 12 | linear_unknown_strand                    |
    # | 13 | spliced_unknown_strand                   |
    circE['circRNA_id']=circE['#chrom'].astype(str)+"##"+circE['start'].astype(str)+"##"+circE['end'].astype(str)+"##"+circE['strand'].astype(str)
    circE.rename({'known_novel' : 'circExplorer_annotation',
                'expected_BSJ_reads' : 'circExplorer_read_count',
                'found_BSJ_reads' : 'circExplorer_found_BSJcounts',
                'linear_same_strand' : 'circExplorer_found_linear_BSJ_same_strand_counts',
                'spliced_same_strand' : 'circExplorer_found_linear_spliced_BSJ_same_strand_counts',
                'linear_opposite_strand' : 'circExplorer_found_linear_BSJ_opposite_strand_counts',
                'spliced_opposite_strand' : 'circExplorer_found_linear_spliced_BSJ_opposite_strand_counts',
                'linear_unknown_strand' : 'circExplorer_found_linear_BSJ_unknown_strand_counts',
                'spliced_unknown_strand' : 'circExplorer_found_linear_spliced_BSJ_unknown_strand_counts'}, axis=1, inplace=True)
    circE.drop(['#chrom','start', 'end','strand'], axis = 1,inplace=True)
    circE.set_index(['circRNA_id'],inplace=True)
    
    circE.fillna(value=-1,inplace=True)

    intcols = [ 'circExplorer_read_count', 
                'circExplorer_found_BSJcounts', 
                'circExplorer_found_linear_BSJ_same_strand_counts', 
                'circExplorer_found_linear_spliced_BSJ_same_strand_counts', 
                'circExplorer_found_linear_BSJ_opposite_strand_counts', 
                'circExplorer_found_linear_spliced_BSJ_opposite_strand_counts', 
                'circExplorer_found_linear_BSJ_unknown_strand_counts', 
                'circExplorer_found_linear_spliced_BSJ_unknown_strand_counts' ]
    strcols = list ( set(circE.columns) - set(intcols) )
    circE = _df_setcol_as_int(circE,intcols)
    circE = _df_setcol_as_str(circE,strcols)

    dfs.append(circE)

    # load ciri
    ciri=pandas.read_csv(args.ciri,sep="\t",header=0,usecols=['chr', 'circRNA_start', 'circRNA_end', '#junction_reads', '#non_junction_reads', 'circRNA_type', 'strand'])
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
    # | 9  | circRNA_type         | type of a circRNA according to positions of its two ends on chromosome (exon, intron or intergenic_region; only available when annotation file is provided)                                                                                                                                                                                                                                                                                                                                             |
    # | 10 | gene_id              | ID of the gene(s) where an exonic or intronic circRNA locates                                                                                                                                                                                                                                                                                                                                                                                                                                           |
    # | 11 | strand               | strand info of a predicted circRNAs (new in CIRI2)                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
    # | 12 | junction_reads_ID    | all of the circular junction read IDs (split by ",")   
    ciri["circRNA_start"]=ciri["circRNA_start"].astype(int)-1
    ciri['circRNA_id']=ciri['chr'].astype(str)+"##"+ciri['circRNA_start'].astype(str)+"##"+ciri['circRNA_end'].astype(str)+"##"+ciri['strand'].astype(str)
    ciri.rename({   '#junction_reads': 'ciri_read_count', 
                    '#non_junction_reads' : 'ciri_linear_read_count', 
                    'circRNA_type': 'ciri_annotation'}, axis=1, inplace=True)
    ciri.drop(['chr','circRNA_start', 'circRNA_end','strand'], axis = 1,inplace=True)
    ciri.set_index(['circRNA_id'],inplace=True)

    ciri.fillna(value=-1,inplace=True)

    intcols = [ 'ciri_read_count',
                'ciri_linear_read_count' ]
    strcols = list ( set(ciri.columns) - set(intcols) )
    ciri = _df_setcol_as_int(ciri,intcols)
    if len(strcols) > 0: ciri = _df_setcol_as_str(ciri,strcols)
    

    dfs.append(ciri)

    # load dcc
    if args.dcc:
        dcc=pandas.read_csv(args.dcc,sep="\t",header=0)
        # output dcc.counts_table.tsv has the following columns:
        # | # | ColName        |
        # |---|----------------|
        # | 1 | chr            |
        # | 2 | start          |
        # | 3 | end            |
        # | 4 | strand         |
        # | 5 | read_count     |
        # | 6 | linear_read_count|
        # | 7 | dcc_annotation | --> this is gene##JunctionType##Start-End Region from CircCoordinates file
        dcc["start"]=dcc["start"].astype(int)-1
        dcc['circRNA_id']=dcc['chr'].astype(str)+"##"+dcc['start'].astype(str)+"##"+dcc['end'].astype(str)+"##"+dcc['strand'].astype(str)
        dcc.rename({'read_count': 'dcc_read_count'}, axis=1, inplace=True)
        dcc.rename({'linear_read_count': 'dcc_linear_read_count'}, axis=1, inplace=True)
        dcc[['dcc_gene', 'dcc_junction_type', 'dcc_annotation2']] = dcc['dcc_annotation'].apply(lambda x: pandas.Series(str(x).split("##")))
        dcc.drop(['chr','start', 'end','strand','dcc_annotation'], axis = 1,inplace=True)
        dcc.rename({'dcc_annotation2': 'dcc_annotation'}, axis=1, inplace=True)
        dcc.set_index(['circRNA_id'],inplace=True)

        dcc.fillna(value=-1,inplace=True)

        intcols = [ 'dcc_read_count',
                    'dcc_linear_read_count' ]
        strcols = list ( set(dcc.columns) - set(intcols) )
        dcc = _df_setcol_as_int(dcc,intcols)
        if len(strcols) > 0: dcc = _df_setcol_as_str(dcc,strcols)    

        dfs.append(dcc)


    # load mapsplice
    if args.mapsplice:
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
        mapsplice.rename({'read_count': 'mapsplice_read_count'}, axis=1, inplace=True)
        mapsplice[['mapsplice_annotation2', 'mapsplice_entropy']] = mapsplice['mapsplice_annotation'].apply(lambda x: pandas.Series(str(x).split("##")))
        mapsplice.drop(['chrom','start', 'end','strand','mapsplice_annotation'], axis = 1,inplace=True)
        mapsplice.rename({'mapsplice_annotation2': 'mapsplice_annotation'}, axis=1, inplace=True)
        mapsplice.set_index(['circRNA_id'],inplace=True)

        mapsplice.fillna(value=-1,inplace=True)

        intcols = [ 'mapsplice_read_count' ]
        mapsplice = _df_setcol_as_int(mapsplice,intcols)
        floatcols = [ 'mapsplice_entropy' ]
        mapsplice = _df_setcol_as_float(mapsplice,floatcols)
        strcols = list ( ( set(mapsplice.columns) - set(intcols) ) - set(floatcols) )
        if len(strcols) > 0: mapsplice = _df_setcol_as_str(mapsplice,strcols) 

        dfs.append(mapsplice)

    # load nclscan
    if args.nclscan:
        nclscan=pandas.read_csv(args.nclscan,sep="\t",header=0)
        # output nslscan table has the following columns:
        # | # | ColName              | Eg.              |
        # |---|----------------------|------------------|
        # | 1 | chrom                | chr1             |
        # | 2 | start                | 1223244          |
        # | 3 | end                  | 1223968          |
        # | 4 | strand               | -                |
        # | 5 | read_count           | 26               |
        # | 6 | nclscan_annotation   | 1                | <--1 for intragenic 0 for intergenic
        includenclscan=True
        if nclscan.shape[0]==0: includenclscan=False
        if includenclscan:
            nclscan["start"]=nclscan["start"].astype(int)-1
            nclscan['circRNA_id']=nclscan['chrom'].astype(str)+"##"+nclscan['start'].astype(str)+"##"+nclscan['end'].astype(str)+"##"+nclscan['strand'].astype(str)
            nclscan.rename({'read_count': 'nclscan_read_count'}, axis=1, inplace=True)
            nclscan.drop(['chrom','start', 'end','strand'], axis = 1,inplace=True)
            nclscan = _df_setcol_as_str(nclscan,['nclscan_annotation'])
            nclscan.loc[nclscan['nclscan_annotation']=="1", 'nclscan_annotation'] = "Intragenic"
            nclscan.loc[nclscan['nclscan_annotation']=="0", 'nclscan_annotation'] = "Intergenic"
            # nclscan.loc[nclscan['nclscan_annotation']!="0" and nclscan['nclscan_annotation']!="1" , 'nclscan_annotation'] = "Unknown"
            nclscan.set_index(['circRNA_id'],inplace=True)

            nclscan.fillna(value=-1,inplace=True)

            intcols = [ 'nclscan_read_count' ]
            strcols = list ( set(nclscan.columns) - set(intcols) )
            nclscan = _df_setcol_as_int(nclscan,intcols)
            if len(strcols) > 0: nclscan = _df_setcol_as_str(nclscan,strcols)  

        dfs.append(nclscan)

    if args.circrnafinder:
        circrnafinder=pandas.read_csv(args.circrnafinder,sep="\t",header=0)
        # output circrnafinder table has the following columns:
        # | # | ColName              | Eg.              |
        # |---|----------------------|------------------|
        # | 1 | chr                  | chr1             |
        # | 2 | start                | 1223244          |
        # | 3 | end                  | 1223968          |
        # | 4 | strand               | -                |
        # | 5 | read_count           | 26               |
        circrnafinder['circRNA_id']=circrnafinder['chr'].astype(str)+"##"+circrnafinder['start'].astype(str)+"##"+circrnafinder['end'].astype(str)+"##"+circrnafinder['strand'].astype(str)
        circrnafinder.rename({'read_count': 'circrnafinder_read_count'}, axis=1, inplace=True)
        circrnafinder.set_index(['circRNA_id'],inplace=True)

        circrnafinder.fillna(value=-1,inplace=True)

        intcols = [ 'circrnafinder_read_count' ]
        strcols = list ( set(circrnafinder.columns) - set(intcols) )
        circrnafinder = _df_setcol_as_int(circrnafinder,intcols)
        if len(strcols) > 0: circrnafinder = _df_setcol_as_str(circrnafinder,strcols)    

        dfs.append(circrnafinder)


    # for df in dfs:
    #     print(df.columns)

    merged_counts=pandas.concat(dfs,axis=1,join="outer",sort=False)
    merged_counts['circRNA_id']=merged_counts.index
    merged_counts.fillna(-1,inplace=True)
    # print(merged_counts.columns)
    merged_counts[ 'ntools'] = 0

    annotation_cols=['circExplorer_annotation','ciri_annotation']
    floatcols = []
    intcols = [ 'circExplorer_read_count', 
                'circExplorer_found_BSJcounts', 
                'circExplorer_found_linear_BSJ_same_strand_counts', 
                'circExplorer_found_linear_spliced_BSJ_same_strand_counts', 
                'circExplorer_found_linear_BSJ_opposite_strand_counts', 
                'circExplorer_found_linear_spliced_BSJ_opposite_strand_counts', 
                'circExplorer_found_linear_BSJ_unknown_strand_counts', 
                'circExplorer_found_linear_spliced_BSJ_unknown_strand_counts' ]

    intcols.extend([ 'ciri_read_count',
                'ciri_linear_read_count' ])
    
    if args.dcc:
        intcols.extend([ 'dcc_read_count',
                    'dcc_linear_read_count' ])
        annotation_cols.extend(['dcc_gene','dcc_junction_type','dcc_annotation'])
    
    if args.mapsplice:
        intcols.extend([ 'mapsplice_read_count' ])
        floatcols.extend([ 'mapsplice_entropy' ])
        annotation_cols.extend(['mapsplice_annotation'])

    if args.nclscan and includenclscan:
        intcols.extend([ 'nclscan_read_count' ])
        annotation_cols.extend(['nclscan_annotation'])
    
    if args.circrnafinder:
        intcols.extend(['circrnafinder_read_count'])

    intcols.extend(['ntools'])
    strcols = list ( ( set(merged_counts.columns) - set(intcols) ) - set(floatcols) )
    merged_counts = _df_setcol_as_int(merged_counts,intcols)
    if len(floatcols)>0: merged_counts = _df_setcol_as_float(merged_counts,floatcols)
    merged_counts = _df_setcol_as_str(merged_counts,strcols)

    # fix annotations == -1
    for c in annotation_cols:
        merged_counts.loc[merged_counts[c]=="-1" , c] = "Unknown"

    merged_counts.loc[merged_counts['circExplorer_read_count'] >= args.minreads, 'ntools'] += 1
    merged_counts.loc[merged_counts['ciri_read_count'] >= args.minreads, 'ntools'] += 1
    if args.dcc: merged_counts.loc[merged_counts['dcc_read_count'] >= args.minreads, 'ntools'] += 1
    if args.mapsplice: merged_counts.loc[merged_counts['mapsplice_read_count'] >= args.minreads, 'ntools'] += 1
    if args.nclscan and includenclscan: merged_counts.loc[merged_counts['nclscan_read_count'] >= args.minreads, 'ntools'] += 1
    if args.circrnafinder: merged_counts.loc[merged_counts['circrnafinder_read_count'] >= args.minreads, 'ntools'] += 1
    merged_counts[['chrom', 'start', 'end', 'strand']] = merged_counts['circRNA_id'].str.split('##', expand=True)
 
    merged_counts=_df_setcol_as_int(merged_counts,['start','end','ntools'])
    merged_counts=_df_setcol_as_str(merged_counts,['chrom','strand'])
    # merged_counts.drop(['circRNA_id'],axis=1,inplace=True)
    # merged_counts['circRNA_id']=merged_counts['chrom'].astype(str)+":"+merged_counts['start'].astype(str)+"-"+merged_counts['end'].astype(str)

    # adding flanking sites
    merged_counts['flanking_sites']="-1"

    sequences = dict((s[1], s[0]) for s in HTSeq.FastaReader(args.reffa, raw_iterator=True))
    for index, row in merged_counts.iterrows():
        # print(index,row)
        bsj = BSJ(chrom=row['chrom'],start=row['start'],end=row['end'],strand=row['strand'])
        bsj.add_flanks(sequences)
        # print(bsj.get_flanks())
        # merged_counts[index]['flanking_sites'] = bsj.get_flanks()
        merged_counts.loc[index, 'flanking_sites'] = bsj.get_flanks()

    # add samplename
    merged_counts['sample_name'] = args.samplename
    merged_counts=_df_setcol_as_str(merged_counts,['sample_name','flanking_sites'])

    # prepare output
    outcols=['chrom', 'start', 'end', 'strand', 'flanking_sites', 'sample_name', 'ntools']
    # add circExplorer columns
    outcols.extend(['circExplorer_read_count',
                    'circExplorer_found_BSJcounts', 
                    'circExplorer_found_linear_BSJ_same_strand_counts', 
                    'circExplorer_found_linear_spliced_BSJ_same_strand_counts', 
                    'circExplorer_found_linear_BSJ_opposite_strand_counts', 
                    'circExplorer_found_linear_spliced_BSJ_opposite_strand_counts', 
                    'circExplorer_found_linear_BSJ_unknown_strand_counts', 
                    'circExplorer_found_linear_spliced_BSJ_unknown_strand_counts'])
    # add ciri columns
    outcols.extend(['ciri_read_count',
                    'ciri_linear_read_count'])
    # add DCC columns
    if args.dcc: 
        outcols.extend(['dcc_read_count',
                        'dcc_linear_read_count'])
    # add MapSplice columns
    if args.mapsplice: outcols.append('mapsplice_read_count')
    # add NCLscan columns
    if args.nclscan and includenclscan: outcols.append('nclscan_read_count')
    # add circRNAfinder columns
    if args.circrnafinder: outcols.append('circrnafinder_read_count')
    # add annotation columns
    outcols.extend(annotation_cols)
    merged_counts = merged_counts[outcols]
    merged_counts.to_csv(args.outfile,sep="\t",header=True,index=False,compression='gzip')


if __name__ == "__main__":
    main()