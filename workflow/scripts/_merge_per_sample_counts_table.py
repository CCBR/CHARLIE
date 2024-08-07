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

def _rev_comp(seq):
    seq = seq.upper()
    seq = seq.replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c")
    seq = seq.upper()[::-1]
    return seq

class BSJ:
    def __init__(self,chrom,start,end,strand="+"):
        self.chrom=chrom
        self.start=int(start)
        self.end=int(end)
        self.strand=strand
        self.splice_site_flank_5="" #donor
        self.splice_site_flank_3="" #acceptor
    
    def add_flanks(self,sequences): # adds flanking assuming + strand
        coord = int(self.end)
        seq = sequences[self.chrom][coord:coord+2]
        self.splice_site_flank_5 = seq.upper()
        coord = int(self.start)
        seq = sequences[self.chrom][coord-2:coord]
        self.splice_site_flank_3 = seq.upper()
    
    def get_flanks(self): # returns + and - strand flanks
        plus_strand = self.splice_site_flank_5+"##"+self.splice_site_flank_3
        minus_strand = _rev_comp(self.splice_site_flank_5)+"##"+_rev_comp(self.splice_site_flank_3)
        return plus_strand,minus_strand

def main() :
    parser = argparse.ArgumentParser(description='Merge per sample Counts from different circRNA detection tools')
    parser.add_argument('--circExplorer', dest='circE', type=str, required=True,
        help='circExplorer2 per-sample counts table')
    parser.add_argument('--circExplorerbwa', dest='circEbwa', type=str, required=True,
        help='circExplorer2_bwa per-sample counts table')
    parser.add_argument('--ciri', dest='ciri', type=str, required=True,
        help='ciri2 per-sample output')
    parser.add_argument('--findcirc', dest='findcirc', type=str, required=False,
        help='findcirc per-sample counts table')
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
    parser.add_argument('--hqcc', dest='hqcc', type=str, required=False, default="circExplorer,circExplorer_bwa",
        help='Comma separated list of high confidence core callers (default="circExplorer,circExplorer_bwa")')
    parser.add_argument('--hqccpn', dest='hqccpn', type=int, required=False, default=1,
        help='Define n:high confidence core callers plus n callers are required to call this circRNA HQ (default 1)')
    parser.add_argument('-o',dest='outfile',required=True,help='merged table')
    args = parser.parse_args()

    sn=args.samplename
    hqcc=args.hqcc
    hqcc=hqcc.strip().lower().split(",")
    hqcclen=len(hqcc)
    required_hqcols=[]
    not_required_hqcols=[]
    dfs=[]

    # load circExplorer
    circE=pandas.read_csv(args.circE,sep="\t",header=0)
    print(circE.columns)
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
    # | 8  | linear_+                       |
    # | 9  | spliced_+                      |
    # | 10 | linear_-                   |
    # | 11 | spliced_-                  |
    # | 12 | linear_.                    |
    # | 13 | spliced_.                   |
    circE['circRNA_id']=circE['#chrom'].astype(str)+"##"+circE['start'].astype(str)+"##"+circE['end'].astype(str)
    circE.rename({'strand'              : 'circExplorer_strand',
                  'known_novel'         : 'circExplorer_annotation',
                'expected_BSJ_reads'    : 'circExplorer_read_count',
                'found_BSJ_reads'       : 'circExplorer_found_BSJcounts',
                'linear_+'              : 'circExplorer_found_linear_BSJ_+_counts',
                'spliced_+'             : 'circExplorer_found_linear_spliced_BSJ_+_counts',
                'linear_-'              : 'circExplorer_found_linear_BSJ_-_counts',
                'spliced_-'             : 'circExplorer_found_linear_spliced_BSJ_-_counts',
                'linear_.'              : 'circExplorer_found_linear_BSJ_._counts',
                'spliced_.'             : 'circExplorer_found_linear_spliced_BSJ_._counts'}, axis=1, inplace=True)
    circE.drop(['#chrom','start', 'end'], axis = 1,inplace=True)
    circE.set_index(['circRNA_id'],inplace=True,drop=True)
    
    circE.fillna(value=-1,inplace=True)
    print(circE.columns)

    intcols = [ 'circExplorer_read_count', 
                'circExplorer_found_BSJcounts', 
                'circExplorer_found_linear_BSJ_+_counts', 
                'circExplorer_found_linear_spliced_BSJ_+_counts', 
                'circExplorer_found_linear_BSJ_-_counts', 
                'circExplorer_found_linear_spliced_BSJ_-_counts', 
                'circExplorer_found_linear_BSJ_._counts', 
                'circExplorer_found_linear_spliced_BSJ_._counts' ]
    strcols = list ( set(circE.columns) - set(intcols) )
    circE = _df_setcol_as_int(circE,intcols)
    circE = _df_setcol_as_str(circE,strcols)

    dfs.append(circE)
    if "circExplorer".lower() in hqcc: 
        required_hqcols.append("circExplorer_read_count")
    else:
        not_required_hqcols.append("circExplorer_read_count")

    #chrom  start   end     strand  read_count      known_novel
    # circExplorer2 with BWA

    circEbwa=pandas.read_csv(args.circEbwa,sep="\t",header=0)
    circEbwa['circRNA_id']=circEbwa['#chrom'].astype(str)+"##"+circEbwa['start'].astype(str)+"##"+circEbwa['end'].astype(str)
    circEbwa.rename({'strand'       : 'circExplorer_bwa_strand',
                     'known_novel'  : 'circExplorer_bwa_annotation',
                     'read_count'   : 'circExplorer_bwa_read_count'}, axis=1, inplace=True)
    circEbwa.drop(['#chrom','start', 'end'], axis = 1,inplace=True)
    circEbwa.set_index(['circRNA_id'],inplace=True,drop=True)
    
    circEbwa.fillna(value=-1,inplace=True)

    intcols = [ 'circExplorer_bwa_read_count' ]
    strcols = list ( set(circEbwa.columns) - set(intcols) )
    circEbwa = _df_setcol_as_int(circEbwa,intcols)
    circEbwa = _df_setcol_as_str(circEbwa,strcols)

    dfs.append(circEbwa)
    if "circExplorer_bwa".lower() in hqcc: 
        required_hqcols.append("circExplorer_bwa_read_count")
    else:
        not_required_hqcols.append("circExplorer_bwa_read_count")

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
    ciri['circRNA_id']=ciri['chr'].astype(str)+"##"+ciri['circRNA_start'].astype(str)+"##"+ciri['circRNA_end'].astype(str)
    ciri.rename({   'strand'                : 'ciri_strand',
                    '#junction_reads'       : 'ciri_read_count', 
                    '#non_junction_reads'   : 'ciri_linear_read_count', 
                    'circRNA_type'          : 'ciri_annotation'}, axis=1, inplace=True)
    ciri.drop(['chr','circRNA_start', 'circRNA_end'], axis = 1,inplace=True)
    ciri.set_index(['circRNA_id'],inplace=True,drop=True)

    ciri.fillna(value=-1,inplace=True)

    intcols = [ 'ciri_read_count',
                'ciri_linear_read_count' ]
    strcols = list ( set(ciri.columns) - set(intcols) )
    ciri = _df_setcol_as_int(ciri,intcols)
    if len(strcols) > 0: ciri = _df_setcol_as_str(ciri,strcols)
    
    dfs.append(ciri)
    if "ciri".lower() in hqcc: 
        required_hqcols.append("ciri_read_count")
    else:
        not_required_hqcols.append("ciri_read_count")

    if args.findcirc:
        findcirc=pandas.read_csv(args.findcirc,sep="\t",header=0)
        print(findcirc.columns)
# add find_circ
# | #  | short_name      | description
# | -- | --------------- | ---------------------------------------------------------------------------------------------------------------- |
# | 1  | chrom           | chromosome/contig name                                                                                           |
# | 2  | start           | left splice site (zero-based)                                                                                    |
# | 3  | end             | right splice site (zero-based). (Always: end > start. 5' 3' depends on strand)                                   |
# | 4  | name            | (provisional) running number/name assigned to junction                                                           |
# | 5  | n_reads         | number of reads supporting the junction (BED 'score')                                                            |
# | 6  | strand          | genomic strand (+ or -)                                                                                          |
# | 7  | n_uniq          | number of distinct read sequences supporting the junction                                                        |
# | 8  | uniq_bridges    | number of reads with both anchors aligning uniquely                                                              |
# | 9  | best_qual_left  | alignment score margin of the best anchor alignment supporting the left splice junction (max=2 \* anchor_length) |
# | 10 | best_qual_right | same for the right splice site                                                                                   |
# | 11 | tissues         | comma-separated, alphabetically sorted list of tissues/samples with this junction                                |
# | 12 | tiss_counts     | comma-separated list of corresponding read-counts                                                                |
# | 13 | edits           | number of mismatches in the anchor extension process                                                             |
# | 14 | anchor_overlap  | number of nucleotides the breakpoint resides within one anchor                                                   |
# | 15 | breakpoints     | number of alternative ways to break the read with flanking GT/AG                                                 |
# | 16 | signal          | flanking dinucleotide splice signal (normally GT/AG)                                                             |
# | 17 | strandmatch     | 'MATCH', 'MISMATCH' or 'NA' for non-stranded analysis                                                            |
# | 18 | category        | list of keywords describing the junction. Useful for quick grep filtering                                        |
        findcirc['circRNA_id']=findcirc['chrom'].astype(str)+"##"+findcirc['start'].astype(str)+"##"+findcirc['end'].astype(str)
        findcirc = findcirc.loc[:, ['circRNA_id', 'n_reads', 'strand']]
        findcirc.rename({ 'strand'      : 'findcirc_strand',
                          'n_reads'     : 'findcirc_read_count'}, axis=1, inplace=True)
        findcirc.set_index(['circRNA_id'],inplace=True,drop=True)

        findcirc.fillna(value=-1,inplace=True)

        intcols = [ 'findcirc_read_count' ]
        strcols = list ( set(findcirc.columns) - set(intcols) )
        findcirc = _df_setcol_as_int(findcirc,intcols)
        if len(strcols) > 0: findcirc = _df_setcol_as_str(findcirc,strcols)    
        dfs.append(findcirc)
        if "findcirc".lower() in hqcc: 
            required_hqcols.append("findcirc_read_count")
        else:
            not_required_hqcols.append("findcirc_read_count")

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
        dcc['circRNA_id']=dcc['chr'].astype(str)+"##"+dcc['start'].astype(str)+"##"+dcc['end'].astype(str)
        dcc.rename({'strand': 'dcc_strand'}, axis=1, inplace=True)
        dcc.rename({'read_count': 'dcc_read_count'}, axis=1, inplace=True)
        dcc.rename({'linear_read_count': 'dcc_linear_read_count'}, axis=1, inplace=True)
        dcc[['dcc_gene', 'dcc_junction_type', 'dcc_annotation2']] = dcc['dcc_annotation'].apply(lambda x: pandas.Series(str(x).split("##")))
        dcc.drop(['chr','start', 'end','dcc_annotation'], axis = 1,inplace=True)
        dcc.rename({'dcc_annotation2': 'dcc_annotation'}, axis=1, inplace=True)
        dcc.set_index(['circRNA_id'],inplace=True,drop=True)

        dcc.fillna(value=-1,inplace=True)

        intcols = [ 'dcc_read_count',
                    'dcc_linear_read_count' ]
        strcols = list ( set(dcc.columns) - set(intcols) )
        dcc = _df_setcol_as_int(dcc,intcols)
        if len(strcols) > 0: dcc = _df_setcol_as_str(dcc,strcols)    

        dfs.append(dcc)
        if "DCC".lower() in hqcc: 
            required_hqcols.append("dcc_read_count")
        else:
            not_required_hqcols.append("dcc_read_count")

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
        mapsplice['circRNA_id']=mapsplice['chrom'].astype(str)+"##"+mapsplice['start'].astype(str)+"##"+mapsplice['end'].astype(str)
        mapsplice.rename({'strand': 'mapsplice_strand'}, axis=1, inplace=True)
        mapsplice.rename({'read_count': 'mapsplice_read_count'}, axis=1, inplace=True)
        mapsplice[['mapsplice_annotation2', 'mapsplice_entropy']] = mapsplice['mapsplice_annotation'].apply(lambda x: pandas.Series(str(x).split("##")))
        mapsplice.drop(['chrom','start', 'end','mapsplice_annotation'], axis = 1,inplace=True)
        mapsplice.rename({'mapsplice_annotation2': 'mapsplice_annotation'}, axis=1, inplace=True)
        mapsplice.set_index(['circRNA_id'],inplace=True,drop=True)

        mapsplice.fillna(value=-1,inplace=True)

        intcols = [ 'mapsplice_read_count' ]
        mapsplice = _df_setcol_as_int(mapsplice,intcols)
        floatcols = [ 'mapsplice_entropy' ]
        mapsplice = _df_setcol_as_float(mapsplice,floatcols)
        strcols = list ( ( set(mapsplice.columns) - set(intcols) ) - set(floatcols) )
        if len(strcols) > 0: mapsplice = _df_setcol_as_str(mapsplice,strcols) 

        dfs.append(mapsplice)
        if "MapSplice".lower() in hqcc: 
            required_hqcols.append("mapsplice_read_count")
        else:
            not_required_hqcols.append("mapsplice_read_count")

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
            nclscan['circRNA_id']=nclscan['chrom'].astype(str)+"##"+nclscan['start'].astype(str)+"##"+nclscan['end'].astype(str)
            nclscan.rename({'strand': 'nclscan_strand'}, axis=1, inplace=True)
            nclscan.rename({'read_count': 'nclscan_read_count'}, axis=1, inplace=True)
            nclscan.drop(['chrom','start', 'end'], axis = 1,inplace=True)
            nclscan = _df_setcol_as_str(nclscan,['nclscan_annotation'])
            nclscan.loc[nclscan['nclscan_annotation']=="1", 'nclscan_annotation'] = "Intragenic"
            nclscan.loc[nclscan['nclscan_annotation']=="0", 'nclscan_annotation'] = "Intergenic"
            # nclscan.loc[nclscan['nclscan_annotation']!="0" and nclscan['nclscan_annotation']!="1" , 'nclscan_annotation'] = "Unknown"
            nclscan.set_index(['circRNA_id'],inplace=True,drop=True)

            nclscan.fillna(value=-1,inplace=True)

            intcols = [ 'nclscan_read_count' ]
            strcols = list ( set(nclscan.columns) - set(intcols) )
            nclscan = _df_setcol_as_int(nclscan,intcols)
            if len(strcols) > 0: nclscan = _df_setcol_as_str(nclscan,strcols)  

        dfs.append(nclscan)
        if "NCLscan".lower() in hqcc: 
            required_hqcols.append("nclscan_read_count")
        else:
            not_required_hqcols.append("nclscan_read_count")

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
        circrnafinder['circRNA_id']=circrnafinder['chr'].astype(str)+"##"+circrnafinder['start'].astype(str)+"##"+circrnafinder['end'].astype(str)
        circrnafinder.rename({'strand': 'circrnafinder_strand'}, axis=1, inplace=True)
        circrnafinder.rename({'read_count': 'circrnafinder_read_count'}, axis=1, inplace=True)
        circrnafinder.drop(['chr','start', 'end'], axis = 1,inplace=True)
        circrnafinder.set_index(['circRNA_id'],inplace=True,drop=True)

        circrnafinder.fillna(value=-1,inplace=True)

        intcols = [ 'circrnafinder_read_count' ]
        strcols = list ( set(circrnafinder.columns) - set(intcols) )
        circrnafinder = _df_setcol_as_int(circrnafinder,intcols)
        if len(strcols) > 0: circrnafinder = _df_setcol_as_str(circrnafinder,strcols)    

        dfs.append(circrnafinder)
        if "circRNAFinder".lower() in hqcc: 
            required_hqcols.append("circrnafinder_read_count")
        else:
            not_required_hqcols.append("circrnafinder_read_count")

    # for df in dfs:
    #     print(df.columns)


    # merged_counts=pandas.concat(dfs,axis=1,join="outer",sort=False)
    # merged_counts['circRNA_id']=merged_counts.index

# above pandas.concat not working as expected
# giving error
#   File "/vf/users/Ziegelbauer_lab/Pipelines/circRNA/230406_activeDev_20284a3/workflow/scripts/_merge_per_sample_counts_table.py", line 396, in <module>
#     main()
#   File "/vf/users/Ziegelbauer_lab/Pipelines/circRNA/230406_activeDev_20284a3/workflow/scripts/_merge_per_sample_counts_table.py", line 289, in main
#     merged_counts=pandas.concat(dfs,axis=1,join="outer",sort=False)
#   File "/usr/local/Anaconda/envs/py3.7/lib/python3.7/site-packages/pandas/util/_decorators.py", line 311, in wrapper
#     return func(*args, **kwargs)
#   File "/usr/local/Anaconda/envs/py3.7/lib/python3.7/site-packages/pandas/core/reshape/concat.py", line 307, in concat
#     return op.get_result()
#   File "/usr/local/Anaconda/envs/py3.7/lib/python3.7/site-packages/pandas/core/reshape/concat.py", line 528, in get_result
#     indexers[ax] = obj_labels.get_indexer(new_labels)
#   File "/usr/local/Anaconda/envs/py3.7/lib/python3.7/site-packages/pandas/core/indexes/base.py", line 3442, in get_indexer
#     raise InvalidIndexError(self._requires_unique_msg)
# pandas.errors.InvalidIndexError: Reindexing only valid with uniquely valued Index objects
# HENCE, replacing concat with this:

    for i,df in enumerate(dfs):
        if i==0:
            merged_counts=df
            merged_counts['circRNA_id']=merged_counts.index
            merged_counts.reset_index(inplace=True,drop=True)
        else:
            df['circRNA_id']=df.index
            df.reset_index(inplace=True,drop=True)
            merged_counts=pandas.merge(merged_counts,df,how='outer',on=['circRNA_id'])
    
    print(merged_counts.columns)
    
    # merged_counts.set_index(['circRNA_id'],inplace=True,drop=True)

    merged_counts.fillna(-1,inplace=True)
    merged_counts[ 'ntools'] = 0
    merged_counts[ 'HQ' ] = "N"
    merged_counts[ 'hqcounts' ] = 0
    merged_counts[ 'nonhqcounts' ] = 0

    annotation_cols=['circExplorer_annotation','ciri_annotation']
    floatcols = []
    strand_cols = ['circExplorer_strand','circExplorer_bwa_strand','ciri_strand']

    intcols = [ 'circExplorer_read_count', 
                'circExplorer_found_BSJcounts', 
                'circExplorer_found_linear_BSJ_+_counts', 
                'circExplorer_found_linear_spliced_BSJ_+_counts', 
                'circExplorer_found_linear_BSJ_-_counts', 
                'circExplorer_found_linear_spliced_BSJ_-_counts', 
                'circExplorer_found_linear_BSJ_._counts', 
                'circExplorer_found_linear_spliced_BSJ_._counts' ]

    intcols.extend([ 'ciri_read_count',
                'ciri_linear_read_count' ])
    
    intcols.extend(['circExplorer_bwa_read_count'])
    annotation_cols.extend(['circExplorer_bwa_annotation'])

    if args.findcirc:
        intcols.extend(['findcirc_read_count'])
        strand_cols.append('findcirc_strand')
    
    if args.dcc:
        intcols.extend([ 'dcc_read_count',
                    'dcc_linear_read_count' ])
        annotation_cols.extend(['dcc_gene','dcc_junction_type','dcc_annotation'])
        strand_cols.append('dcc_strand')
    
    if args.mapsplice:
        intcols.extend([ 'mapsplice_read_count' ])
        floatcols.extend([ 'mapsplice_entropy' ])
        annotation_cols.extend(['mapsplice_annotation'])
        strand_cols.append('mapsplice_strand')

    if args.nclscan and includenclscan:
        intcols.extend([ 'nclscan_read_count' ])
        annotation_cols.extend(['nclscan_annotation'])
        strand_cols.append('nclscan_strand')
    
    if args.circrnafinder:
        intcols.extend(['circrnafinder_read_count'])
        strand_cols.append('circrnafinder_strand')

    intcols.extend(['ntools'])
    intcols.extend(['hqcounts','nonhqcounts'])
    strcols = list ( ( set(merged_counts.columns) - set(intcols) ) - set(floatcols) )
    strcols.append('HQ')
    merged_counts = _df_setcol_as_int(merged_counts,intcols)
    if len(floatcols)>0: merged_counts = _df_setcol_as_float(merged_counts,floatcols)
    merged_counts = _df_setcol_as_str(merged_counts,strcols)

    # fix annotations == -1
    for c in annotation_cols:
        merged_counts.loc[merged_counts[c]=="-1" , c] = "Unknown"
    
    for c in required_hqcols:
        merged_counts.loc[merged_counts[c] >= args.minreads, 'hqcounts'] += 1
    for c in not_required_hqcols:
        merged_counts.loc[merged_counts[c] >= args.minreads, 'nonhqcounts'] += 1
    
    merged_counts.loc[merged_counts['hqcounts'] == hqcclen, 'HQ'] = "Y"
    merged_counts.loc[merged_counts['nonhqcounts'] < args.hqccpn, 'HQ'] = "N"

    merged_counts.loc[merged_counts['circExplorer_read_count'] >= args.minreads, 'ntools'] += 1
    merged_counts.loc[merged_counts['ciri_read_count'] >= args.minreads, 'ntools'] += 1
    merged_counts.loc[merged_counts['circExplorer_bwa_read_count'] >= args.minreads, 'ntools'] += 1
    if args.findcirc: merged_counts.loc[merged_counts['findcirc_read_count'] >= args.minreads, 'ntools'] += 1
    if args.dcc: merged_counts.loc[merged_counts['dcc_read_count'] >= args.minreads, 'ntools'] += 1
    if args.mapsplice: merged_counts.loc[merged_counts['mapsplice_read_count'] >= args.minreads, 'ntools'] += 1
    if args.nclscan and includenclscan: merged_counts.loc[merged_counts['nclscan_read_count'] >= args.minreads, 'ntools'] += 1
    if args.circrnafinder: merged_counts.loc[merged_counts['circrnafinder_read_count'] >= args.minreads, 'ntools'] += 1
    merged_counts[['chrom', 'start', 'end']] = merged_counts['circRNA_id'].str.split('##', expand=True)
 
    merged_counts=_df_setcol_as_int(merged_counts,['start','end','ntools'])
    merged_counts=_df_setcol_as_str(merged_counts,['chrom'])

    # adding flanking sites
    merged_counts['flanking_sites_+']="-1"
    merged_counts['flanking_sites_-']="-1"

    sequences = dict((s[1], s[0]) for s in HTSeq.FastaReader(args.reffa, raw_iterator=True))
    for index, row in merged_counts.iterrows():
        bsj = BSJ(chrom=row['chrom'],start=row['start'],end=row['end'])
        bsj.add_flanks(sequences)
        plus_flank, minus_flank = bsj.get_flanks()
        merged_counts.loc[index, 'flanking_sites_+'] = plus_flank
        merged_counts.loc[index, 'flanking_sites_-'] = minus_flank

    # add samplename
    merged_counts['sample_name'] = args.samplename
    merged_counts=_df_setcol_as_str(merged_counts,['sample_name','flanking_sites_+','flanking_sites_-'])
    print(merged_counts.columns)

    # prepare output ... reorder columns
    outcols=['chrom', 'start', 'end']
    outcols.extend(strand_cols)
    outcols.extend(['flanking_sites_+','flanking_sites_-', 'sample_name', 'ntools', 'HQ'])
    # add circExplorer columns
    outcols.extend(['circExplorer_read_count',
                    'circExplorer_found_BSJcounts', 
                    'circExplorer_found_linear_BSJ_+_counts', 
                    'circExplorer_found_linear_spliced_BSJ_+_counts', 
                    'circExplorer_found_linear_BSJ_-_counts', 
                    'circExplorer_found_linear_spliced_BSJ_-_counts', 
                    'circExplorer_found_linear_BSJ_._counts', 
                    'circExplorer_found_linear_spliced_BSJ_._counts'])
    # add ciri columns
    outcols.extend(['ciri_read_count',
                    'ciri_linear_read_count'])
    # add circExplorer_BWA columns
    outcols.extend(['circExplorer_bwa_read_count'])
    # add find_circ columns
    if args.findcirc:
        outcols.extend(['findcirc_read_count'])
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

    outcols.extend(['hqcounts','nonhqcounts'])
    # add annotation columns
    outcols.extend(annotation_cols)
    merged_counts = merged_counts[outcols]
    merged_counts.to_csv(args.outfile,sep="\t",header=True,index=False,compression='gzip')


if __name__ == "__main__":
    main()
