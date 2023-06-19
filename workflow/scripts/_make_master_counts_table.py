import pandas as pd
import argparse

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

def main() :
    parser = argparse.ArgumentParser(description='Make Master Counts Table with circExplorer_BWA fixes')
    parser.add_argument('--counttablelist', dest='counttablelist', type=str, required=True,
        help='comma separted list of per sample counts tables to merge')
    parser.add_argument('--minreads', dest='minreads', type=int, required=False, default=3,
        help='min read filter')
    parser.add_argument('-o',dest='outfile',required=True,help='master counts table')
    args = parser.parse_args()

    infiles = args.counttablelist
    infiles = infiles.split(",")
    count = 0
    for f in infiles:
        count += 1
        if count==1:
            outdf = pd.read_csv(f,sep="\t",header=0,compression='gzip')
            outdf.set_index(['chrom', 'start', 'end', 'strand', 'flanking_sites', 'sample_name'])
        else:
            tmpdf = pd.read_csv(f,sep="\t",header=0,compression='gzip')
            tmpdf.set_index(['chrom', 'start', 'end', 'strand', 'flanking_sites', 'sample_name'])
            outdf = pd.concat([outdf , tmpdf],axis=0,join="outer",sort=False)
    outdf.reset_index(drop=True,inplace=True)
    outdf.fillna(-1,inplace=True)
    # print(outdf.columns)
    intcols=['start','end','ntools']
    for c in outdf.columns:
        if "count" in c:
            intcols.append(c)
    # print(intcols)
    strcols=list(set(outdf.columns)-set(intcols))
    # print(strcols)
    outdf = _df_setcol_as_int(outdf,intcols)
    outdf = _df_setcol_as_str(outdf,strcols)
    outdf = outdf.sort_values(by=['chrom','start','end','strand'])

    # Applying circExplorer_BWA fixes
    indf=outdf
    non_circExplorer_BWA_df=indf.loc[(indf['strand']!=".")]
    ambiguous_strand_df=indf.loc[(indf['strand']==".") & (indf['circExplorer_bwa_read_count']>=args.minreads)]
    # print(ambiguous_strand_df.head())

    countlist = set(['circExplorer_read_count','ciri_read_count','findcirc_read_count','dcc_read_count','circrnafinder_read_count','mapsplice_read_count','nclscan_read_count'])
    colnames = set(indf.columns)
    excludelist = countlist - colnames
    newcountlist = countlist - excludelist
    # print("-2")
    # non_circExplorer_BWA_df = non_circExplorer_BWA_df.reset_index()
    lookup = dict()
    for index,row in non_circExplorer_BWA_df.iterrows():
        circID = row['chrom'] + "##" + str(row['start']) + "##" + str(row['end'])
        if not circID in lookup:
            lookup[circID] = dict()
            lookup[circID]["+"] = dict()
            lookup[circID]["+"]['sum'] = 0
            lookup[circID]["-"] = dict()
            lookup[circID]["-"]['sum'] = 0
        for c in newcountlist:
            if int(row[c]) != -1: 
                lookup[circID][row['strand']]['sum'] += int(row[c])
                lookup[circID][row['strand']]['flanking_sites'] = row['flanking_sites']
    # print("-1")
    # find winner
    for k,v in lookup.items():
        if v["+"]['sum'] >= v["-"]['sum']: 
            lookup[k]['strand'] = "+"
        else:
            lookup[k]['strand'] = "-"
    # print("0")
    # ambiguous_strand_df = ambiguous_strand_df.reset_index()
    for index,row in ambiguous_strand_df.iterrows():
        circID = row['chrom'] + "##" + str(row['start']) + "##" + str(row['end'])
        if not circID in lookup: continue
        ambiguous_strand_df.at[index,'strand']=lookup[circID]['strand']
        ambiguous_strand_df.at[index,'flanking_sites']=lookup[circID][lookup[circID]['strand']]['flanking_sites']

    non_ambiguous_strand_df = ambiguous_strand_df.loc[ambiguous_strand_df['strand']!="."]
    # print("1")
    bwa_counts = dict()
    circID2index = dict()
    circIDlist=list()
    for index2,row2 in non_ambiguous_strand_df.iterrows():
        circID2 = row2['chrom'] + "##" + str(row2['start']) + "##" + str(row2['end']) + "##" + row2['sample_name']
        circIDlist.append(circID2)
        bwa_counts[circID2] = row2['circExplorer_bwa_read_count']
        circID2index[circID2] = index2

    non_ambiguous_strand_df.insert(0,'circID',circIDlist)

    # print("2")
    claimed=dict()
    for index,row in non_circExplorer_BWA_df.iterrows():
        circID = row['chrom'] + "##" + str(row['start']) + "##" + str(row['end']) + "##" + row['sample_name']
        if circID in bwa_counts:
            claimed[circID]=1
            non_circExplorer_BWA_df.at[index,'circExplorer_bwa_read_count']=bwa_counts[circID]
    # print("3")
    set1=set(bwa_counts.keys())
    set2=set(claimed.keys())
    unclaimed=list(set1-set2)
    # print("4")
    # unclaimed_index=list()
    # for c in unclaimed:
    #     unclaimed_index.append(circID2index[c])
    if len(unclaimed) > 0:
        # unclaimed_non_ambiguous_strand_df = non_ambiguous_strand_df.iloc[unclaimed_index]
        unclaimed_non_ambiguous_strand_df = non_ambiguous_strand_df.loc[non_ambiguous_strand_df['circID'].isin(unclaimed)]
        unclaimed_non_ambiguous_strand_df_wo_circID = unclaimed_non_ambiguous_strand_df.drop(['circID'],axis=1)
        outdf=pd.concat([non_circExplorer_BWA_df,unclaimed_non_ambiguous_strand_df_wo_circID])
    else:
        outdf=non_circExplorer_BWA_df
    # print("5")

    # countlist = set(['circExplorer_read_count','ciri_read_count','findcirc_read_count','dcc_read_count','circrnafinder_read_count','mapsplice_read_count','nclscan_read_count'])
    # newcountlist
    hqcountlist = set(['circExplorer_read_count','circExplorer_bwa_read_count'])
    nonhqcountlist = newcountlist-hqcountlist

    for index,row in outdf.iterrows():
        hqcount=0
        nonhqcount=0
        for c in hqcountlist:
            if row[c] >= args.minreads: hqcount+=1
        for c in nonhqcountlist:
            if row[c] >= args.minreads: nonhqcount+=1
        outdf.at[index,'hqcounts']=hqcount
        outdf.at[index,'nonhqcounts']=nonhqcount
        outdf.at[index,'ntools']=hqcount+nonhqcount
        if hqcount==len(hqcountlist) and nonhqcount >= 1: 
            outdf.at[index,'HQ']="Y"
        else:
            outdf.at[index,'HQ']="N"

    intcols=['start','end','ntools']
    for c in outdf.columns:
        if "count" in c:
            intcols.append(c)
    strcols=list(set(outdf.columns)-set(intcols))
    outdf = _df_setcol_as_int(outdf,intcols)
    outdf = _df_setcol_as_str(outdf,strcols)
    outdf = outdf.sort_values(by=['chrom','start','end','strand'])
    outdf.to_csv(args.outfile,sep="\t",header=True,index=False,compression='gzip')

if __name__ == "__main__":
    main()