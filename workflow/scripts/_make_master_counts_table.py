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
    parser = argparse.ArgumentParser(description='Make Master Counts Table')
    parser.add_argument('--counttablelist', dest='counttablelist', type=str, required=True,
        help='comma separted list of per sample counts tables to merge')
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
    print(outdf.columns)
    intcols=['start','end','ntools']
    for c in outdf.columns:
        if "count" in c:
            intcols.append(c)
    print(intcols)
    strcols=list(set(outdf.columns)-set(intcols))
    print(strcols)
    outdf = _df_setcol_as_int(outdf,intcols)
    outdf = _df_setcol_as_str(outdf,strcols)
    outdf = outdf.sort_values(by=['chrom','start','end','strand'])
    outdf.to_csv(args.outfile,sep="\t",header=True,index=False,compression='gzip')

if __name__ == "__main__":
    main()