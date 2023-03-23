import pandas as pd
import argparse


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
    outdf.reset_index(inplace=True)
    outdf.fillna(-1,inplace=True)
    print(outdf.columns)
    outdf = outdf.sort_values(by=['chrom','start','end','strand'])
    outdf.to_csv(args.outfile,sep="\t",header=True,index=False,compression='gzip')

if __name__ == "__main__":
    main()