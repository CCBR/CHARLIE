import argparse
import sys
import pandas

# output file columns
# 1       #chrom
# 2       start
# 3       end
# 4       strand
# 5       known_novel
# 6       expected_BSJ_reads
# 7       found_BSJ_reads
# 8       linear_BSJ_reads_same_strand
# 9       linear_spliced_BSJ_reads_same_strand
# 10      linear_BSJ_reads_opposite_strand
# 11      linear_spliced_BSJ_reads_opposite_strand


def _df_setcol_as_int(df,collist):
    for c in collist:
        df[[c]]=df[[c]].astype(int)
    return df

def _df_setcol_as_str(df,collist):
    for c in collist:
        df[[c]]=df[[c]].astype(str)
    return df

def main():
    # debug = True
    debug = False
    parser = argparse.ArgumentParser(
    )
    parser.add_argument("--annotationcounts",dest="annotationcounts",required=True,type=str,
        help="annotated_counts.tsv counts file")
    parser.add_argument("--allfoundcounts",dest="allfoundcounts",required=True,type=str,
        help="readcounts.tsv")
    parser.add_argument("--countstable",dest="mergedcounts",required=True,type=str,
        help="merged counts_table.tsv file")
    args = parser.parse_args()

    bcounts = pandas.read_csv(args.annotationcounts,header=0,sep="\t")
    lcounts = pandas.read_csv(args.allfoundcounts,header=0,sep="\t")
    # print(bcounts.head())
    # print(lcounts.head())
    mcounts = bcounts.merge(lcounts,how='outer',on=["#chrom","start","end","strand"])
    mcounts.fillna(value=0,inplace=True)
    strcols = [ '#chrom', 'strand', 'known_novel' ]
    intcols = list ( set(mcounts.columns) - set(strcols) )
    mcounts = _df_setcol_as_str(mcounts,strcols)
    mcounts = _df_setcol_as_int(mcounts,intcols)
    mcounts.drop(["read_count"],axis=1,inplace=True)
    mcounts.to_csv(args.mergedcounts,index=False,doublequote=False,sep="\t")

if __name__ == "__main__":
    main()