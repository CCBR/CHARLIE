#!/bin/bash
# absolute path of input FASTA
fasta=$1
# absolute path of output folder
outdir=$2
if [ -d $outdir ];then rm -rf $outdir;fi
mkdir -p $outdir && cd $outdir
cat $fasta | awk '{
        if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fa")}
                        print $0 >> filename
                                close(filename)
                        }'