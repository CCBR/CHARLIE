#!/bin/bash
fasta=$1
outdir=$2
mkdir -p $outdir && cd $outdir
cat $fasta | awk '{
        if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fasta")}
                        print $0 >> filename
                                close(filename)
                        }'