#!/usr/bin/env bash

# requires
# 1. bedtools
# 2. ucsc
# 3. samtools

# INPUT
bam=$1

bw="${bam%.*}.bw"
bdg="${bam%.*}.bdg"
sizes="${bam%.*}.sizes"
bedtools genomecov -bga -split -ibam $bam > $bdg
bedSort $bdg $bdg
if [ "$(wc -l $bdg|awk '{print $1}')" != "0" ];then
    bedSort $bdg $bdg
    samtools view -H $bam|grep ^@SQ|cut -f2,3|sed "s/SN://g"|sed "s/LN://g" > $sizes
    bedGraphToBigWig $bdg $sizes $bw
fi
rm -f $bdg $sizes