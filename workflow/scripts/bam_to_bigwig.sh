#!/usr/bin/env bash

# requires
# 1. bedtools
# 2. ucsc
# 3. samtools

# INPUT
bam=$1
if [ "$#" == 2 ];then
tmpdir=$2
else
tmpdir=$(dirname $bam)
fi
echo $tmpdir

bam_bn=$(basename $bam)
bw="${bam%.*}.bw"
bdg="${bam_bn%.*}.bdg"
sizes="${bam_bn%.*}.sizes"
bedtools genomecov -bga -split -ibam $bam > ${tmpdir}/${bdg}

if [ "$(wc -l ${tmpdir}/${bdg}|awk '{print $1}')" != "0" ];then
    bedSort ${tmpdir}/${bdg} ${tmpdir}/${bdg}
    samtools view -H $bam | grep ^@SQ | cut -f2,3 | sed "s/SN://g" | sed "s/LN://g" > ${tmpdir}/${sizes}
    bedGraphToBigWig ${tmpdir}/${bdg} ${tmpdir}/${sizes} $bw
fi
rm -f ${tmpdir}/${bdg} ${tmpdir}/${sizes}