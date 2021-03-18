#!/bin/bash
bam=$1
region=$2
samtools view -@4 $bam $region 2>/dev/null|cut -f1|sort -T /dev/shm|uniq|wc -l