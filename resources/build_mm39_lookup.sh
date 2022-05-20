#!/usr/bin/bash

# download files from circBase and chain files for liftover

#wget http://circbase.org/download/mouse_mm9_circRNAs_putative_spliced_sequence.fa.gz
#wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz
#wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm39.over.chain.gz
#wget http://circbase.org/download/mm9_circID_to_name.txt
#wget http://circbase.org/download/mmu_mm9_circRNA.bed
#wget http://circbase.org/download/mmu_mm9_circRNA.txt

# create mID to genename lookup

zcat mouse_mm9_circRNAs_putative_spliced_sequence.fa.gz |grep "^>"|awk -F"|" -v OFS="\t" '{print $1,$3,$4}'|sed "s/>//g" > mID_to_gene.tsv

# create bed6 file from putative circRNA sequences

zcat mouse_mm9_circRNAs_putative_spliced_sequence.fa.gz |grep "^>"|awk -F"|" -v OFS="\t" '{print $2,$1}'|sed "s/>//g" > tmp
awk -v OFS="\t" '{x=substr($1,1,length($1)-1);y=substr($1,length($1));print x,$2,"0",y}' tmp > tmp2
cut -f1 tmp2|sed "s/:/\t/g"|sed "s/-/\t/g" > tmp3
cut -f2- tmp2 > tmp4
paste tmp3 tmp4 > mm9.circBase.bed
rm -f tmp tmp2 tmp3 tmp4

# crossmap the mm9 coordinates to mm10 and then to mm39

crossmap bed mm9ToMm10.over.chain.gz mm9.circBase.bed mm10.circBase.bed
crossmap bed mm10ToMm39.over.chain.gz mm10.circBase.bed mm39.circBase.bed

# collapse the mm39 bed file by name
# idential regions may have different names in column4... we need to collapse them into a comma-separated list in column 4

python collapse_bed_by_names.py mm39.circBase.bed mm39.circBase.collapsed.bed

# add mm9ID as a column

Rscript merge_dataframes.R

# this will produce "mm39_circBase_annotation_lookup.txt"
