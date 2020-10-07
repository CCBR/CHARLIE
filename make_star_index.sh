module load STAR/2.7.0f
mkdir -p ./STAR_index_no_GTF &&
STAR \
--runThreadN 56 \
--runMode genomeGenerate \
--genomeDir ./STAR_index_no_GTF \
--genomeFastaFiles ./ref.fa