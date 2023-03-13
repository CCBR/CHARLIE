#!/usr/bin/env bash
# requires the following scripts
    # filter_bam.py
    # bam_get_max_readlen.py
    # _bedintersect_to_rid2jid.py
# also requires these tools
    # bedtools
    # python3 with pysam library
    # samtools
    # bedSort from ucsc tools

if [ "$#" != 5 ];then
    echo "5 arguments expected!"
    echo "#1 --> path to non-chimeric BAM file"
    echo "#2 --> sample name"
    echo "#3 --> PE or SE"
    echo "#4 --> known BSJs in bed.gz format"
    echo "#5 --> tmpdir"
    echo "OUTPUT: ${sample_name}.rid2jid.tsv.gz is created!"
    exit
fi
set -exo pipefail


#INPUTs
non_chimeric_star2p_bam=$1
sample_name=$2
peorse=$3
bsjbedgz=$4
tmpdir=$5
#OUTPUTs
# ${tmpdir}/${sample_name}.rid2jid.tsv.gz is created!

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )


# Filter all alignments other than the primary alignment

# if [ "$peorse" == "PE" ];then

# python3 ${SCRIPT_DIR}/filter_bam.py \
#     --inbam $non_chimeric_star2p_bam \
#     --outbam ${tmpdir}/${sample_name}.bam \
#     --pe

# else

# python3 ${SCRIPT_DIR}/filter_bam.py \
#     --inbam $non_chimeric_star2p_bam \
#     --outbam ${tmpdir}/${sample_name}.bam \

# fi

# sort primary alignments by query name

# samtools sort -n -l 1 -T $tmpdir -@4 -O BAM -o ${tmpdir}/${sample_name}.qsorted.bam ${tmpdir}/${sample_name}.bam

# if [ "$peorse" == "PE" ];then

# # convert primary alignments to BEDPE format

# bedtools bamtobed -bedpe -i ${tmpdir}/${sample_name}.qsorted.bam > ${tmpdir}/${sample_name}.bedpe
# # grab the START and END coordinates of the entire read-pair aligned
# python3 ${SCRIPT_DIR}/_bedpe2bed.py -i ${tmpdir}/${sample_name}.bedpe -o ${tmpdir}/${sample_name}.bed

# else

# bedtools bamtobed -i ${tmpdir}/${sample_name}.qsorted.bam > ${tmpdir}/${sample_name}.bed

# fi

# bedSort ${tmpdir}/${sample_name}.bed ${tmpdir}/${sample_name}.bed

# fix the max readlength in all alignments

python3 ${SCRIPT_DIR}/bam_get_max_readlen.py -i $non_chimeric_star2p_bam > ${tmpdir}/${sample_name}.maxrl
maxrl=$(cat ${tmpdir}/${sample_name}.maxrl)
two_maxrl=$((maxrl*2))
# two_maxrl="302"

## interect known BSJs with primary alignments and then keep alignments with some alignment within the "inclusion zone"

bedtools intersect -nonamecheck -wa -wb -a $bsjbedgz -b ${tmpdir}/${sample_name}.bed | \
python3 ${SCRIPT_DIR}/_bedintersect_to_rid2jid.py -i - -o ${sample_name}.rid2jid.tsv.gz -m $two_maxrl