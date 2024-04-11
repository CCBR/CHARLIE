#!/usr/bin/env bash
set -e -x -o pipefail

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# SCRIPT_DIR="/vf/users/Ziegelbauer_lab/Pipelines/circRNA/activeDev/workflow/scripts"
SCRIPT_NAME=$( basename -- "${BASH_SOURCE[0]}" )

ARGPARSE_DESCRIPTION="create linear and splice BAMs (and BIGWIGs)"
    source /opt2/argparse.bash || \
    source ${SCRIPT_DIR}/../../resources/argparse.bash || \
    source <(curl -s https://raw.githubusercontent.com/CCBR/Tools/master/scripts/argparse.bash) || \
    exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--nonchimericbam',required=True, help='path to non-chimeric BAM file')
parser.add_argument('--samplename',required=True, help='sample name')
parser.add_argument('--peorse',required=True, help='PE or SE')
parser.add_argument('--bsjbed',required=True, help='known BSJs in bed.gz format')
parser.add_argument('--tmpdir',required=False, default="/tmp", help='temp dir')
parser.add_argument('--rid2jid',required=True, help='path to output file eg. <sample_name>.rid2jid.tsv.gz')
parser.add_argument('--filteredbam',required=True, help='filtered output BAM')
parser.add_argument('--linearbsjlist',required=True, help='gzip-ed list of linear BSJ readids')
parser.add_argument('--splicedbsjlist',required=True, help='gzip-ed list of linear spliced BSJ readids')
parser.add_argument('--jidcounts',required=True, help='path to jid (linear and spliced) counts per jid or BSJ')
parser.add_argument('--linearbsjbam',required=True, help='path to linear BAM with only near BSJ reads')
parser.add_argument('--splicedbsjbam',required=True, help='path to spliced BAM with only near BSJ reads')
parser.add_argument('--linearbam',required=True, help='path to linear BAM with all linear reads')
parser.add_argument('--splicedbam',required=True, help='path to spliced BAM with all linear reads that are spliced')
parser.add_argument('--regions',required=True, help='path to .regions file eg. ref/ref.fa.regions')
parser.add_argument('--host',required=False, default="", help='comma separated list of HOST')
parser.add_argument('--additives',required=False, default="", help='comma separated list of ADDITIVES')
parser.add_argument('--viruses',required=False, default="", help='comma separated list of VIRUSES')
parser.add_argument('--threads',required=False, default=56, help='number of threads')
EOF

threads=$THREADS
if [ "$SLURM_JOB_ID" != "" ];then
alloccpu=$(sacct -j $SLURM_JOB_ID --format "AllocCPUS"|tail -n1|awk '{print $1}')
if [ "$alloccpu" -lt "$threads" ];then
    threads=$alloccpu
fi
fi

start0=$(date +%s.%N)

# requires the following scripts
    # filter_bam.py
    # bam_get_max_readlen.py
    # _bedintersect_to_rid2jid.py
# also requires these tools
    # bedtools
    # python3 with pysam library
    # samtools
    # bedSort from ucsc tools
    # ucsc

#if [ "$#" != 18 ];then
#    echo "$#"
#    echo "18 arguments expected!"
#    echo "#1 --> path to non-chimeric BAM file"
#    echo "#2 --> sample name"
#    echo "#3 --> PE or SE"
#    echo "#4 --> known BSJs in bed.gz format"
#    echo "#5 --> tmpdir"
#    echo "#6 --> gzipped outputfilename eg.${sample_name}.rid2jid.tsv.gz"
#    echo "#7 --> output filtered sample BAM"
#	echo "#8 --> gzip-ed list of linear BSJ readids"
#	echo "#9 --> gzip-ed list of linear spliced BSJ readids"
#	echo "#10 --> jid counts (linear and linear-spliced) per jid or BSJ"
#	echo "#11 --> linear BSJ readids in BAM"
#	echo "#12 --> linear-spliced BSJ readids in BAM"
#	echo "#13 --> .regions file eg. ref/ref.fa.regions"
#	echo "#14 --> host list comma separated .. no spaces"
#	echo "#15 --> additives list comma separated .. no spaces"
#       echo "#16 --> viruses list comma separated .. no spaces"
#	echo "#17 --> all linear readids in BAM"
#	echo "#18 --> all spliced readids in BAM"
#    exit
#fi

set -exo pipefail

function printtime() {
    scriptname=$1
    start0=$2
    start=$3
    msg=$4
    end=$(date +%s.%N)    
    runtime0=$(python -c "print(${end} - ${start0})")
    runtime0=${runtime0%.*}
    runtime=$(python -c "print(${end} - ${start})")
    runtime=${runtime%.*}
    echo "$scriptname | $runtime0 | $runtime | $msg"
}

###################################################################################################
# set inputs and outputs
###################################################################################################

#INPUTs
non_chimeric_star2p_bam=$NONCHIMERICBAM
sample_name=$SAMPLENAME
peorse=$PEORSE
bsjbedgz=$BSJBED
tmpdir=$TMPDIR

#OUTPUTs
rid2jidgzip=$RID2JID
filtered_bam=$FILTEREDBAM
linearrids=$LINEARBSJLIST
splicedrids=$SPLICEDBSJLIST
jidcounts=$JIDCOUNTS
linearbam=$LINEARBSJBAM
splicedbam=$SPLICEDBSJBAM
regions=$REGIONS
host=$HOST
additives=$ADDITIVES
viruses=$VIRUSES
linearbam_all=$LINEARBAM
splicedbam_all=$SPLICEDBAM
# ${tmpdir}/${sample_name}.rid2jid.tsv.gz is created!


outdir=$(dirname $linearbam)
###################################################################################################
## Filter all alignments other than the primary alignment
###################################################################################################

printtime $SCRIPT_NAME $start0 $start0 "filtering bam now!"
start=$(date +%s.%N)

# if [ "1" == "0" ];then

if [ "$peorse" == "PE" ];then

python3 ${SCRIPT_DIR}/filter_bam.py \
    --inbam $non_chimeric_star2p_bam \
    --outbam $filtered_bam \
    --pe

else

python3 ${SCRIPT_DIR}/filter_bam.py \
    --inbam $non_chimeric_star2p_bam \
    --outbam $filtered_bam \

fi

nreads=$(samtools view /gpfs/gsfs8/users/CBLCCBR/kopardevn_tmp/issue57_testing_2/results/GI1_T/circExplorer/GI1_T.bam|wc -l)
if [ "$nreads" == "0" ];then
    echo "SOMETHING WENT WRONG ...filtered BAM is empty!!"
    exit
fi
samtools index -@4 $filtered_bam


###################################################################################################
# convert BAM to read ends BED format
###################################################################################################

printtime $SCRIPT_NAME $start0 $start "converting filtered bam to bed format"
start=$(date +%s.%N)

bedtools bamtobed -split -i $filtered_bam > ${tmpdir}/${sample_name}.bed

python ${SCRIPT_DIR}/_process_bamtobed.py \
    --inbed ${tmpdir}/${sample_name}.bed \
    --outbed ${tmpdir}/${sample_name}.readends.bed \
    --linear ${tmpdir}/${sample_name}.linear.readids.gz \
    --spliced ${tmpdir}/${sample_name}.spliced.readids.gz \

bedSort  ${tmpdir}/${sample_name}.readends.bed  ${tmpdir}/${sample_name}.readends.bed

###################################################################################################
# get BSJ ends BED
###################################################################################################

printtime $SCRIPT_NAME $start0 $start "get bed BSJ ends from BSJ bed file"
start=$(date +%s.%N)

zcat $bsjbedgz |awk -F"\t" -v OFS="\t" '{print $1, $2, $2, $1"##"$2"##"$3"##"$6, $5, $6}' - > ${tmpdir}/BSJ.ends.bed
zcat $bsjbedgz |awk -F"\t" -v OFS="\t" '{print $1, $3, $3, $1"##"$2"##"$3"##"$6, $5, $6}' - >> ${tmpdir}/BSJ.ends.bed
bedSort ${tmpdir}/BSJ.ends.bed ${tmpdir}/BSJ.ends.bed


###################################################################################################
# finding max readlength
###################################################################################################

printtime $SCRIPT_NAME $start0 $start "finding max readlength"
start=$(date +%s.%N)

python3 ${SCRIPT_DIR}/bam_get_max_readlen.py -i $filtered_bam > ${tmpdir}/${sample_name}.maxrl
maxrl=$(cat ${tmpdir}/${sample_name}.maxrl)
two_maxrl=$((maxrl*2))


###################################################################################################
# find closest known BSJs with primary alignments and then keep alignments with some alignment within the "inclusion zone"
# this is run in parallel
###################################################################################################

printtime $SCRIPT_NAME $start0 $start "creating rid2jidgzip"
start=$(date +%s.%N)

shuf ${tmpdir}/${sample_name}.readends.bed | \
    split --lines=10000 --suffix-length=5 --numeric-suffixes=1 - ${tmpdir}/${sample_name}.readends.part.
for f in `ls ${tmpdir}/${sample_name}.readends.part.?????`;do
echo """
bedSort $f $f && \
bedtools closest -nonamecheck \
    -a $f \
    -b ${tmpdir}/BSJ.ends.bed  -d | \
    awk -F\"\\t\" -v OFS=\"\\t\" -v limit=$two_maxrl '{if (\$NF > 0 && \$NF <= limit) {print \$4,\$10}}' | \
    sort --buffer-size=2G --unique --parallel=2 > ${f}.rid2jid.closest.txt
"""
done > ${tmpdir}/para.tmp
grep -v "^$" ${tmpdir}/para.tmp > ${tmpdir}/para
parallel -j $threads < ${tmpdir}/para

#cat ${tmpdir}/${sample_name}.readends.part.?????.rid2jid.closest.txt > ${rid2jidgzip%.*}.tmp
for f in $(ls ${tmpdir}/${sample_name}.readends.part.?????.rid2jid.closest.txt);do
       cat $f
done > ${rid2jidgzip%.*}.tmp
sort --buffer-size=20G --parallel=$threads --unique ${rid2jidgzip%.*}.tmp > ${rid2jidgzip%.*}

pigz -p4 -f ${rid2jidgzip%.*} && rm -f ${rid2jidgzip%.*}.tmp

###################################################################################################
# filter linear and spliced readids to those near BSJs only
###################################################################################################

printtime $SCRIPT_NAME $start0 $start "filter linear and spliced readids to those near BSJs only; creating counts table"
start=$(date +%s.%N)

python3 ${SCRIPT_DIR}/_filter_linear_spliced_readids_w_rid2jid.py \
    --linearin ${tmpdir}/${sample_name}.linear.readids.gz \
    --splicedin ${tmpdir}/${sample_name}.spliced.readids.gz \
    --rid2jid ${rid2jidgzip} \
    --linearout ${linearrids} \
    --splicedout ${splicedrids} \
    --jidcounts ${jidcounts}


# rm -rf ${tmpdir}/${sample_name}.readends.part.*
# the above statement leads to /usr/bin/rm: Argument list too long
# hence,
find ${tmpdir} -maxdepth 1 -name "${sample_name}.readends.part.*" -print0 | xargs -0 rm -f

###################################################################################################
# create BAMS from readids
###################################################################################################

printtime $SCRIPT_NAME $start0 $start "creating linear and spliced BAMs for near BSJ reads and all reads"
start=$(date +%s.%N)

if [ -f ${tmpdir}/para2 ];then rm -f ${tmpdir}/para2;fi

linearbam_bn=$(basename $linearbam)
splicedbam_bn=$(basename $splicedbam)
linearbam_all_bn=$(basename $linearbam_all)
splicedbam_all_bn=$(basename $splicedbam_all)

echo "python3 ${SCRIPT_DIR}/filter_bam_by_readids.py --inputBAM $filtered_bam --outputBAM ${tmpdir}/${linearbam_bn} --readids $linearrids" >> ${tmpdir}/para2
echo "python3 ${SCRIPT_DIR}/filter_bam_by_readids.py --inputBAM $filtered_bam --outputBAM ${tmpdir}/${splicedbam_bn} --readids $splicedrids" >> ${tmpdir}/para2
echo "python3 ${SCRIPT_DIR}/filter_bam_by_readids.py --inputBAM $filtered_bam --outputBAM ${tmpdir}/${linearbam_all_bn} --readids ${tmpdir}/${sample_name}.linear.readids.gz" >> ${tmpdir}/para2
echo "python3 ${SCRIPT_DIR}/filter_bam_by_readids.py --inputBAM $filtered_bam --outputBAM ${tmpdir}/${splicedbam_all_bn} --readids ${tmpdir}/${sample_name}.spliced.readids.gz" >> ${tmpdir}/para2

parallel -j 4 < ${tmpdir}/para2

printtime $SCRIPT_NAME $start0 $start "sorting linear and spliced BAMs for near BSJ reads and all reads"
start=$(date +%s.%N)

if [ -d "${tmpdir}/sorttmp" ];then rm -rf ${tmpdir}/sorttmp;fi && mkdir -p ${tmpdir}/sorttmp
samtools sort -l 9 -T ${tmpdir}/sorttmp --write-index -@${threads} -O BAM -o ${linearbam} ${tmpdir}/${linearbam_bn}

if [ -d "${tmpdir}/sorttmp" ];then rm -rf ${tmpdir}/sorttmp;fi && mkdir -p ${tmpdir}/sorttmp
samtools sort -l 9 -T ${tmpdir}/sorttmp --write-index -@${threads} -O BAM -o ${splicedbam} ${tmpdir}/${splicedbam_bn}

if [ -d "${tmpdir}/sorttmp" ];then rm -rf ${tmpdir}/sorttmp;fi && mkdir -p ${tmpdir}/sorttmp
samtools sort -l 9 -T ${tmpdir}/sorttmp --write-index -@${threads} -O BAM -o ${linearbam_all} ${tmpdir}/${linearbam_all_bn}

if [ -d "${tmpdir}/sorttmp" ];then rm -rf ${tmpdir}/sorttmp;fi && mkdir -p ${tmpdir}/sorttmp
samtools sort -l 9 -T ${tmpdir}/sorttmp --write-index -@${threads} -O BAM -o ${splicedbam_all} ${tmpdir}/${splicedbam_all_bn}

###################################################################################################
# split BAMs by regions
###################################################################################################


if [ -f ${tmpdir}/para3 ];then rm -f ${tmpdir}/para3;fi

echo "python3 ${SCRIPT_DIR}/bam_split_by_regions.py --inbam $linearbam --sample_name $sample_name --regions $regions --prefix linear_BSJ --outdir $outdir --host $host --additives $additives --viruses $viruses" >> ${tmpdir}/para3
echo "python3 ${SCRIPT_DIR}/bam_split_by_regions.py --inbam $splicedbam --sample_name $sample_name --regions $regions --prefix spliced_BSJ --outdir $outdir --host $host --additives $additives --viruses $viruses" >> ${tmpdir}/para3
echo "python3 ${SCRIPT_DIR}/bam_split_by_regions.py --inbam $linearbam_all --sample_name $sample_name --regions $regions --prefix linear --outdir $outdir --host $host --additives $additives --viruses $viruses" >> ${tmpdir}/para3
echo "python3 ${SCRIPT_DIR}/bam_split_by_regions.py --inbam $splicedbam_all --sample_name $sample_name --regions $regions --prefix spliced --outdir $outdir --host $host --additives $additives --viruses $viruses" >> ${tmpdir}/para3

parallel -j 4 < ${tmpdir}/para3

# fi
# two_maxrl="302"

###################################################################################################
# convert BAMs to Bigwigs
###################################################################################################

printtime $SCRIPT_NAME $start0 $start "creating linear and spliced BIGWIGs"
start=$(date +%s.%N)

for folder in $(ls ${tmpdir}/b2b*);do rm -rf $folder;done

if [ -f ${tmpdir}/para3 ];then rm -f ${tmpdir}/para4;fi

echo "mkdir -p ${tmpdir}/b2b_1 && bash ${SCRIPT_DIR}/bam_to_bigwig.sh $linearbam ${tmpdir}/b2b_1" >> ${tmpdir}/para4
echo "mkdir -p ${tmpdir}/b2b_2 && bash ${SCRIPT_DIR}/bam_to_bigwig.sh $splicedbam ${tmpdir}/b2b_2" >> ${tmpdir}/para4
echo "mkdir -p ${tmpdir}/b2b_3 && bash ${SCRIPT_DIR}/bam_to_bigwig.sh $linearbam_all ${tmpdir}/b2b_3" >> ${tmpdir}/para4
echo "mkdir -p ${tmpdir}/b2b_4 && bash ${SCRIPT_DIR}/bam_to_bigwig.sh $splicedbam_all ${tmpdir}/b2b_4" >> ${tmpdir}/para4
count=4
for prefix in "linear_BSJ" "spliced_BSJ" "linear" "spliced";do
for b in $(echo "$host $viruses"|tr ',' ' ');do
    count=$((count+1))
    bam="${outdir}/${sample_name}.${prefix}.${b}.bam"
    echo "mkdir -p ${tmpdir}/b2b_${count} && bash ${SCRIPT_DIR}/bam_to_bigwig.sh $bam ${tmpdir}/b2b_${count}" >> ${tmpdir}/para4
done
done

parallel -j 12 < ${tmpdir}/para4

#rm -rf ${tmpdir}/*
printtime $SCRIPT_NAME $start0 $start "Done!"
