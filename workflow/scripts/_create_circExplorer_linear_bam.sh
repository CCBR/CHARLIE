#!/usr/bin/env bash
# module load parallel
# module load python/3.7
# module load bedtools
# module load ucsc
# module load samtools

threads=56
alloccpu=$(sacct -j $SLURM_JOB_ID --format "AllocCPUS"|tail -n1|awk '{print $1}')
if [ "$alloccpu" -lt "$threads" ];then
	threads=$alloccpu
fi

start0=$(date +%s.%N)

# HERE BE CODE

#end=$(date +%s.%N)    
#runtime=$(python -c "print(${end} - ${start})")

#echo "Runtime was $runtime"

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

if [ "$#" != 16 ];then
    echo "16 arguments expected!"
    echo "#1 --> path to non-chimeric BAM file"
    echo "#2 --> sample name"
    echo "#3 --> PE or SE"
    echo "#4 --> known BSJs in bed.gz format"
    echo "#5 --> tmpdir"
    echo "#6 --> gzipped outputfilename eg.${sample_name}.rid2jid.tsv.gz"
    echo "#7 --> output filtered sample BAM"
	echo "#8 --> gzip-ed list of linear BSJ readids"
	echo "#9 --> gzip-ed list of linear spliced BSJ readids"
	echo "#10 --> jid counts (linear and linear-spliced) per jid or BSJ"
	echo "#11 --> linear BSJ readids in BAM"
	echo "#12 --> linear-spliced BSJ readids in BAM"
	echo "#13 --> .regions file eg. ref/ref.fa.regions"
	echo "#14 --> host list comma separated .. no spaces"
	echo "#15 --> additives list comma separated .. no spaces"
	echo "#16 --> viruses list comma separated .. no spaces"
    exit
fi
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
non_chimeric_star2p_bam=$1
sample_name=$2
peorse=$3
bsjbedgz=$4
tmpdir=$5
#OUTPUTs
rid2jidgzip=$6
filtered_bam=$7
linearrids=$8
splicedrids=$9
jidcounts=${10}
linearbam=${11}
splicedbam=${12}
regions=${13}
host=${14}
additives=${15}
viruses=${16}
# ${tmpdir}/${sample_name}.rid2jid.tsv.gz is created!

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
SCRIPT_NAME=$( basename -- "${BASH_SOURCE[0]}" )

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

# fi

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

# two_maxrl="302"

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

cat ${tmpdir}/${sample_name}.readends.part.?????.rid2jid.closest.txt > ${rid2jidgzip%.*}.tmp
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

printtime $SCRIPT_NAME $start0 $start "creating linear and spliced BAMs"
start=$(date +%s.%N)

linearbam_bn=$(basename $linearbam)
python3 ${SCRIPT_DIR}/filter_bam_by_readids.py \
	--inputBAM $filtered_bam \
	--outputBAM ${tmpdir}/${linearbam_bn} \
	--readids $linearrids

if [ -d "${tmpdir}/sorttmp" ];then rm -rf ${tmpdir}/sorttmp;fi
mkdir -p ${tmpdir}/sorttmp
samtools sort -l 9 -T ${tmpdir}/sorttmp --write-index -@${threads} -O BAM -o ${linearbam} ${tmpdir}/${linearbam_bn}

splicedbam_bn=$(basename $splicedbam)
python3 ${SCRIPT_DIR}/filter_bam_by_readids.py \
	--inputBAM $filtered_bam \
	--outputBAM ${tmpdir}/${splicedbam_bn} \
	--readids $splicedrids

if [ -d "${tmpdir}/sorttmp" ];then rm -rf ${tmpdir}/sorttmp;fi
mkdir -p ${tmpdir}/sorttmp
samtools sort -l 9 -T ${tmpdir}/sorttmp --write-index -@${threads} -O BAM -o ${splicedbam} ${tmpdir}/${splicedbam_bn}

###################################################################################################
# split BAMs by regions
###################################################################################################

outdir=$(dirname $linearbam)
python3 ${SCRIPT_DIR}/bam_split_by_regions.py \
	--inbam $linearbam \
	--sample_name $sample_name \
	--regions $regions \
	--prefix linear \
	--outdir $outdir \
	--host $host \
	--additives $additives \
	--viruses $viruses

python3 ${SCRIPT_DIR}/bam_split_by_regions.py \
	--inbam $splicedbam \
	--sample_name $sample_name \
	--regions $regions \
	--prefix linear_spliced \
	--outdir $outdir \
	--host $host \
	--additives $additives \
	--viruses $viruses

###################################################################################################
# convert BAMs to Bigwigs
###################################################################################################

printtime $SCRIPT_NAME $start0 $start "creating linear and spliced BIGWIGs"
start=$(date +%s.%N)

bash ${SCRIPT_DIR}/bam_to_bigwig.sh $linearbam $tmpdir
bash ${SCRIPT_DIR}/bam_to_bigwig.sh $splicedbam $tmpdir
for prefix in "linear" "linear_spliced";do
for b in $(echo "$host $viruses"|tr ',' ' ');do
	bam="${outdir}/${sample_name}.${prefix}.${b}.BSJ.bam"
	bash ${SCRIPT_DIR}/bam_to_bigwig.sh $bam $tmpdir
done
done

rm -rf ${tmpdir}/*
printtime $SCRIPT_NAME $start0 $start "Done!"
