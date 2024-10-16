#!/usr/bin/env bash

set -e -x -o pipefail

SCRIPTDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
SCRIPTNAME=$( basename -- "${BASH_SOURCE[0]}" )

ARGPARSE_DESCRIPTION="run CircExplorer with STAR junctions file as input"
    source /opt2/argparse.bash || \
    source ${SCRIPTDIR}/../../resources/argparse.bash || \
    source <(curl -s https://raw.githubusercontent.com/CCBR/Tools/master/scripts/argparse.bash) || \
    exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--junctionfile',required=True, help='STAR 2p junction file')
parser.add_argument('--tmpdir',required=True, help='temp folder')
parser.add_argument('--outdir',required=True, help='output folder')
parser.add_argument('--samplename',required=True, help='sample name')
parser.add_argument('--genepred',required=True, help='path to genepred_w_geneid file')
parser.add_argument('--reffa',required=True, help='path to ref.fa')
parser.add_argument('--minreads',required=True, help='filter circRNAs with read support less than this number')
parser.add_argument('--hostminfilter',required=False, type=int, default=150, help='min host width of circRNA')
parser.add_argument('--hostmaxfilter',required=False, type=int, default=1000000000, help='max host width of circRNA')
parser.add_argument('--virusminfilter',required=False, type=int, default=150, help='min virus width of circRNA')
parser.add_argument('--virusmaxfilter',required=False, type=int, default=1000000000, help='max virus width of circRNA')
parser.add_argument('--regions',required=True, help='path to .regions file eg. ref/ref.fa.regions')
parser.add_argument('--host',required=True, help='comma separated list of HOST')
parser.add_argument('--additives',required=True, help='comma separated list of ADDITIVES')
parser.add_argument('--viruses',required=True, help='comma separated list of VIRUSES')
parser.add_argument('--outcount',required=True, help='output count file')
parser.add_argument('--outbsj',required=True, help='output BSJ file')
parser.add_argument('--outannotation',required=True, help='output known.txt file')
EOF


tmpstr=$(basename `mktemp`)
tmpstr=${tmpstr##*.}
TMPDIR="${TMPDIR}/${tmpstr}"
if [ ! -d $TMPDIR ];then mkdir -p $TMPDIR;fi

PARSELOG="${OUTDIR}/${SAMPLENAME}_circexplorer_parse.log"
ORIGINALBSJBED="${OUTBSJ}"
STRANDFIXEDBSJBED="${OUTDIR}/${SAMPLENAME}.back_spliced_junction.strand_fixed.bed"
KNOWNTXT="${OUTANNOTATION}"
FILTEREDKNOWNTXT="${OUTDIR}/${SAMPLENAME}.circRNA_known.filter1.txt"
LOWCONF="${OUTDIR}/low_conf_${SAMPLENAME}.circRNA_known.txt"
FILTEREDLOWCONF="${OUTDIR}/low_conf_${SAMPLENAME}.circRNA_known.filter1.txt"

cd $TMPDIR

# create junction file
grep -v junction_type ${JUNCTIONFILE} > junction

# run CircExplorer2 parse
CIRCexplorer2 parse -t STAR junction > $PARSELOG 2>&1

# copy back original back_spliced BED file
cp back_spliced_junction.bed $ORIGINALBSJBED

python -E ${SCRIPTDIR}/_circExplorer_BSJ_get_strand.py ${JUNCTIONFILE} back_spliced_junction.bed ${MINREADS} > back_spliced_junction.strand_fixed.bed

# copy back strand_fixed BSJ BED
cp back_spliced_junction.strand_fixed.bed $STRANDFIXEDBSJBED

# run CIRCexplorer2 annotate ... also generate low-conf circRNAs
CIRCexplorer2 annotate -r $GENEPRED -g $REFFA -b back_spliced_junction.bed -o circRNA_known.txt --low-confidence

# copy back known.txt and low-conf files
cp circRNA_known.txt $KNOWNTXT
cp low_conf_circRNA_known.txt $LOWCONF

# filter known.txt and low-conf files to exclude calls with low read support
cat $KNOWNTXT |tr '/' '\t'|cut -f1-3,5- |awk -v m=$MINREADS '$4>=m' > $FILTEREDKNOWNTXT
cat $LOWCONF |tr '/' '\t'|cut -f1-3,5- |awk -v m=$MINREADS '$4>=m' > $FILTEREDLOWCONF

python -E ${SCRIPTDIR}/circExplorer_get_annotated_counts_per_sample.py \
	--back_spliced_bed $STRANDFIXEDBSJBED \
	--back_spliced_min_reads $MINREADS \
	--circularRNA_known $FILTEREDKNOWNTXT \
	--low_conf $FILTEREDLOWCONF \
	--host "$HOST" \
	--additives "$ADDITIVES" \
	--viruses "$VIRUSES" \
	--virus_filter_min $VIRUSMINFILTER \
	--host_filter_min $HOSTMINFILTER \
	--host_filter_max $HOSTMAXFILTER \
	--virus_filter_max $VIRUSMAXFILTER \
	--regions $REGIONS \
	-o $OUTCOUNT
