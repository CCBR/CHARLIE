
rule create_circExplorer_BSJ_bam:
    input:
        countstable=rules.circExplorer.output.counts_table,
        chimericbam=rules.star2p.output.chimeric_bam
    output:
        BSJbam=join(WORKDIR,"results","{sample}","circExplorer","{sample}.BSJ.bam"),
        plusBSJbam=join(WORKDIR,"results","{sample}","circExplorer","{sample}.BSJ.plus.bam"),
        minusBSJbam=join(WORKDIR,"results","{sample}","circExplorer","{sample}.BSJ.minus.bam"),
        BSJbed=join(WORKDIR,"results","{sample}","circExplorer","{sample}.BSJ.bed.gz"),
        BSJfoundcounts=join(WORKDIR,"results","{sample}","circExplorer","{sample}.BSJ.foundcounts.tsv"),
        plusBSJbw=join(WORKDIR,"results","{sample}","circExplorer","{sample}.BSJ.plus.bw"),
        minusBSJbw=join(WORKDIR,"results","{sample}","circExplorer","{sample}.BSJ.minus.bw"),
    params:
        sample="{sample}",
        refregions=REF_REGIONS,
        reffa=REF_FA,
        host=HOST,
        additives=ADDITIVES,
        viruses=VIRUSES,
        peorse=get_peorse,
        memG=getmemG("create_circExplorer_BSJ_bam"),
        scriptpe=join(SCRIPTS_DIR,"_create_circExplorer_BSJ_bam_pe.py"),
        scriptse=join(SCRIPTS_DIR,"_create_circExplorer_BSJ_bam_se.py"),
        flankscript=join(SCRIPTS_DIR,"_append_splice_site_flanks_to_BSJs.py"),
        bam2bwscript=join(SCRIPTS_DIR,"bam_to_bigwig.sh"),
        randomstr=str(uuid.uuid4())
    envmodules: TOOLS["python37"]["version"],TOOLS["samtools"]["version"],TOOLS["bedtools"]["version"],TOOLS["ucsc"]["version"]
    threads : getthreads("create_circExplorer_BSJ_bam")
    shell:"""
set -exo pipefail
if [ -d /lscratch/${{SLURM_JOB_ID}} ];then
    TMPDIR="/lscratch/${{SLURM_JOB_ID}}/{params.randomstr}"
else
    TMPDIR="/dev/shm/{params.randomstr}"
fi
if [ ! -d $TMPDIR ];then mkdir -p $TMPDIR;fi

outdir=$(dirname {output.BSJbam})
BSJbedbn=$(basename {output.BSJbed})

if [ "{params.peorse}" == "PE" ];then

python3 {params.scriptpe} \\
	--inbam {input.chimericbam} \\
	--sample_counts_table {input.countstable} \\
	--plusbam ${{TMPDIR}}/{params.sample}.BSJ.plus.bam \\
	--minusbam ${{TMPDIR}}/{params.sample}.BSJ.minus.bam \\
	--outbam ${{TMPDIR}}/{params.sample}.BSJ.bam \\
	--bed ${{TMPDIR}}/${{BSJbedbn}} \\
	--sample_name {params.sample} \\
	--junctionsfound {output.BSJfoundcounts} \\
	--regions {params.refregions} \\
	--host {params.host} \\
	--additives {params.additives} \\
	--viruses {params.viruses} \\
	--outputhostbams --outputvirusbams --outdir $TMPDIR

else

python3 {params.scriptse} \\
	--inbam {input.chimericbam} \\
	--sample_counts_table {input.countstable} \\
	--plusbam ${{TMPDIR}}/{params.sample}.BSJ.plus.bam \\
	--minusbam ${{TMPDIR}}/{params.sample}.BSJ.minus.bam \\
	--outbam ${{TMPDIR}}/{params.sample}.BSJ.bam \\
	--bed ${{TMPDIR}}/${{BSJbedbn}} \\
	--sample_name {params.sample} \\
	--junctionsfound {output.BSJfoundcounts} \\
	--regions {params.refregions} \\
	--host {params.host} \\
	--additives {params.additives} \\
	--viruses {params.viruses} \\
	--outputhostbams --outputvirusbams --outdir $TMPDIR

fi

samtools sort -l 9 -T $TMPDIR --write-index -@{threads} -O BAM -o {output.plusBSJbam} ${{TMPDIR}}/{params.sample}.BSJ.plus.bam 
samtools sort -l 9 -T $TMPDIR --write-index -@{threads} -O BAM -o {output.minusBSJbam} ${{TMPDIR}}/{params.sample}.BSJ.minus.bam 
samtools sort -l 9 -T $TMPDIR --write-index -@{threads} -O BAM -o {output.BSJbam} ${{TMPDIR}}/{params.sample}.BSJ.bam

# for b in {output.plusBSJbam} {output.minusBSJbam} {output.BSJbam} 
for b in {output.plusBSJbam} {output.minusBSJbam}
do
    bash {params.bam2bwscript} $b
done

for i in $(echo {params.host}|tr ',' ' ');do 
    samtools sort -l 9 -T $TMPDIR --write-index -@{threads} -O BAM -o ${{outdir}}/{params.sample}.${{i}}.BSJ.bam ${{TMPDIR}}/{params.sample}.${{i}}.BSJ.bam
    # bash {params.bam2bwscript} ${{outdir}}/{params.sample}.${{i}}.BSJ.bam
done
for i in $(echo {params.viruses}|tr ',' ' ');do 
    samtools sort -l 9 -T $TMPDIR --write-index -@{threads} -O BAM -o ${{outdir}}/{params.sample}.${{i}}.BSJ.bam ${{TMPDIR}}/{params.sample}.${{i}}.BSJ.bam
    # bash {params.bam2bwscript} ${{outdir}}/{params.sample}.${{i}}.BSJ.bam
done

python3 {params.flankscript} --reffa {params.reffa} --inbsjbedgz ${{TMPDIR}}/${{BSJbedbn}} --outbsjbedgz {output.BSJbed}



rm -rf $TMPDIR
"""

rule create_circExplorer_linear_spliced_bams:
    input:
        junction=rules.star2p.output.junction, # Chimeric junctions file
        bam=rules.star2p.output.bam,
        nonchimericbam=rules.star2p.output.non_chimeric_bam,
        bsjbedgz=rules.create_circExplorer_BSJ_bam.output.BSJbed,
        countstable=rules.circExplorer.output.counts_table,
    output:
        rid2jid=join(WORKDIR,"results","{sample}","circExplorer","{sample}.rid2jid.tsv.gz"),
        filtered_bam=join(WORKDIR,"results","{sample}","circExplorer","{sample}.bam"),
        linear_bam=join(WORKDIR,"results","{sample}","circExplorer","{sample}.linear.bam"),                         # linear reads in BSJ inclusion zone
        spliced_bam=join(WORKDIR,"results","{sample}","circExplorer","{sample}.spliced.bam"),                       # spliced-only alignments in the sample
        linear_spliced_bam=join(WORKDIR,"results","{sample}","circExplorer","{sample}.linear.spliced.bam"),         # spliced-only linear alignments in BSJ inclusion zone
        linear_foundcounts=join(WORKDIR,"results","{sample}","circExplorer","{sample}.linear.foundcounts.tsv"),
# strand specific bams
        linear_plus_bam=join(WORKDIR,"results","{sample}","circExplorer","{sample}.linear.plus.bam"),
        linear_minus_bam=join(WORKDIR,"results","{sample}","circExplorer","{sample}.linear.minus.bam"),
        linear_spliced_plus_bam=join(WORKDIR,"results","{sample}","circExplorer","{sample}.linear.spliced.plus.bam"), 
        linear_spliced_minus_bam=join(WORKDIR,"results","{sample}","circExplorer","{sample}.linear.spliced.minus.bam"),
# bigwigs
        # linear_bw=join(WORKDIR,"results","{sample}","circExplorer","{sample}.linear.bw"),
        linear_plus_bw=join(WORKDIR,"results","{sample}","circExplorer","{sample}.linear.plus.bw"),
        linear_minus_bw=join(WORKDIR,"results","{sample}","circExplorer","{sample}.linear.minus.bw"),

        # spliced_bw=join(WORKDIR,"results","{sample}","circExplorer","{sample}.spliced.bw"),

        # linear_spliced_bw=join(WORKDIR,"results","{sample}","circExplorer","{sample}.linear.spliced.bw"),
        linear_spliced_plus_bw=join(WORKDIR,"results","{sample}","circExplorer","{sample}.linear.spliced.plus.bw"), 
        linear_spliced_minus_bw=join(WORKDIR,"results","{sample}","circExplorer","{sample}.linear.spliced.minus.bw"),

    params:
        sample="{sample}",
        memG=getmemG("create_circExplorer_linear_spliced_bams"),
        refregions=REF_REGIONS,
        refregions_host=REF_REGIONS_HOST,
        refregions_viruses=REF_REGIONS_VIRUSES,
        host=HOST,
        additives=ADDITIVES,
        viruses=VIRUSES,
        peorse=get_peorse,
        bashscript=join(SCRIPTS_DIR,"_create_circExplorer_linear_bam_preprocess.sh"),
        pythonscript=join(SCRIPTS_DIR,"_extract_circExplorer_linear_reads.py"),
        bam2bwscript=join(SCRIPTS_DIR,"bam_to_bigwig.sh"),
        outdir=join(WORKDIR,"results","{sample}","circExplorer"),
        randomstr=str(uuid.uuid4())
    envmodules: TOOLS["python37"]["version"],TOOLS["samtools"]["version"],TOOLS["bedtools"]["version"],TOOLS["ucsc"]["version"]
    threads : getthreads("create_circExplorer_linear_spliced_bams")
    shell:"""
set -exo pipefail
if [ -d /lscratch/${{SLURM_JOB_ID}} ];then
    TMPDIR="/lscratch/${{SLURM_JOB_ID}}"
else
    TMPDIR="/dev/shm/{params.randomstr}"
    mkdir -p $TMPDIR
fi
cd {params.outdir}

# get filtered bam (remove secondary/supplementary/etc.) and the rid2jid lookup files

bash {params.bashscript} {input.nonchimericbam} {params.sample} {params.peorse} {input.bsjbedgz} $TMPDIR {output.rid2jid} {output.filtered_bam}
samtools index -@{threads} {output.filtered_bam}

# extract linear, spliced, linear_sliced reads into individual BAMs using the filtered bam and rid2jid from above

linear_bam_bn=$(basename {output.linear_bam})
linear_plus_bam_bn=$(basename {output.linear_plus_bam})
linear_minus_bam_bn=$(basename {output.linear_minus_bam})
spliced_bam_bn=$(basename {output.spliced_bam})
linear_spliced_bam_bn=$(basename {output.linear_spliced_bam})
linear_spliced_plus_bam_bn=$(basename {output.linear_spliced_plus_bam})
linear_spliced_minus_bam_bn=$(basename {output.linear_spliced_minus_bam})

python3 {params.pythonscript} \\
    --inbam {output.filtered_bam} \\
    --rid2jid {output.rid2jid} \\
    --sample_counts_table {input.countstable} \\
    --sample_name {params.sample} \\
    --regions {params.refregions} \\
    --host {params.host} \\
    --additives {params.additives} \\
    --viruses {params.viruses} \\
    --outbam ${{TMPDIR}}/$linear_bam_bn \\
    --outplusbam ${{TMPDIR}}/$linear_plus_bam_bn \\
    --outminusbam ${{TMPDIR}}/$linear_minus_bam_bn \\
    --splicedbam ${{TMPDIR}}/$spliced_bam_bn \\
    --splicedbsjbam ${{TMPDIR}}/$linear_spliced_bam_bn \\
    --splicedbsjplusbam ${{TMPDIR}}/$linear_spliced_plus_bam_bn \\
    --splicedbsjminusbam ${{TMPDIR}}/$linear_spliced_minus_bam_bn \\
    --outputhostbams --outputvirusbams \\
    --outdir {params.outdir} \\
    --countsfound {output.linear_foundcounts}

samtools sort -l 9 -T $TMPDIR --write-index -@{threads} -O BAM -o {output.linear_bam} ${{TMPDIR}}/$linear_bam_bn
samtools sort -l 9 -T $TMPDIR --write-index -@{threads} -O BAM -o {output.linear_plus_bam} ${{TMPDIR}}/$linear_plus_bam_bn
samtools sort -l 9 -T $TMPDIR --write-index -@{threads} -O BAM -o {output.linear_minus_bam} ${{TMPDIR}}/$linear_minus_bam_bn
samtools sort -l 9 -T $TMPDIR --write-index -@{threads} -O BAM -o {output.spliced_bam} ${{TMPDIR}}/$spliced_bam_bn
samtools sort -l 9 -T $TMPDIR --write-index -@{threads} -O BAM -o {output.linear_spliced_bam} ${{TMPDIR}}/$linear_spliced_bam_bn
samtools sort -l 9 -T $TMPDIR --write-index -@{threads} -O BAM -o {output.linear_spliced_plus_bam} ${{TMPDIR}}/$linear_spliced_plus_bam_bn
samtools sort -l 9 -T $TMPDIR --write-index -@{threads} -O BAM -o {output.linear_spliced_minus_bam} ${{TMPDIR}}/$linear_spliced_minus_bam_bn

# for b in {output.linear_bam} {output.linear_plus_bam} {output.linear_minus_bam} {output.spliced_bam} {output.linear_spliced_bam} {output.linear_spliced_plus_bam} {output.linear_spliced_minus_bam}
for b in {output.linear_plus_bam} {output.linear_minus_bam} {output.linear_spliced_plus_bam} {output.linear_spliced_minus_bam}
do
    bash {params.bam2bwscript} $b
done

for regions in {params.refregions_host} {params.refregions_viruses};do
while read a b;do 
    samtools view -@{threads} --write-index -O BAM -o {params.sample}.${{a}}.linear.bam {output.linear_bam} $b 
    samtools view -@{threads} --write-index -O BAM -o {params.sample}.${{a}}.spliced.bam {output.spliced_bam} $b
    samtools view -@{threads} --write-index -O BAM -o {params.sample}.${{a}}.linear_spliced.bam {output.linear_spliced_bam} $b
done < $regions
done

cd $TMPDIR && rm -f *
"""

# rule create_circExplorer_merged_found_counts_table:
# This rule  creates an output file with the following columns:
# chrom
# start
# end
# strand
# expected_BSJ_reads
# found_BSJ_reads
# linear_BSJ_reads_same_strand
# linear_spliced_BSJ_reads_same_strand
# linear_BSJ_reads_opposite_strand
# linear_spliced_BSJ_reads_opposite_strand
localrules: create_circExplorer_merged_found_counts_table
rule create_circExplorer_merged_found_counts_table:
    input:
        bsj_found_counts=rules.create_circExplorer_BSJ_bam.output.BSJfoundcounts,
        linear_found_counts=rules.create_circExplorer_spliced_bams.output.linear_foundcounts
    output:
        count_counts_table=join(WORKDIR,"results","{sample}","circExplorer","{sample}.readcounts.tsv")
    params:
        sample="{sample}",
        pythonscript=join(SCRIPTS_DIR,"_merge_circExplorer_found_counts.py"),
        outdir=join(WORKDIR,"results","{sample}","circExplorer"),
        randomstr=str(uuid.uuid4())    
    shell:"""
set -exo pipefail
if [ -d /lscratch/${{SLURM_JOB_ID}} ];then
    TMPDIR="/lscratch/${{SLURM_JOB_ID}}"
else
    TMPDIR="/dev/shm/{params.randomstr}"
    mkdir -p $TMPDIR
fi
python3 {params.pythonscript} \\
    -b {input.bsj_found_counts} \\
    -l {input.linear_found_counts} \\
    -o {output.count_counts_table}
"""


# rule alignment_stats:
#     input:
#         star2pbam=rules.star2p.output.bam,
#         splicedreadsbam=rules.create_spliced_reads_bam.output.bam,
#         BSJbam=rules.create_circExplorer_BSJ_bam.output.BSJbam,
#     output:
#         alignmentstats=join(WORKDIR,"results","{sample}","STAR2p","alignmentstats.txt")
#     params:
#         regions=REF_REGIONS,
#         bam2nreads_script=join(SCRIPTS_DIR,"  h"),
#         outdir=join(WORKDIR,"results","{sample}","STAR2p")
#     threads: getthreads("alignment_stats")
#     envmodules: TOOLS["samtools"]["version"]
#     shell:"""
# set -exo pipefail
# while read a b;do echo $a;done < {params.regions} > {params.outdir}/tmp1
# while read a b;do bash {params.bam2nreads_script} {input.star2pbam} "$b";done < {params.regions} > {params.outdir}/tmp2
# while read a b;do bash {params.bam2nreads_script} {input.splicedreadsbam} "$b";done < {params.regions} > {params.outdir}/tmp3
# while read a b;do bash {params.bam2nreads_script} {input.BSJbam} "$b";done < {params.regions} > {params.outdir}/tmp4
# echo -ne "#region\taligned_fragments\tspliced_fragments\tBSJ_fragments\n" > {output.alignmentstats}
# paste {params.outdir}/tmp1 {params.outdir}/tmp2 {params.outdir}/tmp3 {params.outdir}/tmp4 >> {output.alignmentstats}
# rm -f {params.outdir}/tmp1 {params.outdir}/tmp2 {params.outdir}/tmp3 {params.outdir}/tmp4 
# """

# localrules: venn
# rule venn:
#     input:
#         circexplorerout=rules.circExplorer.output.annotations,
#         ciriout=rules.ciri.output.cirioutfiltered
#     output:
#         png=join(WORKDIR,"results","{sample}","{sample}.venn_mqc.png"),
#         cirionly=join(WORKDIR,"results","{sample}","{sample}.cirionly.lst"),
#         circexploreronly=join(WORKDIR,"results","{sample}","{sample}.circexploreronly.lst"),
#         common=join(WORKDIR,"results","{sample}","{sample}.common.lst")
#     params:
#         script1=join(SCRIPTS_DIR,"venn.R"),
#         sample="{sample}",
#         outdir=join(WORKDIR,"results","{sample}")
#     container: "docker://nciccbr/ccbr_venn:latest"
#     threads: getthreads("venn")
#     shell:"""
# set -exo pipefail
# cut -f1 {input.ciriout}|grep -v circRNA_ID > /dev/shm/{params.sample}.ciri.lst
# cut -f1-3 {input.circexplorerout}|awk -F"\\t" '{{print $1\":\"$2+1\"|\"$3}}' > /dev/shm/{params.sample}.circExplorer.lst
# if [[ "$(cat /dev/shm/{params.sample}.ciri.lst | wc -l)" != "0" ]];then
#     if [[ "$(cat /dev/shm/{params.sample}.circExplorer.lst | wc -l)" != "0" ]];then
#         2set_venn.R \\
#             -l /dev/shm/{params.sample}.ciri.lst \\
#             -r /dev/shm/{params.sample}.circExplorer.lst  \\
#             -p {output.png} \\
#             -m {output.cirionly} \\
#             -s {output.circexploreronly} \\
#             -c1 "CIRI2" \\
#             -c2 "CircExplorer2" \\
#             -c {output.common} \\
#             -t {params.sample}
#     else
#         for o in {output}
#         do
#             touch $o
#         done
#     fi
# fi
# rm -f /dev/shm{params.sample}*
# """


# rule add_strand_to_circExplorer:
# 	input:
# 		backsplicedjunctions=rules.circExplorer.output.backsplicedjunctions,
# 		bed=rules.split_BAM_create_BW.output.bed
# 	output:
# 		backsplicedjunctions=join(WORKDIR,"results","{sample}","circExplorer","{sample}.back_spliced_junction.bed"),
# 	params:
# 		sample="{sample}",
# 		workdir=WORKDIR,
# 		pythonscript=join(SCRIPTS_DIR,"copy_strand_info.py"),
# 	threads: 2
# 	# envmodules:
# 	shell:"""
# 	cp {input.backsplicedjunctions} /dev/shm
# 	bsj_basename="$(basename {input.backsplicedjunctions})"
# 	python {params.pythonscript} --from {input.bed} --to /dev/shm/${{bsj_basename}} --output {output.backsplicedjunctions}
# """

# rule recall_valid_BSJ_split_BAM_by_strand_create_BW:
#     input:
#         bam=rules.create_circExplorer_BSJ_bam.output.BSJbam,
#     output:
#         bam=join(WORKDIR,"results","{sample}","STAR2p","BSJs","{sample}.BSJ."+HOST+".bam"),
#         plusbam=join(WORKDIR,"results","{sample}","STAR2p","BSJs","{sample}.BSJ."+HOST+".plus.bam"),
#         minusbam=join(WORKDIR,"results","{sample}","STAR2p","BSJs","{sample}.BSJ."+HOST+".minus.bam"),
#         bed=join(WORKDIR,"results","{sample}","STAR2p","BSJs","{sample}.BSJ."+HOST+".bed"),
#         bw=join(WORKDIR,"results","{sample}","STAR2p","BSJs","{sample}.BSJ."+HOST+".bw"),
#         plusbw=join(WORKDIR,"results","{sample}","STAR2p","BSJs","{sample}.BSJ."+HOST+".plus.bw"),
#         minusbw=join(WORKDIR,"results","{sample}","STAR2p","BSJs","{sample}.BSJ."+HOST+".minus.bw"),
#         validBSJbed=join(WORKDIR,"results","{sample}","customBSJs","{sample}.valid_STAR_BSJ.bed")
#     params:
#         sample="{sample}",
#         workdir=WORKDIR,
#         memG=getmemG("recall_valid_BSJ_split_BAM_by_strand_create_BW"),
#         outdir=join(WORKDIR,"results","{sample}","STAR2p","BSJs"),
#         pythonscript=join(SCRIPTS_DIR,"validate_BSJ_reads_and_split_BSJ_bam_by_strand.py"),
#         regions=REF_REGIONS,
#         randomstr=str(uuid.uuid4())
#     threads: getthreads("recall_valid_BSJ_split_BAM_by_strand_create_BW")
#     envmodules: TOOLS["samtools"]["version"],TOOLS["bedtools"]["version"],TOOLS["ucsc"]["version"],TOOLS["sambamba"]["version"],TOOLS["python37"]["version"]
#     shell:"""
# set -exo pipefail
# if [ -d /lscratch/${{SLURM_JOB_ID}} ];then
#     TMPDIR="/lscratch/${{SLURM_JOB_ID}}"
# else
#     TMPDIR="/dev/shm/{params.randomstr}"
#     mkdir -p $TMPDIR
# fi
# cd {params.outdir}
# bam_basename="$(basename {input.bam})"
# while read a b;do
# bam="${{bam_basename%.*}}.${{a}}.bam"
# plusbam="${{bam_basename%.*}}.${{a}}.plus.bam"
# minusbam="${{bam_basename%.*}}.${{a}}.minus.bam"
# bed="${{bam_basename%.*}}.${{a}}.bed"
# samtools view {input.bam} $b -b > ${{TMPDIR}}/${{bam%.*}}.tmp.bam
# sambamba sort --memory-limit={params.memG} --tmpdir=${{TMPDIR}} --nthreads={threads} --out=$bam ${{TMPDIR}}/${{bam%.*}}.tmp.bam
# bw="${{bam%.*}}.bw"
# bdg="${{bam%.*}}.bdg"
# plusbw="${{bam%.*}}.plus.bw"
# plusbdg="${{bam%.*}}.plus.bdg"
# minusbw="${{bam%.*}}.minus.bw"
# minusbdg="${{bam%.*}}.minus.bdg"
# sizes="${{bam%.*}}.sizes"
# # create strand specific bams
# python {params.pythonscript} \
# -i $bam \
# -p $plusbam \
# -m $minusbam \
# -b $bed
# bedSort $bed $bed
# # create bedgraphs
# bedtools genomecov -bga -split -ibam $bam > $bdg
# bedSort $bdg $bdg
# bedtools genomecov -bga -split -ibam $plusbam > $plusbdg
# bedSort $plusbdg $plusbdg
# bedtools genomecov -bga -split -ibam $minusbam > $minusbdg
# bedSort $minusbdg $minusbdg
# # create bigwigs
# if [ "$(wc -l $bdg|awk '{{print $1}}')" != "0" ];then
#     samtools view -H $bam|grep ^@SQ|cut -f2,3|sed "s/SN://g"|sed "s/LN://g" > $sizes
#     bedGraphToBigWig $bdg $sizes $bw
#     if [ "$(wc -l $plusbdg|awk '{{print $1}}')" != "0" ];then
#         bedGraphToBigWig $plusbdg $sizes $plusbw
#     else
#         touch $plusbw
#     fi
#     if [ "$(wc -l $minusbdg|awk '{{print $1}}')" != "0" ];then
#         bedGraphToBigWig $minusbdg $sizes $minusbw
#     else
#         touch $minusbw
#     fi
# else
#     touch $bw
# fi
# rm -f $bdg $sizes

# done < {params.regions}
# cat {params.sample}.BSJ.*.bed |cut -f1-6 > ${{TMPDIR}}/{sample}.tmp.bed
# bedSort ${{TMPDIR}}/{sample}.tmp.bed ${{TMPDIR}}/{sample}.tmp.bed
# awk -F"\\t" -v OFS="\\t" '{{$4="S"NR;print}}' ${{TMPDIR}}/{sample}.tmp.bed > {output.validBSJbed}
# cd $TMPDIR && rm -f *
# """	

# rule add_novel_ciri_BSJs_to_customBSJ:
#     input:
#         ciriout=rules.ciri.output.ciriout,
#         hg38bed=rules.recall_valid_BSJ_split_BAM_by_strand_create_BW.output.bed,
#         validBSJbed=rules.recall_valid_BSJ_split_BAM_by_strand_create_BW.output.validBSJbed,
#     output:
#         novelBSJ=join(WORKDIR,"results","{sample}","customBSJs","{sample}.novel_CIRI_BSJ.bed"),
#         cirireadids=temp(join(WORKDIR,"results","{sample}","customBSJs","{sample}.ciri.readids")),
#         combinedBSJ=join(WORKDIR,"results","{sample}","customBSJs","{sample}.BSJ.bed")
#     params:
#         sample="{sample}",
#         outdir=join(WORKDIR,"results","{sample}","customBSJs"),
#         list_compare_script=join(SCRIPTS_DIR,"_compare_lists.py"),
#         randomstr=str(uuid.uuid4())
#     envmodules: TOOLS["ucsc"]["version"],TOOLS["python37"]["version"]
#     shell:"""
# set -exo pipefail
# if [ -d /lscratch/${{SLURM_JOB_ID}} ];then
#     TMPDIR="/lscratch/${{SLURM_JOB_ID}}"
# else
#     TMPDIR="/dev/shm/{params.randomstr}"
#     mkdir -p $TMPDIR
# fi
# cd {params.outdir}
# grep -v junction_reads_ID {input.ciriout}|awk -F"\\t" '{{print substr($NF,1,length($NF)-1)}}'|tr ',' '\\n'|sort|uniq > {params.sample}.ciri.readids
# stardir=$(dirname {input.hg38bed})
# cat ${{stardir}}/{params.sample}.BSJ.*.bed|awk -F"\\t" '{{print $NF}}'|tr ',' '\\n'|awk -F"##" '{{print $1}}'|sort|uniq > {params.sample}.circExplorer.readids
# python {params.list_compare_script} {params.sample}.ciri.readids {params.sample}.circExplorer.readids 1
# mv a_only.lst {params.sample}.ciri.novel.readids && rm -f a_* b_*
# head -n1 {input.ciriout} > {params.sample}.ciri.novel.out
# while read a;do grep $a {input.ciriout};done < {params.sample}.ciri.novel.readids |sort|uniq >> {params.sample}.ciri.novel.out
# grep -v circRNA_ID {params.sample}.ciri.novel.out | awk -F"\\t" -v OFS="\\t" '{{print $2,$3-1,$4,$2":"$3-1"-"$4,$5,$11}}'  > ${{TMPDIR}}/{params.sample}.novel_CIRI_BSJ.tmp.bed
# awk -F"\\t" '{{print $4}}' ${{TMPDIR}}/{params.sample}.novel_CIRI_BSJ.tmp.bed > ${{TMPDIR}}/{params.sample}.novel_CIRI_BSJ.tmp.lst
# awk -F"\\t" '{{print $1":"$2"-"$3}}' {input.validBSJbed} > ${{TMPDIR}}/{params.sample}.circExplorer.BSJ.lst
# python {params.list_compare_script} ${{TMPDIR}}/{params.sample}.novel_CIRI_BSJ.tmp.lst ${{TMPDIR}}/{params.sample}.circExplorer.BSJ.lst 1
# mv a_only.lst {params.sample}.ciri.novel.BSJ.lst && rm -f a_* b_*
# while read a;do grep $a ${{TMPDIR}}/{params.sample}.novel_CIRI_BSJ.tmp.bed;done < {params.sample}.ciri.novel.BSJ.lst > ${{TMPDIR}}/{params.sample}.ciri.novel.BSJ.bed
# bedSort ${{TMPDIR}}/{params.sample}.ciri.novel.BSJ.bed ${{TMPDIR}}/{params.sample}.ciri.novel.BSJ.bed
# awk -F"\\t" -v OFS="\\t" '{{$4="C"NR;print}}' ${{TMPDIR}}/{params.sample}.ciri.novel.BSJ.bed > {output.novelBSJ}
# cat {input.validBSJbed} {output.novelBSJ} > {output.combinedBSJ}
# bedSort {output.combinedBSJ} {output.combinedBSJ}
# rm -f {params.sample}.circExplorer.readids {params.sample}.ciri.novel.*
# cd $TMPDIR && rm -f *
# """

rule filter_ciri_bam_for_BSJs:
    input:
        bam=rules.ciri.output.ciribam,
        ciriout=rules.ciri.output.cirioutfiltered,
        readids=rules.add_novel_ciri_BSJs_to_customBSJ.output.cirireadids
    output:
        bam=join(WORKDIR,"results","{sample}","ciri","{sample}.bwa.BSJ.bam")
    params:
        sample="{sample}",
        memG=getmemG("filter_ciri_bam_for_BSJs"),
        script=join(SCRIPTS_DIR,"filter_bam_by_readids.py"),
        randomstr=str(uuid.uuid4())
    threads: getthreads("filter_ciri_bam_for_BSJs")
    envmodules: TOOLS["python37"]["version"],TOOLS["sambamba"]["version"]
    shell:"""
set -exo pipefail
if [ -d /lscratch/${{SLURM_JOB_ID}} ];then
    TMPDIR="/lscratch/${{SLURM_JOB_ID}}"
else
    TMPDIR="/dev/shm/{params.randomstr}"
    mkdir -p $TMPDIR
fi
sambamba sort --memory-limit={params.memG} --tmpdir=${{TMPDIR}} --nthreads={threads} --out=${{TMPDIR}}/{params.sample}.ciri.sorted.bam {input.bam}
python {params.script} --inputBAM ${{TMPDIR}}/{params.sample}.ciri.sorted.bam --outputBAM ${{TMPDIR}}/{params.sample}.tmp.bam --readids {input.readids}
sambamba sort --memory-limit={params.memG} --tmpdir=${{TMPDIR}} --nthreads={threads} --out={output.bam} ${{TMPDIR}}/{params.sample}.tmp.bam
cd $TMPDIR && rm -f *
"""