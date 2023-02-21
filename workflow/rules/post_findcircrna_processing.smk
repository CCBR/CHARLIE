
rule create_circExplorer_BSJ_bam:
    input:
        countstable=rules.circExplorer.output.counts_table,
        chimericbam=rules.star2p.output.chimeric_bam
    output:
        BSJbam=join(WORKDIR,"results","{sample}","circExplorer","{sample}.BSJ.bam"),
        plusBSJbam=join(WORKDIR,"results","{sample}","circExplorer","{sample}.BSJ.plus.bam"),
        minusBSJbam=join(WORKDIR,"results","{sample}","circExplorer","{sample}.BSJ.minus.bam"),
        BSJbed=join(WORKDIR,"results","{sample}","circExplorer","{sample}.BSJ.bed")
    params:
        sample="{sample}",
        peorse=get_peorse,
        memG=getmemG("create_circExplorer_BSJ_bam"),
        script1=join(SCRIPTS_DIR,"create_circExplorer_BSJ_bam_pe.py"),
        script2=join(SCRIPTS_DIR,"create_circExplorer_BSJ_bam_se.py"),
        randomstr=str(uuid.uuid4())
    envmodules: TOOLS["python37"]["version"],TOOLS["samtools"]["version"]
    threads : getthreads("create_circExplorer_BSJ_bam")
    shell:"""
set -exo pipefail
if [ -d /lscratch/${{SLURM_JOB_ID}} ];then
    TMPDIR="/lscratch/${{SLURM_JOB_ID}}/{params.randomstr}"
else
    TMPDIR="/dev/shm/{params.randomstr}"
fi
if [ ! -d $TMPDIR ];then mkdir -p $TMPDIR;fi

if [ "{params.peorse}" == "PE" ];then

python3 {params.script1} \\
    --inbam {input.chimericbam} \\
    --sample_counts_table {input.countstable} \\
    -p ${{TMPDIR}}/{params.sample}.plus.bam \\
    -m ${{TMPDIR}}/{params.sample}.minus.bam \\
    -b {output.BSJbed}


else

python3 {params.script2} \\
    --inbam {input.chimericbam} \\
    --sample_counts_table {input.countstable} \\
    -p ${{TMPDIR}}/{params.sample}.plus.bam \\
    -m ${{TMPDIR}}/{params.sample}.minus.bam \\
    -b {output.BSJbed}

fi

samtools sort -l 9 -T $TMPDIR --write-index -@{threads} -O BAM -o {output.plusBSJbam} ${{TMPDIR}}/{params.sample}.plus.bam 
samtools sort -l 9 -T $TMPDIR --write-index -@{threads} -O BAM -o {output.minusBSJbam} ${{TMPDIR}}/{params.sample}.minus.bam 
samtools merge -@{threads} -u -o ${{TMPDIR}}/{params.sample}.BSJ.bam ${{TMPDIR}}/{params.sample}.plus.bam ${{TMPDIR}}/{params.sample}.minus.bam 
samtools sort -l 9 -T $TMPDIR --write-index -@{threads} -O BAM -o {output.BSJbam} ${{TMPDIR}}/{params.sample}.BSJ.bam

rm -rf $TMPDIR
"""

rule alignment_stats:
    input:
        star2pbam=rules.star2p.output.bam,
        splicedreadsbam=rules.create_spliced_reads_bam.output.bam,
        BSJbam=rules.create_circExplorer_BSJ_bam.output.BSJbam,
    output:
        alignmentstats=join(WORKDIR,"results","{sample}","STAR2p","alignmentstats.txt")
    params:
        regions=REF_REGIONS,
        bam2nreads_script=join(SCRIPTS_DIR,"bam_get_uniquely_aligned_fragments.bash"),
        outdir=join(WORKDIR,"results","{sample}","STAR2p")
    threads: getthreads("alignment_stats")
    envmodules: TOOLS["samtools"]["version"]
    shell:"""
set -exo pipefail
while read a b;do echo $a;done < {params.regions} > {params.outdir}/tmp1
while read a b;do bash {params.bam2nreads_script} {input.star2pbam} "$b";done < {params.regions} > {params.outdir}/tmp2
while read a b;do bash {params.bam2nreads_script} {input.splicedreadsbam} "$b";done < {params.regions} > {params.outdir}/tmp3
while read a b;do bash {params.bam2nreads_script} {input.BSJbam} "$b";done < {params.regions} > {params.outdir}/tmp4
echo -ne "#region\taligned_fragments\tspliced_fragments\tBSJ_fragments\n" > {output.alignmentstats}
paste {params.outdir}/tmp1 {params.outdir}/tmp2 {params.outdir}/tmp3 {params.outdir}/tmp4 >> {output.alignmentstats}
rm -f {params.outdir}/tmp1 {params.outdir}/tmp2 {params.outdir}/tmp3 {params.outdir}/tmp4 
"""

localrules: venn
rule venn:
    input:
        circexplorerout=rules.circExplorer.output.annotations,
        ciriout=rules.ciri.output.cirioutfiltered
    output:
        png=join(WORKDIR,"results","{sample}","{sample}.venn_mqc.png"),
        cirionly=join(WORKDIR,"results","{sample}","{sample}.cirionly.lst"),
        circexploreronly=join(WORKDIR,"results","{sample}","{sample}.circexploreronly.lst"),
        common=join(WORKDIR,"results","{sample}","{sample}.common.lst")
    params:
        script1=join(SCRIPTS_DIR,"venn.R"),
        sample="{sample}",
        outdir=join(WORKDIR,"results","{sample}")
    container: "docker://nciccbr/ccbr_venn:latest"
    threads: getthreads("venn")
    shell:"""
set -exo pipefail
cut -f1 {input.ciriout}|grep -v circRNA_ID > /dev/shm/{params.sample}.ciri.lst
cut -f1-3 {input.circexplorerout}|awk -F"\\t" '{{print $1\":\"$2+1\"|\"$3}}' > /dev/shm/{params.sample}.circExplorer.lst
if [[ "$(cat /dev/shm/{params.sample}.ciri.lst | wc -l)" != "0" ]];then
    if [[ "$(cat /dev/shm/{params.sample}.circExplorer.lst | wc -l)" != "0" ]];then
        2set_venn.R \\
            -l /dev/shm/{params.sample}.ciri.lst \\
            -r /dev/shm/{params.sample}.circExplorer.lst  \\
            -p {output.png} \\
            -m {output.cirionly} \\
            -s {output.circexploreronly} \\
            -c1 "CIRI2" \\
            -c2 "CircExplorer2" \\
            -c {output.common} \\
            -t {params.sample}
    else
        for o in {output}
        do
            touch $o
        done
    fi
fi
rm -f /dev/shm{params.sample}*
"""


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

rule recall_valid_BSJ_split_BAM_by_strand_create_BW:
    input:
        bam=rules.create_circExplorer_BSJ_bam.output.BSJbam,
    output:
        bam=join(WORKDIR,"results","{sample}","STAR2p","BSJs","{sample}.BSJ."+HOST+".bam"),
        plusbam=join(WORKDIR,"results","{sample}","STAR2p","BSJs","{sample}.BSJ."+HOST+".plus.bam"),
        minusbam=join(WORKDIR,"results","{sample}","STAR2p","BSJs","{sample}.BSJ."+HOST+".minus.bam"),
        bed=join(WORKDIR,"results","{sample}","STAR2p","BSJs","{sample}.BSJ."+HOST+".bed"),
        bw=join(WORKDIR,"results","{sample}","STAR2p","BSJs","{sample}.BSJ."+HOST+".bw"),
        plusbw=join(WORKDIR,"results","{sample}","STAR2p","BSJs","{sample}.BSJ."+HOST+".plus.bw"),
        minusbw=join(WORKDIR,"results","{sample}","STAR2p","BSJs","{sample}.BSJ."+HOST+".minus.bw"),
        validBSJbed=join(WORKDIR,"results","{sample}","customBSJs","{sample}.valid_STAR_BSJ.bed")
    params:
        sample="{sample}",
        workdir=WORKDIR,
        memG=getmemG("recall_valid_BSJ_split_BAM_by_strand_create_BW"),
        outdir=join(WORKDIR,"results","{sample}","STAR2p","BSJs"),
        pythonscript=join(SCRIPTS_DIR,"validate_BSJ_reads_and_split_BSJ_bam_by_strand.py"),
        regions=REF_REGIONS,
        randomstr=str(uuid.uuid4())
    threads: getthreads("recall_valid_BSJ_split_BAM_by_strand_create_BW")
    envmodules: TOOLS["samtools"]["version"],TOOLS["bedtools"]["version"],TOOLS["ucsc"]["version"],TOOLS["sambamba"]["version"],TOOLS["python37"]["version"]
    shell:"""
set -exo pipefail
if [ -d /lscratch/${{SLURM_JOB_ID}} ];then
    TMPDIR="/lscratch/${{SLURM_JOB_ID}}"
else
    TMPDIR="/dev/shm/{params.randomstr}"
    mkdir -p $TMPDIR
fi
cd {params.outdir}
bam_basename="$(basename {input.bam})"
while read a b;do
bam="${{bam_basename%.*}}.${{a}}.bam"
plusbam="${{bam_basename%.*}}.${{a}}.plus.bam"
minusbam="${{bam_basename%.*}}.${{a}}.minus.bam"
bed="${{bam_basename%.*}}.${{a}}.bed"
samtools view {input.bam} $b -b > ${{TMPDIR}}/${{bam%.*}}.tmp.bam
sambamba sort --memory-limit={params.memG} --tmpdir=${{TMPDIR}} --nthreads={threads} --out=$bam ${{TMPDIR}}/${{bam%.*}}.tmp.bam
bw="${{bam%.*}}.bw"
bdg="${{bam%.*}}.bdg"
plusbw="${{bam%.*}}.plus.bw"
plusbdg="${{bam%.*}}.plus.bdg"
minusbw="${{bam%.*}}.minus.bw"
minusbdg="${{bam%.*}}.minus.bdg"
sizes="${{bam%.*}}.sizes"
# create strand specific bams
python {params.pythonscript} \
-i $bam \
-p $plusbam \
-m $minusbam \
-b $bed
bedSort $bed $bed
# create bedgraphs
bedtools genomecov -bga -split -ibam $bam > $bdg
bedSort $bdg $bdg
bedtools genomecov -bga -split -ibam $plusbam > $plusbdg
bedSort $plusbdg $plusbdg
bedtools genomecov -bga -split -ibam $minusbam > $minusbdg
bedSort $minusbdg $minusbdg
# create bigwigs
if [ "$(wc -l $bdg|awk '{{print $1}}')" != "0" ];then
    samtools view -H $bam|grep ^@SQ|cut -f2,3|sed "s/SN://g"|sed "s/LN://g" > $sizes
    bedGraphToBigWig $bdg $sizes $bw
    if [ "$(wc -l $plusbdg|awk '{{print $1}}')" != "0" ];then
        bedGraphToBigWig $plusbdg $sizes $plusbw
    else
        touch $plusbw
    fi
    if [ "$(wc -l $minusbdg|awk '{{print $1}}')" != "0" ];then
        bedGraphToBigWig $minusbdg $sizes $minusbw
    else
        touch $minusbw
    fi
else
    touch $bw
fi
rm -f $bdg $sizes

done < {params.regions}
cat {params.sample}.BSJ.*.bed |cut -f1-6 > ${{TMPDIR}}/{sample}.tmp.bed
bedSort ${{TMPDIR}}/{sample}.tmp.bed ${{TMPDIR}}/{sample}.tmp.bed
awk -F"\\t" -v OFS="\\t" '{{$4="S"NR;print}}' ${{TMPDIR}}/{sample}.tmp.bed > {output.validBSJbed}
cd $TMPDIR && rm -f *
"""	

rule add_novel_ciri_BSJs_to_customBSJ:
    input:
        ciriout=rules.ciri.output.ciriout,
        hg38bed=rules.recall_valid_BSJ_split_BAM_by_strand_create_BW.output.bed,
        validBSJbed=rules.recall_valid_BSJ_split_BAM_by_strand_create_BW.output.validBSJbed,
    output:
        novelBSJ=join(WORKDIR,"results","{sample}","customBSJs","{sample}.novel_CIRI_BSJ.bed"),
        cirireadids=temp(join(WORKDIR,"results","{sample}","customBSJs","{sample}.ciri.readids")),
        combinedBSJ=join(WORKDIR,"results","{sample}","customBSJs","{sample}.BSJ.bed")
    params:
        sample="{sample}",
        outdir=join(WORKDIR,"results","{sample}","customBSJs"),
        list_compare_script=join(SCRIPTS_DIR,"_compare_lists.py"),
        randomstr=str(uuid.uuid4())
    envmodules: TOOLS["ucsc"]["version"],TOOLS["python37"]["version"]
    shell:"""
set -exo pipefail
if [ -d /lscratch/${{SLURM_JOB_ID}} ];then
    TMPDIR="/lscratch/${{SLURM_JOB_ID}}"
else
    TMPDIR="/dev/shm/{params.randomstr}"
    mkdir -p $TMPDIR
fi
cd {params.outdir}
grep -v junction_reads_ID {input.ciriout}|awk -F"\\t" '{{print substr($NF,1,length($NF)-1)}}'|tr ',' '\\n'|sort|uniq > {params.sample}.ciri.readids
stardir=$(dirname {input.hg38bed})
cat ${{stardir}}/{params.sample}.BSJ.*.bed|awk -F"\\t" '{{print $NF}}'|tr ',' '\\n'|awk -F"##" '{{print $1}}'|sort|uniq > {params.sample}.circExplorer.readids
python {params.list_compare_script} {params.sample}.ciri.readids {params.sample}.circExplorer.readids 1
mv a_only.lst {params.sample}.ciri.novel.readids && rm -f a_* b_*
head -n1 {input.ciriout} > {params.sample}.ciri.novel.out
while read a;do grep $a {input.ciriout};done < {params.sample}.ciri.novel.readids |sort|uniq >> {params.sample}.ciri.novel.out
grep -v circRNA_ID {params.sample}.ciri.novel.out | awk -F"\\t" -v OFS="\\t" '{{print $2,$3-1,$4,$2":"$3-1"-"$4,$5,$11}}'  > ${{TMPDIR}}/{params.sample}.novel_CIRI_BSJ.tmp.bed
awk -F"\\t" '{{print $4}}' ${{TMPDIR}}/{params.sample}.novel_CIRI_BSJ.tmp.bed > ${{TMPDIR}}/{params.sample}.novel_CIRI_BSJ.tmp.lst
awk -F"\\t" '{{print $1":"$2"-"$3}}' {input.validBSJbed} > ${{TMPDIR}}/{params.sample}.circExplorer.BSJ.lst
python {params.list_compare_script} ${{TMPDIR}}/{params.sample}.novel_CIRI_BSJ.tmp.lst ${{TMPDIR}}/{params.sample}.circExplorer.BSJ.lst 1
mv a_only.lst {params.sample}.ciri.novel.BSJ.lst && rm -f a_* b_*
while read a;do grep $a ${{TMPDIR}}/{params.sample}.novel_CIRI_BSJ.tmp.bed;done < {params.sample}.ciri.novel.BSJ.lst > ${{TMPDIR}}/{params.sample}.ciri.novel.BSJ.bed
bedSort ${{TMPDIR}}/{params.sample}.ciri.novel.BSJ.bed ${{TMPDIR}}/{params.sample}.ciri.novel.BSJ.bed
awk -F"\\t" -v OFS="\\t" '{{$4="C"NR;print}}' ${{TMPDIR}}/{params.sample}.ciri.novel.BSJ.bed > {output.novelBSJ}
cat {input.validBSJbed} {output.novelBSJ} > {output.combinedBSJ}
bedSort {output.combinedBSJ} {output.combinedBSJ}
rm -f {params.sample}.circExplorer.readids {params.sample}.ciri.novel.*
cd $TMPDIR && rm -f *
"""

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