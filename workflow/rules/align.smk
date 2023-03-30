## rules
#  --outSJfilterOverhangMin 15 15 15 15
#  --alignSJoverhangMin 15
#  --alignSJDBoverhangMin 15
#  --outFilterScoreMin 1
#  --outFilterMatchNmin 1
#  --outFilterMismatchNmax 2
#  --chimSegmentMin 15 
#  --chimScoreMin 15
#  --chimJunctionOverhangMin 15 --> --alignSJoverhangMin and --chimJunctionOverhangMin should use the same value to make the circRNA expression and linear gene expression level comparable.
#  --seedSearchStartLmax 30 ... for SE only
#  ref: https://github.com/dieterich-lab/DCC
rule star1p:
    input:
        sa=rules.create_index.output.sa,
        R1=rules.cutadapt.output.of1,
        R2=rules.cutadapt.output.of2,
        gtf=rules.create_index.output.fixed_gtf,
    output:
        junction=join(WORKDIR,"results","{sample}","STAR1p","{sample}_p1.SJ.out.tab"),
        chimeric_junctions=join(WORKDIR,"results","{sample}","STAR1p","{sample}_p1.Chimeric.out.junction"),
        mate1_chimeric_junctions=join(WORKDIR,"results","{sample}","STAR1p","mate1","{sample}"+"_mate1.Chimeric.out.junction"),
        mate2_chimeric_junctions=join(WORKDIR,"results","{sample}","STAR1p","mate2","{sample}"+"_mate2.Chimeric.out.junction"),
        # get_mate_outputs
    params:
        sample="{sample}",
        peorse=get_peorse,
        workdir=WORKDIR,
        outdir=join(WORKDIR,"results","{sample}","STAR1p"),
        starindexdir=STAR_INDEX_DIR,
        alignTranscriptsPerReadNmax=TOOLS["star"]["alignTranscriptsPerReadNmax"],
        randomstr=str(uuid.uuid4())
    envmodules: TOOLS["star"]["version"]
    threads: getthreads("star1p")
    shell:"""
set -exo pipefail
if [ -d /lscratch/${{SLURM_JOB_ID}} ];then
    TMPDIR="/lscratch/${{SLURM_JOB_ID}}/{params.randomstr}"
else
    TMPDIR="/dev/shm/{params.randomstr}"
fi

if [ ! -d {params.outdir} ];then mkdir {params.outdir};fi
if [ "{params.peorse}" == "PE" ];then
# paired-end
    overhang=$(zcat {input.R1} {input.R2} | awk -v maxlen=100 'NR%4==2 {{if (length($1) > maxlen+0) maxlen=length($1)}}; END {{print maxlen-1}}')
    echo "sjdbOverhang for STAR: ${{overhang}}"
    cd {params.outdir}
    STAR --genomeDir {params.starindexdir} \\
    --outSAMstrandField None  \\
    --outFilterMultimapNmax 20 \\
    --outSJfilterOverhangMin 15 15 15 15 \\
    --alignSJoverhangMin 15 \\
    --alignSJDBoverhangMin 15 \\
    --outFilterScoreMin 1 \\
    --outFilterMatchNmin 1 \\
    --outFilterMismatchNmax 2 \\
    --outFilterMismatchNoverLmax 0.3  \\
    --alignIntronMin 20 \\
    --alignIntronMax 1000000 \\
    --alignMatesGapMax 1000000 \\
    --readFilesIn {input.R1} {input.R2} \\
    --readFilesCommand zcat \\
    --runThreadN {threads} \\
    --outFileNamePrefix {params.sample}_p1. \\
    --chimSegmentMin 15 \\
    --chimScoreMin 15 \\
    --chimJunctionOverhangMin 15 \\
    --chimMultimapNmax 10 \\
    --chimOutType Junctions \\
    --alignTranscriptsPerReadNmax {params.alignTranscriptsPerReadNmax} \\
    --outSAMtype None \\
    --alignEndsProtrude 10 ConcordantPair \\
    --outFilterIntronMotifs None \\
    --sjdbGTFfile {input.gtf} \\
    --outTmpDir ${{TMPDIR}} \\
    --sjdbOverhang $overhang

    # mate1
    mkdir -p "{params.outdir}/mate1" && cd {params.outdir}/mate1
    STAR --genomeDir {params.starindexdir} \\
    --outSAMstrandField None  \\
    --outFilterMultimapNmax 20 \\
    --outSJfilterOverhangMin 15 15 15 15 \\
    --alignSJoverhangMin 15 \\
    --alignSJDBoverhangMin 15 \\
    --seedSearchStartLmax 30 \\
    --outFilterScoreMin 1 \\
    --outFilterMatchNmin 1 \\
    --outFilterMismatchNmax 2 \\
    --outFilterMismatchNoverLmax 0.3  \\
    --alignIntronMin 20 \\
    --alignIntronMax 1000000 \\
    --alignMatesGapMax 1000000 \\
    --readFilesIn {input.R1} \\
    --readFilesCommand zcat \\
    --runThreadN {threads} \\
    --outFileNamePrefix {params.sample}_mate1. \\
    --chimSegmentMin 15 \\
    --chimScoreMin 15 \\
    --chimJunctionOverhangMin 15 \\
    --chimMultimapNmax 10 \\
    --chimOutType Junctions \\
    --alignTranscriptsPerReadNmax {params.alignTranscriptsPerReadNmax} \\
    --outSAMtype None \\
    --alignEndsProtrude 10 ConcordantPair \\
    --outFilterIntronMotifs None \\
    --sjdbGTFfile {input.gtf} \\
    --outTmpDir ${{TMPDIR}} \\
    --sjdbOverhang $overhang

    # mate2
    mkdir -p "{params.outdir}/mate2" && cd {params.outdir}/mate2
    STAR --genomeDir {params.starindexdir} \\
    --outSAMstrandField None  \\
    --outFilterMultimapNmax 20 \\
    --outSJfilterOverhangMin 15 15 15 15 \\
    --alignSJoverhangMin 15 \\
    --alignSJDBoverhangMin 15 \\
    --seedSearchStartLmax 30 \\
    --outFilterScoreMin 1 \\
    --outFilterMatchNmin 1 \\
    --outFilterMismatchNmax 2 \\
    --outFilterMismatchNoverLmax 0.3  \\
    --alignIntronMin 20 \\
    --alignIntronMax 1000000 \\
    --alignMatesGapMax 1000000 \\
    --readFilesIn {input.R2} \\
    --readFilesCommand zcat \\
    --runThreadN {threads} \\
    --outFileNamePrefix {params.sample}_mate2. \\
    --chimSegmentMin 15 \\
    --chimScoreMin 15 \\
    --chimJunctionOverhangMin 15 \\
    --chimMultimapNmax 10 \\
    --chimOutType Junctions \\
    --alignTranscriptsPerReadNmax {params.alignTranscriptsPerReadNmax} \\
    --outSAMtype None \\
    --alignEndsProtrude 10 ConcordantPair \\
    --outFilterIntronMotifs None \\
    --sjdbGTFfile {input.gtf} \\
    --outTmpDir ${{TMPDIR}} \\
    --sjdbOverhang $overhang

else

#single-end
    overhang=$(zcat {input.R1} | awk -v maxlen=100 'NR%4==2 {{if (length($1) > maxlen+0) maxlen=length($1)}}; END {{print maxlen-1}}')
    echo "sjdbOverhang for STAR: ${{overhang}}"
    cd {params.outdir}
    STAR --genomeDir {params.starindexdir} \\
    --outSAMstrandField None  \\
    --outFilterMultimapNmax 20 \\
    --outSJfilterOverhangMin 15 15 15 15 \\
    --alignSJoverhangMin 15 \\
    --alignSJDBoverhangMin 15 \\
    --seedSearchStartLmax 30 \\
    --outFilterScoreMin 1 \\
    --outFilterMatchNmin 1 \\
    --outFilterMismatchNmax 2 \\
    --outFilterMismatchNoverLmax 0.3  \\
    --alignIntronMin 20 \\
    --alignIntronMax 1000000 \\
    --alignMatesGapMax 1000000 \\
    --readFilesIn {input.R1} \\
    --readFilesCommand zcat \\
    --runThreadN {threads} \\
    --outFileNamePrefix {params.sample}_p1. \\
    --chimSegmentMin 15 \\
    --chimScoreMin 15 \\
    --chimJunctionOverhangMin 15 \\
    --chimMultimapNmax 10 \\
    --chimOutType Junctions \\
    --alignTranscriptsPerReadNmax {params.alignTranscriptsPerReadNmax} \\
    --outSAMtype None \\
    --alignEndsProtrude 10 ConcordantPair \\
    --outFilterIntronMotifs None \\
    --sjdbGTFfile {input.gtf} \\
    --outTmpDir ${{TMPDIR}} \\
    --sjdbOverhang $overhang
    mkdir -p $(dirname {output.mate1_chimeric_junctions})
    touch {output.mate1_chimeric_junctions}
    mkdir -p $(dirname {output.mate2_chimeric_junctions})
    touch {output.mate2_chimeric_junctions}
fi
rm -rf {params.outdir}/{params.sample}_p1.Aligned.out.bam
"""

rule merge_SJ_tabs:
    input:
        expand(join(WORKDIR,"results","{sample}","STAR1p","{sample}_p1.SJ.out.tab"), sample=SAMPLES)
    output:
        pass1sjtab=join(WORKDIR,"results","pass1.out.tab")
    params:
        script1=join(SCRIPTS_DIR,"apply_junction_filters.py"),
        regions=REF_REGIONS,
        filter1regions=HOST_ADDITIVES,
        filter1_noncanonical=config['star_1pass_filter_host_noncanonical'],
        filter1_unannotated=config['star_1pass_filter_host_unannotated'],
        filter2_noncanonical=config['star_1pass_filter_viruses_noncanonical'],
        filter2_unannotated=config['star_1pass_filter_viruses_unannotated']
    threads: getthreads("merge_SJ_tabs")
    shell:"""
set -exo pipefail
cat {input} | \\
python {params.script1} \\
    --regions {params.regions} \\
    --filter1regions {params.filter1regions} \\
    --filter1_noncanonical {params.filter1_noncanonical} \\
    --filter1_unannotated {params.filter1_unannotated} \\
    --filter2_noncanonical {params.filter2_noncanonical} \\
    --filter2_unannotated {params.filter1_unannotated} | \\
cut -f1-4 | sort -k1,1 -k2,2n | uniq > {output.pass1sjtab}
"""

rule star2p:
    input:
        R1=rules.cutadapt.output.of1,
        R2=rules.cutadapt.output.of2,
        gtf=rules.create_index.output.fixed_gtf,
        pass1sjtab=rules.merge_SJ_tabs.output.pass1sjtab
    output:
        junction=join(WORKDIR,"results","{sample}","STAR2p","{sample}_p2.Chimeric.out.junction"),
        unsortedbam=temp(join(WORKDIR,"results","{sample}","STAR2p","{sample}_p2.Aligned.out.bam")),
        bam=join(WORKDIR,"results","{sample}","STAR2p","{sample}_p2.bam"),
        chimeric_bam=join(WORKDIR,"results","{sample}","STAR2p","{sample}_p2.chimeric.bam"),
        non_chimeric_bam=join(WORKDIR,"results","{sample}","STAR2p","{sample}_p2.non_chimeric.bam"),
        genecounts=join(WORKDIR,"results","{sample}","STAR2p","{sample}_p2.ReadsPerGene.out.tab")
    params:
        sample="{sample}",
        memG=getmemG("star2p"),
        peorse=get_peorse,
        workdir=WORKDIR,
        outdir=join(WORKDIR,"results","{sample}","STAR2p"),
        starindexdir=STAR_INDEX_DIR,
        alignTranscriptsPerReadNmax=TOOLS["star"]["alignTranscriptsPerReadNmax"],
        randomstr=str(uuid.uuid4())
    envmodules: TOOLS["star"]["version"],TOOLS["sambamba"]["version"], TOOLS["samtools"]["version"]
    threads: getthreads("star2p")
    shell:"""
set -exo pipefail
if [ -d /lscratch/${{SLURM_JOB_ID}} ];then
    TMPDIR="/lscratch/${{SLURM_JOB_ID}}/{params.randomstr}"
else
    TMPDIR="/dev/shm/{params.randomstr}"
fi

if [ ! -d {params.outdir} ];then mkdir {params.outdir};fi
limitSjdbInsertNsj=$(wc -l {input.pass1sjtab}|awk '{{print $1+1}}')
if [ "$limitSjdbInsertNsj" -lt "400000" ];then limitSjdbInsertNsj="400000";fi
if [ "{params.peorse}" == "PE" ];then
# paired-end
    overhang=$(zcat {input.R1} {input.R2} | awk -v maxlen=100 'NR%4==2 {{if (length($1) > maxlen+0) maxlen=length($1)}}; END {{print maxlen-1}}')
    echo "sjdbOverhang for STAR: ${{overhang}}"
    cd {params.outdir}
    STAR --genomeDir {params.starindexdir} \\
    --outSAMstrandField None  \\
    --outFilterType BySJout \\
    --outFilterMultimapNmax 20 \\
    --outSJfilterOverhangMin 15 15 15 15 \\
    --alignSJoverhangMin 15 \\
    --alignSJDBoverhangMin 15 \\
    --outFilterScoreMin 1 \\
    --outFilterMatchNmin 1 \\
    --outFilterMismatchNmax 2 \\
    --outFilterMismatchNoverLmax 0.3  \\
    --alignIntronMin 20 \\
    --alignIntronMax 2000000 \\
    --alignMatesGapMax 2000000 \\
    --readFilesIn {input.R1} {input.R2} \\
    --readFilesCommand  zcat \\
    --runThreadN {threads} \\
    --outFileNamePrefix {params.sample}_p2. \\
    --sjdbFileChrStartEnd {input.pass1sjtab} \\
    --chimSegmentMin 15 \\
    --chimScoreMin 15 \\
    --chimJunctionOverhangMin 15 \\
    --chimOutType Junctions WithinBAM \\
    --chimMultimapNmax 10 \\
    --limitSjdbInsertNsj $limitSjdbInsertNsj \\
    --alignTranscriptsPerReadNmax {params.alignTranscriptsPerReadNmax} \\
    --outSAMtype BAM Unsorted \\
    --alignEndsProtrude 10 ConcordantPair \\
    --outFilterIntronMotifs None \\
    --sjdbGTFfile {input.gtf} \\
    --quantMode GeneCounts \\
    --outTmpDir ${{TMPDIR}} \\
    --sjdbOverhang $overhang \\
    --outBAMcompression 0 \\
    --outSAMattributes All
else
#single-end
    overhang=$(zcat {input.R1} | awk -v maxlen=100 'NR%4==2 {{if (length($1) > maxlen+0) maxlen=length($1)}}; END {{print maxlen-1}}')
    echo "sjdbOverhang for STAR: ${{overhang}}"
    cd {params.outdir}
    STAR --genomeDir {params.starindexdir} \\
    --outSAMstrandField None  \\
    --outFilterType BySJout \\
    --outFilterMultimapNmax 20 \\
    --outSJfilterOverhangMin 15 15 15 15 \\
    --alignSJoverhangMin 15 \\
    --alignSJDBoverhangMin 15 \\
    --seedSearchStartLmax 30 \\
    --outFilterScoreMin 1 \\
    --outFilterMatchNmin 1 \\
    --outFilterMismatchNmax 2 \\
    --outFilterMismatchNoverLmax 0.3  \\
    --alignIntronMin 20 \\
    --alignIntronMax 2000000 \\
    --alignMatesGapMax 2000000 \\
    --readFilesIn {input.R1} \\
    --readFilesCommand  zcat \\
    --runThreadN {threads} \\
    --outFileNamePrefix {params.sample}_p2. \\
    --sjdbFileChrStartEnd {input.pass1sjtab} \\
    --chimSegmentMin 15 \\
    --chimScoreMin 15 \\
    --chimJunctionOverhangMin 15 \\
    --chimOutType Junctions WithinBAM \\
    --chimMultimapNmax 10 \\
    --limitSjdbInsertNsj $limitSjdbInsertNsj \\
    --alignTranscriptsPerReadNmax {params.alignTranscriptsPerReadNmax} \\
    --outSAMtype BAM Unsorted \\
    --alignEndsProtrude 10 ConcordantPair \\
    --outFilterIntronMotifs None \\
    --sjdbGTFfile {input.gtf} \\
    --quantMode GeneCounts \\
    --outTmpDir ${{TMPDIR}} \\
    --sjdbOverhang $overhang \\
    --outBAMcompression 0 \\
    --outSAMattributes All
fi
sleep 120
if [ ! -d $TMPDIR ];then mkdir -p $TMPDIR;fi
samtools view -H {output.unsortedbam} > ${{TMPDIR}}/{params.sample}_p2.non_chimeric.sam
cp ${{TMPDIR}}/{params.sample}_p2.non_chimeric.sam ${{TMPDIR}}/{params.sample}_p2.chimeric.sam
# ref https://github.com/alexdobin/STAR/issues/678
samtools view -@ {threads} {output.unsortedbam} | grep "ch:A:1" >> ${{TMPDIR}}/{params.sample}_p2.chimeric.sam
samtools view -@ {threads} {output.unsortedbam} | grep -v "ch:A:1" >> ${{TMPDIR}}/{params.sample}_p2.non_chimeric.sam
ls -alrth
for i in 1 2 3;do
    if [ ! -d ${{TMPDIR}}/{params.randomstr}_${{i}} ];then mkdir -p ${{TMPDIR}}/{params.randomstr}_${{i}};fi
done
samtools view -@ {threads} -b -S ${{TMPDIR}}/{params.sample}_p2.chimeric.sam | \\
samtools sort \\
    -l 9 \\
    -T ${{TMPDIR}}/{params.randomstr}_1 \\
    --write-index \\
    -@ {threads} \\
    --output-fmt BAM \\
    -o {output.chimeric_bam} -
samtools view -@ {threads} -b -S ${{TMPDIR}}/{params.sample}_p2.non_chimeric.sam | \\
samtools sort \\
    -l 9 \\
    -T ${{TMPDIR}}/{params.randomstr}_2 \\
    --write-index \\
    -@ {threads} \\
    --output-fmt BAM \\
    -o {output.non_chimeric_bam} -
samtools sort \\
    -l 9 \\
    -T ${{TMPDIR}}/{params.randomstr}_3 \\
    --write-index \\
    -@ {threads} \\
    --output-fmt BAM \\
    -o {output.bam} {output.unsortedbam}
"""


rule star_circrnafinder:
    input:
        R1=rules.cutadapt.output.of1,
        R2=rules.cutadapt.output.of2,
        gtf=rules.create_index.output.fixed_gtf,
    output:
        chimericsam=join(WORKDIR,"results","{sample}","STAR_circRNAFinder","{sample}.Chimeric.out.sam"),
        chimericjunction=join(WORKDIR,"results","{sample}","STAR_circRNAFinder","{sample}.Chimeric.out.junction"),
        sjouttab=join(WORKDIR,"results","{sample}","STAR_circRNAFinder","{sample}.SJ.out.tab"),
# -rw-rw----  1 kopardevn kopardevn 2.0M Mar 30 00:16 GI1_T.Chimeric.out.sam
# -rw-rw----  1 kopardevn kopardevn 289K Mar 30 00:16 GI1_T.Chimeric.out.junction
# -rw-rw----  1 kopardevn kopardevn  87K Mar 30 00:16 GI1_T.SJ.out.tab
        bam=temp(join(WORKDIR,"results","{sample}","STAR_circRNAFinder","{sample}.Aligned.out.bam"))
    params:
        sample="{sample}",
        memG=getmemG("star2p"),
        peorse=get_peorse,
        workdir=WORKDIR,
        starindexdir=STAR_INDEX_DIR,
        alignTranscriptsPerReadNmax=TOOLS["star"]["alignTranscriptsPerReadNmax"],
        randomstr=str(uuid.uuid4())
    envmodules: TOOLS["star"]["version"],TOOLS["sambamba"]["version"], TOOLS["samtools"]["version"]
    threads: getthreads("star_circrnafinder")
    shell:"""
set -exo pipefail
if [ -d /lscratch/${{SLURM_JOB_ID}} ];then
    TMPDIR="/lscratch/${{SLURM_JOB_ID}}/{params.randomstr}"
else
    TMPDIR="/dev/shm/{params.randomstr}"
fi

outdir=$(dirname {output.bam})
if [ ! -d $outdir ];then mkdir -p $outdir;fi
cd $outdir

if [ "{params.peorse}" == "PE" ];then
# paired-end
    STAR --genomeDir {params.starindexdir} \\
    --readFilesIn {input.R1} {input.R2} \\
    --readFilesCommand  zcat \\
    --runThreadN {threads} \\
	--chimSegmentMin 20 \\
	--chimScoreMin 1 \\
	--alignIntronMax 1000000 \\
	--outFilterMismatchNoverReadLmax 0.02 \\
	--alignTranscriptsPerReadNmax 100000 \\
	--twopassMode Basic \\
	--outSAMtype BAM Unsorted \\
	--chimOutType Junctions SeparateSAMold \\
	--outFilterMultimapNmax 2 \\
	--outFileNamePrefix {params.sample}. \\
    --outBAMcompression 0 \\
    --outTmpDir $TMPDIR \\
	--sjdbGTFfile {input.gtf}

else
#single-end
    STAR --genomeDir {params.starindexdir} \\
    --readFilesIn {input.R1} \\
    --readFilesCommand  zcat \\
    --runThreadN {threads} \\
	--chimSegmentMin 20 \\
	--chimScoreMin 1 \\
	--alignIntronMax 1000000 \\
	--outFilterMismatchNoverReadLmax 0.02 \\
	--alignTranscriptsPerReadNmax 100000 \\
	--twopassMode Basic \\
	--outSAMtype BAM Unsorted \\
	--chimOutType Junctions SeparateSAMold \\
	--outFilterMultimapNmax 2 \\
	--outFileNamePrefix {params.sample}. \\
    --outBAMcompression 0 \\
    --outTmpDir $TMPDIR \\
	--sjdbGTFfile {input.gtf}

fi

sleep 120

# Used to sort the BAM file after effect , but realized that it is not required by circRNA_Finder scripts
# Hence deleting it to save digital footprint by making it temp in output block 

# if [ ! -d $TMPDIR ];then mkdir -p $TMPDIR;fi

# mv {output.bam} {params.sample}.tmp.Aligned.out.bam

# samtools sort \\
#     -l 9 \\
#     -T ${{TMPDIR}}/{params.randomstr}_1 \\
#     --write-index \\
#     -@ {threads} \\
#     --output-fmt BAM \\
#     -o {output.bam} {params.sample}.tmp.Aligned.out.bam && \\
# rm -f {params.sample}.tmp.Aligned.out.bam

# rm -rf $TMPDIR

"""

rule estimate_duplication:
    input:
        bam=rules.star2p.output.bam
    output:
        metrics=join(WORKDIR,"qc","picard_MarkDuplicates","{sample}.MarkDuplicates.metrics.txt")
    params:
        sample="{sample}",
        memG=getmemG("estimate_duplication"),
    envmodules: TOOLS["picard"]["version"], TOOLS["java"]["version"]
    shell:"""
set -exo pipefail
java -Xmx{params.memG} -jar ${{PICARD_JARPATH}}/picard.jar MarkDuplicates I={input.bam} O=/dev/shm/{params.sample}.mark_dup.bam M={output.metrics}
rm -f /dev/shm/{params.sample}*
"""


