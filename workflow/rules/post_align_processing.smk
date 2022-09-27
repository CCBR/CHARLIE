rule create_BSJ_bam:
    input:
        junction=rules.star2p.output.junction,
        bam=rules.star2p.output.bam
    output:
        readids=join(WORKDIR,"results","{sample}","STAR2p","{sample}.BSJ.readids"),
        bam=join(WORKDIR,"results","{sample}","STAR2p","{sample}.BSJ.bam")
    params:
        sample="{sample}",
        memG=getmemG("create_BSJ_bam"),
        script1=join(SCRIPTS_DIR,"junctions2readids.py"),
        script2=join(SCRIPTS_DIR,"filter_bam_by_readids.py"),
        script3=join(SCRIPTS_DIR,"filter_bam_for_BSJs.py"),
        outdir=join(WORKDIR,"results","{sample}","STAR2p"),
        randomstr=str(uuid.uuid4())
    envmodules: TOOLS["python37"]["version"],TOOLS["sambamba"]["version"],TOOLS["samtools"]["version"]
    threads : getthreads("create_BSJ_bam")
    shell:"""
set -exo pipefail
if [ -d /lscratch/${{SLURM_JOB_ID}} ];then
    TMPDIR="/lscratch/${{SLURM_JOB_ID}}"
else
    TMPDIR="/dev/shm/{params.randomstr}"
    mkdir -p $TMPDIR
fi
cd {params.outdir}

## get BSJ readids along with chrom,site,cigar etc.
python {params.script1} -j {input.junction} > {output.readids}

## extract only the uniq readids
cut -f1 {output.readids} | sort | uniq > ${{TMPDIR}}/{params.sample}.readids

# INVALID_INDEX_FILE_POINTER 1 issues give errors hence re-indexing
samtools index {input.bam}
## downsize the star2p bam file to a new bam file with only BSJ reads ... these may still contain alignments which are chimeric but not BSJ
## note the argument --readids here is just a list of readids
python {params.script2} --inputBAM {input.bam} --outputBAM ${{TMPDIR}}/{params.sample}.chimeric.bam --readids ${{TMPDIR}}/{params.sample}.readids
sambamba sort --memory-limit={params.memG} --tmpdir=${{TMPDIR}} --nthreads={threads} --out=${{TMPDIR}}/{params.sample}.chimeric.sorted.bam ${{TMPDIR}}/{params.sample}.chimeric.bam
rm -f ${{TMPDIR}}/{params.sample}.chimeric.bam*

## using the downsized star2p bam file containing chimeric alignments ...included all the BSJs... we now extract only the BSJs
## note the argument --readids here is a tab delimited file created by junctions2readids.py ... reaids,chrom,strand,sites,cigars,etc.
python {params.script3} --inputBAM ${{TMPDIR}}/{params.sample}.chimeric.sorted.bam --outputBAM ${{TMPDIR}}/{params.sample}.BSJs.tmp.bam --readids {output.readids}
sambamba sort --memory-limit={params.memG} --tmpdir=${{TMPDIR}} --nthreads={threads} --out=${{TMPDIR}}/{params.sample}.BSJs.tmp.sorted.bam ${{TMPDIR}}/{params.sample}.BSJs.tmp.bam
rm -f ${{TMPDIR}}/{params.sample}.BSJs.tmp.bam*

## some alignments are repeated/duplicated in the output for some reason ... hence deduplicating
samtools view -H ${{TMPDIR}}/{params.sample}.BSJs.tmp.sorted.bam > ${{TMPDIR}}/{params.sample}.BSJs.tmp.dedup.sam
samtools view ${{TMPDIR}}/{params.sample}.BSJs.tmp.sorted.bam | sort | uniq >> ${{TMPDIR}}/{params.sample}.BSJs.tmp.dedup.sam
samtools view -bS ${{TMPDIR}}/{params.sample}.BSJs.tmp.dedup.sam > ${{TMPDIR}}/{params.sample}.BSJs.tmp.dedup.bam
sambamba sort --memory-limit={params.memG} --tmpdir=${{TMPDIR}} --nthreads={threads} --out={output.bam} ${{TMPDIR}}/{params.sample}.BSJs.tmp.dedup.bam
cd ${{TMPDIR}} && rm -f *
"""

rule create_spliced_reads_bam:
    input:
        bam=rules.star2p.output.bam,
        tab=rules.merge_SJ_tabs.output.pass1sjtab
    output:
        bam=join(WORKDIR,"results","{sample}","STAR2p","{sample}.spliced_reads.bam")
    params:
        sample="{sample}",
        memG=getmemG("create_spliced_reads_bam"),
        script1=join(SCRIPTS_DIR,"filter_bam_for_splice_reads.py"),
        outdir=join(WORKDIR,"results","{sample}","STAR2p"),
        randomstr=str(uuid.uuid4())
    envmodules: TOOLS["python37"]["version"],TOOLS["sambamba"]["version"],TOOLS["samtools"]["version"]
    threads : getthreads("create_spliced_reads_bam")
    shell:"""
set -exo pipefail
if [ -d /lscratch/${{SLURM_JOB_ID}} ];then
    TMPDIR="/lscratch/${{SLURM_JOB_ID}}"
else
    TMPDIR="/dev/shm/{params.randomstr}"
    mkdir -p $TMPDIR
fi
cd {params.outdir}
python {params.script1} --inbam {input.bam} --outbam ${{TMPDIR}}/{params.sample}.SR.bam --tab {input.tab}
sambamba sort --memory-limit={params.memG} --tmpdir=${{TMPDIR}} --nthreads={threads} --out={output.bam} ${{TMPDIR}}/{params.sample}.SR.bam
cd $TMPDIR && rm -f *
"""

rule create_linear_reads_bam:
    input:
        junction=rules.star2p.output.junction, # Chimeric junctions file
        bam=rules.star2p.output.bam
    output:
        bam=temp(join(WORKDIR,"results","{sample}","STAR2p","{sample}.linear_reads.bam"))
    params:
        sample="{sample}",
        memG=getmemG("create_linear_reads_bam"),
        peorse=get_peorse,
        script1=join(SCRIPTS_DIR,"filter_bam_for_linear_reads.py"),
        outdir=join(WORKDIR,"results","{sample}","STAR2p"),
        randomstr=str(uuid.uuid4())
    envmodules: TOOLS["python37"]["version"],TOOLS["sambamba"]["version"],TOOLS["samtools"]["version"]
    threads : getthreads("create_linear_reads_bam")
    shell:"""
set -exo pipefail
if [ -d /lscratch/${{SLURM_JOB_ID}} ];then
    TMPDIR="/lscratch/${{SLURM_JOB_ID}}"
else
    TMPDIR="/dev/shm/{params.randomstr}"
    mkdir -p $TMPDIR
fi
cd {params.outdir}
if [ "{params.peorse}" == "PE" ];then
python {params.script1} --inputBAM {input.bam} --outputBAM ${{TMPDIR}}/{params.sample}.tmp.bam -j {input.junction} -p
else 
python {params.script1} --inputBAM {input.bam} --outputBAM ${{TMPDIR}}/{params.sample}.tmp.bam -j {input.junction}
fi
sambamba sort --memory-limit={params.memG} --tmpdir=${{TMPDIR}} --nthreads={threads} --out={output.bam} ${{TMPDIR}}/{params.sample}.tmp.bam
cd $TMPDIR && rm -f *
"""



rule split_splice_reads_BAM_create_BW:
    input:
        bam=rules.create_spliced_reads_bam.output.bam
    output:
        bam=join(WORKDIR,"results","{sample}","STAR2p","{sample}.spliced_reads."+HOST+".bam"),
        bw=join(WORKDIR,"results","{sample}","STAR2p","{sample}.spliced_reads."+HOST+".bw")
    params:
        sample="{sample}",
        workdir=WORKDIR,
        memG=getmemG("split_splice_reads_BAM_create_BW"),
        outdir=join(WORKDIR,"results","{sample}","STAR2p"),
        regions=REF_REGIONS,
        randomstr=str(uuid.uuid4())
    threads: getthreads("split_splice_reads_BAM_create_BW")
    envmodules: TOOLS["samtools"]["version"],TOOLS["bedtools"]["version"],TOOLS["ucsc"]["version"],TOOLS["sambamba"]["version"]
    shell:"""
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
samtools view {input.bam} $b -b > ${{TMPDIR}}/${{bam%.*}}.tmp.bam
sambamba sort --memory-limit={params.memG} --tmpdir=${{TMPDIR}} --nthreads={threads} --out=$bam ${{TMPDIR}}/${{bam%.*}}.tmp.bam
bw="${{bam%.*}}.bw"
bdg="${{bam%.*}}.bdg"
sizes="${{bam%.*}}.sizes"
bedtools genomecov -bga -split -ibam $bam > $bdg
bedSort $bdg $bdg
if [ "$(wc -l $bdg|awk '{{print $1}}')" != "0" ];then
samtools view -H $bam|grep ^@SQ|cut -f2,3|sed "s/SN://g"|sed "s/LN://g" > $sizes
bedGraphToBigWig $bdg $sizes $bw
else
touch $bw
fi
rm -f $bdg $sizes
done < {params.regions}
cd $TMPDIR && rm -f *
"""	

rule split_linear_reads_BAM_create_BW:
# This rule is identical to split_splice_reads_BAM_create_BW, "spliced_reads" is replaced by "linear_reads"
    input:
        bam=rules.create_linear_reads_bam.output.bam
    output:
        bam=join(WORKDIR,"results","{sample}","STAR2p","{sample}.linear_reads."+HOST+".bam"),
        bw=join(WORKDIR,"results","{sample}","STAR2p","{sample}.linear_reads."+HOST+".bw")
    params:
        sample="{sample}",
        workdir=WORKDIR,
        memG=getmemG("split_linear_reads_BAM_create_BW"),
        outdir=join(WORKDIR,"results","{sample}","STAR2p"),
        regions=REF_REGIONS,
        randomstr=str(uuid.uuid4())
    threads: getthreads("split_linear_reads_BAM_create_BW")
    envmodules: TOOLS["samtools"]["version"],TOOLS["bedtools"]["version"],TOOLS["ucsc"]["version"],TOOLS["sambamba"]["version"]
    shell:"""
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
samtools view {input.bam} $b -b > ${{TMPDIR}}/${{bam%.*}}.tmp.bam
sambamba sort --memory-limit={params.memG} --tmpdir=${{TMPDIR}} --nthreads={threads} --out=$bam ${{TMPDIR}}/${{bam%.*}}.tmp.bam
bw="${{bam%.*}}.bw"
bdg="${{bam%.*}}.bdg"
sizes="${{bam%.*}}.sizes"
bedtools genomecov -bga -split -ibam $bam > $bdg
bedSort $bdg $bdg
if [ "$(wc -l $bdg|awk '{{print $1}}')" != "0" ];then
samtools view -H $bam|grep ^@SQ|cut -f2,3|sed "s/SN://g"|sed "s/LN://g" > $sizes
bedGraphToBigWig $bdg $sizes $bw
else
touch $bw
fi
rm -f $bdg $sizes
done < {params.regions}
cd $TMPDIR && rm -f *
"""	

rule alignment_stats:
    input:
        star2pbam=rules.star2p.output.bam,
        splicedreadsbam=rules.create_spliced_reads_bam.output.bam,
        BSJbam=rules.create_BSJ_bam.output.bam,
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

localrules: merge_genecounts
rule merge_genecounts:
    input:
        expand(join(WORKDIR,"results","{sample}","STAR2p","{sample}_p2.ReadsPerGene.out.tab"),sample=SAMPLES)
    output:
        join(WORKDIR,"results","unstranded_STAR_GeneCounts.tsv"),
        join(WORKDIR,"results","stranded_STAR_GeneCounts.tsv"),
        join(WORKDIR,"results","revstranded_STAR_GeneCounts.tsv")
    params:
        outdir=join(WORKDIR,"results"),
        rscript=join(SCRIPTS_DIR,"merge_ReadsPerGene_counts.R")
    envmodules: TOOLS["R"]["version"]
    shell:"""
set -exo pipefail
cd {params.outdir}
Rscript {params.rscript}
"""