
# rule split_splice_reads_BAM_create_BW:
#     input:
#         bam=rules.create_spliced_reads_bam.output.bam
#     output:
#         bam=join(WORKDIR,"results","{sample}","STAR2p","{sample}.spliced_reads."+HOST+".bam"),
#         bw=join(WORKDIR,"results","{sample}","STAR2p","{sample}.spliced_reads."+HOST+".bw")
#     params:
#         sample="{sample}",
#         workdir=WORKDIR,
#         memG=getmemG("split_splice_reads_BAM_create_BW"),
#         outdir=join(WORKDIR,"results","{sample}","STAR2p"),
#         regions=REF_REGIONS,
#         randomstr=str(uuid.uuid4())
#     threads: getthreads("split_splice_reads_BAM_create_BW")
#     envmodules: TOOLS["samtools"]["version"],TOOLS["bedtools"]["version"],TOOLS["ucsc"]["version"],TOOLS["sambamba"]["version"]
#     shell:"""
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
# samtools view {input.bam} $b -b > ${{TMPDIR}}/${{bam%.*}}.tmp.bam
# sambamba sort --memory-limit={params.memG} --tmpdir=${{TMPDIR}} --nthreads={threads} --out=$bam ${{TMPDIR}}/${{bam%.*}}.tmp.bam
# bw="${{bam%.*}}.bw"
# bdg="${{bam%.*}}.bdg"
# sizes="${{bam%.*}}.sizes"
# bedtools genomecov -bga -split -ibam $bam > $bdg
# bedSort $bdg $bdg
# if [ "$(wc -l $bdg|awk '{{print $1}}')" != "0" ];then
# samtools view -H $bam|grep ^@SQ|cut -f2,3|sed "s/SN://g"|sed "s/LN://g" > $sizes
# bedGraphToBigWig $bdg $sizes $bw
# else
# touch $bw
# fi
# rm -f $bdg $sizes
# done < {params.regions}
# cd $TMPDIR && rm -f *
# """	

# rule split_linear_reads_BAM_create_BW:
# # This rule is identical to split_splice_reads_BAM_create_BW, "spliced_reads" is replaced by "linear_reads"
#     input:
#         bam=rules.create_linear_reads_bam.output.bam
#     output:
#         bam=join(WORKDIR,"results","{sample}","STAR2p","{sample}.linear_reads."+HOST+".bam"),
#         bw=join(WORKDIR,"results","{sample}","STAR2p","{sample}.linear_reads."+HOST+".bw")
#     params:
#         sample="{sample}",
#         workdir=WORKDIR,
#         memG=getmemG("split_linear_reads_BAM_create_BW"),
#         outdir=join(WORKDIR,"results","{sample}","STAR2p"),
#         regions=REF_REGIONS,
#         randomstr=str(uuid.uuid4())
#     threads: getthreads("split_linear_reads_BAM_create_BW")
#     envmodules: TOOLS["samtools"]["version"],TOOLS["bedtools"]["version"],TOOLS["ucsc"]["version"],TOOLS["sambamba"]["version"]
#     shell:"""
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
# samtools view {input.bam} $b -b > ${{TMPDIR}}/${{bam%.*}}.tmp.bam
# sambamba sort --memory-limit={params.memG} --tmpdir=${{TMPDIR}} --nthreads={threads} --out=$bam ${{TMPDIR}}/${{bam%.*}}.tmp.bam
# bw="${{bam%.*}}.bw"
# bdg="${{bam%.*}}.bdg"
# sizes="${{bam%.*}}.sizes"
# bedtools genomecov -bga -split -ibam $bam > $bdg
# bedSort $bdg $bdg
# if [ "$(wc -l $bdg|awk '{{print $1}}')" != "0" ];then
# samtools view -H $bam|grep ^@SQ|cut -f2,3|sed "s/SN://g"|sed "s/LN://g" > $sizes
# bedGraphToBigWig $bdg $sizes $bw
# else
# touch $bw
# fi
# rm -f $bdg $sizes
# done < {params.regions}
# cd $TMPDIR && rm -f *
# """	


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