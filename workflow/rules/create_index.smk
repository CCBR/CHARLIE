rule create_index:
    input:
        FASTAS_REGIONS_GTFS
    output:
        genepred_w_geneid=join(REF_DIR,"ref.genes.genepred_w_geneid"),
        sa=join(REF_DIR,"STAR_no_GTF","SA"),
        bwt=join(REF_DIR,"ref.sa"),
        fixed_gtf=join(REF_DIR,"ref.fixed.gtf"),
        transcripts_fa=join(REF_DIR,"ref.transcripts.fa"),
        lncRNA_transcripts_fa=join(REF_DIR,"ref.dummy.fa"),
    params:
        reffa=REF_FA,
        refgtf=REF_GTF,
        refdir=REF_DIR,
        script1=join(SCRIPTS_DIR,"_add_geneid2genepred.py"),
        script2=join(SCRIPTS_DIR,"_multifasta2separatefastas.sh"),
        script3=join(SCRIPTS_DIR,"fix_gtfs.py"),
        randomstr=str(uuid.uuid4())
    envmodules: TOOLS["star"]["version"], TOOLS["bwa"]["version"], TOOLS["samtools"]["version"], TOOLS["ucsc"]["version"], TOOLS["cufflinks"]["version"]
    threads: getthreads("create_index")
    shell:"""
set -exo pipefail
cd {params.refdir}
samtools faidx {params.reffa} && \
    cut -f1-2 {params.reffa}.fai > {params.reffa}.sizes
bwa index -p ref {params.reffa} > bwa_index.log
gtfToGenePred -ignoreGroupsWithoutExons {params.refgtf} ref.genes.genepred && \
    python {params.script1} {params.refgtf} ref.genes.genepred > {output.genepred_w_geneid}
stardir=$(dirname {output.sa})
mkdir -p STAR_no_GTF && \
    STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir $stardir --genomeFastaFiles {params.reffa} > STAR-build.log
# MapSplice requires the {params.reffa} multifasta to be split into separate fastas
bash {params.script2} {params.reffa} {params.refdir}/separate_fastas
# may have to create bowtie1 index here

# NCLscan files
python {params.script3} --ingtf {params.refgtf} --outgtf {output.fixed_gtf}

gffread -w {output.transcripts_fa} -g {params.reffa} {output.fixed_gtf}

touch {output.lncRNA_transcripts_fa}
"""