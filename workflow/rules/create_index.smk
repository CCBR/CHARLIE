rule create_index:
    input:
        FASTAS_REGIONS_GTFS
    output:
        genepred_w_geneid=join(REF_DIR,"ref.genes.genepred_w_geneid"),
        sa=join(REF_DIR,"STAR_no_GTF","SA"),
        bwt=join(REF_DIR,"ref.sa")
    params:
        refdir=REF_DIR,
        script1=join(SCRIPTS_DIR,"_add_geneid2genepred.py"),
        randomstr=str(uuid.uuid4())
    envmodules: TOOLS["star"]["version"], TOOLS["bwa"]["version"], TOOLS["samtools"]["version"], TOOLS["ucsc"]["version"]
    threads: getthreads("create_index")
    shell:"""
set -exo pipefail
cd {params.refdir}
samtools faidx ref.fa && \
    cut -f1-2 ref.fa.fai > ref.fa.sizes
bwa index -p ref ref.fa > bwa_index.log
gtfToGenePred -ignoreGroupsWithoutExons ref.gtf ref.genes.genepred && \
    python {params.script1} ref.gtf ref.genes.genepred > ref.genes.genepred_w_geneid
mkdir -p STAR_no_GTF && \
    STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir ./STAR_no_GTF --genomeFastaFiles ref.fa > STAR-build.log 
"""