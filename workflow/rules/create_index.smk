AWK1 = r"""-F"/" '{print $NF}'"""
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
        fastalst=join(REF_DIR,"separate_fastas","separate_fastas.lst"),
        ndx=join(REF_DIR,"NCLscan_index","AllRef.ndx"),
    params:
        reffa=REF_FA,
        refgtf=REF_GTF,
        refdir=REF_DIR,
        script1=join(SCRIPTS_DIR,"_add_geneid2genepred.py"),
        script2=join(SCRIPTS_DIR,"_multifasta2separatefastas.sh"),
        script3=join(SCRIPTS_DIR,"fix_gtfs.py"),
        randomstr=str(uuid.uuid4()),
        nclscan_dir=config['nclscan_dir'],
        nclscan_config=config['nclscan_config'],
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
ls {params.refdir}/separate_fastas/*.fa | awk {AWK1} > {output.fastalst}
# may have to create bowtie1 index here.. has to be a separate rule ... see below


# NCLscan files
python {params.script3} --ingtf {params.refgtf} --outgtf {output.fixed_gtf}
gffread -w {output.transcripts_fa} -g {params.reffa} {output.fixed_gtf}
touch {output.lncRNA_transcripts_fa}
{params.nclscan_dir}/bin/create_reference.py -c {params.nclscan_config}
"""

TRSED = r"""tr '\n' ',' | sed 's/.$//g'"""
rule create_mapsplice_index:
    input:
        fastalst=rules.create_index.output.fastalst
    output:
        rev1ebwt=join(REF_DIR,"separate_fastas_index.rev.1.ebwt"),
    params:
        separate_fastas=join(REF_DIR,"separate_fastas"),
        ebwt=join(REF_DIR,"separate_fastas_index"),
    threads: getthreads("create_mapsplice_index")
    container: "docker://cgrlab/mapsplice2:latest"
    shell:"""
set -exo pipefail
fastas=$(ls {params.separate_fastas}/*.fa| {TRSED})
/opt/MapSplice2/bin/bowtie-build \
$fastas \
{params.ebwt}
"""
