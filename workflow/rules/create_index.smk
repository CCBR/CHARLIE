AWK1 = r"""-F"/" '{print $NF}'"""


rule create_index:
    input:
        # FASTAS_REGIONS_GTFS
        list(map(lambda x: ancient(x), FASTAS_REGIONS_GTFS)),
    output:
        genepred_w_geneid=join(REF_DIR, "ref.genes.genepred_w_geneid"),
        sa=join(REF_DIR, "STAR_no_GTF", "SA"),
        fixed_gtf=join(REF_DIR, "ref.fixed.gtf"),
        transcripts_fa=join(REF_DIR, "ref.transcripts.fa"),
        lncRNA_transcripts_fa=join(REF_DIR, "ref.dummy.fa"),
        fastalst=join(REF_DIR, "separate_fastas", "separate_fastas.lst"),
        ndx=join(REF_DIR, "NCLscan_index", "AllRef.ndx"),
    params:
        reffa=REF_FA,
        refgtf=REF_GTF,
        refdir=REF_DIR,
        script1=join(SCRIPTS_DIR, "_add_geneid2genepred.py"),
        script2=join(SCRIPTS_DIR, "_multifasta2separatefastas.sh"),
        script3=join(SCRIPTS_DIR, "fix_gtfs.py"),
        nclscan_config=config["nclscan_config"],
    container: config['containers']['star_ucsc_cufflinks']
    threads: getthreads("create_index")
    shell:
        """
set -exo pipefail
cd {params.refdir}
samtools faidx {params.reffa} && \\
    cut -f1-2 {params.reffa}.fai > {params.reffa}.sizes

# bwa index -p ref {params.reffa} > bwa_index.log ... created in a separate rule

# NCLscan files
python -E {params.script3} --ingtf {params.refgtf} --outgtf {output.fixed_gtf}
gffread -w {output.transcripts_fa} -g {params.reffa} {output.fixed_gtf}
touch {output.lncRNA_transcripts_fa}
create_reference.py -c {params.nclscan_config}

gtfToGenePred -ignoreGroupsWithoutExons {output.fixed_gtf} ref.genes.genepred && \\
    python -E {params.script1} {output.fixed_gtf} ref.genes.genepred > {output.genepred_w_geneid}

stardir=$(dirname {output.sa})
mkdir -p $stardir && \\
STAR \\
    --runThreadN {threads} \\
    --runMode genomeGenerate \\
    --genomeDir $stardir \\
    --genomeFastaFiles {params.reffa}

# MapSplice requires the {params.reffa} multifasta to be split into separate fastas
bash {params.script2} {params.reffa} {params.refdir}/separate_fastas
ls {params.refdir}/separate_fastas/*.fa | awk {AWK1} > {output.fastalst}
# may have to create bowtie1 index here.. has to be a separate rule ... see below
"""


TRSED = r"""tr '\n' ',' | sed 's/.$//g'"""


rule create_mapsplice_index:
    input:
        fastalst=rules.create_index.output.fastalst,
    output:
        rev1ebwt=join(REF_DIR, "separate_fastas_index.rev.1.ebwt"),
    params:
        separate_fastas=join(REF_DIR, "separate_fastas"),
        ebwt=join(REF_DIR, "separate_fastas_index"),
    threads: getthreads("create_mapsplice_index")
    container: "docker://cgrlab/mapsplice2:latest"
    shell:
        """
set -exo pipefail
fastas=$(ls {params.separate_fastas}/*.fa| {TRSED})
/opt/MapSplice2/bin/bowtie-build \
$fastas \
{params.ebwt}
"""

rule create_bwa_index:
    input:
        # FASTAS_REGIONS_GTFS
        list(map(lambda x: ancient(x), FASTAS_REGIONS_GTFS)),
    output:
        bwt=join(REF_DIR,"ref.bwt"),
        log=join(REF_DIR,"bwa_index.log")
    params:
        reffa=REF_FA
    container: config['containers']["base"]
    shell:"""
set -exo pipefail
refdir=$(dirname {params.reffa})
cd $refdir
bwa index -p ref {params.reffa} > bwa_index.log
"""

rule create_bowtie2_index:
    input:
        # FASTAS_REGIONS_GTFS
        list(map(lambda x: ancient(x), FASTAS_REGIONS_GTFS)),
    output:
        bt2=join(REF_DIR,"ref.1.bt2")
    params:
        reffa=REF_FA
    container: config['containers']["base"]
    shell:"""
set -exo pipefail
refdir=$(dirname {params.reffa})
cd $refdir
bowtie2-build {params.reffa} ref
"""

rule create_bowtie1_index:
    input:
        # FASTAS_REGIONS_GTFS
        list(map(lambda x: ancient(x), FASTAS_REGIONS_GTFS)),
    output:
        bt2=join(REF_DIR,"ref.1.ebwt")
    params:
        reffa=REF_FA
    container: config['containers']["bowtie1"]
    shell:"""
set -exo pipefail
refdir=$(dirname {params.reffa})
cd $refdir
bowtie-build {params.reffa} ref
"""
