AWK1 = r"""-F"/" '{print $NF}'"""


rule create_index:
    input:
        # FASTAS_REGIONS_GTFS
        list(map(lambda x: ancient(x), FASTAS_REGIONS_GTFS)),
    output:
        genepred_w_geneid=join(REF_DIR, "ref.genes.genepred_w_geneid"),
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
        randomstr=str(uuid.uuid4()),
        nclscan_dir=config["nclscan_dir"],
        nclscan_docker_dir=config["nclscan_docker_dir"],
        nclscan_config=config["nclscan_config"],
        nclscan_docker_config=config["nclscan_docker_config"],
    envmodules:
        TOOLS["bwa"]["version"],
        TOOLS["samtools"]["version"],
        TOOLS["ucsc"]["version"],
        TOOLS["cufflinks"]["version"],
    container: config['containers']['star']
    threads: getthreads("create_index")
    shell:
        """
set -exo pipefail
cd {params.refdir}
if [[ -f "/.dockerenv" || -f /singularity ]];then
    # you are inside a docker or singularity container
    ncldir="{params.nclscan_docker_dir}"
    nclconf="{params.nclscan_docker_config}"
else
    ncldir="{params.nclscan_dir}"
    nclconf="{params.nclscan_config}"
fi
samtools faidx {params.reffa} && \
    cut -f1-2 {params.reffa}.fai > {params.reffa}.sizes

# NCLscan files
python {params.script3} --ingtf {params.refgtf} --outgtf {output.fixed_gtf}
gffread -w {output.transcripts_fa} -g {params.reffa} {output.fixed_gtf}
touch {output.lncRNA_transcripts_fa}
${{ncldir}}/bin/create_reference.py -c ${{nclconf}}
# the above script internally runs bwa and that takes a while!

gtfToGenePred -ignoreGroupsWithoutExons {output.fixed_gtf} ref.genes.genepred && \
    python {params.script1} {output.fixed_gtf} ref.genes.genepred > {output.genepred_w_geneid}

# MapSplice requires the {params.reffa} multifasta to be split into separate fastas
bash {params.script2} {params.reffa} {params.refdir}/separate_fastas
ls {params.refdir}/separate_fastas/*.fa | awk {AWK1} > {output.fastalst}
# may have to create bowtie1 index here.. has to be a separate rule ... see below
"""


TRSED = r"""tr '\n' ',' | sed 's/.$//g'"""

rule create_star_index:
    input:
        # FASTAS_REGIONS_GTFS
        list(map(lambda x: ancient(x), FASTAS_REGIONS_GTFS)),
    output:
        sa=join(REF_DIR, "STAR_no_GTF", "SA"),
    params:
        reffa=REF_FA
    envmodules: TOOLS["star"]["version"]
    container: config['containers']['star']
    threads: getthreads("create_star_index")
    shell:
        """
set -exo pipefail
stardir=$(dirname {output.sa})
parentdir=$(dirname ${{stardir}}) && cd $parentdir
mkdir -p $stardir && \\
STAR \\
    --runThreadN {threads} \\
    --runMode genomeGenerate \\
    --genomeDir $stardir \\
    --genomeFastaFiles {params.reffa}
"""

rule create_mapsplice_index:
    input:
        fastalst=rules.create_index.output.fastalst,
    output:
        rev1ebwt=join(REF_DIR, "separate_fastas_index.rev.1.ebwt"),
    params:
        separate_fastas=join(REF_DIR, "separate_fastas"),
        ebwt=join(REF_DIR, "separate_fastas_index"),
    threads: getthreads("create_mapsplice_index")
    container: config['containers']['mapsplice']
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
    params:
        reffa=REF_FA
    envmodules: TOOLS["bwa"]["version"]
    container: config['containers']['star']
    shell:"""
set -exo pipefail
refdir=$(dirname {params.reffa})
cd $refdir
# increasing -b from default ... based off of https://github.com/lh3/bwa/issues/104
bwa index -b 10000000000 -p ref {params.reffa}
"""

rule create_bowtie2_index:
    input:
        # FASTAS_REGIONS_GTFS
        list(map(lambda x: ancient(x), FASTAS_REGIONS_GTFS)),
    output:
        bt2=join(REF_DIR,"ref.1.bt2")
    params:
        reffa=REF_FA
    envmodules: TOOLS["bowtie2"]["version"]
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
    envmodules: TOOLS["bowtie1"]["version"]
    shell:"""
set -exo pipefail
refdir=$(dirname {params.reffa})
cd $refdir
bowtie-build {params.reffa} ref
"""  