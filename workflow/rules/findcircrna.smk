# Find circRNAs using:
# ------------------------------------------------------------------------------------------------
# TOOL          |    Main output file
# ------------------------------------------------------------------------------------------------
# circExplorer2 | "results","{sample}","circExplorer","{sample}.circExplorer.counts_table.tsv"
# ciri2         | "results","{sample}","ciri","{sample}.ciri.out"
# CLEAR         | "results","{sample}","CLEAR","quant.txt.annotated" ... not used as this is just filtered circExplorer output
# DCC           | "results","{sample}","DCC","{sample}.dcc.counts_table.tsv"
# MapSplice     | "results","{sample}","MapSplice","circular_RNAs.txt"
# NCLscan       | "results","{sample}","NCLscan","{sample}.nclscan.counts_table.tsv"
# and annotate them


## function
# def get_dcc_inputs(wildcards):
#     filelist=[]
#     for s in SAMPLES:
#         filelist.append(join(WORKDIR,"results",s,"STAR1p",s+"_p1.Chimeric.out.junction"))
#         filelist.append(join(WORKDIR,"results",s,"STAR1p","mate1",s+"_mate1.Chimeric.out.junction"))
#         filelist.append(join(WORKDIR,"results",s,"STAR1p","mate2",s+"_mate2.Chimeric.out.junction"))
#     return filelist

def get_nclscan_target_files_per_sample(wildcards):
    targetfiles=dict()
    s=wildcards.sample
    if SAMPLESDF.loc[[s],"PEorSE"][0]=="PE": # SE is already take care of by function get_nclscan_target_files
        targetfiles['fixed_gtf']=join(REF_DIR,"ref.fixed.gtf")
        targetfiles['ndx']=join(REF_DIR,"NCLscan_index","AllRef.ndx")
        targetfiles['R1']=join(WORKDIR,"results",s,"trim",s+".R1.trim.fastq.gz")
        targetfiles['R2']=join(WORKDIR,"results",s,"trim",s+".R2.trim.fastq.gz")
    return targetfiles # empty if SE and will not run the rule at all!

def get_per_sample_files_to_merge(wildcards):
    filedict={}
    s=wildcards.sample
    filedict['circExplorer']=join(WORKDIR,"results",s,"circExplorer",s+".circExplorer.counts_table.tsv")
    filedict['CIRI']=join(WORKDIR,"results",s,"ciri",s+".ciri.out.filtered")
    # # if RUN_CLEAR:
    # #     filedict['CLEAR']=join(WORKDIR,"results","{sample}","CLEAR","quant.txt.annotated")
    if RUN_DCC:
        filedict['DCC']=join(WORKDIR,"results","{sample}","DCC","{sample}.dcc.counts_table.tsv.filtered")
    if RUN_MAPSPLICE:
        filedict['MapSplice']=join(WORKDIR,"results","{sample}","MapSplice","{sample}.mapsplice.counts_table.tsv.filtered")
    if RUN_NCLSCAN:
        filedict['NCLscan']=join(WORKDIR,"results","{sample}","NCLscan","{sample}.nclscan.counts_table.tsv.filtered")
    return(filedict)


## rules
# rule circExplorer:
# find circRNAs using circExplorer2 and then annotate it with known gene-annotations
# "annotate" requires GENCODE genes in the following columns:
# | #  |    ColName  | Description                  |
# |----|-------------|------------------------------|
# | 1  | geneName    | Name of gene                 |
# | 2  | isoformName | Name of isoform              |
# | 3  | chrom       | Reference sequence           |
# | 4  | strand      | + or - for strand            |
# | 5  | txStart     | Transcription start position |
# | 6  | txEnd       | Transcription end position   |
# | 7  | cdsStart    | Coding region start          |
# | 8  | cdsEnd      | Coding region end            |
# | 9  | exonCount   | Number of exons              |
# | 10 | exonStarts  | Exon start positions         |
# | 11 | exonEnds    | Exon end positions           |

# outout "known" file has the following columns:
# | #  | ColName     |  Description                        |
# |----|-------------|-------------------------------------|
# | 1  | chrom       | Chromosome                          |
# | 2  | start       | Start of circular RNA               |
# | 3  | end         | End of circular RNA                 |
# | 4  | name        | Circular RNA/Junction reads         |
# | 5  | score       | Flag of fusion junction realignment |
# | 6  | strand      | + or - for strand                   |
# | 7  | thickStart  | No meaning                          |
# | 8  | thickEnd    | No meaning                          |
# | 9  | itemRgb     | 0,0,0                               |
# | 10 | exonCount   | Number of exons                     |
# | 11 | exonSizes   | Exon sizes                          |
# | 12 | exonOffsets | Exon offsets                        |
# | 13 | readNumber  | Number of junction reads            |
# | 14 | circType    | Type of circular RNA                |
# | 15 | geneName    | Name of gene                        |
# | 16 | isoformName | Name of isoform                     |
# | 17 | index       | Index of exon or intron             |
# | 18 | flankIntron | Left intron/Right intron            |

# output low confidence file columns
# | # | ColName   |  Description                        |
# |---|-----------|-------------------------------------|
# | 1 | chrom     | Chromosome                          |
# | 2 | start     | Start of circular RNA               |
# | 3 | end       | End of circular RNA                 |
# | 4 | name      | Circular RNA/Junction reads         |
# | 5 | score     | Flag of fusion junction realignment |
# | 6 | strand    | + or - for strand                   |
# | 7 | leftInfo  | Gene:Isoform:Index of left exon     |
# | 8 | rightInfo | Gene:Isoform:Index of right exon    |

# STEPS:
# 1. parse the chimeric junctions file from STAR to CircExplorer2 'parse' to generate the back_spliced_junction BED file
# 2. parse the back_spliced_junction BED from above along with known splicing annotations to CircExplorer2 'parse' to create
#       a. circularRNA_known.txt ... circRNAs around known gene exons
#       b. low_conf_circularRNA_known.txt .... circRNAs with low confidence
# 3. parse back_spliced_junction BED along with circularRNA_known.txt and low_conf_circularRNA_known.txt to custom python script
# to create an aggregated list of BSJs with following columns:
# | # | ColName     |
# |---|-------------|
# | 1 | chrom       |
# | 2 | start       |
# | 3 | end         |
# | 4 | strand      |
# | 5 | read_count  |
# | 6 | known_novel |
# known_novel can have 3 different values:
# a. known ... BSJ is around a known gene-exon
# b. novel ... BSJ is not around a known gene-exon ... it is absent in circularRNA_known.txt or low_conf_circularRNA_known.txt files 
#       but present in ack_spliced_junction BED
# c. low_conf ... BSJ is around a known gene-exon but circExplorer called it with low-confidence.
# ref: https://circexplorer2.readthedocs.io/en/latest/
rule circExplorer:
    input:
        junctionfile=rules.star2p.output.junction
    output:
        backsplicedjunctions=join(WORKDIR,"results","{sample}","circExplorer","{sample}.back_spliced_junction.bed"),
        annotations=join(WORKDIR,"results","{sample}","circExplorer","{sample}.circularRNA_known.txt"),
        annotation_counts_table=join(WORKDIR,"results","{sample}","circExplorer","{sample}.circExplorer.annotation_counts.tsv")
    params:
        sample="{sample}",
        bsj_min_nreads=config['circexplorer_bsj_circRNA_min_reads'], # in addition to "known" and "low-conf" circRNAs identified by circexplorer, we also include those found in back_spliced.bed file but not classified as known/low-conf only if the number of reads supporting the BSJ call is greater than this number
        outdir=join(WORKDIR,"results","{sample}","circExplorer"),
        genepred=rules.create_index.output.genepred_w_geneid,
        reffa=REF_FA,
        refregions=REF_REGIONS,
        host=HOST,
        additives=ADDITIVES,
        viruses=VIRUSES,
        minsize_host=config["minsize_host"],
        maxsize_host=config["maxsize_host"],
        minsize_virus=config["minsize_virus"],
        maxsize_virus=config["maxsize_virus"],
        script=join(SCRIPTS_DIR,"circExplorer_get_annotated_counts_per_sample.py") # this produces an annotated counts table to which counts found in BAMs need to be appended
    threads: getthreads("circExplorer")
    envmodules: TOOLS["circexplorer"]["version"]
    shell:"""
set -exo pipefail
if [ ! -d {params.outdir} ];then mkdir {params.outdir};fi
cd {params.outdir}
mv {input.junctionfile} {input.junctionfile}.original
grep -v junction_type {input.junctionfile}.original > {input.junctionfile}
CIRCexplorer2 parse \\
    -t STAR \\
    {input.junctionfile} > {params.sample}_circexplorer_parse.log 2>&1
mv back_spliced_junction.bed {output.backsplicedjunctions}
mv {input.junctionfile}.original {input.junctionfile}
CIRCexplorer2 annotate \\
-r {params.genepred} \\
-g {params.reffa} \\
-b {output.backsplicedjunctions} \\
-o $(basename {output.annotations}) \\
--low-confidence

python {params.script} \\
    --back_spliced_bed {output.backsplicedjunctions} \\
    --back_spliced_min_reads {params.bsj_min_nreads} \\
    --circularRNA_known {output.annotations} \\
    --low_conf low_conf_$(basename {output.annotations}) \\
    --host {params.host} \\
    --additives {params.additives} \\
    --viruses {params.viruses} \\
    --regions {params.refregions} \\
    --host_filter_min {params.minsize_host} \\
    --host_filter_max {params.maxsize_host} \\
    --virus_filter_min {params.minsize_virus} \\
    --virus_filter_max {params.maxsize_virus} \\
    -o {output.annotation_counts_table}
"""


# rule ciri:
# call circRNAs using CIRI2. The output file has following columns:
# | #  | colName              | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
# |----|----------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
# | 1  | circRNA_ID           | ID of a predicted circRNA in the pattern of "chr:start|end"                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
# | 2  | chr                  | chromosome of a predicted circRNA                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
# | 3  | circRNA_start        | start loci of a predicted circRNA on the chromosome                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
# | 4  | circRNA_end          | end loci of a predicted circRNA on the chromosome                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
# | 5  | #junction_reads      | circular junction read (also called as back-spliced junction read) count of a predicted circRNA                                                                                                                                                                                                                                                                                                                                                                                                         |
# | 6  | SM_MS_SMS            | unique CIGAR types of a predicted circRNA. For example, a circRNAs have three junction reads: read A (80M20S, 80S20M), read B (80M20S, 80S20M), read C (40M60S, 40S30M30S, 70S30M), then its has two SM types (80S20M, 70S30M), two MS types (80M20S, 70M30S) and one SMS type (40S30M30S). Thus its SM_MS_SMS should be 2_2_1.                                                                                                                                                                         |
# | 7  | #non_junction_reads  | non-junction read count of a predicted circRNA that mapped across the circular junction but consistent with linear RNA instead of being back-spliced                                                                                                                                                                                                                                                                                                                                                    |
# | 8  | junction_reads_ratio | ratio of circular junction reads calculated by 2 * #junction_reads/(2 * #junction_reads+#non_junction_reads). #junction_reads is multiplied by two because a junction read is generated from two ends of circular junction but only counted once while a non-junction read is from one end. It has to be mentioned that the non-junction reads are still possibly from another larger circRNA, so the junction_reads_ratio based on it may be an inaccurate estimation of relative expression of the circRNA. |
# | 9  | circRNA_type         | type of a circRNA according to positions of its two ends on chromosome (exon, intron or intergenic_region; only available when annotation file is provided)                                                                                                                                                                                                                                                                                                                                             |
# | 10 | gene_id              | ID of the gene(s) where an exonic or intronic circRNA locates                                                                                                                                                                                                                                                                                                                                                                                                                                           |
# | 11 | strand               | strand info of a predicted circRNAs (new in CIRI2)                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
# | 12 | junction_reads_ID    | all of the circular junction read IDs (split by ",")                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
# ref: https://ciri-cookbook.readthedocs.io/en/latest/CIRI2.html#an-example-of-running-ciri2
rule ciri:
    input:
        bwt=rules.create_index.output.bwt,
        R1=rules.cutadapt.output.of1,
        R2=rules.cutadapt.output.of2,
        gtf=rules.create_index.output.fixed_gtf,
    output:
        cirilog=join(WORKDIR,"results","{sample}","ciri","{sample}.ciri.log"),
        bwalog=join(WORKDIR,"results","{sample}","ciri","{sample}.bwa.log"),
        ciribam=join(WORKDIR,"results","{sample}","ciri","{sample}.ciri.cram"),
        ciriout=join(WORKDIR,"results","{sample}","ciri","{sample}.ciri.out"),
        cirioutfiltered=join(WORKDIR,"results","{sample}","ciri","{sample}.ciri.out.filtered"),
    params:
        sample="{sample}",
        memG=getmemG("ciri"),
        outdir=join(WORKDIR,"results","{sample}","ciri"),
        peorse=get_peorse,
        genepred=rules.create_index.output.genepred_w_geneid,
        reffa=REF_FA,
        bwaindex=BWA_INDEX,
        ciripl=config['ciri_perl_script'],
        bsj_min_nreads=config['minreadcount'],
        refregions=REF_REGIONS,
        host=HOST,
        additives=ADDITIVES,
        viruses=VIRUSES,
        minsize_host=config["minsize_host"],
        maxsize_host=config["maxsize_host"],
        minsize_virus=config["minsize_virus"],
        maxsize_virus=config["maxsize_virus"],
        script=join(SCRIPTS_DIR,"filter_ciriout.py"),
        randomstr=str(uuid.uuid4()),
    threads: getthreads("ciri")
    envmodules: TOOLS["bwa"]["version"], TOOLS["samtools"]["version"]
    shell:"""
set -exo pipefail
if [ -d /lscratch/${{SLURM_JOB_ID}} ];then
    TMPDIR="/lscratch/${{SLURM_JOB_ID}}/{params.randomstr}"
else
    TMPDIR="/dev/shm/{params.randomstr}"
fi
if [ ! -d $TMPDIR ];then mkdir -p $TMPDIR;fi
cd {params.outdir}
if [ "{params.peorse}" == "PE" ];then
    ## paired-end
    bwa mem -t {threads} -T 19 \\
    {params.bwaindex} \\
    {input.R1} {input.R2} \\
    > {params.sample}.bwa.sam 2> {output.bwalog}
else
    ## single-end
    bwa mem -t {threads} -T 19 \\
    {params.bwaindex} \\
    {input.R1} \\
    > {params.sample}.bwa.sam 2> {output.bwalog}
fi
perl {params.ciripl} \\
-I {params.sample}.bwa.sam \\
-O {output.ciriout} \\
-F {params.reffa} \\
-A {input.gtf} \\
-G {output.cirilog} -T {threads}
samtools view -@{threads} -T {params.reffa} -CS {params.sample}.bwa.sam | samtools sort -l 9 -T $TMPDIR --write-index -@{threads} -O CRAM -o {output.ciribam} -
rm -rf {params.sample}.bwa.sam
python {params.script} \\
    --ciriout {output.ciriout} \\
    --back_spliced_min_reads {params.bsj_min_nreads} \\
    --host {params.host} \\
    --additives {params.additives} \\
    --viruses {params.viruses} \\
    --regions {params.refregions} \\
    --host_filter_min {params.minsize_host} \\
    --host_filter_max {params.maxsize_host} \\
    --virus_filter_min {params.minsize_virus} \\
    --virus_filter_max {params.maxsize_virus} \\
    -o {output.cirioutfiltered}
"""



rule create_ciri_count_matrix:
# DEPRECATED
    input:
        expand(join(WORKDIR,"results","{sample}","ciri","{sample}.ciri.out"),sample=SAMPLES)
    output:
        matrix=join(WORKDIR,"results","ciri_count_matrix.txt")
    params:
        script=join(SCRIPTS_DIR,"Create_ciri_count_matrix.py"),
        lookup=ANNOTATION_LOOKUP,
        outdir=join(WORKDIR,"results"),
        hostID=HOST+"ID"
    envmodules: TOOLS["python37"]["version"]
    shell:"""
set -exo pipefail
cd {params.outdir}
python {params.script} {params.lookup} {params.hostID}
"""

rule create_circexplorer_count_matrix:
# DEPRECATED
    input:
        expand(join(WORKDIR,"results","{sample}","circExplorer","{sample}.circularRNA_known.txt"),sample=SAMPLES)
    output:
        matrix=join(WORKDIR,"results","circExplorer_count_matrix.txt"),
        matrix2=join(WORKDIR,"results","circExplorer_BSJ_count_matrix.txt")
    params:
        script=join(SCRIPTS_DIR,"Create_circExplorer_count_matrix.py"),
        script2=join(SCRIPTS_DIR,"Create_circExplorer_BSJ_count_matrix.py"),
        lookup=ANNOTATION_LOOKUP,
        outdir=join(WORKDIR,"results"),
        hostID=HOST+"ID"
    envmodules: TOOLS["python37"]["version"]
    shell:"""
cd {params.outdir}
python {params.script} {params.lookup} {params.hostID}
python {params.script2} {params.lookup} {params.hostID}
"""

# rule clear:
# quantify circRNAs using CLEAR
# This uses the "known" circRNAs from circExplorer2
# Hence, not considered as a separate method for detecting circRNAs
# CLEAR (aka. CircExplorer3) is run for completeness of the circExplorer pipeline
# and to extract "Relative expression of circRNA" for downstream purposes
# CLEAR does not quantify "Relative expression of circRNA" for novel circRNA, ie., 
# circRNAs not labeled as "known" possible due to poor genome annotation.
# circRNA is labled as "known" if its coordinates match with exons of known genes! 
# quant.txt is a TSV with the following columns:
# | #  | ColName     | Description                         |
# |----|-------------|-------------------------------------|
# | 1  | chrom       | Chromosome                          |
# | 2  | start       | Start of circular RNA               |
# | 3  | end         | End of circular RNA                 |
# | 4  | name        | Circular RNA/Junction reads         |
# | 5  | score       | Flag of fusion junction realignment |
# | 6  | strand      | + or - for strand                   |
# | 7  | thickStart  | No meaning                          |
# | 8  | thickEnd    | No meaning                          |
# | 9  | itemRgb     | 0,0,0                               |
# | 10 | exonCount   | Number of exons                     |
# | 11 | exonSizes   | Exon sizes                          |
# | 12 | exonOffsets | Exon offsets                        |
# | 13 | readNumber  | Number of junction reads            |
# | 14 | circType    | Type of circular RNA                |
# | 15 | geneName    | Name of gene                        |
# | 16 | isoformName | Name of isoform                     |
# | 17 | index       | Index of exon or intron             |
# | 18 | flankIntron | Left intron/Right intron            |
# | 19 | FPBcirc     | Expression of circRNA               |
# | 20 | FPBlinear   | Expression of cognate linear RNA    |
# | 21 | CIRCscore   | Relative expression of circRNA      |
rule clear:
    input:
        bam=rules.star2p.output.bam,
        circexplorerout=rules.circExplorer.output.annotations,
    output:
        quantfile=join(WORKDIR,"results","{sample}","CLEAR","quant.txt")
    params:
        genepred=rules.create_index.output.genepred_w_geneid,
    container: "docker://nciccbr/ccbr_clear:latest"
    threads: getthreads("clear")
    shell:"""
set -exo pipefail
circ_quant \\
-c {input.circexplorerout} \\
-b {input.bam} \\
-t \\
-r {params.genepred} \\
-o {output.quantfile}
"""

# rule annotate_clear_output:
# annotate CLEAR output with circRNA databases
# the ".annotated" file columns are:
# | #  | ColName            |
# |----|--------------------|
# | 1  | hg38ID             |
# | 2  | quant_chrom        |
# | 3  | quant_start        |
# | 4  | quant_end          |
# | 5  | quant_name         |
# | 6  | quant_score        |
# | 7  | quant_quant_strand |
# | 8  | quant_thickStart   |
# | 9  | quant_thickEnd     |
# | 10 | quant_itemRgb      |
# | 11 | quant_exonCount    |
# | 12 | quant_exonSizes    |
# | 13 | quant_exonOffsets  |
# | 14 | quant_readNumber   |
# | 15 | quant_circType     |
# | 16 | quant_geneName     |
# | 17 | quant_isoformName  |
# | 18 | quant_index        |
# | 19 | quant_flankIntron  |
# | 20 | quant_FPBcirc      |
# | 21 | quant_FPBlinear    |
# | 22 | quant_CIRCscore    |
# | 23 | hg19ID             |
# | 24 | strand             |
# | 25 | circRNA.ID         |
# | 26 | genomic.length     |
# | 27 | spliced.seq.length |
# | 28 | samples            |
# | 29 | repeats            |
# | 30 | annotation         |
# | 31 | best.transcript    |
# | 32 | gene.symbol        |
# | 33 | circRNA.study      |
localrules: annotate_clear_output
rule annotate_clear_output:
    input:
        quantfile=rules.clear.output.quantfile
    output:
        annotatedquantfile=join(WORKDIR,"results","{sample}","CLEAR","quant.txt.annotated")
    params:
        script=join(SCRIPTS_DIR,"annotate_clear_quant.py"),
        lookup=ANNOTATION_LOOKUP,
        cleardir=join(WORKDIR,"results","{sample}","CLEAR"),
        hostID=HOST+"ID"
    shell:"""
set -exo pipefail
## cleanup quant.txt* dirs before annotation
find {params.cleardir} -maxdepth 1 -type d -name "quant.txt*" -exec rm -rf {{}} \;
if [[ "$(cat {input.quantfile} | wc -l)" != "0" ]]
then
python {params.script} {params.lookup} {input.quantfile} {params.hostID}
else
touch {output.annotatedquantfile}
fi
"""		

localrules: dcc_create_samplesheets
rule dcc_create_samplesheets:
    input:
        f1=join(WORKDIR,"results","{sample}","STAR1p","{sample}"+"_p1.Chimeric.out.junction"),
        f2=join(WORKDIR,"results","{sample}","STAR1p","mate1","{sample}"+"_mate1.Chimeric.out.junction"),
        f3=join(WORKDIR,"results","{sample}","STAR1p","mate2","{sample}"+"_mate2.Chimeric.out.junction"),
    output:
        ss=join(WORKDIR,"results","{sample}","DCC","samplesheet.txt"),
        m1=join(WORKDIR,"results","{sample}","DCC","mate1.txt"),
        m2=join(WORKDIR,"results","{sample}","DCC","mate2.txt"),
    shell:"""
set -exo pipefail
outdir=$(dirname {output.ss})
if [ ! -d $outdir ];then mkdir -p $outdir;fi
echo "{input.f1}" > {output.ss}
echo "{input.f2}" > {output.m1}
echo "{input.f3}" > {output.m2}
"""

# rule dcc:
# output files
# CircRNACount columns are:
# | # | ColName                        |
# |---|--------------------------------|
# | 1 | Chr                            |
# | 2 | Start                          |
# | 3 | End                            |
# | 4 | Strand                         |
# | 5 | <sample_junction_filename>     | <-- circRNA read counts for this sample
#
# CircCoordinates columns are:
# | # | ColName          | Description                                                                  |
# |---|------------------|------------------------------------------------------------------------------|
# | 1 | Chr              |                                                                              |
# | 2 | Start            |                                                                              |
# | 3 | End              |                                                                              |
# | 4 | Gene             |                                                                              |
# | 5 | JunctionType     | 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT |
# | 6 | Strand           |                                                                              |
# | 7 | Start-End Region | eg. intron-intergenic, exon-exon, intergenic-intron, etc.                    |
# | 8 | OverallRegion    | the genomic features circRNA coordinates interval covers                     |

# output dcc.counts_table.tsv has the following columns:
# | # | ColName        |
# |---|----------------|
# | 1 | chr            |
# | 2 | start          |
# | 3 | end            |
# | 4 | strand         |
# | 5 | read_count     |
# | 6 | dcc_annotation | --> this is JunctionType##Start-End Region from CircCoordinates file
rule dcc:
    input:
        ss=rules.dcc_create_samplesheets.output.ss,
        m1=rules.dcc_create_samplesheets.output.m1,
        m2=rules.dcc_create_samplesheets.output.m2,
        bam=rules.star2p.output.bam,
        gtf=rules.create_index.output.fixed_gtf,
    output:
        cr=join(WORKDIR,"results","{sample}","DCC","CircRNACount"),
        cc=join(WORKDIR,"results","{sample}","DCC","CircCoordinates"),
        linear=join(WORKDIR,"results","{sample}","DCC","LinearCount"),
        ct=join(WORKDIR,"results","{sample}","DCC","{sample}.dcc.counts_table.tsv"),
        ctf=join(WORKDIR,"results","{sample}","DCC","{sample}.dcc.counts_table.tsv.filtered"),
    threads: getthreads("dcc")
    envmodules: TOOLS["python27"]["version"]
    params:
        peorse=get_peorse,
        dcc_strandedness=config['dcc_strandedness'],
        rep=REPEATS_GTF,
        fa=REF_FA,
        randomstr=str(uuid.uuid4()),
        script=join(SCRIPTS_DIR,"create_dcc_per_sample_counts_table.py"),
        bsj_min_nreads=config['minreadcount'],
        refregions=REF_REGIONS,
        host=HOST,
        additives=ADDITIVES,
        viruses=VIRUSES,
        minsize_host=config["minsize_host"],
        maxsize_host=config["maxsize_host"],
        minsize_virus=config["minsize_virus"],
        maxsize_virus=config["maxsize_virus"],
        script2=join(SCRIPTS_DIR,"filter_dcc.py")
    shell:"""
set -exo pipefail
if [ -d /lscratch/${{SLURM_JOB_ID}} ];then
    TMPDIR="/lscratch/${{SLURM_JOB_ID}}/{params.randomstr}"
else
    TMPDIR="/dev/shm/{params.randomstr}"
fi
if [ ! -d $TMPDIR ];then mkdir -p $TMPDIR;fi

. "/data/CCBR_Pipeliner/db/PipeDB/Conda/etc/profile.d/conda.sh"
conda activate DCC
cd $(dirname {output.cr})
if [ "{params.peorse}" == "PE" ];then
DCC @{input.ss} \\
    --temp ${{TMPDIR}}/DCC \\
    --threads {threads} \\
    --detect --gene \\
    --bam {input.bam} \\
    {params.dcc_strandedness} \\
    --annotation {input.gtf} \\
    --chrM -G \\
    --rep_file {params.rep} \\
    --refseq {params.fa} \\
    --PE-independent \\
    -mt1 @{input.m1} \\
    -mt2 @{input.m2}
else
DCC @{input.ss} \\
    --temp ${{TMPDIR}}/DCC \\
    --threads {threads} \\
    --detect --gene \\
    --bam {input.bam} \\
    {params.dcc_strandedness} \\
    --annotation {input.gtf} \\
    --chrM -G \\
    --rep_file {params.rep} \\
    --refseq {params.fa} 
fi

ls -alrth ${{TMPDIR}}

paste {output.cr} {output.linear} | cut -f1-5,9 > ${{TMPDIR}}/CircRNALinearCount

python {params.script} \\
  --CircCoordinates {output.cc} --CircRNALinearCount ${{TMPDIR}}/CircRNALinearCount -o {output.ct}

python {params.script2} \\
    --in_dcc_counts_table {output.ct} \\
    --out_dcc_filtered_counts_table {output.ctf} \\
    --back_spliced_min_reads {params.bsj_min_nreads} \\
    --host {params.host} \\
    --additives {params.additives} \\
    --viruses {params.viruses} \\
    --regions {params.refregions} \\
    --host_filter_min {params.minsize_host} \\
    --host_filter_max {params.maxsize_host} \\
    --virus_filter_min {params.minsize_virus} \\
    --virus_filter_max {params.maxsize_virus}
"""


# rule mapsplice:
# output "circular_RNA.txt" has following columns
# ref: https://github.com/Aufiero/circRNAprofiler/blob/master/R/importFilesPredictionTool.R
# | #  | ColName                         | Example                   |
# |----|---------------------------------|---------------------------|
# | 1  | chrom                           | chr1~chr1                 |
# | 2  | doner_end                       | 1223244                   |
# | 3  | acceptor_start                  | 1223968                   |
# | 4  | id                              | FUSIONJUNC_427            |
# | 5  | coverage                        | 26                        |
# | 6  | strand                          | --                        |
# | 7  | rgb                             | 255,0,0                   |
# | 8  | block_count                     | 2                         |
# | 9  | block_size                      | 147,130,147,138,          |
# | 10 | block_distance                  | 0,855,                    |
# | 11 | entropy                         | 2.811419                  |
# | 12 | flank_case                      | 5                         |
# | 13 | flank_string                    | GTAG                      |
# | 14 | min_mismatch                    | 1                         |
# | 15 | max_mismatch                    | 2                         |
# | 16 | ave_mismatch                    | 1.038462                  |
# | 17 | max_min_suffix                  | 72                        |
# | 18 | max_min_prefix                  | 71                        |
# | 19 | min_anchor_difference           | 7                         |
# | 20 | unique_read_count               | 26                        |
# | 21 | multi_read_count                | 0                         |
# | 22 | paired_read_count               | 6                         |
# | 23 | left_paired_read_count          | 3                         |
# | 24 | right_paired_read_count         | 3                         |
# | 25 | multiple_paired_read_count      | 0                         |
# | 26 | unique_paired_read_count        | 6                         |
# | 27 | single_read_count               | 20                        |
# | 28 | encompassing_read               | 0                         |
# | 29 | doner_start                     | 1223391                   |
# | 30 | acceptor_end                    | 1223838                   |
# | 31 | doner_iosforms                  | 1223244,114M474N113M|     |
# | 32 | acceptor_isoforms               | 1223246,112M474N137M|     |
# | 33 | obsolete1                       | 0                         |
# | 34 | obsolete2                       | 0                         |
# | 35 | obsolete3                       | 0.404444                  |
# | 36 | obsolete4                       | 0.455556                  |
# | 37 | minimal_doner_isoform_length    | 227                       |
# | 38 | maximal_doner_isoform_length    | 227                       |
# | 39 | minimal_acceptor_isoform_length | 249                       |
# | 40 | maximal_acceptor_isoform_length | 249                       |
# | 41 | paired_reads_entropy            | 1.56071                   |
# | 42 | mismatch_per_bp                 | 0.00687723                |
# | 43 | anchor_score                    | 1                         |
# | 44 | max_doner_fragment              | 227                       |
# | 45 | max_acceptor_fragment           | 249                       |
# | 46 | max_cur_fragment                | 396                       |
# | 47 | min_cur_fragment                | 323                       |
# | 48 | ave_cur_fragment                | 348.167                   |
# | 49 | doner_encompass_unique          | 0                         |
# | 50 | doner_encompass_multiple        | 0                         |
# | 51 | acceptor_encompass_unique       | 0                         |
# | 52 | acceptor_encompass_multiple     | 0                         |
# | 53 | doner_match_to_normal           | doner_exact_matched       |
# | 54 | acceptor_match_to_normal        | acceptor_exact_matched    |
# | 55 | doner_seq                       | GAGGAACTCAAAGTGGATGAGGAAA |
# | 56 | acceptor_seq                    | CTTCCGGTCAGTGTTCACATCCACC |
# | 57 | match_gene_strand               | 0                         |
# | 58 | annotated_type                  | from_fusion               |
# | 59 | fusion_type                     | normal                    |
# | 60 | gene_strand                     | *                         |
# | 61 | annotated_gene_donor            | SDF4,                     |
# | 62 | annotated_gene_acceptor         | SDF4,                     |
rule mapsplice:
    input:
        bwt=rules.create_mapsplice_index.output.rev1ebwt, 
        R1=rules.cutadapt.output.of1,
        R2=rules.cutadapt.output.of2,
        gtf=rules.create_index.output.fixed_gtf,
    output:
        # rev1ebwt=join(REF_DIR,"separate_fastas_index.rev.1.ebwt"),
        sam=temp(join(WORKDIR,"results","{sample}","MapSplice","alignments.sam")),
        circRNAs=join(WORKDIR,"results","{sample}","MapSplice","circular_RNAs.txt"),
    params:
        peorse=get_peorse,
        minmaplen=MAPSPLICE_MIN_MAP_LEN,
        filtering=MAPSPLICE_FILTERING,
        separate_fastas=join(REF_DIR,"separate_fastas"),
        ebwt=join(REF_DIR,"separate_fastas_index"),
        outdir=join(WORKDIR,"results","{sample}","MapSplice"),
        randomstr=str(uuid.uuid4()),
    threads: getthreads("mapsplice")
    container: "docker://cgrlab/mapsplice2:latest"
    shell:"""
set -exo pipefail
if [ -d /lscratch/${{SLURM_JOB_ID}} ];then
    TMPDIR="/lscratch/${{SLURM_JOB_ID}}/{params.randomstr}"
else
    TMPDIR="/dev/shm/{params.randomstr}"
fi
if [ ! -d $TMPDIR ];then mkdir -p $TMPDIR;fi

MSHOME="/opt/MapSplice2"
# singularity exec -B /data/Ziegelbauer_lab,/data/kopardevn \
#     /data/kopardevn/SandBox/MapSplice/mapsplice2.sif \

if [ "{params.peorse}" == "PE" ];then

R1fn=$(basename {input.R1})
R2fn=$(basename {input.R2})
zcat {input.R1} > ${{TMPDIR}}/${{R1fn%.*}}
zcat {input.R2} > ${{TMPDIR}}/${{R2fn%.*}}

python $MSHOME/mapsplice.py \\
 -1 ${{TMPDIR}}/${{R1fn%.*}} \\
 -2 ${{TMPDIR}}/${{R2fn%.*}} \\
 -c {params.separate_fastas} \\
 -p {threads} \\
 --min-map-len {params.minmaplen} \\
 -x {params.ebwt} \\
 --non-canonical-double-anchor \\
 --non-canonical-single-anchor \\
 --filtering {params.filtering} \\
 --fusion-non-canonical --min-fusion-distance 200 \\
 --gene-gtf {input.gtf} \\
 -o {params.outdir}

else

R1fn=$(basename {input.R1})
zcat {input.R1} > ${{TMPDIR}}/${{R1fn%.*}}

python $MSHOME/mapsplice.py \
 -1 ${{TMPDIR}}/${{R1fn%.*}} \
 -c {params.separate_fastas} \
 -p {threads} \
 -x {params.ebwt} \
 --non-canonical-double-anchor \
 --non-canonical-single-anchor \
 --filtering 1 \
 --fusion-non-canonical --min-fusion-distance 200 \
 --gene-gtf {input.gtf} \
 -o {params.outdir}

fi
"""

# rule mapsplice_postprocess:
# the above file is filtered to only include the following columns:
# | # | ColName              | Eg.              |
# |---|----------------------|------------------|
# | 1 | chrom                | chr1             |
# | 2 | start                | 1223244          |
# | 3 | end                  | 1223968          |
# | 4 | strand               | -                |
# | 5 | read_count           | 26               |
# | 6 | mapsplice_annotation | normal##2.811419 | <--"fusion_type"##"entropy" 
# "fusion_type" is either "normal" or "overlapping" ... higher "entropy" values are better!
# mapsplice output contains an alignment.sam file which can be really large. Hence converting it to a sorted bam
# to save space
rule mapsplice_postprocess:
    input:
        sam=rules.mapsplice.output.sam,
        circRNAs=rules.mapsplice.output.circRNAs
    output:
        ct=join(WORKDIR,"results","{sample}","MapSplice","{sample}.mapsplice.counts_table.tsv"),
        ctf=join(WORKDIR,"results","{sample}","MapSplice","{sample}.mapsplice.counts_table.tsv.filtered"),
        bam=join(WORKDIR,"results","{sample}","MapSplice","{sample}.mapsplice.cram"),
        bai=join(WORKDIR,"results","{sample}","MapSplice","{sample}.mapsplice.cram.crai"),
    envmodules: TOOLS["samtools"]["version"], TOOLS["python27"]["version"]
    params:
        script=join(SCRIPTS_DIR,"create_mapsplice_per_sample_counts_table.py"),
        memG=getmemG("mapsplice_postprocess"),
        randomstr=str(uuid.uuid4()),
        bsj_min_nreads=config['minreadcount'],
        refregions=REF_REGIONS,
        reffa=REF_FA,
        host=HOST,
        additives=ADDITIVES,
        viruses=VIRUSES,
        minsize_host=config["minsize_host"],
        maxsize_host=config["maxsize_host"],
        minsize_virus=config["minsize_virus"],
        maxsize_virus=config["maxsize_virus"],
    threads: getthreads("mapsplice_postprocess")
    shell:"""
set -exo pipefail
if [ -d /lscratch/${{SLURM_JOB_ID}} ];then
    TMPDIR="/lscratch/${{SLURM_JOB_ID}}/{params.randomstr}"
else
    TMPDIR="/dev/shm/{params.randomstr}"
fi
if [ ! -d $TMPDIR ];then mkdir -p $TMPDIR;fi
python {params.script} \\
    --circularRNAstxt {input.circRNAs} \\
    -o {output.ct} \\
    -fo {output.ctf} \\
    --back_spliced_min_reads {params.bsj_min_nreads} \\
    --host {params.host} \\
    --additives {params.additives} \\
    --viruses {params.viruses} \\
    --regions {params.refregions} \\
    --host_filter_min {params.minsize_host} \\
    --host_filter_max {params.maxsize_host} \\
    --virus_filter_min {params.minsize_virus} \\
    --virus_filter_max {params.maxsize_virus}
cd $TMPDIR
samtools view -@{threads} -T {params.reffa} -CS {input.sam} | samtools sort -l 9 -T $TMPDIR --write-index -@{threads} -O CRAM -o alignments.cram -
rsync -az --progress alignments.cram {output.bam}
rsync -az --progress alignments.cram.crai {output.bai}
"""


# rule nclscan:
# result.txt output file has the following columns:
# ref: https://github.com/TreesLab/NCLscan
# | #  | Description                                | ColName
# |----|--------------------------------------------|------------ 
# | 1  | Chromosome name of the donor side (5'ss)   | chrd
# | 2  | Junction coordinate of the donor side      | coordd
# | 3  | Strand of the donor side                   | strandd
# | 4  | Chromosome name of the acceptor side (3'ss)| chra
# | 5  | Junction coordinate of the acceptor side   | coorda
# | 6  | Strand of the acceptor side                | stranda
# | 7  | Gene name of the donor side                | gened
# | 8  | Gene name of the acceptor side             | genea
# | 9  | Intragenic (1) or intergenic (0) case      | case
# | 10 | Total number of all supporting reads       | reads
# | 11 | Total number of junc-reads                 | jreads
# | 12 | Total number of span-reads                 | sreads
# the above file is filtered to only include the following columns:
# | # | ColName              | Eg.              |
# |---|----------------------|------------------|
# | 1 | chrom                | chr1             |
# | 2 | start                | 1223244          |
# | 3 | end                  | 1223968          |
# | 4 | strand               | -                |
# | 5 | read_count           | 26               |
# | 6 | nclscan_annotation   | 1                | <--1+1 for intragenic 0+1 for intergenic
rule nclscan:
    input:
        unpack(get_nclscan_target_files_per_sample)
    output:
        result=join(WORKDIR,"results","{sample}","NCLscan","{sample}.result"),
        ct=join(WORKDIR,"results","{sample}","NCLscan","{sample}.nclscan.counts_table.tsv"),
        ctf=join(WORKDIR,"results","{sample}","NCLscan","{sample}.nclscan.counts_table.tsv.filtered"),
    envmodules:
        TOOLS["ncl_required_modules"]
    threads: getthreads("nclscan")
    params:
        workdir=WORKDIR,
        sample="{sample}",
        peorse=get_peorse,
        nclscan_dir=config['nclscan_dir'],
        nclscan_config=config['nclscan_config'],
        script=join(SCRIPTS_DIR,"create_nclscan_per_sample_counts_table.py"),
        randomstr=str(uuid.uuid4()),
        bsj_min_nreads=config['minreadcount'],
        refregions=REF_REGIONS,
        host=HOST,
        additives=ADDITIVES,
        viruses=VIRUSES,
        minsize_host=config["minsize_host"],
        maxsize_host=config["maxsize_host"],
        minsize_virus=config["minsize_virus"],
        maxsize_virus=config["maxsize_virus"],
    shell:"""
set -exo pipefail
if [ -d /lscratch/${{SLURM_JOB_ID}} ];then
    TMPDIR="/lscratch/${{SLURM_JOB_ID}}/{params.randomstr}"
else
    TMPDIR="/dev/shm/{params.randomstr}"
fi
if [ ! -d $TMPDIR ];then mkdir -p $TMPDIR;fi
outdir=$(dirname {output.result})
results_bn=$(basename {output.result})

if [ "{params.peorse}" == "PE" ];then
{params.nclscan_dir}/NCLscan.py -c {params.nclscan_config} -pj {params.sample} -o $TMPDIR --fq1 {input.R1} --fq2 {input.R2}
rsync -az --progress ${{TMPDIR}}/${{results_bn}} {output.result}
python {params.script} \\
    --result {output.result} \\
    -o {output.ct} \\
    -fo {output.ctf} \\
    --back_spliced_min_reads {params.bsj_min_nreads} \\
    --host {params.host} \\
    --additives {params.additives} \\
    --viruses {params.viruses} \\
    --regions {params.refregions} \\
    --host_filter_min {params.minsize_host} \\
    --host_filter_max {params.maxsize_host} \\
    --virus_filter_min {params.minsize_virus} \\
    --virus_filter_max {params.maxsize_virus}
# else
#     outdir=$(dirname {output.result})
#     if [ ! -d $outdir ];then
#         mkdir -p $outdir
#     fi
#     touch {output.result}
#     touch {output.ct}
# This part is redundant as it is already taken care of by get_nclscan_target_files function!
fi
"""

rule circrnafinder:
    input:
        chimericsam=join(WORKDIR,"results","{sample}","STAR_circRNAFinder","{sample}.Chimeric.out.sam"),
        chimericjunction=join(WORKDIR,"results","{sample}","STAR_circRNAFinder","{sample}.Chimeric.out.junction"),
        sjouttab=join(WORKDIR,"results","{sample}","STAR_circRNAFinder","{sample}.SJ.out.tab"),
    output:
        bed=join(WORKDIR,"results","{sample}","circRNA_finder","{sample}.filteredJunctions.bed"),
        ctf=join(WORKDIR,"results","{sample}","circRNA_finder","{sample}.circRNA_finder.counts_table.tsv.filtered"),
        chimericbam=join(WORKDIR,"results","{sample}","circRNA_finder","{sample}.Chimeric.out.sorted.bam")
    params:
        postProcessStarAlignment_script=join(config['circrnafinder_dir'],"postProcessStarAlignment.pl"),
        bsj_min_nreads=config['minreadcount'],
        randomstr=str(uuid.uuid4())
    envmodules: TOOLS["perl"]["version"],TOOLS["samtools"]["version"]
    shell:"""
set -exo pipefail
if [ -d /lscratch/${{SLURM_JOB_ID}} ];then
    TMPDIR="/lscratch/${{SLURM_JOB_ID}}/{params.randomstr}"
else
    TMPDIR="/dev/shm/{params.randomstr}"
fi
if [ ! -d $TMPDIR ];then mkdir -p $TMPDIR;fi

starDir=$(dirname {input.chimericsam})
outDir=$(dirname {output.bed})

if [ -d $outDir ];then rm -rf $outDir;fi
if [ ! -d $outDir ];then mkdir -p $outDir;fi

{params.postProcessStarAlignment_script} \\
    --starDir ${{starDir}}/ \\
    --outDir ${{outDir}}/

echo -ne "chr\tstart\tend\tstrand\tread_count" > {output.ctf}
awk -F"\\t" -v OFS="\\t" -v minreads={params.bsj_min_nreads} '{{if ($5>=minreads) {{print $1,$2,$3,$6,$5}}}}' {output.bed} >> {output.ctf}

"""

def _boolean2str(x): # "1" for True and "0" for False
    if x==True:
        return "1"
    else:
        return "0"

# rule merge_per_sample:
# The output file looks like this:
# | Col# | ColName                              | Description                                                                                                                                                                                                        |
# |------|--------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
# | 1    | circRNA_id                           | chrom:start-end                                                                                                                                                                                                    |
# | 2    | strand                               | "+ or -"                                                                                                                                                                                                           |
# | 3    | <samplename>_ntools                  | number of tools which have this circRNA_id detected                                                                                                                                                                |
# | 4    | <samplename>_circExplorer_read_count |                                                                                                                                                                                                                    |
# | 5    | <samplename>_ciri_read_count         |                                                                                                                                                                                                                    |
# | 6    | <samplename>_dcc_read_count          |                                                                                                                                                                                                                    |
# | 7    | <samplename>_mapsplice_read_count    |                                                                                                                                                                                                                    |
# | 8    | <samplename>_nclscan_read_count      |                                                                                                                                                                                                                    |
# | 9    | circExplorer_annotation              | options are known, low_conf, novel                                                                                                                                                                                 |
# | 10   | ciri_annotation                      | options are exon, intron, intergenic_region                                                                                                                                                                        |
# | 11   | dcc_annotation                       | JunctionType##Start-End Region from CircCoordinates file; 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT;Start-End Region eg. intron-intergenic, exon-exon, intergenic-intron, etc.  |
# | 12   | mapsplice_annotation                 | "fusion_type"##"entropy"; "fusion_type" is either "normal" or "overlapping" ... higher "entropy" values are better!                                                                                                |
# | 13   | nclscan_annotation                   |  1+1 for intragenic 0+1 for intergenic                                                                                                                                                                             |  


if RUN_DCC and RUN_MAPSPLICE and RUN_NCLSCAN:
    rule merge_per_sample:
        input:
            unpack(get_per_sample_files_to_merge)
        output:
            merged_counts=join(WORKDIR,"results","{sample}","{sample}.circRNA_counts.txt.gz"),
        params:
            script=join(SCRIPTS_DIR,"_merge_per_sample_counts_table.py"),
            samplename="{sample}",
            peorse=get_peorse,
            reffa=REF_FA,
            minreadcount=config['minreadcount'] # this filter is redundant as inputs are already pre-filtered.
        envmodules:
            TOOLS["python37"]["version"]
        shell:"""
    set -exo pipefail
    outdir=$(dirname {output.merged_counts})

    parameters=" --circExplorer {input.circExplorer}"
    parameters="$parameters --ciri {input.CIRI}"
    parameters="$parameters --dcc {input.DCC}"
    parameters="$parameters --mapsplice {input.MapSplice}"
    if [[ "{params.peorse}" == "PE" ]]; then    # NCLscan is run only if the sample is PE
        parameters="$parameters --nclscan {input.NCLscan}"
    fi
    parameters="$parameters --min_read_count_reqd {params.minreadcount}"
    parameters="$parameters --reffa {params.reffa}"
    parameters="$parameters --samplename {params.samplename} -o {output.merged_counts}"

    echo "python {params.script} $parameters" 
    python {params.script} $parameters
    """
elif RUN_DCC and not RUN_MAPSPLICE and not RUN_NCLSCAN:
    rule merge_per_sample:
        input:
            unpack(get_per_sample_files_to_merge)
        output:
            merged_counts=join(WORKDIR,"results","{sample}","{sample}.circRNA_counts.txt.gz"),
        params:
            script=join(SCRIPTS_DIR,"_merge_per_sample_counts_table.py"),
            samplename="{sample}",
            peorse=get_peorse,
            reffa=REF_FA,
            minreadcount=config['minreadcount'] # this filter is redundant as inputs are already pre-filtered.
        envmodules:
            TOOLS["python37"]["version"]
        shell:"""
    set -exo pipefail
    outdir=$(dirname {output.merged_counts})

    parameters=" --circExplorer {input.circExplorer}"
    parameters="$parameters --ciri {input.CIRI}"
    parameters="$parameters --dcc {input.DCC}"
    parameters="$parameters --min_read_count_reqd {params.minreadcount}"
    parameters="$parameters --reffa {params.reffa}"
    parameters="$parameters --samplename {params.samplename} -o {output.merged_counts}"

    echo "python {params.script} $parameters" 
    python {params.script} $parameters
    """
elif RUN_DCC and RUN_MAPSPLICE and not RUN_NCLSCAN:
    rule merge_per_sample:
        input:
            unpack(get_per_sample_files_to_merge)
        output:
            merged_counts=join(WORKDIR,"results","{sample}","{sample}.circRNA_counts.txt.gz"),
        params:
            script=join(SCRIPTS_DIR,"_merge_per_sample_counts_table.py"),
            samplename="{sample}",
            peorse=get_peorse,
            reffa=REF_FA,
            minreadcount=config['minreadcount'] # this filter is redundant as inputs are already pre-filtered.
        envmodules:
            TOOLS["python37"]["version"]
        shell:"""
    set -exo pipefail
    outdir=$(dirname {output.merged_counts})

    parameters=" --circExplorer {input.circExplorer}"
    parameters="$parameters --ciri {input.CIRI}"
    parameters="$parameters --dcc {input.DCC}"
    parameters="$parameters --mapsplice {input.MapSplice}"
    parameters="$parameters --min_read_count_reqd {params.minreadcount}"
    parameters="$parameters --reffa {params.reffa}"
    parameters="$parameters --samplename {params.samplename} -o {output.merged_counts}"

    echo "python {params.script} $parameters" 
    python {params.script} $parameters
    """
elif not RUN_DCC and RUN_MAPSPLICE and not RUN_NCLSCAN:
    rule merge_per_sample:
        input:
            unpack(get_per_sample_files_to_merge)
        output:
            merged_counts=join(WORKDIR,"results","{sample}","{sample}.circRNA_counts.txt.gz"),
        params:
            script=join(SCRIPTS_DIR,"_merge_per_sample_counts_table.py"),
            samplename="{sample}",
            peorse=get_peorse,
            reffa=REF_FA,
            minreadcount=config['minreadcount'] # this filter is redundant as inputs are already pre-filtered.
        envmodules:
            TOOLS["python37"]["version"]
        shell:"""
    set -exo pipefail
    outdir=$(dirname {output.merged_counts})

    parameters=" --circExplorer {input.circExplorer}"
    parameters="$parameters --ciri {input.CIRI}"
    parameters="$parameters --mapsplice {input.MapSplice}"
    parameters="$parameters --min_read_count_reqd {params.minreadcount}"
    parameters="$parameters --reffa {params.reffa}"
    parameters="$parameters --samplename {params.samplename} -o {output.merged_counts}"

    echo "python {params.script} $parameters" 
    python {params.script} $parameters
    """
elif not RUN_DCC and not RUN_MAPSPLICE and RUN_NCLSCAN:
    rule merge_per_sample:
        input:
            unpack(get_per_sample_files_to_merge)
        output:
            merged_counts=join(WORKDIR,"results","{sample}","{sample}.circRNA_counts.txt.gz"),
        params:
            script=join(SCRIPTS_DIR,"_merge_per_sample_counts_table.py"),
            samplename="{sample}",
            peorse=get_peorse,
            reffa=REF_FA,
            minreadcount=config['minreadcount'] # this filter is redundant as inputs are already pre-filtered.
        envmodules:
            TOOLS["python37"]["version"]
        shell:"""
    set -exo pipefail
    outdir=$(dirname {output.merged_counts})

    parameters=" --circExplorer {input.circExplorer}"
    parameters="$parameters --ciri {input.CIRI}"
    if [[ "{params.peorse}" == "PE" ]]; then    # NCLscan is run only if the sample is PE
        parameters="$parameters --nclscan {input.NCLscan}"
    fi
    parameters="$parameters --min_read_count_reqd {params.minreadcount}"
    parameters="$parameters --reffa {params.reffa}"
    parameters="$parameters --samplename {params.samplename} -o {output.merged_counts}"

    echo "python {params.script} $parameters" 
    python {params.script} $parameters
    """
else: # DCC, MapSplice and NCLScan are all off!
    rule merge_per_sample:
        input:
            unpack(get_per_sample_files_to_merge)
        output:
            merged_counts=join(WORKDIR,"results","{sample}","{sample}.circRNA_counts.txt.gz"),
        params:
            script=join(SCRIPTS_DIR,"_merge_per_sample_counts_table.py"),
            samplename="{sample}",
            peorse=get_peorse,
            reffa=REF_FA,
            minreadcount=config['minreadcount'] # this filter is redundant as inputs are already pre-filtered.
        envmodules:
            TOOLS["python37"]["version"]
        shell:"""
    set -exo pipefail
    outdir=$(dirname {output.merged_counts})

    parameters=" --circExplorer {input.circExplorer}"
    parameters="$parameters --ciri {input.CIRI}"
    parameters="$parameters --min_read_count_reqd {params.minreadcount}"
    parameters="$parameters --reffa {params.reffa}"
    parameters="$parameters --samplename {params.samplename} -o {output.merged_counts}"

    echo "python {params.script} $parameters" 
    python {params.script} $parameters
    """

localrules: create_master_counts_file
# rule create_master_counts_file:
# merge all per-sample counts tables into a single giant counts matrix and annotate it with known circRNA databases
rule create_master_counts_file:
    input:
        expand(join(WORKDIR,"results","{sample}","{sample}.circRNA_counts.txt.gz"),sample=SAMPLES),
    output:
        matrix=join(WORKDIR,"results","circRNA_master_counts.tsv.gz")
    params:
        script=join(SCRIPTS_DIR,"_make_master_counts_table.py"),
        resultsdir=join(WORKDIR,"results"),
        lookup_table=ANNOTATION_LOOKUP
    envmodules:
        TOOLS['python37']['version']
    shell:"""
set -exo pipefail
count=0
for f in {input};do
    count=$((count+1))
    if [ "$count" == "1" ];then
        infiles="$f"
    else
        infiles="$infiles,$f"
    fi
done

python {params.script} \\
    --counttablelist $infiles \\
    -o {output.matrix}
"""
