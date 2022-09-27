# Find circRNAs using:
# circExplorer2
# ciri2
# CLEAR
# and annotate them



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
        counts_table=join(WORKDIR,"results","{sample}","circExplorer","{sample}.circExplorer.counts_table.tsv")
    params:
        sample="{sample}",
        outdir=join(WORKDIR,"results","{sample}","circExplorer"),
        genepred=rules.create_index.output.genepred_w_geneid,
        reffa=REF_FA,
        script=join(SCRIPTS_DIR,"create_circExplorer_per_sample_counts_table.py")
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

python {params.script} \
    --back_spliced_bed {output.backsplicedjunctions} \
    --circularRNA_known {output.annotations} \
    --low_conf low_conf_$(basename {output.annotations}) \
    -o {output.counts_table}
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
# | 8  | junction_reads_ratio | ratio of circular junction reads calculated by 2#junction_reads/(2#junction_reads+#non_junction_reads). #junction_reads is multiplied by two because a junction read is generated from two ends of circular junction but only counted once while a non-junction read is from one end. It has to be mentioned that the non-junction reads are still possibly from another larger circRNA, so the junction_reads_ratio based on it may be an inaccurate estimation of relative expression of the circRNA. |
# | 9  | x         | type of a circRNA according to positions of its two ends on chromosome (exon, intron or intergenic_region; only available when annotation file is provided)                                                                                                                                                                                                                                                                                                                                             |
# | 10 | gene_id              | ID of the gene(s) where an exonic or intronic circRNA locates                                                                                                                                                                                                                                                                                                                                                                                                                                           |
# | 11 | strand               | strand info of a predicted circRNAs (new in CIRI2)                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
# | 12 | junction_reads_ID    | all of the circular junction read IDs (split by ",")                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
# ref: https://ciri-cookbook.readthedocs.io/en/latest/CIRI2.html#an-example-of-running-ciri2
rule ciri:
    input:
        R1=rules.cutadapt.output.of1,
        R2=rules.cutadapt.output.of2
    output:
        cirilog=join(WORKDIR,"results","{sample}","ciri","{sample}.ciri.log"),
        bwalog=join(WORKDIR,"results","{sample}","ciri","{sample}.bwa.log"),
        ciribam=temp(join(WORKDIR,"results","{sample}","ciri","{sample}.bwa.bam")),
        ciriout=join(WORKDIR,"results","{sample}","ciri","{sample}.ciri.out")
    params:
        sample="{sample}",
        outdir=join(WORKDIR,"results","{sample}","ciri"),
        peorse=get_peorse,
        genepred=rules.create_index.output.genepred_w_geneid,
        reffa=REF_FA,
        bwaindex=BWA_INDEX,
        gtf=REF_GTF,
        ciripl=config['ciri_perl_script']
    threads: getthreads("ciri")
    envmodules: TOOLS["bwa"]["version"], TOOLS["samtools"]["version"]
    shell:"""
set -exo pipefail
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
-A {params.gtf} \\
-G {output.cirilog} -T {threads}
samtools view -@{threads} -bS {params.sample}.bwa.sam > {output.ciribam}
rm -rf {params.sample}.bwa.sam
"""

localrules: merge_per_sample_circRNA_counts
# rule merge_per_sample_circRNA_counts:
# merges counts from circExplorer2 and CIRI2 for all identified circRNAs. 
# The output file columns are:
# | # | ColName                               |
# |---|---------------------------------------|
# | 1 | circRNA_id                            |
# | 2 | strand                                |
# | 3 | <samplename>_circExplorer_read_count  |
# | 4 | <samplename>_ciri_read_count          |
# | 5 | <samplename>_circExplorer_known_novel | --> options are known, low_conf, novel
# | 6 | <samplename>_circRNA_type             | --> options are exon, intron, intergenic_region
# | 7 | <samplename>_ntools                   | --> number of tools calling this BSJ/circRNA
rule merge_per_sample_circRNA_counts:
    input:
        circExplorer_table=rules.circExplorer.output.counts_table,
        ciri_table=rules.ciri.output.ciriout
    output:
        merged_counts=join(WORKDIR,"results","{sample}","circRNA_counts.txt")
    params:
        script=join(SCRIPTS_DIR,"merge_per_sample_counts_table.py"),
        samplename="{sample}"
    shell:"""
set -exo pipefail
python {params.script} \
    --circExplorer {input.circExplorer_table} \
    --ciri {input.ciri_table} \
    --samplename {params.samplename} \
    -o {output.merged_counts}
"""

localrules: create_counts_matrix
# rule create_counts_matrix:
# merge all per-sample counts tables into a single giant counts matrix and annotate it with known circRNA databases
rule create_counts_matrix:
    input:
        expand(join(WORKDIR,"results","{sample}","circRNA_counts.txt"),sample=SAMPLES),
    output:
        matrix=join(WORKDIR,"results","circRNA_counts_matrix.tsv")
    params:
        script=join(SCRIPTS_DIR,"merge_counts_tables_2_counts_matrix.py"),
        resultsdir=join(WORKDIR,"results"),
        lookup_table=ANNOTATION_LOOKUP
    shell:"""
set -exo pipefail
python {params.script} \
    --results_folder {params.resultsdir} \
    --lookup_table {params.lookup_table} \
    -o {output.matrix}
"""

rule create_ciri_count_matrix:
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
