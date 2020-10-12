from os.path import join
import os

#tool versions
CUTADAPT_VERSION="1.18"
STAR_VERSION="2.7.0f"
BWA_VERSION="0.7.17"
CIRCEXPLORER_VERSION="2.3.5"

#resouce absolute paths
SCRIPTS_DIR="/data/Ziegelbauer_lab/circRNADetection/scripts"
STAR_INDEX_DIR="/data/Ziegelbauer_lab/circRNADetection/resources/STAR_index_no_GTF"
HG38_PLUS_VIRUSES_FA="/data/Ziegelbauer_lab/circRNADetection/resources/hg38_plus_viruses.fa"
HG38_PLUS_VIRUSES_GTF="/data/Ziegelbauer_lab/circRNADetection/resources/hg38_gencodeV35_plus_viruses.gtf"
HG38_PLUS_VIRUSES_BWA_INDEX="/data/Ziegelbauer_lab/circRNADetection/resources/hg38_plus_viruses.fa"
CIRI_PERL_SCRIPT="/data/Ziegelbauer_lab/circRNADetection/resources/CIRI_v2.0.6/CIRI2.pl"
GENEPRED_W_GENEID="/data/Ziegelbauer_lab/circRNADetection/resources/genes.genepred_w_geneid"
REFFA="/data/Ziegelbauer_lab/circRNADetection/resources/hg38_plus_viruses.fa"

#test

#SAMPLES=["C14_2_1","C14_2_2","C14_2_3","C911_1","C911_2","C911_3","GI1_N","GI1_T","GI2_N","GI2_T","GI3_N","GI3_T","GI4_N","GI4_T","GI5_N","GI5_T","GI6_N","GI6_T","LN1_T","NT_1","NT_2","NT_3","Skin1_N","Skin1_T","Skin2_N","Skin2_T","Skin3_N","Skin3_T","Skin4_N","Skin4_T","Skin5_N","Skin5_T","Skin6_N","Skin6_T","Skin7_N","Skin7_T","Skin8_N","Skin8_T","SLK_17_C","SLK_17_G","SLK_20_C","SLK_20_G","SLK_26_C","SLK_26_G","sampleinput_1","sampleinput_2","ORF57IPsample_1","ORF57IPsample_2"]
# SAMPLES=['C14_2_1', 'C14_2_2', 'C14_2_3', 'C911_1', 'C911_2', 'C911_3', 'GI1_N', 'GI1_T', 'GI2_N', 'GI2_T', 'GI3_N', 'GI3_T', 'GI4_N', 'GI4_T', 'GI5_N', 'GI5_T', 'GI6_N', 'GI6_T', 'HK1M', 'HK1R', 'HK3M', 'HK3R', 'Infected25HC_1', 'Infected25HC_2', 'Infected25HC_3', 'Infected25HC_4', 'InfectedVehicle_1', 'InfectedVehicle_2', 'InfectedVehicle_3', 'InfectedVehicle_4', 'KS041411_RNase', 'KS041411_total', 'LN1_T', 'Mock25HC_1', 'Mock25HC_2', 'Mock25HC_3', 'Mock25HC_4', 'MockVehicle_1', 'MockVehicle_2', 'MockVehicle_3', 'MockVehicle_4', 'NT_1', 'NT_2', 'NT_3', 'Skin1_N', 'Skin1_T', 'Skin2_N', 'Skin2_T', 'Skin3_N', 'Skin3_T', 'Skin4_N', 'Skin4_T', 'Skin5_N', 'Skin5_T', 'Skin6_N', 'Skin6_T', 'Skin7_N', 'Skin7_T', 'Skin8_N', 'Skin8_T', 'SLK_17_C', 'SLK_17_G', 'SLK_20_C', 'SLK_20_G', 'SLK_26_C', 'SLK_26_G', 'BCBL1_induced', 'BCBL1_uninduced', "sampleinput_1", "sampleinput_2", "ORF57IPsample_1", "ORF57IPsample_2" ]
SAMPLES=["GI1_N", "GI1_T"]
WORKDIR="/data/Ziegelbauer_lab/circRNADetection/ccbr983_v3"

rule all:
	input:
		expand(join(WORKDIR,"STAR1p","{sample}_p1.SJ.out.tab"), sample=SAMPLES),
		join(WORKDIR,"STAR1p","pass1.out.tab"),
		expand(join(WORKDIR,"STAR2p","{sample}_p2.Chimeric.out.junction"), sample=SAMPLES),
		# expand(join(WORKDIR,"{sample}","{sample}_chrKSHV_only.back_spliced_junction.bed"), sample=SAMPLES),	
		# expand(join(WORKDIR,"{sample}","{sample}_human_only.back_spliced_junction.bed"), sample=SAMPLES),
		# expand(join(WORKDIR,"{sample}","{sample}_human_only.circularRNA_known.txt"), sample=SAMPLES),
		expand(join(WORKDIR,"{sample}","{sample}.circularRNA_known.txt"),sample=SAMPLES),
		expand(join(WORKDIR,"{sample}","{sample}.ciri.out"),sample=SAMPLES)

# rule index_rl:
# 	input:
# 		fq=join(WORKDIR,"fastqs","{sample}.R1.fastq.gz")
# 	output:
# 		rl=join(WORKDIR,"trim","{sample}.index_rl")
# 	params:
# 		sample="{sample}",
# 		workdir=WORKDIR
# 	shell:"""
# python {params.workdir}/get_index_rl.py {input.fq} > {output.rl}
# """

rule cutadapt:
	input:
		f1=join(WORKDIR,"fastqs","{sample}.R1.fastq.gz"),
		f2=join(WORKDIR,"fastqs","{sample}.R2.fastq.gz")
	output:
		of1=join(WORKDIR,"trim","{sample}.R1.trim.fastq.gz"),
		of2=join(WORKDIR,"trim","{sample}.R2.trim.fastq.gz")
	params:
		sample="{sample}",
		cutadaptversion=CUTADAPT_VERSION
	threads: 56
	shell:"""
module load cutadapt/{params.cutadaptversion}
cutadapt --pair-filter=any \
--nextseq-trim=2 \
--trim-n \
-n 5 -O 5 \
-q 10,10 -m 35:35 \
-b file:/data/CCBR_Pipeliner/db/PipeDB/dev/TruSeq_and_nextera_adapters.consolidated.fa \
-B file:/data/CCBR_Pipeliner/db/PipeDB/dev/TruSeq_and_nextera_adapters.consolidated.fa \
-j {threads} \
-o {output.of1} -p {output.of2} \
{input.f1} {input.f2}
"""

rule star1p:
	input:
		fastqr1=join(WORKDIR,"trim","{sample}.R1.trim.fastq.gz"),
		fastqr2=join(WORKDIR,"trim","{sample}.R2.trim.fastq.gz")
	output:
		junction=join(WORKDIR,"STAR1p","{sample}_p1.SJ.out.tab"),
		bam=temp(join(WORKDIR,"STAR1p","{sample}.Aligned.out.bam"))
	params:
		sample="{sample}",
		workdir=WORKDIR,
		starversion=STAR_VERSION,
		starindexdir=STAR_INDEX_DIR,
		gtf=HG38_PLUS_VIRUSES_GTF
	threads: 56
	shell:"""
module load STAR/{params.starversion};
meanrl=$(zcat {input.fastqr1} {input.fastqr2}|grep "^G\|^A\|^C\|^T\|^N"|awk "{{sum=sum+length(\$1);count=count+1}}END{{printf(\\"%d\\",sum/count)}}")
echo $meanrl
overhang=$((meanrl-1))
echo $overhang
cd {params.workdir}/STAR1p
STAR --genomeDir {params.starindexdir} \
--outSAMstrandField None  \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.3  \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--readFilesIn {input.fastqr1} {input.fastqr2} \
--readFilesCommand zcat \
--runThreadN {threads} \
--outFileNamePrefix {params.sample}_p1. \
--chimSegmentMin 20 \
--chimMultimapNmax 10 \
--chimOutType Junctions \
--alignTranscriptsPerReadNmax 1200000 \
--outSAMtype BAM Unsorted \
--alignEndsProtrude 10 ConcordantPair \
--outFilterIntronMotifs None \
--sjdbGTFfile {params.gtf} \
--outTmpDir=/lscratch/$SLURM_JOB_ID/{params.sample} \
--sjdbOverhang $overhang
"""

rule merge_SJ_tabs:
	input:
		expand(join(WORKDIR,"STAR1p","{sample}_p1.SJ.out.tab"), sample=SAMPLES)
	output:
		join(WORKDIR,"STAR1p","pass1.out.tab")
	threads: 1
	shell:"""
cat {input}|sort -k1,1 -k2,2n|uniq > {output}
"""

rule star2p:
	input:
		fastqr1=join(WORKDIR,"trim","{sample}.R1.trim.fastq.gz"),
		fastqr2=join(WORKDIR,"trim","{sample}.R2.trim.fastq.gz"),
		pass1sjtab=join(WORKDIR,"STAR1p","pass1.out.tab")
	output:
		junction=join(WORKDIR,"STAR2p","{sample}_p2.Chimeric.out.junction"),
		bam=temp(join(WORKDIR,"STAR2p","{sample}_p2.Aligned.out.bam"))
	params:
		sample="{sample}",
		workdir=WORKDIR,
		starversion=STAR_VERSION,
		starindexdir=STAR_INDEX_DIR,
		gtf=HG38_PLUS_VIRUSES_GTF
	threads: 56
	shell:"""
module load STAR/{params.starversion}
meanrl=$(zcat {input.fastqr1} {input.fastqr2}|grep "^G\|^A\|^C\|^T\|^N"|awk "{{sum=sum+length(\$1);count=count+1}}END{{printf(\\"%d\\",sum/count)}}")
echo $meanrl
overhang=$((meanrl-1))
echo $overhang
cd {params.workdir}/STAR2p
STAR --genomeDir {params.starindexdir} \
--outSAMstrandField None  \
--outFilterType BySJout \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.3  \
--alignIntronMin 20 \
--alignIntronMax 2000000 \
--alignMatesGapMax 2000000 \
--readFilesIn {input.fastqr1} {input.fastqr2} \
--readFilesCommand  zcat \
--runThreadN 56 \
--outFileNamePrefix {params.sample}_p2. \
--sjdbFileChrStartEnd {input.pass1sjtab} \
--chimSegmentMin 20 \
--chimOutType Junctions \
--chimMultimapNmax 10 \
--limitSjdbInsertNsj 3000000 \
--alignTranscriptsPerReadNmax 2000000 \
--outSAMtype BAM Unsorted \
--alignEndsProtrude 10 ConcordantPair \
--outFilterIntronMotifs None \
--sjdbGTFfile {params.gtf} \
--outTmpDir=/lscratch/$SLURM_JOB_ID/{params.sample} \
--sjdbOverhang $overhang
"""

# rule kshv_only_back_spliced_bed:
# 	input:
# 		join(WORKDIR,"STAR2p","{sample}_p2.Chimeric.out.junction")
# 	output:
# 		join(WORKDIR,"{sample}","{sample}_chrKSHV_only.back_spliced_junction.bed")
# 	params:
# 		sample="{sample}",
# 		workdir=WORKDIR,
# 		circexplorerversion=CIRCEXPLORER_VERSION,
# 		scriptsdir=SCRIPTS_DIR
# 	threads: 1
# 	shell:"""
# module load circexplorer2/{params.circexplorerversion}
# cd {params.workdir}
# python {params.scriptsdir}/filter_junction.py {input} > {params.sample}/{params.sample}_chrKSHV_only_Chimeric.out.junction
# cd {params.sample}
# CIRCexplorer2 parse -t STAR {params.sample}_chrKSHV_only_Chimeric.out.junction > {params.sample}_chrKSHV_only_Chimeric.out.junction_parse.log 2>&1
# mv back_spliced_junction.bed {params.sample}_chrKSHV_only.back_spliced_junction.bed
# """

# rule human_only_back_spliced_bed:
# 	input:
# 		join(WORKDIR,"STAR2p","{sample}_p2.Chimeric.out.junction")
# 	output:
# 		join(WORKDIR,"{sample}","{sample}_human_only.back_spliced_junction.bed")
# 	params:
# 		sample="{sample}",
# 		workdir=WORKDIR,
# 		circexplorerversion=CIRCEXPLORER_VERSION,
# 		scriptsdir=SCRIPTS_DIR
# 	threads: 1
# 	shell:"""
# module load circexplorer2/{params.circexplorerversion}
# cd {params.workdir}
# python {params.scriptsdir}/filter_junction_human.py {input} |sort -k1,1 -k2,2n > {params.sample}/{params.sample}_human_only_Chimeric.out.junction
# cd {params.sample}
# CIRCexplorer2 parse -t STAR {params.sample}_human_only_Chimeric.out.junction > {params.sample}_human_only_Chimeric.out.junction_parse.log 2>&1
# mv back_spliced_junction.bed {params.sample}_human_only.back_spliced_junction.bed
# """

rule annotate_human_circRNA:
	input:
		join(WORKDIR,"{sample}","{sample}_human_only.back_spliced_junction.bed")
	output:
		join(WORKDIR,"{sample}","{sample}_human_only.circularRNA_known.txt")
	params:
		sample="{sample}",
		workdir=WORKDIR,
		circexplorerversion=CIRCEXPLORER_VERSION,
		genepred=GENEPRED_W_GENEID,
		reffa=REFFA
	threads: 1
	shell:"""
module load circexplorer2/{params.circexplorerversion}
cd {params.workdir}
CIRCexplorer2 annotate \
-r {params.genepred} \
-g {params.reffa} \
-b {input} \
-o {params.sample}/{params.sample}_human_only.circularRNA_known.txt
"""

rule annotate_circRNA:
	input:
		junctionfile=join(WORKDIR,"STAR2p","{sample}_p2.Chimeric.out.junction")
	output:
		backsplicedjunctions=join(WORKDIR,"{sample}","{sample}.back_spliced_junction.bed"),
		annotations=join(WORKDIR,"{sample}","{sample}.circularRNA_known.txt")
	params:
		sample="{sample}",
		workdir=WORKDIR,
		circexplorerversion=CIRCEXPLORER_VERSION,
		genepred=GENEPRED_W_GENEID,
		reffa=REFFA
	threads: 1
	shell:"""
module load circexplorer2/{params.circexplorerversion}
cd {params.workdir}/{params.sample}
CIRCexplorer2 parse \
	-t STAR \
	{input.junctionfile} > {params.sample}_circexplorer_parse.log 2>&1
mv back_spliced_junction.bed {output.backsplicedjunctions}
CIRCexplorer2 annotate \
-r {params.genepred} \
-g {params.reffa} \
-b {output.backsplicedjunctions} \
-o {output.annotations}
"""


rule ciri:
	input:
		fastqr1=join(WORKDIR,"trim","{sample}.R1.trim.fastq.gz"),
		fastqr2=join(WORKDIR,"trim","{sample}.R2.trim.fastq.gz")
	output:
		cirilog=join(WORKDIR,"{sample}","{sample}.ciri.log"),
		bwalog=join(WORKDIR,"{sample}","{sample}.bwa.log"),
		cirisam=temp(join(WORKDIR,"{sample}","{sample}.sam")),
		ciriout=join(WORKDIR,"{sample}","{sample}.ciri.out")
	params:
		sample="{sample}",
		workdir=WORKDIR,
		reffa=HG38_PLUS_VIRUSES_FA,
		gtf=HG38_PLUS_VIRUSES_GTF,
		bwaindex=HG38_PLUS_VIRUSES_BWA_INDEX,
		ciripl=CIRI_PERL_SCRIPT,
		bwaversion=BWA_VERSION
	threads: 56
	shell:"""
module load bwa/{params.bwaversion}
cd {params.workdir}
bwa mem -t {threads} -T 19 \
{params.bwaindex} \
{input.fastqr1} {input.fastqr2} \
> {output.cirisam} 2> {output.bwalog}
perl {params.ciripl} \
-I {output.cirisam} \
-O {output.ciriout} \
-F {params.reffa} \
-A {params.gtf} \
-G {output.cirilog} -T {threads}
"""
