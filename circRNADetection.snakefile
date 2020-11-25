from os.path import join
import sys
import os
import pandas as pd
import yaml
import pprint

def get_fastqs(wildcards):
	d=dict()
	d["R1"]=SAMPLESDF["R1"][wildcards.sample]
	d["R2"]=SAMPLESDF["R2"][wildcards.sample]
	return d

def get_peorse(wildcards):
	return SAMPLESDF["PEorSE"][wildcards.sample]

def check_existence(filename):
	"""Checks if file exists on filesystem
	:param filename <str>: Name of file to check
	"""
	filename=filename.strip()
	if not os.path.exists(filename):
		sys.exit("File: {} does not exists!".format(filename))
	return True


def check_readaccess(filename):
	"""Checks permissions to see if user can read a file
	:param filename <str>: Name of file to check
	"""
	filename=filename.strip()
	check_existence(filename)
	if not os.access(filename,os.R_OK):
		sys.exit("File: {} exists, but user cannot read from file due to permissions!".format(filename))
	return True


def check_writeaccess(filename):
	"""Checks permissions to see if user can write to a file
	:param filename <str>: Name of file to check
	"""
	filename=filename.strip()
	check_existence(filename)
	if not os.access(filename,os.W_OK):
		sys.exit("File: {} exists, but user cannot write to file due to permissions!".format(filename))
	return True

##### load config and sample sheets #####

check_readaccess("config/config.yaml")
configfile: "config/config.yaml"

#resouce absolute path
WORKDIR=config['workdir']
SCRIPTS_DIR=join(WORKDIR,"scripts")
RESOURCES_DIR=join(WORKDIR,"resources")
if not os.path.exists(join(WORKDIR,"fastqs")):
	os.mkdir(join(WORKDIR,"fastqs"))
for f in ["samples", "tools", "cluster"]:
	check_readaccess(config[f])

SAMPLESDF = pd.read_csv(config["samples"],sep="\t",header=0,index_col="sampleName")
SAMPLES = list(SAMPLESDF.index)
SAMPLESDF["R1"]=join(RESOURCES_DIR,"dummy")
SAMPLESDF["R2"]=join(RESOURCES_DIR,"dummy")
SAMPLESDF["PEorSE"]="PE"

for sample in SAMPLES:
	R1file=SAMPLESDF["path_to_R1_fastq"][sample]
	R2file=SAMPLESDF["path_to_R2_fastq"][sample]
	# print(sample,R1file,R2file)
	check_readaccess(R1file)
	R1filenewname=join(WORKDIR,"fastqs",sample+".R1.fastq.gz")
	if os.path.exists(R1filenewname):
		os.remove(R1filenewname)
	os.symlink(R1file,R1filenewname)
	SAMPLESDF["R1"][sample]=R1filenewname
	if str(R2file)!='nan':
		check_readaccess(R2file)
		R2filenewname=join(WORKDIR,"fastqs",sample+".R2.fastq.gz")
		if os.path.exists(R2filenewname):
			os.remove(R2filenewname)
		os.symlink(R2file,R2filenewname)
		SAMPLESDF["R2"][sample]=R2filenewname
	else:
		SAMPLESDF["PEorSE"][sample]="SE"
pprint.pprint(SAMPLESDF["PEorSE"])
exit

## Load tools from YAML file
with open(config["tools"]) as f:
	TOOLS = yaml.safe_load(f)

rule all:
	input:
		expand(join(WORKDIR,"trim","{sample}.R1.trim.fastq.gz"),sample=SAMPLES),
		expand(join(WORKDIR,"trim","{sample}.R2.trim.fastq.gz"),sample=SAMPLES)
		

		# expand(join(WORKDIR,"STAR1p","{sample}_p1.SJ.out.tab"), sample=SAMPLES),
		# join(WORKDIR,"STAR1p","pass1.out.tab"),
		# expand(join(WORKDIR,"STAR2p","{sample}_p2.Chimeric.out.junction"), sample=SAMPLES),
		# expand(join(WORKDIR,"{sample}","{sample}.circularRNA_known.txt"),sample=SAMPLES),
		# expand(join(WORKDIR,"{sample}","{sample}.ciri.out"),sample=SAMPLES),
		# join(WORKDIR,"ciri_count_matrix.txt"),
		# join(WORKDIR,"ciri_count_matrix_with_annotations.txt"),
		# join(WORKDIR,"circExplorer_count_matrix.txt"),
		# join(WORKDIR,"circExplorer_count_matrix_with_annotations.txt")


rule cutadapt:
	input:
		unpack(get_fastqs)
	output:
		of1=join(WORKDIR,"trim","{sample}.R1.trim.fastq.gz"),
		of2=join(WORKDIR,"trim","{sample}.R2.trim.fastq.gz")
	params:
		sample="{sample}",
		peorse=get_peorse,
		adapters=join(RESOURCES_DIR,"TruSeq_and_nextera_adapters.consolidated.fa")
	envmodules: TOOLS["cutadapt"]["version"]
	threads: 56
	shell:"""
if [ "{params.peorse}" == "PE" ];then
	## Paired-end
	cutadapt --pair-filter=any \
	--nextseq-trim=2 \
	--trim-n \
	-n 5 -O 5 \
	-q 10,10 -m 35:35 \
	-b file:{params.adapters} \
	-B file:{params.adapters} \
	-j {threads} \
	-o {output.of1} -p {output.of2} \
	{input.R1} {input.R2}
else
	## Single-end
	cutadapt \
	--nextseq-trim=2 \
	--trim-n \
	-n 5 -O 5 \
	-q 10,10 -m 35 \
	-b file:{params.adapters} \
	-j {threads} \
	-o {output.of1} \
	{input.R1}
	touch {output.of2}
fi
"""

# rule star1p:
# 	input:
# 		unpack(get)
# 		fastqr1=join(WORKDIR,"trim","{sample}.R1.trim.fastq.gz"),
# 		fastqr2=join(WORKDIR,"trim","{sample}.R2.trim.fastq.gz")
# 	output:
# 		junction=join(WORKDIR,"STAR1p","{sample}_p1.SJ.out.tab"),
# 		bam=temp(join(WORKDIR,"STAR1p","{sample}_p1.Aligned.out.bam"))
# 	params:
# 		sample="{sample}",
# 		workdir=WORKDIR,
# 		starversion=STAR_VERSION,
# 		starindexdir=STAR_INDEX_DIR,
# 		gtf=HG38_PLUS_VIRUSES_GTF
# 	threads: 56
# 	shell:"""
# module load STAR/{params.starversion};
# meanrl=$(zcat {input.fastqr1} {input.fastqr2}|grep "^G\|^A\|^C\|^T\|^N"|awk "{{sum=sum+length(\$1);count=count+1}}END{{printf(\\"%d\\",sum/count)}}")
# echo $meanrl
# overhang=$((meanrl-1))
# echo $overhang
# cd {params.workdir}/STAR1p
# STAR --genomeDir {params.starindexdir} \
# --outSAMstrandField None  \
# --outFilterMultimapNmax 20 \
# --alignSJoverhangMin 8 \
# --alignSJDBoverhangMin 1 \
# --outFilterMismatchNmax 999 \
# --outFilterMismatchNoverLmax 0.3  \
# --alignIntronMin 20 \
# --alignIntronMax 1000000 \
# --alignMatesGapMax 1000000 \
# --readFilesIn {input.fastqr1} {input.fastqr2} \
# --readFilesCommand zcat \
# --runThreadN {threads} \
# --outFileNamePrefix {params.sample}_p1. \
# --chimSegmentMin 20 \
# --chimMultimapNmax 10 \
# --chimOutType Junctions \
# --alignTranscriptsPerReadNmax 1200000 \
# --outSAMtype BAM Unsorted \
# --alignEndsProtrude 10 ConcordantPair \
# --outFilterIntronMotifs None \
# --sjdbGTFfile {params.gtf} \
# --outTmpDir=/lscratch/$SLURM_JOB_ID/{params.sample} \
# --sjdbOverhang $overhang
# """

# rule merge_SJ_tabs:
# 	input:
# 		expand(join(WORKDIR,"STAR1p","{sample}_p1.SJ.out.tab"), sample=SAMPLES)
# 	output:
# 		join(WORKDIR,"STAR1p","pass1.out.tab")
# 	threads: 1
# 	shell:"""
# cat {input}|sort -k1,1 -k2,2n|uniq > {output}
# """

# rule star2p:
# 	input:
# 		fastqr1=join(WORKDIR,"trim","{sample}.R1.trim.fastq.gz"),
# 		fastqr2=join(WORKDIR,"trim","{sample}.R2.trim.fastq.gz"),
# 		pass1sjtab=join(WORKDIR,"STAR1p","pass1.out.tab")
# 	output:
# 		junction=join(WORKDIR,"STAR2p","{sample}_p2.Chimeric.out.junction"),
# 		bam=temp(join(WORKDIR,"STAR2p","{sample}_p2.Aligned.out.bam"))
# 	params:
# 		sample="{sample}",
# 		workdir=WORKDIR,
# 		starversion=STAR_VERSION,
# 		starindexdir=STAR_INDEX_DIR,
# 		gtf=HG38_PLUS_VIRUSES_GTF
# 	threads: 56
# 	shell:"""
# module load STAR/{params.starversion}
# meanrl=$(zcat {input.fastqr1} {input.fastqr2}|grep "^G\|^A\|^C\|^T\|^N"|awk "{{sum=sum+length(\$1);count=count+1}}END{{printf(\\"%d\\",sum/count)}}")
# echo $meanrl
# overhang=$((meanrl-1))
# echo $overhang
# cd {params.workdir}/STAR2p
# STAR --genomeDir {params.starindexdir} \
# --outSAMstrandField None  \
# --outFilterType BySJout \
# --outFilterMultimapNmax 20 \
# --alignSJoverhangMin 8 \
# --alignSJDBoverhangMin 1 \
# --outFilterMismatchNmax 999 \
# --outFilterMismatchNoverLmax 0.3  \
# --alignIntronMin 20 \
# --alignIntronMax 2000000 \
# --alignMatesGapMax 2000000 \
# --readFilesIn {input.fastqr1} {input.fastqr2} \
# --readFilesCommand  zcat \
# --runThreadN 56 \
# --outFileNamePrefix {params.sample}_p2. \
# --sjdbFileChrStartEnd {input.pass1sjtab} \
# --chimSegmentMin 20 \
# --chimOutType Junctions \
# --chimMultimapNmax 10 \
# --limitSjdbInsertNsj 5000000 \
# --alignTranscriptsPerReadNmax 2000000 \
# --outSAMtype BAM Unsorted \
# --alignEndsProtrude 10 ConcordantPair \
# --outFilterIntronMotifs None \
# --sjdbGTFfile {params.gtf} \
# --outTmpDir=/lscratch/$SLURM_JOB_ID/{params.sample} \
# --sjdbOverhang $overhang
# """

# rule annotate_circRNA:
# 	input:
# 		junctionfile=join(WORKDIR,"STAR2p","{sample}_p2.Chimeric.out.junction")
# 	output:
# 		backsplicedjunctions=join(WORKDIR,"{sample}","{sample}.back_spliced_junction.bed"),
# 		annotations=join(WORKDIR,"{sample}","{sample}.circularRNA_known.txt")
# 	params:
# 		sample="{sample}",
# 		workdir=WORKDIR,
# 		circexplorerversion=CIRCEXPLORER_VERSION,
# 		genepred=GENEPRED_W_GENEID,
# 		reffa=REFFA
# 	threads: 1
# 	shell:"""
# module load circexplorer2/{params.circexplorerversion}
# cd {params.workdir}/{params.sample}
# CIRCexplorer2 parse \
# 	-t STAR \
# 	{input.junctionfile} > {params.sample}_circexplorer_parse.log 2>&1
# mv back_spliced_junction.bed {output.backsplicedjunctions}
# CIRCexplorer2 annotate \
# -r {params.genepred} \
# -g {params.reffa} \
# -b {output.backsplicedjunctions} \
# -o {output.annotations}
# """


# rule ciri:
# 	input:
# 		fastqr1=join(WORKDIR,"trim","{sample}.R1.trim.fastq.gz"),
# 		fastqr2=join(WORKDIR,"trim","{sample}.R2.trim.fastq.gz")
# 	output:
# 		cirilog=join(WORKDIR,"{sample}","{sample}.ciri.log"),
# 		bwalog=join(WORKDIR,"{sample}","{sample}.bwa.log"),
# 		cirisam=temp(join(WORKDIR,"{sample}","{sample}.sam")),
# 		ciriout=join(WORKDIR,"{sample}","{sample}.ciri.out")
# 	params:
# 		sample="{sample}",
# 		workdir=WORKDIR,
# 		reffa=HG38_PLUS_VIRUSES_FA,
# 		gtf=HG38_PLUS_VIRUSES_GTF,
# 		bwaindex=HG38_PLUS_VIRUSES_BWA_INDEX,
# 		ciripl=CIRI_PERL_SCRIPT,
# 		bwaversion=BWA_VERSION
# 	threads: 56
# 	shell:"""
# module load bwa/{params.bwaversion}
# cd {params.workdir}
# bwa mem -t {threads} -T 19 \
# {params.bwaindex} \
# {input.fastqr1} {input.fastqr2} \
# > {output.cirisam} 2> {output.bwalog}
# perl {params.ciripl} \
# -I {output.cirisam} \
# -O {output.ciriout} \
# -F {params.reffa} \
# -A {params.gtf} \
# -G {output.cirilog} -T {threads}
# """

# rule create_ciri_count_matrix:
# 	input:
# 		expand(join(WORKDIR,"{sample}","{sample}.ciri.out"),sample=SAMPLES)
# 	output:
# 		join(WORKDIR,"ciri_count_matrix.txt"),
# 		join(WORKDIR,"ciri_count_matrix_with_annotations.txt")
# 	params:
# 		script=join(SCRIPTS_DIR,"Create_ciri_count_matrix.py"),
# 		lookup=join(SCRIPTS_DIR,"hg19_hg38_annotated_lookup.txt")
# 	shell:"""
# module load python/3.8
# python {params.script} {params.lookup}
# """

# rule create_circexplorer_count_matrix:
# 	input:
# 		expand(join(WORKDIR,"{sample}","{sample}.circularRNA_known.txt"),sample=SAMPLES)
# 	output:
# 		join(WORKDIR,"circExplorer_count_matrix.txt"),
# 		join(WORKDIR,"circExplorer_count_matrix_with_annotations.txt")
# 	params:
# 		script=join(SCRIPTS_DIR,"Create_circExplorer_count_matrix.py"),
# 		lookup=join(SCRIPTS_DIR,"hg19_hg38_annotated_lookup.txt")
# 	shell:"""
# module load python/3.8
# python {params.script} {params.lookup}
# """
