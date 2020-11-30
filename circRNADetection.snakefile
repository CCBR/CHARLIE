from os.path import join
import sys
import os
import pandas as pd
import yaml
import pprint

def get_file_size(filename):
	filename=filename.strip()
	if check_readaccess(filename):
		return os.stat(filename).st_size

def get_fastqs(wildcards):
	d=dict()
	d["R1"]=SAMPLESDF["R1"][wildcards.sample]
	d["R2"]=SAMPLESDF["R2"][wildcards.sample]
	return d

def get_raw_and_trim_fastqs(wildcards):
	d=dict()
	d["R1"]=SAMPLESDF["R1"][wildcards.sample]
	d["R2"]=SAMPLESDF["R2"][wildcards.sample]
	d["R1trim"]=join(WORKDIR,"results",wildcards.sample,wildcards.sample+".R1.trim.fastq.gz")
	d["R2trim"]=join(WORKDIR,"results",wildcards.sample,wildcards.sample+".R2.trim.fastq.gz")	
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
if not os.path.exists(join(WORKDIR,"results")):
	os.mkdir(join(WORKDIR,"results"))
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
	if not os.path.exists(R1filenewname):
		os.symlink(R1file,R1filenewname)
	SAMPLESDF["R1"][sample]=R1filenewname
	if str(R2file)!='nan':
		check_readaccess(R2file)
		R2filenewname=join(WORKDIR,"fastqs",sample+".R2.fastq.gz")
		if not os.path.exists(R2filenewname):
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
		## cutadapt
		expand(join(WORKDIR,"results","{sample}","{sample}.R1.trim.fastq.gz"),sample=SAMPLES),
		expand(join(WORKDIR,"results","{sample}","{sample}.R2.trim.fastq.gz"),sample=SAMPLES),
		## fastqc
		expand(join(WORKDIR,"qc","fastqc","{sample}.R1.trim_fastqc.zip"),sample=SAMPLES),
		## star1p
		expand(join(WORKDIR,"results","{sample}","STAR1p","{sample}_p1.SJ.out.tab"),sample=SAMPLES),
		## merge junctions
		join(WORKDIR,"results","pass1.out.tab"),
		## star2p
		expand(join(WORKDIR,"results","{sample}","STAR2p","{sample}_p2.Chimeric.out.junction"),sample=SAMPLES),
		expand(join(WORKDIR,"results","{sample}","STAR2p","{sample}_p2.Aligned.out.bam"),sample=SAMPLES),
		## circExplorer
		expand(join(WORKDIR,"results","{sample}","STAR2p","{sample}_p2.Chimeric.out.junction"), sample=SAMPLES),
		expand(join(WORKDIR,"results","{sample}","circExplorer","{sample}.circularRNA_known.txt"),sample=SAMPLES),
		## ciri
		expand(join(WORKDIR,"results","{sample}","ciri","{sample}.ciri.out"),sample=SAMPLES),
		## ciri aggregate count matrix
		join(WORKDIR,"results","ciri_count_matrix.txt"),
		join(WORKDIR,"results","ciri_count_matrix_with_annotations.txt"),
		## circExplorer aggregate count matrix
		join(WORKDIR,"results","circExplorer_count_matrix.txt"),
		join(WORKDIR,"results","circExplorer_count_matrix_with_annotations.txt")


rule cutadapt:
	input:
		unpack(get_fastqs)
	output:
		of1=join(WORKDIR,"results","{sample}","{sample}.R1.trim.fastq.gz"),
		of2=join(WORKDIR,"results","{sample}","{sample}.R2.trim.fastq.gz")
	params:
		sample="{sample}",
		workdir=WORKDIR,
		outdir=join(WORKDIR,"results","{sample}"),
		peorse=get_peorse,
		adapters=join(RESOURCES_DIR,"TruSeq_and_nextera_adapters.consolidated.fa")
	envmodules: TOOLS["cutadapt"]["version"]
	threads: 56
	shell:"""
if [ ! -d {params.outdir} ];then mkdir {params.outdir};fi
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

rule fastqc:
	input:
		unpack(get_raw_and_trim_fastqs)
	output:
		join(WORKDIR,"qc","fastqc","{sample}.R1.trim_fastqc.zip")
	params:
		outdir=join(WORKDIR,"qc","fastqc")
	threads: 16
	envmodules: TOOLS["fastqc"]["version"]
	shell:"""
files=""
for f in {input};do
if [ "$(wc $f|awk '{{print $1}}')" != "0" ];then
files=$(echo -ne "$files $f")
fi
done
fastqc $files -t {threads} -o {params.outdir}
"""

rule star1p:
	input:
		R1=rules.cutadapt.output.of1,
		R2=rules.cutadapt.output.of2
	output:
		junction=join(WORKDIR,"results","{sample}","STAR1p","{sample}_p1.SJ.out.tab")
	params:
		sample="{sample}",
		peorse=get_peorse,
		workdir=WORKDIR,
		outdir=join(WORKDIR,"results","{sample}","STAR1p"),
		starindexdir=config['star_index_dir'],
		alignTranscriptsPerReadNmax=TOOLS["star"]["alignTranscriptsPerReadNmax"],
		gtf=config['hg38_plus_virusus_gtf']
	envmodules: TOOLS["star"]["version"]
	threads: 56
	shell:"""
if [ -d /dev/shm/{params.sample} ];then rm -rf /dev/shm/{params.sample};fi
if [ ! -d {params.outdir} ];then mkdir {params.outdir};fi
if [ "{params.peorse}" == "PE" ];then
# paired-end
	overhang=$(zcat {input} | awk -v maxlen=100 'NR%4==2 {{if (length($1) > maxlen+0) maxlen=length($1)}}; END {{print maxlen-1}}')
	echo "sjdbOverhang for STAR: ${{overhang}}"
	cd {params.outdir}
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
	--readFilesIn {input.R1} {input.R2} \
	--readFilesCommand zcat \
	--runThreadN {threads} \
	--outFileNamePrefix {params.sample}_p1. \
	--chimSegmentMin 20 \
	--chimMultimapNmax 10 \
	--chimOutType Junctions \
	--alignTranscriptsPerReadNmax {params.alignTranscriptsPerReadNmax} \
	--outSAMtype BAM Unsorted \
	--alignEndsProtrude 10 ConcordantPair \
	--outFilterIntronMotifs None \
	--sjdbGTFfile {params.gtf} \
	--outTmpDir=/dev/shm/{params.sample} \
	--sjdbOverhang $overhang

else

#single-end
	overhang=$(zcat {input.R1} | awk -v maxlen=100 'NR%4==2 {{if (length($1) > maxlen+0) maxlen=length($1)}}; END {{print maxlen-1}}')
	echo "sjdbOverhang for STAR: ${{overhang}}"
	cd {params.outdir}
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
	--readFilesIn {input.R1} \
	--readFilesCommand zcat \
	--runThreadN {threads} \
	--outFileNamePrefix {params.sample}_p1. \
	--chimSegmentMin 20 \
	--chimMultimapNmax 10 \
	--chimOutType Junctions \
	--alignTranscriptsPerReadNmax {params.alignTranscriptsPerReadNmax} \
	--outSAMtype BAM Unsorted \
	--alignEndsProtrude 10 ConcordantPair \
	--outFilterIntronMotifs None \
	--sjdbGTFfile {params.gtf} \
	--outTmpDir=/dev/shm/{params.sample} \
	--sjdbOverhang $overhang

fi
rm -rf {params.outdir}/{params.sample}_p1._STARgenome
rm -rf {params.outdir}/{params.sample}_p1.Aligned.out.bam
"""

rule merge_SJ_tabs:
	input:
		expand(join(WORKDIR,"results","{sample}","STAR1p","{sample}_p1.SJ.out.tab"), sample=SAMPLES)
	output:
		pass1sjtab=join(WORKDIR,"results","pass1.out.tab")
	threads: 1
	shell:"""
cat {input} |sort|uniq|awk -F \"\\t\" '{{if ($5>0 && $6==1) {{print}}}}'|cut -f1-4|sort -k1,1 -k2,2n|uniq > {output.pass1sjtab}
"""

rule star2p:
	input:
		R1=rules.cutadapt.output.of1,
		R2=rules.cutadapt.output.of2,
		pass1sjtab=rules.merge_SJ_tabs.output.pass1sjtab
	output:
		junction=join(WORKDIR,"results","{sample}","STAR2p","{sample}_p2.Chimeric.out.junction"),
		bam=join(WORKDIR,"results","{sample}","STAR2p","{sample}_p2.Aligned.out.bam")
	params:
		sample="{sample}",
		peorse=get_peorse,
		workdir=WORKDIR,
		outdir=join(WORKDIR,"results","{sample}","STAR2p"),
		starindexdir=config['star_index_dir'],
		alignTranscriptsPerReadNmax=TOOLS["star"]["alignTranscriptsPerReadNmax"],
		gtf=config['hg38_plus_virusus_gtf']
	envmodules: TOOLS["star"]["version"]
	threads: 56
	shell:"""
if [ -d /dev/shm/{params.sample} ];then rm -rf /dev/shm/{params.sample};fi
if [ ! -d {params.outdir} ];then mkdir {params.outdir};fi
limitSjdbInsertNsj=$(wc -l {input.pass1sjtab}|awk '{{print $1+1}}')
if [ "$limitSjdbInsertNsj" -lt "400000" ];then limitSjdbInsertNsj="400000";fi
if [ "{params.peorse}" == "PE" ];then
# paired-end
	overhang=$(zcat {input.R1} {input.R2} | awk -v maxlen=100 'NR%4==2 {{if (length($1) > maxlen+0) maxlen=length($1)}}; END {{print maxlen-1}}')
	echo "sjdbOverhang for STAR: ${{overhang}}"
	cd {params.outdir}
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
	--readFilesIn {input.R1} {input.R2} \
	--readFilesCommand  zcat \
	--runThreadN 56 \
	--outFileNamePrefix {params.sample}_p2. \
	--sjdbFileChrStartEnd {input.pass1sjtab} \
	--chimSegmentMin 20 \
	--chimOutType Junctions \
	--chimMultimapNmax 10 \
	--limitSjdbInsertNsj $limitSjdbInsertNsj \
	--alignTranscriptsPerReadNmax {params.alignTranscriptsPerReadNmax} \
	--outSAMtype BAM Unsorted \
	--alignEndsProtrude 10 ConcordantPair \
	--outFilterIntronMotifs None \
	--sjdbGTFfile {params.gtf} \
	--outTmpDir=/dev/shm/{params.sample} \
	--sjdbOverhang $overhang

else
#single-end
	overhang=$(zcat {input.R1} | awk -v maxlen=100 'NR%4==2 {{if (length($1) > maxlen+0) maxlen=length($1)}}; END {{print maxlen-1}}')
	echo "sjdbOverhang for STAR: ${{overhang}}"
	cd {params.outdir}
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
	--readFilesIn {input.R1} \
	--readFilesCommand  zcat \
	--runThreadN 56 \
	--outFileNamePrefix {params.sample}_p2. \
	--sjdbFileChrStartEnd {input.pass1sjtab} \
	--chimSegmentMin 20 \
	--chimOutType Junctions \
	--chimMultimapNmax 10 \
	--limitSjdbInsertNsj $limitSjdbInsertNsj \
	--alignTranscriptsPerReadNmax {params.alignTranscriptsPerReadNmax} \
	--outSAMtype BAM Unsorted \
	--alignEndsProtrude 10 ConcordantPair \
	--outFilterIntronMotifs None \
	--sjdbGTFfile {params.gtf} \
	--outTmpDir=/dev/shm/{params.sample} \
	--sjdbOverhang $overhang
fi
rm -rf {params.outdir}/{params.sample}_p2._STARgenome
"""

rule annotate_circRNA:
	input:
		junctionfile=rules.star2p.output.junction
	output:
		backsplicedjunctions=join(WORKDIR,"results","{sample}","circExplorer","{sample}.back_spliced_junction.bed"),
		annotations=join(WORKDIR,"results","{sample}","circExplorer","{sample}.circularRNA_known.txt")
	params:
		sample="{sample}",
		workdir=WORKDIR,
		outdir=join(WORKDIR,"results","{sample}","circExplorer"),
		genepred=config['genepred_w_geneid'],
		reffa=config['reffa']
	threads: 1
	envmodules: TOOLS["circexplorer"]["version"]
	shell:"""
if [ ! -d {params.outdir} ];then mkdir {params.outdir};fi
cd {params.outdir}
mv {input.junctionfile} {input.junctionfile}.original
grep -v junction_type {input.junctionfile}.original > {input.junctionfile}
CIRCexplorer2 parse \
	-t STAR \
	{input.junctionfile} > {params.sample}_circexplorer_parse.log 2>&1
mv back_spliced_junction.bed {output.backsplicedjunctions}
mv {input.junctionfile}.original {input.junctionfile}
CIRCexplorer2 annotate \
-r {params.genepred} \
-g {params.reffa} \
-b {output.backsplicedjunctions} \
-o $(basename {output.annotations}) \
--low-confidence
"""

rule ciri:
	input:
		R1=rules.cutadapt.output.of1,
		R2=rules.cutadapt.output.of2
	output:
		cirilog=join(WORKDIR,"results","{sample}","ciri","{sample}.ciri.log"),
		bwalog=join(WORKDIR,"results","{sample}","ciri","{sample}.bwa.log"),
		ciribam=join(WORKDIR,"results","{sample}","ciri","{sample}.bwa.bam"),
		ciriout=join(WORKDIR,"results","{sample}","ciri","{sample}.ciri.out")
	params:
		sample="{sample}",
		workdir=WORKDIR,
		outdir=join(WORKDIR,"results","{sample}","ciri"),
		peorse=get_peorse,
		genepred=config['genepred_w_geneid'],
		reffa=config['reffa'],
		bwaindex=config['hg38_plus_virusus_bwa_index'],
		gtf=config['hg38_plus_virusus_gtf'],
		ciripl=config['ciri_perl_script']
	threads: 56
	envmodules: TOOLS["bwa"]["version"], TOOLS["samtools"]["version"]
	shell:"""
cd {params.outdir}
if [ "{params.peorse}" == "PE" ];then
	## paired-end
	bwa mem -t {threads} -T 19 \
	{params.bwaindex} \
	{input.R1} {input.R2} \
	> {params.sample}.bwa.sam 2> {output.bwalog}
else
	## single-end
	bwa mem -t {threads} -T 19 \
	{params.bwaindex} \
	{input.R1} \
	> {params.sample}.bwa.sam 2> {output.bwalog}
fi
perl {params.ciripl} \
-I {params.sample}.bwa.sam \
-O {output.ciriout} \
-F {params.reffa} \
-A {params.gtf} \
-G {output.cirilog} -T {threads}
samtools view -@{threads} -bS {params.sample}.bwa.sam > {output.ciribam}
rm -rf {params.sample}.bwa.sam
"""

rule create_ciri_count_matrix:
	input:
		expand(join(WORKDIR,"results","{sample}","ciri","{sample}.ciri.out"),sample=SAMPLES)
	output:
		matrix=join(WORKDIR,"results","ciri_count_matrix.txt"),
		matrix_w_annotations=join(WORKDIR,"results","ciri_count_matrix_with_annotations.txt")
	params:
		script=join(SCRIPTS_DIR,"Create_ciri_count_matrix.py"),
		lookup=join(RESOURCES_DIR,"hg19_hg38_annotated_lookup.txt"),
		outdir=join(WORKDIR,"results")
	envmodules: "python/3.7"
	shell:"""
cd {params.outdir}
python {params.script} {params.lookup}
# mv $(basename {output.matrix}) {output.matrix}
# mv $(basename {output.matrix_w_annotations}) {output.matrix_w_annotations}
"""

rule create_circexplorer_count_matrix:
	input:
		expand(join(WORKDIR,"results","{sample}","circExplorer","{sample}.circularRNA_known.txt"),sample=SAMPLES)
	output:
		matrix=join(WORKDIR,"results","circExplorer_count_matrix.txt"),
		matrix_w_annotations=join(WORKDIR,"results","circExplorer_count_matrix_with_annotations.txt")
	params:
		script=join(SCRIPTS_DIR,"Create_circExplorer_count_matrix.py"),
		lookup=join(RESOURCES_DIR,"hg19_hg38_annotated_lookup.txt"),
		outdir=join(WORKDIR,"results")
	envmodules: "python/3.7"
	shell:"""
cd {params.outdir}
python {params.script} {params.lookup}
# mv $(basename {output.matrix}) {output.matrix}
# mv $(basename {output.matrix_w_annotations}) {output.matrix_w_annotations} 
"""
