from os.path import join
import sys
import os
import pandas as pd
import yaml


# no truncations during print pandas data frames
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)

def get_clear_target_files(SAMPLES,runclear):
	targetfiles=[]
	if runclear==True or runclear=="True":
		for s in SAMPLES:
			targetfiles.append(join(WORKDIR,"results",s,"CLEAR","quant.txt"))
	return targetfiles


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
	d["R1trim"]=join(WORKDIR,"results",wildcards.sample,"trim",wildcards.sample+".R1.trim.fastq.gz")
	d["R2trim"]=join(WORKDIR,"results",wildcards.sample,"trim",wildcards.sample+".R2.trim.fastq.gz")	
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

# check_readaccess("config/config.yaml")
# configfile: "config/config.yaml"

## set memory limit 
## used for sambamba sort, etc
MEMORYG="100G"

#resouce absolute path
WORKDIR=config['workdir']
SCRIPTS_DIR=config['scriptsdir']
RESOURCES_DIR=config['resourcesdir']
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
	SAMPLESDF.loc[[sample],"R1"]=R1filenewname
	if str(R2file)!='nan':
		check_readaccess(R2file)
		R2filenewname=join(WORKDIR,"fastqs",sample+".R2.fastq.gz")
		if not os.path.exists(R2filenewname):
			os.symlink(R2file,R2filenewname)
		SAMPLESDF.loc[[sample],"R2"]=R2filenewname
	else:
		SAMPLESDF.loc[[sample],"PEorSE"]="SE"
# print(SAMPLESDF)

## Load tools from YAML file
with open(config["tools"]) as f:
	TOOLS = yaml.safe_load(f)

localrules: multiqc

rule all:
	input:
		## cutadapt
		expand(join(WORKDIR,"results","{sample}","trim","{sample}.R1.trim.fastq.gz"),sample=SAMPLES),
		expand(join(WORKDIR,"results","{sample}","trim","{sample}.R2.trim.fastq.gz"),sample=SAMPLES),
		## fastqc
		expand(join(WORKDIR,"qc","fastqc","{sample}.R1.trim_fastqc.zip"),sample=SAMPLES),
		## star1p
		expand(join(WORKDIR,"results","{sample}","STAR1p","{sample}_p1.SJ.out.tab"),sample=SAMPLES),
		## merge junctions
		join(WORKDIR,"results","pass1.out.tab"),
		## star2p
		expand(join(WORKDIR,"results","{sample}","STAR2p","{sample}_p2.Chimeric.out.junction"),sample=SAMPLES),
		expand(join(WORKDIR,"results","{sample}","STAR2p","{sample}_p2.Aligned.sortedByCoord.out.bam"),sample=SAMPLES),
		## picard MarkDuplicates metrics
		expand(join(WORKDIR,"qc","picard_MarkDuplicates","{sample}.MarkDuplicates.metrics.txt"),sample=SAMPLES),
		## BSJ bam
		expand(join(WORKDIR,"results","{sample}","STAR2p","{sample}.BSJ.bam"),sample=SAMPLES),
		## circExplorer
		expand(join(WORKDIR,"results","{sample}","circExplorer","{sample}.circularRNA_known.txt"),sample=SAMPLES),
		## ciri
		expand(join(WORKDIR,"results","{sample}","ciri","{sample}.ciri.out"),sample=SAMPLES),
		## ciri aggregate count matrix
		join(WORKDIR,"results","ciri_count_matrix.txt"),
		## circExplorer aggregate count matrix
		join(WORKDIR,"results","circExplorer_count_matrix.txt"),
		join(WORKDIR,"results","circExplorer_BSJ_count_matrix.txt"),
		## BSJ bams and bigwigs
		expand(join(WORKDIR,"results","{sample}","STAR2p","{sample}.BSJ.hg38.bam"),sample=SAMPLES),
		expand(join(WORKDIR,"results","{sample}","STAR2p","{sample}.BSJ.hg38.bw"),sample=SAMPLES),
		## Spliced reads bams and bigwigs
		expand(join(WORKDIR,"results","{sample}","STAR2p","{sample}.spliced_reads.hg38.bam"),sample=SAMPLES),
		expand(join(WORKDIR,"results","{sample}","STAR2p","{sample}.spliced_reads.hg38.bw"),sample=SAMPLES),
		## CLEAR quant output
		get_clear_target_files(SAMPLES,config['run_clear']),
		## CIRI2 and CIRCEXPLORER2 VENN DIAGRAM
		expand(join(WORKDIR,"results","{sample}","{sample}.venn.png"),sample=SAMPLES),
		expand(join(WORKDIR,"results","{sample}","{sample}.cirionly.lst"),sample=SAMPLES),
		expand(join(WORKDIR,"results","{sample}","{sample}.circexploreronly.lst"),sample=SAMPLES),
		expand(join(WORKDIR,"results","{sample}","{sample}.common.lst"),sample=SAMPLES),
		## multiqc report
		join(WORKDIR,"multiqc_report.html")

rule cutadapt:
	input:
		unpack(get_fastqs)
	output:
		of1=join(WORKDIR,"results","{sample}","trim","{sample}.R1.trim.fastq.gz"),
		of2=join(WORKDIR,"results","{sample}","trim","{sample}.R2.trim.fastq.gz")
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
		gtf=config['ref_gtf']
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
	--outSAMtype None \
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
	--outSAMtype None \
	--alignEndsProtrude 10 ConcordantPair \
	--outFilterIntronMotifs None \
	--sjdbGTFfile {params.gtf} \
	--outTmpDir=/dev/shm/{params.sample} \
	--sjdbOverhang $overhang

fi
rm -rf {params.outdir}/{params.sample}_p1._STARgenome
# rm -rf {params.outdir}/{params.sample}_p1.Aligned.out.bam
"""

rule merge_SJ_tabs:
	input:
		expand(join(WORKDIR,"results","{sample}","STAR1p","{sample}_p1.SJ.out.tab"), sample=SAMPLES)
	output:
		pass1sjtab=join(WORKDIR,"results","pass1.out.tab")
	params:
		script1=join(SCRIPTS_DIR,"apply_junction_filters.py"),
		regions=config['regions'],
		filter1regions=config['star1passfilter1regions'],
		filter1_noncanonical=config['star1passfilter1_noncanonical'],
		filter1_unannotated=config['star1passfilter1_unannotated'],
		filter2_noncanonical=config['star1passfilter2_noncanonical'],
		filter2_unannotated=config['star1passfilter2_unannotated']
	threads: 1
	shell:"""
cat {input} | \
python {params.script1} \
	--regions {params.regions} \
	--filter1regions {params.filter1regions} \
	--filter1_noncanonical {params.filter1_noncanonical} \
	--filter1_unannotated {params.filter1_unannotated} \
	--filter2_noncanonical {params.filter2_noncanonical} \
	--filter2_unannotated {params.filter1_unannotated} | \
cut -f1-4 | sort -k1,1 -k2,2n | uniq > {output.pass1sjtab}
"""

rule star2p:
	input:
		R1=rules.cutadapt.output.of1,
		R2=rules.cutadapt.output.of2,
		pass1sjtab=rules.merge_SJ_tabs.output.pass1sjtab
	output:
		junction=join(WORKDIR,"results","{sample}","STAR2p","{sample}_p2.Chimeric.out.junction"),
		bam=join(WORKDIR,"results","{sample}","STAR2p","{sample}_p2.Aligned.sortedByCoord.out.bam")
	params:
		sample="{sample}",
		peorse=get_peorse,
		workdir=WORKDIR,
		outdir=join(WORKDIR,"results","{sample}","STAR2p"),
		starindexdir=config['star_index_dir'],
		alignTranscriptsPerReadNmax=TOOLS["star"]["alignTranscriptsPerReadNmax"],
		gtf=config['ref_gtf']
	envmodules: TOOLS["star"]["version"],TOOLS["sambamba"]["version"]
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
	--chimOutType Junctions WithinBAM \
	--chimMultimapNmax 10 \
	--limitSjdbInsertNsj $limitSjdbInsertNsj \
	--alignTranscriptsPerReadNmax {params.alignTranscriptsPerReadNmax} \
	--outSAMtype BAM SortedByCoordinate \
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
	--chimOutType Junctions WithinBAM \
	--chimMultimapNmax 10 \
	--limitSjdbInsertNsj $limitSjdbInsertNsj \
	--alignTranscriptsPerReadNmax {params.alignTranscriptsPerReadNmax} \
	--outSAMtype BAM SortedByCoordinate \
	--alignEndsProtrude 10 ConcordantPair \
	--outFilterIntronMotifs None \
	--sjdbGTFfile {params.gtf} \
	--outTmpDir=/dev/shm/{params.sample} \
	--sjdbOverhang $overhang
fi
rm -rf {params.outdir}/{params.sample}_p2._STARgenome
## ensure the star2p file is indexed ... is should already be sorted by STAR
if [ ! -f {output.bam}.bai ];then
sambamba index {output.bam}
fi
"""

rule estimate_duplication:
	input:
		bam=rules.star2p.output.bam
	output:
		metrics=join(WORKDIR,"qc","picard_MarkDuplicates","{sample}.MarkDuplicates.metrics.txt")
	params:
		sample="{sample}",
		memG=MEMORYG,
	envmodules: TOOLS["picard"]["version"]
	shell:"""
java -Xmx{params.memG} -jar ${{PICARD_JARPATH}}/picard.jar MarkDuplicates I={input.bam} O=/dev/shm/{params.sample}.mark_dup.bam M={output.metrics}
"""

rule create_BSJ_bam:
	input:
		junction=rules.star2p.output.junction,
		bam=rules.star2p.output.bam
	output:
		readids=join(WORKDIR,"results","{sample}","STAR2p","{sample}.BSJ.readids"),
		bam=join(WORKDIR,"results","{sample}","STAR2p","{sample}.BSJ.bam")
	params:
		sample="{sample}",
		memG=MEMORYG,
		script1=join(SCRIPTS_DIR,"junctions2readids.py"),
		script2=join(SCRIPTS_DIR,"filter_bam_by_readids.py"),
		script3=join(SCRIPTS_DIR,"filter_bam_for_BSJs.py"),
		outdir=join(WORKDIR,"results","{sample}","STAR2p")
	envmodules: TOOLS["python37"]["version"],TOOLS["sambamba"]["version"],TOOLS["samtools"]["version"]
	threads : 4
	shell:"""
cd {params.outdir}

## get BSJ readids along with chrom,site,cigar etc.
python {params.script1} -j {input.junction} > {output.readids}

## extract only the uniq readids
cut -f1 {output.readids} | sort | uniq > /dev/shm/{params.sample}.readids

## downsize the star2p bam file to a new bam file with only BSJ reads ... these may still contain alignments which are chimeric but not BSJ
## note the argument --readids here is just a list of readids
python {params.script2} --inputBAM {input.bam} --outputBAM /dev/shm/{params.sample}.chimeric.bam --readids /dev/shm/{params.sample}.readids
sambamba sort --memory-limit={params.memG} --tmpdir=/dev/shm --nthreads={threads} --out=/dev/shm/{params.sample}.chimeric.sorted.bam /dev/shm/{params.sample}.chimeric.bam
rm -f /dev/shm/{params.sample}.chimeric.bam*

## using the downsized star2p bam file containing chimeric alignments ...included all the BSJs... we now extract only the BSJs
## note the argument --readids here is a tab delimited file created by junctions2readids.py ... reaids,chrom,strand,sites,cigars,etc.
python {params.script3} --inputBAM /dev/shm/{params.sample}.chimeric.sorted.bam --outputBAM /dev/shm/{params.sample}.BSJs.tmp.bam --readids {output.readids}
sambamba sort --memory-limit={params.memG} --tmpdir=/dev/shm --nthreads={threads} --out=/dev/shm/{params.sample}.BSJs.tmp.sorted.bam /dev/shm/{params.sample}.BSJs.tmp.bam
rm -f /dev/shm/{params.sample}.BSJs.tmp.bam*

## some alignments are repeated/duplicated in the output for some reason ... hence deduplicating
samtools view -H /dev/shm/{params.sample}.BSJs.tmp.sorted.bam > /dev/shm/{params.sample}.BSJs.tmp.dedup.sam
samtools view /dev/shm/{params.sample}.BSJs.tmp.sorted.bam | sort | uniq >> /dev/shm/{params.sample}.BSJs.tmp.dedup.sam
samtools view -bS /dev/shm/{params.sample}.BSJs.tmp.dedup.sam > /dev/shm/{params.sample}.BSJs.tmp.dedup.bam
sambamba sort --memory-limit={params.memG} --tmpdir=/dev/shm --nthreads={threads} --out={output.bam} /dev/shm/{params.sample}.BSJs.tmp.dedup.bam
rm -f /dev/shm/{params.sample}.BSJs.tmp.dedup.bam*
"""

rule create_spliced_reads_bam:
	input:
		bam=rules.star2p.output.bam,
		tab=rules.merge_SJ_tabs.output.pass1sjtab
	output:
		bam=join(WORKDIR,"results","{sample}","STAR2p","{sample}.spliced_reads.bam")
	params:
		sample="{sample}",
		memG=MEMORYG,
		script1=join(SCRIPTS_DIR,"filter_bam_for_splice_reads.py"),
		outdir=join(WORKDIR,"results","{sample}","STAR2p")
	envmodules: TOOLS["python37"]["version"],TOOLS["sambamba"]["version"],TOOLS["samtools"]["version"]
	threads : 4
	shell:"""
cd {params.outdir}
python {params.script1} --inbam {input.bam} --outbam /dev/shm/{params.sample}.SR.bam --tab {input.tab}
sambamba sort --memory-limit={params.memG} --tmpdir=/dev/shm --nthreads={threads} --out={output.bam} /dev/shm/{params.sample}.SR.bam
rm -f /dev/shm/{params.sample}.SR.bam*
"""

rule split_splice_reads_BAM_create_BW:
	input:
		bam=rules.create_spliced_reads_bam.output.bam
	output:
		bam=join(WORKDIR,"results","{sample}","STAR2p","{sample}.spliced_reads.hg38.bam"),
		bw=join(WORKDIR,"results","{sample}","STAR2p","{sample}.spliced_reads.hg38.bw")
	params:
		sample="{sample}",
		workdir=WORKDIR,
		memG=MEMORYG,
		outdir=join(WORKDIR,"results","{sample}","STAR2p"),
		regions=config['regions']
	threads: 2
	envmodules: TOOLS["samtools"]["version"],TOOLS["bedtools"]["version"],TOOLS["ucsc"]["version"],TOOLS["sambamba"]["version"]
	shell:"""
cd {params.outdir}
bam_basename="$(basename {input.bam})"
while read a b;do
bam="${{bam_basename%.*}}.${{a}}.bam"
samtools view {input.bam} $b -b > /dev/shm/${{bam%.*}}.tmp.bam
sambamba sort --memory-limit={params.memG} --tmpdir=/dev/shm --nthreads={threads} --out=$bam /dev/shm/${{bam%.*}}.tmp.bam
bw="${{bam%.*}}.bw"
bdg="${{bam%.*}}.bdg"
sizes="${{bam%.*}}.sizes"
bedtools genomecov -bg -ibam $bam > $bdg
bedSort $bdg $bdg
if [ "$(wc -l $bdg|awk '{{print $1}}')" != "0" ];then
samtools view -H $bam|grep ^@SQ|cut -f2,3|sed "s/SN://g"|sed "s/LN://g" > $sizes
bedGraphToBigWig $bdg $sizes $bw
else
touch $bw
fi
rm -f $bdg $sizes
done < {params.regions}
"""	


rule split_BAM_create_BW:
	input:
		bam=rules.create_BSJ_bam.output.bam
	output:
		bam=join(WORKDIR,"results","{sample}","STAR2p","{sample}.BSJ.hg38.bam"),
		bw=join(WORKDIR,"results","{sample}","STAR2p","{sample}.BSJ.hg38.bw")
	params:
		sample="{sample}",
		workdir=WORKDIR,
		memG=MEMORYG,
		outdir=join(WORKDIR,"results","{sample}","STAR2p"),
		regions=config['regions']
	threads: 2
	envmodules: TOOLS["samtools"]["version"],TOOLS["bedtools"]["version"],TOOLS["ucsc"]["version"],TOOLS["sambamba"]["version"]
	shell:"""
cd {params.outdir}
bam_basename="$(basename {input.bam})"
while read a b;do
bam="${{bam_basename%.*}}.${{a}}.bam"
samtools view {input.bam} $b -b > /dev/shm/${{bam%.*}}.tmp.bam
sambamba sort --memory-limit={params.memG} --tmpdir=/dev/shm --nthreads={threads} --out=$bam /dev/shm/${{bam%.*}}.tmp.bam
bw="${{bam%.*}}.bw"
bdg="${{bam%.*}}.bdg"
sizes="${{bam%.*}}.sizes"
bedtools genomecov -bg -ibam $bam > $bdg
bedSort $bdg $bdg
if [ "$(wc -l $bdg|awk '{{print $1}}')" != "0" ];then
samtools view -H $bam|grep ^@SQ|cut -f2,3|sed "s/SN://g"|sed "s/LN://g" > $sizes
bedGraphToBigWig $bdg $sizes $bw
else
touch $bw
fi
rm -f $bdg $sizes
done < {params.regions}
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
		reffa=config['ref_fa']
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
		reffa=config['ref_fa'],
		bwaindex=config['ref_bwa_index'],
		gtf=config['ref_gtf'],
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
		matrix=join(WORKDIR,"results","ciri_count_matrix.txt")
	params:
		script=join(SCRIPTS_DIR,"Create_ciri_count_matrix.py"),
		lookup=join(RESOURCES_DIR,"hg38_2_hg19_lookup.txt"),
		outdir=join(WORKDIR,"results")
	envmodules: TOOLS["python37"]["version"]
	shell:"""
cd {params.outdir}
python {params.script} {params.lookup}
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
		lookup=join(RESOURCES_DIR,"hg38_2_hg19_lookup.txt"),
		outdir=join(WORKDIR,"results")
	envmodules: TOOLS["python37"]["version"]
	shell:"""
cd {params.outdir}
python {params.script} {params.lookup}
python {params.script2} {params.lookup}
"""

rule clear:
	input:
		bam=rules.star2p.output.bam,
		circexplorerout=rules.annotate_circRNA.output.annotations,
	output:
		outf1=join(WORKDIR,"results","{sample}","CLEAR","quant.txt")
	params:
		genepred=config['genepred_w_geneid'],
	container: "docker://nciccbr/ccbr_clear:latest"
	threads: 4
	shell:"""
circ_quant \
-c {input.circexplorerout} \
-b {input.bam} \
-t \
-r {params.genepred} \
-o {output.outf1}
"""

rule venn:
	input:
		circexplorerout=rules.annotate_circRNA.output.annotations,
		ciriout=rules.ciri.output.ciriout
	output:
		png=join(WORKDIR,"results","{sample}","{sample}.venn.png"),
		cirionly=join(WORKDIR,"results","{sample}","{sample}.cirionly.lst"),
		circexploreronly=join(WORKDIR,"results","{sample}","{sample}.circexploreronly.lst"),
		common=join(WORKDIR,"results","{sample}","{sample}.common.lst")
	params:
		script1=join(SCRIPTS_DIR,"venn.R"),
		sample="{sample}",
		outdir=join(WORKDIR,"results","{sample}")
	container: "docker://nciccbr/ccbr_venn:latest"
	threads: 2
	shell:"""
cut -f1 {input.ciriout}|grep -v circRNA_ID > /dev/shm/{params.sample}.ciri.lst
cut -f1-3 {input.circexplorerout}|awk -F"\\t" '{{print $1\":\"$2+1\"|\"$3}}' > /dev/shm/{params.sample}.circExplorer.lst
2set_venn.R \
	-l /dev/shm/{params.sample}.ciri.lst \
	-r /dev/shm/{params.sample}.circExplorer.lst  \
	-p {output.png} \
	-m {output.cirionly} \
	-s {output.circexploreronly} \
	-c1 "CIRI2" \
	-c2 "CircExplorer2" \
	-c {output.common} \
	-t {params.sample}
"""

rule multiqc:
	input:
		expand(join(WORKDIR,"results","{sample}","STAR2p","{sample}_p2.Aligned.sortedByCoord.out.bam"),sample=SAMPLES),
		expand(join(WORKDIR,"qc","picard_MarkDuplicates","{sample}.MarkDuplicates.metrics.txt"),sample=SAMPLES),
		expand(join(WORKDIR,"qc","fastqc","{sample}.R1.trim_fastqc.zip"),sample=SAMPLES),
		expand(join(WORKDIR,"results","{sample}","trim","{sample}.R1.trim.fastq.gz"),sample=SAMPLES)
	output:
		html=join(WORKDIR,"multiqc_report.html")
	params:
		outdir=WORKDIR
	envmodules: TOOLS["multiqc"]["version"]
	shell:"""
cd {params.outdir}
multiqc .
"""
