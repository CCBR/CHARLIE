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
		starindexdir=STAR_INDEX_DIR,
		alignTranscriptsPerReadNmax=TOOLS["star"]["alignTranscriptsPerReadNmax"],
		gtf=REF_GTF,
		randomstr=str(uuid.uuid4())
	envmodules: TOOLS["star"]["version"]
	threads: getthreads("star1p")
	shell:"""
set -exo pipefail
if [ -d /lscratch/${{SLURM_JOB_ID}} ];then
	TMPDIR="/lscratch/${{SLURM_JOB_ID}}/{params.randomstr}"
else
	TMPDIR="/dev/shm/{params.randomstr}"
fi

if [ ! -d {params.outdir} ];then mkdir {params.outdir};fi
if [ "{params.peorse}" == "PE" ];then
# paired-end
	overhang=$(zcat {input} | awk -v maxlen=100 'NR%4==2 {{if (length($1) > maxlen+0) maxlen=length($1)}}; END {{print maxlen-1}}')
	echo "sjdbOverhang for STAR: ${{overhang}}"
	cd {params.outdir}
	STAR --genomeDir {params.starindexdir} \\
	--outSAMstrandField None  \\
	--outFilterMultimapNmax 20 \\
	--alignSJoverhangMin 8 \\
	--alignSJDBoverhangMin 1 \\
	--outFilterMismatchNmax 999 \\
	--outFilterMismatchNoverLmax 0.3  \\
	--alignIntronMin 20 \\
	--alignIntronMax 1000000 \\
	--alignMatesGapMax 1000000 \\
	--readFilesIn {input.R1} {input.R2} \\
	--readFilesCommand zcat \\
	--runThreadN {threads} \\
	--outFileNamePrefix {params.sample}_p1. \\
	--chimSegmentMin 20 \\
	--chimMultimapNmax 10 \\
	--chimOutType Junctions \\
	--alignTranscriptsPerReadNmax {params.alignTranscriptsPerReadNmax} \\
	--outSAMtype None \\
	--alignEndsProtrude 10 ConcordantPair \\
	--outFilterIntronMotifs None \\
	--sjdbGTFfile {params.gtf} \\
	--outTmpDir=${{TMPDIR}} \\
	--sjdbOverhang $overhang

else

#single-end
	overhang=$(zcat {input.R1} | awk -v maxlen=100 'NR%4==2 {{if (length($1) > maxlen+0) maxlen=length($1)}}; END {{print maxlen-1}}')
	echo "sjdbOverhang for STAR: ${{overhang}}"
	cd {params.outdir}
	STAR --genomeDir {params.starindexdir} \\
	--outSAMstrandField None  \\
	--outFilterMultimapNmax 20 \\
	--alignSJoverhangMin 8 \\
	--alignSJDBoverhangMin 1 \\
	--outFilterMismatchNmax 999 \\
	--outFilterMismatchNoverLmax 0.3  \\
	--alignIntronMin 20 \\
	--alignIntronMax 1000000 \\
	--alignMatesGapMax 1000000 \\
	--readFilesIn {input.R1} \\
	--readFilesCommand zcat \\
	--runThreadN {threads} \\
	--outFileNamePrefix {params.sample}_p1. \\
	--chimSegmentMin 20 \\
	--chimMultimapNmax 10 \\
	--chimOutType Junctions \\
	--alignTranscriptsPerReadNmax {params.alignTranscriptsPerReadNmax} \\
	--outSAMtype None \\
	--alignEndsProtrude 10 ConcordantPair \\
	--outFilterIntronMotifs None \\
	--sjdbGTFfile {params.gtf} \\
	--outTmpDir=${{TMPDIR}} \\
	--sjdbOverhang $overhang

fi
rm -rf {params.outdir}/{params.sample}_p1.Aligned.out.bam
"""

rule merge_SJ_tabs:
	input:
		expand(join(WORKDIR,"results","{sample}","STAR1p","{sample}_p1.SJ.out.tab"), sample=SAMPLES)
	output:
		pass1sjtab=join(WORKDIR,"results","pass1.out.tab")
	params:
		script1=join(SCRIPTS_DIR,"apply_junction_filters.py"),
		regions=REF_REGIONS,
		filter1regions=HOST_ADDITIVES,
		filter1_noncanonical=config['star_1pass_filter_host_noncanonical'],
		filter1_unannotated=config['star_1pass_filter_host_unannotated'],
		filter2_noncanonical=config['star_1pass_filter_viruses_noncanonical'],
		filter2_unannotated=config['star_1pass_filter_viruses_unannotated']
	threads: getthreads("merge_SJ_tabs")
	shell:"""
set -exo pipefail
cat {input} | \\
python {params.script1} \\
	--regions {params.regions} \\
	--filter1regions {params.filter1regions} \\
	--filter1_noncanonical {params.filter1_noncanonical} \\
	--filter1_unannotated {params.filter1_unannotated} \\
	--filter2_noncanonical {params.filter2_noncanonical} \\
	--filter2_unannotated {params.filter1_unannotated} | \\
cut -f1-4 | sort -k1,1 -k2,2n | uniq > {output.pass1sjtab}
"""

rule star2p:
	input:
		R1=rules.cutadapt.output.of1,
		R2=rules.cutadapt.output.of2,
		pass1sjtab=rules.merge_SJ_tabs.output.pass1sjtab
	output:
		junction=join(WORKDIR,"results","{sample}","STAR2p","{sample}_p2.Chimeric.out.junction"),
		bam=join(WORKDIR,"results","{sample}","STAR2p","{sample}_p2.Aligned.sortedByCoord.out.bam"),
		genecounts=join(WORKDIR,"results","{sample}","STAR2p","{sample}_p2.ReadsPerGene.out.tab")
	params:
		sample="{sample}",
		peorse=get_peorse,
		workdir=WORKDIR,
		outdir=join(WORKDIR,"results","{sample}","STAR2p"),
		starindexdir=STAR_INDEX_DIR,
		alignTranscriptsPerReadNmax=TOOLS["star"]["alignTranscriptsPerReadNmax"],
		gtf=REF_GTF,
		randomstr=str(uuid.uuid4())
	envmodules: TOOLS["star"]["version"],TOOLS["sambamba"]["version"], TOOLS["samtools"]["version"]
	threads: getthreads("star2p")
	shell:"""
set -exo pipefail
if [ -d /lscratch/${{SLURM_JOB_ID}} ];then
	TMPDIR="/lscratch/${{SLURM_JOB_ID}}/{params.randomstr}"
else
	TMPDIR="/dev/shm/{params.randomstr}"
fi

if [ ! -d {params.outdir} ];then mkdir {params.outdir};fi
limitSjdbInsertNsj=$(wc -l {input.pass1sjtab}|awk '{{print $1+1}}')
if [ "$limitSjdbInsertNsj" -lt "400000" ];then limitSjdbInsertNsj="400000";fi
if [ "{params.peorse}" == "PE" ];then
# paired-end
	overhang=$(zcat {input.R1} {input.R2} | awk -v maxlen=100 'NR%4==2 {{if (length($1) > maxlen+0) maxlen=length($1)}}; END {{print maxlen-1}}')
	echo "sjdbOverhang for STAR: ${{overhang}}"
	cd {params.outdir}
	STAR --genomeDir {params.starindexdir} \\
	--outSAMstrandField None  \\
	--outFilterType BySJout \\
	--outFilterMultimapNmax 20 \\
	--alignSJoverhangMin 8 \\
	--alignSJDBoverhangMin 1 \\
	--outFilterMismatchNmax 999 \\
	--outFilterMismatchNoverLmax 0.3  \\
	--alignIntronMin 20 \\
	--alignIntronMax 2000000 \\
	--alignMatesGapMax 2000000 \\
	--readFilesIn {input.R1} {input.R2} \\
	--readFilesCommand  zcat \\
	--runThreadN 56 \\
	--outFileNamePrefix {params.sample}_p2. \\
	--sjdbFileChrStartEnd {input.pass1sjtab} \\
	--chimSegmentMin 20 \\
	--chimOutType Junctions WithinBAM \\
	--chimMultimapNmax 10 \\
	--limitSjdbInsertNsj $limitSjdbInsertNsj \\
	--alignTranscriptsPerReadNmax {params.alignTranscriptsPerReadNmax} \\
	--outSAMtype BAM SortedByCoordinate \\
	--alignEndsProtrude 10 ConcordantPair \\
	--outFilterIntronMotifs None \\
	--sjdbGTFfile {params.gtf} \\
	--quantMode GeneCounts \\
	--outTmpDir=${{TMPDIR}} \\
	--sjdbOverhang $overhang

else
#single-end
	overhang=$(zcat {input.R1} | awk -v maxlen=100 'NR%4==2 {{if (length($1) > maxlen+0) maxlen=length($1)}}; END {{print maxlen-1}}')
	echo "sjdbOverhang for STAR: ${{overhang}}"
	cd {params.outdir}
	STAR --genomeDir {params.starindexdir} \\
	--outSAMstrandField None  \\
	--outFilterType BySJout \\
	--outFilterMultimapNmax 20 \\
	--alignSJoverhangMin 8 \\
	--alignSJDBoverhangMin 1 \\
	--outFilterMismatchNmax 999 \\
	--outFilterMismatchNoverLmax 0.3  \\
	--alignIntronMin 20 \\
	--alignIntronMax 2000000 \\
	--alignMatesGapMax 2000000 \\
	--readFilesIn {input.R1} \\
	--readFilesCommand  zcat \\
	--runThreadN 56 \\
	--outFileNamePrefix {params.sample}_p2. \\
	--sjdbFileChrStartEnd {input.pass1sjtab} \\
	--chimSegmentMin 20 \\
	--chimOutType Junctions WithinBAM \\
	--chimMultimapNmax 10 \\
	--limitSjdbInsertNsj $limitSjdbInsertNsj \\
	--alignTranscriptsPerReadNmax {params.alignTranscriptsPerReadNmax} \\
	--outSAMtype BAM SortedByCoordinate \\
	--alignEndsProtrude 10 ConcordantPair \\
	--outFilterIntronMotifs None \\
	--sjdbGTFfile {params.gtf} \\
	--quantMode GeneCounts \\
	--outTmpDir=${{TMPDIR}} \\
	--sjdbOverhang $overhang
fi
## ensure the star2p file is indexed ... is should already be sorted by STAR
sleep 120
samtools index {output.bam}
"""

rule estimate_duplication:
	input:
		bam=rules.star2p.output.bam
	output:
		metrics=join(WORKDIR,"qc","picard_MarkDuplicates","{sample}.MarkDuplicates.metrics.txt")
	params:
		sample="{sample}",
		memG=getmemG("estimate_duplication"),
	envmodules: TOOLS["picard"]["version"]
	shell:"""
set -exo pipefail
java -Xmx{params.memG} -jar ${{PICARD_JARPATH}}/picard.jar MarkDuplicates I={input.bam} O=/dev/shm/{params.sample}.mark_dup.bam M={output.metrics}
rm -f /dev/shm/{params.sample}*
"""