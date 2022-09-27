## functions

def get_fastqs(wildcards):
	d=dict()
	d["R1"]=SAMPLESDF["R1"][wildcards.sample]
	d["R2"]=SAMPLESDF["R2"][wildcards.sample]
	return d

## rules

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
		cutadapt_min_length=config["cutadapt_min_length"],
		adapters=join(RESOURCES_DIR,"TruSeq_and_nextera_adapters.consolidated.fa")
	envmodules: TOOLS["cutadapt"]["version"]
	threads: getthreads("cutadapt")
	shell:"""
set -exo pipefail
if [ ! -d {params.outdir} ];then mkdir {params.outdir};fi
if [ "{params.peorse}" == "PE" ];then
	## Paired-end
	cutadapt --pair-filter=any \\
	--nextseq-trim=2 \\
	--trim-n \\
	-n 5 -O 5 \\
	-q 10,10 -m {params.cutadapt_min_length}:{params.cutadapt_min_length} \\
	-b file:{params.adapters} \\
	-B file:{params.adapters} \\
	-j {threads} \\
	-o {output.of1} -p {output.of2} \\
	{input.R1} {input.R2}
else
	## Single-end
	cutadapt \\
	--nextseq-trim=2 \\
	--trim-n \\
	-n 5 -O 5 \\
	-q 10,10 -m {params.cutadapt_min_length} \\
	-b file:{params.adapters} \\
	-j {threads} \\
	-o {output.of1} \\
	{input.R1}
	touch {output.of2}
fi
"""