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
		cutadapt_n=config['cutadapt_n'],
		cutadapt_max_n=config['cutadapt_max_n'],
		cutadapt_O=config['cutadapt_O'],
		cutadapt_q=config['cutadapt_q'],
		adapters=join(RESOURCES_DIR,"TruSeq_and_nextera_adapters.consolidated.fa"),
		randomstr=str(uuid.uuid4()),
	envmodules: TOOLS["cutadapt"]["version"]
	threads: getthreads("cutadapt")
	shell:"""
set -exo pipefail

# set TMPDIR
if [ -d /lscratch/${{SLURM_JOB_ID}} ];then
    TMPDIR="/lscratch/${{SLURM_JOB_ID}}/{params.randomstr}"
else
    TMPDIR="/dev/shm/{params.randomstr}"
fi
if [ ! -d $TMPDIR ];then mkdir -p $TMPDIR;fi

# set conda environment
. "/data/CCBR_Pipeliner/db/PipeDB/Conda/etc/profile.d/conda.sh"
conda activate fastqfilter

if [ ! -d {params.outdir} ];then mkdir {params.outdir};fi

of1bn=$(basename {output.of1})
of2bn=$(basename {output.of2})

if [ "{params.peorse}" == "PE" ];then
	## Paired-end
	cutadapt --pair-filter=any \\
	--nextseq-trim=2 \\
	--trim-n \\
	--max-n {params.cutadapt_max_n} \\
	-n {params.cutadapt_n} -O {params.cutadapt_O} \\
	-q {params.cutadapt_q},{params.cutadapt_q} -m {params.cutadapt_min_length}:{params.cutadapt_min_length} \\
	-b file:{params.adapters} \\
	-B file:{params.adapters} \\
	-j {threads} \\
	-o ${{TMPDIR}}/${{of1bn}} -p ${{TMPDIR}}/${{of2bn}} \\
	{input.R1} {input.R2}

	fastq-filter \\
		-q {params.cutadapt_q} \\
		-o {output.of1} -o {output.of2} \\
		${{TMPDIR}}/${{of1bn}} ${{TMPDIR}}/${{of2bn}}

else
	## Single-end
	cutadapt \\
	--nextseq-trim=2 \\
	--trim-n \\
	--max-n {params.cutadapt_max_n} \\
	-n {params.cutadapt_n} -O {params.cutadapt_O} \\
	-q {params.cutadapt_q},{params.cutadapt_q} -m {params.cutadapt_min_length} \\
	-b file:{params.adapters} \\
	-j {threads} \\
	-o ${{TMPDIR}}/${{of1bn}} \\
	{input.R1}
	
	touch {output.of2}

	fastq-filter \\
		-q {params.cutadapt_q} \\
		-o {output.of1} \\
		${{TMPDIR}}/${{of1bn}}

fi

# filter for average read quality


"""