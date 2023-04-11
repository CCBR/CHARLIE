## functions


def get_raw_and_trim_fastqs(wildcards):
    d = dict()
    d["R1"] = SAMPLESDF["R1"][wildcards.sample]
    d["R2"] = SAMPLESDF["R2"][wildcards.sample]
    d["R1trim"] = join(
        WORKDIR,
        "results",
        wildcards.sample,
        "trim",
        wildcards.sample + ".R1.trim.fastq.gz",
    )
    d["R2trim"] = join(
        WORKDIR,
        "results",
        wildcards.sample,
        "trim",
        wildcards.sample + ".R2.trim.fastq.gz",
    )
    return d


## rules


rule fastqc:
    input:
        unpack(get_raw_and_trim_fastqs),
    output:
        join(WORKDIR, "qc", "fastqc", "{sample}.R1.trim_fastqc.zip"),
    params:
        outdir=join(WORKDIR, "qc", "fastqc"),
    threads: getthreads("fastqc")
    envmodules:
        TOOLS["fastqc"]["version"],
    shell:
        """
set -exo pipefail
files=""
for f in {input};do
if [ "$(wc $f|awk '{{print $1}}')" != "0" ];then
files=$(echo -ne "$files $f")
fi
done
fastqc $files -t {threads} -o {params.outdir}
"""


localrules:
    multiqc,


rule multiqc:
    input:
        expand(
            join(
                WORKDIR,
                "results",
                "{sample}",
                "STAR2p",
                "{sample}_p2.Aligned.sortedByCoord.out.bam",
            ),
            sample=SAMPLES,
        ),
        expand(
            join(
                WORKDIR,
                "qc",
                "picard_MarkDuplicates",
                "{sample}.MarkDuplicates.metrics.txt",
            ),
            sample=SAMPLES,
        ),
        expand(
            join(WORKDIR, "qc", "fastqc", "{sample}.R1.trim_fastqc.zip"),
            sample=SAMPLES,
        ),
        expand(
            join(WORKDIR, "results", "{sample}", "trim", "{sample}.R1.trim.fastq.gz"),
            sample=SAMPLES,
        ),
    output:
        html=join(WORKDIR, "multiqc_report.html"),
    params:
        outdir=WORKDIR,
    envmodules:
        TOOLS["multiqc"]["version"],
    shell:
        """
cd {params.outdir}
multiqc --force --ignore *.snakemake --ignore *multiqc* --verbose .
"""
