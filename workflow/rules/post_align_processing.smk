
localrules:
    merge_genecounts,


rule merge_genecounts:
    input:
        expand(
            join(
                WORKDIR,
                "results",
                "{sample}",
                "STAR2p",
                "{sample}_p2.ReadsPerGene.out.tab",
            ),
            sample=SAMPLES,
        ),
    output:
        join(WORKDIR, "results", "unstranded_STAR_GeneCounts.tsv"),
        join(WORKDIR, "results", "stranded_STAR_GeneCounts.tsv"),
        join(WORKDIR, "results", "revstranded_STAR_GeneCounts.tsv"),
    params:
        outdir=join(WORKDIR, "results"),
        rscript=join(SCRIPTS_DIR, "merge_ReadsPerGene_counts.R"),
    container: config['containers']["R"]
    shell:
        """
set -exo pipefail
cd {params.outdir}
Rscript {params.rscript}
"""
