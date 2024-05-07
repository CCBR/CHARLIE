# AWK_SUM_READCOUNTS_CMD=r"""{sum=sum+$1}END{print sum}"""


def get_alignment_stats_input(wildcards):
    sample = wildcards.sample
    d = dict()
    d['star2bam']=join(WORKDIR,"results",sample,"STAR2p",sample+"_p2.bam")
    d['star2bam_chimeric']=join(WORKDIR,"results",sample,"STAR2p",sample+"_p2.chimeric.bam")
    d['star2bam_non_chimeric']=join(WORKDIR,"results",sample,"STAR2p",sample+"_p2.non_chimeric.bam")
    d['filtered_bam']=join(WORKDIR,"results",sample,"circExplorer",sample+".bam")  
    d['linearbam']=join(WORKDIR,"results",sample,"circExplorer",sample+".linear.bam")
    d['splicedbam']=join(WORKDIR,"results",sample,"circExplorer",sample+".spliced.bam")
    d['BSJbam']=join(WORKDIR,"results",sample,"circExplorer",sample+".BSJ.bam")
    d['ciribam']=join(WORKDIR,"results",sample,"ciri",sample+".ciri.bam")
    if RUN_MAPSPLICE:
        d["mapsplicebam"] = join(
            WORKDIR, "results", "{sample}", "MapSplice", "{sample}.mapsplice.cram"
        )
    return d


rule create_circExplorer_BSJ_bam:
    input:
        countstable=rules.circExplorer.output.annotation_counts_table,
        chimericbam=rules.star2p.output.chimeric_bam,
    output:
        BSJbam=join(WORKDIR, "results", "{sample}", "circExplorer", "{sample}.BSJ.bam"),
        plusBSJbam=join(
            WORKDIR, "results", "{sample}", "circExplorer", "{sample}.BSJ.plus.bam"
        ),
        minusBSJbam=join(
            WORKDIR, "results", "{sample}", "circExplorer", "{sample}.BSJ.minus.bam"
        ),
        BSJbed=join(
            WORKDIR, "results", "{sample}", "circExplorer", "{sample}.BSJ.bed.gz"
        ),
        BSJfoundcounts=join(
            WORKDIR,
            "results",
            "{sample}",
            "circExplorer",
            "{sample}.BSJ.foundcounts.tsv",
        ),
        BSJbw=join(WORKDIR, "results", "{sample}", "circExplorer", "{sample}.BSJ.bw"),
        plusBSJbw=join(
            WORKDIR, "results", "{sample}", "circExplorer", "{sample}.BSJ.plus.bw"
        ),
        minusBSJbw=join(
            WORKDIR, "results", "{sample}", "circExplorer", "{sample}.BSJ.minus.bw"
        ),
    params:
        sample="{sample}",
        refregions=REF_REGIONS,
        reffa=REF_FA,
        host=HOST,
        additives=ADDITIVES,
        viruses=VIRUSES,
        peorse=get_peorse,
        memG=getmemG("create_circExplorer_BSJ_bam"),
        scriptpe=join(SCRIPTS_DIR, "_create_circExplorer_BSJ_bam_pe.py"),
        scriptse=join(SCRIPTS_DIR, "_create_circExplorer_BSJ_bam_se.py"),
        flankscript=join(SCRIPTS_DIR, "_append_splice_site_flanks_to_BSJs.py"),
        bam2bwscript=join(SCRIPTS_DIR, "bam_to_bigwig.sh"),
        randomstr=str(uuid.uuid4()),
    container: config['containers']["star_ucsc_cufflinks"]
    threads: getthreads("create_circExplorer_BSJ_bam")
    shell:
        """
set -exo pipefail
if [ -d /lscratch/${{SLURM_JOB_ID}} ];then
    TMPDIR="/lscratch/${{SLURM_JOB_ID}}/{params.randomstr}"
else
    TMPDIR="/dev/shm/{params.randomstr}"
fi
if [ ! -d $TMPDIR ];then mkdir -p $TMPDIR;fi

outdir=$(dirname {output.BSJbam})
BSJbedbn=$(basename {output.BSJbed})

if [ "{params.peorse}" == "PE" ];then

python3 {params.scriptpe} \\
    --inbam {input.chimericbam} \\
    --sample_counts_table {input.countstable} \\
    --plusbam ${{TMPDIR}}/{params.sample}.BSJ.plus.bam \\
    --minusbam ${{TMPDIR}}/{params.sample}.BSJ.minus.bam \\
    --outbam ${{TMPDIR}}/{params.sample}.BSJ.bam \\
    --bed ${{TMPDIR}}/${{BSJbedbn}} \\
    --sample_name {params.sample} \\
    --junctionsfound {output.BSJfoundcounts} \\
    --regions {params.refregions} \\
    --host "{params.host}" \\
    --additives "{params.additives}" \\
    --viruses "{params.viruses}" \\
    --outputhostbams --outputvirusbams --outdir $TMPDIR

else

python3 {params.scriptse} \\
    --inbam {input.chimericbam} \\
    --sample_counts_table {input.countstable} \\
    --plusbam ${{TMPDIR}}/{params.sample}.BSJ.plus.bam \\
    --minusbam ${{TMPDIR}}/{params.sample}.BSJ.minus.bam \\
    --outbam ${{TMPDIR}}/{params.sample}.BSJ.bam \\
    --bed ${{TMPDIR}}/${{BSJbedbn}} \\
    --sample_name {params.sample} \\
    --junctionsfound {output.BSJfoundcounts} \\
    --regions {params.refregions} \\
    --host "{params.host}" \\
    --additives "{params.additives}" \\
    --viruses "{params.viruses}" \\
    --outputhostbams --outputvirusbams --outdir $TMPDIR

fi

samtools sort -l 9 -T $TMPDIR --write-index -@{threads} -O BAM -o {output.plusBSJbam} ${{TMPDIR}}/{params.sample}.BSJ.plus.bam 
samtools sort -l 9 -T $TMPDIR --write-index -@{threads} -O BAM -o {output.minusBSJbam} ${{TMPDIR}}/{params.sample}.BSJ.minus.bam 
samtools sort -l 9 -T $TMPDIR --write-index -@{threads} -O BAM -o {output.BSJbam} ${{TMPDIR}}/{params.sample}.BSJ.bam

for b in {output.plusBSJbam} {output.minusBSJbam} {output.BSJbam} 
# for b in {output.plusBSJbam} {output.minusBSJbam}
do
    bash {params.bam2bwscript} $b $TMPDIR
done

for i in $(echo {params.host}|tr ',' ' ');do 
    samtools sort -l 9 -T $TMPDIR --write-index -@{threads} -O BAM -o ${{outdir}}/{params.sample}.${{i}}.BSJ.bam ${{TMPDIR}}/{params.sample}.${{i}}.BSJ.bam
    bash {params.bam2bwscript} ${{outdir}}/{params.sample}.${{i}}.BSJ.bam $TMPDIR
done
for i in $(echo {params.viruses}|tr ',' ' ');do 
    samtools sort -l 9 -T $TMPDIR --write-index -@{threads} -O BAM -o ${{outdir}}/{params.sample}.${{i}}.BSJ.bam ${{TMPDIR}}/{params.sample}.${{i}}.BSJ.bam
    bash {params.bam2bwscript} ${{outdir}}/{params.sample}.${{i}}.BSJ.bam $TMPDIR
done

python3 {params.flankscript} --reffa {params.reffa} --inbsjbedgz ${{TMPDIR}}/${{BSJbedbn}} --outbsjbedgz {output.BSJbed}



rm -rf $TMPDIR
"""


rule create_circExplorer_linear_spliced_bams:
    input:
        junction=rules.star2p.output.junction,  # Chimeric junctions file
        bam=rules.star2p.output.bam,
        nonchimericbam=rules.star2p.output.non_chimeric_bam,
        bsjbedgz=rules.create_circExplorer_BSJ_bam.output.BSJbed,
        countstable=rules.circExplorer.output.annotation_counts_table,
    output:
        rid2jid=join(
            WORKDIR, "results", "{sample}", "circExplorer", "{sample}.rid2jid.tsv.gz"
        ),
        filtered_bam=join(
            WORKDIR, "results", "{sample}", "circExplorer", "{sample}.bam"
        ),
        linear_readids=join(
            WORKDIR,
            "results",
            "{sample}",
            "circExplorer",
            "{sample}.linear.BSJ.readids.gz",
        ),
        spliced_readids=join(
            WORKDIR,
            "results",
            "{sample}",
            "circExplorer",
            "{sample}.spliced.BSJ.readids.gz",
        ),
        linear_spliced_counts=join(
            WORKDIR,
            "results",
            "{sample}",
            "circExplorer",
            "{sample}.linear_spliced.counts.tsv",
        ),
        linear_BSJ_bam=join(
            WORKDIR, "results", "{sample}", "circExplorer", "{sample}.linear_BSJ.bam"
        ),  # linear reads in BSJ inclusion zone
        spliced_BSJ_bam=join(
            WORKDIR, "results", "{sample}", "circExplorer", "{sample}.spliced_BSJ.bam"
        ),  # linear spliced-only alignments in the sample
        linear_BSJ_bw=join(
            WORKDIR, "results", "{sample}", "circExplorer", "{sample}.linear_BSJ.bw"
        ),
        spliced_BSJ_bw=join(
            WORKDIR, "results", "{sample}", "circExplorer", "{sample}.spliced_BSJ.bw"
        ),
        linear_bam=join(
            WORKDIR, "results", "{sample}", "circExplorer", "{sample}.linear.bam"
        ),  # linear reads in BSJ inclusion zone
        spliced_bam=join(
            WORKDIR, "results", "{sample}", "circExplorer", "{sample}.spliced.bam"
        ),  # linear spliced-only alignments in the sample
        linear_bw=join(
            WORKDIR, "results", "{sample}", "circExplorer", "{sample}.linear.bw"
        ),
        spliced_bw=join(
            WORKDIR, "results", "{sample}", "circExplorer", "{sample}.spliced.bw"
        ),
    params:
        sample="{sample}",
        memG=getmemG("create_circExplorer_linear_spliced_bams"),
        refregions=REF_REGIONS,
        refregions_host=REF_REGIONS_HOST,
        refregions_viruses=REF_REGIONS_VIRUSES,
        host=HOST,
        additives=ADDITIVES,
        viruses=VIRUSES,
        peorse=get_peorse,
        bashscript=join(SCRIPTS_DIR, "_create_circExplorer_linear_bam.v2.sh"),
        outdir=join(WORKDIR, "results", "{sample}", "circExplorer"),
        randomstr=str(uuid.uuid4()),
    container: config['containers']["star_ucsc_cufflinks"]
    threads: getthreads("create_circExplorer_linear_spliced_bams")
    shell:
        """
set -exo pipefail
if [ -d /lscratch/${{SLURM_JOB_ID}} ];then
    TMPDIR="/lscratch/${{SLURM_JOB_ID}}/{params.randomstr}"
else
    TMPDIR="/dev/shm/{params.randomstr}"
fi
if [ -d $TMPDIR ];then rm -rf $TMPDIR;fi
mkdir -p $TMPDIR

cd {params.outdir}

# get filtered bam (remove secondary/supplementary/etc.) and the rid2jid lookup files

bash {params.bashscript} \\
    --nonchimericbam {input.nonchimericbam} \\
    --samplename {params.sample} \\
    --peorse {params.peorse} \\
    --bsjbed {input.bsjbedgz} \\
    --tmpdir $TMPDIR \\
    --rid2jid {output.rid2jid} \\
    --filteredbam {output.filtered_bam} \\
    --linearbsjlist {output.linear_readids} \\
    --splicedbsjlist {output.spliced_readids} \\
    --jidcounts {output.linear_spliced_counts} \\
    --linearbsjbam {output.linear_BSJ_bam} \\
    --splicedbsjbam {output.spliced_BSJ_bam} \\
    --regions {params.refregions} \\
    --host "{params.host}" \\
    --additives "{params.additives}" \\
    --viruses "{params.viruses}" \\
    --linearbam {output.linear_bam} \\
    --splicedbam {output.spliced_bam} \\
    --threads {threads}
rm -rf $TMPDIR
"""


# rule create_circExplorer_merged_found_counts_table:
# This rule  creates an output file with the following columns:
# chrom
# start
# end
# strand
# expected_BSJ_reads
# found_BSJ_reads
# linear_BSJ_reads_same_strand
# linear_spliced_BSJ_reads_same_strand
# linear_BSJ_reads_opposite_strand
# linear_spliced_BSJ_reads_opposite_strand
localrules:
    create_circExplorer_merged_found_counts_table,


rule create_circExplorer_merged_found_counts_table:
    input:
        annotation_counts=rules.circExplorer.output.annotation_counts_table,
        bsj_found_counts=rules.create_circExplorer_BSJ_bam.output.BSJfoundcounts,
        linear_spliced_counts=rules.create_circExplorer_linear_spliced_bams.output.linear_spliced_counts,
    output:
        found_counts_table=join(
            WORKDIR, "results", "{sample}", "circExplorer", "{sample}.readcounts.tsv"
        ),
        count_counts_table=join(
            WORKDIR,
            "results",
            "{sample}",
            "circExplorer",
            "{sample}.circExplorer.counts_table.tsv",
        ),
    params:
        sample="{sample}",
        pythonscript=join(SCRIPTS_DIR, "_merge_circExplorer_found_counts.py"),
        pythonscript2=join(
            SCRIPTS_DIR, "create_circExplorer_per_sample_counts_table.py"
        ),
        outdir=join(WORKDIR, "results", "{sample}", "circExplorer"),
        randomstr=str(uuid.uuid4()),
    shell:
        """
set -exo pipefail
if [ -d /lscratch/${{SLURM_JOB_ID}} ];then
    TMPDIR="/lscratch/${{SLURM_JOB_ID}}"
else
    TMPDIR="/dev/shm/{params.randomstr}"
    mkdir -p $TMPDIR
fi
python3 {params.pythonscript} \\
    -b {input.bsj_found_counts} \\
    -l {input.linear_spliced_counts} \\
    -o {output.found_counts_table}

python3 {params.pythonscript2} \\
    --annotationcounts {input.annotation_counts} \\
    --allfoundcounts {output.found_counts_table} \\
    --countstable {output.count_counts_table}
"""


# localrules: create_circExplorer_per_sample_counts_table
# rule create_circExplorer_per_sample_counts_table:
#     input:
#         annotation_counts=rules.circExplorer.output.annotation_counts_table,
#         found_counts=rules.create_circExplorer_merged_found_counts_table.output.found_counts_table
#     output:
#         count_counts_table=join(WORKDIR,"results","{sample}","circExplorer","{sample}.circExplorer.counts_table.tsv")
#     params:
#         sample="{sample}",
#         pythonscript=join(SCRIPTS_DIR,"create_circExplorer_per_sample_counts_table.py"),
#         outdir=join(WORKDIR,"results","{sample}","circExplorer"),
#         randomstr=str(uuid.uuid4())
#     shell:"""
# set -exo pipefail
# if [ -d /lscratch/${{SLURM_JOB_ID}} ];then
#     TMPDIR="/lscratch/${{SLURM_JOB_ID}}"
# else
#     TMPDIR="/dev/shm/{params.randomstr}"
#     mkdir -p $TMPDIR
# fi
# python3 {params.pythonscript} \\
#     --annotationcounts {input.annotation_counts} \\
#     --allfoundcounts {input.found_counts} \\
#     --countstable {output.count_counts_table}
# """

if RUN_MAPSPLICE:

    rule alignment_stats:
        input:
            unpack(get_alignment_stats_input),
            # star2bam=rules.star2p.output.bam,
            # splicedbam=rules.create_circExplorer_linear_spliced_bams.output.spliced_bam,
            # linearbam=rules.create_circExplorer_linear_spliced_bams.output.linear_bam,
            # linearsplicedbam=rules.create_circExplorer_linear_spliced_bams.output.linear_spliced_bam,
            # BSJbam=rules.create_circExplorer_BSJ_bam.output.BSJbam,
            # ciribam=rules.ciri.output.ciribam,
        output:
            alignmentstats=join(WORKDIR, "results", "{sample}", "alignmentstats.txt"),
        params:
            sample="{sample}",
            regions=REF_REGIONS,
            peorse=get_peorse,
            run_mapsplice=N_RUN_MAPSPLICE,
            bash2nreads_pyscript=join(SCRIPTS_DIR, "_bam_get_alignment_stats.py"),
            randomstr=str(uuid.uuid4()),
        threads: getthreads("alignment_stats")
        container: config['containers']["base"]
        shell:
            """
    set -exo pipefail
    if [ -d /lscratch/${{SLURM_JOB_ID}} ];then
        TMPDIR="/lscratch/${{SLURM_JOB_ID}}"
    else
        TMPDIR="/dev/shm/{params.randomstr}"
        mkdir -p $TMPDIR
    fi
    for bamfile in {input};do
        bamfile_bn=$(basename $bamfile)
        if [ "{params.peorse}" == "PE" ];then
        echo "python3 {params.bash2nreads_pyscript} --inbam $bamfile --regions {params.regions} --pe > ${{TMPDIR}}/${{bamfile_bn}}.counts"
        else
        echo "python3 {params.bash2nreads_pyscript} --inbam $bamfile --regions {params.regions} > ${{TMPDIR}}/${{bamfile_bn}}.counts"
        fi
    done > ${{TMPDIR}}/do_bamstats
    parallel -j 2 < ${{TMPDIR}}/do_bamstats

    print_bam_results () {{
        bamfile=$1
        bamfile_bn=$(basename $bamfile)
        stats_file=${{TMPDIR}}/${{bamfile_bn}}.counts
        prefix=$2
        while read b a;do echo -ne "${{prefix}}_${{a}}\\t${{b}}\\n";done < $stats_file
    }}

    echo -ne "sample\\t{params.sample}\\n" > {output.alignmentstats}
    print_bam_results {input.star2bam} "STAR" >> {output.alignmentstats}
    print_bam_results {input.star2bam_chimeric} "STAR_chimeric" >> {output.alignmentstats}
    print_bam_results {input.star2bam_non_chimeric} "STAR_non_chimeric" >> {output.alignmentstats}
    # print_bam_results {input.filtered_bam} "STAR_non_chimeric_filtered" >> {output.alignmentstats}
    print_bam_results {input.linearbam} "CircExplorer_linear" >> {output.alignmentstats}
    print_bam_results {input.splicedbam} "CircExplorer_spliced" >> {output.alignmentstats}
    print_bam_results {input.BSJbam} "CircExplorer_BSJ" >> {output.alignmentstats}
    print_bam_results {input.ciribam} "CIRI" >> {output.alignmentstats}
    print_bam_results {input.mapsplicebam} "MapSplice" >> {output.alignmentstats}
    """

else:

    rule alignment_stats:
        input:
            unpack(get_alignment_stats_input),
        output:
            alignmentstats=join(WORKDIR, "results", "{sample}", "alignmentstats.txt"),
        params:
            sample="{sample}",
            regions=REF_REGIONS,
            peorse=get_peorse,
            run_mapsplice=N_RUN_MAPSPLICE,
            bash2nreads_pyscript=join(SCRIPTS_DIR, "_bam_get_alignment_stats.py"),
            randomstr=str(uuid.uuid4()),
        threads: getthreads("alignment_stats")
        container: config['containers']["base"]
        shell:
            """
    set -exo pipefail
    if [ -d /lscratch/${{SLURM_JOB_ID}} ];then
        TMPDIR="/lscratch/${{SLURM_JOB_ID}}"
    else
        TMPDIR="/dev/shm/{params.randomstr}"
        mkdir -p $TMPDIR
    fi
    for bamfile in {input};do
        bamfile_bn=$(basename $bamfile)
        if [ "{params.peorse}" == "PE" ];then
        echo "python3 {params.bash2nreads_pyscript} --inbam $bamfile --regions {params.regions} --pe > ${{TMPDIR}}/${{bamfile_bn}}.counts"
        else
        echo "python3 {params.bash2nreads_pyscript} --inbam $bamfile --regions {params.regions} > ${{TMPDIR}}/${{bamfile_bn}}.counts"
        fi
    done > ${{TMPDIR}}/do_bamstats
    parallel -j 2 < ${{TMPDIR}}/do_bamstats

    print_bam_results () {{
        bamfile=$1
        bamfile_bn=$(basename $bamfile)
        stats_file=${{TMPDIR}}/${{bamfile_bn}}.counts
        prefix=$2
        while read b a;do echo -ne "${{prefix}}_${{a}}\\t${{b}}\\n";done < $stats_file
    }}

    echo -ne "sample\\t{params.sample}\\n" > {output.alignmentstats}
    print_bam_results {input.star2bam} "STAR" >> {output.alignmentstats}
    print_bam_results {input.star2bam_chimeric} "STAR_chimeric" >> {output.alignmentstats}
    print_bam_results {input.star2bam_non_chimeric} "STAR_non_chimeric" >> {output.alignmentstats}
    print_bam_results {input.filtered_bam} "STAR_filtered" >> {output.alignmentstats}
    print_bam_results {input.linearbam} "CircExplorer_linear" >> {output.alignmentstats}
    print_bam_results {input.splicedbam} "CircExplorer_spliced" >> {output.alignmentstats}
    print_bam_results {input.BSJbam} "CircExplorer_BSJ" >> {output.alignmentstats}
    print_bam_results {input.ciribam} "CIRI" >> {output.alignmentstats}

    """


localrules:
    merge_alignment_stats,


rule merge_alignment_stats:
    input:
        expand(
            join(WORKDIR, "results", "{sample}", "alignmentstats.txt"), sample=SAMPLES
        ),
    output:
        join(WORKDIR, "results", "alignmentstats.txt"),
    params:
        randomstr=str(uuid.uuid4()),
    shell:
        """
set -exo pipefail
if [ -d /lscratch/${{SLURM_JOB_ID}} ];then
    TMPDIR="/lscratch/${{SLURM_JOB_ID}}"
else
    TMPDIR="/dev/shm/{params.randomstr}"
    mkdir -p $TMPDIR
fi
count=0
for f in {input};do
    count=$((count+1))
    if [ "$count" == "1" ];then
        cp $f {output}
    else
        cut -f2 $f > ${{TMPDIR}}/${{count}}
        paste {output} ${{TMPDIR}}/${{count}} > ${{TMPDIR}}/${{count}}.tmp
        mv ${{TMPDIR}}/${{count}}.tmp {output}
    fi
done 
"""


rule create_hq_bams:
    input:
        inbam=rules.create_circExplorer_BSJ_bam.output.BSJbam,
        countstable=rules.create_master_counts_file.output.matrix,
    output:
        outbam=join(WORKDIR, "results", "HQ_BSJ_bams","{sample}.HQ_only.BSJ.bam")
    params:
        script=join(SCRIPTS_DIR,"_bam_filter_BSJ_for_HQonly.py"),
        samplename="{sample}",
        regions=REF_REGIONS,
        host=HOST,
        additives=ADDITIVES,
        viruses=VIRUSES,
    container: config['containers']["base"]
    shell:"""
set -exo pipefail
outdir=$(dirname {output.outbam})
if [ ! -d $outdir ];then mkdir -p $outdir;fi
cd $outdir
python3 {params.script} \\
    -i {input.inbam} \\
    -t {input.countstable} \\
    -o {output.outbam} \\
    --regions {params.regions} \\
    --host "{params.host}" \\
    --additives "{params.additives}" \\
    --viruses "{params.viruses}" \\
    --sample_name {params.samplename}
samtools index {output.outbam}
for bam in $(ls {params.samplename}.*.HQ_only.BSJ.bam);do
    if [ ! -f "${{bam}}.bai" ];then
        samtools index $bam
    fi
done
"""
