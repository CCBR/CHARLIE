#!/usr/bin/env bash
#$ -N charlie
#$ -cwd
#$ -l h_rt=48:00:00
#$ -l h_vmem=40g
#$ -pe sharedmem 2
#$ -o logs/$JOB_NAME-$JOB_ID-$HOSTNAME.out
#$ -e logs/$JOB_NAME-$JOB_ID-$HOSTNAME.err

$MODULE_LOAD
$EXPORT_SING_CACHE_DIR_CMD

snakemake -s $SNAKEFILE \
    --directory $WORKDIR \
    --use-singularity \
    --singularity-args "$SINGULARITY_BINDS" \
    --use-envmodules \
    --printshellcmds \
    --latency-wait 300 \
    --configfile $CONFIGFILE \
    --profile $CLUSTER_PROFILE \
    -j 500 \
    --rerun-incomplete \
    --rerun-triggers $trigger \
    --retries 2 \
    --keep-going \
    --stats ${WORKDIR}/snakemake.stats \
    2>&1 | tee ${WORKDIR}/snakemake.log

if [ "$?" -eq "0" ];then
  snakemake -s $SNAKEFILE \
  --directory $WORKDIR \
  --report ${WORKDIR}/runqsub_snakemake_report.html \
  --configfile $CONFIGFILE
fi
