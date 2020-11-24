#!/bin/bash
module load python/3.7
module load snakemake/5.24.1
module load singularity

PIPELINE_HOME="/data/RBL_NCI/ONT/ont_rnaseq/workflow/ONT_RNAseq"
# make current folder as working directory
a=$(readlink -f $0)
echo $a
for i in `seq 1 20`;do
b=$(echo $a|sed "s/\/gpfs\/gsfs${i}\/users/\/data/g")
a=$b
done
echo $a

WORKDIR=$(dirname $a)
if [ ! -d ${WORKDIR}/.singularity_cache ];then mkdir ${WORKDIR}/.singularity_cache;fi
export SINGULARITY_CACHEDIR="${WORKDIR}/.singularity_cache"
echo "Working Dir: $WORKDIR"
cd $WORKDIR

# make sure that config folder exists in the current folder
# if [ -d $WORKDIR/config ]; then rm -rf $WORKDIR/config;fi
# cp -r $PIPELINE_HOME/config $WORKDIR/
# cp -r $PIPELINE_HOME/workflow/envs $WORKDIR/envs
if [ ! -d $WORKDIR/logs/cluster ]; then mkdir -p $WORKDIR/logs/cluster;fi

if [ ! -d slurmfiles ];then mkdir -p slurmfiles ;fi

if [ -f ${WORKDIR}/snakemake.log ];then 
  modtime=$(stat ${WORKDIR}/snakemake.log |grep Modify|awk '{print $2,$3}'|awk -F"." '{print $1}'|sed "s/ //g"|sed "s/-//g"|sed "s/://g")
  mv ${WORKDIR}/snakemake.log ${WORKDIR}/slurmfiles/snakemake.log.${modtime} && gzip -n ${WORKDIR}/slurmfiles/snakemake.log.${modtime}
fi
if [ -f ${WORKDIR}/snakemake.stats ];then 
  modtime=$(stat ${WORKDIR}/snakemake.stats |grep Modify|awk '{print $2,$3}'|awk -F"." '{print $1}'|sed "s/ //g"|sed "s/-//g"|sed "s/://g")
  mv ${WORKDIR}/snakemake.stats ${WORKDIR}/slurmfiles/snakemake.stats.${modtime} && gzip -n ${WORKDIR}/slurmfiles/snakemake.stats.${modtime}
fi

# --use-conda \

snakemake $@ -s ${PIPELINE_HOME}/workflow/Snakefile \
--directory $WORKDIR \
--use-envmodules \
--printshellcmds \
--use-singularity \
--singularity-args "-B ${WORKDIR}" \
--latency-wait 120 \
--configfile ${WORKDIR}/config/config.yaml \
--cluster-config ${PIPELINE_HOME}/config/cluster.json \
--cluster "sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name {cluster.name} --output {cluster.output} --error {cluster.error}" \
-j 500 \
--rerun-incomplete \
--keep-going \
--stats ${WORKDIR}/snakemake.stats \
2>&1|tee ${WORKDIR}/snakemake.log

mv slurm*out slurmfiles && for a in $(ls slurmfiles/slurm*out);do gzip -n $a;done
