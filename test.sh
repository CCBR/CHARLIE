#!/usr/bin/env bash
set -euo pipefail

# ## setting PIPELINE_HOME
# ## clone the pipeline to a folder
# ## git clone https://github.com/kopardev/circRNA.git
# ## and set that as the pipeline home
PIPELINE_HOME="/data/Ziegelbauer_lab/circRNADetection/scripts/circRNA"
# ## make current folder as working directory
a=$(readlink -f $0)
# echo $a
for i in `seq 1 20`;do
b=$(echo $a|sed "s/\/gpfs\/gsfs${i}\/users/\/data/g")
a=$b
done
# echo $a
WORKDIR=$(dirname $a)

function usage() { cat << EOF
test.sh: test the workflow.
USAGE:
  bash test.sh <MODE> 
Required Positional Argument:
  MODE: [Type: Str] Valid options:
    a) init: initial workdir
    b) run: run test data
EOF
}

function err() { cat <<< "$@" 1>&2; }

function init() {
echo "Working Dir: $WORKDIR"
cd $WORKDIR

# make sure that config folder exists in the current folder
if [ ! -d $WORKDIR/config ]; then cp -r $PIPELINE_HOME/config $WORKDIR/;echo "Config Dir: $WORKDIR/config";fi

#create log and stats folders
if [ ! -d $WORKDIR/logs ]; then mkdir -p $WORKDIR/logs;echo "Logs Dir: $WORKDIR/logs";fi
if [ ! -d $WORKDIR/stats ];then mkdir -p $WORKDIR/stats;echo "Stats Dir: $WORKDIR/stats";fi
}

function run () {
  ## initialize if not already done
  if [ -d $WORKDIR/config ]; then err "Error: config folder not found ... initialize first!";usage && exit 1;fi
  ## Archive previous run files
  if [ -f ${WORKDIR}/snakemake.log ];then 
    modtime=$(stat ${WORKDIR}/snakemake.log |grep Modify|awk '{print $2,$3}'|awk -F"." '{print $1}'|sed "s/ //g"|sed "s/-//g"|sed "s/://g")
    mv ${WORKDIR}/snakemake.log ${WORKDIR}/stats/snakemake.log.${modtime} && gzip -n ${WORKDIR}/stats/snakemake.log.${modtime}
  fi
  if [ -f ${WORKDIR}/snakemake.stats ];then 
    modtime=$(stat ${WORKDIR}/snakemake.stats |grep Modify|awk '{print $2,$3}'|awk -F"." '{print $1}'|sed "s/ //g"|sed "s/-//g"|sed "s/://g")
    mv ${WORKDIR}/snakemake.stats ${WORKDIR}/stats/snakemake.stats.${modtime} && gzip -n ${WORKDIR}/stats/snakemake.stats.${modtime}
  fi


  module load python/3.7
  module load snakemake/5.24.1
  module load singularity

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

  mv slurm*out stats && for a in $(ls stats/slurm*out);do gzip -n $a;done
}

function main(){

  if [ $# -eq 0 ]; then usage; exit 1; fi

  case $1 in
    init) init && exit 0;;
    run) run && exit 0;;
    -h | --help | help) usage && exit 0;;
    -* | --*) err "Error: Failed to provide mode: <init|run>."; usage && exit 1;;
    *) err "Error: Failed to provide mode: <init|run>. '${1}' is not supported."; usage && exit 1;;
  esac
}

main "$@"





