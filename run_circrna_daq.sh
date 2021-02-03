#!/usr/bin/env bash
# cirna_daq
# d=detection
# a=annotation
# q=quantification
set -euo pipefail
module purge

function get_git_commitid_tag() {
  cd $1
  gid=$(git rev-parse HEAD)
  tag=$(git describe --tags $gid)
  echo -ne "$gid\t$tag"
}

# ## setting PIPELINE_HOME
# ## clone the pipeline to a folder
# ## git clone https://github.com/kopardev/circRNA.git
# ## and set that as the pipeline home
PIPELINE_HOME="/data/Ziegelbauer_lab/circRNADetection/scripts/circRNA"
echo "Pipeline Dir: $PIPELINE_HOME"
# ## make current folder as working directory
a=$(readlink -f $0)
# echo $a
for i in `seq 1 20`;do
b=$(echo $a|sed "s/\/gpfs\/gsfs${i}\/users/\/data/g")
a=$b
done
# echo $a
WORKDIR=$(dirname $a)
GIT_COMMIT_TAG=$(get_git_commitid_tag $PIPELINE_HOME)
echo "Git Commit/Tag: $GIT_COMMIT_TAG"

function usage() { cat << EOF
test.sh: test the workflow.
USAGE:
  bash test.sh <MODE> 
Required Positional Argument:
  MODE: [Type: Str] Valid options:
    a) init: initial workdir
    b) run: run test data
    c) cleanup: delete folders to get ready for re-init
    d) reset: cleanup followed by init
    e) dryrun: snakemake --dry-run
    f) unlock: snakemake --unlock
    g) runlocal: run without submitting to sbatch
EOF
}

function err() { cat <<< "$@" 1>&2; }

function init() {
echo "Working Dir: $WORKDIR"
cd $WORKDIR

# make sure that config folder exists in the current folder
if [ ! -d $WORKDIR/config ]; then cp -r $PIPELINE_HOME/config $WORKDIR/;echo "Config Dir: $WORKDIR/config";fi
if [ ! -d $WORKDIR/scripts ]; then cp -r $PIPELINE_HOME/scripts $WORKDIR/;echo "Scripts Dir: $WORKDIR/scripts";fi
if [ ! -d $WORKDIR/resources ]; then cp -r $PIPELINE_HOME/resources $WORKDIR/;echo "Resources Dir: $WORKDIR/resources";fi

#create log and stats folders
if [ ! -d $WORKDIR/logs ]; then mkdir -p $WORKDIR/logs;echo "Logs Dir: $WORKDIR/logs";fi
if [ ! -d $WORKDIR/stats ];then mkdir -p $WORKDIR/stats;echo "Stats Dir: $WORKDIR/stats";fi

}

function run () {
  ## initialize if not already done
  echo "Working Dir: $WORKDIR"
  if [ ! -d $WORKDIR/config ]; then err "Error: config folder not found ... initialize first!";usage && exit 1;fi
  for f in config.yaml tools.yaml cluster.json samples.tsv; do
    if [ ! -f $WORKDIR/config/$f ]; then err "Error: '${f}' file not found in config folder ... initialize first!";usage && exit 1;fi
  done
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

  # --use-singularity \
  # --singularity-args " -B ${WORKDIR}:${WORKDIR} -B /data/Ziegelbauer_lab/resources/:/data/Ziegelbauer_lab/resources/" \

  # --use-conda \
  # --use-envmodules \


  if [ "$1" == "local" ];then

  snakemake -s ${PIPELINE_HOME}/circRNADetection.snakefile \
  --directory $WORKDIR \
  --printshellcmds \
  --use-singularity \
  --singularity-args " -B ${WORKDIR}:${WORKDIR} -B /data/Ziegelbauer_lab/resources/:/data/Ziegelbauer_lab/resources/" \
  --use-envmodules \
  --latency-wait 120 \
  --configfile ${WORKDIR}/config/config.yaml \
  --cores all \
  --stats ${WORKDIR}/snakemake.stats \
  2>&1|tee ${WORKDIR}/snakemake.log

  else

  snakemake $1 -s ${PIPELINE_HOME}/circRNADetection.snakefile \
  --directory $WORKDIR \
  --use-singularity \
  --singularity-args " -B ${WORKDIR}:${WORKDIR} -B /data/Ziegelbauer_lab/resources/:/data/Ziegelbauer_lab/resources/" \
  --use-envmodules \
  --printshellcmds \
  --latency-wait 120 \
  --configfile ${WORKDIR}/config/config.yaml \
  --cluster-config ${WORKDIR}/config/cluster.json \
  --cluster "sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name {cluster.name} --output {cluster.output} --error {cluster.error}" \
  -j 500 \
  --rerun-incomplete \
  --keep-going \
  --stats ${WORKDIR}/snakemake.stats \
  2>&1|tee ${WORKDIR}/snakemake.log

  fi

  mv slurm*out stats && for a in $(ls stats/slurm*out);do gzip -n $a;done
}

function cleanup() {
  rm -rf ${WORKDIR}/config
  rm -rf ${WORKDIR}/logs
  rm -rf ${WORKDIR}/stats
  rm -rf ${WORKDIR}/scripts
  rm -rf ${WORKDIR}/resources
  rm -rf *snakemake*
}


function main(){

  if [ $# -eq 0 ]; then usage; exit 1; fi

  case $1 in
    init) init && exit 0;;
    dryrun) run --dry-run && exit 0;;
    unlock) run --unlock && exit 0;;
    run) run "" && exit 0;;
    runlocal) run local && exit 0;;
    cleanup) cleanup && exit 0;;
    reset) cleanup && init && exit 0;;
    -h | --help | help) usage && exit 0;;
    -* | --*) err "Error: Failed to provide mode: <init|run>."; usage && exit 1;;
    *) err "Error: Failed to provide mode: <init|run>. '${1}' is not supported."; usage && exit 1;;
  esac
}

main "$@"





