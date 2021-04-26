#!/usr/bin/env bash
# Author: Vishal Koparde, Ph.D.
# CCBR, NCI
# (c) 2021
#
#
# cirna_daq
# d=detection
# a=annotation
# q=quantification
set -eo pipefail
module purge

# you may want to re-set these
#######
EXTRA_SINGULARITY_BINDS="/lscratch"
PYTHONVERSION="3.7"
SNAKEMAKEVERSION="5.24.1"
#######


SCRIPTNAME="$0"
SCRIPTDIRNAME=$(readlink -f $(dirname $0))
SCRIPTBASENAME=$(readlink -f $(basename $0))


function get_git_commitid_tag() {
  cd $1
  gid=$(git rev-parse HEAD)
  tag=$(git describe --tags $gid)
  echo -ne "$gid\t$tag"
}

# ## setting PIPELINE_HOME
PIPELINE_HOME=$(readlink -f $(dirname "$0"))
echo "Pipeline Dir: $PIPELINE_HOME"
# set snakefile
SNAKEFILE="${PIPELINE_HOME}/workflow/Snakefile"

# get github commit tag
GIT_COMMIT_TAG=$(get_git_commitid_tag $PIPELINE_HOME)
echo "Git Commit/Tag: $GIT_COMMIT_TAG"

function usage() { cat << EOF

${SCRIPTBASENAME}
--> run circRNA Detection Annotation Quantification Pipeline

USAGE:
  bash ${SCRIPTNAME} -m/--runmode=<RUNMODE> -w/--workdir=<WORKDIR>
Required Arguments:
1.  RUNMODE: [Type: String] Valid options:
    *) init : initialize workdir
    *) run : run with slurm
    *) reset : DELETE workdir dir and re-init it
    *) dryrun : dry run snakemake to generate DAG
    *) unlock : unlock workdir if locked by snakemake
    *) runlocal : run without submitting to sbatch
2.  WORKDIR: [Type: String]: Absolute or relative path to the output folder with write permissions.
EOF
}

function err() { cat <<< "
#
#
#
  $@
#
#
#
" && usage && exit 1 1>&2; }

# function check_arguments() {
#   if [ "$#" -eq "1" ]; then err "init needs an absolute path to the working dir"; usage; exit 1; fi
#   if [ "$#" -gt "2" ]; then err "init takes only one more argument"; usage; exit 1;fi
#   WORKDIR=$2
# # echo $WORKDIR
# # x=$(echo $WORKDIR|awk '{print substr($1,1,1)}')
# # if [ "$x" != "/" ]; then err "working dir should be supplied as an absolute path"; usage; exit 1; fi
#   WORKDIR=$(readlink -f "$WORKDIR")
#   echo "Working Dir: $WORKDIR"
#   # exit 1
# }

function init() {

# create output folder
if [ -d $WORKDIR ];then err "Folder $WORKDIR already exists!"; fi
mkdir -p $WORKDIR

# copy config and samples files
sed -e "s/PIPELINE_HOME/${PIPELINE_HOME//\//\\/}/g" -e "s/WORKDIR/${WORKDIR//\//\\/}/g" ${PIPELINE_HOME}/config/config.yaml > $WORKDIR/config.yaml
cp ${PIPELINE_HOME}/config/samples.tsv $WORKDIR/

#create log and stats folders
if [ ! -d $WORKDIR/logs ]; then mkdir -p $WORKDIR/logs;echo "Logs Dir: $WORKDIR/logs";fi
if [ ! -d $WORKDIR/stats ];then mkdir -p $WORKDIR/stats;echo "Stats Dir: $WORKDIR/stats";fi

echo "Done Initializing $WORKDIR. You can now edit $WORKDIR/config.yaml and $WORKDIR/samples.tsv"

}

function check_essential_files() {
  if [ ! -d $WORKDIR ];then err "Folder $WORKDIR does not exist!"; fi
  for f in config.yaml samples.tsv; do
    if [ ! -f $WORKDIR/$f ]; then err "Error: '${f}' file not found in workdir ... initialize first!";fi
  done
}

function reconfig(){
  # rebuild config file and replace the config.yaml in the WORKDIR
  # this is only for dev purposes when new key-value pairs are being added to the config file
  check_essential_files
  sed -e "s/PIPELINE_HOME/${PIPELINE_HOME//\//\\/}/g" -e "s/WORKDIR/${WORKDIR//\//\\/}/g" ${PIPELINE_HOME}/config/config.yaml > $WORKDIR/config.yaml
  echo "$WORKDIR/config.yaml has been updated!"
}

function runcheck(){
  check_essential_files
  module load python/$PYTHONVERSION
  module load snakemake/$SNAKEMAKEVERSION
}

function dryrun() {
  runcheck
  run "--dry-run"
}

function unlock() {
  runcheck
  run "--unlock"  
}

function set_singularity_binds(){
# this functions tries find what folders to bind
# biowulf specific
  echo "$PIPELINE_HOME" > ${WORKDIR}/tmp1
  echo "$WORKDIR" >> ${WORKDIR}/tmp1
  grep -o '\/.*' <(cat ${WORKDIR}/config.yaml ${WORKDIR}/samples.tsv)|tr '\t' '\n'|grep -v ' \|\/\/'|sort|uniq >> ${WORKDIR}/tmp1
  grep gpfs ${WORKDIR}/tmp1|awk -F'/' -v OFS='/' '{print $1,$2,$3,$4,$5}' |sort|uniq > ${WORKDIR}/tmp2
  grep -v gpfs ${WORKDIR}/tmp1|awk -F'/' -v OFS='/' '{print $1,$2,$3}'|sort|uniq > ${WORKDIR}/tmp3
  while read a;do readlink -f $a;done < ${WORKDIR}/tmp3 > ${WORKDIR}/tmp4
  binds=$(cat ${WORKDIR}/tmp2 ${WORKDIR}/tmp3 ${WORKDIR}/tmp4|sort|uniq |tr '\n' ',')
  rm -f ${WORKDIR}/tmp?
  binds=$(echo $binds|awk '{print substr($1,1,length($1)-1)}')
  SINGULARITY_BINDS="-B $EXTRA_SINGULARITY_BINDS,$binds"
}

function printbinds(){
  set_singularity_binds
  echo $SINGULARITY_BINDS
}

function runlocal() {
  runcheck
  set_singularity_binds
  if [ "$SLURM_JOB_ID" == "" ];then err "runlocal can only be done on an interactive node"; exit 1; fi
  module load singularity
  run "local"
}

function runslurm() {
  runcheck
  set_singularity_binds
  run "slurm"
}

function create_runinfo {
  if [ -f ${WORKDIR}/runinfo.yaml ];then
    modtime=$(stat ${WORKDIR}/runinfo.yaml|grep Modify|awk '{print $2,$3}'|awk -F"." '{print $1}'|sed "s/ //g"|sed "s/-//g"|sed "s/://g")
    mv ${WORKDIR}/runinfo.yaml ${WORKDIR}/runinfo.yaml.${modtime}
  fi
  echo "Pipeline Dir: $PIPELINE_HOME" > ${WORKDIR}/runinfo.yaml
  echo "Git Commit/Tag: $GIT_COMMIT_TAG" >> ${WORKDIR}/runinfo.yaml
  userlogin=$(whoami)
  username=$(finger $userlogin|grep ^Login|awk -F"Name: " '{print $2}')
  echo "Login: $userlogin" >> ${WORKDIR}/runinfo.yaml
  echo "Name: $username" >> ${WORKDIR}/runinfo.yaml
  g=$(groups)
  echo "Groups: $g" >> ${WORKDIR}/runinfo.yaml
  d=$(date)
  echo "Date/Time: $d" >> ${WORKDIR}/runinfo.yaml
}

function preruncleanup() {
  echo "Running..."

  # check initialization
  check_essential_files 

  cd $WORKDIR
  ## Archive previous run files
  if [ -f ${WORKDIR}/snakemake.log ];then 
    modtime=$(stat ${WORKDIR}/snakemake.log |grep Modify|awk '{print $2,$3}'|awk -F"." '{print $1}'|sed "s/ //g"|sed "s/-//g"|sed "s/://g")
    mv ${WORKDIR}/snakemake.log ${WORKDIR}/stats/snakemake.${modtime}.log
    if [ -f ${WORKDIR}/snakemake.log.HPC_summary.txt ];then 
      mv ${WORKDIR}/snakemake.log.HPC_summary.txt ${WORKDIR}/stats/snakemake.${modtime}.log.HPC_summary.txt
    fi
    if [ -f ${WORKDIR}/snakemake.stats ];then 
      mv ${WORKDIR}/snakemake.stats ${WORKDIR}/stats/snakemake.${modtime}.stats
    fi
  fi
  nslurmouts=$(find ${WORKDIR} -maxdepth 1 -name "slurm-*.out" |wc -l)
  if [ "$nslurmouts" != "0" ];then
    for f in $(ls ${WORKDIR}/slurm-*.out);do mv ${f} ${WORKDIR}/logs/;done
  fi

  create_runinfo

}

function run() {


  if [ "$1" == "local" ];then

  preruncleanup

  snakemake -s $SNAKEFILE\
  --directory $WORKDIR \
  --printshellcmds \
  --use-singularity \
  --singularity-args "$SINGULARITY_BINDS" \
  --use-envmodules \
  --latency-wait 120 \
  --configfile ${WORKDIR}/config.yaml \
  --cores all \
  --stats ${WORKDIR}/snakemake.stats \
  2>&1|tee ${WORKDIR}/snakemake.log

  if [ "$?" -eq "0" ];then
    snakemake -s ${PIPELINE_HOME}/circRNADetection.snakefile \
    --report ${WORKDIR}/runlocal_snakemake_report.html \
    --directory $WORKDIR \
    --configfile ${WORKDIR}/config.yaml
  fi

  elif [ "$1" == "slurm" ];then
  
  preruncleanup

  cat > ${WORKDIR}/submit_script.sbatch << EOF
#!/bin/bash
#SBATCH --job-name="circRNA_DAQ"
#SBATCH --mem=10g
#SBATCH --partition="ccr,norm"
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=2

module load python/$PYTHONVERSION
module load snakemake/$SNAKEMAKEVERSION
module load singularity

cd \$SLURM_SUBMIT_DIR

snakemake -s $SNAKEFILE \
--directory $WORKDIR \
--use-singularity \
--singularity-args "$SINGULARITY_BINDS" \
--use-envmodules \
--printshellcmds \
--latency-wait 120 \
--configfile ${WORKDIR}/config.yaml \
--cluster-config ${PIPELINE_HOME}/resources/cluster.json \
--cluster "sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name {cluster.name} --output {cluster.output} --error {cluster.error} --qos {cluster.qos}" \
-j 500 \
--rerun-incomplete \
--keep-going \
--stats ${WORKDIR}/snakemake.stats \
2>&1|tee ${WORKDIR}/snakemake.log

if [ "\$?" -eq "0" ];then
  snakemake -s ${PIPELINE_HOME}/circRNADetection.snakefile \
  --directory $WORKDIR \
  --report ${WORKDIR}/runslurm_snakemake_report.html \
  --configfile ${WORKDIR}/config.yaml 
fi

bash <(curl https://raw.githubusercontent.com/CCBR/Tools/master/Biowulf/gather_cluster_stats_biowulf.sh 2>/dev/null) ${WORKDIR}/snakemake.log > ${WORKDIR}/snakemake.log.HPC_summary.txt

EOF

  sbatch ${WORKDIR}/submit_script.sbatch

  else # dry-run and unlock

snakemake $1 -s $SNAKEFILE \
--directory $WORKDIR \
--use-envmodules \
--printshellcmds \
--latency-wait 120 \
--configfile ${WORKDIR}/config.yaml \
--cluster-config ${PIPELINE_HOME}/resources/cluster.json \
--cluster "sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name {cluster.name} --output {cluster.output} --error {cluster.error}" \
-j 500 \
--rerun-incomplete \
--keep-going \
--reason \
--stats ${WORKDIR}/snakemake.stats

  fi

}

function reset() {
  #delete the workdir and re-initialize it
  echo "Working Dir: $WORKDIR"
  if [ ! -d $WORKDIR ];then err "Folder $WORKDIR does not exist!";fi
  echo "Deleting $WORKDIR"
  rm -rf $WORKDIR
  echo "Re-Initializing $WORKDIR"
  init
}


function main(){

  if [ $# -eq 0 ]; then usage; exit 1; fi

  for i in "$@"
  do
  case $i in
      -m=*|--runmode=*)
        RUNMODE="${i#*=}"
      ;;
      -w=*|--workdir=*)
        WORKDIR="${i#*=}"
      ;;
      *)
        err "Unknown argument $i!"    # unknown option
      ;;
  esac
  done
  WORKDIR=$(readlink -f "$WORKDIR")
  echo "Working Dir: $WORKDIR"
  # echo SCRIPTNAME = ${SCRIPTNAME}
  # echo RUNMODE = ${RUNMODE}
  # echo WORKDIR = ${WORKDIR}
  # exit;


  case $RUNMODE in
    init) init && exit 0;;
    dryrun) dryrun && exit 0;;
    unlock) unlock && exit 0;;
    run) runslurm && exit 0;;
    runlocal) runlocal && exit 0;;
    reset) reset && exit 0;;
    dry) dryrun && exit 0;;                      # hidden option
    local) runlocal && exit 0;;                  # hidden option
    reconfig) reconfig && exit 0;;               # hidden option for debugging
    printbinds) printbinds && exit 0;;           # hidden option
    *) err "Unknown RUNMODE \"$RUNMODE\"";;
  esac
}


main "$@"





