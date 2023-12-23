#!/usr/bin/env bash
# Author: Vishal Koparde, Ph.D.
# CCBR, NCI
# (c) 2021
# CHARLIE

set -eo pipefail
module purge

# decide trigger
trigger="mtime"
# trigger="input"
# trigger="code"

##########################################################################################
# functions
##########################################################################################

function get_git_commitid_tag() {
  cd $1
  gid=$(git rev-parse HEAD)
  tag=$(git describe --tags $gid)
  echo -ne "$gid\t$tag"
}

##########################################################################################
# initial setup
##########################################################################################

# set PIPELINE_HOME
PIPELINE_HOME=$(readlink -f $(dirname "$0"))

# set snakefile
SNAKEFILE="${PIPELINE_HOME}/workflow/Snakefile"

# get github commit tag
GIT_COMMIT_TAG=$(get_git_commitid_tag $PIPELINE_HOME)

##########################################################################################
# Some more set up
##########################################################################################

CLUSTER_SBATCH_CMD="sbatch --parsable --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name {cluster.name} --output {cluster.output} --error {cluster.error}"
# if [ "$HOSTNAME" == "biowulf.nih.gov" ];then
# if [ "$SLURM_CLUSTER_NAME" == "biowulf" ];then
  EXTRA_SINGULARITY_BINDS="/lscratch"
  CLUSTER_SBATCH_CMD="$CLUSTER_SBATCH_CMD --gres {cluster.gres}"
# fi
PYTHONVERSION="3.7"
SNAKEMAKEVERSION="7.19.1"
# SNAKEMAKEVERSION="5.24.1"

# set defaults
HOST="hg38"
ADDITIVES="ERCC"
VIRUSES="NC_009333.1"
MANIFEST="${PIPELINE_HOME}/config/samples.tsv"

# set variables
SCRIPTNAME="$0"
SCRIPTDIRNAME=$(readlink -f $(dirname $0))
SCRIPTBASENAME=$(basename $0)

##########################################################################################
# USAGE
##########################################################################################

function usage() { cat << EOF

##########################################################################################

Welcome to charlie(v0.10.1)
 _______  __   __  _______  ______    ___      ___   _______ 
|       ||  | |  ||   _   ||    _ |  |   |    |   | |       |
|       ||  |_|  ||  |_|  ||   | ||  |   |    |   | |    ___|
|       ||       ||       ||   |_||_ |   |    |   | |   |___ 
|      _||       ||       ||    __  ||   |___ |   | |    ___|
|     |_ |   _   ||   _   ||   |  | ||       ||   | |   |___ 
|_______||__| |__||__| |__||___|  |_||_______||___| |_______|

C_ircrnas in H_ost A_nd vi_R_uses ana_L_ysis p_I_p_E_line

##########################################################################################

This pipeline was built by CCBR (https://bioinformatics.ccr.cancer.gov/ccbr)
Please contact Vishal Koparde for comments/questions (vishal.koparde@nih.gov)

##########################################################################################

CHARLIE can be used to DAQ(Detect/Annotate/Quantify) circRNAs in hosts and viruses.

Here is the list of hosts and viruses that are currently supported:

HOSTS:
  * hg38          [Human]
  * mm39          [Mouse]

ADDITIVES:
  * ERCC          [External RNA Control Consortium sequences]
  * BAC16Insert   [insert from rKSHV.219-derived BAC clone of the full-length KSHV genome]

VIRUSES:
  * NC_007605.1   [Human gammaherpesvirus 4 (Epstein-Barr virus)]
  * NC_006273.2   [Human betaherpesvirus 5 (Cytomegalovirus )]
  * NC_001664.4   [Human betaherpesvirus 6A (HHV-6A)]
  * NC_000898.1   [Human betaherpesvirus 6B (HHV-6B)]
  * NC_001716.2   [Human betaherpesvirus 7 (HHV-7)]
  * NC_009333.1   [Human gammaherpesvirus 8 (KSHV)]
  * NC_045512.2   [Severe acute respiratory syndrome(SARS)-related coronavirus]
  * MN485971.1    [HIV from Belgium]
  * NC_001806.2   [Human alphaherpesvirus 1 (Herpes simplex virus type 1)](strain 17) (HSV-1)]
  * KT899744.1    [HSV-1 strain KOS]
  * MH636806.1    [MHV68 (Murine herpesvirus 68 strain WUMS)]

##########################################################################################

USAGE:
  bash ${SCRIPTNAME} -w/--workdir=<WORKDIR> -m/--runmode=<RUNMODE>

Required Arguments:
1.  WORKDIR     : [Type: String]: Absolute or relative path to the output folder with write permissions.

2.  RUNMODE     : [Type: String] Valid options:
    * init      : initialize workdir
    * dryrun    : dry run snakemake to generate DAG
    * run       : run with slurm
    * runlocal  : run without submitting to sbatch
    ADVANCED RUNMODES (use with caution!!)
    * unlock    : unlock WORKDIR if locked by snakemake NEVER UNLOCK WORKDIR WHERE PIPELINE IS CURRENTLY RUNNING!
    * reconfig  : recreate config file in WORKDIR (debugging option) EDITS TO config.yaml WILL BE LOST!
    * reset     : DELETE workdir dir and re-init it (debugging option) EDITS TO ALL FILES IN WORKDIR WILL BE LOST!
    * printbinds: print singularity binds (paths)
    * local     : same as runlocal

Optional Arguments:

--host|-g       : supply host at command line. hg38 or mm39.                                            (--runmode=init only)
--additives|-a  : supply comma-separated list of additives at command line. ERCC or BAC16Insert or both (--runmode=init only)
--viruses|-v    : supply comma-separated list of viruses at command line                                (--runmode=init only)
--manifest|-s   : absolute path to samples.tsv. This will be copied to output folder                    (--runmode=init only)
--changegrp|-z  : change group to "Ziegelbauer_lab" before running anything. Biowulf-only. Useful for correctly setting permissions.
--help|-h       : print this help


Example commands:
  bash ${SCRIPTNAME} -w=/my/ouput/folder -m=init
  bash ${SCRIPTNAME} -w=/my/ouput/folder -m=dryrun
  bash ${SCRIPTNAME} -w=/my/ouput/folder -m=run

##########################################################################################

VersionInfo: 
  python          : $PYTHONVERSION 
  snakemake       : $SNAKEMAKEVERSION 
  pipeline_home   : $PIPELINE_HOME
  git commit/tag  : $GIT_COMMIT_TAG

##########################################################################################
EOF
}

##########################################################################################
# ERR
##########################################################################################

function err() { usage && cat <<< "
#
# ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR
#
  $@
#
# ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR
#
" && exit 1 1>&2; }

##########################################################################################
# INIT
##########################################################################################

function init() {

# create output folder
if [ -d $WORKDIR ];then err "Folder $WORKDIR already exists!"; fi
mkdir -p $WORKDIR

# copy config and samples files
if [ ! -f $CONFIGFILE ];then
sed -e "s/PIPELINE_HOME/${PIPELINE_HOME//\//\\/}/g" \
    -e "s/WORKDIR/${WORKDIR//\//\\/}/g" \
    -e "s/HOST/${HOST}/g" \
    -e "s/ADDITIVES/${ADDITIVES}/g" \
    -e "s/VIRUSES/${VIRUSES}/g" \
    ${PIPELINE_HOME}/config/config.yaml > $CONFIGFILE
fi
if [ ! -f $WORKDIR/nclscan.config ];then
sed -e "s/PIPELINE_HOME/${PIPELINE_HOME//\//\\/}/g" -e "s/WORKDIR/${WORKDIR//\//\\/}/g" ${PIPELINE_HOME}/resources/NCLscan.config.template > $WORKDIR/nclscan.config
fi
if [ ! -f $CLUSTERFILE ];then
cp ${PIPELINE_HOME}/resources/cluster.json $CLUSTERFILE
fi
if [ ! -f $WORKDIR/samples.tsv ];then
cp $MANIFEST $WORKDIR/samples.tsv
fi
if [ ! -f $JOBBY ];then
cp ${PIPELINE_HOME}/workflow/scripts/jobby $JOBBY && chmod a+x $JOBBY
fi
if [ ! -f $JOBBY2 ];then
cp ${PIPELINE_HOME}/workflow/scripts/run_jobby_on_snakemake_log $WORKDIR && chmod a+x $JOBBY2
fi

#create log and stats folders
if [ ! -d $WORKDIR/logs ]; then mkdir -p $WORKDIR/logs;echo "Logs Dir: $WORKDIR/logs";fi
if [ ! -d $WORKDIR/stats ];then mkdir -p $WORKDIR/stats;echo "Stats Dir: $WORKDIR/stats";fi

echo "Done Initializing $WORKDIR. You can now edit $WORKDIR/config.yaml and $WORKDIR/samples.tsv"

}

##########################################################################################
# CHECK ESSENTIAL FILES
##########################################################################################

function check_essential_files() {
  if [ ! -d $WORKDIR ];then err "Folder $WORKDIR does not exist!"; fi
  for f in config.yaml samples.tsv nclscan.config jobby cluster.json; do
    if [ ! -f $WORKDIR/$f ]; then err "Error: '${f}' file not found in workdir ... initialize first!";fi
  done
}

##########################################################################################
# CHANGE GROUP to Ziegelbauer_lab
##########################################################################################

function change_grp() {
  grps=$(groups)
  found=0
  current_grp=""
  count=0
  for g in $grps
  do
    count=$((count+1))
    if [ "$count" == "1" ];then
      current_grp=$g
    fi
    if [ "Ziegelbauer_lab" == $g ]
    then
      found=1
      break
    fi
  done
  if [ "$found" == "1" ]
  then
    if [ "$current_grp" != "Ziegelbauer_lab" ]
    then
      echo "Current group     : $current_grp"
      echo "Changing to group : Ziegelbauer_lab"
      # newgrp Ziegelbauer_lab
      exec sg Ziegelbauer_lab "$0 $*"
    fi
  fi
  # main "$@"
}

##########################################################################################
# RECONFIG ... recreate config.yaml and overwrite old version
##########################################################################################

function reconfig(){
  # rebuild config file and replace the config.yaml in the WORKDIR
  # this is only for dev purposes when new key-value pairs are being added to the config file
  check_essential_files
  sed -e "s/PIPELINE_HOME/${PIPELINE_HOME//\//\\/}/g" -e "s/WORKDIR/${WORKDIR//\//\\/}/g" ${PIPELINE_HOME}/config/config.yaml > $WORKDIR/config.yaml
  echo "$WORKDIR/config.yaml has been updated!"
}

##########################################################################################
# RUNCHECK ... check essential files and load required packages
##########################################################################################

function runcheck(){
  check_essential_files
  module load python/$PYTHONVERSION
  module load snakemake/$SNAKEMAKEVERSION
}

##########################################################################################
# DRYRUN ... also run automatically before actual run
##########################################################################################

function dryrun() {
  runcheck
  timestamp=$(date +"%y%m%d%H%M%S")
  nfiles=$(find ${WORKDIR} -maxdepth 1 -name "dryrun.*.log"|wc -l)
  if [ "$nfiles" != "0" ];then
    for f in $(ls ${WORKDIR}/dryrun.*.log);do
      mv $f ${WORKDIR}/stats/
    done
  fi
  run "--dry-run" | tee ${WORKDIR}/dryrun.${timestamp}.log
}

function touch() {
  runcheck
  timestamp=$(date +"%y%m%d%H%M%S")
  run "--touch" | tee ${WORKDIR}/touch.${timestamp}.log
}

##########################################################################################
# UNLOCK
##########################################################################################

function unlock() {
  runcheck
  run "--unlock"  
}

##########################################################################################
# SET SINGULARITY BINDS ... bind required singularity folders appropriately
##########################################################################################

function set_singularity_binds(){
# this functions tries find what folders to bind
# biowulf specific
  echo "$PIPELINE_HOME" > ${WORKDIR}/tmp1
  echo "$WORKDIR" >> ${WORKDIR}/tmp1
  grep -o '\/.*' <(cat ${WORKDIR}/config.yaml ${WORKDIR}/samples.tsv)|dos2unix|tr '\t' '\n'|grep -v ' \|\/\/'|sort|uniq >> ${WORKDIR}/tmp1
  grep gpfs ${WORKDIR}/tmp1|awk -F'/' -v OFS='/' '{print $1,$2,$3,$4,$5}'| grep "[a-zA-Z0-9]" |sort|uniq > ${WORKDIR}/tmp2
  grep -v gpfs ${WORKDIR}/tmp1|awk -F'/' -v OFS='/' '{print $1,$2,$3}'| grep "[a-zA-Z0-9]"|sort|uniq > ${WORKDIR}/tmp3
  while read a;do readlink -f $a;done < ${WORKDIR}/tmp3 | grep "[a-zA-Z0-9]"> ${WORKDIR}/tmp4
  binds=$(cat ${WORKDIR}/tmp2 ${WORKDIR}/tmp3 ${WORKDIR}/tmp4|sort|uniq |tr '\n' ',')
  rm -f ${WORKDIR}/tmp?
  binds=$(echo $binds|awk '{print substr($1,1,length($1)-1)}')
  SINGULARITY_BINDS="-B $EXTRA_SINGULARITY_BINDS,$binds"
}

##########################################################################################
# PRINT SINGULARITY BINDS ... print bound singularity folders for debugging
##########################################################################################

function printbinds(){
  set_singularity_binds
  echo $SINGULARITY_BINDS
}

##########################################################################################
# RUNLOCAL ... run directly on local interactive node ... no submission to SLURM
##########################################################################################

function runlocal() {
  runcheck
  set_singularity_binds
  if [ "$SLURM_JOB_ID" == "" ];then err "runlocal can only be done on an interactive node"; exit 1; fi
  module load singularity
  run "local"
}

##########################################################################################
# RUNSLURM ... submit head job to slurm which will spawn other jobs on SLURM
##########################################################################################

function runslurm() {
  runcheck
  set_singularity_binds
  run "--dry-run " && run "slurm"
}

##########################################################################################
# CREATE RUNINFO ... create runinfo.yaml in workdir
##########################################################################################

function create_runinfo {
  modtime=$1
  if [ "$modtime" == "" ];then
   modtime=$(stat ${WORKDIR}/runinfo.yaml|grep Modify|awk '{print $2,$3}'|awk -F"." '{print $1}'|sed "s/ //g"|sed "s/-//g"|sed "s/://g")
  fi
  if [ -f ${WORKDIR}/runinfo.yaml ];then
    mv ${WORKDIR}/runinfo.yaml ${WORKDIR}/stats/runinfo.${modtime}.yaml
  fi
  echo "Pipeline Dir: $PIPELINE_HOME" > ${WORKDIR}/runinfo.yaml
  echo "Git Commit/Tag: $GIT_COMMIT_TAG" >> ${WORKDIR}/runinfo.yaml
  userlogin=$(whoami)
  if [[ `which finger 2>/dev/null` ]];then 
	  username=$(finger $userlogin |grep ^Login | awk -F"Name: " '{print $2}'); 
  elif [[ `which lslogins 2>/dev/null` ]];then 
	  username=$(lslogins -u $userlogin | grep ^Geco | awk -F": " '{print $2}' | awk '{$1=$1;print}'); 
  else username="";fi
  echo "Login: $userlogin" >> ${WORKDIR}/runinfo.yaml
  echo "Name: $username" >> ${WORKDIR}/runinfo.yaml
  g=$(groups)
  echo "Groups: $g" >> ${WORKDIR}/runinfo.yaml
  d=$(date)
  echo "Date/Time: $d" >> ${WORKDIR}/runinfo.yaml
}

##########################################################################################
# PRERUN CLEANUP ... get ready to run .. park old logs/stats etc.
##########################################################################################

function preruncleanup() {
  echo "Running..."

  # check initialization
  check_essential_files 

  cd $WORKDIR
  modtime=""
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
    if [ -f ${WORKDIR}/snakemake.log.jobinfo ];then 
      mv ${WORKDIR}/snakemake.log.jobinfo ${WORKDIR}/stats/snakemake.${modtime}.log.jobinfo
    fi
  fi
  nslurmouts=$(find ${WORKDIR} -maxdepth 1 -name "slurm-*.out" |wc -l)
  if [ "$nslurmouts" != "0" ];then
    for f in $(ls ${WORKDIR}/slurm-*.out);do mv ${f} ${WORKDIR}/logs/;done
  fi

  create_runinfo modtime

}

##########################################################################################
# RUN wrapper for all possible run's
# a. dryrun
# b. unlock
# c. local run
# d. slurm run, etc.
##########################################################################################

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
  --configfile $CONFIGFILE \
  --cores all \
  --rerun-incomplete \
  --rerun-triggers $trigger \
  --retries 2 \
  --keep-going \
  --stats ${WORKDIR}/snakemake.stats \
  2>&1|tee ${WORKDIR}/snakemake.log

  if [ "$?" -eq "0" ];then
    snakemake -s $SNAKEFILE \
    --report ${WORKDIR}/runlocal_snakemake_report.html \
    --directory $WORKDIR \
    --configfile $CONFIGFILE
  fi

  elif [ "$1" == "slurm" ];then
  
  preruncleanup

  cat > ${WORKDIR}/submit_script.sbatch << EOF
#!/bin/bash
#SBATCH --job-name="charlie"
#SBATCH --mem=40g
#SBATCH --partition="ccr,norm"
#SBATCH --time=48:00:00
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
--configfile $CONFIGFILE \
--cluster-config $CLUSTERFILE \
--cluster "$CLUSTER_SBATCH_CMD" \
--cluster-status $CLUSTERSTATUSCMD \
-j 500 \
--rerun-incomplete \
--rerun-triggers $trigger \
--retries 2 \
--keep-going \
--stats ${WORKDIR}/snakemake.stats \
2>&1|tee ${WORKDIR}/snakemake.log

if [ "\$?" -eq "0" ];then
  snakemake -s $SNAKEFILE \
  --directory $WORKDIR \
  --report ${WORKDIR}/runslurm_snakemake_report.html \
  --configfile $CONFIGFILE 
fi

EOF

  sbatch ${WORKDIR}/submit_script.sbatch

  elif [ "$1" == "--touch" ];then

snakemake $1 -s $SNAKEFILE \
--directory $WORKDIR \
--configfile $CONFIGFILE \
--cores 1

  else # dry-run and unlock

echo $CLUSTER_SBATCH_CMD

snakemake $1 -s $SNAKEFILE \
--directory $WORKDIR \
--use-envmodules \
--printshellcmds \
--latency-wait 120 \
--configfile $CONFIGFILE \
--cluster-config $CLUSTERFILE \
--cluster "$CLUSTER_SBATCH_CMD" \
-j 500 \
--rerun-incomplete \
--rerun-triggers $trigger \
--keep-going \
--reason \
--stats ${WORKDIR}/snakemake.stats

  fi

}

##########################################################################################
# RESET ... delete workdir and then initialize
##########################################################################################

function reset() {
  #delete the workdir and re-initialize it
  echo "Working Dir: $WORKDIR"
  if [ ! -d $WORKDIR ];then err "Folder $WORKDIR does not exist!";fi
  echo "Deleting $WORKDIR"
  rm -rf $WORKDIR
  echo "Re-Initializing $WORKDIR"
  init
}

##########################################################################################
# MAIN ... command line argument parsing
##########################################################################################

function main(){

  CHANGEGRP=0

  if [ $# -eq 0 ]; then usage; exit 1; fi

  allargs="$@"

  for i in "$@"
  do
  case $i in
      -m=*|--runmode=*)
        RUNMODE="${i#*=}"
      ;;
      -w=*|--workdir=*)
        WORKDIR="${i#*=}"
      ;;
      -z|--changegrp)
        CHANGEGRP=1
      ;;
      -g=*|--host=*)
        HOST="${i#*=}"
      ;;
      -a=*|--additives=*)
        ADDITIVES="${i#*=}"
      ;;
      -v=*|--viruses=*)
        VIRUSES="${i#*=}"
      ;;      
      -s=*|--manifest=*)
        MANIFEST="${i#*=}"
        if [ ! -f $MANIFEST ];then err "File $MANIFEST does NOT exist!";fi
      ;;      
      -h|--help)
        usage && exit 0;
      ;;
      *)
        err "Unknown argument $i!"    # unknown option
      ;;
  esac
  done
  WORKDIR=$(readlink -f "$WORKDIR")
  echo "Working Dir: $WORKDIR"

  # required files
  CONFIGFILE="${WORKDIR}/config.yaml"
  CLUSTERFILE="${WORKDIR}/cluster.json"
  JOBBY="${WORKDIR}/jobby"
  JOBBY2="${WORKDIR}/run_jobby_on_snakemake_log"

  CLUSTERSTATUSCMD="${PIPELINE_HOME}/resources/cluster_status.sh"

  # change group to Ziegelbauer_lab before doing anything
  if [ "$CHANGEGRP" == "1" ]; then change_grp "$allargs"; fi

  case $RUNMODE in
    init) init && exit 0;;
    dryrun) dryrun && exit 0;;
    unlock) unlock && exit 0;;
    run) runslurm && exit 0;;
    runlocal) runlocal && exit 0;;
    reset) reset && exit 0;;
    touch) touch && exit 0;;
    dry) dryrun && exit 0;;                      # hidden option
    local) runlocal && exit 0;;                  # hidden option
    reconfig) reconfig && exit 0;;               # hidden option for debugging
    printbinds) printbinds && exit 0;;           # hidden option
    help) usage && exit 0;;                      # print help
    *) err "Unknown RUNMODE \"$RUNMODE\"";;
  esac
}

##########################################################################################
# run main!
##########################################################################################

main "$@"





