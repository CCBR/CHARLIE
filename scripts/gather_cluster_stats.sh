#!/bin/bash
## AUTHOR : Vishal N. Koparde, Ph.D., CCBR, NCI
## DATE : Feb 2021
## This scripts gathers cluster related statistics for jobs run on Biowulf using Snakemake by:
##   > extracting "external" jobids from snakemake.log files
##   > gather cluster stats for each job using "jobdata" and "jobhist" commands
##   > sorts output by job submission time
##   > output TSV file
##

function get_jobid_stats {
jobid=$1
declare -A jobdataarray
notaccountablejob=$(jobhist $jid|grep "No accounting"|wc -l)
if [ "$notaccountablejob" == "1" ];then
	jobdataarray["submit_time"]="JOBNOTACCOUNTABLE"
	jobdataarray["jobid"]="$jobid"
else
	jobdata $jobid > ${jobid}.tmp
	awk -F"\t" '{if (NF==2) {print}}' ${jobid}.tmp > ${jobid}.data && rm -f ${jobid}.tmp
	while read a b;do
		jobdataarray["$a"]="$b"
	done < ${jobid}.data
	rm -f ${jobid}.data
	st=${jobdataarray["submit_time"]}
	jobdataarray["human_submit_time"]=$(date -d @$st|sed "s/ /_/g")
fi
echo -ne "${jobdataarray["submit_time"]}\t"
echo -ne "${jobdataarray["human_submit_time"]}\t"
echo -ne "${jobdataarray["jobid"]}:${jobdataarray["state"]}:${jobdataarray["job_name"]}\t"
echo -ne "${jobdataarray["alloc_node"]}\t"
echo -ne "${jobdataarray["queued"]}:${jobdataarray["elapsed"]}:${jobdataarray["time_limit"]}\t"
echo -ne "${jobdataarray["avg_cpus"]}:${jobdataarray["max_cpu_used"]}:${jobdataarray["cpus_per_task"]}\t"
echo -ne "${jobdataarray["avg_mem"]}:${jobdataarray["max_mem_used"]}:${jobdataarray["total_mem"]}\t"
echo -ne "${jobdataarray["partition"]}:${jobdataarray["qos"]}\t"
echo -ne "${jobdataarray["username"]}:${jobdataarray["groupname"]}:${jobdataarray["account"]}\t"
echo -ne "${jobdataarray["work_dir"]}\t"
echo -ne "${jobdataarray["std_out"]}\t"
echo -ne "${jobdataarray["std_err"]}\n"
}

workdir=$(readlink -f $(dirname "$0"))
find $workdir -name "snakemake*log" -exec grep "with external jobid" {} \; |awk '{print $NF}' |sed "s/['.]//g" > jobids.lst
echo -ne "##SubmitTime\tHumanSubmitTime\tJobID:JobState:JobName\tNode\tQueueTime:RunTime:TimeLimit\tAvgCPU:MaxCPU:CPULimit\tAvgMEM:MaxMEM:MEMLimit\tPartition:QOS\tUsername:Group:Account\tWorkdir\tStdOut\tStdErr\n" > HPC_summary.txt
while read jid;do
get_jobid_stats $jid
# done < jobids.lst
done < jobids.lst |sort -k1,1n >> HPC_summary.txt
rm -f jobids.lst