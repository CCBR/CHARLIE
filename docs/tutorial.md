# circRNA DAQ Pipeline

![img](https://img.shields.io/github/issues/kopardev/circRNA?style=for-the-badge)![img](https://img.shields.io/github/forks/kopardev/circRNA?style=for-the-badge)![img](https://img.shields.io/github/stars/kopardev/circRNA?style=for-the-badge)![img](https://img.shields.io/github/license/kopardev/circRNA?style=for-the-badge)

### Tutorial

#### Prerequisites

* [Biowulf](https://hpc.nih.gov/) account: Biowulf account can be requested [here](https://hpc.nih.gov/docs/accounts.html)

* Membership to Ziegelbauer user group on Biowulf. You can check this by typing the following command:

  ```bash
  % groups
  ```
  

output:

  ```bash
  CCBR kopardevn Ziegelbauer_lab
  ```

  If `Ziegelbauer_lab` is not listed then you can email a request to be added to the groups [here](mailto:staff@hpc.nih.gov)

#### Location

Different versions of circRNA DAQ pipeline have been parked at `/data/Ziegelbauer_lab/Pipelines/circRNA`

```bash
% ls -la /data/Ziegelbauer_lab/Pipelines/circRNA
```

output:

```bash
total 135
drwxrws--T. 6 kopardevn Ziegelbauer_lab 4096 Feb  4 18:27 .
drwxrws--T. 3 kopardevn Ziegelbauer_lab 4096 Jan 12 09:05 ..
lrwxrwxrwx. 1 kopardevn Ziegelbauer_lab   54 Feb  3 14:29 dev -> /data/Ziegelbauer_lab/circRNADetection/scripts/circRNA
drwxrws---. 3 kopardevn Ziegelbauer_lab 4096 Jan 12 09:07 v0.1.0
drwxrws---. 6 kopardevn Ziegelbauer_lab 4096 Jan 12 17:32 v0.2.1
drwxrws---. 6 kopardevn Ziegelbauer_lab 4096 Jan 14 10:24 v0.3.3
drwxrws---. 7 kopardevn Ziegelbauer_lab 4096 Feb  4 18:26 v0.4.0
```

The exacts versions listed here may changed as newer versions are added. Also, the `dev` version is pointing to the most recent untagged version of the pipeline (use at own risk!)

#### Init

To get help about the pipeline you can run:

```bash
% bash /data/Ziegelbauer_lab/Pipelines/circRNA/v0.4.0/run_circrna_daq.sh --help
```

output:

```bash
Pipeline Dir: /data/Ziegelbauer_lab/circRNADetection/scripts/circRNA
Git Commit/Tag: b2f387c1b6854646d12974cd16da1168d93bb43b	v0.4.0-14-gb2f387c
run_circrna_daq.sh: run the workflow to DAQ (detect, annotate and quantify circRNAs)
USAGE:
  bash run_circrna_daq.sh <MODE>
Required Positional Argument:
  MODE: [Type: Str] Valid options:
    a) init <path_to_workdir> : initialize workdir
    b) run <path_to_workdir>: run with slurm
    c) reset <path_to_workdir> : DELETE workdir dir and re-init it
    e) dryrun <path_to_workdir> : dry run snakemake to generate DAG
    f) unlock <path_to_workdir> : unlock workdir if locked by snakemake
    g) runlocal <path_to_workdir>: run without submitting to sbatch
```

You can replace `v0.4.0` in the above command with the latest version to use a newer version. `run_circrna_daq.sh` was called `test.sh` in versions older than `v0.4.0`.

To initial the working directory run:

```bash
% bash /data/Ziegelbauer_lab/Pipelines/circRNA/v0.4.0/run_circrna_daq.sh init /scratch/circRNA_daq_test
```

> **NOTE**
>
> With version 0.6.0 or newer the arguments to the wrapper script need to be provided with -w/-m prefixes
>
> ```bash
> bash /data/Ziegelbauer_lab/Pipelines/circRNA/v0.6.0/run_circrna_daq.sh -m=init -w=/scratch/circRNA_daq_test
> ```
>
> 

This assumes that `/scratch/circRNA_daq_test` does not exist before running this command and is at a location where write permissions are available.

```bash
% bash /data/Ziegelbauer_lab/Pipelines/circRNA/v0.4.0/run_circrna_daq.sh init /scratch/circRNA_daq_test
```

​	output:

```bash
Pipeline Dir: /data/Ziegelbauer_lab/circRNADetection/scripts/circRNA
Git Commit/Tag: b2f387c1b6854646d12974cd16da1168d93bb43b	v0.4.0-14-gb2f387c
Working Dir: /scratch/circRNA_daq_test
/data/Ziegelbauer_lab/circRNADetection/scripts/circRNA
/scratch/circRNA_daq_test
Logs Dir: /scratch/circRNA_daq_test/logs
Stats Dir: /scratch/circRNA_daq_test/stats
Done Initializing /scratch/circRNA_daq_test. You can now edit /scratch/circRNA_daq_test/config.yaml and /scratch/circRNA_daq_test/samples.tsv
```

The above command creates `/scratch/circRNA_daq_test` folder and creates 2 subfolders `logs` and `stats` inside that folder along with `config.yaml` and `samples.tsv` files.

```bash
% tree /scratch/circRNA_daq_test
```

​	output:

```bash
/scratch/circRNA_daq_test
├── config.yaml
├── logs
├── samples.tsv
└── stats

2 directories, 2 files
```

##### config.yaml

This file is used to fine tune the execution of the pipeline by setting:

* sample sheet location ... aka `samples.tsv`
* whether to run CLEAR pipeline or not by setting **run_clear** to `True` or `False`
* describes the location of other resources/indexes/tools etc. Generally, these do NOT need to be changed.

##### samples.tsv

Tab delimited definition of sample sheet. The header is fixed and each row represents a sample. It has 3 columns:

1. sampleName = Name of the sample. This has to be unique.
2. path_to_R1_fastq = absolute path to the read1 fastq.gz file.
3. path_to_R2_fastq = absolute path to the read2 fastq.gz file. If the sample was sequenced in single-end mode, then leave this blank.

Running **init** will put the following example file in the workdir supplied and it looks like this:

```bash
% more /scratch/circRNA_daq_test/samples.tsv
```

output:

```bash
sampleName	path_to_R1_fastq	path_to_R2_fastq
GI1_N	/data/Ziegelbauer_lab/circRNADetection/rawdata/ccbr983/fastq2/5_GI112118_norm_S4_R1_001.fastq.gz	/data/Ziegelbauer_lab/circRNADetection/rawdata/ccbr983/fastq2/5_GI112118_norm_S4_R2_001.fastq.gz
GI1_T	/data/Ziegelbauer_lab/circRNADetection/rawdata/ccbr983/fastq2/6_GI112118_tum_S5_R1_001.fastq.gz
```

#### Dryrun

Once the `samples.tsv` file has been edited appropriately to include the desired samples, it is a good idea to **dryrun** the pipeline to ensure that everything will work as desired. Dryrun can be run as follows:

```bash
% bash /data/Ziegelbauer_lab/Pipelines/circRNA/v0.4.0/run_circrna_daq.sh dryrun /scratch/circRNA_daq_test
```

[Here](dryrun_example.txt) is the output of the above command.

#### Run

Upon verifying that dryrun is successful. You can then submit the job to the cluster using the following command:

```bash
% bash /data/Ziegelbauer_lab/Pipelines/circRNA/v0.4.0/run_circrna_daq.sh run /scratch/circRNA_daq_test
```

> NOTE:
>
> With v0.6.0 or later the above command will be:
>
> ```bash
> bash /data/Ziegelbauer_lab/Pipelines/circRNA/v0.4.0/run_circrna_daq.sh -m=run -w=/scratch/circRNA_daq_test
> ```
>
> 

output:

```bash
Pipeline Dir: /data/Ziegelbauer_lab/circRNADetection/scripts/circRNA
Git Commit/Tag: 37419bc0eb196fb1e137849ebeb2739a1c12126c	v0.4.0-16-g37419bc
Working Dir: /scratch/circRNA_daq_test
[+] Loading python 3.7  ...
[+] Loading snakemake  5.24.1
Running...
ls: cannot access /scratch/circRNA_daq_test/slurm-*.out: No such file or directory
7930909
```

`7930909` is the jobid returned by the slurm job scheduler on biowulf. This means that the job was successfully submitted, it will spawn off other subjobs which in-turn will be run and outputs will be moved to the `results` folder created inside the working directory supplied at command line. You can check the status of your queue of jobs in biowulf running:

```bash
% squeue -u `whoami`
```

output:

```bash
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
           7930909  ccr,norm  circRNA kopardev PD       0:00      1 (None)
```

`ST` in the above results stands for Status and `PD` means Pending. The status will change from pending(`PD`) to running(`R`) to completed as jobs are run on the cluster.

Next, just sit tight until the pipeline finishes. You can keep monitoring the queue as shown above. If there are no jobs running on biowulf, then your pipeline has finished (or errored out!)

Once completed the output should something like this:

```bash
% tree /scratch/circRNA_daq_test
```

output:

```bash
/scratch/circRNA_daq_test
├── config.yaml
├── fastqs
│   ├── GI1_N.R1.fastq.gz -> /data/Ziegelbauer_lab/circRNADetection/rawdata/ccbr983/fastq2/5_GI112118_norm_S4_R1_001.fastq.gz
│   ├── GI1_N.R2.fastq.gz -> /data/Ziegelbauer_lab/circRNADetection/rawdata/ccbr983/fastq2/5_GI112118_norm_S4_R2_001.fastq.gz
│   └── GI1_T.R1.fastq.gz -> /data/Ziegelbauer_lab/circRNADetection/rawdata/ccbr983/fastq2/6_GI112118_tum_S5_R1_001.fastq.gz
├── logs
│   ├── 7658773.7658778.cutadapt.sample=GI1_T.err
│   ├── 7658773.7658778.cutadapt.sample=GI1_T.out
│   ├── 7658773.7658779.cutadapt.sample=GI1_N.err
│   ├── 7658773.7658779.cutadapt.sample=GI1_N.out
│   ├── 7658773.7659728.fastqc.sample=GI1_T.err
│   ├── 7658773.7659728.fastqc.sample=GI1_T.out
│   ├── 7658773.7659729.star1p.sample=GI1_T.err
│   ├── 7658773.7659729.star1p.sample=GI1_T.out
│   ├── 7658773.7659730.ciri.sample=GI1_T.err
│   ├── 7658773.7659730.ciri.sample=GI1_T.out
│   ├── 7658773.7659731.clear.sample=GI1_T.err
│   ├── 7658773.7659731.clear.sample=GI1_T.out
│   ├── 7658773.7660529.ciri.sample=GI1_N.err
│   ├── 7658773.7660529.ciri.sample=GI1_N.out
│   ├── 7658773.7660738.clear.sample=GI1_N.err
│   ├── 7658773.7660738.clear.sample=GI1_N.out
│   ├── 7658773.7660753.fastqc.sample=GI1_N.err
│   ├── 7658773.7660753.fastqc.sample=GI1_N.out
│   ├── 7658773.7660754.star1p.sample=GI1_N.err
│   ├── 7658773.7660754.star1p.sample=GI1_N.out
│   ├── 7658773.7661100.merge_SJ_tabs..err
│   ├── 7658773.7661100.merge_SJ_tabs..out
│   ├── 7658773.7661134.star2p.sample=GI1_T.err
│   ├── 7658773.7661134.star2p.sample=GI1_T.out
│   ├── 7658773.7661135.star2p.sample=GI1_N.err
│   ├── 7658773.7661135.star2p.sample=GI1_N.out
│   ├── 7658773.7662317.create_BSJ_bam.sample=GI1_T.err
│   ├── 7658773.7662317.create_BSJ_bam.sample=GI1_T.out
│   ├── 7658773.7662319.annotate_circRNA.sample=GI1_T.err
│   ├── 7658773.7662319.annotate_circRNA.sample=GI1_T.out
│   ├── 7658773.7662563.annotate_circRNA.sample=GI1_N.err
│   ├── 7658773.7662563.annotate_circRNA.sample=GI1_N.out
│   ├── 7658773.7662564.create_BSJ_bam.sample=GI1_N.err
│   ├── 7658773.7662564.create_BSJ_bam.sample=GI1_N.out
│   ├── 7658773.7662567.split_BAM_create_BW.sample=GI1_T.err
│   ├── 7658773.7662567.split_BAM_create_BW.sample=GI1_T.out
│   ├── 7658773.7663000.create_circexplorer_count_matrix..err
│   ├── 7658773.7663000.create_circexplorer_count_matrix..out
│   ├── 7658773.7663642.create_ciri_count_matrix..err
│   ├── 7658773.7663642.create_ciri_count_matrix..out
│   ├── 7658773.7663643.split_BAM_create_BW.sample=GI1_N.err
│   ├── 7658773.7663643.split_BAM_create_BW.sample=GI1_N.out
│   └── slurm-7658773.out.gz
├── qc
│   └── fastqc
│       ├── GI1_N.R1_fastqc.html
│       ├── GI1_N.R1_fastqc.zip
│       ├── GI1_N.R1.trim_fastqc.html
│       ├── GI1_N.R1.trim_fastqc.zip
│       ├── GI1_N.R2_fastqc.html
│       ├── GI1_N.R2_fastqc.zip
│       ├── GI1_N.R2.trim_fastqc.html
│       ├── GI1_N.R2.trim_fastqc.zip
│       ├── GI1_T.R1_fastqc.html
│       ├── GI1_T.R1_fastqc.zip
│       ├── GI1_T.R1.trim_fastqc.html
│       └── GI1_T.R1.trim_fastqc.zip
├── report.html
├── results
│   ├── circExplorer_BSJ_count_matrix.txt
│   ├── circExplorer_BSJ_count_matrix_with_annotations.txt
│   ├── circExplorer_count_matrix.txt
│   ├── circExplorer_count_matrix_with_annotations.txt
│   ├── ciri_count_matrix.txt
│   ├── ciri_count_matrix_with_annotations.txt
│   ├── GI1_N
│   │   ├── circExplorer
│   │   │   ├── GI1_N.back_spliced_junction.bed
│   │   │   ├── GI1_N_circexplorer_parse.log
│   │   │   ├── GI1_N.circularRNA_known.txt
│   │   │   └── low_conf_GI1_N.circularRNA_known.txt
│   │   ├── ciri
│   │   │   ├── CIRIerror.log
│   │   │   ├── GI1_N.bwa.bam
│   │   │   ├── GI1_N.bwa.log
│   │   │   ├── GI1_N.ciri.log
│   │   │   └── GI1_N.ciri.out
│   │   ├── CLEAR
│   │   │   ├── circ
│   │   │   │   └── bsj.bed
│   │   │   ├── fusion
│   │   │   │   └── junctions.bed
│   │   │   ├── hisat
│   │   │   │   └── sp.txt
│   │   │   └── quant
│   │   │       └── quant.txt
│   │   ├── STAR1p
│   │   │   ├── GI1_N_p1.Chimeric.out.junction
│   │   │   ├── GI1_N_p1.Log.final.out
│   │   │   ├── GI1_N_p1.Log.out
│   │   │   ├── GI1_N_p1.Log.progress.out
│   │   │   └── GI1_N_p1.SJ.out.tab
│   │   ├── STAR2p
│   │   │   ├── GI1_N.BSJ.bam
│   │   │   ├── GI1_N.BSJ.bam.bai
│   │   │   ├── GI1_N.BSJ.ERCC.bam
│   │   │   ├── GI1_N.BSJ.ERCC.bam.bai
│   │   │   ├── GI1_N.BSJ.ERCC.bw
│   │   │   ├── GI1_N.BSJ.hg38.bam
│   │   │   ├── GI1_N.BSJ.hg38.bam.bai
│   │   │   ├── GI1_N.BSJ.hg38.bw
│   │   │   ├── GI1_N.BSJ.MN485971.1.bam
│   │   │   ├── GI1_N.BSJ.MN485971.1.bam.bai
│   │   │   ├── GI1_N.BSJ.MN485971.1.bw
│   │   │   ├── GI1_N.BSJ.NC_000898.1.bam
│   │   │   ├── GI1_N.BSJ.NC_000898.1.bam.bai
│   │   │   ├── GI1_N.BSJ.NC_000898.1.bw
│   │   │   ├── GI1_N.BSJ.NC_001664.4.bam
│   │   │   ├── GI1_N.BSJ.NC_001664.4.bam.bai
│   │   │   ├── GI1_N.BSJ.NC_001664.4.bw
│   │   │   ├── GI1_N.BSJ.NC_001716.2.bam
│   │   │   ├── GI1_N.BSJ.NC_001716.2.bam.bai
│   │   │   ├── GI1_N.BSJ.NC_001716.2.bw
│   │   │   ├── GI1_N.BSJ.NC_006273.2.bam
│   │   │   ├── GI1_N.BSJ.NC_006273.2.bam.bai
│   │   │   ├── GI1_N.BSJ.NC_006273.2.bw
│   │   │   ├── GI1_N.BSJ.NC_007605.1.bam
│   │   │   ├── GI1_N.BSJ.NC_007605.1.bam.bai
│   │   │   ├── GI1_N.BSJ.NC_007605.1.bw
│   │   │   ├── GI1_N.BSJ.NC_009333.1.bam
│   │   │   ├── GI1_N.BSJ.NC_009333.1.bam.bai
│   │   │   ├── GI1_N.BSJ.NC_009333.1.bw
│   │   │   ├── GI1_N.BSJ.NC_045512.2.bam
│   │   │   ├── GI1_N.BSJ.NC_045512.2.bam.bai
│   │   │   ├── GI1_N.BSJ.NC_045512.2.bw
│   │   │   ├── GI1_N.BSJ.readids
│   │   │   ├── GI1_N.BSJ.rRNA.bam
│   │   │   ├── GI1_N.BSJ.rRNA.bam.bai
│   │   │   ├── GI1_N.BSJ.rRNA.bw
│   │   │   ├── GI1_N_p2.Aligned.sortedByCoord.out.bam
│   │   │   ├── GI1_N_p2.Aligned.sortedByCoord.out.bam.bai
│   │   │   ├── GI1_N_p2.Chimeric.out.junction
│   │   │   ├── GI1_N_p2.Log.final.out
│   │   │   ├── GI1_N_p2.Log.out
│   │   │   ├── GI1_N_p2.Log.progress.out
│   │   │   └── GI1_N_p2.SJ.out.tab
│   │   └── trim
│   │       ├── GI1_N.R1.trim.fastq.gz
│   │       └── GI1_N.R2.trim.fastq.gz
│   ├── GI1_T
│   │   ├── circExplorer
│   │   │   ├── GI1_T.back_spliced_junction.bed
│   │   │   ├── GI1_T_circexplorer_parse.log
│   │   │   ├── GI1_T.circularRNA_known.txt
│   │   │   └── low_conf_GI1_T.circularRNA_known.txt
│   │   ├── ciri
│   │   │   ├── CIRIerror.log
│   │   │   ├── GI1_T.bwa.bam
│   │   │   ├── GI1_T.bwa.log
│   │   │   ├── GI1_T.ciri.log
│   │   │   └── GI1_T.ciri.out
│   │   ├── CLEAR
│   │   │   ├── circ
│   │   │   │   └── bsj.bed
│   │   │   ├── fusion
│   │   │   │   └── junctions.bed
│   │   │   ├── hisat
│   │   │   │   └── sp.txt
│   │   │   └── quant
│   │   │       └── quant.txt
│   │   ├── STAR1p
│   │   │   ├── GI1_T_p1.Chimeric.out.junction
│   │   │   ├── GI1_T_p1.Log.final.out
│   │   │   ├── GI1_T_p1.Log.out
│   │   │   ├── GI1_T_p1.Log.progress.out
│   │   │   └── GI1_T_p1.SJ.out.tab
│   │   ├── STAR2p
│   │   │   ├── GI1_T.BSJ.bam
│   │   │   ├── GI1_T.BSJ.bam.bai
│   │   │   ├── GI1_T.BSJ.ERCC.bam
│   │   │   ├── GI1_T.BSJ.ERCC.bam.bai
│   │   │   ├── GI1_T.BSJ.ERCC.bw
│   │   │   ├── GI1_T.BSJ.hg38.bam
│   │   │   ├── GI1_T.BSJ.hg38.bam.bai
│   │   │   ├── GI1_T.BSJ.hg38.bw
│   │   │   ├── GI1_T.BSJ.MN485971.1.bam
│   │   │   ├── GI1_T.BSJ.MN485971.1.bam.bai
│   │   │   ├── GI1_T.BSJ.MN485971.1.bw
│   │   │   ├── GI1_T.BSJ.NC_000898.1.bam
│   │   │   ├── GI1_T.BSJ.NC_000898.1.bam.bai
│   │   │   ├── GI1_T.BSJ.NC_000898.1.bw
│   │   │   ├── GI1_T.BSJ.NC_001664.4.bam
│   │   │   ├── GI1_T.BSJ.NC_001664.4.bam.bai
│   │   │   ├── GI1_T.BSJ.NC_001664.4.bw
│   │   │   ├── GI1_T.BSJ.NC_001716.2.bam
│   │   │   ├── GI1_T.BSJ.NC_001716.2.bam.bai
│   │   │   ├── GI1_T.BSJ.NC_001716.2.bw
│   │   │   ├── GI1_T.BSJ.NC_006273.2.bam
│   │   │   ├── GI1_T.BSJ.NC_006273.2.bam.bai
│   │   │   ├── GI1_T.BSJ.NC_006273.2.bw
│   │   │   ├── GI1_T.BSJ.NC_007605.1.bam
│   │   │   ├── GI1_T.BSJ.NC_007605.1.bam.bai
│   │   │   ├── GI1_T.BSJ.NC_007605.1.bw
│   │   │   ├── GI1_T.BSJ.NC_009333.1.bam
│   │   │   ├── GI1_T.BSJ.NC_009333.1.bam.bai
│   │   │   ├── GI1_T.BSJ.NC_009333.1.bw
│   │   │   ├── GI1_T.BSJ.NC_045512.2.bam
│   │   │   ├── GI1_T.BSJ.NC_045512.2.bam.bai
│   │   │   ├── GI1_T.BSJ.NC_045512.2.bw
│   │   │   ├── GI1_T.BSJ.readids
│   │   │   ├── GI1_T.BSJ.rRNA.bam
│   │   │   ├── GI1_T.BSJ.rRNA.bam.bai
│   │   │   ├── GI1_T.BSJ.rRNA.bw
│   │   │   ├── GI1_T_p2.Aligned.sortedByCoord.out.bam
│   │   │   ├── GI1_T_p2.Aligned.sortedByCoord.out.bam.bai
│   │   │   ├── GI1_T_p2.Chimeric.out.junction
│   │   │   ├── GI1_T_p2.Log.final.out
│   │   │   ├── GI1_T_p2.Log.out
│   │   │   ├── GI1_T_p2.Log.progress.out
│   │   │   └── GI1_T_p2.SJ.out.tab
│   │   └── trim
│   │       ├── GI1_T.R1.trim.fastq.gz
│   │       └── GI1_T.R2.trim.fastq.gz
│   └── pass1.out.tab
├── samples.tsv
├── stats
│   ├── snakemake.log.20210204103112.gz
│   └── snakemake.stats.20210204103112.gz
├── submit_script.sbatch
```

#### Runlocal

If you have grabbed an interactive node using `sinteractive` and have a small test dataset in the `samples.tsv` and simply want to quickly check if everything works as expected. You can run **locally**, i.e. directly on the interactive node without submitting to the cluster using:

```bash
% bash /data/Ziegelbauer_lab/Pipelines/circRNA/v0.4.0/run_circrna_daq.sh runlocal /scratch/circRNA_daq_test
```

This is only for testing purposes, do not use it for running 10s of samples as you will be soon timed out of the interactive node.

