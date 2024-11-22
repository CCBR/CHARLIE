# Tutorial

#### Prerequisites

- [Biowulf](https://hpc.nih.gov/) account: Biowulf account can be requested [here](https://hpc.nih.gov/docs/accounts.html).

#### Installation

CHARLIE is already installed on biowulf.
It is included in the ccbrpipeliner module from release 7 onward.
To load the module run:

```bash
module load ccbrpipeliner/7
```

#### Init

To get help about the pipeline you can run:

```bash
charlie --help
```

output:

```bash

##########################################################################################

Welcome to charlie(v0.10.0-dev)
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

CHARLIE can be used to DAQ (Detect/Annotate/Quantify) circRNAs in hosts and viruses.

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
  charlie -w/--workdir=<WORKDIR> -m/--runmode=<RUNMODE>

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

--singcache|-c  : singularity cache directory. Default is `/data/${USER}/.singularity` if available, or falls back to `${WORKDIR}/.singularity`. Use this flag to specify a different singularity cache directory.
--host|-g       : supply host at command line. hg38 or mm39.                                            (--runmode=init only)
--additives|-a  : supply comma-separated list of additives at command line. ERCC or BAC16Insert or both (--runmode=init only)
--viruses|-v    : supply comma-separated list of viruses at command line                                (--runmode=init only)
--manifest|-s   : absolute path to samples.tsv. This will be copied to output folder                    (--runmode=init only)
--changegrp|-z  : change group to "Ziegelbauer_lab" before running anything. Biowulf-only. Useful for correctly setting permissions.
--help|-h       : print this help


Example commands:
  charlie -w=/my/output/folder -m=init
  charlie -w=/my/output/folder -m=dryrun
  charlie -w=/my/output/folder -m=run

##########################################################################################

VersionInfo:
  python          : 3
  snakemake       : 7
  pipeline_home   : /gpfs/gsfs10/users/CCBR_Pipeliner/Pipelines/CHARLIE/.v0.11.1
  git commit/tag  : 613fb617f1ed426fb8900f98e599ca0497a67cc0    v0.11.0-49-g613fb61

##########################################################################################
```

> NOTE:
> You can replace `v0.10.0` in the above command with the latest version to use a newer version. `run_circrna_daq.sh` was called `test.sh` in versions older than `v0.4.0`.

To initial the working directory run:

```bash
charlie -w=<path to output dir> -m=init
```

This assumes that `<path to output dir>` does not exist before running the above command and is at a location where write permissions are available.

The above command creates `<path to output dir>` folder and creates 2 subfolders `logs` and `stats` inside that folder along with `config.yaml` and `samples.tsv` files.

```bash
tree <path to output dir>
```

##### config.yaml

This file is used to fine tune the execution of the pipeline by setting:

- sample sheet location ... aka `samples.tsv`
- the temporary directory -- make sure this is correct for your computing environment.
- which circRNA finding tools to use by editing these:
  - run_clear: True
  - run_dcc: True
  - run_mapsplice: False
  - run_circRNAFinder: True
  - run_nclscan: False
  - run_findcirc: False
- describes the location of other resources/indexes/tools etc. Generally, these do NOT need to be changed.

##### samples.tsv

Tab delimited definition of sample sheet. The header is fixed and each row represents a sample. It has 3 columns:

1. sampleName = Name of the sample. This has to be unique.
2. path_to_R1_fastq = absolute path to the read1 fastq.gz file.
3. path_to_R2_fastq = absolute path to the read2 fastq.gz file. If the sample was sequenced in single-end mode, then leave this blank.

The `/data/CCBR_Pipeliner/testdata/circRNA/humans` folder in the repo has test dataset:

```bash
tree /data/CCBR_Pipeliner/testdata/circRNA/humans
/data/CCBR_Pipeliner/testdata/circRNA/humans
├── GI1_N_ss.R1.fastq.gz
├── GI1_N_ss.R2.fastq.gz
├── GI1_T_ss.R1.fastq.gz
└── samples.tsv
```

`GI1_N` is a PE sample while `GI1_T` is a SE sample.

#### Dryrun

Once the `samples.tsv` file has been edited appropriately to include the desired samples, it is a good idea to **dryrun** the pipeline to ensure that everything will work as desired. Dryrun can be run as follows:

```bash
charlie -w=<path to output dir> -m=dryrun
```

This will create the reference fasta and gtf file based on the selections made in the `config.yaml`. Hence, can take a few minutes to run.

#### Run

Upon verifying that dryrun is successful. You can then submit the job to the cluster using the following command:

```bash
charlie -w=<path to output dir> -m=run
```

which will produce something like this:

```
... ... skipping ~1000 lines
...
...
Job stats:
job                                              count    min threads    max threads
---------------------------------------------  -------  ----------
all                                                  1              1              1
annotate_clear_output                                2              1              1
circExplorer                                         2              2              2
circExplorer_bwa                                     2              2              2
circrnafinder                                        2              1              1
ciri                                                 2             56             56
clear                                                2              2              2
create_bowtie2_index                                 1              1              1
create_bwa_index                                     1              1              1
create_circExplorer_BSJ_bam                          2              4              4
create_circExplorer_linear_spliced_bams              2             56             56
create_circExplorer_merged_found_counts_table        2              1              1
create_hq_bams                                       2              1              1
create_index                                         1             56             56
create_master_counts_file                            1              1              1
cutadapt                                             2             56             56
dcc                                                  2              4              4
dcc_create_samplesheets                              2              1              1
estimate_duplication                                 2              1              1
fastqc                                               2              4              4
find_circ                                            2             56             56
find_circ_align                                      2             56             56
merge_SJ_tabs                                        1              2              2
merge_alignment_stats                                1              1              1
merge_genecounts                                     1              1              1
merge_per_sample                                     2              1              1
star1p                                               2             56             56
star2p                                               2             56             56
star_circrnafinder                                   2             56             56
total                                               52              1             56

Reasons:
    (check individual jobs above for details)
    input files updated by another job:
        alignment_stats, all, annotate_clear_output, circExplorer, circExplorer_bwa, circrnafinder, ciri, clear, create_circExplorer_BSJ_bam, create_circExplorer_linear_spliced_bams, create_circExplorer_merged_found_counts_table, create_hq_bams, create_master_counts_file, dcc, dcc_create_samplesheets, estimate_duplication, fastqc, find_circ, find_circ_align, merge_SJ_tabs, merge_alignment_stats, merge_genecounts, merge_per_sample, star1p, star2p, star_circrnafinder
    missing output files:
        alignment_stats, annotate_clear_output, circExplorer, circExplorer_bwa, circrnafinder, ciri, clear, create_bowtie2_index, create_bwa_index, create_circExplorer_BSJ_bam, create_circExplorer_linear_spliced_bams, create_circExplorer_merged_found_counts_table, create_hq_bams, create_index, create_master_counts_file, cutadapt, dcc, dcc_create_samplesheets, estimate_duplication, fastqc, find_circ, find_circ_align, merge_SJ_tabs, merge_alignment_stats, merge_genecounts, merge_per_sample, star1p, star2p, star_circrnafinder

This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.
Running...
14743440
```

In this example, `14743440` is the jobid returned by the slurm job scheduler on biowulf. This means that the job was successfully submitted, it will spawn off other subjobs which in-turn will be run and outputs will be moved to the `results` folder created inside the working directory supplied at command line. You can check the status of your queue of jobs in biowulf running:

```bash
squeue -u `whoami`
```

output:

```bash
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
           14743440  ccr,norm  circRNA kopardev PD       0:00      1 (None)
```

`ST` in the above results stands for Status and `PD` means Pending. The status will change from pending(`PD`) to running(`R`) to completed as jobs are run on the cluster.

Next, just sit tight until the pipeline finishes. You can keep monitoring the queue as shown above. If there are no jobs running on biowulf, then your pipeline has finished (or errored out!)

Once completed the output should something like this:

```bash
tree <path to output dir>
```

output:

```bash
.
├── cluster.json
├── config.yaml
├── dryrun.231222103505.log
├── fastqs
│   ├── GI1_N.R1.fastq.gz -> /data/CCBR_Pipeliner/testdata/circRNA/human/GI1_N_ss.R1.fastq.gz
│   ├── GI1_N.R2.fastq.gz -> /data/CCBR_Pipeliner/testdata/circRNA/human/GI1_N_ss.R2.fastq.gz
│   └── GI1_T.R1.fastq.gz -> /data/CCBR_Pipeliner/testdata/circRNA/human/GI1_T_ss.R1.fastq.gz
├── logs/snakemake.log.jobby.short
├── logs/snakemake.log.jobby.txt
├── logs
│   ... log files ...
│   ... skipping ...
├── nclscan.config
├── qc
│   ├── fastqc
│   │   ├── GI1_N.R1_fastqc.html
│   │   ├── GI1_N.R1_fastqc.zip
│   │   ├── GI1_N.R1.trim_fastqc.html
│   │   ├── GI1_N.R1.trim_fastqc.zip
│   │   ├── GI1_N.R2_fastqc.html
│   │   ├── GI1_N.R2_fastqc.zip
│   │   ├── GI1_N.R2.trim_fastqc.html
│   │   ├── GI1_N.R2.trim_fastqc.zip
│   │   ├── GI1_T.R1_fastqc.html
│   │   ├── GI1_T.R1_fastqc.zip
│   │   ├── GI1_T.R1.trim_fastqc.html
│   │   └── GI1_T.R1.trim_fastqc.zip
│   └── picard_MarkDuplicates
│       ├── GI1_N.MarkDuplicates.metrics.txt
│       └── GI1_T.MarkDuplicates.metrics.txt
├── ref
│   ├── bwa_index.log
│   ├── gene_id_2_gene_name.tsv
│   ├── NCLscan_index
│   │   ├── AllRef.fa
│   │   ├── AllRef.fa.amb
│   │   ├── AllRef.fa.ann
│   │   ├── AllRef.fa.bwt
│   │   ├── AllRef.fa.pac
│   │   ├── AllRef.fa.sa
│   │   ├── AllRef.ndx
│   │   └── RepChrM.fa
│   ├── ref.1.bt2
│   ├── ref.2.bt2
│   ├── ref.3.bt2
│   ├── ref.4.bt2
│   ├── ref.amb
│   ├── ref.ann
│   ├── ref.bwt
│   ├── ref.dummy.fa
│   ├── ref.fa
│   ├── ref.fa.byo_index
│   ├── ref.fa.fai
│   ├── ref.fa.regions
│   ├── ref.fa.regions.host
│   ├── ref.fa.regions.viruses
│   ├── ref.fa.sizes
│   ├── ref.fixed.gtf
│   ├── ref.genes.genepred
│   ├── ref.genes.genepred_w_geneid
│   ├── ref.gtf
│   ├── ref.pac
│   ├── ref.rev.1.bt2
│   ├── ref.rev.2.bt2
│   ├── ref.sa
│   ├── ref.transcripts.fa
│   ├── separate_fastas
│   │   ├── chr10.fa
│   │   ├── chr11.fa
│   │   ├── chr12.fa
│   │   ├── chr13.fa
│   │   ├── chr14.fa
│   │   ├── chr15.fa
│   │   ├── chr16.fa
│   │   ├── chr17.fa
│   │   ├── chr18.fa
│   │   ├── chr19.fa
│   │   ├── chr1.fa
│   │   ├── chr20.fa
│   │   ├── chr21.fa
│   │   ├── chr22.fa
│   │   ├── chr2.fa
│   │   ├── chr3.fa
│   │   ├── chr45S.fa
│   │   ├── chr4.fa
│   │   ├── chr5.fa
│   │   ├── chr5S.fa
│   │   ├── chr6.fa
│   │   ├── chr7.fa
│   │   ├── chr8.fa
│   │   ├── chr9.fa
│   │   ├── chrM.fa
│   │   ├── chrX.fa
│   │   ├── chrY.fa
│   │   ├── ERCC_00002.fa
│   │   ├── ERCC_00003.fa
│   │   ├── ERCC_00004.fa
│   │   ├── ERCC_00009.fa
│   │   ├── ERCC_00012.fa
│   │   ├── ERCC_00013.fa
│   │   ├── ERCC_00014.fa
│   │   ├── ERCC_00016.fa
│   │   ├── ERCC_00017.fa
│   │   ├── ERCC_00019.fa
│   │   ├── ERCC_00022.fa
│   │   ├── ERCC_00024.fa
│   │   ├── ERCC_00025.fa
│   │   ├── ERCC_00028.fa
│   │   ├── ERCC_00031.fa
│   │   ├── ERCC_00033.fa
│   │   ├── ERCC_00034.fa
│   │   ├── ERCC_00035.fa
│   │   ├── ERCC_00039.fa
│   │   ├── ERCC_00040.fa
│   │   ├── ERCC_00041.fa
│   │   ├── ERCC_00042.fa
│   │   ├── ERCC_00043.fa
│   │   ├── ERCC_00044.fa
│   │   ├── ERCC_00046.fa
│   │   ├── ERCC_00048.fa
│   │   ├── ERCC_00051.fa
│   │   ├── ERCC_00053.fa
│   │   ├── ERCC_00054.fa
│   │   ├── ERCC_00057.fa
│   │   ├── ERCC_00058.fa
│   │   ├── ERCC_00059.fa
│   │   ├── ERCC_00060.fa
│   │   ├── ERCC_00061.fa
│   │   ├── ERCC_00062.fa
│   │   ├── ERCC_00067.fa
│   │   ├── ERCC_00069.fa
│   │   ├── ERCC_00071.fa
│   │   ├── ERCC_00073.fa
│   │   ├── ERCC_00074.fa
│   │   ├── ERCC_00075.fa
│   │   ├── ERCC_00076.fa
│   │   ├── ERCC_00077.fa
│   │   ├── ERCC_00078.fa
│   │   ├── ERCC_00079.fa
│   │   ├── ERCC_00081.fa
│   │   ├── ERCC_00083.fa
│   │   ├── ERCC_00084.fa
│   │   ├── ERCC_00085.fa
│   │   ├── ERCC_00086.fa
│   │   ├── ERCC_00092.fa
│   │   ├── ERCC_00095.fa
│   │   ├── ERCC_00096.fa
│   │   ├── ERCC_00097.fa
│   │   ├── ERCC_00098.fa
│   │   ├── ERCC_00099.fa
│   │   ├── ERCC_00104.fa
│   │   ├── ERCC_00108.fa
│   │   ├── ERCC_00109.fa
│   │   ├── ERCC_00111.fa
│   │   ├── ERCC_00112.fa
│   │   ├── ERCC_00113.fa
│   │   ├── ERCC_00116.fa
│   │   ├── ERCC_00117.fa
│   │   ├── ERCC_00120.fa
│   │   ├── ERCC_00123.fa
│   │   ├── ERCC_00126.fa
│   │   ├── ERCC_00130.fa
│   │   ├── ERCC_00131.fa
│   │   ├── ERCC_00134.fa
│   │   ├── ERCC_00136.fa
│   │   ├── ERCC_00137.fa
│   │   ├── ERCC_00138.fa
│   │   ├── ERCC_00142.fa
│   │   ├── ERCC_00143.fa
│   │   ├── ERCC_00144.fa
│   │   ├── ERCC_00145.fa
│   │   ├── ERCC_00147.fa
│   │   ├── ERCC_00148.fa
│   │   ├── ERCC_00150.fa
│   │   ├── ERCC_00154.fa
│   │   ├── ERCC_00156.fa
│   │   ├── ERCC_00157.fa
│   │   ├── ERCC_00158.fa
│   │   ├── ERCC_00160.fa
│   │   ├── ERCC_00162.fa
│   │   ├── ERCC_00163.fa
│   │   ├── ERCC_00164.fa
│   │   ├── ERCC_00165.fa
│   │   ├── ERCC_00168.fa
│   │   ├── ERCC_00170.fa
│   │   ├── ERCC_00171.fa
│   │   ├── GL000008.2.fa
│   │   ├── GL000009.2.fa
│   │   ├── GL000194.1.fa
│   │   ├── GL000195.1.fa
│   │   ├── GL000205.2.fa
│   │   ├── GL000208.1.fa
│   │   ├── GL000213.1.fa
│   │   ├── GL000214.1.fa
│   │   ├── GL000216.2.fa
│   │   ├── GL000218.1.fa
│   │   ├── GL000219.1.fa
│   │   ├── GL000220.1.fa
│   │   ├── GL000221.1.fa
│   │   ├── GL000224.1.fa
│   │   ├── GL000225.1.fa
│   │   ├── GL000226.1.fa
│   │   ├── KI270302.1.fa
│   │   ├── KI270303.1.fa
│   │   ├── KI270304.1.fa
│   │   ├── KI270305.1.fa
│   │   ├── KI270310.1.fa
│   │   ├── KI270311.1.fa
│   │   ├── KI270312.1.fa
│   │   ├── KI270315.1.fa
│   │   ├── KI270316.1.fa
│   │   ├── KI270317.1.fa
│   │   ├── KI270320.1.fa
│   │   ├── KI270322.1.fa
│   │   ├── KI270329.1.fa
│   │   ├── KI270330.1.fa
│   │   ├── KI270333.1.fa
│   │   ├── KI270334.1.fa
│   │   ├── KI270335.1.fa
│   │   ├── KI270336.1.fa
│   │   ├── KI270337.1.fa
│   │   ├── KI270338.1.fa
│   │   ├── KI270340.1.fa
│   │   ├── KI270362.1.fa
│   │   ├── KI270363.1.fa
│   │   ├── KI270364.1.fa
│   │   ├── KI270366.1.fa
│   │   ├── KI270371.1.fa
│   │   ├── KI270372.1.fa
│   │   ├── KI270373.1.fa
│   │   ├── KI270374.1.fa
│   │   ├── KI270375.1.fa
│   │   ├── KI270376.1.fa
│   │   ├── KI270378.1.fa
│   │   ├── KI270379.1.fa
│   │   ├── KI270381.1.fa
│   │   ├── KI270382.1.fa
│   │   ├── KI270383.1.fa
│   │   ├── KI270384.1.fa
│   │   ├── KI270385.1.fa
│   │   ├── KI270386.1.fa
│   │   ├── KI270387.1.fa
│   │   ├── KI270388.1.fa
│   │   ├── KI270389.1.fa
│   │   ├── KI270390.1.fa
│   │   ├── KI270391.1.fa
│   │   ├── KI270392.1.fa
│   │   ├── KI270393.1.fa
│   │   ├── KI270394.1.fa
│   │   ├── KI270395.1.fa
│   │   ├── KI270396.1.fa
│   │   ├── KI270411.1.fa
│   │   ├── KI270412.1.fa
│   │   ├── KI270414.1.fa
│   │   ├── KI270417.1.fa
│   │   ├── KI270418.1.fa
│   │   ├── KI270419.1.fa
│   │   ├── KI270420.1.fa
│   │   ├── KI270422.1.fa
│   │   ├── KI270423.1.fa
│   │   ├── KI270424.1.fa
│   │   ├── KI270425.1.fa
│   │   ├── KI270429.1.fa
│   │   ├── KI270435.1.fa
│   │   ├── KI270438.1.fa
│   │   ├── KI270442.1.fa
│   │   ├── KI270448.1.fa
│   │   ├── KI270465.1.fa
│   │   ├── KI270466.1.fa
│   │   ├── KI270467.1.fa
│   │   ├── KI270468.1.fa
│   │   ├── KI270507.1.fa
│   │   ├── KI270508.1.fa
│   │   ├── KI270509.1.fa
│   │   ├── KI270510.1.fa
│   │   ├── KI270511.1.fa
│   │   ├── KI270512.1.fa
│   │   ├── KI270515.1.fa
│   │   ├── KI270516.1.fa
│   │   ├── KI270517.1.fa
│   │   ├── KI270518.1.fa
│   │   ├── KI270519.1.fa
│   │   ├── KI270521.1.fa
│   │   ├── KI270522.1.fa
│   │   ├── KI270528.1.fa
│   │   ├── KI270529.1.fa
│   │   ├── KI270530.1.fa
│   │   ├── KI270538.1.fa
│   │   ├── KI270539.1.fa
│   │   ├── KI270544.1.fa
│   │   ├── KI270548.1.fa
│   │   ├── KI270579.1.fa
│   │   ├── KI270580.1.fa
│   │   ├── KI270581.1.fa
│   │   ├── KI270582.1.fa
│   │   ├── KI270583.1.fa
│   │   ├── KI270584.1.fa
│   │   ├── KI270587.1.fa
│   │   ├── KI270588.1.fa
│   │   ├── KI270589.1.fa
│   │   ├── KI270590.1.fa
│   │   ├── KI270591.1.fa
│   │   ├── KI270593.1.fa
│   │   ├── KI270706.1.fa
│   │   ├── KI270707.1.fa
│   │   ├── KI270708.1.fa
│   │   ├── KI270709.1.fa
│   │   ├── KI270710.1.fa
│   │   ├── KI270711.1.fa
│   │   ├── KI270712.1.fa
│   │   ├── KI270713.1.fa
│   │   ├── KI270714.1.fa
│   │   ├── KI270715.1.fa
│   │   ├── KI270716.1.fa
│   │   ├── KI270717.1.fa
│   │   ├── KI270718.1.fa
│   │   ├── KI270719.1.fa
│   │   ├── KI270720.1.fa
│   │   ├── KI270721.1.fa
│   │   ├── KI270722.1.fa
│   │   ├── KI270723.1.fa
│   │   ├── KI270724.1.fa
│   │   ├── KI270725.1.fa
│   │   ├── KI270726.1.fa
│   │   ├── KI270727.1.fa
│   │   ├── KI270728.1.fa
│   │   ├── KI270729.1.fa
│   │   ├── KI270730.1.fa
│   │   ├── KI270731.1.fa
│   │   ├── KI270732.1.fa
│   │   ├── KI270733.1.fa
│   │   ├── KI270734.1.fa
│   │   ├── KI270735.1.fa
│   │   ├── KI270736.1.fa
│   │   ├── KI270737.1.fa
│   │   ├── KI270738.1.fa
│   │   ├── KI270739.1.fa
│   │   ├── KI270740.1.fa
│   │   ├── KI270741.1.fa
│   │   ├── KI270742.1.fa
│   │   ├── KI270743.1.fa
│   │   ├── KI270744.1.fa
│   │   ├── KI270745.1.fa
│   │   ├── KI270746.1.fa
│   │   ├── KI270747.1.fa
│   │   ├── KI270748.1.fa
│   │   ├── KI270749.1.fa
│   │   ├── KI270750.1.fa
│   │   ├── KI270751.1.fa
│   │   ├── KI270752.1.fa
│   │   ├── KI270753.1.fa
│   │   ├── KI270754.1.fa
│   │   ├── KI270755.1.fa
│   │   ├── KI270756.1.fa
│   │   ├── KI270757.1.fa
│   │   ├── NC_009333.1.fa
│   │   └── separate_fastas.lst
│   └── STAR_no_GTF
│       ├── chrLength.txt
│       ├── chrNameLength.txt
│       ├── chrName.txt
│       ├── chrStart.txt
│       ├── Genome
│       ├── genomeParameters.txt
│       ├── Log.out
│       ├── SA
│       └── SAindex
├── results
│   ├── alignmentstats.txt
│   ├── circRNA_master_counts.tsv.gz
│   ├── GI1_N
│   │   ├── alignmentstats.txt
│   │   ├── circExplorer
│   │   │   ├── GI1_N.back_spliced_junction.bed
│   │   │   ├── GI1_N.back_spliced_junction.strand_fixed.bed
│   │   │   ├── GI1_N.bam
│   │   │   ├── GI1_N.bam.bai
│   │   │   ├── GI1_N.BSJ.bam
│   │   │   ├── GI1_N.BSJ.bam.csi
│   │   │   ├── GI1_N.BSJ.bed.gz
│   │   │   ├── GI1_N.BSJ.bw
│   │   │   ├── GI1_N.BSJ.foundcounts.tsv
│   │   │   ├── GI1_N.BSJ.minus.bam
│   │   │   ├── GI1_N.BSJ.minus.bam.csi
│   │   │   ├── GI1_N.BSJ.minus.bw
│   │   │   ├── GI1_N.BSJ.plus.bam
│   │   │   ├── GI1_N.BSJ.plus.bam.csi
│   │   │   ├── GI1_N.BSJ.plus.bw
│   │   │   ├── GI1_N.circExplorer.annotation_counts.tsv
│   │   │   ├── GI1_N.circExplorer.counts_table.tsv
│   │   │   ├── GI1_N_circexplorer_parse.log
│   │   │   ├── GI1_N.circRNA_known.filter1.txt
│   │   │   ├── GI1_N.circularRNA_known.txt
│   │   │   ├── GI1_N.hg38.BSJ.bam
│   │   │   ├── GI1_N.hg38.BSJ.bam.csi
│   │   │   ├── GI1_N.hg38.BSJ.bw
│   │   │   ├── GI1_N.linear.bam
│   │   │   ├── GI1_N.linear.bam.csi
│   │   │   ├── GI1_N.linear_BSJ.bam
│   │   │   ├── GI1_N.linear_BSJ.bam.csi
│   │   │   ├── GI1_N.linear_BSJ.bw
│   │   │   ├── GI1_N.linear_BSJ.hg38.bam
│   │   │   ├── GI1_N.linear_BSJ.hg38.bw
│   │   │   ├── GI1_N.linear_BSJ.NC_009333.1.bam
│   │   │   ├── GI1_N.linear_BSJ.NC_009333.1.bw
│   │   │   ├── GI1_N.linear.BSJ.readids.gz
│   │   │   ├── GI1_N.linear.bw
│   │   │   ├── GI1_N.linear.hg38.bam
│   │   │   ├── GI1_N.linear.hg38.bw
│   │   │   ├── GI1_N.linear.NC_009333.1.bam
│   │   │   ├── GI1_N.linear.NC_009333.1.bw
│   │   │   ├── GI1_N.linear_spliced.counts.tsv
│   │   │   ├── GI1_N.NC_009333.1.BSJ.bam
│   │   │   ├── GI1_N.NC_009333.1.BSJ.bam.csi
│   │   │   ├── GI1_N.NC_009333.1.BSJ.bw
│   │   │   ├── GI1_N.readcounts.tsv
│   │   │   ├── GI1_N.rid2jid.tsv.gz
│   │   │   ├── GI1_N.spliced.bam
│   │   │   ├── GI1_N.spliced.bam.csi
│   │   │   ├── GI1_N.spliced_BSJ.bam
│   │   │   ├── GI1_N.spliced_BSJ.bam.csi
│   │   │   ├── GI1_N.spliced_BSJ.bw
│   │   │   ├── GI1_N.spliced_BSJ.hg38.bam
│   │   │   ├── GI1_N.spliced_BSJ.hg38.bw
│   │   │   ├── GI1_N.spliced_BSJ.NC_009333.1.bam
│   │   │   ├── GI1_N.spliced_BSJ.NC_009333.1.bw
│   │   │   ├── GI1_N.spliced.BSJ.readids.gz
│   │   │   ├── GI1_N.spliced.bw
│   │   │   ├── GI1_N.spliced.hg38.bam
│   │   │   ├── GI1_N.spliced.hg38.bw
│   │   │   ├── GI1_N.spliced.NC_009333.1.bam
│   │   │   ├── GI1_N.spliced.NC_009333.1.bw
│   │   │   ├── low_conf_GI1_N.circRNA_known.filter1.txt
│   │   │   └── low_conf_GI1_N.circRNA_known.txt
│   │   ├── circExplorer_BWA
│   │   │   ├── GI1_N.back_spliced_junction.bed
│   │   │   ├── GI1_N.back_spliced_junction.strand_fixed.bed
│   │   │   ├── GI1_N.circExplorer_bwa.annotation_counts.tsv
│   │   │   ├── GI1_N_circexplorer_parse.log
│   │   │   ├── GI1_N.circRNA_known.filter1.txt
│   │   │   ├── GI1_N.circularRNA_known.txt
│   │   │   ├── low_conf_GI1_N.circRNA_known.filter1.txt
│   │   │   └── low_conf_GI1_N.circRNA_known.txt
│   │   ├── circRNA_finder
│   │   │   ├── GI1_N.Chimeric.out.sorted.bam
│   │   │   ├── GI1_N.Chimeric.out.sorted.bam.bai
│   │   │   ├── GI1_N.circRNA_finder.counts_table.tsv.filtered
│   │   │   ├── GI1_N.filteredJunctions.bed
│   │   │   ├── GI1_N.s_filteredJunctions.bed
│   │   │   └── GI1_N.s_filteredJunctions_fw.bed
│   │   ├── ciri
│   │   │   ├── CIRIerror.log
│   │   │   ├── GI1_N.bwa.log
│   │   │   ├── GI1_N.ciri.bam
│   │   │   ├── GI1_N.ciri.bam.csi
│   │   │   ├── GI1_N.ciri.log
│   │   │   ├── GI1_N.ciri.out
│   │   │   └── GI1_N.ciri.out.filtered
│   │   ├── CLEAR
│   │   │   ├── quant.txt
│   │   │   └── quant.txt.annotated
│   │   ├── DCC
│   │   │   ├── CircCoordinates
│   │   │   ├── CircRNACount
│   │   │   ├── CircSkipJunctions
│   │   │   ├── DCC-2023-12-20_1149.log
│   │   │   ├── GI1_N.dcc.counts_table.tsv
│   │   │   ├── GI1_N.dcc.counts_table.tsv.filtered
│   │   │   ├── LinearCount
│   │   │   ├── mate1.txt
│   │   │   ├── mate2.txt
│   │   │   └── samplesheet.txt
│   │   ├── find_circ
│   │   │   ├── GI1_N_anchors.fastq.gz
│   │   │   ├── GI1_N.unmapped.bam
│   │   │   └── GI1_N.unmapped.bam.csi
│   │   ├── GI1_N.circRNA_counts.txt.gz
│   │   ├── merge_per_sample.sh
│   │   ├── STAR1p
│   │   │   ├── GI1_N_p1.Chimeric.out.junction
│   │   │   ├── GI1_N_p1.Log.final.out
│   │   │   ├── GI1_N_p1.Log.out
│   │   │   ├── GI1_N_p1.Log.progress.out
│   │   │   ├── GI1_N_p1.SJ.out.tab
│   │   │   ├── mate1
│   │   │   │   ├── GI1_N_mate1.Chimeric.out.junction
│   │   │   │   ├── GI1_N_mate1.Log.final.out
│   │   │   │   ├── GI1_N_mate1.Log.out
│   │   │   │   ├── GI1_N_mate1.Log.progress.out
│   │   │   │   └── GI1_N_mate1.SJ.out.tab
│   │   │   └── mate2
│   │   │       ├── GI1_N_mate2.Chimeric.out.junction
│   │   │       ├── GI1_N_mate2.Chimeric.out.junction.fixed
│   │   │       ├── GI1_N_mate2.Log.final.out
│   │   │       ├── GI1_N_mate2.Log.out
│   │   │       ├── GI1_N_mate2.Log.progress.out
│   │   │       └── GI1_N_mate2.SJ.out.tab
│   │   ├── STAR2p
│   │   │   ├── GI1_N_p2.bam
│   │   │   ├── GI1_N_p2.bam.csi
│   │   │   ├── GI1_N_p2.chimeric.bam
│   │   │   ├── GI1_N_p2.chimeric.bam.csi
│   │   │   ├── GI1_N_p2.Chimeric.out.junction
│   │   │   ├── GI1_N_p2.Log.final.out
│   │   │   ├── GI1_N_p2.Log.out
│   │   │   ├── GI1_N_p2.Log.progress.out
│   │   │   ├── GI1_N_p2.non_chimeric.bam
│   │   │   ├── GI1_N_p2.non_chimeric.bam.csi
│   │   │   ├── GI1_N_p2.ReadsPerGene.out.tab
│   │   │   └── GI1_N_p2.SJ.out.tab
│   │   ├── STAR_circRNAFinder
│   │   │   ├── GI1_N.Chimeric.out.junction
│   │   │   ├── GI1_N.Chimeric.out.sam
│   │   │   ├── GI1_N.Log.final.out
│   │   │   ├── GI1_N.Log.out
│   │   │   ├── GI1_N.Log.progress.out
│   │   │   └── GI1_N.SJ.out.tab
│   │   └── trim
│   │       ├── GI1_N.R1.trim.fastq.gz
│   │       └── GI1_N.R2.trim.fastq.gz
│   ├── GI1_T
│   │   ├── alignmentstats.txt
│   │   ├── circExplorer
│   │   │   ├── GI1_T.back_spliced_junction.bed
│   │   │   ├── GI1_T.back_spliced_junction.strand_fixed.bed
│   │   │   ├── GI1_T.bam
│   │   │   ├── GI1_T.bam.bai
│   │   │   ├── GI1_T.BSJ.bam
│   │   │   ├── GI1_T.BSJ.bam.csi
│   │   │   ├── GI1_T.BSJ.bed.gz
│   │   │   ├── GI1_T.BSJ.bw
│   │   │   ├── GI1_T.BSJ.foundcounts.tsv
│   │   │   ├── GI1_T.BSJ.minus.bam
│   │   │   ├── GI1_T.BSJ.minus.bam.csi
│   │   │   ├── GI1_T.BSJ.minus.bw
│   │   │   ├── GI1_T.BSJ.plus.bam
│   │   │   ├── GI1_T.BSJ.plus.bam.csi
│   │   │   ├── GI1_T.BSJ.plus.bw
│   │   │   ├── GI1_T.circExplorer.annotation_counts.tsv
│   │   │   ├── GI1_T.circExplorer.counts_table.tsv
│   │   │   ├── GI1_T_circexplorer_parse.log
│   │   │   ├── GI1_T.circRNA_known.filter1.txt
│   │   │   ├── GI1_T.circularRNA_known.txt
│   │   │   ├── GI1_T.hg38.BSJ.bam
│   │   │   ├── GI1_T.hg38.BSJ.bam.csi
│   │   │   ├── GI1_T.hg38.BSJ.bw
│   │   │   ├── GI1_T.linear.bam
│   │   │   ├── GI1_T.linear.bam.csi
│   │   │   ├── GI1_T.linear_BSJ.bam
│   │   │   ├── GI1_T.linear_BSJ.bam.csi
│   │   │   ├── GI1_T.linear_BSJ.bw
│   │   │   ├── GI1_T.linear_BSJ.hg38.bam
│   │   │   ├── GI1_T.linear_BSJ.hg38.bw
│   │   │   ├── GI1_T.linear_BSJ.NC_009333.1.bam
│   │   │   ├── GI1_T.linear_BSJ.NC_009333.1.bw
│   │   │   ├── GI1_T.linear.BSJ.readids.gz
│   │   │   ├── GI1_T.linear.bw
│   │   │   ├── GI1_T.linear.hg38.bam
│   │   │   ├── GI1_T.linear.hg38.bw
│   │   │   ├── GI1_T.linear.NC_009333.1.bam
│   │   │   ├── GI1_T.linear.NC_009333.1.bw
│   │   │   ├── GI1_T.linear_spliced.counts.tsv
│   │   │   ├── GI1_T.NC_009333.1.BSJ.bam
│   │   │   ├── GI1_T.NC_009333.1.BSJ.bam.csi
│   │   │   ├── GI1_T.NC_009333.1.BSJ.bw
│   │   │   ├── GI1_T.readcounts.tsv
│   │   │   ├── GI1_T.rid2jid.tsv.gz
│   │   │   ├── GI1_T.spliced.bam
│   │   │   ├── GI1_T.spliced.bam.csi
│   │   │   ├── GI1_T.spliced_BSJ.bam
│   │   │   ├── GI1_T.spliced_BSJ.bam.csi
│   │   │   ├── GI1_T.spliced_BSJ.bw
│   │   │   ├── GI1_T.spliced_BSJ.hg38.bam
│   │   │   ├── GI1_T.spliced_BSJ.hg38.bw
│   │   │   ├── GI1_T.spliced_BSJ.NC_009333.1.bam
│   │   │   ├── GI1_T.spliced_BSJ.NC_009333.1.bw
│   │   │   ├── GI1_T.spliced.BSJ.readids.gz
│   │   │   ├── GI1_T.spliced.bw
│   │   │   ├── GI1_T.spliced.hg38.bam
│   │   │   ├── GI1_T.spliced.hg38.bw
│   │   │   ├── GI1_T.spliced.NC_009333.1.bam
│   │   │   ├── GI1_T.spliced.NC_009333.1.bw
│   │   │   ├── low_conf_GI1_T.circRNA_known.filter1.txt
│   │   │   └── low_conf_GI1_T.circRNA_known.txt
│   │   ├── circExplorer_BWA
│   │   │   ├── GI1_T.back_spliced_junction.bed
│   │   │   ├── GI1_T.back_spliced_junction.strand_fixed.bed
│   │   │   ├── GI1_T.circExplorer_bwa.annotation_counts.tsv
│   │   │   ├── GI1_T_circexplorer_parse.log
│   │   │   ├── GI1_T.circRNA_known.filter1.txt
│   │   │   ├── GI1_T.circularRNA_known.txt
│   │   │   ├── low_conf_GI1_T.circRNA_known.filter1.txt
│   │   │   └── low_conf_GI1_T.circRNA_known.txt
│   │   ├── circRNA_finder
│   │   │   ├── GI1_T.Chimeric.out.sorted.bam
│   │   │   ├── GI1_T.Chimeric.out.sorted.bam.bai
│   │   │   ├── GI1_T.circRNA_finder.counts_table.tsv.filtered
│   │   │   ├── GI1_T.filteredJunctions.bed
│   │   │   ├── GI1_T.s_filteredJunctions.bed
│   │   │   └── GI1_T.s_filteredJunctions_fw.bed
│   │   ├── ciri
│   │   │   ├── CIRIerror.log
│   │   │   ├── GI1_T.bwa.log
│   │   │   ├── GI1_T.ciri.bam
│   │   │   ├── GI1_T.ciri.bam.csi
│   │   │   ├── GI1_T.ciri.log
│   │   │   ├── GI1_T.ciri.out
│   │   │   └── GI1_T.ciri.out.filtered
│   │   ├── CLEAR
│   │   │   ├── quant.txt
│   │   │   └── quant.txt.annotated
│   │   ├── DCC
│   │   │   ├── CircCoordinates
│   │   │   ├── CircRNACount
│   │   │   ├── CircSkipJunctions
│   │   │   ├── DCC-2023-12-20_1206.log
│   │   │   ├── GI1_T.dcc.counts_table.tsv
│   │   │   ├── GI1_T.dcc.counts_table.tsv.filtered
│   │   │   ├── LinearCount
│   │   │   ├── mate1.txt
│   │   │   ├── mate2.txt
│   │   │   └── samplesheet.txt
│   │   ├── find_circ
│   │   │   ├── GI1_T_anchors.fastq.gz
│   │   │   ├── GI1_T.unmapped.bam
│   │   │   └── GI1_T.unmapped.bam.csi
│   │   ├── GI1_T.circRNA_counts.txt.gz
│   │   ├── merge_per_sample.sh
│   │   ├── STAR1p
│   │   │   ├── GI1_T_p1.Chimeric.out.junction
│   │   │   ├── GI1_T_p1.Chimeric.out.junction.circRNA
│   │   │   ├── GI1_T_p1.Chimeric.out.junction.circRNAmapped
│   │   │   ├── GI1_T_p1.Log.final.out
│   │   │   ├── GI1_T_p1.Log.out
│   │   │   ├── GI1_T_p1.Log.progress.out
│   │   │   ├── GI1_T_p1.SJ.out.tab
│   │   │   ├── mate1
│   │   │   │   └── GI1_T_mate1.Chimeric.out.junction
│   │   │   └── mate2
│   │   │       └── GI1_T_mate2.Chimeric.out.junction
│   │   ├── STAR2p
│   │   │   ├── GI1_T_p2.bam
│   │   │   ├── GI1_T_p2.bam.csi
│   │   │   ├── GI1_T_p2.chimeric.bam
│   │   │   ├── GI1_T_p2.chimeric.bam.csi
│   │   │   ├── GI1_T_p2.Chimeric.out.junction
│   │   │   ├── GI1_T_p2.Log.final.out
│   │   │   ├── GI1_T_p2.Log.out
│   │   │   ├── GI1_T_p2.Log.progress.out
│   │   │   ├── GI1_T_p2.non_chimeric.bam
│   │   │   ├── GI1_T_p2.non_chimeric.bam.csi
│   │   │   ├── GI1_T_p2.ReadsPerGene.out.tab
│   │   │   └── GI1_T_p2.SJ.out.tab
│   │   ├── STAR_circRNAFinder
│   │   │   ├── GI1_T.Chimeric.out.junction
│   │   │   ├── GI1_T.Chimeric.out.sam
│   │   │   ├── GI1_T.Log.final.out
│   │   │   ├── GI1_T.Log.out
│   │   │   ├── GI1_T.Log.progress.out
│   │   │   └── GI1_T.SJ.out.tab
│   │   └── trim
│   │       ├── GI1_T.R1.trim.fastq.gz
│   │       └── GI1_T.R2.trim.fastq.gz
│   ├── HQ_BSJ_bams
│   │   ├── GI1_N.hg38.HQ_only.BSJ.bam
│   │   ├── GI1_N.hg38.HQ_only.BSJ.bam.bai
│   │   ├── GI1_N.HQ_only.BSJ.bam
│   │   ├── GI1_N.HQ_only.BSJ.bam.bai
│   │   ├── GI1_N.NC_009333.1.HQ_only.BSJ.bam
│   │   ├── GI1_N.NC_009333.1.HQ_only.BSJ.bam.bai
│   │   ├── GI1_T.hg38.HQ_only.BSJ.bam
│   │   ├── GI1_T.hg38.HQ_only.BSJ.bam.bai
│   │   ├── GI1_T.HQ_only.BSJ.bam
│   │   ├── GI1_T.HQ_only.BSJ.bam.bai
│   │   ├── GI1_T.NC_009333.1.HQ_only.BSJ.bam
│   │   └── GI1_T.NC_009333.1.HQ_only.BSJ.bam.bai
│   ├── pass1.out.tab
│   ├── revstranded_STAR_GeneCounts.tsv
│   ├── stranded_STAR_GeneCounts.tsv
│   └── unstranded_STAR_GeneCounts.tsv
├── runinfo.yaml
├── runslurm_snakemake_report.html
├── samples.tsv
├── snakemake.log
├── snakemake.stats
├── stats
│   ├── dryrun.231222102114.log
│   ├── runinfo.modtime.yaml
│   └── snakemake.20231222104214.log
├── submit_script.sbatch
└── tools.txt
```

#### Error tracking

Using the following command to find FAILED jobs:

```bash
grep FAIL logs/snakemake.log.jobby.short
```

The above command also gives `.err` and `.out` log files which can give further insights on reasons for failure and changes required to be made for a successful run.

#### Expected output:

The main output file is `results/circRNA_master_counts.tsv.gz`. Here are the top 3 tiles from an example output:

| Column_number | Column_title                                    | Example_1  | Example_2  | Example_3  |
| ------------- | ----------------------------------------------- | ---------- | ---------- | ---------- |
| 1             | chrom                                           | GL000220.1 | GL000220.1 | GL000220.1 |
| 2             | start                                           | 107635     | 112482     | 118578     |
| 3             | end                                             | 151634     | 156427     | 118759     |
| 4             | circExplorer_strand                             | \-1        | \-1        | \-1        |
| 5             | circExplorer_bwa_strand                         | .          | .          | .          |
| 6             | ciri_strand                                     | \-1        | \-1        | \-1        |
| 7             | dcc_strand                                      | \-1        | \-1        | \-1        |
| 8             | circrnafinder_strand                            | \-1        | \-1        | \-1        |
| 9             | flanking*sites*+                                | CC##GC     | GC##CC     | CC##GC     |
| 10            | flanking*sites*-                                | GG##GC     | GC##GG     | GG##GC     |
| 11            | sample_name                                     | GI1_N      | GI1_N      | GI1_N      |
| 12            | ntools                                          | 1          | 1          | 1          |
| 13            | HQ                                              | N          | N          | N          |
| 14            | circExplorer_read_count                         | \-1        | \-1        | \-1        |
| 15            | circExplorer_found_BSJcounts                    | \-1        | \-1        | \-1        |
| 16            | circExplorer*found_linear_BSJ*+\_counts         | \-1        | \-1        | \-1        |
| 17            | circExplorer*found_linear_spliced_BSJ*+\_counts | \-1        | \-1        | \-1        |
| 18            | circExplorer*found_linear_BSJ*-\_counts         | \-1        | \-1        | \-1        |
| 19            | circExplorer*found_linear_spliced_BSJ*-\_counts | \-1        | \-1        | \-1        |
| 20            | circExplorer*found_linear_BSJ*.\_counts         | \-1        | \-1        | \-1        |
| 21            | circExplorer*found_linear_spliced_BSJ*.\_counts | \-1        | \-1        | \-1        |
| 22            | ciri_read_count                                 | \-1        | \-1        | \-1        |
| 23            | ciri_linear_read_count                          | \-1        | \-1        | \-1        |
| 24            | circExplorer_bwa_read_count                     | 3          | 7          | 3          |
| 25            | dcc_read_count                                  | \-1        | \-1        | \-1        |
| 26            | dcc_linear_read_count                           | \-1        | \-1        | \-1        |
| 27            | circrnafinder_read_count                        | \-1        | \-1        | \-1        |
| 28            | hqcounts                                        | 1          | 1          | 1          |
| 29            | nonhqcounts                                     | 0          | 0          | 0          |
| 30            | circExplorer_annotation                         | Unknown    | Unknown    | Unknown    |
| 31            | ciri_annotation                                 | Unknown    | Unknown    | Unknown    |
| 32            | circExplorer_bwa_annotation                     | novel      | novel      | novel      |
| 33            | dcc_gene                                        | Unknown    | Unknown    | Unknown    |
| 34            | dcc_junction_type                               | Unknown    | Unknown    | Unknown    |
| 35            | dcc_annotation                                  | Unknown    | Unknown    | Unknown    |

Expected output from the sample data is stored under `.tests/expected_output`.
