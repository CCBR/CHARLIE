# CHARLIE
![img](https://img.shields.io/github/issues/CCBR/CHARLIE?style=for-the-badge)![img](https://img.shields.io/github/forks/CCBR/CHARLIE?style=for-the-badge)![img](https://img.shields.io/github/stars/CCBR/CHARLIE?style=for-the-badge)![img](https://img.shields.io/github/license/CCBR/CHARLIE?style=for-the-badge)


### Table of Contents
- [CHARLIE - **C**ircrnas in **H**ost **A**nd vi**R**uses ana**L**ysis p**I**p**E**line](#charlie)
    - [Table of Contents](#table-of-contents)
    - [1. Introduction](#1-introduction)
    - [2. Flowchart](#2-flowchart)
    - [3. Software Dependencies](#3-software-dependencies)
    - [4. Usage](#4-usages)

### 1. Introduction

 **C**ircrnas in **H**ost **A**nd vi**R**uses ana**L**ysis p**I**p**E**line

Things to know about CHARLIE:

- Snakemake workflow to detect, annotate and quantify (DAQ) host and viral circular RNAs.
- Primirarily developed to run on [BIOWULF](https://hpc.nih.gov/)
- Reach out to [Vishal Koparde](mailto:vishal.koparde@nihgov) for questions/comments/requests.


This circularRNA detection pipeline uses CIRCExplorer2, CIRI2 and many other tools in parallel to detect, quantify and annotate circRNAs. Here is a list of tools that can be run using CHARLIE:

| circRNA Detection Tool | Aligner(s) | Run by default |
| ---------------------- | ---------- | -------------- |
| [CIRCExplorer2](https://github.com/YangLab/CIRCexplorer2)          | STAR<sup>1</sup>       | Yes            |
| [CIRI2](https://sourceforge.net/projects/ciri/files/CIRI2/)                  | BWA<sup>1</sup>        | Yes            |
| [CIRCExplorer2](https://github.com/YangLab/CIRCexplorer2)          | BWA<sup>1</sup>        | Yes            |
| [CLEAR](https://github.com/YangLab/CLEAR)                  | STAR<sup>1</sup>       | Yes            |
| [DCC](https://github.com/dieterich-lab/DCC)                    | STAR<sup>2</sup>       | Yes            |
| [circRNAFinder](https://github.com/bioxfu/circRNAFinder)          | STAR<sup>3</sup>       | Yes            |
| [find_circ](https://github.com/marvin-jens/find_circ)              | Bowtie2    | Yes            |
| [MapSplice](https://github.com/merckey/MapSplice2)              | BWA<sup>2</sup>        | No             |
| [NCLScan](https://github.com/TreesLab/NCLscan)                | NovoAlign  | No             |

> Note: STAR<sup>1</sup>, STAR<sup>2</sup>, STAR<sup>3</sup> denote 3 different sets of alignment parameters, etc.

> Note: BWA<sup>1</sup>, BWA<sup>2</sup> denote 2 different alignment parameters, etc.

### 2. Flowchart
![](docs/images/CHARLIE_v0.8.x.png)

For complete documentation with tutorial go [here](https://CCBR.github.io/CHARLIE/).

> DISCLAIMER: New circRNA tools have been added CHARLIE and the documentation is currently out of date!

### 3. Software Dependencies

The following version of various bioinformatics tools are using within CHARLIE:

| tool          | version   |
| ------------- | --------- |
| blat          | 3.5       |
| bedtools      | 2.30.0    |
| bowtie        | 2-2.5.1   |
| bowtie        | 1.3.1     |
| bwa           | 0.7.17    |
| circexplorer2 | 2.3.8     |
| cufflinks     | 2.2.1     |
| cutadapt      | 1.18      |
| fastqc        | 0.11.9    |
| hisat         | 2.2.1.0   |
| java          | 1.8.0_181 |
| multiqc       | 1.9       |
| parallel      | 20231122  |
| perl          | 5.24      |
| picard        | 2.25.0    |
| python        | 2.7       |
| python        | 3.7       |
| sambamba      | 0.8.0     |
| samtools      | 1.16.1    |
| STAR          | 2.7.6a    |
| stringtie     | 2.1.4     |
| ucsc          | 450       |
| R             | 4.0.3     |
| novocraft     | 4.03.05   |


### 4. Usage

```bash
 % ./charlie


##########################################################################################

Welcome to
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
  bash /data/Ziegelbauer_lab/Pipelines/circRNA/activeDev/charlie -w/--workdir=<WORKDIR> -m/--runmode=<RUNMODE>

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
  bash /data/Ziegelbauer_lab/Pipelines/circRNA/activeDev/charlie -w=/my/ouput/folder -m=init
  bash /data/Ziegelbauer_lab/Pipelines/circRNA/activeDev/charlie -w=/my/ouput/folder -m=dryrun
  bash /data/Ziegelbauer_lab/Pipelines/circRNA/activeDev/charlie -w=/my/ouput/folder -m=run

##########################################################################################

VersionInfo:
  python          : 3.7
  snakemake       : 7.19.1
  pipeline_home   : /vf/users/Ziegelbauer_lab/Pipelines/circRNA/activeDev
  git commit/tag  : 1ae5ca091976364369784f67adffbbbf1dcdb7d5    v0.8-197-g1ae5ca0

##########################################################################################
```

> DISCLAIMER:
> 
> CHARLIE is built to be run only on [BIOWULF](https://hpc.nih.gov). A newer HPC-agnostic version of CHARLIE is planned for 2024.
