# CHARLIE
![img](https://img.shields.io/github/issues/CCBR/CHARLIE?style=for-the-badge)![img](https://img.shields.io/github/forks/CCBR/CHARLIE?style=for-the-badge)![img](https://img.shields.io/github/stars/CCBR/CHARLIE?style=for-the-badge)![img](https://img.shields.io/github/license/CCBR/CHARLIE?style=for-the-badge)

### **C**ircrnas in **H**ost **A**nd vi**R**uses ana**L**ysis p**I**p**E**line


Things to know about CHARLIE:

- Snakemake workflow to detect, annotate and quantify (DAQ) host and viral circular RNAs.
- Primirarily developed to run on [BIOWULF](https://hpc.nih.gov/)
- Reach out to [Vishal Koparde](mailto:vishal.koparde@nihgov) for questions/comments/requests.


This circularRNA detection pipeline uses CIRCExplorer2, CIRI2 and many other tools in parallel to detect, quantify and annotate circRNAs. Here is a flowchart of v0.3.3:
![img](https://github.com/CCBR/CHARLIE/blob/master/circRNA_v0.3.3.png)

For complete documentation with tutorial go [here](https://ccbr.github.io/CCBR_circRNA_DAQ/)


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
  * KT899744      [HSV-1 strain KOS]
  * MH636806.1    [MHV68 (Murine herpesvirus 68 strain WUMS)]

##########################################################################################

USAGE:
  bash ./charlie -w/--workdir=<WORKDIR> -m/--runmode=<RUNMODE>

Required Arguments:
1.  WORKDIR     : [Type: String]: Absolute or relative path to the output folder with write permissions.

2.  RUNMODE     : [Type: String] Valid options:
    * init      : initialize workdir
    * dryrun    : dry run snakemake to generate DAG
    * run       : run with slurm
    * runlocal  : run without submitting to sbatch
    * help      : print this help
    ADVANCED RUNMODES (use with caution!!)
    * unlock    : unlock WORKDIR if locked by snakemake NEVER UNLOCK WORKDIR WHERE PIPELINE IS CURRENTLY RUNNING!
    * reconfig  : recreate config file in WORKDIR (debugging option) EDITS TO config.yaml WILL BE LOST!
    * reset     : DELETE workdir dir and re-init it (debugging option) EDITS TO ALL FILES IN WORKDIR WILL BE LOST!
    * printbinds: print singularity binds (paths)
    * local     : same as runlocal

Example commands:
  bash ./charlie -w=/my/ouput/folder -m=init
  bash ./charlie -w=/my/ouput/folder -m=dryrun
  bash ./charlie -w=/my/ouput/folder -m=run

##########################################################################################

VersionInfo:
  python          : 3.7
  snakemake       : 7.19.1
  pipeline_home   : /blah/blah/blah
  git commit/tag  : b1f7695dffba79a245900257b34ebd2d0d91abc7    v0.8-148-gb1f7695

##########################################################################################
```