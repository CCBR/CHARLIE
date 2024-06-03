# CHARLIE development version

# CHARLIE 0.10.1

- strand are reported together, strand from all callers are reported,
- both + and - flanking sites are reported,
- rev-comp function updated,
- updated versions of tools to match available tools on BIOWULF.

# CHARLIE 0.9.0

Significant upgrades since the last release:

- updates to wrapper script, many new arguments/options added
- new per-sample counts table format
- new all-sample master counts matrix with min-nreads filtering and ntools column to show number of tools supporting the circRNA call
- new version of Snakemake
- cluster_status script added for forced completion of pipeline upon TIMEOUTs
- updated flowchart from lucid charts
- added circRNAfinder, find_circ, circExplorer2_bwa and other tools
- optimized execution and resource requirements
- updated viral annotations (Thanks Sara!)
- new method to extract linear counts, create linear BAMs using circExplorer2 outputs
- new job reporting using jobby and its derivatives
- separated creation of BWA and BOWTIE2 index from creation of STAR index to speed things up
- parallelized find_circ
- better cleanup (eg. deleting _STARgenome folders, etc.) for much smaller digital footprint
- multitude of comments throughout the snakefiles including listing of output file column descriptions
- preliminary GH actions added

# CHARLIE 0.7.0

- 5 circRNA callers
- all-sample counts matrix with annotations

# CHARLIE 0.6.9

- Optimized pysam scripts
- fixed premature completion of singularity rules

# CHARLIE 0.6.5

- updated config.yaml to use the latest HSV-1 annotations received from Sarah (050421)

# CHARLIE 0.6.4

- create linear reads BAM file
- create linear reads BigWigs for each region in the .regions file.

# CHARLIE 0.6.3

- QOS not working for Taka... removed from cluster.json
- recall rule requires python/3.7 ... env module updated

# CHARLIE 0.6.2

- BSJ files are in BSJ subfolder... bug fix for v0.6.1

# CHARLIE 0.6.1

- customBSJs recalled from STAR alignments
    - only for PE
    - removes erroneously called CircExplorer BSJs
- create sense and anti-sense BSJ BAMs and BW for each reference (host+viruses)
- find reads which contribute to CIRI BSJs but not on the STAR list of BSJ reads, see if they contribute to novel (not called by STAR) BSJs and append novel BSJs to customBSJ list

# CHARLIE 0.6.0

cutadapt_min_length to cutadapt rule... setting it to 15 in config (for miRNAs, Biot and short viral features)
