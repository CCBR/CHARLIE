# CHARLIE

![img](https://img.shields.io/github/issues/kopardev/circRNA?style=for-the-badge)![img](https://img.shields.io/github/forks/kopardev/circRNA?style=for-the-badge)![img](https://img.shields.io/github/stars/kopardev/circRNA?style=for-the-badge)![img](https://img.shields.io/github/license/kopardev/circRNA?style=for-the-badge)

### Version/Release highlights

#### v0.1.0

* base version
* PE only support

#### v0.2.x

* SE support added .. PE/SE samples handled concurrently
* `envmodules` used in Snakemake in place of `module load` statements

#### v0.3.x

* Lookup table for hg38 to hg19 circRNA annotations is updated... this eliminate one-to-many hits from the previous version
* BSJs extracted as different bam file.
* flowchart added
* adding slurmjobid to log/err file names
* v0.3.1 has significant (>10X) performance improvements at BSJ bam creation
* v0.3.3 splits BSJ bams into human and viral bams, and also converts them to bigwigs
* v0.3.4 adds hg38_rRNA_masked_plus_rRNA_plus_viruses_plus_ERCC reference (source:Sarah)

#### v0.4.x 

* [CLEAR](https://github.com/YangLab/CLEAR) added.
* wrapper script (`run_circrna_daq.sh`) added for local and cluster execution.

* "spliced reads only" bam created and split by regions

#### v0.5.x

* `run_clear` is now set to True (as default)
* `circ_quant` replaces `clear_quant` in the CLEAR rule. In order words, we are reusing the STAR alignment file and the circExplorer2 output file for running CLEAR. No need to run HISAT2 and TopHat (fusion-search with Bowtie1). This is much quicker.
* Using picard to estimate duplicates using[ *MarkDuplicates*](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-)
* Generating a per-run [multiqc](https://multiqc.info/) HTML report
* Using [*eulerr*](https://www.rdocumentation.org/packages/eulerr/versions/6.1.0) R package to generate CIRI-CircExplorer circRNA Venn diagrams and include them in the mulitqc report
* Gather per job cluster metadata like queue time, run time, job state etc. Stats are compiled in **HPC_summary** file
* CLEAR pipeline *quant.txt* file is annotated for known circRNAs
* `WORKDIR` can now be a relative path
* bam2bw conversion fix for BSJ and spliced_reads. [Issue](https://github.com/kopardev/circRNA/issues/17) closed!


#### v0.6.x

* cutadapt_min_length to cutadapt rule... setting it to 15 in config.
* customBSJs recalled from STAR alignments
  * only for PE
  * removes erroneously called CircExplorer BSJs
* create sense and anti-sense BSJ BAMs and BW for each reference (host+viruses)
* find reads which contribute to CIRI BSJs but not on the STAR list of BSJ reads, see if they contribute to novel (not called by STAR) BSJs and append novel BSJs to customBSJ list
* create linear reads BAM file
* create linear reads BigWigs for each region in the .regions file.
* Optimized pysam scripts
* fixed premature completion of singularity rules

#### v0.7.x

* 5 circRNA callers
* all-sample counts matrix with annotations

#### v0.9.x

* updates to wrapper script, many new arguments/options added
* new per-sample counts table format
* new all-sample master counts matrix with min-nreads filtering and ntools column to show number of tools supporting the circRNA call
* new version of Snakemake
* cluster_status script added for forced completion of pipeline upon TIMEOUTs
updated flowchart from lucid charts
* added circRNAfinder, find_circ, circExplorer2_bwa and other tools
* optimized execution and resource requirements
* updated viral annotations (Thanks Sara!)
* new method to extract linear counts, create linear BAMs using circExplorer2 outputs
* new job reporting using jobby and its derivatives
separated creation of BWA and BOWTIE2 index from creation of STAR index to speed things up
* parallelized find_circ
better cleanup (eg. deleting _STARgenome folders, etc.) for much smaller digital footprint
* multitude of comments throughout the snakefiles including listing of output file column descriptions
* preliminary GH actions added

#### v0.10.x

* strand are reported together, strand from all callers are reported, 
* both + and - flanking sites are reported, 
* rev-comp function updated,
* updated versions of tools to match available tools on BIOWULF.