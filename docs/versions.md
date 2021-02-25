# circRNA DAQ Pipeline

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