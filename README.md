# circRNA

![img](https://img.shields.io/github/issues/kopardev/circRNA?style=for-the-badge)![img](https://img.shields.io/github/forks/kopardev/circRNA?style=for-the-badge)![img](https://img.shields.io/github/stars/kopardev/circRNA?style=for-the-badge)![img](https://img.shields.io/github/license/kopardev/circRNA?style=for-the-badge)

This circularRNA detection pipeline uses CIRCExplorer2 and CIRI2 in parallel to detect, quantify and annotate circRNAs. Here is a flowchart of v0.3.3:
![img](https://github.com/kopardev/circRNA/blob/master/circRNA_v0.3.3.png)

# Version/Release highlights
## v0.1.0
* base version
* PE only support
## v0.2.x
* SE support added .. PE/SE samples handled concurrently
* `envmodules` used in Snakemake in place of `module load` statements
## v0.3.x
* Lookup table for hg38 to hg19 circRNA annotations is updated... this eliminate one-to-many hits from the previous version
* BSJs extracted as different bam file.
* flowchart added
* adding slurmjobid to log/err file names
* v0.3.1 has significant (>10X) performance improvements at BSJ bam creation
* v0.3.3 splits BSJ bams into human and viral bams, and also converts them to bigwigs
* v0.3.4 adds hg38_rRNA_masked_plus_rRNA_plus_viruses_plus_ERCC reference (source:Sarah)
