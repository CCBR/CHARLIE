# CHARLIE

**C**ircrnas in **H**ost **A**nd vi**R**uses ana**L**ysis p**I**p**E**line

This is [CCBR](https://bioinformatics.ccr.cancer.gov/ccbr/)'s pipeline for DAQ (**D**etection, **A**nnotation & **Q**uantification) of cirRNAs in host and viruses using multiple circRNA detection tools and custom scripts. The input fastqs to this pipeline are generally RNAseq Fastqs.

This pipeline us built using the [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipelining framework.
This is a snakemake workflow to detect, annotate and quantify host and viral circular RNAs.

![img](https://img.shields.io/github/issues/kopardev/circRNA?style=for-the-badge)![img](https://img.shields.io/github/forks/kopardev/circRNA?style=for-the-badge)![img](https://img.shields.io/github/stars/kopardev/circRNA?style=for-the-badge)![img](https://img.shields.io/github/license/kopardev/circRNA?style=for-the-badge)

This circularRNA detection pipeline uses CIRCExplorer2 and CIRI2 in parallel to detect, quantify and annotate circRNAs. Here is a flowchart of v0.3.3:
![img](https://github.com/kopardev/circRNA/blob/master/circRNA_v0.3.3.png)

For complete documentation with tutorial go [here](https://kopardev.github.io/circRNA/)

> DISCLAIMER: New circRNA tools have been added CHARLIE and the documentation is currently out of date!

> Please contact [Vishal Koparde](mailto:vishal.koparde@nih.gov) for questions,requests,comments,etc. about CHARLIE.
