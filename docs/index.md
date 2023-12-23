## CHARLIE

![img](https://img.shields.io/github/issues/kopardev/circRNA?style=for-the-badge)![img](https://img.shields.io/github/forks/kopardev/circRNA?style=for-the-badge)![img](https://img.shields.io/github/stars/kopardev/circRNA?style=for-the-badge)![img](https://img.shields.io/github/license/kopardev/circRNA?style=for-the-badge)

CHARLIE=**C**ircrnas in **H**ost **A**nd vi**R**uses ana**L**ysis p**I**p**E**line

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


<!-- ### References

<sup>1</sup>Zhang XO*, Dong R*, Zhang Y, Zhang JL, Luo Z, Zhang J, Chen LL, Yang L. Diverse alternative back-splicing and alternative splicing landscape of circular RNAs. *Genome Res*, 2016, 26:1277-1287 doi:10.1101/gr.202895.115

<sup>2</sup>Yuan Gao, Jinyang Zhang and Fangqing Zhao. Circular RNA identification based on multiple seed matching. Briefings in Bioinformatics (2017) doi: 10.1093/bib/bbx014.

<sup>3</sup>Ma XK, Wang MR, Liu CX, Dong R, Carmichael GG, Chen LL and Yang L. A CLEAR pipeline for direct comparison of circular and linear RNA expression. 2019, bioRxiv doi: 10.1101/668657 -->

