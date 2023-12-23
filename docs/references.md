## CHARLIE

![img](https://img.shields.io/github/issues/kopardev/circRNA?style=for-the-badge)![img](https://img.shields.io/github/forks/kopardev/circRNA?style=for-the-badge)![img](https://img.shields.io/github/stars/kopardev/circRNA?style=for-the-badge)![img](https://img.shields.io/github/license/kopardev/circRNA?style=for-the-badge)

The reference sequences comprises of the host genome and the viral genomes.

### Fasta

*hg38* and *mm39* genome builds are chosen to represent hosts. Ribosomal sequences (*45S, 5S*) are downloaded from NCBI. *hg38* and *mm39* were masked for rRNA sequence and *45S* and *5S* sequences from NCBI are appended as separate chromosomes. The following viral sequences were appended to the rRNA masked *hg38* reference:

```bash
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
```

Location: The entire resource bundle is available at `/data/CCBR_Pipeliner/db/PipeDB/charlie/fastas_gtfs` on BIOWULF. This location also have additional bash scritpts required for aggregating annotations and building indices required by different aligners. 

When `-m=dryrun` is run for the first time after initialization (`-m=init`), the appropriate host+additives+viruses fasta and gtf files are created on the fly, which are then used to build aligner reference indexes automatically. 
