## References

The reference sequences comprises of the host genome and the viral genomes.

### Fasta

_hg38_ and _mm39_ genome builds are chosen to represent hosts. Ribosomal sequences (_45S, 5S_) are downloaded from NCBI. _hg38_ and _mm39_ were masked for rRNA sequence and _45S_ and _5S_ sequences from NCBI are appended as separate chromosomes. The following viral sequences were appended to the rRNA masked _hg38_ reference:

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
