## circRNA DAQ  Pipeline

![img](https://img.shields.io/github/issues/kopardev/circRNA?style=for-the-badge)![img](https://img.shields.io/github/forks/kopardev/circRNA?style=for-the-badge)![img](https://img.shields.io/github/stars/kopardev/circRNA?style=for-the-badge)![img](https://img.shields.io/github/license/kopardev/circRNA?style=for-the-badge)

The reference sequences comprises of the host genome and the viral genomes.

### Fasta

*hg38* genome build is chosen to represent human (host). Human ribosomal sequences (*45S, 5S*) are downloaded from NCBI. *hg38* was masked for rRNA sequence and *45S* and *5S* sequences from NCBI are appended as separate chromosomes. The following viral sequences were appended to the rRNA masked *hg38* reference:

| RefSeq  Sequence | RefSeq assembly accession | Notes                                                 |
| ---------------- | ------------------------- | ----------------------------------------------------- |
| NC_007605.1      | GCF_002402265.1           | Human gammaherpesvirus 4 (Epstein-Barr virus)         |
| NC_000898.1      | GCF_000846365.1           | Human betaherpesvirus 6B                              |
| NC_001664.4      | GCF_000845685.2           | Human betaherpesvirus 6A                              |
| NC_001716.2      | GCF_000848125.1           | Human betaherpesvirus 7                               |
| NC_006273.2      | GCF_000845245.1           | Human betaherpesvirus 5                               |
| NC_009333.1      | GCF_000838265.1           | Human gammaherpesvirus 8                              |
| NC_045512.2      | GCF_009858895.2           | Severe acute respiratory syndrome-related coronavirus |
| MN485971.1       | xx                        | HIV from Belgium ... GTF is hand curated              |

Location: The entire resource bundle is available at `/data/Ziegelbauer_lab/resources/hg38_rRNA_masked_plus_rRNA_plus_viruses_plus_ERCC` on biowulf. This location also have additional bash scritpts required for aggregating annotations and building indices required by different aligners.

<u>**Update** **(02/10/21)**</u>

The following viral sequence has also been appended to the reference:

| RefSeq  Sequence | RefSeq assembly accession | Notes                                                        |
| ---------------- | ------------------------- | ------------------------------------------------------------ |
| NC_001806.2      | GCF_000859985.2           | [Human alphaherpesvirus 1 (Herpes simplex virus type 1)](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=10298&lvl=3&lin=f&keep=1&srchmode=1&unlock) (strain 17) |

Location: The updated resource bundle is at `/data/Ziegelbauer_lab/resources/hg38_rRNA_masked_plus_rRNA_plus_viruses_plus_ERCC.v2` on biowulf

### Annotations

Gencode release 36 is used to annotate the human transcripts. Customized RefSeq annotations are used for annotating the viral sequences. 

#### Viruses

To properly annotate viral sequences the following process was followed:

**Downloading annotations from NCBI:**

For example, to download gtf annotations for KSHV genome search "KSHV" at [NCBI](https://www.ncbi.nlm.nih.gov/assembly/)

![image-20200807133058218](https://tva1.sinaimg.cn/large/007S8ZIlgy1ghirj71jv5j30qn01taa5.jpg)

On the results page (assembly page), click on "FTP directory for RefSeq assembly"

![image-20200807133320457](https://tva1.sinaimg.cn/large/007S8ZIlgy1ghirln4fc7j30s60sawj6.jpg) 

The FTP page has the relevant sequence and annotations files:

![image-20200807133421026](https://tva1.sinaimg.cn/large/007S8ZIlgy1ghirmr2a4dj30t30ajq66.jpg)

Download the file ending with *.gtf.gz*:

```bash
> curl -L "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/838/265/GCF_000838265.1_ViralProj14158/GCF_000838265.1_ViralProj14158_genomic.gtf.gz" > GCF_000838265.1_ViralProj14158_genomic.gtf.gz

> gzip -d GCF_000838265.1_ViralProj14158_genomic.gtf.gz
```

**What needs to be fixed?**

The GTF file downloaded from NBCI (RefSeq annotations) needs fixing in order to work seemlessly with STAR and CCBR_Pipeliner:

1. GTF is missing lines for "transcript" in column 3. Transcript coordinates can be different from "gene" coordinates. I am going to pick all non-"gene" lines for a particular gene, filter them for a particular "transcript", ~~find the leftmost and rightmost coordinate and use them as coordinates for the new "transcript" line in the new GTF file~~. <font color='blue'>Using leftmost and rightmost script eliminates the UTR regions and hence if a gene has only one transcript (most common scenario), the transcript coordinates mimic the gene coordinates.</font>
   eg. for **HHV8GK18_gp01**  gene these are the lines in the original file

   ```bash
   NC_009333.1	RefSeq	gene	4	1083	.	+	.	gene_id "HHV8GK18_gp01"; db_xref "GeneID:4961511"; gbkey "Gene"; gene "K1"; gene_biotype "protein_coding"; locus_tag "HHV8GK18_gp01";
   NC_009333.1	RefSeq	exon	29	1083	.	+	.	gene_id "HHV8GK18_gp01"; transcript_id "HHV8GK18_gp01"; gbkey "mRNA"; gene "K1"; locus_tag "HHV8GK18_gp01"; exon_number "1";
   NC_009333.1	RefSeq	CDS	105	941	.	+	0	gene_id "HHV8GK18_gp01"; transcript_id "HHV8GK18_gp01"; gbkey "CDS"; gene "K1"; locus_tag "HHV8GK18_gp01"; note "contains an Ig domain"; product "K1"; protein_id "YP_001129350.1"; exon_number "1";
   NC_009333.1	RefSeq	start_codon	105	107	.	+	0	gene_id "HHV8GK18_gp01"; transcript_id "HHV8GK18_gp01"; gbkey "CDS"; gene "K1"; locus_tag "HHV8GK18_gp01"; note "contains an Ig domain"; product "K1"; protein_id "YP_001129350.1"; exon_number "1";
   NC_009333.1	RefSeq	stop_codon	942	944	.	+	0	gene_id "HHV8GK18_gp01"; transcript_id "HHV8GK18_gp01"; gbkey "CDS"; gene "K1"; locus_tag "HHV8GK18_gp01"; note "contains an Ig domain"; product "K1"; protein_id "YP_001129350.1"; exon_number "1";
   ```

2. Some transcripts have ids as "unknown_transcript_1". New name will be of the format "**gene_id**+_transcript_1". All child features for the transcript in question need  the **transcript_id** to be changed to the new name.

   ```bash
   NC_009333.1	RefSeq	gene	3179	17026	.	+	.	gene_id "HHV8GK18_gp03"; db_xref "GeneID:4961521"; gbkey "Gene"; gene "ORF6"; gene_biotype "protein_coding"; locus_tag "HHV8GK18_gp03";
   NC_009333.1	RefSeq	CDS	3179	6574	.	+	0	gene_id "HHV8GK18_gp03"; transcript_id "unknown_transcript_1"; gbkey "CDS"; gene "ORF6"; locus_tag "HHV8GK18_gp03"; note "herpesvirus core gene UL29 family"; product "ORF6"; protein_id "YP_001129352.1";
   NC_009333.1	RefSeq	start_codon	3179	3181	.	+	0	gene_id "HHV8GK18_gp03"; transcript_id "unknown_transcript_1"; gbkey "CDS"; gene "ORF6"; locus_tag "HHV8GK18_gp03"; note "herpesvirus core gene UL29 family"; product "ORF6"; protein_id "YP_001129352.1";
   NC_009333.1	RefSeq	stop_codon	6575	6577	.	+	0	gene_id "HHV8GK18_gp03"; transcript_id "unknown_transcript_1"; gbkey "CDS"; gene "ORF6"; locus_tag "HHV8GK18_gp03"; note "herpesvirus core gene UL29 family"; product "ORF6"; protein_id "YP_001129352.1";
   ```

3. Some **gene_id** 's are empty. I am planning on reporting these in a separate file ending with *.extralines*. These need to be edited manually and appended to the output GTF file. I will be using "chromosome_name:start_coordinate-end_coordinate" as the format for generating a unique **gene_id** for replacement GTF entries. Eg.

   ```bash
   NC_009333.1	RefSeq	exon	118075	118097	.	-	.	gene_id ""; transcript_id "unknown_transcript_1"; gbkey "misc_RNA"; product "miR-K10"; exon_number "1";
   NC_009333.1	RefSeq	exon	127997	129368	.	+	.	gene_id ""; transcript_id "unknown_transcript_1"; gbkey "mRNA"; note "K14, ORF74; polyA_site not reported"; exon_number "1";
   NC_009333.1	RefSeq	exon	129517	130671	.	+	.	gene_id ""; transcript_id "unknown_transcript_1"; gbkey "mRNA"; note "K14, ORF74; polyA_site not reported"; exon_number "2";
   ```

   will be changed to something like this:

   ```bash
   NC_009333.1     RefSeq  gene    118075  118097  .       -       .       gene_id "NC_009333.1:118075-118097"; gene_name "NC_009333.1:118075-118097"; gene_biotype "miRNA"; gbkey "misc_RNA"; product "miR-K10";
   NC_009333.1     RefSeq  transcript      118075  118097  .       -       .       gene_id "NC_009333.1:118075-118097"; gene_name "NC_009333.1:118075-118097"; gene_biotype "miRNA"; transcript_id "NC_009333.1:118075-118097_
   transcript_1"; transcript_name "NC_009333.1:118075-118097_transcript_1"; transcript_type "miRNA"; gbkey "misc_RNA"; product "miR-K10";
   NC_009333.1     RefSeq  exon    118075  118097  .       -       .       gene_id "NC_009333.1:118075-118097"; transcript_id "NC_009333.1:118075-118097_transcript_1"; gbkey "misc_RNA"; product "miR-K10"; exon_number "1";
   NC_009333.1     RefSeq  gene    127997  130671  .       +       .       gene_id "NC_009333.1:127997_130671"; gene_name "NC_009333.1:127997_130671"; gene_biotype "mRNA"; gbkey "mRNA"; note "K14, ORF74; polyA_site not reporte
   d";
   NC_009333.1     RefSeq  transcript    127997  130671  .       +       .       gene_id "NC_009333.1:127997_130671"; gene_name "NC_009333.1:127997_130671"; gene_biotype "mRNA"; transcript_id "NC_009333.1:127997_130671_transcr
   ipt_1"; transcript_name "NC_009333.1:127997_130671_transcript_1"; transcript_type "NC_009333.1:127997_130671_transcript_1";gbkey "mRNA"; note "K14, ORF74; polyA_site not reported";
   NC_009333.1     RefSeq  exon    127997  129368  .       +       .       gene_id "NC_009333.1:127997_130671"; transcript_id "NC_009333.1:127997_130671_transcript_1"; gbkey "mRNA"; note "K14, ORF74; polyA_site not reported";
   exon_number "1";
   NC_009333.1     RefSeq  exon    129517  130671  .       +       .       gene_id "NC_009333.1:127997_130671"; transcript_id "NC_009333.1:127997_130671_transcript_1"; gbkey "mRNA"; note "K14, ORF74; polyA_site not reported";
   exon_number "2";
   ```

4. Only few genes have **exon** features. Most of them have **CDS** only. The **exon** line is added to have coordinates same as the **CDS**, if it is missing.

   eg. This

   ```bash
   NC_009333.1	RefSeq	gene	15756	17026	.	+	.	gene_id "HHV8GK18_gp08"; gene_name "ORF11"; gene_biotype "protein_coding"; db_xref "GeneID:4961439"; gbkey "Gene"; locus_tag "HHV8GK18_gp08";
   NC_009333.1	RefSeq	transcript	15756	16979	.	+	.	gene_id "HHV8GK18_gp08"; gene_name "ORF11"; gene_biotype "protein_coding"; transcript_id "HHV8GK18_gp08_transcript_1"; transcript_name "HHV8GK18_gp08_transcript_1"; transcript_type "protein_coding";
   NC_009333.1	RefSeq	CDS	15756	16976	.	+	0	gene_id "HHV8GK18_gp08"; transcript_id "HHV8GK18_gp08_transcript_1"; gbkey "CDS"; gene "ORF11"; locus_tag "HHV8GK18_gp08"; note "derived from herpesvirus dUTPase; related to HHV-5 UL82, UL83 and UL84"; product "ORF11"; protein_id "YP_001129357.1";
   NC_009333.1	RefSeq	start_codon	15756	15758	.	+	0	gene_id "HHV8GK18_gp08"; transcript_id "HHV8GK18_gp08_transcript_1"; gbkey "CDS"; gene "ORF11"; locus_tag "HHV8GK18_gp08"; note "derived from herpesvirus dUTPase; related to HHV-5 UL82, UL83 and UL84"; product "ORF11"; protein_id "YP_001129357.1";
   NC_009333.1	RefSeq	stop_codon	16977	16979	.	+	0	gene_id "HHV8GK18_gp08"; transcript_id "HHV8GK18_gp08_transcript_1"; gbkey "CDS"; gene "ORF11"; locus_tag "HHV8GK18_gp08"; note "derived from herpesvirus dUTPase; related to HHV-5 UL82, UL83 and UL84"; product "ORF11"; protein_id "YP_001129357.1";
   ```

   becomes

   ```bash
   NC_009333.1	RefSeq	gene	15756	17026	.	+	.	gene_id "HHV8GK18_gp08"; gene_name "ORF11"; gene_biotype "protein_coding"; db_xref "GeneID:4961439"; gbkey "Gene"; locus_tag "HHV8GK18_gp08";
   NC_009333.1	RefSeq	transcript	15756	16979	.	+	.	gene_id "HHV8GK18_gp08"; gene_name "ORF11"; gene_biotype "protein_coding"; transcript_id "HHV8GK18_gp08_transcript_1"; transcript_name "HHV8GK18_gp08_transcript_1"; transcript_type "protein_coding";
   NC_009333.1	RefSeq	exon	15756	16976	.	+	.	gene_id "HHV8GK18_gp08"; transcript_id "HHV8GK18_gp08_transcript_1"; gbkey "CDS"; gene "ORF11"; locus_tag "HHV8GK18_gp08"; note "derived from herpesvirus dUTPase; related to HHV-5 UL82, UL83 and UL84"; product "ORF11"; protein_id "YP_001129357.1";
   NC_009333.1	RefSeq	CDS	15756	16976	.	+	0	gene_id "HHV8GK18_gp08"; transcript_id "HHV8GK18_gp08_transcript_1"; gbkey "CDS"; gene "ORF11"; locus_tag "HHV8GK18_gp08"; note "derived from herpesvirus dUTPase; related to HHV-5 UL82, UL83 and UL84"; product "ORF11"; protein_id "YP_001129357.1";
   NC_009333.1	RefSeq	start_codon	15756	15758	.	+	0	gene_id "HHV8GK18_gp08"; transcript_id "HHV8GK18_gp08_transcript_1"; gbkey "CDS"; gene "ORF11"; locus_tag "HHV8GK18_gp08"; note "derived from herpesvirus dUTPase; related to HHV-5 UL82, UL83 and UL84"; product "ORF11"; protein_id "YP_001129357.1";
   NC_009333.1	RefSeq	stop_codon	16977	16979	.	+	0	gene_id "HHV8GK18_gp08"; transcript_id "HHV8GK18_gp08_transcript_1"; gbkey "CDS"; gene "ORF11"; locus_tag "HHV8GK18_gp08"; note "derived from herpesvirus dUTPase; related to HHV-5 UL82, UL83 and UL84"; product "ORF11"; protein_id "YP_001129357.1";
   ```

**Python script to automate:**

Simple python script take in RefSeq's GTF and spits out new GTF. This script (*fix_refseq_gtf.py*) is provided in the scripts folder of the repo.:

```python
# Author: Vishal N. Koparde
# CCBR NCI
# Date: Aug, 2020


import sys,copy,argparse


parser = argparse.ArgumentParser()
parser.add_argument('-i',dest='ingtf', required=True, type=str, help="Input RefSeq GTF ..downloaded from NCBI ftp server")
parser.add_argument('-o',dest='outgtf', required=True, type=str, help="Modified Output RefSeq GTF")
args = parser.parse_args()

def get_gene_id(column9):
    x=column9.strip().split()
    for i,value in enumerate(x):
        if value=="gene_id":
            gene_id_index=i+1
            break
    gene_id=x[gene_id_index]
    return gene_id

def get_gene_biotype(column9):
    x=column9.strip().split()
    found=0
    for i,value in enumerate(x):
        if value=="gene_type" or value=="gene_biotype":
            gene_biotype_index=i+1
            found=1
            break
    if found==0:
        return '"unknown";'
    gene_biotype=x[gene_biotype_index]
    return gene_biotype

def get_gene_name(column9):
    x=column9.strip().split()
    found=0
    for i,value in enumerate(x):
        if value=="gene" or value=="gene_name":
            gene_index=i+1
            found=1
            break
    if found==0:
        return ""
    gene_name=x[gene_index]
    return gene_name

def get_transcript_id(column9):
    x=column9.strip().split()
    found=0
    for i,value in enumerate(x):
        if value=="transcript_id":
            transcript_id_index=i+1
            found=1
            break
    if found==0:
        return '"transcript_id_unknown";'
    transcript_id=x[transcript_id_index]
    return transcript_id

def fix_transcript_id(column9,g):
    x=column9.strip().split()
    found=0
    for i,value in enumerate(x):
        if value=="transcript_id":
            transcript_id_index=i+1
            found=1
            break
    x[transcript_id_index]=g
    if found==0:
        x.append("transcript_id")
        x.append(g)
    x=" ".join(x)
    return x   

def create_new_transript_id(g,i):
    n=g.split('"')
    n[-2]+="_transcript_"+str(i)
    n='"'.join(n)
    return n

def are_exons_present(transcript_lines):
    for l in transcript_lines:
        l_split=l.strip().split("\t")
        if l_split[2]=="exon":
            return True
    else:
        return False

#create genelist
genelist=[]
gene_coords=dict()
all_gtflines=list(filter(lambda x:not x.startswith("#"),open(args.ingtf).readlines()))
blank_gene_id_lines=[]
for f in all_gtflines:
    its_a_gene=0
    if f.strip().split("\t")[2]=="gene":
        its_a_gene=1
    gene_id=get_gene_id(f.strip().split("\t")[8])
    if gene_id=='"";':
        blank_gene_id_lines.append(f)
        continue
    genelist.append(gene_id)
    if its_a_gene==1 and not gene_id in gene_coords:
        gene_coords[gene_id]=(int(f.strip().split("\t")[3]),int(f.strip().split("\t")[4]))
genelist=list(set(genelist))
# print(genelist)
# print(len(blank_gene_id_lines))

#get genes2transcripts ... this is only for verifying that every gene has only 1 transript... this is the assumption
gene_id_2_transcript_ids=dict()
for g in genelist:
    if not g in gene_id_2_transcript_ids:
        gene_id_2_transcript_ids[g]=list()
    lines_with_gene_id=list(filter(lambda x: g in x,all_gtflines))
    non_gene_lines=list(filter(lambda x:x.split("\t")[2]!="gene",lines_with_gene_id))
    for l in non_gene_lines:
        t_id=get_transcript_id(l.strip().split("\t")[8])
        if t_id!='"transcript_id_unknown";':
            gene_id_2_transcript_ids[g].append(t_id)
            gene_id_2_transcript_ids[g]=list(set(gene_id_2_transcript_ids[g]))

geneid2transcriptidfile=open(args.ingtf+".geneid2transcriptid",'w')
for k,v in gene_id_2_transcript_ids.items():
    geneid2transcriptidfile.write("%s\t%s\n"%(k,v))
geneid2transcriptidfile.close()

#get genenames
gene_id_2_gene_name=dict()
for g in genelist:
    if not g in gene_id_2_gene_name:
        gene_id_2_gene_name[g]=list()
    lines_with_gene_id=list(filter(lambda x: g in x,all_gtflines))
    gene_line=list(filter(lambda x:x.split("\t")[2]=="gene",lines_with_gene_id))
    # if len(gene_line)==0:
    #     for l in lines_with_gene_id:
    #         print(l,)
    gene_line=gene_line[0]
    gene_name=get_gene_name(gene_line.split("\t")[8])
    if gene_name=="":
        gene_name=g
    gene_id_2_gene_name[g]=gene_name
# for k,v in gene_id_2_gene_name.items():
#     print(k,v)
    
#get transcript coordinates
gene_id_2_transcript_coordinates=dict()
for g in genelist:
    # print("gene=",g)
    if not g in gene_id_2_transcript_coordinates:
        gene_id_2_transcript_coordinates[g]=list()
    if len(gene_id_2_transcript_ids[g])==1:
        gene_id_2_transcript_coordinates[g].append(gene_coords[g])
    else:
        lines_with_gene_id=list(filter(lambda x: g in x,all_gtflines))
        non_gene_lines=list(filter(lambda x:x.split("\t")[2]!="gene",lines_with_gene_id))
        for t in gene_id_2_transcript_ids[g]:
            # print("transcript=",t)
            transcript_lines=list(filter(lambda x:t in x,non_gene_lines))
            coords=[]
            for l in transcript_lines:
                # print(l.strip())
                l_split=l.split("\t")
                coords.append(int(l_split[3]))
                coords.append(int(l_split[4]))
            # print()
            gene_id_2_transcript_coordinates[g].append((min(coords),max(coords)))
    # print(gene_id_2_transcript_coordinates[g])
# for k,v in gene_id_2_transcript_coordinates.items():
    # print(k,v)
# exit()

#get gene biotype\
gene_id_2_gene_biotype=dict()
for g in genelist:
    lines_with_gene_id=list(filter(lambda x: g in x,all_gtflines))
    gene_line=list(filter(lambda x:x.split("\t")[2]=="gene",lines_with_gene_id))
    gene_line=gene_line[0]
    gene_biotype=get_gene_biotype(gene_line.split("\t")[8])
    gene_id_2_gene_biotype[g]=gene_biotype
# for k,v in gene_id_2_gene_biotype.items():
#     print(k,v)

out=open(args.outgtf,'w')    
for g in genelist:
    lines_with_gene_id=list(filter(lambda x: g in x,all_gtflines))
    gene_line=list(filter(lambda x:x.split("\t")[2]=="gene",lines_with_gene_id))
    gene_line=gene_line[0]
    gene_line=gene_line.split("\t")
    others=gene_line.pop(-1)
    gene_line_copy=copy.copy(gene_line)
    # other key value pairs to add in the gene_line(col9)
    others_to_add=[]
    # print("others=",others)
    for o in others.strip().split("; "):
        # print("o=",o)
        o2=o.split(" ")
        # print("o2=",o2)
        key=o2[0]
        value=o2[1:]
        value=" ".join(value)
        # print("key=",key)
        # print("value=",value)
        if key in ["gene_id","gene","gene_name","gene_type","gene_biotype"]:
            continue
        else:
            others_to_add.append(key)
            if not ";" in value:
                others_to_add.append(value+";")
            else:
                others_to_add.append(value)

    col9=[]
    col9.append("gene_id")
    col9.append(g)
    col9.append("gene_name")
    col9.append(gene_id_2_gene_name[g])
    col9.append("gene_biotype")
    col9.append(gene_id_2_gene_biotype[g])
    col9plus=copy.copy(col9)
    col9plus.extend(others_to_add)
    gene_col9=" ".join(col9plus)
    gene_line.append(gene_col9)
    gene_line="\t".join(gene_line)
    out.write("%s\n"%(gene_line))

    non_gene_lines=list(filter(lambda x:x.split("\t")[2]!="gene",lines_with_gene_id))
    for i,t in enumerate(gene_id_2_transcript_ids[g]):
        transcript_line=copy.copy(gene_line_copy)
        transcript_line[2]="transcript"
        transcript_line[3]=str(gene_id_2_transcript_coordinates[g][i][0])
        transcript_line[4]=str(gene_id_2_transcript_coordinates[g][i][1])
        new_trascript_id=create_new_transript_id(g,i+1)
        transcript_col9=copy.copy(col9)
        transcript_col9.append("transcript_id")
        transcript_col9.append(new_trascript_id)
        transcript_col9.append("transcript_name")
        transcript_col9.append(new_trascript_id)
        transcript_col9.append("transcript_type")
        transcript_col9.append(gene_id_2_gene_biotype[g])
        transcript_col9=" ".join(transcript_col9)
        transcript_line.append(transcript_col9)
        transcript_line="\t".join(transcript_line)   
        out.write("%s\n"%(transcript_line))

        transcript_lines=list(filter(lambda x:t in x,non_gene_lines))
        have_exons=are_exons_present(transcript_lines)
        for l in transcript_lines:
            # print(l)
            l=l.strip().split("\t")
            tofix=l.pop(-1)
            l.append(fix_transcript_id(tofix,new_trascript_id))
            if l[2]=="CDS" and have_exons==False:
                l2=copy.copy(l)
                l2[7]="."
                l2[2]="exon"
                l2="\t".join(l2)
                out.write("%s\n"%(l2))
            l="\t".join(l)
            out.write("%s\n"%(l))
            # print(l)
out.close()

out=open(args.ingtf+".extralines",'w')
for b in blank_gene_id_lines:
    out.write(b)
out.close()
```

Final versions of the GTF as located here:

`/data/Ziegelbauer_lab/circRNADetection/viral_db/GTF/final_versions`