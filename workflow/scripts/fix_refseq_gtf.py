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
