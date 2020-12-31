f=open("hg19_hg38_annotated_lookup.txt")
hg38_2_hg19=dict()
for l in f.readlines():
	l=l.strip().split("\t")
	hg19ID=l[0]
	hg38ID=l[1]
	strand=l[2]
	circRNA_ID=l[3]
	genomic_length=l[4]
	spliced_seq_length=l[5]
	samples=l[6].split(",")	
	repeats=l[7]
	annotation=l[8].split(",")	
	best_transcript=l[9]
	gene_symbol=l[10]
	circRNA_study=l[11].split(",")
	if not hg38ID in hg38_2_hg19:
		hg38_2_hg19[hg38ID]=dict()
		hg38_2_hg19[hg38ID]['hg19ID']=list()
		hg38_2_hg19[hg38ID]['circRNA_ID']=list()
		hg38_2_hg19[hg38ID]['samples']=list()
		hg38_2_hg19[hg38ID]['annotation']=list()
		hg38_2_hg19[hg38ID]['circRNA_study']=list()
	hg38_2_hg19[hg38ID]['hg19ID'].append(hg19ID)
	hg38_2_hg19[hg38ID]['strand']=strand
	hg38_2_hg19[hg38ID]['circRNA_ID'].append(circRNA_ID)
	hg38_2_hg19[hg38ID]['genomic_length']=genomic_length
	hg38_2_hg19[hg38ID]['spliced_seq_length']=spliced_seq_length
	hg38_2_hg19[hg38ID]['samples'].extend(samples)
	hg38_2_hg19[hg38ID]['repeats']=repeats
	hg38_2_hg19[hg38ID]['annotation'].extend(annotation)
	hg38_2_hg19[hg38ID]['best_transcript']=best_transcript
	hg38_2_hg19[hg38ID]['gene_symbol']=gene_symbol
	hg38_2_hg19[hg38ID]['circRNA_study'].extend(circRNA_study)

#print("\t".join(["hg38ID","hg19ID","strand","circRNA.ID","genomic.length","spliced.seq.length","samples","repeats","annotation","best.transcript","gene.symbol","circRNA.study"]),)
for k,v in hg38_2_hg19.items():
	l=list()
	l.append(k)
	l.append(",".join(set(v['hg19ID'])))
	l.append(v['strand'])
	l.append(",".join(set(v['circRNA_ID'])))
	l.append(v['genomic_length'])
	l.append(v['spliced_seq_length'])
	l.append(",".join(set(v['samples'])))
	l.append(v['repeats'])
	l.append(",".join(set(v['annotation'])))
	l.append(v['best_transcript'])
	l.append(v['gene_symbol'])
	l.append(",".join(set(v['circRNA_study'])))
	print("\t".join(l),)


