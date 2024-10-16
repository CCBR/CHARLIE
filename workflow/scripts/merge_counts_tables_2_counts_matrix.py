#!/usr/bin/env python
# coding: utf-8


import pandas
import argparse
import re
from pathlib import Path
import os
import numpy

debug=False

# no truncations during if debug: print pandas data frames
pandas.set_option('display.max_rows', None)
pandas.set_option('display.max_columns', None)
pandas.set_option('display.width', None)
pandas.set_option('display.max_colwidth', None)


parser = argparse.ArgumentParser(description='Merge per sample counts tables to a single annotated counts matrix')
parser.add_argument('--per_sample_tables', nargs='+', dest='ctables', type=argparse.FileType('r'), required=True,
                    help='space separated list of input per-sample count tables')
parser.add_argument('--lookup_table', dest='lookup', type=argparse.FileType('r'), required=True,
                    help='annotation lookup table (host-only)')
parser.add_argument('-o',dest='outfile',required=True,type=argparse.FileType('w'),help='merged countsmatrix')
args = parser.parse_args()

if debug:
    print(args)

def prefix_counts(colname):
    # returns true if the col needs to be an int
    if colname.endswith("_read_count"): return True
    if colname.endswith("_ntools"): return True
    if colname.endswith(".length"): return True
    return False

def prefix_annotations(colname):
    if colname.endswith("_annotation"): return True
    return False


def atof(text):
    try:
        retval = float(text)
    except ValueError:
        retval = text
    return retval

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    float regex comes from https://stackoverflow.com/a/12643073/190597
    '''
    return [ atof(c) for c in re.split(r'[+-]?([0-9]+(?:[.][0-9]*)?|[.][0-9]+)', str(text)) ]

def get_count_and_annotation_columns(df):
    count_cols=['circRNA_id2']
    annotation_cols=['circRNA_id2']
    for col in df.columns:
        if prefix_counts(col):
            count_cols.append(col)
        if prefix_annotations(col):
            annotation_cols.append(col)
    return count_cols,annotation_cols

def readin_counts_file(f):
    intable=pandas.read_csv(f,sep="\t",header=0)
    intable['circRNA_id2']=intable['circRNA_id'].astype(str)+"##"+intable['strand'].astype(str)
    intable.drop(['circRNA_id','strand'],axis=1,inplace=True)
    count_cols,annotation_cols=get_count_and_annotation_columns(intable)
    count_table=intable[count_cols]
    annotation_table=intable[annotation_cols]
    count_table.set_index(['circRNA_id2'],inplace=True)
    annotation_table.set_index(['circRNA_id2'],inplace=True)
    return(count_table,annotation_table)


# per_sample_files=list(Path(args.folder).rglob("*.circRNA_counts.txt"))
# per_sample_files=list(filter(lambda x: os.stat(x).st_size !=0, per_sample_files))
# per_sample_files.sort(key=natural_keys)

per_sample_files=args.ctables

if debug: print(per_sample_files)

annotation_tables=list()
f=per_sample_files[0]
if debug: print("Currently reading file:"+str(f))
ctable,atable=readin_counts_file(f)
annotation_tables.append(atable)
count_matrix=ctable.copy()
# count_matrix.set_index(['circRNA_id2'],inplace=True)
if debug: print("Head of this file looks like this:")
if debug: print(count_matrix.head())
for i in range(1,len(per_sample_files)):
    f=per_sample_files[i]
    if debug: print("Currently reading file:"+str(f))
    ctable,atable=readin_counts_file(f)
    # ctable.set_index(['circRNA_id2'],inplace=True)
    if debug: print("Head of this file looks like this:")
    if debug: print(ctable.head())
    count_matrix=pandas.concat([count_matrix,ctable],axis=1,join="outer",sort=False)
    count_matrix.fillna(0,inplace=True)
    annotation_tables.append(atable)

for i,a in enumerate(annotation_tables):
    if i==0:
        amatrix=a.copy()
    else:
        oldi=set(list(amatrix.index))
        newi=set(list(a.index))
        toadd=newi-oldi
        suba=a.loc[list(toadd)]
        amatrix=pandas.concat([amatrix,suba])


if debug: print(count_matrix.head())
if debug: print(annotation_tables[0].head())
if debug: print(count_matrix.shape)
if debug: print(annotation_tables[0].shape)
if debug: print(annotation_tables[1].shape)
if debug: print(amatrix.shape)


annotations=pandas.read_csv(args.lookup,sep="\t",header=0)
annotations_cols=annotations.columns
# annotations.set_index([annotations_cols[0]],inplace=True)
annotations['circRNA_id2']=annotations[annotations_cols[0]].astype(str)+"##"+annotations['strand'].astype(str)
annotations.set_index(annotations['circRNA_id2'],inplace=True)
# annotations.drop(['strand'],axis=1,inplace=True)
if debug: print(annotations.head())
if debug: print(annotations.shape)



# count_matrix=pandas.concat([count_matrix,annotations],axis=1,join="outer",sort=False)
cmatrix = pandas.merge(amatrix,annotations,left_index=True,right_index=True,sort=False,how='left')
cmatrix['circRNA_id2']=cmatrix.index
count_matrix = pandas.merge(count_matrix,cmatrix,left_index=True,right_index=True,sort=False,how='left')
count_matrix.replace('.',numpy.nan,inplace=True)
count_matrix.fillna(0,inplace=True)
# count_matrix.replace(re.compile('\.'),'0', regex=True,inplace=True)

if debug: print(count_matrix.head())
if debug: print(count_matrix.shape)

coltypes=dict()
for col in count_matrix.columns:
    coltypes[col]=str
    if prefix_counts(col):
        # count_matrix[[col]].replace(re.compile('\.'),'0', regex=True,inplace=True)
        coltypes[col]=int
count_matrix = count_matrix.astype(coltypes)
count_matrix[['circRNA_coord', 'circRNA_strand']] = count_matrix['circRNA_id2'].str.split('##', expand=True)
count_matrix.drop(['circRNA_id2'],axis=1,inplace=True)
cols=list(count_matrix.columns)
col1index=cols.index('circRNA_coord')
col2index=cols.index('circRNA_strand')
other_indices=list(set(range(len(cols)))-set([col1index,col2index]))
new_order=['circRNA_coord','circRNA_strand']
for i in other_indices:
    new_order.append(cols[i])
count_matrix=count_matrix[new_order]
if debug: print(count_matrix.head())

df2 = count_matrix[list(filter(lambda x:x.endswith("_read_count"),list(count_matrix.columns)))]
df2 = df2.astype('int')
count_matrix['sum_of_all_counts'] = df2.sum(axis=1)
df3 = count_matrix[list(filter(lambda x:x.endswith("_ntools"),list(count_matrix.columns)))]
df3 = df3.astype('int')
count_matrix['sum_of_all_ntools'] = df3.sum(axis=1)

count_matrix = count_matrix.sort_values(by=['sum_of_all_ntools','sum_of_all_counts'], ascending=False)
count_matrix.drop(['sum_of_all_ntools','sum_of_all_counts'],axis=1,inplace=True)
count_matrix.to_csv(args.outfile,sep="\t",header=True,index=False)


