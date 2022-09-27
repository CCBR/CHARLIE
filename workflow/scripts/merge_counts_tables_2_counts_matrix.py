#!/usr/bin/env python
# coding: utf-8


import pandas
import argparse
import re
from pathlib import Path
import os

debug=False

# no truncations during if debug: print pandas data frames
pandas.set_option('display.max_rows', None)
pandas.set_option('display.max_columns', None)
pandas.set_option('display.width', None)
pandas.set_option('display.max_colwidth', None)


parser = argparse.ArgumentParser(description='Merge per sample counts tables to a single annotated counts matrix')
parser.add_argument('--results_folder', dest='folder', type=str, required=True,
                    help='folder whose subfolders contain all the per sample counts tables')
parser.add_argument('--lookup_table', dest='lookup', type=str, required=True,
                    help='annotation lookup table')
parser.add_argument('-o',dest='outfile',required=True,help='merged countsmatrix')
args = parser.parse_args()

def prefix_test(colname):
    # returns true if the col needs to be an int
    if colname.endswith("_circExplorer_read_count"): return True
    if colname.endswith("_ciri_read_count"): return True
    if colname.endswith("_ntools"): return True
    if colname.endswith(".length"): return True
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


per_sample_files=list(Path(args.folder).rglob("circRNA_counts.txt"))
per_sample_files=list(filter(lambda x: os.stat(x).st_size !=0, per_sample_files))
per_sample_files.sort(key=natural_keys)

if len(per_sample_files)==0:
    if os.path.exists(args.outfile):
        os.remove(args.outfile)
    exit()

f=per_sample_files[0]
if debug: print("Currently reading file:"+str(f))
count_matrix=pandas.read_csv(f,sep="\t",header=0)
count_matrix['circRNA_id2']=count_matrix['circRNA_id'].astype(str)+"##"+count_matrix['strand'].astype(str)
count_matrix.drop(['circRNA_id','strand'],axis=1,inplace=True)
count_matrix.set_index(['circRNA_id2'],inplace=True)
if debug: print("Head of this file looks like this:")
if debug: print(count_matrix.head())
for i in range(1,len(per_sample_files)):
    f=per_sample_files[i]
    if debug: print("Currently reading file:"+str(f))
    x=pandas.read_csv(f,sep="\t",header=0)
    x['circRNA_id2']=x['circRNA_id'].astype(str)+"##"+x['strand'].astype(str)
    x.set_index(['circRNA_id2'],inplace=True)
    x.drop(['circRNA_id','strand'],axis=1,inplace=True)
    if debug: print("Head of this file looks like this:")
    if debug: print(x.head())
    count_matrix=pandas.concat([count_matrix,x],axis=1,join="outer",sort=False)
    count_matrix.fillna(0,inplace=True)

if debug: print(count_matrix.head())


annotations=pandas.read_csv(args.lookup,sep="\t",header=0)
annotations_cols=annotations.columns
# annotations.set_index([annotations_cols[0]],inplace=True)
if debug: print(annotations.head())
annotations['circRNA_id2']=annotations[annotations_cols[0]].astype(str)+"##"+annotations['strand'].astype(str)
annotations.set_index(annotations['circRNA_id2'],inplace=True)
annotations.drop(['strand'],axis=1,inplace=True)
if debug: print(annotations.head())
if debug: print(annotations.shape)

# count_matrix=pandas.concat([count_matrix,annotations],axis=1,join="outer",sort=False)
count_matrix = pandas.merge(count_matrix,annotations,left_index=True,right_index=True,sort=False,how='left')
count_matrix['circRNA_id2']=count_matrix.index
count_matrix.fillna(0,inplace=True)
count_matrix.replace('\.','0', regex=True,inplace=True)
coltypes=dict()
for col in count_matrix.columns:
    # if debug: print('column', col,':', type(col[0]))
    # if debug: print('column', col,': prefix test result: ', prefix_test(col))
    coltypes[col]=str
    if prefix_test(col):
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

count_matrix.to_csv(args.outfile,sep="\t",header=True,index=False)


