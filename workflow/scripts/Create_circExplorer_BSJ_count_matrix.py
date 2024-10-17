#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas
import re
import numpy as np
from pathlib import Path
import os
import matplotlib.pyplot as plt
import sys
lookupfile=sys.argv[1]
hostID=sys.argv[2]
# In[27]:


def readthefile(f):
    sampleName=f.name.replace(".back_spliced_junction.bed","")
    x=pandas.read_csv(f,sep="\t",header=None)
    x.columns=["chr","start","end","name_count","score","strand"]
    x['id'] = x["chr"]+":"+x["start"].map(str)+"-"+x["end"].map(str)
    x[['name',sampleName]] = x.name_count.str.split("/",expand=True)
    x=x.loc[:,["id",sampleName]]
    x.set_index(["id"],inplace=True)
    return(x)
        


# In[2]:


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


# In[3]:
outfilename1="circExplorer_BSJ_count_matrix.txt"
outfilename="circExplorer_BSJ_count_matrix_with_annotations.txt"

files_circExplorer=list(Path(os.getcwd()).rglob("*.back_spliced_junction.bed"))
files_circExplorer=list(filter(lambda x:"_only.back" not in str(x),files_circExplorer))
files_circExplorer=list(filter(lambda x: os.stat(x).st_size !=0, files_circExplorer))
files_circExplorer.sort(key=natural_keys)
if len(files_circExplorer)==0:
    for f in [outfilename1,outfilename]:
        if os.path.exists(f):
            os.remove(f)
        os.mknod(f)
    exit()


# In[35]:


circE_count_matrix=readthefile(files_circExplorer[0])
print(circE_count_matrix.head())


# In[36]:


for j in range(1,len(files_circExplorer)):
    x=readthefile(files_circExplorer[j])
    circE_count_matrix=pandas.concat([circE_count_matrix,x],axis=1,join="outer",sort=False)
circE_count_matrix=circE_count_matrix.sort_index()
print(circE_count_matrix.head())


# In[37]:


circE_count_matrix.fillna(0,inplace=True)
circE_count_matrix.head()
circE_count_matrix.to_csv(outfilename1,sep="\t",header=True)


# In[38]:


annotations=pandas.read_csv(lookupfile,sep="\t",header=0)
annotations.set_index([hostID],inplace=True)
annotations.head()


# In[11]:


annotations.shape


# In[39]:


x=circE_count_matrix.join(annotations)
x.to_csv(outfilename,sep="\t",header=True)


# In[14]:


print(circE_count_matrix.shape)
print(x.shape)


