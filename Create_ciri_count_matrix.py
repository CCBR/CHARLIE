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
#get_ipython().run_line_magic('matplotlib', 'inline')

lookupfile=sys.argv[1]
hostID=sys.argv[2]
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


#files_circExplorer=list(Path(os.getcwd()).rglob("*_human_only.circularRNA_known.txt"))
files_ciri=list(Path(os.getcwd()).rglob("*.ciri.out"))
#filter out files in the "old" folder
#files_circExplorer=list(filter(lambda x: not re.search('/old/', str(x)), files_circExplorer))
files_ciri=list(filter(lambda x: not re.search('/old/', str(x)), files_ciri))
#files_circExplorer.sort(key=natural_keys)
files_ciri.sort(key=natural_keys)


# In[4]:


f=files_ciri[0]
sampleName=f.name.replace(".ciri.out","")
x=pandas.read_csv(f,sep="\t",header=0,usecols=["chr","circRNA_start","circRNA_end","#junction_reads"])
print(x.head())
x["circRNA_start"]=x["circRNA_start"].astype(int)-1
x[hostID]=x["chr"].astype(str)+":"+x["circRNA_start"].astype(str)+"-"+x["circRNA_end"].astype(str)
x[sampleName+"_ciri"]=x["#junction_reads"].astype(str)
x.drop(["chr","circRNA_start","circRNA_end","#junction_reads"],inplace=True,axis=1)
x.set_index([hostID],inplace=True)
ciri_count_matrix=x
print(ciri_count_matrix.head())


# In[5]:


for f in files_ciri[1:]:
    sampleName=f.name.replace(".ciri.out","")
    print(f,sampleName)
    x=pandas.read_csv(f,sep="\t",header=0,usecols=["chr","circRNA_start","circRNA_end","#junction_reads"])
    x["circRNA_start"]=x["circRNA_start"].astype(int)-1
    x[hostID]=x["chr"].astype(str)+":"+x["circRNA_start"].astype(str)+"-"+x["circRNA_end"].astype(str)
    x[sampleName+"_ciri"]=x["#junction_reads"].astype(str)
    x.drop(["chr","circRNA_start","circRNA_end","#junction_reads"],inplace=True,axis=1)
    x.set_index([hostID],inplace=True)
    ciri_count_matrix=pandas.concat([ciri_count_matrix,x],axis=1,join="outer",sort=False)
ciri_count_matrix.head()


# In[6]:


ciri_count_matrix.fillna(0,inplace=True)
ciri_count_matrix.head()
ciri_count_matrix.to_csv("ciri_count_matrix.txt",sep="\t",header=True)


# In[7]:


annotations=pandas.read_csv(lookupfile,sep="\t",header=0)
annotations.set_index([hostID],inplace=True)
annotations.head()


# In[8]:


x=ciri_count_matrix.join(annotations)
x.to_csv("ciri_count_matrix_with_annotations.txt",sep="\t",header=True)


# In[ ]:




