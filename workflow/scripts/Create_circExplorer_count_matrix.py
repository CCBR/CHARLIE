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

outfilename1="circExplorer_count_matrix.txt"
outfilename="circExplorer_count_matrix_with_annotations.txt"

files_circExplorer=list(Path(os.getcwd()).rglob("*.circularRNA_known.txt"))
files_circExplorer=list(filter(lambda x: False if str(x).find("low_conf")!=-1 else True, files_circExplorer))
files_circExplorer=list(filter(lambda x: os.stat(x).st_size !=0, files_circExplorer))
files_circExplorer.sort(key=natural_keys)
print(files_circExplorer)
if len(files_circExplorer)==0:
	for f in [outfilename1,outfilename]:
		if os.path.exists(f):
			os.remove(f)
		os.mknod(f)
	exit()

# In[12]:


f=files_circExplorer[0]
sampleName=f.name.replace(".circularRNA_known.txt","")
print("Reading file:",f)
print("Sample Name:",sampleName)
x=pandas.read_csv(f,sep="\t",header=None,usecols=[0,1,2,12])
x[hostID]=x[0].astype(str)+":"+x[1].astype(str)+"-"+x[2].astype(str)
x[sampleName+"_circE"]=x[12].astype(str)
x.drop([0,1,2,12],inplace=True,axis=1)
x.set_index([hostID],inplace=True)
circE_count_matrix=x


# In[8]:



print(circE_count_matrix.head(),circE_count_matrix.tail())
print(circE_count_matrix.shape)


# In[13]:

for i in range(1,len(files_circExplorer)):
	f=files_circExplorer[i]
	print("Currently reading file:"+str(f))
	x=pandas.read_csv(f,sep="\t",header=None,usecols=[0,1,2,12])
	print("Head of this file looks like this:")
	print(x.head())
	sampleName=f.name.replace(".circularRNA_known.txt","")
	# x=pandas.read_csv(f,sep="\t",header=None,usecols=[0,1,2,12])
	print("SampleName is:"+sampleName)
	x[hostID]=x[0].astype(str)+":"+x[1].astype(str)+"-"+x[2].astype(str)
	x[sampleName+"_circE"]=x[12].astype(str)
	print(x.head())
	x.drop([0,1,2,12],inplace=True,axis=1)
	x.set_index([hostID],inplace=True)
	print(x.head())
	print("Before concat")
	print(circE_count_matrix.head())


# In[14]:


	circE_count_matrix = circE_count_matrix.loc[~circE_count_matrix.index.duplicated(keep='first')]
	x = x.loc[~x.index.duplicated(keep='first')]
	circE_count_matrix=pandas.concat([circE_count_matrix,x],axis=1,join="outer",sort=False)
	print("After concat")
	print(circE_count_matrix.head())


# In[9]:


circE_count_matrix.fillna(0,inplace=True)
print(circE_count_matrix.head())
circE_count_matrix.to_csv(outfilename1,sep="\t",header=True)


# In[10]:


annotations=pandas.read_csv(lookupfile,sep="\t",header=0)
annotations.set_index([hostID],inplace=True)
annotations.head()


# In[11]:


annotations.shape


# In[12]:


x=circE_count_matrix.join(annotations)
x.to_csv(outfilename,sep="\t",header=True)


# In[14]:


print(circE_count_matrix.shape)
print(x.shape)

