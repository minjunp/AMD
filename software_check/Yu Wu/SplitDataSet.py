#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import manifold
from sklearn.manifold import TSNE
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import NullFormatter
from sklearn.decomposition import PCA

from time import time


# In[25]:


"""
By doing this file, we can finally seperate original dataset into test and train dataset 

"""


# In[2]:


df = pd.read_csv('genelevel_rsem_expectedcounts_byrid_nooutliers.counts.matrix.tsv', sep='\t',index_col=0)
df.head() #(58051, 453)
#df.shape 


# In[3]:


anno = pd.read_csv('AnnotationFile.tsv', sep = '\t', index_col = [0])
anno.head()


# In[4]:


#Removing X&Y chromosome
anno = anno[anno.chromosome_name != 'X']
anno.head()
anno = anno[anno.chromosome_name != 'Y']
anno.head()


# In[5]:


both = df.join(anno)#
both= both.dropna(axis=0,how='any')
both.shape # (55169, 460)
both.head()


# In[6]:


cleandf = both.drop(['external_gene_name','chromosome_name','start_position','end_position','strand','gene_length','gene_biotype'],axis = 1)
cleandf.head()


# In[7]:


ndf = cleandf.T
ndf['MGS_LEVEL'] = 0
ndf.head()


# In[8]:


new = pd.DataFrame()
new = ndf.reset_index()
new.head()


# In[9]:


mlist = new ['index'].apply(lambda x:x[-1]).tolist()
ndf['MGS_LEVEL'] = mlist
level = ndf.pop('MGS_LEVEL')
ndf.insert(0,'MGS_LEVEL',level)
ndf.head()


# In[10]:


from sklearn.model_selection import train_test_split

#split dataset into train and test dataset
X, y = ndf.iloc[:,1:], ndf.iloc[:,0]
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=0)

#shape of train and test
X_train.shape, y_train.shape # ((407, 55169), (407,))
X_test.shape, y_test.shape #((46, 55169), (46,))


# In[11]:


X_test.head()


# In[12]:


one = X_train[X_train.index.str.contains('_4')]
one.shape # 97 / 160 / 98 / 52


# In[13]:


stage4 = X_test[X_test.index.str.contains('_4')]
#stage1- 8
#stage2 - 15 
#stage3 -14
#stage4 - 9 


# In[14]:


result = X_train.T
result.head()


# In[22]:


result.to_csv("./train_genelevel_90%.tsv", sep='\t')


# In[23]:


result2= X_test.T
result2.shape


# In[24]:


result2.to_csv("./test_genelevel_10%.tsv", sep='\t')


# In[ ]:




