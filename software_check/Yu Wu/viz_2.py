#!/usr/bin/env python
# coding: utf-8

# In[1]:


"""
This viz_2 file is to visualize the data to determine distributions, and patterns related to the AMD stages

Section1 - import data
Section2 - show distributions plot (age_amd stage/gender_amd stage/genes expression)
Section3 - use dimension reduction method t-SNE 
"""


# In[2]:


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


# # Import data

# In[3]:


df = pd.read_csv('train_genelevel_90%.tsv', sep='\t',index_col=0)
df.head()


# In[4]:


retina = pd.read_csv('meta_retina.csv', index_col = [0])
retina.shape # 523


# In[5]:


sub_retina = pd.DataFrame(retina,columns = ['age','sex','mgs_level'])
sub_retina.head()


# In[6]:


df = sub_retina.join(df.T)
df.head()


# In[7]:


df.shape #523..


# In[8]:


fclean = df.dropna(axis=0,how='any')
fclean.shape #407


# In[9]:


fclean.head()


# # Distributions

# In[10]:


sns.set(color_codes = True)
sns.set_style("white")
sns.violinplot(x = fclean["mgs_level"], y = fclean["age"] )


# In[11]:


fig = sns.countplot(x = "mgs_level", data=fclean)
#fig.figure.savefig("count.png",dpi = 300) 
plt.show()


# In[12]:


ct = pd.crosstab(fclean.mgs_level, fclean.sex)
ct.plot.bar(stacked=True)
plt.legend(title='sex')
plt.show()


# In[13]:


dfn = pd.read_csv('train_genelevel_90%.tsv', sep='\t',index_col=0)
dfn.head()


# In[14]:


dfn['mean'] = dfn.mean(axis = 1)
dfn.head()


# In[15]:


dfn['mean'].describe()


# In[16]:


bins = [0,25,50,100,300,600,900,1200,1500,2000,4000,8000,12000,16000,21000]
cats = pd.cut(dfn['mean'], bins)


# In[17]:


def get_stats(group):
    return {'count': group.count()}

grouped = dfn['mean'].groupby(cats)
bin_counts = grouped.apply(get_stats).unstack()
#print (bin_counts)
bin_counts


# In[18]:


bin_counts.index = ['25', '50', '100', '300', '600', '900', '1200',
                    '1500', '2000', '4000','8000','12000','16000','21000']
bin_counts.index.name = 'mean'
bin_counts.plot(kind='bar', alpha=0.5, rot=0)


# # Dimension Reduction

# In[19]:


Y = TSNE(n_components = 2).fit_transform(fclean.drop(['age','sex','mgs_level'],axis=1))
Y.shape


# In[20]:


fclean['tsne-2d-one'] = Y[:,0]
fclean['tsne-2d-two'] = Y[:,1]
fclean.head()


# In[21]:


plt.figure(figsize = (10,10))
fig = sns.scatterplot(
    x = "tsne-2d-one", y ="tsne-2d-two", hue = "sex",
    palette = sns.color_palette("hls", 2),
    data = fclean
)
fig.figure.savefig("label by x and y.png",dpi = 300)


# In[22]:


plt.figure(figsize = (10,10))
fig = sns.scatterplot(
    x = "tsne-2d-one", y ="tsne-2d-two", hue = "mgs_level",
    palette = sns.color_palette("hls", 4),
    data = fclean
)
fig.figure.savefig("label mgs_level.png",dpi = 300)


# In[23]:


fclean['age'].describe()


# In[24]:


b = [50,60,70,80,90,110]
c = pd.cut(fclean['age'], b)


# In[25]:


def get_stats(group):
    return {'count': group.count()}

grouped = fclean['age'].groupby(c)
bin_counts = grouped.apply(get_stats).unstack()
#print (bin_counts)
bin_counts


# In[26]:


from math import floor, ceil


# In[27]:


fclean['AGE'] = 0
fclean.head()


# In[28]:


#fclean['AGE'] = floor(fclean['age'].values / 10)
#fclean.head()

fclean['AGE'] = fclean['age'].map(lambda x: floor(x/10) * 10)
fclean.head()


# In[29]:


plt.figure(figsize = (10,10))
fig = sns.scatterplot(
    x = "tsne-2d-one", y ="tsne-2d-two", hue = "AGE",
    palette = sns.color_palette("hls", 6),
    data = fclean
)
fig.figure.savefig("age.png",dpi = 300)


# # Dimension reduction for each stage

# In[30]:


stage1 = fclean[fclean.index.str.contains('_1')]
stage1.head()


# In[31]:


Y1 = TSNE(n_components = 2).fit_transform(stage1.drop(['tsne-2d-one','tsne-2d-two','AGE','age','sex','mgs_level'],axis=1))
Y1.shape


# In[32]:


stage1['t1_1'] = Y1[:,0]
stage1['t2_1'] = Y1[:,1]
stage1.head()


# In[33]:


plt.figure(figsize = (10,10))
fig = sns.scatterplot(
    x = "tsne-2d-one", y ="tsne-2d-two", hue = "mgs_level",
    palette = sns.color_palette("hls", 1),
    data = stage1
)
# fig.figure.savefig("label x and y.png",dpi = 300)


# In[34]:


stage4 = fclean[fclean.index.str.contains('_4')]
stage4.head()


# In[35]:


Y4 = TSNE(n_components = 2).fit_transform(stage4.drop(['tsne-2d-one','tsne-2d-two','AGE','age','sex','mgs_level'],axis=1))
Y4.shape


# In[36]:


stage4['t1_4'] = Y4[:,0]
stage4['t2_4'] = Y4[:,1]
stage4.head()


# In[37]:


plt.figure(figsize = (10,10))
fig = sns.scatterplot(
    x = "tsne-2d-one", y ="tsne-2d-two", hue = "mgs_level",
    palette = sns.color_palette("hls", 1),
    data = stage4
)
# fig.figure.savefig("label x and y.png",dpi = 300)


# In[ ]:




