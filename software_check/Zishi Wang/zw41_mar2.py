#!/usr/bin/env python
# coding: utf-8

# # Data Visualization & Analysis

# ## Load Data

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import umap.umap_ as umap # install "umap-learn" instead of "umap"
#from IPython.display import display


# In[2]:


# Read raw data
original = pd.read_csv('Data/genelevel_rsem_expectedcounts_byrid_nooutliers.counts.matrix.tsv', sep = '\t', index_col = [0])
original = original.transpose()
original.index.name = 'r_id'
original.reset_index(inplace = True)

#display(original.head())
#display(original.shape)


# In[3]:


# Read batch-corrected data
cleaned = pd.read_csv('Data/batchcorrected_log2cpm.tsv', sep = '\t', index_col = [0])
cleaned = cleaned.transpose()
cleaned.index.name = 'r_id'
cleaned.reset_index(inplace = True)

#display(cleaned.head())
#display(cleaned.shape)


# In[4]:


# Read trainning set of data
train = pd.read_csv('Data/train_genelevel_90%.tsv', sep = '\t', index_col = [0])
train = train.transpose()
train.index.name = 'r_id'
train.reset_index(inplace = True)

#display(train.head())
#display(train.shape)


# ## Relation with Batches

# In[5]:


# Merge gene data in training set with batches to take into account
meta_retina = pd.read_csv('Data/meta_retina.csv', encoding = 'mac_roman')
batches = pd.DataFrame(meta_retina, columns = ["rna_isolation_batch", "library_prepper", "sex", "age"])
df = batches.join(train, how = 'inner')
df = df.dropna(axis = 0)

#display(df.head())
#display(df.shape)


# In[6]:


# Functions that shows the relation between expression and batch for a specific gene
def map_to_age(df, gene_id):
    
    df_stage = []
    for i in range(4):
        df_stage.append(df[df.r_id.str.contains('_' + str(i + 1))])
    labels = ['No AMD / normal', 'Early stage AMD', 'Intermediate AMD', 'Advanced AMD']
    fig, ax = plt.subplots(figsize = (16, 12), dpi = 80, facecolor = 'w', edgecolor = 'k')
    for i in range(4):
        ax.scatter(df_stage[i]['age'], df_stage[i][gene_id], label = labels[i])
    ax.legend()
    plt.title('Gene Expression v.s. Age (id: ' + gene_id +')', fontsize = 18)
    plt.xlabel('Age', fontsize = 18)
    plt.ylabel('The Degree of Gene Expression', fontsize = 18)
    plt.show()

def map_to_gender(df, gene_id):
    
    df_stage = []
    for i in range(4):
        df_stage.append(df[df.r_id.str.contains('_' + str(i + 1))])
    gender1 = ['F', 'M']
    gender2 = ['Female', 'Male']
    labels = ['No AMD / normal', 'Early stage AMD', 'Intermediate AMD', 'Advanced AMD']
    colors = ['blue', 'orange', 'green', 'red']
    tick = round((max(df[gene_id]) - min(df[gene_id]))) / 5
    ticks = []
    for i in range(-2, 5):
        ticks.append(round(min(df[gene_id])) + i * tick)
    fig, ax = plt.subplots(nrows = 2, ncols = 4, figsize = (16, 6))
    for i in range(2):
        for j in range(4):
            ax[i][j].scatter(df_stage[j][df_stage[j].sex == gender1[i]]['r_id'], 
                             df_stage[j][df_stage[j].sex == gender1[i]][gene_id], c = colors[j])
            ax[i][j].get_xaxis().set_visible(False)
            ax[i][j].set_yticks(ticks)
            #ax[i][j].set_ylabel('The Degree of Gene Expression')
            ax[i][j].set_title(gender2[i] + ", " + labels[j])
    fig.suptitle('Gene Expression v.s. Gender (id: ' + gene_id + ')', fontsize = 18)
    plt.show()

def map_to_lib_prep(df, gene_id):
    
    df_stage = []
    for i in range(4):
        df_stage.append(df[df.r_id.str.contains('_' + str(i + 1))])
    lib_prep = ['MRS', 'RRP']
    labels = ['No AMD / normal', 'Early stage AMD', 'Intermediate AMD', 'Advanced AMD']
    colors = ['blue', 'orange', 'green', 'red']
    tick = round((max(df[gene_id]) - min(df[gene_id]))) / 5
    ticks = []
    for i in range(-2, 5):
        ticks.append(round(min(df[gene_id])) + i * tick)
    fig, ax = plt.subplots(nrows = 2, ncols = 4, figsize = (16, 6))
    for i in range(2):
        for j in range(4):
            ax[i][j].scatter(df_stage[j][df_stage[j].library_prepper == lib_prep[i]]['r_id'], 
                             df_stage[j][df_stage[j].library_prepper == lib_prep[i]][gene_id], c = colors[j])
            ax[i][j].get_xaxis().set_visible(False)
            ax[i][j].set_yticks(ticks)
            #ax[i][j].set_ylabel('The Degree of Gene Expression')
            ax[i][j].set_title(lib_prep[i] + ", " + labels[j])
    fig.suptitle('Gene Expression v.s. Library Prepper (id: ' + gene_id + ')', fontsize = 18)
    plt.show()
    
def map_to_ria(df, gene_id):
    
    df_stage = []
    for i in range(4):
        df_stage.append(df[df.r_id.str.contains('_' + str(i + 1))])
    ria = ['isobatch1', 'isobatch2']
    labels = ['No AMD / normal', 'Early stage AMD', 'Intermediate AMD', 'Advanced AMD']
    colors = ['blue', 'orange', 'green', 'red']
    tick = round((max(df[gene_id]) - min(df[gene_id]))) / 5
    ticks = []
    for i in range(-2, 5):
        ticks.append(round(min(df[gene_id])) + i * tick)
    fig, ax = plt.subplots(nrows = 2, ncols = 4, figsize = (16, 6))
    for i in range(2):
        for j in range(4):
            ax[i][j].scatter(df_stage[j][df_stage[j].rna_isolation_batch == ria[i]]['r_id'], 
                             df_stage[j][df_stage[j].rna_isolation_batch == ria[i]][gene_id], c = colors[j])
            ax[i][j].get_xaxis().set_visible(False)
            ax[i][j].set_yticks(ticks)
            #ax[i][j].set_ylabel('The Degree of Gene Expression')
            ax[i][j].set_title(ria[i] + ", " + labels[j])
    fig.suptitle('Gene Expression v.s. RNA Isolation Batch (id: ' + gene_id + ')', fontsize = 18)
    plt.show()

def map_to_all(df, gene_id):
    
    map_to_age(df, gene_id)
    map_to_gender(df, gene_id)
    map_to_lib_prep(df, gene_id)
    map_to_ria(df, gene_id)


# In[7]:


# Usage example of methods above
gene_id = 'ENSG00000000419'
map_to_all(df, gene_id)


# ## Dimentional Reduction

# In[8]:


# Default parameters of umap.UMAP: 
#     n_components = 2 (The dimension of the space to embed into)
#     n_neighbors = 15 (The size of local neighborhood, global v.s. local)
#     min_dist = 0.1 (The effective minimum distance between embedded points, degree of tightness)
#     metric = 'euclidean' (The metric to use to compute distances in high dimensional space)
params = [2, 30, 0.2, 'euclidean']


# In[9]:


# Divide original data into 4 group according to AMD stages
original_stage1 = original[original.r_id.str.contains('_1')]
original_stage2 = original[original.r_id.str.contains('_2')]
original_stage3 = original[original.r_id.str.contains('_3')]
original_stage4 = original[original.r_id.str.contains('_4')]
X = np.vstack([original_stage1, original_stage2, original_stage3, original_stage4])

label_stage1 = np.zeros(len(original_stage1))
label_stage2 = np.ones(len(original_stage2))
label_stage3 = np.full(len(original_stage3), 2)
label_stage4 = np.full(len(original_stage4), 3)
y = np.concatenate([label_stage1, label_stage2, label_stage3, label_stage4])


# In[10]:


# Ideal case where correct classification (target array) is provided
X_embedded = umap.UMAP(
    n_components = params[0],
    n_neighbors = params[1],
    min_dist = params[2],
    metric = params[3]
).fit_transform(X, y)

fig, ax = plt.subplots(figsize = (8, 6), dpi = 80, facecolor = 'w', edgecolor = 'k')
pc1 = ax.scatter(X_embedded[y == 0, 0], X_embedded[y == 0, 1], label = 'No AMD / normal')
pc2 = ax.scatter(X_embedded[y == 1, 0], X_embedded[y == 1, 1], label = 'Early stage AMD')
pc3 = ax.scatter(X_embedded[y == 2, 0], X_embedded[y == 2, 1], label = 'Intermediate AMD')
pc4 = ax.scatter(X_embedded[y == 3, 0], X_embedded[y == 3, 1], label = 'Advanced AMD')
ax.set_title("Original Data with Target Stage")
ax.legend()
plt.show()


# In[11]:


# Actual case where no classification is provided
X_embedded = umap.UMAP(
    n_components = params[0],
    n_neighbors = params[1],
    min_dist = params[2],
    metric = params[3]
).fit_transform(X)

fig, ax = plt.subplots(figsize = (8, 6), dpi = 80, facecolor = 'w', edgecolor = 'k')
pc1 = ax.scatter(X_embedded[y == 0, 0], X_embedded[y == 0, 1], label = 'No AMD / normal')
pc2 = ax.scatter(X_embedded[y == 1, 0], X_embedded[y == 1, 1], label = 'Early stage AMD')
pc3 = ax.scatter(X_embedded[y == 2, 0], X_embedded[y == 2, 1], label = 'Intermediate AMD')
pc4 = ax.scatter(X_embedded[y == 3, 0], X_embedded[y == 3, 1], label = 'Advanced AMD')
ax.set_title("Original Data without Target Stage")
ax.legend()
plt.show()


# In[12]:


# Divide cleaned data into 4 group according to AMD stages
cleaned_stage1 = cleaned[cleaned.r_id.str.contains('_1')]
cleaned_stage2 = cleaned[cleaned.r_id.str.contains('_2')]
cleaned_stage3 = cleaned[cleaned.r_id.str.contains('_3')]
cleaned_stage4 = cleaned[cleaned.r_id.str.contains('_4')]
X = np.vstack([cleaned_stage1, cleaned_stage2, cleaned_stage3, cleaned_stage4])

label_stage1 = np.zeros(len(cleaned_stage1))
label_stage2 = np.ones(len(cleaned_stage2))
label_stage3 = np.full(len(cleaned_stage3), 2)
label_stage4 = np.full(len(cleaned_stage4), 3)
y = np.concatenate([label_stage1, label_stage2, label_stage3, label_stage4])


# In[13]:


# Ideal case where correct classification (target array) is provided
X_embedded = umap.UMAP(
    n_components = params[0],
    n_neighbors = params[1],
    min_dist = params[2],
    metric = params[3]
).fit_transform(X, y)

fig, ax = plt.subplots(figsize = (8, 6), dpi = 80, facecolor = 'w', edgecolor = 'k')
pc1 = ax.scatter(X_embedded[y == 0, 0], X_embedded[y == 0, 1], label = 'No AMD / normal')
pc2 = ax.scatter(X_embedded[y == 1, 0], X_embedded[y == 1, 1], label = 'Early stage AMD')
pc3 = ax.scatter(X_embedded[y == 2, 0], X_embedded[y == 2, 1], label = 'Intermediate AMD')
pc4 = ax.scatter(X_embedded[y == 3, 0], X_embedded[y == 3, 1], label = 'Advanced AMD')
ax.set_title("Batch-corrected Data with Target Stage")
ax.legend()
plt.show()


# In[14]:


# Actual case where no classification is provided
X_embedded = umap.UMAP(
    n_components = params[0],
    n_neighbors = params[1],
    min_dist = params[2],
    metric = params[3]
).fit_transform(X)

fig, ax = plt.subplots(figsize = (8, 6), dpi = 80, facecolor = 'w', edgecolor = 'k')
pc1 = ax.scatter(X_embedded[y == 0, 0], X_embedded[y == 0, 1], label = 'No AMD / normal')
pc2 = ax.scatter(X_embedded[y == 1, 0], X_embedded[y == 1, 1], label = 'Early stage AMD')
pc3 = ax.scatter(X_embedded[y == 2, 0], X_embedded[y == 2, 1], label = 'Intermediate AMD')
pc4 = ax.scatter(X_embedded[y == 3, 0], X_embedded[y == 3, 1], label = 'Advanced AMD')
ax.set_title("Batch-corrected Data without Target Stage")
ax.legend()
plt.show()


# ## Dimentional Reduction for Data in Training Set

# In[15]:


# Divide data in training set into 4 group according to AMD stages
train_stage1 = train[train.r_id.str.contains('_1')]
train_stage2 = train[train.r_id.str.contains('_2')]
train_stage3 = train[train.r_id.str.contains('_3')]
train_stage4 = train[train.r_id.str.contains('_4')]
X = np.vstack([train_stage1, train_stage2, train_stage3, train_stage4])

label_stage1 = np.zeros(len(train_stage1))
label_stage2 = np.ones(len(train_stage2))
label_stage3 = np.full(len(train_stage3), 2)
label_stage4 = np.full(len(train_stage4), 3)
y = np.concatenate([label_stage1, label_stage2, label_stage3, label_stage4])


# In[16]:


# Ideal case where correct classification (target array) is provided
X_embedded = umap.UMAP(
    n_components = params[0],
    n_neighbors = params[1],
    min_dist = params[2],
    metric = params[3]
).fit_transform(X, y)

fig, ax = plt.subplots(figsize = (8, 6), dpi = 80, facecolor = 'w', edgecolor = 'k')
pc1 = ax.scatter(X_embedded[y == 0, 0], X_embedded[y == 0, 1], label = 'No AMD / normal')
pc2 = ax.scatter(X_embedded[y == 1, 0], X_embedded[y == 1, 1], label = 'Early stage AMD')
pc3 = ax.scatter(X_embedded[y == 2, 0], X_embedded[y == 2, 1], label = 'Intermediate AMD')
pc4 = ax.scatter(X_embedded[y == 3, 0], X_embedded[y == 3, 1], label = 'Advanced AMD')
ax.set_title("Training Data with Target (AMD Stage)")
ax.legend()
plt.show()


# In[17]:


# Actual case where no classification is provided
X_embedded = umap.UMAP(
    n_components = params[0],
    n_neighbors = params[1],
    min_dist = params[2],
    metric = params[3]
).fit_transform(X)

fig, ax = plt.subplots(figsize = (8, 6), dpi = 80, facecolor = 'w', edgecolor = 'k')
pc1 = ax.scatter(X_embedded[y == 0, 0], X_embedded[y == 0, 1], label = 'No AMD / normal')
pc2 = ax.scatter(X_embedded[y == 1, 0], X_embedded[y == 1, 1], label = 'Early stage AMD')
pc3 = ax.scatter(X_embedded[y == 2, 0], X_embedded[y == 2, 1], label = 'Intermediate AMD')
pc4 = ax.scatter(X_embedded[y == 3, 0], X_embedded[y == 3, 1], label = 'Advanced AMD')
ax.set_title("Training Data without Target (AMD Stage)")
ax.legend()
plt.show()


# In[18]:


# Divide data in training set into 4 group according to AMD stages
df_F = df[df.sex.str.contains('F')].drop(['rna_isolation_batch', 'library_prepper', 'sex', 'age'], axis = 1)
df_M = df[df.sex.str.contains('M')].drop(['rna_isolation_batch', 'library_prepper', 'sex', 'age'], axis = 1)
X = np.vstack([df_F, df_M])

label_F = np.zeros(len(df_F))
label_M = np.ones(len(df_M))
y = np.concatenate([label_F, label_M])


# In[19]:


# Ideal case where correct classification (target array) is provided
X_embedded = umap.UMAP(
    n_components = params[0],
    n_neighbors = params[1],
    min_dist = params[2],
    metric = params[3]
).fit_transform(X, y)

fig, ax = plt.subplots(figsize = (8, 6), dpi = 80, facecolor = 'w', edgecolor = 'k')
pc1 = ax.scatter(X_embedded[y == 0, 0], X_embedded[y == 0, 1], label = 'Female')
pc2 = ax.scatter(X_embedded[y == 1, 0], X_embedded[y == 1, 1], label = 'Male')
ax.set_title("Training Data with Target (Gender)")
ax.legend()
plt.show()


# In[20]:


# Actual case where no classification is provided
X_embedded = umap.UMAP(
    n_components = params[0],
    n_neighbors = params[1],
    min_dist = params[2],
    metric = params[3]
).fit_transform(X)

fig, ax = plt.subplots(figsize = (8, 6), dpi = 80, facecolor = 'w', edgecolor = 'k')
pc1 = ax.scatter(X_embedded[y == 0, 0], X_embedded[y == 0, 1], label = 'Female')
pc2 = ax.scatter(X_embedded[y == 1, 0], X_embedded[y == 1, 1], label = 'Male')
ax.set_title("Training Data without Target (Gender)")
ax.legend()
plt.show()


# In[21]:


# Divide data in training set into 4 group according to AMD stages
df_1 = df[df.rna_isolation_batch.str.contains('1')].drop(['rna_isolation_batch', 'library_prepper', 'sex', 'age'], axis = 1)
df_2 = df[df.rna_isolation_batch.str.contains('2')].drop(['rna_isolation_batch', 'library_prepper', 'sex', 'age'], axis = 1)
X = np.vstack([df_1, df_2])

label_1 = np.zeros(len(df_1))
label_2 = np.ones(len(df_2))
y = np.concatenate([label_1, label_2])


# In[22]:


# Ideal case where correct classification (target array) is provided
X_embedded = umap.UMAP(
    n_components = params[0],
    n_neighbors = params[1],
    min_dist = params[2],
    metric = params[3]
).fit_transform(X, y)

fig, ax = plt.subplots(figsize = (8, 6), dpi = 80, facecolor = 'w', edgecolor = 'k')
pc1 = ax.scatter(X_embedded[y == 0, 0], X_embedded[y == 0, 1], label = 'isobatch1')
pc2 = ax.scatter(X_embedded[y == 1, 0], X_embedded[y == 1, 1], label = 'isobatch2')
ax.set_title("Training Data with Target (RNA Isolation Batch)")
ax.legend()
plt.show()


# In[23]:


# Actual case where no classification is provided
X_embedded = umap.UMAP(
    n_components = params[0],
    n_neighbors = params[1],
    min_dist = params[2],
    metric = params[3]
).fit_transform(X)

fig, ax = plt.subplots(figsize = (8, 6), dpi = 80, facecolor = 'w', edgecolor = 'k')
pc1 = ax.scatter(X_embedded[y == 0, 0], X_embedded[y == 0, 1], label = 'isobatch1')
pc2 = ax.scatter(X_embedded[y == 1, 0], X_embedded[y == 1, 1], label = 'isobatch2')
ax.set_title("Training Data without Target (RNA Isolation Batch)")
ax.legend()
plt.show()


# In[24]:


# Divide data in training set into 4 group according to AMD stages
df_MRS = df[df.library_prepper.str.contains('MRS')].drop(['rna_isolation_batch', 'library_prepper', 'sex', 'age'], axis = 1)
df_RRP = df[df.library_prepper.str.contains('RRP')].drop(['rna_isolation_batch', 'library_prepper', 'sex', 'age'], axis = 1)
X = np.vstack([df_MRS, df_RRP])

label_MRS = np.zeros(len(df_MRS))
label_RRP = np.ones(len(df_RRP))
y = np.concatenate([label_MRS, label_RRP])


# In[25]:


# Ideal case where correct classification (target array) is provided
X_embedded = umap.UMAP(
    n_components = params[0],
    n_neighbors = params[1],
    min_dist = params[2],
    metric = params[3]
).fit_transform(X, y)

fig, ax = plt.subplots(figsize = (8, 6), dpi = 80, facecolor = 'w', edgecolor = 'k')
pc1 = ax.scatter(X_embedded[y == 0, 0], X_embedded[y == 0, 1], label = 'MRS')
pc2 = ax.scatter(X_embedded[y == 1, 0], X_embedded[y == 1, 1], label = 'RRP')
ax.set_title("Training Data with Target (Library Prepper)")
ax.legend()
plt.show()


# In[26]:


# Actual case where no classification is provided
X_embedded = umap.UMAP(
    n_components = params[0],
    n_neighbors = params[1],
    min_dist = params[2],
    metric = params[3]
).fit_transform(X)

fig, ax = plt.subplots(figsize = (8, 6), dpi = 80, facecolor = 'w', edgecolor = 'k')
pc1 = ax.scatter(X_embedded[y == 0, 0], X_embedded[y == 0, 1], label = 'MRS')
pc2 = ax.scatter(X_embedded[y == 1, 0], X_embedded[y == 1, 1], label = 'RRP')
ax.set_title("Training Data without Target (Libaray Prepper)")
ax.legend()
plt.show()


# In[ ]:




