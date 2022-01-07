#!/usr/bin/env python
# coding: utf-8

# ### For futher info on BorutaPy package see: <br><br>https://github.com/scikit-learn-contrib/boruta_py 
# ##### Articles:
# #### BorutaPy: <br> https://www.jstatsoft.org/article/view/v036i11
# #### Robustness of RF-based feature selection: <br> https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-8 <br>
# 

# In[39]:


import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


# In[40]:


os.chdir("/home/pchabets/Dropbox/STRESS_INDEX/")


# ### Load in top-variance-selected transcriptomics data

# In[41]:


expr_train = pd.read_csv("scripts/VM/Python/output/transcriptomics_variance_selection_top_5000_TRAIN.csv")
expr_train = expr_train.iloc[:,1:]
expr_test = pd.read_csv("scripts/VM/Python/output/transcriptomics_variance_selection_top_5000_TEST.csv")
expr_test = expr_test.iloc[:,1:]


# ### Use Boruta for variabe selection

# #### Turn labels into 0's and 1's 
# non-remitted = 0, remitted = 1

# In[42]:


#not_remitted = 0, remitted = 1

from sklearn.preprocessing import LabelEncoder
lbl = LabelEncoder() 
expr_train['Remitted_depression'] = lbl.fit_transform(expr_train['Remitted_depression'])


# In[43]:


# Base Boruta feature selection only on train set - X and y are from trainset only
X = expr_train.drop(['Remitted_depression', 'pident'], inplace=False, axis=1)
y = expr_train['Remitted_depression']


# #### initialize Boruta

# In[44]:


# set max depth and percentage
max_depth = 3
percentage = 95


# In[45]:


from boruta import BorutaPy
from sklearn.ensemble import RandomForestClassifier


# In[46]:


rf = RandomForestClassifier(
    n_jobs=14,
    max_depth = max_depth,
    class_weight='balanced',
    verbose=1
)


# In[47]:


feat_selector = BorutaPy(
    estimator = rf,
    n_estimators='auto',
    max_iter=1000,
    perc=percentage,
#     random_state=101,
    verbose=2
)


# #### fit Boruta

# In[48]:


feat_selector.fit(np.array(X), np.array(y))


# ### Check results

# In[49]:


# check selected features
feat_selector.support_


# In[50]:


# check ranking of features
feat_selector.ranking_


# In[51]:


strong = X.columns[feat_selector.support_].to_list()
weak = X.columns[feat_selector.support_weak_].to_list()
print('Selected features:', strong)
print('Potentially irrelevant features:', weak)


# ### Transform X to selected features X

# In[52]:


# only keep confirmed features, discard tentative features 
X_filtered = X.loc[:,feat_selector.support_]


# In[53]:


X_filtered.shape


# ### Combine into new dataframe with labels and write to file

# #### Train set

# In[54]:


# participant IDs of train set
pid_train = expr_train['pident']


# In[55]:


# concatenate into 1 df
selected_set_TRAIN = pd.concat([pid_train, y, X_filtered], axis=1)


# In[56]:


# transform labels back to characters, not_remitted = 0, remitted = 1
selected_set_TRAIN['Remitted_depression'] = selected_set_TRAIN['Remitted_depression'].transform(lambda x: 'not_remitted' if x == 0 else 'remitted')


# In[57]:


# save to file
selected_set_TRAIN.to_csv("scripts/VM/Python/output/transcriptomics_TopVariance_boruta_selection_{}%_TRAIN_maxdepth{}.csv".format(percentage, max_depth))


# #### Test set

# In[58]:


# participaexpr_testof test set
pid_test = expr_test['pident'] 

# labels of test set
y_test = expr_test['Remitted_depression']

# features of test set
X_test = expr_test.drop(['pident', 'Remitted_depression'], inplace=False, axis=1) 


# In[59]:


# filtered features of test set
X_test_filtered = X_test.loc[:,feat_selector.support_] 


# In[60]:


# concatenate into 1 df
selected_set_TEST = pd.concat([pid_test, y_test, X_test_filtered], axis=1) 


# In[61]:


# save to file
selected_set_TEST.to_csv("scripts/VM/Python/output/transcriptomics_TopVariance_boruta_selection_{}%_TEST_maxdepth{}.csv".format(percentage, max_depth))


# In[62]:


print("done")


# In[ ]:




