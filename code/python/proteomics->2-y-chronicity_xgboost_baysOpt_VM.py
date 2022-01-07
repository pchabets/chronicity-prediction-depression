#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import pandas as pd
import numpy as np

from datetime import datetime

import pyreadstat as prs

import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import uniform, randint

from sklearn import preprocessing
from sklearn.metrics import auc, roc_curve, roc_auc_score, accuracy_score, confusion_matrix, plot_confusion_matrix
from sklearn.model_selection import RepeatedKFold, RepeatedStratifiedKFold, cross_val_score, GridSearchCV, KFold, RandomizedSearchCV, train_test_split

from skopt import BayesSearchCV

import xgboost as xgb


# In[2]:


import sklearn
sklearn.__version__


# In[3]:


os.chdir("/home/pchabets/Dropbox/STRESS_INDEX/")


# In[4]:


# timer function, copied from https://www.kaggle.com/tilii7/hyperparameter-grid-search-with-xgboost
def timer(start_time=None):
    if not start_time:
        start_time = datetime.now()
        return start_time
    elif start_time:
        thour, temp_sec = divmod((datetime.now() - start_time).total_seconds(), 3600)
        tmin, tsec = divmod(temp_sec, 60)
        print('\n Time taken: %i hours %i minutes and %s seconds.' % (thour, tmin, round(tsec, 2)))


# In[5]:


#proteomics
prtm = pd.read_csv("data/blood_and_saliva_variables/W1/proteomics/output/proteomics_replaced_outliers.csv")
prtm = prtm.drop('applate', axis=1)


# In[6]:


#label data
y = pd.read_spss("data/outcome_groups/DSM_groups.sav", convert_categoricals=True)
y = y.drop('Remitted_comorbid', axis=1)
y = y.rename(columns={'pident':'Pident'})
y['Pident'] = y['Pident'].astype(int)


# In[7]:


whole_set = pd.merge(y, prtm, how='inner', on='Pident')


# In[8]:


whole_set


# In[9]:


#Turn labels into 0's and 1's: non-remitted = 0, remitted = 1
from sklearn.preprocessing import LabelEncoder
lbl = LabelEncoder() 
whole_set['Remitted_depression'] = lbl.fit_transform(whole_set['Remitted_depression']) 


# In[10]:


X = whole_set.drop(['Pident','Remitted_depression'], axis=1)
y = whole_set.pop('Remitted_depression')


# ### log10 transform data, can be done before train/test split because pointwise operation

# In[11]:


X = np.log10(X) 


# ##### optional:  binarize analytes with > 90% lowest values (below detection limit values), using whole dataframe

# In[12]:


# X = X.apply(
#     lambda x: (x.apply(lambda y: 0 if y==min(x) else 1)) if sum(x==min(x))/whole_set.shape[0] > 0.9 else x
# )


# In[13]:


X


# ### Corr plot

# In[14]:


plt.figure(figsize = (18,12))
sns.heatmap(X.corr(), cmap='coolwarm', vmin=-1, vmax=1) 


# ### train-test split (80-20)

# In[14]:


from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=101, shuffle=True)


# ### Preprocessing: remove zero-variance variables and impute NaN's with median

# In[15]:


from sklearn.feature_selection import VarianceThreshold
from sklearn.impute import SimpleImputer


# #### Select ZV columns from trainings data and remove them from train + test X

# In[16]:


selector = VarianceThreshold()
selector.fit(X_train)


# In[17]:


#Store variable names for later conversion back to dataframe
X_names = list(X_train.columns[selector.get_support()])


# In[18]:


X_train = selector.transform(X_train)
X_test = selector.transform(X_test)


# In[19]:


print(str(X.shape[1]-X_train.shape[1]) + " analytes with zero variance")


# #### Calculate median for each column in train X, and replace NaNs with that value in both train and test X

# In[20]:


imputer = SimpleImputer(strategy="median")


# In[21]:


imputer.fit(X_train)


# In[22]:


X_train = imputer.transform(X_train)
X_test = imputer.transform(X_test)


# In[23]:


# Convert back to dataframe
X_train = pd.DataFrame(data=X_train,columns=X_names) 
X_test = pd.DataFrame(data=X_test, columns=X_names)


# ### Model
# Bayesian optimization of hyperparameter tuning is used. Code was adapted from https://www.kaggle.com/nanomathias/bayesian-optimization-of-xgboost-lb-0-9769

# In[24]:


#set number of parallel jobs to execute with the Bayes Search
N_JOBS = 14


# In[25]:


#set number of iterations of search job
ITERATIONS = 1000


# For more info on hyperparameters: see https://xgboost.readthedocs.io/en/latest/parameter.html

# In[26]:


# Classifier
bayes_cv_tuner = BayesSearchCV(
    estimator = xgb.XGBClassifier(
        n_jobs = 1,
        objective = 'binary:logistic',
        booster = 'gbtree',
        eval_metric = 'auc', 
        use_label_encoder=False,
        verbosity=1
    ),
    search_spaces = {
        'learning_rate': (0.01, 1.0, 'log-uniform'),
        'min_child_weight': (0, 16, 'uniform'),
        'max_depth': (1, 6, 'uniform'),
#         'max_delta_step': (0, 1),
        'subsample': (0.5, 1.0, 'uniform'),
        'colsample_bytree': (0.01, 1.0, 'uniform'),
#         'colsample_bylevel': (0.01, 1.0, 'uniform'),
#         'colsample_bynode': (0.01, 1.0, 'uniform'),
#         'reg_lambda': (1e-9, 1000.0, 'log-uniform'),
#         'reg_alpha': (1e-9, 1.0, 'log-uniform'),
        'gamma': (0.01, 3.0, 'log-uniform'),
        'n_estimators': (50, 1000, 'uniform')
    },    
    cv = RepeatedKFold(
        n_splits=10, 
        n_repeats=10, 
        random_state=42
    ),
    n_jobs = N_JOBS,
    n_iter = ITERATIONS,   
    verbose= 3,
    random_state = 42,
    refit = True    
)


# In[27]:


def status_print(optim_result):
    """Status callback during bayesian hyperparameter search"""
    
    # Get all the models tested so far in DataFrame format
    all_models = pd.DataFrame(bayes_cv_tuner.cv_results_)    
    
    # Get current parameters and the best parameters    
    best_params = pd.Series(bayes_cv_tuner.best_params_)
    print('Model #{}\nBest ROC-AUC: {}\nBest params: {}\n'.format(
        len(all_models),
        np.round(bayes_cv_tuner.best_score_, 4),
        bayes_cv_tuner.best_params_
    ))


# ### Fit the model

# In[ ]:


#fit the model
start_time = timer(None) 
fitted_model = bayes_cv_tuner.fit(X_train, y_train, callback=status_print)
timer(start_time)


# In[ ]:


print('\n Best estimator:')
print(fitted_model.best_estimator_)
print('\n Best hyperparameters:')
print(fitted_model.best_params_)

results = pd.DataFrame(fitted_model.cv_results_)


# In[ ]:


#Make predictions with best fitted model
y_pred = fitted_model.best_estimator_.predict(X_test)
y_scores = fitted_model.best_estimator_.predict(X_test, output_margin=True)


# In[ ]:


#Plot confusion matrix
plot_confusion_matrix(fitted_model, X_test, y_test, normalize='true', cmap='viridis')


# In[ ]:


print('\nAccuracy: \n\n'+str(accuracy_score(y_test, y_pred)))


# ### Save model

# In[28]:


import joblib


# In[29]:


# joblib.dump(fitted_model.best_estimator_, "scripts/VM/Python/output/BaysOpt_model_best_estimator_1000.pkl")
# joblib.dump(fitted_model.best_params_,"scripts/VM/Python/output/BaysOpt_model_best_params_1000.pkl")
# joblib.dump(fitted_model, "scripts/VM/Python/output/BaysOpt_model_1000.sav")


# In[30]:


#load model
best_estimator = joblib.load("scripts/VM/Python/output/BaysOpt_model_best_estimator_1000.pkl")
best_params = joblib.load("scripts/VM/Python/output/BaysOpt_model_best_params_1000.pkl")
model = joblib.load("scripts/VM/Python/output/BaysOpt_model_1000.sav")


# In[ ]:





# In[31]:


best_estimator


# ### Predictions

# In[32]:


#Make predictions with best fitted model
y_pred = best_estimator.predict(X_test)
y_scores = best_estimator.predict(X_test, output_margin=True)


# In[33]:


#Confusion matrix
cm = confusion_matrix(y_test, y_pred, normalize='pred')
print('Confusion Matrix: \n\n'+str(cm))
print('\nAccuracy: \n\n'+str(accuracy_score(y_test, y_pred)))


# In[34]:


#Plot confusion matrix
import seaborn as sns
plt.figure(figsize = (8,6))
ax = plt.axes()
sns.set(font_scale=1.5)
sns.set_style("white")
sns.heatmap(cm, 
            annot=True,
            vmin=0,
            vmax=1,
            cmap='coolwarm',
            xticklabels=['non-remitted', 'remitted'], 
            yticklabels=['non-remitted', 'remitted'])
ax.set(xlabel='PREDICTED', ylabel = 'ACTUAL')


# ### Plot ROC curve

# In[35]:


from sklearn.metrics import plot_roc_curve
plt.clf()
plt.figure(figsize=(7,7))
axes = plt.axes()
plot_roc_curve(best_estimator, X_test, y_test, ax=axes)  
plt.plot([0,1],[0,1],'r--') 


# In[37]:


roc_auc_score(y_test, y_scores)


# In[37]:


from sklearn.metrics import f1_score, accuracy_score, balanced_accuracy_score, recall_score, precision_score, auc

print("Precision: {}".format(precision_score(y_test,y_pred).round(2)))
print("Recall: {}".format(recall_score(y_test,y_pred).round(2)))
print("F1 Score: {}".format(f1_score(y_test,y_pred).round(2)))
print("Accuracy: {}".format(accuracy_score(y_test,y_pred).round(2)))
print("Balanced accuracy: {}".format(balanced_accuracy_score(y_test,y_pred).round(2)))


# ### Model explanation

# In[42]:


import shap
import statistics
# non-remitted = 0, remitted = 1


# In[43]:


# explain the model's predictions using SHAP
explainer = shap.Explainer(best_estimator)
shap_values = explainer(X_train)


# In[44]:


# dataframe of shapley values per case per analyte
shap_df = pd.DataFrame(shap_values.values, columns=shap_values.feature_names)


# In[45]:


# average absolute shap value per analyte
avg_shap_abs = shap_df.apply(lambda x: statistics.mean(abs(x)),axis=0)
avg_shap_abs = avg_shap_abs.sort_values(ascending=False)
avg_shap_abs = pd.DataFrame(avg_shap_abs, columns=["Avg Absolute Shapley Value"])
avg_shap_abs.reset_index(inplace=True)
avg_shap_abs = avg_shap_abs.rename(columns = {'index':'Protein'})
avg_shap_abs['Protein'] = avg_shap_abs['Protein'].apply(lambda x: x[2:])
avg_shap_abs


# In[46]:


#write to file
avg_shap_abs.to_csv("scripts/VM/Python/output/avg_shap_abs_VM.csv")


# In[47]:


# average shap value per analyte, keeping 'direction'
avg_shap = shap_df.apply(lambda x: statistics.mean(x),axis=0)
avg_shap = avg_shap.reindex(avg_shap.abs().sort_values(ascending=False).index)
avg_shap = pd.DataFrame(avg_shap, columns=["Average Shapley Value"])
avg_shap.reset_index(inplace=True)
avg_shap = avg_shap.rename(columns = {'index':'Protein'})
avg_shap['Protein'] = avg_shap['Protein'].apply(lambda x: x[2:])
avg_shap


# In[48]:


#write to file
avg_shap.to_csv("scripts/VM/Python/output/avg_shap_direction_VM.csv")


# In[49]:


# visualize shap values for each analyte across instances
shap.initjs()
shap.plots.heatmap(shap_values, max_display=10)


# In[50]:


# visualize top 10 proteins with highest mean absolute shap values
shap.initjs()
shap.plots.bar(shap_values, max_display=10)


# In[51]:


# visualize all prediction's explanation
shap.initjs()
shap.plots.force(explainer.expected_value, shap_values.values, features=X_train, feature_names=X_names)


# In[52]:


# visualize the first prediction's explanation
shap.initjs()
shap.plots.force(shap_values[0],feature_names=X_names)


# In[53]:


# visualize the first prediction's explanation as waterfall plot
shap.initjs()
shap.plots.waterfall(shap_values[0], max_display=10)


# ### Alternatively, build exact explainer using masker

# In[54]:


# build a clustering of the features based on shared information about y
clustering = shap.utils.hclust(X_train, y_train)


# In[ ]:


# above we implicitly used shap.maskers.Independent by passing a raw dataframe as the masker
# now we explicitly use a Partition masker that uses the clustering we just computed
masker = shap.maskers.Partition(X_train, clustering=clustering)

# build an Exact explainer and explain the model predictions on the given dataset
explainer = shap.explainers.Exact(best_estimator.predict_proba, masker)
shap_values2 = explainer(X_train[:232])

# get just the explanations for the positive class
shap_values2 = shap_values2[...,1]


# In[ ]:


shap_values2


# In[ ]:




