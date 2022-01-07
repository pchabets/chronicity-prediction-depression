#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import tensorflow as tf

import autokeras as ak

os.chdir("/home/pchabets/Dropbox/STRESS_INDEX/")


# In[2]:


#proteomics
prtm = pd.read_csv("data/blood_and_saliva_variables/W1/proteomics/output/proteomics_replaced_outliers.csv")
prtm = prtm.drop('applate', axis=1)


# In[3]:


#label data
y = pd.read_spss("data/outcome_groups/DSM_groups.sav", convert_categoricals=True)
y = y.drop('Remitted_comorbid', axis=1)
y = y.rename(columns={'pident':'Pident'})
y['Pident'] = y['Pident'].astype(int)


# In[4]:


whole_set = pd.merge(y, prtm, how='inner', on='Pident')


# In[5]:


whole_set


# In[6]:


#Turn labels into 0's and 1's: non-remitted = 0, remitted = 1
from sklearn.preprocessing import LabelEncoder
lbl = LabelEncoder() 
whole_set['Remitted_depression'] = lbl.fit_transform(whole_set['Remitted_depression']) 


# In[7]:


X = whole_set.drop(['Pident','Remitted_depression'], axis=1)
y = whole_set.pop('Remitted_depression')


# ### log10 transform data, can be done before train/test split because pointwise operation

# In[8]:


X = np.log10(X) 


# ### Train-test split (80-20)

# In[9]:


from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=101, shuffle=True)


# ### Preprocessing: remove zero-variance variables and impute NaN's with median

# In[10]:


from sklearn.feature_selection import VarianceThreshold
from sklearn.impute import SimpleImputer


# #### Select ZV columns from trainings data and remove them from train + test X

# In[11]:


selector = VarianceThreshold()
selector.fit(X_train)


# In[12]:


X_names = list(X_train.columns[selector.get_support()])


# In[13]:


X_train = selector.transform(X_train)
X_test = selector.transform(X_test)


# In[14]:


print(str(X.shape[1]-X_train.shape[1]) + " analytes with zero variance")


# #### Calculate median for each column in train X, and replace NaNs with that value in both train and test X

# In[15]:


imputer = SimpleImputer(strategy="median")


# In[16]:


imputer.fit(X_train)


# In[17]:


X_train = imputer.transform(X_train)
X_test = imputer.transform(X_test)


# #### Scale data with MinMax scaling

# In[18]:


from sklearn.preprocessing import MinMaxScaler


# In[19]:


scaler = MinMaxScaler()


# In[20]:


scaler.fit(X_train)


# In[21]:


X_train = scaler.transform(X_train)
X_test = scaler.transform(X_test)


# In[22]:


# Convert back to dataframe
X_train = pd.DataFrame(data=X_train,columns=X_names) 
X_test = pd.DataFrame(data=X_test, columns=X_names)


# ## Creating the model

# In[23]:


clf = ak.StructuredDataClassifier(overwrite=False, 
                                  max_trials=100, 
                                  objective='val_accuracy',
                                  directory="/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/Python.scripts/output"
                                 )


# In[24]:


# Feed the structured data classifier with training data.
clf.fit(X_train, y_train, 
#         epochs=10, 
        validation_split=0.3)


# In[68]:


# model = clf.export_model()


# In[69]:


model.summary()


# ## Model evaluation

# In[76]:


# Predict with the best model.
predicted_y = clf.predict(X_test)


# In[27]:


# Evaluate the best model with testing data.
clf.evaluate(X_train, y_train)


# In[28]:


model.evaluate(X_test, y_test)


# In[75]:


predicted_y


# In[30]:


from sklearn.metrics import classification_report,confusion_matrix


# In[31]:


print(classification_report(y_test,y_pred=predicted_y))
print(confusion_matrix(y_test,predicted_y, normalize="true"))


# ### Save model 

# In[34]:


# try:
#     model.save("/scripts/VM/Python/output/model_autokeras", save_format="tf")
# except Exception:
#     model.save("/scripts/VM/Python/output/model_autokeras.h5")


# ### Load model for further plotting/analysis

# In[37]:


loaded_model = tf.keras.models.load_model("scripts/VM/Python/output/model_autokeras", 
                                          custom_objects=ak.CUSTOM_OBJECTS
                                         )


# In[41]:


loaded_model.evaluate(X_test, y_test)


# ### Plot ROC curve and calculate AUC

# In[70]:


from sklearn.metrics import plot_roc_curve, roc_curve, auc, plot_confusion_matrix


# In[61]:


y_pred_keras = loaded_model.predict(X_test).ravel()
fpr_keras, tpr_keras, thresholds_keras = roc_curve(y_test, y_pred_keras)
auc_keras = auc(fpr_keras, tpr_keras)


# In[67]:


plt.figure(figsize=(6,6))
plt.plot([0, 1], [0, 1], 'k--')
plt.plot(fpr_keras, tpr_keras, label='AK-model (AUC = {:.3f})'.format(auc_keras))
plt.xlabel('False positive rate')
plt.ylabel('True positive rate')
plt.title('ROC curve')
plt.legend(loc='best')
plt.show()


# ### Plot confusion matrix

# In[111]:


cm = confusion_matrix(y_test, predicted_y, normalize="pred", labels=[0, 1])


# In[112]:


#non-remitted = 0, remitted = 1
pd.DataFrame(cm)


# In[174]:


#Plot confusion matrix
import seaborn as sns
plt.figure(figsize = (10,7))
ax = plt.axes()
sns.heatmap(cm, 
            annot=True,
            vmin=0,
            vmax=1,
            cmap='coolwarm',
            xticklabels=['non-remitted', 'remitted'], 
            yticklabels=['non-remitted', 'remitted'])
ax.set(xlabel='PREDICTED', ylabel = 'ACTUAL')
sns.set(font_scale=1.8)


# In[176]:


from sklearn.metrics import f1_score, accuracy_score, balanced_accuracy_score, recall_score, precision_score

print("Precision: {}".format(precision_score(y_test,predicted_y).round(2)))
print("Recall: {}".format(recall_score(y_test,predicted_y).round(2)))
print("F1 Score: {}".format(f1_score(y_test,predicted_y).round(2)))
print("Accuracy: {}".format(accuracy_score(y_test,predicted_y).round(2)))
print("Balanced accuracy: {}".format(balanced_accuracy_score(y_test,predicted_y).round(2)))

