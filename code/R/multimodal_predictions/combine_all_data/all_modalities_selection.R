library(tidyverse)
library(caret)

setwd("/Users/philippehabets/Dropbox/STRESS_INDEX/")

# read in segTrain and segTest (preprocessed in combine_all_data.R)
segTrain <- read_csv("data/multimodal/wave1_omics_clinical_segTrain.csv")
segTest <- read_csv("data/multimodal/wave1_omics_clinical_segTest.csv")

# Proteomic selection
proteomic_shap <- readRDS("scripts/Predictions_explorative/R.scripts/2-yearChronicity/proteomics/data/shap_values_proteomics.RDS")

# Transcriptomic selection
transcriptomic <- readRDS("scripts/Predictions_explorative/R.scripts/2-yearChronicity/transcriptomics/data/trancriptomic_selection_boruta.RDS")

# load data to know what columns to deselect 
prtm <- read_csv("data/blood_and_saliva_variables/W1/proteomics/output/proteomics_replaced_outliers.csv")
prtm <- colnames(prtm)[3:ncol(prtm)]
ge <- read_csv('data/blood_and_saliva_variables/W1/transcriptomics/transcriptomics_2ychronicity.csv')
ge <- colnames(ge)[4:ncol(ge)]

# List of proteomics to exclude
prtm_after_preprocessing <- prtm[prtm %in% colnames(segTrain)]
excl_prot <- prtm_after_preprocessing[!prtm_after_preprocessing %in% proteomic_shap]

# list of transcriptomics to exlcude
excl_GE <- ge[!ge %in% transcriptomic]

# remove abbundant proteomics and transcriptomics
train <- segTrain %>% 
  select(-c(excl_GE, excl_prot))
test <- segTest %>% 
  select(-c(excl_GE, excl_prot))

# load in columns from metabolomic data so we have a list to select for PCA preprocessing
mtbl <- read.csv("data/blood_and_saliva_variables/W1/metabolomics/output/mtbl_measurements.csv")
mtbl <- colnames(mtbl)[-1]

# compute PCA of lipids and metabolites based on train set
pp <- preProcess(train[mtbl], method = 'pca', thresh = 0.95)
transformedTrain <- predict(pp, newdata = train[mtbl])
transformedTest <- predict(pp, newdata = test[mtbl])

# remove mtbl and replace with PCs
train_final <- train %>% 
  select(-mtbl) %>% 
  cbind(transformedTrain)

test_final <- test %>% 
  select(-mtbl) %>% 
  cbind(transformedTest)
  
# write to file
# write_csv(train_final, "data/multimodal/wave1_omics_clinical_segTrain_final_selection.csv")
# write_csv(test_final, "data/multimodal/wave1_omics_clinical_segTest_final_selection.csv")










