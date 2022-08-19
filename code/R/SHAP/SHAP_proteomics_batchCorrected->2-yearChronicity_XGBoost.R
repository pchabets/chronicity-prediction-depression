library(tidyverse)
library(caret)
library(xgboost)
library(SHAPforxgboost)

setwd("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/")

# load in train and test data
segTrain <- read_rds("R.scripts/2-yearChronicity/output/segTrain_bagImputedPreProcessingProteomicsBatchCorrected.RDS") %>% 
  select(-c(pident, applate))
segTest <- read_rds("R.scripts/2-yearChronicity/output/segTest_bagImputedPreProcessingProteomicsBatchCorrected.RDS") %>% 
  select(-c(pident, applate))

## Load model
XGbTune <- readRDS("R.scripts/2-yearChronicity/output/XGbTune_batchCorrected.RDS")

## Gain variable importance
var_imp <- varImp(XGbTune, scale = FALSE)
var_imp$importance
analyte_names <- read_csv("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/analyte_names.csv")

analyte_importance <- var_imp$importance %>% 
  mutate(symbol = rownames(var_imp$importance)) %>% 
  mutate(symbol = ifelse(symbol %in% names(segTrain)[2:11], symbol, str_sub(rownames(var_imp$importance), start = 3))) %>% 
  left_join(y = analyte_names, by = c("symbol" = "Analyte abbreviation")) %>% 
  dplyr::rename(Gain = Overall) %>% 
  mutate(Gain = round(Gain, 5)) %>% 
  select(2,3,1)
View(analyte_importance)

########  #SHAP (Shapley Additive Explanation) of variables ########
# Top 10 in plot with buildin plot function of xgboost package
xgb.plot.shap(data = as.matrix(segTrain[XGbTune$finalModel$feature_names]), model = XGbTune$finalModel, top_n = 10)
  
# With SHAP for xgboost package:
# function to clean feature names
strip_ap <- function(x) {
  #function to strip variable names of "a" or "ap" if protein
  if(str_starts(x, "a")) {
    #if string starts with "a", it is either a protein or clinical variable
    if(str_starts(x, "ap")){
      #if string starts with "ap", it's a protein, strip fist two letters
      x <- str_sub(x, start = 3)
    } else {
      #if string starts with only "a", it's a clinical variable, so strip first letter
      x <- str_sub(x, start = 2)
    }
  } else {
    #if string does not start with "a", it is either 'Age' or "Sex" variable: keep all letters
    x <- x
  }
  return(x)
}

dataX <- as.matrix(segTrain[,XGbTune$finalModel$feature_names])
colnames(dataX) <- sapply(colnames(dataX), strip_ap)

feature_names <- XGbTune$finalModel$feature_names 
clean_names <- sapply(feature_names, strip_ap)
XGbTune$finalModel$feature_names <- clean_names

shap_values <- shap.values(xgb_model = XGbTune$finalModel, X_train = dataX)
shap_values$mean_shap_score

plot <- shap.plot.summary.wrap1(XGbTune$finalModel,dataX, top_n=10, dilute = F)
plot 
# ggsave("/Users/philippehabets/Dropbox/STRESS_INDEX/manuscripts/prediction_2ychronicity/Figures/output/shap_proteomics.pdf",
#        plot,
#        width = 15, height = 15, units = "cm")

analyte_shap <- data.frame(symbol = sapply(names(shap_values$mean_shap_score), strip_ap), 
                           mean_shap_score =  shap_values$mean_shap_score) %>% 
  left_join(y = analyte_names, by = c("symbol" = "Analyte abbreviation")) %>% 
  select(1,3,2)
# write_excel_csv(analyte_shap, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/shap_xgb_proteomicsOnly.csv")


