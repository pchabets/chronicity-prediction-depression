library(tidyverse)
library(caret)
library(xgboost)
library(SHAPforxgboost)

# load in train and test data
segTrain <- read_csv("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/output/segTrain_bagImputedPreProcessingProteomicsBatchCorrected_AND_clinical.csv") %>%
  mutate(Remitted_depression = factor(Remitted_depression, levels=c("remitted", "not_remitted"))) %>% as.data.frame()
segTest <- read_csv("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/output/segTest_bagImputedPreProcessingProteomicsBatchCorrected_AND_clinical.csv") %>%
  mutate(Remitted_depression = factor(Remitted_depression, levels=c("remitted", "not_remitted"))) %>% as.data.frame()

## Load model
XGbTune <- readRDS("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/XGbTune2.RDS")

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
# write_excel_csv(analyte_importance, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/analyte_importance_XGB-Gain.csv")

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
colnames(dataX) <- c('sex', tolower(colnames(dataX)[2:10]), colnames(dataX)[11:length(colnames(dataX))])

feature_names <- XGbTune$finalModel$feature_names 
clean_names <- sapply(feature_names, strip_ap)
clean_names <- c('sex', tolower(clean_names[2:10]), clean_names[11:length(clean_names)])
XGbTune$finalModel$feature_names <- clean_names

shap_values <- shap.values(xgb_model = XGbTune$finalModel, X_train = dataX)
shap_values$mean_shap_score

plot <- shap.plot.summary.wrap1(XGbTune$finalModel,dataX, top_n=10, dilute = F)
plot 
# ggsave("/Users/philippehabets/Dropbox/STRESS_INDEX/manuscripts/prediction_2ychronicity/Figures/output/shap_prot_clin.pdf",
#        plot,
#        width = 15, height = 15, units = "cm")

analyte_shap <- data.frame(symbol = sapply(names(shap_values$mean_shap_score), strip_ap), 
                           mean_shap_score =  shap_values$mean_shap_score) %>% 
  left_join(y = analyte_names, by = c("symbol" = "Analyte abbreviation")) %>% 
  select(1,3,2)
# write_excel_csv(analyte_shap, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/shap_xgb_proteomicAndClinical.csv")



############################################################################################################
############################################################################################################
## Fit logistic model to look at linear relation between FBG and response variable chronicity
train_inference <- segTrain %>% 
  mutate(Remitted_depression = ifelse(Remitted_depression == "remitted", 0, 1))

summary(glm(formula = Remitted_depression ~ apFBG + 
              Sexefemale + 
              Age + 
              aids + aidssev + 
              apIL15 + apIL1a + apACE + apCFHR1 + apFetA + apIL2ra + apB2M, 
            data = train_inference, family = "binomial"))

test_inference <- segTest %>% 
  mutate(Remitted_depression = ifelse(Remitted_depression == "remitted", 0,1))

summary(glm(formula = Remitted_depression ~ apFBG + 
              Sexefemale + 
              Age + 
              aids + aidssev + 
              apIL15 + apIL1a + apACE + apCFHR1 + apFetA + apIL2ra + apB2M, 
            data = test_inference, family = "binomial"))

# Group effect
train_remitted <- segTrain %>% 
  filter(Remitted_depression == "remitted")
train_chronic <- segTrain %>% 
  filter(Remitted_depression == "not_remitted")

mean(train_remitted$apFBG)
mean(train_chronic$apFBG)


