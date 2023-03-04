library(tidyverse)
library(foreign)
library(modEvA)
library(caret)
library(splines)
library(cowplot)


setwd("/Users/philippehabets/Dropbox/STRESS_INDEX/")

######################### functions ########################
as.numeric.factor <- function(x){
  #function for turning factor values (levels) into numeric values
  as.numeric(levels(x))[x]
}

decompose_d2 <- function(m_conf, m_pred, m_conf_pred){
  d2_conf <- Dsquared(m_conf)
  d2_pred <- Dsquared(m_pred)
  d2_conf_pred <- Dsquared(m_conf_pred)
  
  conf_unexplained <- 1 - d2_conf
  pred_unexplained <- 1 - d2_pred
  
  delta_pred <- d2_conf_pred - d2_conf
  delta_conf <- d2_conf_pred - d2_pred
  
  partial_pred <- delta_pred / conf_unexplained
  partial_conf <- delta_conf / pred_unexplained
  
  shared = d2_conf_pred - delta_conf - delta_pred
  
  res <- c('confounds' = d2_conf,
           'predictions' = d2_pred,
           'confounds+predictions' = d2_conf_pred,
           'delta confounds' = delta_conf,
           'delta predictions' = delta_pred,
           'partial confounds' = partial_conf,
           'partial predicitons' = partial_pred,
           'shared' = shared)
  res <- as.data.frame(res)
  res$d2_type <- rownames(res)
  return(res)
}
###########################################################

# load ID's of train and test set to match confounder data
train_ids <- read_csv("scripts/Predictions_explorative/R.scripts/data/xgboost_prot_clin_train.csv") %>% 
  mutate(pident = as.factor(pident))
test_ids <- read_csv("scripts/Predictions_explorative/R.scripts/data/xgboost_prot_clin_test.csv") %>% 
  mutate(pident = as.factor(pident))

# Confounder set to see if ML model predicts over and beyond what is just in the confounder information:
##BMI
bmi <- read.spss("data/biological_variables/W1/HeightWeightBMI/N1_A357D.sav",
                     to.data.frame = TRUE, use.missings = FALSE) %>% 
  mutate(Pident = as.factor(Pident)) %>% 
  mutate(abmi = sapply(abmi, as.numeric.factor))

##Years of education & recruitment type
demographic <- read.spss("data/demographic_and_sociological_variables/W1/DOB_age_gender_nationality_and_education_of_respondents/N1_A100R.sav",
                             to.data.frame = TRUE, use.missings = FALSE) %>% 
  mutate(pident = as.factor(pident),
         female = ifelse(Sexe == "female", 1, 0)) %>% 
  select(pident, female, Age, aedu, aframe02)

##AD use at baseline (TCA, SSRI, or "Other antidepressants")
AD <- read.spss("/Users/philippehabets/Dropbox/STRESS_INDEX/data/biological_variables/W1/medication.current.use/N1_A354D.sav",
                    to.data.frame = TRUE, use.missings = FALSE) %>% 
  mutate(AD = "No") %>% 
  mutate_at(vars(AD), ~if_else((atca_fr != "No" | assri_fr != "No" | aother_ad_fr != "No"), "Yes", .x)) %>%
  mutate(PIDENT = as.factor(PIDENT)) %>% 
  select(PIDENT, AD) 

##Psychopharmaca use past 3 years
PF <- read.spss("/Users/philippehabets/Dropbox/STRESS_INDEX/data/biological_variables/W1/medication.antidepressants.past3yrs/N1_A355R.sav",
                    to.data.frame = TRUE, use.missings = FALSE) %>% 
  select(Pident, appfmuse) %>% 
  droplevels.data.frame()

# confounder dfs
confounders_train <- train_ids %>% 
  select(pident, Remitted_depression) %>% 
  left_join(AD, by = c("pident" = "PIDENT")) %>% 
  left_join(PF, by = c("pident" = "Pident")) %>% 
  left_join(demographic, by = "pident") %>% 
  left_join(bmi, by = c("pident" = "Pident")) %>% 
  select(-pident) %>% 
  mutate(Remitted_depression = ifelse(Remitted_depression == "remitted", 0, 1))

options(na.action = "na.pass")
modelMatrix <- model.matrix(Remitted_depression ~ ., data = confounders_train)[,-1] #transform categorical nominal features into numeric values with dummy variables (k-1)
confounders_train <- cbind(Remitted_depression = confounders_train$Remitted_depression, as.data.frame(modelMatrix))
names(confounders_train) <- str_replace_all(names(confounders_train), " ", "") #remove spaces in columnames

confounders_test <- test_ids %>% 
  select(pident, Remitted_depression) %>% 
  left_join(AD, by = c("pident" = "PIDENT")) %>% 
  left_join(PF, by = c("pident" = "Pident")) %>% 
  left_join(demographic, by = "pident") %>% 
  left_join(bmi, by = c("pident" = "Pident")) %>% 
  select(-pident) %>% 
  mutate(Remitted_depression = ifelse(Remitted_depression == "remitted", 0, 1))

options(na.action = "na.pass")
modelMatrix <- model.matrix(Remitted_depression ~ ., data = confounders_test)[,-1] #transform categorical nominal features into numeric values with dummy variables (k-1)
confounders_test <- cbind(Remitted_depression = confounders_test$Remitted_depression, as.data.frame(modelMatrix))
names(confounders_test) <- str_replace_all(names(confounders_test), " ", "") #remove spaces in columnames


# load train and test data
segTrain <- read_csv("scripts/Predictions_explorative/R.scripts/2-yearChronicity/output/segTrain_bagImputedPreProcessingProteomicsBatchCorrected_AND_clinical.csv") %>%
  mutate(Remitted_depression = factor(Remitted_depression, levels=c("remitted", "not_remitted"))) %>% as.data.frame()
segTest <- read_csv("scripts/Predictions_explorative/R.scripts/2-yearChronicity/output/segTest_bagImputedPreProcessingProteomicsBatchCorrected_AND_clinical.csv") %>%
  mutate(Remitted_depression = factor(Remitted_depression, levels=c("remitted", "not_remitted"))) %>% as.data.frame()

# load model
XGbTune <- readRDS("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/XGbTune2.RDS")

# predictions on test set
accuracies <- list()
for(i in seq(0.35, 0.65, 0.0001)){
  pred_ <- factor(ifelse(predict(XGbTune, segTest[,-1], type="prob")$remitted>i,"remitted","not_remitted"), levels = c("remitted", "not_remitted"))
  cm_ <- confusionMatrix(pred_, segTest[,1])
  accuracies <- append(accuracies, cm_$overall[1])
}
names(accuracies) <- seq(0.35, 0.65, 0.0001)
opt_thresh <- as.numeric(names(accuracies[accuracies == max(sapply(accuracies, max))]))[1]
XGbPred <- factor(ifelse(predict(XGbTune, segTest[,-1], type="prob")$remitted>opt_thresh,"remitted","not_remitted"), levels = c("remitted", "not_remitted"))

### check if model performs on the basis of patterns over and above predictability of confounder variables 
### code adjusted from: https://github.com/dinga92/confounds_paper/blob/master/example-analysis.Rmd

#df with predictions
predTest <- as.data.frame(cbind(Remitted_depression = as.numeric(segTest$Remitted_depression)-1, predicted = as.numeric(XGbPred)-1))

#combined df
combTest <- cbind.data.frame(predTest, confounders_test[,-1])

# fit models we will compare
# 1. model with only confounds as predictors  + (nonlinear expansions using splines)
m_conf <- glm(Remitted_depression ~
                ADYes +
                appfmuseNo + 
                female +
                bs(Age, df=5) +
                bs(aedu, df =5) +
                aframe02Specialisedmentalhealthcare +
                aframe02Generalpopulation +
                bs(abmi, df = 5),
              data=confounders_test, family = binomial)

# 2. model with only ML predictions as predictors
m_pred <- glm(Remitted_depression ~ predicted, data=predTest, family = binomial)

# 3. model with both ML predictions and confounds as covariates
m_conf_pred <- glm(Remitted_depression ~ ., data=combTest, family = binomial)

# get partial and delta r2 for all confounds
decomposed_d2 <- decompose_d2(m_conf, m_pred, m_conf_pred)

# inspect decomposed d2
decomposed_d2
summary(m_conf)
summary(m_pred)
summary(m_conf_pred)

# make a table for  future plotting
t1 <- decomposed_d2[c(1:5,8),]
t1$conf <- 'all'
t1 <- rbind(t1, t2)
t1

t_small <- t1[t1$d2_type %in% 
                c('delta predictions', 'shared', 'delta confounds'),]
t_small$d2_type <- factor(t_small$d2_type,
                          levels=rev(c('delta predictions',
                                       'shared',
                                       'delta confounds')),
                          labels=c('confounders only',
                                   'confounders + predictions', 
                                   'predictions only'
                          ))
p <- ggplot(t_small, aes(y=res, x=conf, fill=d2_type)) +
  geom_bar(stat='identity') +
  coord_flip() +
  labs(y=expression( D^{2}),
       x='Confounders') +
  theme_minimal_vgrid() +
  theme(aspect.ratio = 0.25) +
  scale_fill_manual(values = c("#F0E442", "#9CCA9D", "#56B4E9"),
                    name='Deviance explained by')
p



