library(tidyverse)
library(plyr)
library(foreign)
library(caret)
library(xgboost)
library(SHAPforxgboost)
library(doMC)
library(pROC)
setwd("/Users/philippehabets/Dropbox/STRESS_INDEX/data/")

#Initialize parallel processing
registerDoMC(detectCores()-2)
getDoParWorkers()

###function for visualizing confusionMatrix:
draw_confusion_matrix <- function(cm, ratio=FALSE) {
  
  total <- sum(cm$table)
  if(ratio == FALSE){
    res <- as.numeric(cm$table)
  } else {
    res <- as.numeric(cm$table)
    res <- c(res[1]/(res[1]+res[2]), 1-res[1]/(res[1]+res[2]), res[3]/(res[3]+res[4]), 1-res[3]/(res[3]+res[4]))
    res <- round(res, 2)
  }
  
  # Generate color gradients. Palettes come from RColorBrewer.
  greenPalette <- c("#F7FCF5","#E5F5E0","#C7E9C0","#A1D99B","#74C476","#41AB5D","#238B45","#006D2C","#00441B")
  redPalette <- c("#FFF5F0","#FEE0D2","#FCBBA1","#FC9272","#FB6A4A","#EF3B2C","#CB181D","#A50F15","#67000D")
  getColor <- function (greenOrRed = "green", amount = 0) {
    if (amount == 0)
      return("#FFFFFF")
    palette <- greenPalette
    if (greenOrRed == "red")
      palette <- redPalette
    colorRampPalette(palette)(100)[10 + ceiling(90 * amount / total)]
  }
  
  # set the basic layout
  layout(matrix(c(1,1,2)))
  par(mar=c(2,2,2,2))
  plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  title('CONFUSION MATRIX', cex.main=2)
  
  # create the matrix 
  classes = colnames(cm$table)
  rect(150, 430, 240, 370, col=getColor("green", res[1]))
  text(195, 435, classes[1], cex=1.2)
  rect(250, 430, 340, 370, col=getColor("red", res[3]))
  text(295, 435, classes[2], cex=1.2)
  text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
  text(245, 450, 'Actual', cex=1.3, font=2)
  rect(150, 305, 240, 365, col=getColor("red", res[2]))
  rect(250, 305, 340, 365, col=getColor("green", res[4]))
  text(140, 400, classes[1], cex=1.2, srt=90)
  text(140, 335, classes[2], cex=1.2, srt=90)
  
  # add in the cm results
  text(195, 400, res[1], cex=1.6, font=2, col='black')
  text(195, 335, res[2], cex=1.6, font=2, col='black')
  text(295, 400, res[3], cex=1.6, font=2, col='black')
  text(295, 335, res[4], cex=1.6, font=2, col='black')
  
  # add in the specifics 
  plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
  text(10, 85, names(cm$byClass[1]), cex=1.2, font=2)
  text(10, 70, round(as.numeric(cm$byClass[1]), 3), cex=1.2)
  text(30, 85, names(cm$byClass[2]), cex=1.2, font=2)
  text(30, 70, round(as.numeric(cm$byClass[2]), 3), cex=1.2)
  text(50, 85, names(cm$byClass[5]), cex=1.2, font=2)
  text(50, 70, round(as.numeric(cm$byClass[5]), 3), cex=1.2)
  text(70, 85, names(cm$byClass[6]), cex=1.2, font=2)
  text(70, 70, round(as.numeric(cm$byClass[6]), 3), cex=1.2)
  text(90, 85, names(cm$byClass[7]), cex=1.2, font=2)
  text(90, 70, round(as.numeric(cm$byClass[7]), 3), cex=1.2)
  
  # add in the accuracy information 
  text(10, 35, names(cm$overall[1]), cex=1.5, font=2)
  text(10, 20, round(as.numeric(cm$overall[1]), 3), cex=1.4)
  text(50, 35, names(cm$byClass[11]), cex=1.5, font=2)
  text(50, 20, round(as.numeric(cm$byClass[11]), 3), cex=1.4)
  text(90, 35, names(cm$overall[2]), cex=1.5, font=2)
  text(90, 20, round(as.numeric(cm$overall[2]), 3), cex=1.4)
}
###############

##create dataframe from all data
w3Chronicity <- read.spss("outcome_groups/DSM_groups.sav", to.data.frame = TRUE, use.missings = FALSE)
w3Chronicity <- w3Chronicity[,c(1,3)]

#age, gender, years of education
demographic.raw <- read.spss("demographic_and_sociological_variables/W1/DOB_age_gender_nationality_and_education_of_respondents/N1_A100R.sav", to.data.frame = TRUE, use.missings = FALSE)
demographic <- demographic.raw %>% 
  mutate(pident = as.factor(pident)) %>% 
  select(pident, Sexe, Age, aedu)

#NEO-big five 
bigFive_raw <- read.spss("psychological_and_lifestyle_variables/W1/NEO-FFI_big_five_personality_test/N1_A240D.sav", to.data.frame = TRUE, use.missings = FALSE)
toNA <- function(x){
  # function to replace levels that indicate a missing value with NA
  if(x == "Too many missings" | x == "Q1 not returned"){
    x <- NA
  }
  x
}
as.numeric.factor <- function(x){
  #function for turning factor values (levels) into numeric values
  as.numeric(levels(x))[x]
}
bigFive <- bigFive_raw %>% 
  select(pident, aneurot, aextrave, aopenes, aagreeab, aconscie) %>% 
  apply(2,function(x){sapply(x, toNA)}) %>% 
  as.data.frame() %>% 
  # filter(rowSums(is.na(.)) != ncol(.)-1) %>% #delete all rows with no BigFive data. 
  mutate(pident = as.factor(pident)) %>% 
  mutate_at(vars(-pident), ~sapply(.x, as.numeric))

#IDS
IDS_raw <- read.spss("psychological_and_lifestyle_variables/W1/IDS/N1_A235D.sav", to.data.frame = T, use.missings = F)
as.character.factor <- function(x){
  #function for turning factor values (levels) into character values
  as.character(levels(x))[x]
}
removeAddition <- function(x){
  # turn missings into NA
  # remove additions and keep score
  if(x %in% c("No scale: too many missings","Q1 not returned")){
    x <- NA
  } else if(str_count(x) > 2){
    x <- str_trunc(x, 2, "right", "")
  }
  x
}
IDS <- IDS_raw %>% 
  select(pident, aids, aidssev) %>% 
  mutate(aids = as.character.factor(aids)) %>% 
  mutate(aids = sapply(aids, removeAddition)) %>% 
  mutate(aidssev = na_if(aidssev, "Q1 not returned")) %>% 
  mutate(aidssev = na_if(aidssev, "No scale: too many missings")) %>% 
  mutate(aids = as.numeric(aids))

#merge into dataframe
whole_set <- left_join(w3Chronicity, demographic, by = "pident") %>% 
  left_join(bigFive, by = "pident") %>% 
  left_join(IDS, by = "pident") %>% 
  mutate_if(is.factor, droplevels) %>% 
  mutate(aidssev = as.numeric(aidssev)) #turn ordinal data into numerical data

options(na.action = "na.pass")
modelMatrix <- model.matrix(Remitted_depression ~ ., data = whole_set[,-1])[,-1] #transform categorical nominal features into numeric values
whole_set <- cbind(pident = whole_set$pident, Remitted_depression = whole_set$Remitted_depression, as.data.frame(modelMatrix))

whole_set$Remitted_depression <- revalue(whole_set$Remitted_depression, c("remitted"="remitted", "non-remitted"="not_remitted")) 
table(whole_set$Remitted_depression) # 407 (58%) remitted, 397 (42%) non-remitted

# write to file for use in human vs machine
# write_csv(whole_set, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/human_vs_machine/data/whole_set_clinical.csv")

##############################################################################################################
## Using heldout test set 
##############################################################################################################
## partition data
set.seed(42)
trainIndex <- createDataPartition(y=whole_set$Remitted_depression, times=1, p=0.8, list=F)

# training set
segTrain <- as.data.frame(whole_set[trainIndex,])
prepTrain <- preProcess(segTrain[, 3:12],
                   #used preProcessing arguments only for bigFive and IDS data
                   method = c("medianImpute"))
transformedTrain <- predict(prepTrain, newdata = segTrain[, 3:12])
segTrain[3:12] <- transformedTrain

# testing set (preprocess similarly, i.e. use train imputations on testdata)
segTest <- as.data.frame(whole_set[-trainIndex,])
transformedTest <- predict(prepTrain, newdata = segTest[, 3:12])
segTest[3:12] <- transformedTest

summary(segTrain$Remitted_depression) #remitted/non-remitted = 326/318
summary(segTest$Remitted_depression) #remitted/non-remitted = 81/79

# Export train/test data for human vs machine
# write_csv(segTrain, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/human_vs_machine/data/clinical_train_imputed.csv")
# write_csv(segTest, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/human_vs_machine/data/clinical_test_imputed.csv")

# Remove pident for further analysis and export files for python based analysis on same set
segTrain <- segTrain[,-1]
segTest <- segTest[,-1]
# write_csv(segTrain[,-1], "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/VM/data/clinical_train_imputed.csv")
# write_csv(segTest[,-1], "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/VM/data/clinical_test_imputed.csv")

## Fit model
## set seed list for reproduction
set.seed(42)
seeds <- vector(mode = "list", length = 101)
for(i in 1:100) seeds[[i]] <- sample.int(100, 50)
seeds[[101]] <- sample.int(100,1)

train_ctrl <- trainControl(method="repeatedcv",
                           number = 10,
                           repeats = 10,
                           summaryFunction = twoClassSummary,
                           search = "random",
                           classProbs = TRUE,
                           seeds = seeds,
                           allowParallel = TRUE,
                           verboseIter = TRUE)

# tune_grid <- expand.grid(nrounds = seq(from = 200, to = nrounds, by = 50),
#                          eta = c(0.025, 0.05, 0.1, 0.3),
#                          max_depth = c(2, 3, 4, 5, 6),
#                          gamma = 0,
#                          colsample_bytree = 1,
#                          min_child_weight = 1,
#                          subsample = 1)

XGbTune <- train(x = segTrain[,-1],
                 y = segTrain[,1],
                 method = "xgbTree",
                 tuneLength = 1000,
                 #tuneGrid = grid,
                 metric = "ROC",
                 trControl = train_ctrl,
                 preProcess = NULL,
                 verbose = TRUE)

XGbTune
# saveRDS(XGbTune, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/XGbTune3.RDS")
XGbTune <- readRDS("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/XGbTune3.RDS")
XGbTune$finalModel$tuneValue

# # helper function for the plots
# tuneplot <- function(x, probs = .90) {
#   ggplot(x) +
#     coord_cartesian(ylim = c(quantile(x$results$RMSE, probs = probs), min(x$results$RMSE))) +
#     theme_bw()
# }
# 
# tuneplot(XGbTune)

## predictions on test set:
XGbPred <- predict(XGbTune, segTest[,-1])
cm <- confusionMatrix(XGbPred, segTest$Remitted_depression)
draw_confusion_matrix(cm, ratio = TRUE)

accuracies <- list()
for(i in seq(0.4, 0.6, 0.01)){
  pred_ <- factor(ifelse(predict(XGbTune, segTest[,-1], type="prob")$remitted>i,"remitted","not_remitted"), levels = c("remitted", "not_remitted"))
  cm_ <- confusionMatrix(pred_, segTest$Remitted_depression)
  accuracies <- append(accuracies, cm_$overall[1])
}
names(accuracies) <- seq(0.4, 0.6, 0.01)
opt_thresh <- as.numeric(names(accuracies[accuracies == max(sapply(accuracies, max))]))[1]
XGbPred <- factor(ifelse(predict(XGbTune, segTest[,-1], type="prob")$remitted>opt_thresh,"remitted","not_remitted"), levels = c("remitted", "not_remitted"))
cm <- confusionMatrix(XGbPred, factor(segTest$Remitted_depression, levels = c("remitted", "not_remitted")))
draw_confusion_matrix(cm, ratio = TRUE)
# export in pdf 7 * 6 inches


# variable importance
var_imp <- varImp(XGbTune, scale = FALSE)
var_imp$importance
analyte_names <- read_csv("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/analyte_names.csv")

analyte_importance <- var_imp$importance %>% 
  mutate(symbol = str_sub(rownames(var_imp$importance), start = 2)) %>% 
  left_join(y = analyte_names, by = c("symbol" = "Analyte abbreviation")) %>% 
  dplyr::rename(Gain = Overall) %>% 
  mutate(Gain = round(Gain, 5)) %>% 
  select(2,3,1)
# write_excel_csv(analyte_importance, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/clinical_only_analyte_importance_XGB-Gain.csv")


#Get predicted class probabilities so we can build ROC curve.
XGbProbs <- predict(XGbTune, segTest[,-1], type="prob")
head(XGbProbs)

#build ROC curve
# Setting levels: control = not_remitted, case = remitted --> same as in confusion matrix
# Setting direction: controls < cases
XGbROC <- roc(factor(segTest[,1], levels = c("not_remitted","remitted")), XGbProbs[,"remitted"])

# saveRDS(XGbROC, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/roc_2ychronicity_clinical.RDS")
auc(XGbROC)
# ci.auc(XGbROC, method = "bootstrap", boot.n = 10000, boot.stratified = TRUE)
# roc.test(XGbROC_clin, XGbROC_prot_clin, method = "bootstrap", boot.stratified = TRUE, boot.n = 10000)

#plot ROC curve
closest <- function(ROC, optimal_threshold) {
  # function for determining coordinates of specificity and sensitivitywith optimal cutoff for accuracy
  coord <- which(abs(ROC$thresholds-optimal_threshold)==min(abs(ROC$thresholds-optimal_threshold))) 
  spec <- ROC$specificities[coord]
  sens <- ROC$sensitivities[coord]
  output <- list(spec, sens)
  names(output) <- c("specificity", "sensitivity")
  return(output)
}
cutoff <- closest(XGbROC, opt_thresh)
plot_df <- data.frame(specificity = XGbROC$specificities, sensitivity = XGbROC$sensitivities)

plot <- ggplot(plot_df, aes(specificity, sensitivity)) +
  geom_path() +
  scale_x_reverse() +
  geom_abline(intercept = 1, slope = 1, colour='grey') +
  geom_point(aes(x = cutoff$specificity, y = cutoff$sensitivity), size = 5, color = "#F87474") +
  geom_text(data = subset(plot_df, specificity==cutoff$specificity & sensitivity ==cutoff$sensitivity), 
            aes(label = paste0("balanced \naccuracy = ", as.character(round(cm$byClass[11],2))), 
                hjust = cutoff$specificity-0.7, 
                vjust = cutoff$sensitivity+0.7), 
            colour = "black", 
            size=6) +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
  theme_classic() +
  theme(text = element_text(size=18)) +
  labs(title = paste0("AUROC =", signif(auc(XGbROC), 2))) +
  theme(plot.title = element_text(size=16))
plot

# ggsave("/Users/philippehabets/Dropbox/STRESS_INDEX/manuscripts/prediction_2ychronicity/Figures/output/auc_clinical.pdf",
#        plot,
#        width = 15, height = 15, units = "cm")


#SHAP (Shapley Additive Explanation) of variables
########  #SHAP (Shapley Additive Explanation) of variables ########
# Top 10 in plot with buildin plot function of xgboost package
xgb.plot.shap(data = as.matrix(segTrain[XGbTune$finalModel$feature_names]), model = XGbTune$finalModel, top_n = 10)

# With SHAP for xgboost package:
dataX <- as.matrix(segTrain[,XGbTune$finalModel$feature_names])
shap_values <- shap.values(xgb_model = XGbTune$finalModel, X_train = dataX)
shap_values$mean_shap_score

shap.plot.summary.wrap1(XGbTune$finalModel,dataX, top_n=50, dilute = F)

strip_a <- function(x) {
  #function to strip variable names of "a" if not Age or Sex
  if(str_starts(x, "a")) {
    #if string starts with "a",
    if(!str_starts(x, "Ag")){
      #if string starts nog with "Ag" (age)
      x <- str_sub(x, start = 2)
    }
  } else {
    #if string does not start with "a", it is either 'Age' or "Sex" variable: keep all letters
    x <- x
  }
  return(x)
}

analyte_shap <- data.frame(symbol = sapply(names(shap_values$mean_shap_score), strip_a), 
                           mean_shap_score =  shap_values$mean_shap_score) %>% 
  mutate(symbol = ifelse(symbol == "Sexefemale", "Sex", symbol)) %>% 
  mutate_all(~tolower(.x))

# write_excel_csv(analyte_shap, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/shap_xgb_clinical.csv")



############################################################################################################
# following example on: https://liuyanguu.github.io/post/2018/10/14/shap-visualization-for-xgboost/
############################################################################################################
dataX <- as.matrix(segTrain[,-1])
param_list <- list(objective = "binary:logistic",  
                   eta = 0.06012163,
                   max_depth = 1,
                   gamma = 2.26149,
                   subsample = 0.8548254)

labels <- data.frame(label = segTrain[,1]) %>% 
  mutate(label = ifelse(label == "not_remitted", 1, 0))

mod <- xgboost::xgboost(data = dataX, 
                        label = as.matrix(labels), 
                        params = param_list, nrounds = 897,
                        verbose = TRUE, 
                        nthread = 10,
                        early_stopping_rounds = 8)

shap_values <- shap.values(xgb_model = mod, X_train = dataX)
shap_values$mean_shap_score

analyte_shap <- data.frame(symbol = str_sub(names(shap_values$mean_shap_score), start = 3), 
                           mean_shap_score =  shap_values$mean_shap_score) %>% 
  left_join(y = analyte_names, by = c("symbol" = "Analyte abbreviation")) %>% 
  select(1,3,2)
#write_excel_csv(analyte_shap, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/analyte_shap_xgb.csv")












