library(tidyverse)
library(plyr)
library(foreign)
library(caret)
library(xgboost)
library(tictoc)
library(SHAPforxgboost)
library(doMC)
library(pROC)

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
  text(30, 35, names(cm$overall[1]), cex=1.5, font=2)
  text(30, 20, round(as.numeric(cm$overall[1]), 3), cex=1.4)
  text(70, 35, names(cm$overall[2]), cex=1.5, font=2)
  text(70, 20, round(as.numeric(cm$overall[2]), 3), cex=1.4)
}
###############

##create dataframe: 1170 individual (join of 1860 individuals with defined outcome and 1837 individuals with proteomics measurements)
w3Chronicity <- read.spss("/Users/philippehabets/Dropbox/STRESS_INDEX/data/outcome_groups/DSM_groups.sav", to.data.frame = TRUE, use.missings = FALSE)
w3Chronicity <- w3Chronicity[,c(1,3)]
prtm <- read_csv("/Users/philippehabets/Dropbox/STRESS_INDEX/data/blood_and_saliva_variables/W1/proteomics/output/proteomics_replaced_outliers.csv")
prtm$Pident <- as.factor(prtm$Pident) 

whole_set <- inner_join(w3Chronicity, prtm, by = c("pident" = "Pident"))
whole_set$Remitted_depression <- revalue(whole_set$Remitted_depression, c("remitted"="remitted", "non-remitted"="not_remitted"))
table(whole_set$Remitted_depression) # 312 (51%) remitted, 299 (49%) non-remitted

##############################################################################################################
## Using heldout test set 
##############################################################################################################
## partition data
set.seed(42)

trainIndex <- createDataPartition(y=whole_set$Remitted_depression, times=1, p=0.8, list=F)

## export test set for sample selection for human rating 
# exportTest <-as.data.frame(whole_set[-trainIndex, c(1:2)])
# write_csv(exportTest, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/human_vs_machine/data/additional_test.csv")

segTrain <- as.data.frame(whole_set[trainIndex,-c(1,3)])
segTrain[,1] <- as.factor(segTrain[,1])
segTrain <- segTrain %>% 
  #log10 transform as preprocessing step (not included as option in caret, so done manually here)
  mutate_at(vars(-1), log10)

segTest <- as.data.frame(whole_set[-trainIndex,-c(1,3)])
segTest[,1] <- as.factor(segTest[,1])
segTest <- segTest %>% 
  #log10 transform as preprocessing step (not included as option in caret, so done manually here)
  mutate_at(vars(-1), log10)

summary(segTrain$Remitted_depression) #remitted/non-remitted = 250/240
summary(segTest$Remitted_depression) #remitted/non-remitted = 62/59

## Fit model
## set seed list for reproduction
#function to time progress:
toc_hour <- function() {
  x <- toc()
  hours <- (x$toc - x$tic) / 3600
  cat(paste(c("or ", as.character(round(hours, digits = 2)), " hours", "\n"),collapse = ""))
}
hyperparams <- 1000

set.seed(42)
seeds <- vector(mode = "list", length = 101)
for(i in 1:100) seeds[[i]] <- sample.int(1000, hyperparams)
seeds[[101]] <- sample.int(1000,1)

# train_ctrl <- trainControl(method="repeatedcv",
#                            workers = 10,
#                            number = 10,
#                            repeats = 10,
#                            summaryFunction = twoClassSummary,
#                            search = "random",
#                            classProbs = TRUE,
#                            preProcOptions = list(knnSummary = median, k=5), 
#                            seeds = seeds,
#                            allowParallel = TRUE,
#                            verboseIter = TRUE)

adaptControl <- trainControl(method = "adaptive_cv",
                             number = 10, repeats = 10,
                             adaptive = list(min = 10, alpha = 0.05, method = "gls", complete = TRUE), #using adaptive resampling to have more efficient computation, see http://topepo.github.io/caret/adaptive-resampling.html
                             classProbs = TRUE,
                             summaryFunction = twoClassSummary,
                             preProcOptions = list(knnSummary = median, k=5), 
                             search = "random",
                             allowParallel = TRUE,
                             verboseIter = TRUE)

# tune_grid <- expand.grid(nrounds = seq(from = 200, to = nrounds, by = 50),
#                          eta = c(0.025, 0.05, 0.1, 0.3),
#                          max_depth = c(2, 3, 4, 5, 6),
#                          gamma = 0,
#                          colsample_bytree = 1,
#                          min_child_weight = 1,
#                          subsample = 1)
tic()
XGbTune <- train(x = segTrain[,-1],
                 y = segTrain[,1],
                 method = "xgbTree",
                 tuneLength = hyperparams,
                 #tuneGrid = grid,
                 metric = "ROC",
                 trControl = adaptControl, #train_ctrl, 
                 preProcess = c("zv","knnImpute"),
                 verbose = TRUE)
toc_hour()

# XGbTune
# saveRDS(XGbTune, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/output/XGbTune.RDS")
XGbTune <- readRDS("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/output/XGbTune.RDS")
XGbTune$finalModel
plot(density(XGbTune$resample$ROC))

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
cm <- confusionMatrix(XGbPred, segTest[,1])
draw_confusion_matrix(cm, ratio = TRUE)
var_imp <- varImp(XGbTune, scale = FALSE)
var_imp$importance
analyte_names <- read_csv("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/analyte_names.csv")

analyte_importance <- var_imp$importance %>% 
  mutate(symbol = str_sub(rownames(var_imp$importance), start = 3)) %>% 
  left_join(y = analyte_names, by = c("symbol" = "Analyte abbreviation")) %>% 
  dplyr::rename(Gain = Overall) %>% 
  mutate(Gain = round(Gain, 5)) %>% 
  select(2,3,1)
# write_excel_csv(analyte_importance, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/analyte_importance_XGB-Gain.csv")


#Get predicted class probabilities so we can build ROC curve.
XGbProbs <- predict(XGbTune, segTest[,-1], type="prob")
head(XGbProbs)

#build ROC curve
XGbROC <- roc(segTest[,1], XGbProbs[,"remitted"])
auc(XGbROC)

#plot ROC curve
ggplot(data.frame(specificity = XGbROC$specificities, sensitivity = XGbROC$sensitivities), aes(specificity, sensitivity)) +
  geom_path() +
  scale_x_reverse() +
  geom_abline(intercept = 1, slope = 1, colour='grey') +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
  theme_classic() +
  labs(title = paste0("AUROC =", signif(auc(XGbROC), 3)))


#SHAP (Shapley Additive Explanation) of variables
library(reticulate)
Sys.setenv(RETICULATE_MINICONDA_PATH = "/opt/miniconda3/")
use_miniconda("shap_installed")

library(DALEX)
exp_xgb <- explain(XGbTune, data = segTrain[,-1], y = as.numeric(segTrain[,1])-1)

library(shapper)
# shapper::install_shap()
ive_xgb <- shap(exp_xgb, new_observation = segTest[1,-1])
ive_xgb


############################################################################################################
# following example on: https://liuyanguu.github.io/post/2018/10/14/shap-visualization-for-xgboost/
############################################################################################################
dataX <- as.matrix(segTrain[,-1])
param_list <- list(objective = "binary:logistic",  
                   eta = 0.186079,
                   max_depth = 1,
                   gamma = 1.38372,
                   subsample = 0.691072,
                   colsample_bytree = 0.6961078,
                   min_child_weight = 9)

labels <- data.frame(label = segTrain[,1]) %>% 
  mutate(label = ifelse(label == "not_remitted", 1, 0))

mod <- xgboost::xgboost(data = dataX, 
                        label = as.matrix(labels), 
                        params = param_list, nrounds = 878,
                        verbose = TRUE, 
                        nthread = 10,
                        early_stopping_rounds = 8)

shap_values <- shap.values(xgb_model = mod, X_train = dataX)
shap_values$mean_shap_score

analyte_shap <- data.frame(symbol = str_sub(names(shap_values$mean_shap_score), start = 3), 
                           mean_shap_score =  shap_values$mean_shap_score) %>% 
  left_join(y = analyte_names, by = c("symbol" = "Analyte abbreviation")) %>% 
  select(1,3,2)
# write_excel_csv(analyte_shap, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/analyte_shap_xgb.csv")

shap_values_prtm <- shap_values$shap_score
shap_long_prtm <- shap.prep(shap_contrib = shap_values_prtm, X_train = dataX)
shap.plot.summary(shap_long_prtm)
shap.plot.summary(shap_long_prtm, x_bound = 1.5, dilute = 10)
shap.plot.summary.wrap1(model = mod, dataX, top_n = 50, dilute = FALSE)

shap.plot.dependence(data_long = shap_long_prtm, x="apBDNF",
                     y = "apCortisol", color_feature = "apCortisol")

# **SHAP force plot**
plot_data <- shap.prep.stack.data(shap_contrib = shap_values_prtm,
                                  top_n = 10,
                                  n_groups = 4)
shap.plot.force_plot(plot_data)
shap.plot.force_plot(plot_data, zoom_in_group = 2)

# plot all the clusters:
shap.plot.force_plot_bygroup(plot_data)

############################################################################################################
############################################################################################################























