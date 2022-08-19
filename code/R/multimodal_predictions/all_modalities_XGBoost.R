library(tidyverse)
library(plyr)
library(foreign)
library(doMC)
library(caret)
library(xgboost)
library(pROC)
library(tictoc)
setwd("/Users/philippehabets/Dropbox/STRESS_INDEX/")

segTrain <- read_csv("data/multimodal/wave1_omics_clinical_segTrain_final_selection.csv") %>% 
  select(-c(pident, applate))
segTest <- read_csv("data/multimodal/wave1_omics_clinical_segTest_final_selection.csv") %>% 
  select(-c(pident, applate))

segTrain$Remitted_depression <- as.factor(segTrain$Remitted_depression)
segTest$Remitted_depression <- as.factor(segTest$Remitted_depression)

summary(segTrain$Remitted_depression) #remitted/non-remitted = 209/196
summary(segTest$Remitted_depression) #remitted/non-remitted = 52/49

#Initialize parallel processing
registerDoMC(detectCores()-2)
getDoParWorkers()

## Fit model
#set number of hyperparameter combinations to use in CV
hyperparams <- 1000

## set seed list for reproduction
set.seed(101)
seeds <- vector(mode = "list", length = 101)
for(i in 1:100) seeds[[i]] <- sample.int(1000, hyperparams)
seeds[[101]] <- sample.int(1000,1)

adaptControl <- trainControl(method = "adaptive_cv",
                             number = 10, repeats = 10,
                             adaptive = list(min = 10, alpha = 0.05, method = "gls", complete = TRUE), #using adaptive resampling to have more efficient computation, see http://topepo.github.io/caret/adaptive-resampling.html
                             classProbs = TRUE, 
                             summaryFunction = twoClassSummary,
                             # preProcOptions = list(freqCut = 90/10, uniqueCut=10),
                             search = "random",
                             # seeds = seeds,
                             allowParallel = TRUE,
                             verboseIter = TRUE)

toc_hour <- function() {
  # function to keep track of time passed
  x <- toc()
  hours <- (x$toc - x$tic) / 3600
  print(paste(c("or ", as.character(round(hours, digits = 2)), " hours"),collapse = ""))
}

##train model
tic() 
xgbTune <- train(x = segTrain[,-1],
                 y = segTrain$Remitted_depression,
                 method = "xgbTree",
                 tuneLength = hyperparams,
                 #tuneGrid = grid,
                 metric = "ROC",
                 trControl = adaptControl, #train_ctrl, 
                 preProcess = c("zv"),
                 verbose = TRUE)
toc_hour()

xgbTune
# saveRDS(xgbTune, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/output/XGbTune_allFeatures_final.RDS")
xgbTune <- readRDS("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/output/XGbTune_allFeatures_final.RDS")
xgbTune$finalModel
plot(density(xgbTune$resample$ROC))

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

## predictions on test set:
XGbPred <- predict(xgbTune, segTest[,-1])
cm <- confusionMatrix(XGbPred, segTest$Remitted_depression)
draw_confusion_matrix(cm, ratio = TRUE)

#Variable importance
var_imp <- varImp(xgbTune, scale = FALSE)
var_imp$importance

#Get predicted class probabilities so we can build ROC curve.
XGbProbs <- predict(xgbTune, segTest[,-1], type="prob")
head(XGbProbs)

#build ROC curve
XGbROC <- roc(segTest$Remitted_depression, XGbProbs[,"not_remitted"])
auc(XGbROC)

#plot ROC curve
ggplot(data.frame(specificity = XGbROC$specificities, sensitivity = XGbROC$sensitivities), aes(specificity, sensitivity)) +
  geom_path() +
  scale_x_reverse() +
  geom_abline(intercept = 1, slope = 1, colour='grey') +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
  theme_classic() +
  labs(title = paste0("AUROC =", signif(auc(XGbROC), 3)))

# test significance of difference with best modal (proteomics + clinical)
xgb_prot_clin <- readRDS("scripts/Predictions_explorative/R.scripts/data/roc_2ychronicity_clinical.RDS")
roc.test(xgb_prot_clin, XGbROC, method='delong')

########  #SHAP (Shapley Additive Explanation) of variables ########
# Top 10 in plot with buildin plot function of xgboost package
xgb.plot.shap(data = as.matrix(segTrain[xgbTune$finalModel$feature_names]), model = xgbTune$finalModel, top_n = 10)

# With SHAP for xgboost package:
library(SHAPforxgboost)
dataX <- as.matrix(segTrain[,xgbTune$finalModel$feature_names])
shap_values <- shap.values(xgb_model = xgbTune$finalModel, X_train = dataX)
shap_values$mean_shap_score

shap.plot.summary.wrap1(xgbTune$finalModel,dataX, top_n=50, dilute = F)

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
    #if string does not start with "a", it is either 'Age', "Sex", PRS, metabolite-PCA or transcriptomic variable: keep all letters
    x <- x
  }
  return(x)
}

analyte_names <- read_csv("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/analyte_names.csv")
analyte_shap <- data.frame(symbol = sapply(names(shap_values$mean_shap_score), strip_ap), 
                           mean_shap_score =  shap_values$mean_shap_score) %>% 
  left_join(y = analyte_names, by = c("symbol" = "Analyte abbreviation")) %>% 
  select(1,3,2)
# write_excel_csv(analyte_shap, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/shap_xgb_all_modalities.csv")









