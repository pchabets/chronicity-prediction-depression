library(tidyverse)
library(caret)
library(kernlab)
library(doMC)
library(pROC)
library(tictoc)
setwd("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/")

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

# Load train and test set (same as in XGBoost script)
segTrain <- read_csv("2-yearChronicity/output/segTrain_bagImputedPreProcessingProteomicsBatchCorrected_AND_clinical.csv") %>%
  mutate(Remitted_depression = factor(Remitted_depression, levels=c("remitted", "not_remitted"))) %>% as.data.frame()

segTest <- read_csv("2-yearChronicity/output/segTest_bagImputedPreProcessingProteomicsBatchCorrected_AND_clinical.csv") %>%
  mutate(Remitted_depression = factor(Remitted_depression, levels=c("remitted", "not_remitted"))) %>% as.data.frame()

## check if levels have same order in train and test set:
levels(segTrain$Remitted_depression) == levels(segTest$Remitted_depression)

#################  fit model  ###################
toc_hour <- function() {
  # function to keep track of time passed
  x <- toc()
  hours <- (x$toc - x$tic) / 3600
  print(paste(c("or ", as.character(round(hours, digits = 2)), " hours"),collapse = ""))
}
hyperparams <- 1000 #number of hyperparameter combinations to use in tuning

# ## set seed list for reproduction
# set.seed(42)
# seeds <- vector(mode = "list", length = 101)
# for(i in 1:100) seeds[[i]] <- sample.int(1000, hyperparams)
# seeds[[101]] <- sample.int(1000,1)

adaptControl <- trainControl(method = "adaptive_cv",
                             number = 10, repeats = 10,
                             adaptive = list(min = 10, alpha = 0.05, method = "gls", complete = TRUE), #using adaptive resampling to have more efficient computation, see http://topepo.github.io/caret/adaptive-resampling.html
                             classProbs = TRUE,
                             summaryFunction = twoClassSummary,
                             search = "random", 
                             allowParallel = TRUE,
                             verboseIter = TRUE)

tic()
svmTune <- train(x = segTrain[,-1],
                 y = segTrain[,1],
                 method = "svmLinear",
                 tuneLength = hyperparams,
                 # tuneGrid = tune_grid,
                 metric = "ROC",
                 trControl = adaptControl,
                 preProcess = list("medianImpute" = colnames(segTrain)[2:11],
                                   "center" = colnames(segTrain)[-1], 
                                   "scale" = colnames(segTrain)[-1]), # will only be performed on clinical variables as proteomic variables are already imputed, center + scale for svm
                 verbose = TRUE)

toc_hour()

# saveRDS(svmTune, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/svmLinearTune2.RDS")
svmTune <- readRDS("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/svmLinearTune2.RDS")


## predictions on test set:
svmPred <- predict(svmTune, segTest[,-1])

## With alternative threshold: 
accuracies <- list()
for(i in seq(0.4, 0.6, 0.01)){
  pred_ <- factor(ifelse(predict(svmTune, segTest[,-1], type="prob")$remitted>i,"remitted","not_remitted"), levels = c("remitted", "not_remitted"))
  cm_ <- confusionMatrix(pred_, segTest[,1])
  accuracies <- append(accuracies, cm_$overall[1])
}
names(accuracies) <- seq(0.4, 0.6, 0.01)
opt_thresh <- as.numeric(names(accuracies[accuracies == max(sapply(accuracies, max))]))
svmPred <- factor(ifelse(predict(svmTune, segTest[,-1], type="prob")$remitted>opt_thresh,"remitted","not_remitted"), levels = c("remitted", "not_remitted"))

# Confustion matrix
cm <- confusionMatrix(svmPred, segTest[,1])
draw_confusion_matrix(cm, ratio = TRUE)

var_imp <- varImp(svmTune, scale = FALSE)
var_imp$importance

#Get predicted class probabilities so we can build ROC curve.
svmProbs <- predict(svmTune, segTest[,-1], type="prob")
head(svmProbs)

#build ROC curve
svmROC <- roc(segTest[,1], svmProbs[,"remitted"])
# saveRDS(svmROC, file = "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/roc_2ychronicity_prot_clin_SVM-L.RDS")
auc(svmROC)
ci.auc(svmROC)

#plot ROC curve
ggplot(data.frame(specificity = svmROC$specificities, sensitivity = svmROC$sensitivities), aes(specificity, sensitivity)) +
  geom_path() +
  scale_x_reverse() +
  geom_abline(intercept = 1, slope = 1, colour='grey') +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
  theme_classic() +
  labs(title = paste0("AUROC =", signif(auc(svmROC), 3)))

