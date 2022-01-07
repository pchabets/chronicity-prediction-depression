library(tidyverse)
library(plyr)
library(foreign)
library(doMC)
library(caret)
library(xgboost)
library(pROC)
library(tictoc)
setwd("/Users/philippehabets/Dropbox/STRESS_INDEX/")

# segTrain <- read_csv("scripts/VM/Python/output/transcriptomics_TopVariance_boruta_selection_95%_TRAIN_maxdepth1.csv")
segTrain <- read_csv("scripts/VM/Python/output/transcriptomics_TopVariance_boruta_selection_95%_TRAIN_maxdepth3.csv")
# segTrain <- read_csv("scripts/VM/Python/output/transcriptomics_variance_selection_top_250_TRAIN.csv")
# segTrain <- read_csv("scripts/VM/Python/output/transcriptomics_variance_selection_top_500_TRAIN.csv")
# segTrain <- read_csv("scripts/VM/Python/output/transcriptomics_boruta_selection_95%_TRAIN.csv")
# segTrain <- read_csv("scripts/VM/Python/output/transcriptomics_boruta_selection_98%_TRAIN.csv")
# segTrain <- read_csv("scripts/VM/Python/output/transcriptomics_boruta_selection_95%_TRAIN_maxdepth3.csv")
# segTrain <- read_csv("scripts/VM/Python/output/transcriptomics_boruta_selection_98%_TRAIN_maxdepth3.csv")
# segTrain <- read_csv("scripts/VM/Python/output/transcriptomics_boruta_selection_95%_TRAIN_maxdepth1.csv")
# segTrain <- read_csv("scripts/VM/Python/output/transcriptomics_boruta_selection_98%_TRAIN_maxdepth1.csv")

# segTest <- read_csv("scripts/VM/Python/output/transcriptomics_TopVariance_boruta_selection_95%_TEST_maxdepth1.csv")
segTest <- read_csv("scripts/VM/Python/output/transcriptomics_TopVariance_boruta_selection_95%_TEST_maxdepth3.csv")
# segTest <- read_csv("scripts/VM/Python/output/transcriptomics_variance_selection_top_250_TEST.csv")
# segTest <- read_csv("scripts/VM/Python/output/transcriptomics_variance_selection_top_500_TEST.csv")
# segTest <- read_csv("scripts/VM/Python/output/transcriptomics_boruta_selection_95%_TEST.csv")
# segTest <- read_csv("scripts/VM/Python/output/transcriptomics_boruta_selection_98%_TEST.csv")
# segTest <- read_csv("scripts/VM/Python/output/transcriptomics_boruta_selection_95%_TEST_maxdepth3.csv")
# segTest <- read_csv("scripts/VM/Python/output/transcriptomics_boruta_selection_98%_TEST_maxdepth3.csv")
# segTest <- read_csv("scripts/VM/Python/output/transcriptomics_boruta_selection_95%_TEST_maxdepth1.csv")
# segTest <- read_csv("scripts/VM/Python/output/transcriptomics_boruta_selection_98%_TEST_maxdepth1.csv")

segTrain <- as.data.frame(segTrain[,-c(1:2)])
segTest <- as.data.frame(segTest[,-c(1:2)])

segTrain$Remitted_depression <- as.factor(segTrain$Remitted_depression)
segTest$Remitted_depression <- as.factor(segTest$Remitted_depression)

summary(segTrain$Remitted_depression) #remitted/non-remitted = 276/260
summary(segTest$Remitted_depression) #remitted/non-remitted = 68/65

#Initialize parallel processing
registerDoMC(detectCores()-2)
getDoParWorkers()


## Fit model
#set number of hyperparameter combinations to use in CV
hyperparams <- 100

## set seed list for reproduction
set.seed(101)
seeds <- vector(mode = "list", length = 101)
for(i in 1:100) seeds[[i]] <- sample.int(1000, hyperparams)
seeds[[101]] <- sample.int(1000,1)

# train_ctrl <- trainControl(method="repeatedcv",
#                            number = 10,
#                            repeats = 10,
#                            summaryFunction = twoClassSummary,
#                            search = "random",
#                            classProbs = TRUE,
#                            preProcOptions = list(thresh = 0.95), 
#                            seeds = seeds,
#                            allowParallel = TRUE,
#                            verboseIter = TRUE)

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
                 y = segTrain[,1],
                 method = "xgbTree",
                 tuneLength = hyperparams,
                 #tuneGrid = grid,
                 metric = "ROC",
                 trControl = adaptControl, #train_ctrl, 
                 preProcess = c("zv"),
                 verbose = TRUE)
toc_hour()
print("need to save!")

# saveRDS(xgbTune, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/output/XGbTune_transcriptomics_95%.RDS")
xgbTune <- readRDS("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/output/XGbTune_transcriptomics_TopVariance_boruta_selection_95%_maxdepth3.RDS")
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
  text(30, 35, names(cm$overall[1]), cex=1.5, font=2)
  text(30, 20, round(as.numeric(cm$overall[1]), 3), cex=1.4)
  text(70, 35, names(cm$overall[2]), cex=1.5, font=2)
  text(70, 20, round(as.numeric(cm$overall[2]), 3), cex=1.4)
}
###############

## predictions on test set:
XGbPred <- predict(xgbTune, segTest[,-1])
cm <- confusionMatrix(XGbPred, segTest[,1])
draw_confusion_matrix(cm, ratio = TRUE)
var_imp <- varImp(xgbTune, scale = FALSE)
var_imp$importance

#Get predicted class probabilities so we can build ROC curve.
XGbProbs <- predict(xgbTune, segTest[,-1], type="prob")
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
  labs(title = paste0("AUROC =", signif(auc(XGbROC), 2)))








