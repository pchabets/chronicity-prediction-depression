library(tidyverse)
library(caret)
library(glmnet)
library(tictoc)
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
  text(10, 35, names(cm$overall[1]), cex=1.5, font=2)
  text(10, 20, round(as.numeric(cm$overall[1]), 3), cex=1.4)
  text(50, 35, names(cm$byClass[11]), cex=1.5, font=2)
  text(50, 20, round(as.numeric(cm$byClass[11]), 3), cex=1.4)
  text(90, 35, names(cm$overall[2]), cex=1.5, font=2)
  text(90, 20, round(as.numeric(cm$overall[2]), 3), cex=1.4)
}
###############

## Load in train and test data as preprocessed in XGBoost script
segTrain <- read_rds("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/output/segTrain_bagImputedPreProcessingProteomicsBatchCorrected.RDS")
segTrain <- segTrain %>% select(-c(pident, applate))

segTest <- read_rds("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/output/segTest_bagImputedPreProcessingProteomicsBatchCorrected.RDS")
segTest <- segTest %>% select(-c(pident, applate))

##############################################################################################################
## Model training and validation
##############################################################################################################
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

adaptControl <- trainControl(method = "adaptive_cv",
                             number = 10, repeats = 10,
                             adaptive = list(min = 10, alpha = 0.05, method = "gls", complete = TRUE), #using adaptive resampling to have more efficient computation, see http://topepo.github.io/caret/adaptive-resampling.html
                             classProbs = TRUE,
                             summaryFunction = twoClassSummary,
                             preProcOptions = list(), 
                             search = "random",
                             # seeds = seeds,
                             allowParallel = TRUE,
                             verboseIter = TRUE)

tic()
elTune <- train(x = segTrain[,-1],
                 y = segTrain[,1],
                 method = "glmnet",
                 tuneLength = hyperparams,
                 #tuneGrid = grid,
                 metric = "ROC",
                 trControl = adaptControl, #train_ctrl, 
                 preProcess = c("center", "scale"),
                 verbose = TRUE)
toc_hour()

# saveRDS(elTune, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/output/elTune_batchCorrected.RDS")
elTune <- readRDS("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/output/elTune_batchCorrected.RDS")

## predictions on test set:
elPred <- predict(elTune, segTest[,-1])
cm <- confusionMatrix(elPred, segTest[,1])
draw_confusion_matrix(cm, ratio = TRUE)

## Variable importance (coefficient sorting)
coef <- coef(elTune$finalModel, elTune$bestTune$lambda)
coef_abs <- data.frame(analyte = rownames(coef), coefficient = coef[,1])
coef_abs <- coef_abs %>%
  dplyr::slice(-1) %>% 
  arrange(desc(abs(coefficient)))
coef_abs <- rbind(data.frame(analyte = "intercept", coefficient = coef[1,]), coef_abs) %>% 
  mutate(coefficient = round(coefficient, 3))

analyte_names <- read_csv("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/analyte_names.csv")
analyte_importance <- coef_abs %>% 
  mutate(analyte = str_sub(analyte, start = 3)) %>% 
  left_join(y = analyte_names, by = c("analyte" = "Analyte abbreviation")) %>% 
  select(1,3,2)
# write_excel_csv(analyte_importance, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/analyte_importance_ELASTICNET.csv")

#Get predicted class probabilities so we can build ROC curve.
elProbs <- predict(elTune, segTest[,-1], type="prob")
head(elProbs)

#build ROC curve
elROC <- roc(segTest[,1], elProbs[,"remitted"])
# saveRDS(elROC, file = "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/roc_2ychronicity_proteomics_elnet.RDS")
auc(elROC)

#plot ROC curve
library(data.table)
ggplot(data.table("1-FPR" = elROC$specificities, TPR = elROC$sensitivities), aes(`1-FPR`, TPR)) +
  geom_path() +
  scale_x_reverse() +
  geom_abline(intercept = 1, slope = 1, colour='grey') +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
  theme_classic() +
  labs(title = paste0("AUROC =", signif(auc(elROC), 3)))








