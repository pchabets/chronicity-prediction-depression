library(tidyverse)
library(doMC)
library(caret)
library(xgboost)
library(pROC)
library(tictoc)

setwd("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/additional_analyses_revision/")

# Load data
segTrain <- read_csv("data/transcriptomics_postselection_boruta_95%_TRAIN_maxdepth3.csv") %>% select(-1)
segTest <- read_csv("data/transcriptomics_postselection_boruta_95%_TEST_maxdepth3.csv") %>% select(-1)

segTrain$Remitted_depression <- as.factor(segTrain$Remitted_depression)
segTest$Remitted_depression <- as.factor(segTest$Remitted_depression)

summary(segTrain$Remitted_depression) #remitted/non-remitted = 275/260
summary(segTest$Remitted_depression) #remitted/non-remitted = 69/65

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
                 y = segTrain %>% pull(Remitted_depression),
                 method = "xgbTree",
                 tuneLength = hyperparams,
                 #tuneGrid = grid,
                 metric = "ROC",
                 trControl = adaptControl, #train_ctrl, 
                 preProcess = c("zv"),
                 verbose = TRUE)
toc_hour()
# saveRDS(xgbTune, "data/XGbTune_transcriptomics.RDS")

xgbTune <- readRDS("data/XGbTune_transcriptomics.RDS")
xgbTune$finalModel$tuneValue
plot(density(xgbTune$resample$ROC))

#variable importance check
var_imp <- varImp(xgbTune, scale = FALSE)
var_imp$importance

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
cm <- confusionMatrix(XGbPred, segTest %>% pull(Remitted_depression))
draw_confusion_matrix(cm, ratio = TRUE)

accuracies <- list()
for(i in seq(0.35, 0.65, 0.0001)){
  pred_ <- factor(ifelse(predict(xgbTune, segTest[,-1], type="prob")$remitted>i,"remitted","not_remitted"), levels = c("not_remitted", "remitted"))
  cm_ <- confusionMatrix(pred_, segTest %>% pull(Remitted_depression))
  accuracies <- append(accuracies, cm_$byClass[11])
}
names(accuracies) <- seq(0.35, 0.65, 0.0001)
opt_thresh <- as.numeric(names(accuracies[accuracies == max(sapply(accuracies, max))]))[1]
XGbPred <- factor(ifelse(predict(xgbTune, segTest[,-1], type="prob")$remitted>opt_thresh,"remitted","not_remitted"), levels = c("not_remitted", "remitted"))
cm <- confusionMatrix(XGbPred, factor(segTest$Remitted_depression, levels = c("not_remitted", "remitted")))
draw_confusion_matrix(cm, ratio = TRUE)

#variable importance
var_imp <- varImp(xgbTune, scale = FALSE)
var_imp$importance

#Get predicted class probabilities so we can build ROC curve.
XGbProbs <- predict(xgbTune, segTest[,-1], type="prob")
head(XGbProbs)

#build ROC curve
XGbROC <- roc(segTest %>% pull(Remitted_depression), XGbProbs[,"remitted"])
auc(XGbROC)

#plot ROC curve
ggplot(data.frame(specificity = XGbROC$specificities, sensitivity = XGbROC$sensitivities), aes(specificity, sensitivity)) +
  geom_path() +
  scale_x_reverse() +
  geom_abline(intercept = 1, slope = 1, colour='grey') +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
  theme_classic() +
  labs(title = paste0("AUROC =", signif(auc(XGbROC), 2)))


#plot ROC curve
closest <- function(ROC, optimal_threshold) {
  # function for determining coordinates of specificity and sensitivity with optimal cutoff for accuracy
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
  geom_point(aes(x = cutoff$specificity, y = cutoff$sensitivity), size = 3, color = "#F87474") +
  geom_text(data = subset(plot_df, specificity==cutoff$specificity & sensitivity ==cutoff$sensitivity), aes(label = paste0("balanced \naccuracy =", as.character(round(cm$byClass[11],2))), hjust = cutoff$specificity-0.5, vjust = cutoff$sensitivity+0.7), colour = "black", size=5) +
  coord_fixed(ratio = 0.9, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
  theme_classic() +
  labs(title = paste0("AUROC =", signif(auc(XGbROC), 2))) +
  theme(text = element_text(size=16))
plot

ggsave("/Users/philippehabets/Dropbox/STRESS_INDEX/manuscripts/prediction_2ychronicity/BiolPsy/Figures/output/transcriptomics_kegg.pdf",
       plot,
       width = 23, height = 15, units = "cm")














