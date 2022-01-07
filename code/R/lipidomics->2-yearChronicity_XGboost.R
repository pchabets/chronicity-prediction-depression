library(tidyverse)
library(plyr)
library(foreign)
library(caret)
library(xgboost)
library(doMC)
library(pROC)
library(tictoc)

setwd("/Users/philippehabets/Dropbox/STRESS_INDEX/")
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
w3Chronicity <- read.spss("data/outcome_groups/DSM_groups.sav", to.data.frame = TRUE, use.missings = FALSE)
w3Chronicity <- w3Chronicity[,c(1,3)]

#metabolomics (= lipidomics):
#define some needed functions:
as.character.factor <- function(x) {as.character(levels(x))[x]} #function for turning factor values (levels) into character values
replace <- function(x){if_else(x =="Very low", "-1", if_else(x=="No result", "-2", x))} #replace "Very low" with "-1", replace "No result" with "-2"

mtbl <- read.csv("data/blood_and_saliva_variables/W1/metabolomics/output/mtbl_measurements.csv")
mtbl <- mtbl %>%
  mutate_if(is.factor, as.character.factor) %>% #all columns as character to make mutate_all with replace-function possible
  mutate_if(is.character, replace) %>% #replace all "Very low" and "No result" values with "-1" and "-2" respectively
  mutate_all(as.numeric) %>%  #turn all metabolites into numeric for imputations
  mutate(pident = as.factor(pident))

whole_set <- inner_join(w3Chronicity, mtbl, by = "pident")
whole_set$Remitted_depression <- revalue(whole_set$Remitted_depression, c("remitted"="remitted", "non-remitted"="not_remitted")) #otherwise an error will occur, as non-remitted is not a valeid R variable name.
table(whole_set$Remitted_depression) # 398 (50,4%) remitted, 392 (49,6%) non-remitted

##############################################################################################################
## Using heldout test set 
##############################################################################################################
## partition data
set.seed(42)
trainIndex <- createDataPartition(y=whole_set$Remitted_depression, times=1, p=0.8, list=F)

segTrain <- as.data.frame(whole_set[trainIndex,-1])
table(segTrain$Remitted_depression) # 319 vs 314

segTest <- as.data.frame(whole_set[-trainIndex,-1])
table(segTest$Remitted_depression) # 79 vs 78

##############################################################################################################
## Preprocessing
##############################################################################################################
minEx <- function(x){
  #function to get minimum of feature without taking -1 and -2 into account
  min(x[!x %in% c(-1,-2)], na.rm=TRUE)
}
trainLipidMins <- sapply(segTrain[names(mtbl)[-1]], minEx) #store min values of all lipid variables in trainset

# Some needed functions for setting up replacement function applicable to both train and test set, but only informed on min values of train set
is_any_low <- function(x){any(x == -1, na.rm = TRUE)} #function to check if any value is -1, meaning it is "Very low"
is_any_NR <- function(x){any(x == -2, na.rm = TRUE)} #function to check if any value is -2, meaning it has "No result"
replaceMinCol <- function(x){
  #check for each column if it contains -1. If so, execute replaceMinVal() function on this column
  for(i in c(1:ncol(x))){
    if(is_any_low(x[,i])) {
      for (j in c(1:length(x[,i]))) {
        if(x[j,i] == -1) {
          x[j,i] <- trainLipidMins[i] 
        }
        x[,i] <- x[,i]
      }
    }
    x <- x
  }
  return(x)
}

## Train data:
# Replace -1 with minimum values of that analyte
segTrain[names(trainLipidMins)] <- replaceMinCol(segTrain[names(trainLipidMins)])
# Replace -2 with NA
segTrain[names(trainLipidMins)] <- segTrain[names(trainLipidMins)] %>% 
  mutate_if(is_any_NR, ~na_if(.x, -2)) 

# Check if any lipidomic variable has more than 10% missings
sum(colMeans(is.na(segTrain))>0.1) # nope

#impute lipidomic variables
pp <- preProcess(segTrain[names(trainLipidMins)], method = c("zv", "knnImpute"), k = 5, knnSummary = median)
transformedTrain <- predict(pp, newdata = segTrain[names(trainLipidMins)])
segTrain[names(trainLipidMins)] <- transformedTrain

## Test data:
# Replace -1 with minimum values of that analyte
segTest[names(trainLipidMins)] <- replaceMinCol(segTest[names(trainLipidMins)])
# Replace -2 with NA
segTest[names(trainLipidMins)] <- segTest[names(trainLipidMins)] %>% 
  mutate_if(is_any_NR, ~na_if(.x, -2)) 

#imputetest data based on traindata fit of preprocessing 
transformedTest <- predict(pp, newdata = segTest[names(trainLipidMins)])
segTest[names(trainLipidMins)] <- transformedTest #lipidomics

##############################################################################################################
## Fit model 
##############################################################################################################
toc_hour <- function() {
  # function to keep track of time passed
  x <- toc()
  hours <- (x$toc - x$tic) / 3600
  print(paste(c("or ", as.character(round(hours, digits = 2)), " hours"),collapse = ""))
}
hyperparams <- 1000 #number of hyperparameter combinations to use in tuning

## set seed list for reproduction
set.seed(42)
seeds <- vector(mode = "list", length = 101)
for(i in 1:100) seeds[[i]] <- sample.int(1000, hyperparams)
seeds[[101]] <- sample.int(1000,1)

adaptControl <- trainControl(method = "adaptive_cv",
                             number = 10, repeats = 10,
                             adaptive = list(min = 10, alpha = 0.05, method = "gls", complete = TRUE), #using adaptive resampling to have more efficient computation, see http://topepo.github.io/caret/adaptive-resampling.html
                             classProbs = TRUE,
                             summaryFunction = twoClassSummary,
                             search = "random", 
                             allowParallel = TRUE,
                             preProcOptions = list(thresh = 0.99),
                             # seeds = seeds,
                             verboseIter = TRUE)
tic()
xgbTune <- train(x = segTrain[,-1],
                 y = segTrain[,1],
                 method = "xgbTree",
                 tuneLength = hyperparams,
                 # tuneGrid = tune_grid,
                 metric = "ROC",
                 trControl = adaptControl,
                 preProcess = list(pca = names(trainLipidMins)),
                 verbose = TRUE)
toc_hour()

xgbTune
# saveRDS(xgbTune, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/enTune4_v2.RDS")
xgbTune <- readRDS("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/enTune4_v2.RDS")

xgbTune$finalModel
plot(xgbTune, metric = "ROC", scales = list(x = list(log =2)))
roc_imp <- varImp(enTune, scale = FALSE) # ROC per variable
roc_imp

## predictions on test set:
xgbPred <- predict(xgbTune, segTest[,-1])
cm <- confusionMatrix(xgbPred, segTest[,1])
draw_confusion_matrix(cm, ratio = TRUE)

#Get predicted class probabilities so we can build ROC curve.
xgbProbs <- predict(xgbTune, segTest[,-1], type="prob")
head(xgbProbs)

#build ROC curve
xgbROC <- roc(segTest[,1], xgbProbs[,"remitted"])
auc(xgbROC)

#plot ROC curve
plot(enROC, type = "S")
ggplot(data.frame(specificity = xgbROC$specificities, sensitivity = xgbROC$sensitivities), aes(specificity, sensitivity)) +
  geom_path() +
  scale_x_reverse() +
  geom_abline(intercept = 1, slope = 1, colour='grey') +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
  theme_classic() +
  labs(title = paste0("AUROC =", signif(auc(xgbROC), 3)))




















