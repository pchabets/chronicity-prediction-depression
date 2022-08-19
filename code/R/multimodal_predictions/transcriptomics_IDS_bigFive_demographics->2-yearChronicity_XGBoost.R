library(tidyverse)
library(plyr)
library(foreign)
library(doMC)
library(caret)
library(xgboost)
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
  text(10, 35, names(cm$overall[1]), cex=1.5, font=2)
  text(10, 20, round(as.numeric(cm$overall[1]), 3), cex=1.4)
  text(50, 35, names(cm$byClass[11]), cex=1.5, font=2)
  text(50, 20, round(as.numeric(cm$byClass[11]), 3), cex=1.4)
  text(90, 35, names(cm$overall[2]), cex=1.5, font=2)
  text(90, 20, round(as.numeric(cm$overall[2]), 3), cex=1.4)
}
###############

##create dataframe
transcriptomicsTrain <- read_csv("scripts/VM/Python/output/transcriptomics_TopVariance_boruta_selection_95%_TRAIN_maxdepth3.csv")
transcriptomicsTest <- read_csv("scripts/VM/Python/output/transcriptomics_TopVariance_boruta_selection_95%_TEST_maxdepth3.csv")

transcriptomicsTrain <- as.data.frame(transcriptomicsTrain[,-1]) %>% 
  mutate(pident = as.factor(pident)) %>% 
  mutate(Remitted_depression = as.factor(Remitted_depression))
transcriptomicsTest <- as.data.frame(transcriptomicsTest[,-1]) %>% 
  mutate(pident = as.factor(pident)) %>% 
  mutate(Remitted_depression = as.factor(Remitted_depression))
  
#age, gender, years of education
demographic.raw <- read.spss("data/demographic_and_sociological_variables/W1/DOB_age_gender_nationality_and_education_of_respondents/N1_A100R.sav", to.data.frame = TRUE, use.missings = FALSE)
demographic <- demographic.raw %>% 
  mutate(pident = as.factor(pident)) %>% 
  select(pident, Sexe, Age, aedu)

#NEO-big five 
bigFive_raw <- read.spss("data/psychological_and_lifestyle_variables/W1/NEO-FFI_big_five_personality_test/N1_A240D.sav", to.data.frame = TRUE, use.missings = FALSE)
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
IDS_raw <- read.spss("data/psychological_and_lifestyle_variables/W1/IDS/N1_A235D.sav", to.data.frame = T, use.missings = F)
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

segTrain <- left_join(transcriptomicsTrain, demographic, by = "pident") %>% 
  left_join(bigFive, by = "pident") %>% 
  left_join(IDS, by = "pident") %>% 
  mutate_if(is.factor, droplevels) %>% 
  mutate(aidssev = as.numeric(aidssev)) #turn ordinal data into numerical data

segTest <- left_join(transcriptomicsTest, demographic, by = "pident") %>% 
  left_join(bigFive, by = "pident") %>% 
  left_join(IDS, by = "pident") %>% 
  mutate_if(is.factor, droplevels) %>% 
  mutate(aidssev = as.numeric(aidssev)) #turn ordinal data into numerical data

options(na.action = "na.pass")

modelMatrixTrain <- model.matrix(Remitted_depression ~ ., data = segTrain[,-1])[,-1] #transform categorical nominal features into numeric values
segTrain <- cbind(Remitted_depression = segTrain$Remitted_depression, as.data.frame(modelMatrixTrain))
table(segTrain$Remitted_depression) # 260 (49%) non-remitted, 276 (51%) remitted

modelMatrixTest <- model.matrix(Remitted_depression ~ ., data = segTest[,-1])[,-1] #transform categorical nominal features into numeric values
segTest <- cbind(Remitted_depression = segTest$Remitted_depression, as.data.frame(modelMatrixTest))
table(segTest$Remitted_depression) # 65 (49%) non-remitted, 68 (51%) remitted

which(sapply(segTrain, anyNA)) #no transcriptomics missing, only some of clinical variables

##############################################################################################################
## Using heldout test set 
##############################################################################################################

#################  fit model  ###################
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
                             verboseIter = TRUE)

tic()
XGbTune <- train(x = segTrain[,-1],
                 y = segTrain[,1],
                 method = "xgbTree",
                 tuneLength = hyperparams,
                 metric = "ROC",
                 trControl = adaptControl,
                 preProcess = c("medianImpute", "zv"),
                 verbose = TRUE)

toc_hour()

# saveRDS(XGbTune, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/XGbTune_transcriptomics_clinical.RDS")
XGbTune <- readRDS("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/XGbTune_transcriptomics_clinical.RDS")
XGbTune$finalModel$tuneValue

## predictions on test set:
XGbPred <- predict(XGbTune, segTest[,-1])
cm <- confusionMatrix(XGbPred, segTest[,1])
draw_confusion_matrix(cm, ratio = TRUE)

#Variable importance
var_imp <- varImp(XGbTune, scale = FALSE)
var_imp$importance

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
  labs(title = paste0("AUROC =", signif(auc(XGbROC), 2)))


