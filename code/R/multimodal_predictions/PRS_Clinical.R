library(tidyverse)
library(plyr)
library(naniar)
library(caret)
library(xgboost)
library(glmnet)
library(tictoc)
library(foreign)
library(doMC)
library(pROC)

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
###

##create dataframe from all data
w3Chronicity <- read.spss("data/outcome_groups/DSM_groups.sav", to.data.frame = TRUE, use.missings = FALSE)
w3Chronicity <- w3Chronicity[,c(1,3)]

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

# PRS data
# read in PRS data
prs <- read_tsv("genetic_data/PRS_NESDA_NTR/PRSs_NESDA.tsv") %>% 
  mutate(FAMID = str_trim(FAMID)) 
key <- read.spss("genetic_data/PRS_NESDA_NTR/Key GODOT.sav", to.data.frame = T, use.missings = F) %>% 
  mutate(DNAid = str_trim(DNAid)) 

#transform PRS IDs to NESDA IDs
prs <- prs %>% 
  inner_join(key, by = c("FAMID" = "DNAid")) %>% 
  relocate(pident, .before = FAMID) %>% 
  select(-c(FAMID, INDID, Age, sex, DNAstudy)) %>% 
  mutate(pident = as.factor(pident))

#merge into dataframe
whole_set <- inner_join(w3Chronicity, prs, by = "pident") %>% 
  left_join(demographic, by = "pident") %>% 
  left_join(bigFive, by = "pident") %>% 
  left_join(IDS, by = "pident") %>% 
  mutate_if(is.factor, droplevels) %>% 
  mutate(aidssev = as.numeric(aidssev)) #turn ordinal data into numerical data

options(na.action = "na.pass")
modelMatrix <- model.matrix(Remitted_depression ~ ., data = whole_set[,-1])[,-1] #transform categorical nominal features into numeric values
whole_set <- cbind(pident = whole_set$pident, Remitted_depression = whole_set$Remitted_depression, as.data.frame(modelMatrix))

whole_set$Remitted_depression <- revalue(whole_set$Remitted_depression, c("remitted"="remitted", "non-remitted"="not_remitted")) 
table(whole_set$Remitted_depression) 


##############################################################################################################
## Using heldout test set 
##############################################################################################################
## partition data
set.seed(42)
trainIndex <- createDataPartition(y=whole_set$Remitted_depression, times=1, p=0.8, list=F)

# training set
segTrain <- as.data.frame(whole_set[trainIndex,-1])

# testing set (preprocess similarly, i.e. use train imputations on testdata)
segTest <- as.data.frame(whole_set[-trainIndex,-1])

summary(segTrain$Remitted_depression) #remitted/non-remitted = 289/272
summary(segTest$Remitted_depression) #remitted/non-remitted = 72/68

## Fit model
#function to time progress:
toc_hour <- function() {
  x <- toc()
  hours <- (x$toc - x$tic) / 3600
  cat(paste(c("or ", as.character(round(hours, digits = 2)), " hours", "\n"),collapse = ""))
}

#Set number of combinations of hyperparameters to use
hyperparams <- 1000

## set seed list for reproduction
set.seed(42)
seeds <- vector(mode = "list", length = 101)
for(i in 1:100) seeds[[i]] <- sample.int(1000, hyperparams)
seeds[[101]] <- sample.int(100,1)

#Function to run with different model:
train_model <- function(trainset, method = "glmnet", preProcess = c(), preProcOptions = list()) {
  ### trainset is dataframe with first column = outcomes, other columns are features
  ### method is method to pass to train() function as character, e.g. "xgbTree" or "glmnet"
  ### preProcess = character vector to pass to preProcess argument in train() function
  ### preProcOptions = character vector to pass to preProcOptions argument in train() function
  train_ctrl <- trainControl(method = "adaptive_cv",
                             number = 10, repeats = 10,
                             adaptive = list(min = 10, alpha = 0.05, method = "gls", complete = TRUE), #using adaptive resampling to have more efficient computation, see http://topepo.github.io/caret/adaptive-resampling.html
                             summaryFunction = twoClassSummary,
                             search = "random",
                             classProbs = TRUE,
                             preProcOptions = preProcOptions,
                             # seeds = seeds,
                             allowParallel = TRUE,
                             verboseIter = TRUE)
  
  tic()
  model <- train(x = segTrain[,-1],
                   y = segTrain[,1],
                   method = method,
                   tuneLength = hyperparams,
                   #tuneGrid = grid,
                   metric = "ROC",
                   trControl = train_ctrl,
                   preProcess = preProcess,
                   verbose = TRUE)
  toc_hour()
  return(model)
}

#xgboost:
XGbTune <- train_model(segTrain, method = "xgbTree", preProcess = c("medianImpute"))
# saveRDS(XGbTune, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/output/XGbTune_29PRS_clinical.RDS")
XGbTune <- readRDS("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/output/XGbTune_29PRS_clinical.RDS")
XGbTune$finalModel$tuneValue
plot(density(XGbTune$resample$ROC))

## predictions on test set:
XGbPred <- predict(XGbTune, segTest[,-1])
cm <- confusionMatrix(XGbPred, as.factor(segTest[,1]))
draw_confusion_matrix(cm, ratio = TRUE)

# Opimize accuracy given ROC
accuracies <- list()
for(i in seq(0.4, 0.6, 0.01)){
  pred_ <- factor(ifelse(predict(XGbTune, segTest[,-1], type="prob")$remitted>i,"remitted","not_remitted"), levels = c("remitted", "not_remitted"))
  cm_ <- confusionMatrix(pred_, segTest[,1])
  accuracies <- append(accuracies, cm_$overall[1])
}
names(accuracies) <- seq(0.4, 0.6, 0.01)
opt_thresh <- as.numeric(names(accuracies[accuracies == max(sapply(accuracies, max))]))
XGbPred <- factor(ifelse(predict(XGbTune, segTest[,-1], type="prob")$remitted>opt_thresh,"remitted","not_remitted"), levels = c("remitted", "not_remitted"))
cm <- confusionMatrix(XGbPred, as.factor(segTest[,1]))
draw_confusion_matrix(cm, ratio = TRUE)

# Variable importance
var_imp <- varImp(XGbTune, scale = FALSE)
var_imp$importance

#Get predicted class probabilities so we can build ROC curve.
XGbProbs <- predict(XGbTune, segTest[,-1], type="prob")
head(XGbProbs)

#build ROC curve
XGbROC <- roc(segTest[,1], XGbProbs[,"not_remitted"])
auc(XGbROC)

#plot ROC curve
ggplot(data.frame(specificity = XGbROC$specificities, sensitivity = XGbROC$sensitivities), aes(specificity, sensitivity)) +
  geom_path() +
  scale_x_reverse() +
  geom_abline(intercept = 1, slope = 1, colour='grey') +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
  theme_classic() +
  labs(title = paste0("AUROC =", signif(auc(XGbROC), 3)))


# #elastic net:
# elnet <- train_model(segTrain, method = "glmnet", preProcess = c("medianImpute", "center", "scale")) 
# 
# ## predictions on test set:
# elnetPred <- predict(elnet, segTest[,-1])
# cm <- confusionMatrix(elnetPred, as.factor(segTest[,1]))
# draw_confusion_matrix(cm, ratio = TRUE)
# var_imp <- varImp(elnet, scale = FALSE)
# var_imp
# 
# #Get predicted class probabilities so we can build ROC curve.
# elnetProbs <- predict(elnet, segTest[,-1], type="prob")
# head(elnetProbs)
# 
# #build ROC curve
# elnetROC <- roc(segTest[,1], elnetProbs[,"not_remitted"])
# auc(elnetROC)
# 
# #plot ROC curve
# ggplot(data.frame(specificity = elnetROC$specificities, sensitivity = elnetROC$sensitivities), aes(specificity, sensitivity)) +
#   geom_path() +
#   scale_x_reverse() +
#   geom_abline(intercept = 1, slope = 1, colour='grey') +
#   coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
#   theme_classic() +
#   labs(title = paste0("AUROC =", signif(auc(elnetROC), 3)))








