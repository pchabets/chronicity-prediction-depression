library(tidyverse)
library(plyr)
library(foreign)
library(caret)
library(kernlab)
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
  text(30, 35, names(cm$overall[1]), cex=1.5, font=2)
  text(30, 20, round(as.numeric(cm$overall[1]), 3), cex=1.4)
  text(70, 35, names(cm$overall[2]), cex=1.5, font=2)
  text(70, 20, round(as.numeric(cm$overall[2]), 3), cex=1.4)
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

whole_set <- left_join(w3Chronicity, demographic, by = "pident") %>% 
  left_join(bigFive, by = "pident") %>% 
  left_join(IDS, by = "pident") %>% 
  mutate_if(is.factor, droplevels) %>% 
  mutate(aidssev = as.numeric(aidssev)) #turn ordinal data into numerical data
  
options(na.action = "na.pass")
modelMatrix <- model.matrix(Remitted_depression ~ ., data = whole_set[,-1])[,-1] #transform categorical nominal features into numeric values
whole_set <- cbind(Remitted_depression = whole_set$Remitted_depression, as.data.frame(modelMatrix))

whole_set$Remitted_depression <- revalue(whole_set$Remitted_depression, c("remitted"="remitted", "non-remitted"="not_remitted")) 
table(whole_set$Remitted_depression) # 407 (51%) remitted, 397 (49%) non-remitted

##############################################################################################################
## Using heldout test set 
##############################################################################################################
## partition data
set.seed(42)

trainIndex <- createDataPartition(y=whole_set$Remitted_depression, times=1, p=0.8, list=F)

segTrain <- as.data.frame(whole_set[trainIndex,])
prepTrain <- preProcess(segTrain[, 3:ncol(segTrain)],
                   #used preProcessing arguments only for bigFive and IDS numerical data (ordinal categorical severity data excluded)
                   method = c("center", "scale", "medianImpute"))
transformedTrain <- predict(prepTrain, newdata = segTrain[, 3:ncol(segTrain)])
segTrain[3:ncol(segTrain)] <- transformedTrain

# similarly preprocess testdata
segTest <- as.data.frame(whole_set[-trainIndex,])
transformedTest <- predict(prepTrain, newdata = segTest[, 3:ncol(segTest)])
segTest[3:ncol(segTest)] <- transformedTest

summary(segTrain$Remitted_depression) #remitted/non-remitted = 326/318
summary(segTest$Remitted_depression) #remitted/non-remitted = 81/79

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

# grid <- expand.grid(sigma = seq(0.001, 0.0015, length.out = 25),
#                     C = seq(6, 6.2, length.out = 25))

svmTune <- train(x = segTrain[,-1],
                 y = segTrain[,1],
                 method = "svmLinear",
                 tuneLength = 50,
                 # tuneGrid = grid,
                 metric = "ROC",
                 trControl = train_ctrl,
                 preProcess = NULL,
                 verbose = TRUE)

svmTune
# saveRDS(svmTune, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/svmTune3_linear.RDS")
svmTune <- readRDS("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/svmTune3_linear.RDS")
svmTune$finalModel

plot(svmTune, metric = "ROC", scales = list(x = list(log =2)))
roc_imp <- varImp(svmTune, scale = FALSE) # ROC per variable
roc_imp

## predictions on test set:
svmPred <- predict(svmTune, segTest[,-1])
cm <- confusionMatrix(svmPred, segTest[,1])
draw_confusion_matrix(cm, ratio = TRUE)

#Get predicted class probabilities so we can build ROC curve.
svmProbs <- predict(svmTune, segTest[,-1], type="prob")
head(svmProbs)

#build ROC curve
svmROC <- roc(segTest[,1], svmProbs[,"remitted"])
# saveRDS(svmROC, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/svmL_roc_2ychr_clin.RDS")
auc(svmROC)

#plot ROC curve
plot(svmROC, type = "S")
ggplot(data.frame(specificity = svmROC$specificities, sensitivity = svmROC$sensitivities), aes(specificity, sensitivity)) +
  geom_path() +
  scale_x_reverse() +
  geom_abline(intercept = 1, slope = 1, colour='grey') +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
  theme_classic() +
  labs(title = paste0("AUROC =", signif(auc(svmROC), 3)))









