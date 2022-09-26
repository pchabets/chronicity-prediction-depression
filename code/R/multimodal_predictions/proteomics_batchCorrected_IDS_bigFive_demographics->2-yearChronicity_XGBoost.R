library(tidyverse)
library(plyr)
library(foreign)
library(caret)
library(xgboost)
library(doMC)
library(pROC)
library(tictoc)
library(SHAPforxgboost)

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
  text(10, 70, round(as.numeric(cm$byClass[1]), 2), cex=1.2)
  text(30, 85, names(cm$byClass[2]), cex=1.2, font=2)
  text(30, 70, round(as.numeric(cm$byClass[2]), 2), cex=1.2)
  text(50, 85, names(cm$byClass[5]), cex=1.2, font=2)
  text(50, 70, round(as.numeric(cm$byClass[5]), 2), cex=1.2)
  text(70, 85, names(cm$byClass[6]), cex=1.2, font=2)
  text(70, 70, round(as.numeric(cm$byClass[6]), 2), cex=1.2)
  text(90, 85, names(cm$byClass[7]), cex=1.2, font=2)
  text(90, 70, round(as.numeric(cm$byClass[7]), 2), cex=1.2)
  
  # add in the accuracy information 
  text(10, 35, names(cm$overall[1]), cex=1.5, font=2)
  text(10, 20, round(as.numeric(cm$overall[1]), 2), cex=1.4)
  text(50, 35, names(cm$byClass[11]), cex=1.5, font=2)
  text(50, 20, round(as.numeric(cm$byClass[11]), 2), cex=1.4)
  text(90, 35, names(cm$overall[2]), cex=1.5, font=2)
  text(90, 20, round(as.numeric(cm$overall[2]), 2), cex=1.4)
}
###############

##create dataframe
w3Chronicity <- read.spss("outcome_groups/DSM_groups.sav", to.data.frame = TRUE, use.missings = FALSE)
w3Chronicity <- w3Chronicity[,c(1,3)]

prtm <- read_csv("blood_and_saliva_variables/W1/proteomics/output/proteomics_replaced_outliers.csv")
prtm <- prtm %>% 
  mutate(Pident = as.factor(Pident))
  
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
  inner_join(prtm, by = c("pident" = "Pident")) %>% #keeping only cases with proteomics measured
  mutate_if(is.factor, droplevels) %>% 
  mutate(applate = as.factor(applate)) %>% 
  mutate(aidssev = as.numeric(aidssev)) #turn ordinal data into numerical data

options(na.action = "na.pass")
modelMatrix <- model.matrix(Remitted_depression ~ ., data = whole_set %>% select(-c(pident, applate)))[,-1] #transform categorical nominal features into numeric values
whole_set <- bind_cols(pident = whole_set$pident, Remitted_depression = whole_set$Remitted_depression, applate = whole_set$applate, as.data.frame(modelMatrix))

whole_set$Remitted_depression <- revalue(whole_set$Remitted_depression, c("remitted"="remitted", "non-remitted"="not_remitted")) 
table(whole_set$Remitted_depression) # 312 (51%) remitted, 299 (49%) non-remitted

##############################################################################################################
## Calculate effect of batch on outcome label distribution
##############################################################################################################
contingency <- table(whole_set$Remitted_depression, whole_set$applate)
chisq.test(contingency, simulate.p.value = T)

##############################################################################################################
## Preprocess by removing batch effect using srs-normalisation per batch in train-set, and transform test set 
## based on those srs-normalisations
##############################################################################################################
##############################################################################################################
### partition data ###########################################################################################
set.seed(42)
trainIndex <- createDataPartition(y=whole_set$Remitted_depression, times=1, p=0.8, list=F)

segTrain <- as.data.frame(whole_set[trainIndex,])
segTest <- as.data.frame(whole_set[-trainIndex,])

summary(segTrain$Remitted_depression) #remitted/non-remitted = 250/240
summary(segTest$Remitted_depression) #remitted/non-remitted = 62/59

## Proceed with segTrain (no imputations have occurred yet)
##############################################################################################################
## for each batch (applate), calculate median and IQR for each analyte in train set #############################
plates <- levels(segTrain$applate)

# create list of dataframes for each batch
batch_list <- list()
for(batch in 1:length(plates)) {
  batch_list[[batch]] <- segTrain %>% filter(applate == plates[batch])
}
names(batch_list) <- plates

# for each batch, calculate mean and sd for each analyte
batch_median <- list()
batch_IQR <- list()
for(batch in 1:length(plates)){
  batch_median[[batch]] <- sapply(batch_list[[plates[batch]]][,-c(1:13)], median.default, na.rm=T)
  batch_IQR[[batch]] <- sapply(batch_list[[plates[batch]]][,-c(1:13)], IQR, na.rm=T)
}
names(batch_median) <- plates
names(batch_IQR) <- plates

batch_srs <- vector(mode = "list", length = length(plates))
for(batch in 1:length(plates)) {
  batch_srs[[batch]] <- setNames(data.frame(matrix(ncol = ncol(whole_set), nrow = nrow(batch_list[[batch]]))), colnames(whole_set))
  for(col in 1:(ncol(prtm)-2)) { #-2 because first two columns are no analytes
    column <- batch_list[[batch]][col+13] # +13 because first three columns are no analytes
    if(!all(is.na(column)) && batch_IQR[[batch]][col] !=0) { # keep constants
      column <- sapply(column, function(analyte) {
        analyte <- 1/(1+exp(-(analyte-(batch_median[[batch]][col]))/(batch_IQR[[batch]][col]/1.35)))
      })
    }
    batch_srs[[batch]][col+13] <- column
    batch_srs[[batch]][1:13] <- batch_list[[batch]][1:13]
  }
}
names(batch_srs) <- plates

# rescale to unit interval
batch_srs_min <- list()
batch_srs_max <- list()
for(batch in 1:length(plates)){
  batch_srs_min[[batch]] <- sapply(batch_srs[[plates[batch]]][,-c(1:13)], min, na.rm=T)
  batch_srs_max[[batch]] <- sapply(batch_srs[[plates[batch]]][,-c(1:13)], max, na.rm=T)
}
names(batch_srs_min) <- plates
names(batch_srs_max) <- plates

batch_srs_norm <- vector(mode = "list", length = length(plates))
for(batch in 1:length(plates)) {
  batch_srs_norm[[batch]] <- setNames(data.frame(matrix(ncol = ncol(whole_set), nrow = nrow(batch_srs[[batch]]))), colnames(whole_set))
  for(col in 1:(ncol(prtm)-2)) { #-2 because first two columns are no analytes
    column <- batch_srs[[batch]][col+13] # +13 because first three columns are no analytes
    if(!all(is.na(column)) && !is.na(batch_srs_max[[batch]][col]) && batch_srs_max[[batch]][col] != batch_srs_min[[batch]][col]) { # keep constants
      column <- sapply(column, function(analyte) {
        analyte <- (analyte-batch_srs_min[[batch]][col])/(batch_srs_max[[batch]][col]-batch_srs_min[[batch]][col])
      })
    }
    batch_srs_norm[[batch]][col+13] <- column
    batch_srs_norm[[batch]][1:13] <- batch_srs[[batch]][1:13]
  }
}
names(batch_srs_norm) <- plates

# transform back to single dataframe
segTrain <- bind_rows(batch_srs_norm) # or batch_srs of no rescaling is required
segTrain <- segTrain %>% mutate_at(vars(-c(1:13)), as.numeric)

# Because in some batches there was no variance and/or IQR, and min and max were equal, values were kept original here. Because either strategy to deal with this might result in preservation of batch information, 
# we drop these columns here.
no_IQR <- names(which(sapply(segTrain[,-c(1:13)], function(x){ return(max(x, na.rm=T)>1)})))
segTrain <- segTrain %>% select(-no_IQR)

# store columns so we can select same columns in test set
trainCols <- colnames(segTrain) 

## remove zero variance variables and impute missing variables
pp <- preProcess(segTrain[, -c(1:13)], method = c("zv", "bagImpute"))
transformedTrain <- predict(pp, newdata = segTrain[, -c(1:13)])
segTrain <- bind_cols(segTrain[, c(1:13)], transformedTrain)

## shuffle rows back so no ordering on batch
shuffled <- as.data.frame(whole_set[trainIndex,])
segTrain <- inner_join(shuffled[,1:2], segTrain, by = c("pident","Remitted_depression"))

#remove pident and batch as predictor variables
segTrain <- segTrain %>% select(-c(pident, applate))

# write_csv(segTrain, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/output/segTrain_bagImputedPreProcessingProteomicsBatchCorrected_AND_clinical.csv")
segTrain <- read_csv("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/output/segTrain_bagImputedPreProcessingProteomicsBatchCorrected_AND_clinical.csv") %>%
  mutate(Remitted_depression = factor(Remitted_depression, levels=c("remitted", "not_remitted"))) %>% as.data.frame()

##############################################################################################################
## apply the same procedure for the test set using trainingset information only  #############################

# create list of dataframes for each batch
batch_list_y <- list()
for(batch in 1:length(plates)) {
  batch_list_y[[batch]] <- segTest %>% filter(applate == plates[batch])
}
names(batch_list_y) <- plates

# SRS normalize using median and IQR from analytes in trainingsdata
batch_srs_y <- vector(mode = "list", length = length(plates))
for(batch in 1:length(plates)) {
  batch_srs_y[[batch]] <- setNames(data.frame(matrix(ncol = ncol(whole_set), nrow = nrow(batch_list_y[[batch]]))), colnames(whole_set))
  for(col in 1:(ncol(prtm)-2)) { #-2 because first two columns are no analytes
    column <- batch_list_y[[batch]][col+13] # +13 because first three columns are no analytes
    if(!all(is.na(column)) && batch_IQR[[batch]][col] !=0) { # keep constants
      column <- sapply(column, function(analyte) {
        analyte <- 1/(1+exp(-(analyte-(batch_median[[batch]][col]))/(batch_IQR[[batch]][col]/1.35)))
      })
    }
    batch_srs_y[[batch]][col+13] <- column
    batch_srs_y[[batch]][1:13] <- batch_list_y[[batch]][1:13]
  }
}
names(batch_srs_y) <- plates

# Rescale to unit interval using min and max of analyte from trainings data.
batch_srs_norm_y <- vector(mode = "list", length = length(plates))
for(batch in 1:length(plates)) {
  batch_srs_norm_y[[batch]] <- setNames(data.frame(matrix(ncol = ncol(whole_set), nrow = nrow(batch_srs_y[[batch]]))), colnames(whole_set))
  for(col in 1:(ncol(prtm)-2)) { #-2 because first two columns are no analytes
    column <- batch_srs_y[[batch]][col+13] # +13 because first three columns are no analytes
    if(!all(is.na(column)) && !is.na(batch_srs_max[[batch]][col]) && batch_srs_max[[batch]][col] != batch_srs_min[[batch]][col]) { # keep constants
      column <- sapply(column, function(analyte) {
        analyte <- (analyte-batch_srs_min[[batch]][col])/(batch_srs_max[[batch]][col]-batch_srs_min[[batch]][col])
      })
    }
    batch_srs_norm_y[[batch]][col+13] <- column
    batch_srs_norm_y[[batch]][1:13] <- batch_srs_y[[batch]][1:13]
  }
}
names(batch_srs_norm_y) <- plates

# transform back to single dataframe
segTest <- bind_rows(batch_srs_y)
segTest[,-c(1:13)] <- sapply(segTest[,-c(1:13)], as.numeric)

# remove columns with only NA values in trainset
segTest <- segTest %>% select(all_of(trainCols))

## remove zero variance variables and impute missing variables based on trainingsdata
transformedTest <- predict(pp, newdata = segTest[, -c(1:13)])
segTest <- bind_cols(segTest[, c(1:13)], transformedTest)

# function to truncate x < 0 and x > 1 to 0 and 1 respectively (this is due to scaling with minMax based on traindata)
truncate <- function(x) {
  for(i in 1:length(x)) {
    if(x[i] < 0) {
      x[i] <- 0
    }
    else if(x[i] > 1) {
      x[i] <- 1
    }
  }
  return(x)
}
segTest <- segTest %>% mutate_at(vars(-c(1:13)), truncate) 

## remove pident and applate column
segTest <- segTest %>% select(-c(pident, applate))

# write_csv(segTest, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/output/segTest_bagImputedPreProcessingProteomicsBatchCorrected_AND_clinical.csv")
segTest <- read_csv("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/output/segTest_bagImputedPreProcessingProteomicsBatchCorrected_AND_clinical.csv") %>%
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
XGbTune <- train(x = segTrain[,-1],
                 y = segTrain[,1],
                 method = "xgbTree",
                 tuneLength = hyperparams,
                 # tuneGrid = tune_grid,
                 metric = "ROC",
                 trControl = adaptControl,
                 preProcess = c("medianImpute"), # will only be performed on clinical variables as proteomic variables are already imputed
                 verbose = TRUE)

toc_hour()

XGbTune$finalModel$tuneValue

# saveRDS(XGbTune, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/XGbTune2.RDS")
XGbTune <- readRDS("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/XGbTune2.RDS")
XGbTune$finalModel$tuneValue

## predictions on test set:
XGbPred <- predict(XGbTune, segTest[,-1])

## With alternative threshold: 
accuracies <- list()
for(i in seq(0.4, 0.6, 0.01)){
  pred_ <- factor(ifelse(predict(XGbTune, segTest[,-1], type="prob")$remitted>i,"remitted","not_remitted"), levels = c("remitted", "not_remitted"))
  cm_ <- confusionMatrix(pred_, segTest[,1])
  accuracies <- append(accuracies, cm_$overall[1])
}
names(accuracies) <- seq(0.4, 0.6, 0.01)
opt_thresh <- as.numeric(names(accuracies[accuracies == max(sapply(accuracies, max))]))[1]
XGbPred <- factor(ifelse(predict(XGbTune, segTest[,-1], type="prob")$remitted>opt_thresh,"remitted","not_remitted"), levels = c("remitted", "not_remitted"))

# Confustion matrix
cm <- confusionMatrix(XGbPred, segTest[,1], positive = 'remitted')
draw_confusion_matrix(cm, ratio = TRUE)

#Get predicted class probabilities so we can build ROC curve.
XGbProbs <- predict(XGbTune, segTest[,-1], type="prob")
head(XGbProbs)

#build ROC curve
# Setting levels: control = not_remitted, case = remitted --> same as in confusion matrix
# Setting direction: controls < cases

XGbROC <- roc(factor(segTest[,1], levels = c("not_remitted","remitted")), XGbProbs[,"remitted"])
# saveRDS(XGbROC, file = "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/roc_2ychronicity_prot_clin.RDS")
auc(XGbROC)
ci.auc(XGbROC)

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
  geom_text(data = subset(plot_df, specificity==cutoff$specificity & sensitivity ==cutoff$sensitivity), aes(label = paste0("balanced \naccuracy = ", as.character(round(cm$byClass[11],2))), hjust = cutoff$specificity-1, vjust = cutoff$sensitivity+0.7), colour = "black", size=5) +
  coord_fixed(ratio = 0.9, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
  theme_classic() +
  labs(title = paste0("AUROC =", signif(auc(XGbROC), 2))) +
  theme(text = element_text(size=16))
plot

# ggsave("/Users/philippehabets/Dropbox/STRESS_INDEX/manuscripts/prediction_2ychronicity/Figures/output/auc_prot_clinical.pdf",
#        plot,
#        width = 23, height = 15, units = "cm")

