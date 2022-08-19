library(tidyverse)
library(plyr)
library(foreign)
library(caret)
library(limma)
library(xgboost)
library(tictoc)
library(SHAPforxgboost)
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
  # title('CONFUSION MATRIX', cex.main=2)
  
  # create the matrix 
  classes = colnames(cm$table)
  rect(150, 430, 240, 370, col=getColor("green", res[1]))
  text(195, 435, classes[1], cex=3)
  rect(250, 430, 340, 370, col=getColor("red", res[3]))
  text(295, 435, classes[2], cex=3)
  text(125, 370, 'Predicted', cex=3, srt=90, font=2)
  text(245, 450, 'Actual', cex=3, font=2)
  rect(150, 305, 240, 365, col=getColor("red", res[2]))
  rect(250, 305, 340, 365, col=getColor("green", res[4]))
  text(140, 400, classes[1], cex=3, srt=90)
  text(140, 335, classes[2], cex=3, srt=90)
  
  # add in the cm results
  text(195, 400, res[1], cex=2.5, font=2, col='black')
  text(195, 335, res[2], cex=2.5, font=2, col='black')
  text(295, 400, res[3], cex=2.5, font=2, col='black')
  text(295, 335, res[4], cex=2.5, font=2, col='black')
  
  # add in the specifics 
  plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
  text(10, 85, names(cm$byClass[1]), cex=3, font=2)
  text(10, 70, round(as.numeric(cm$byClass[1]), 3), cex=3)
  text(40, 85, names(cm$byClass[2]), cex=3, font=2)
  text(40, 70, round(as.numeric(cm$byClass[2]), 3), cex=3)
  text(65, 85, names(cm$byClass[5]), cex=3, font=2)
  text(65, 70, round(as.numeric(cm$byClass[5]), 3), cex=3)
  text(90, 85, names(cm$byClass[6]), cex=3, font=2)
  text(90, 70, round(as.numeric(cm$byClass[6]), 3), cex=3)
  
  # add in the accuracy information 
  text(10, 35, names(cm$overall[1]), cex=3, font=2)
  text(10, 20, round(as.numeric(cm$overall[1]), 3), cex=3)
  text(40, 35, names(cm$byClass[11]), cex=3, font=2)
  text(40, 20, round(as.numeric(cm$byClass[11]), 3), cex=3)
  text(65, 35, names(cm$overall[2]), cex=3, font=2)
  text(65, 20, round(as.numeric(cm$overall[2]), 3), cex=3)
  text(90, 35, names(cm$byClass[7]), cex=3, font=2)
  text(90, 20, round(as.numeric(cm$byClass[7]), 3), cex=3)
}
###############

##create dataframe: 1170 individual (join individuals with defined outcome and 1837 individuals with proteomics measurements)
w3Chronicity <- read.spss("/Users/philippehabets/Dropbox/STRESS_INDEX/data/outcome_groups/DSM_groups.sav", to.data.frame = TRUE, use.missings = FALSE)
w3Chronicity <- w3Chronicity[,c(1,3)]
prtm <- read_csv("/Users/philippehabets/Dropbox/STRESS_INDEX/data/blood_and_saliva_variables/W1/proteomics/output/proteomics_replaced_outliers.csv") %>% 
  mutate(Pident = as.factor(Pident)) 

# plot(density(unlist(flatten(prtm[,-c(1,2)])), na.rm=T))
whole_set <- inner_join(w3Chronicity, prtm, by = c("pident" = "Pident")) %>% 
  mutate(applate = as.factor(applate))
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

##############################################################################################################
##############################################################################################################
# store segTrain without batchcorrection for tSNE plotting (deal with NA's for tSNE here by bagImputation and zv)
# train_uncorrected <- segTrain
# pp_uncorr <- preProcess(train_uncorrected[, -c(1:3)], method = c("zv", "bagImpute"))
# transformedTrain_uncorr <- predict(pp_uncorr, newdata = train_uncorrected[, -c(1:3)])
# train_uncorrected[, -c(1:3)] <- transformedTrain_uncorr
# write_csv(train_uncorrected, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/output/prtm_train_uncorrected.csv")
##############################################################################################################
##############################################################################################################

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
  batch_median[[batch]] <- sapply(batch_list[[plates[batch]]][,-c(1:3)], median.default, na.rm=T)
  batch_IQR[[batch]] <- sapply(batch_list[[plates[batch]]][,-c(1:3)], IQR, na.rm=T)
}
names(batch_median) <- plates
names(batch_IQR) <- plates

batch_srs <- vector(mode = "list", length = length(plates))
for(batch in 1:length(plates)) {
  batch_srs[[batch]] <- setNames(data.frame(matrix(ncol = ncol(whole_set), nrow = nrow(batch_list[[batch]]))), colnames(whole_set))
  for(col in 1:(ncol(prtm)-2)) { #-2 because first two columns are no analytes
  column <- batch_list[[batch]][col+3] # +3 because first three columns are no analytes
  if(!all(is.na(column)) && batch_IQR[[batch]][col] !=0) { # keep constants
    column <- sapply(column, function(analyte) {
      analyte <- 1/(1+exp(-(analyte-(batch_median[[batch]][col]))/(batch_IQR[[batch]][col]/1.35)))
    })
  }
  batch_srs[[batch]][col+3] <- column
  batch_srs[[batch]][1:3] <- batch_list[[batch]][1:3]
  }
}
names(batch_srs) <- plates

# rescale to unit interval
batch_srs_min <- list()
batch_srs_max <- list()
for(batch in 1:length(plates)){
  batch_srs_min[[batch]] <- sapply(batch_srs[[plates[batch]]][,-c(1:3)], min, na.rm=T)
  batch_srs_max[[batch]] <- sapply(batch_srs[[plates[batch]]][,-c(1:3)], max, na.rm=T)
}
names(batch_srs_min) <- plates
names(batch_srs_max) <- plates

batch_srs_norm <- vector(mode = "list", length = length(plates))
for(batch in 1:length(plates)) {
  batch_srs_norm[[batch]] <- setNames(data.frame(matrix(ncol = ncol(whole_set), nrow = nrow(batch_srs[[batch]]))), colnames(whole_set))
  for(col in 1:(ncol(prtm)-2)) { #-2 because first two columns are no analytes
    column <- batch_srs[[batch]][col+3] # +3 because first three columns are no analytes
    if(!all(is.na(column)) && !is.na(batch_srs_max[[batch]][col]) && batch_srs_max[[batch]][col] != batch_srs_min[[batch]][col]) { # keep constants
      column <- sapply(column, function(analyte) {
        analyte <- (analyte-batch_srs_min[[batch]][col])/(batch_srs_max[[batch]][col]-batch_srs_min[[batch]][col])
      })
    }
    batch_srs_norm[[batch]][col+3] <- column
    batch_srs_norm[[batch]][1:3] <- batch_srs[[batch]][1:3]
  }
}
names(batch_srs_norm) <- plates

# transform back to single dataframe
segTrain <- bind_rows(batch_srs_norm)
segTrain <- segTrain %>% mutate_at(vars(-c(1:3)), as.numeric) #because columns are now matrices so turn to numerical vectors

# Because in some batches there was no variance and/or IQR, and min and max were equal, values were kept original here. Because either strategy to deal with this might result in preservation of batch information, 
# we drop these columns here.
no_IQR <- names(which(sapply(segTrain[,-c(1:3)], function(x){ return(max(x, na.rm=T)>1)})))
segTrain <- segTrain %>% select(-no_IQR)

# remove columns with only NA values 
# segTrain <- segTrain %>% select_if(~sum(!is.na(.)) > 0)
trainCols <- colnames(segTrain) # store for later test set transformation

# turn all NaN/NAs to NA
to_na <- function(x){
  x[is.nan(x)] <- NA
  return(as.numeric(x))
}
segTrain <- segTrain %>% mutate_at(vars(-c(1:3)), to_na)

## remove zero variance variables and impute missing variables
pp <- preProcess(segTrain[, -c(1:3)], method = c("zv", "bagImpute"), na.remove = F)
transformedTrain <- predict(pp, newdata = segTrain[, -c(1:3)])
segTrain <- bind_cols(segTrain[,c(1:3)], transformedTrain)

## shuffle rows back so no ordering on batch
shuffled <- as.data.frame(whole_set[trainIndex,])
segTrain <- inner_join(shuffled[,1:2], segTrain, by = c("pident","Remitted_depression"))

# write_rds(segTrain, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/output/segTrain_bagImputedPreProcessingProteomicsBatchCorrected.RDS")
segTrain <- read_rds("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/output/segTrain_bagImputedPreProcessingProteomicsBatchCorrected.RDS")
# write_csv(segTrain, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/output/segTrain_bagImputedPreProcessingProteomicsBatchCorrected.csv")
# segTrain <- read_csv("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/output/segTrain_bagImputedPreProcessingProteomicsBatchCorrected.csv") %>%
#   mutate(Remitted_depression = factor(Remitted_depression, levels=c("remitted", "not_remitted"))) %>% as.data.frame()

## store batch corrected train set for tSNE
# train_batchCorrected <- segTrain
# write_csv(train_batchCorrected, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/output/prtm_train_batchCorrected.csv")

#remove pident and batch as predictor variables
segTrain <- segTrain %>% select(-c(pident, applate))

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
    column <- batch_list_y[[batch]][col+3] # +3 because first three columns are no analytes
    if(!all(is.na(column)) && batch_IQR[[batch]][col] !=0) { # keep constants
      column <- sapply(column, function(analyte) {
        analyte <- 1/(1+exp(-(analyte-(batch_median[[batch]][col]))/(batch_IQR[[batch]][col]/1.35)))
      })
    }
    batch_srs_y[[batch]][col+3] <- column
    batch_srs_y[[batch]][1:3] <- batch_list_y[[batch]][1:3]
  }
}
names(batch_srs_y) <- plates

# Rescale to unit interval using min and max of analyte from trainings data.
batch_srs_norm_y <- vector(mode = "list", length = length(plates))
for(batch in 1:length(plates)) {
  batch_srs_norm_y[[batch]] <- setNames(data.frame(matrix(ncol = ncol(whole_set), nrow = nrow(batch_srs_y[[batch]]))), colnames(whole_set))
  for(col in 1:(ncol(prtm)-2)) { #-2 because first two columns are no analytes
    column <- batch_srs_y[[batch]][col+3] # +3 because first three columns are no analytes
    if(!all(is.na(column)) && !is.na(batch_srs_max[[batch]][col]) && batch_srs_max[[batch]][col] != batch_srs_min[[batch]][col]) { # keep constants
      column <- sapply(column, function(analyte) {
        analyte <- (analyte-batch_srs_min[[batch]][col])/(batch_srs_max[[batch]][col]-batch_srs_min[[batch]][col])
      })
    }
    batch_srs_norm_y[[batch]][col+3] <- column
    batch_srs_norm_y[[batch]][1:3] <- batch_srs_y[[batch]][1:3]
  }
}
names(batch_srs_norm_y) <- plates


# transform back to single dataframe
segTest <- bind_rows(batch_srs_norm_y)
segTest[,-c(1:3)] <- sapply(segTest[,-c(1:3)], as.numeric)

# keep only columns that were kept in trainset
segTest <- segTest %>% select(all_of(trainCols))

## remove zero variance variables and impute missing variables based on trainingsdata
transformedTest <- predict(pp, newdata = segTest[, -c(1:3)])
segTest <- bind_cols(segTest[, c(1:3)], transformedTest)

## function to truncate x < 0 and x > 1 to 0 and 1 (this is due to scaling with minMax based on traindata)
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
segTest <- segTest %>% mutate_at(vars(-c(1:3)), truncate) 

# write_rds(segTest, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/output/segTest_bagImputedPreProcessingProteomicsBatchCorrected.RDS")
segTest <- read_rds("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/output/segTest_bagImputedPreProcessingProteomicsBatchCorrected.RDS")
# write_csv(segTest, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/output/segTest_bagImputedPreProcessingProteomicsBatchCorrected.csv")
# segTest <- read_csv("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/output/segTest_bagImputedPreProcessingProteomicsBatchCorrected.csv") %>%
#   mutate(Remitted_depression = factor(Remitted_depression, levels=c("remitted", "not_remitted"))) %>% as.data.frame()


## remove pident and applate column
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
XGbTune <- train(x = segTrain[,-1],
                 y = segTrain[,1],
                 method = "xgbTree",
                 tuneLength = hyperparams,
                 #tuneGrid = grid,
                 metric = "ROC",
                 trControl = adaptControl, #train_ctrl, 
                 preProcess = c(),
                 verbose = TRUE)
toc_hour()

# saveRDS(XGbTune, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/output/XGbTune_batchCorrected.RDS")
XGbTune <- readRDS("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/output/XGbTune_batchCorrected.RDS")

XGbTune$finalModel$tuneValue
plot(hist(XGbTune$resample$ROC))

## predictions on test set:
XGbPred <- predict(XGbTune, segTest[,-1])
cm <- confusionMatrix(XGbPred, segTest[,1])
draw_confusion_matrix(cm, ratio = TRUE)

accuracies <- list()
for(i in seq(0.4, 0.6, 0.01)){
  pred_ <- factor(ifelse(predict(XGbTune, segTest[,-1], type="prob")$remitted>i,"remitted","not_remitted"), levels = c("remitted", "not_remitted"))
  cm_ <- confusionMatrix(pred_, segTest[,1])
  accuracies <- append(accuracies, cm_$overall[1])
}
names(accuracies) <- seq(0.4, 0.6, 0.01)
opt_thresh <- as.numeric(names(accuracies[accuracies == max(sapply(accuracies, max))]))[1]
XGbPred <- factor(ifelse(predict(XGbTune, segTest[,-1], type="prob")$remitted>opt_thresh,"remitted","not_remitted"), levels = c("remitted", "not_remitted"))
cm <- confusionMatrix(XGbPred, as.factor(segTest$Remitted_depression))
draw_confusion_matrix(cm, ratio = TRUE) # export save as pdf with dim 10x10inch

# variable importance
var_imp <- varImp(XGbTune, scale = FALSE)
var_imp$importance

analyte_names <- read_csv("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/analyte_names.csv")

analyte_importance <- var_imp$importance %>% 
  mutate(symbol = str_sub(rownames(var_imp$importance), start = 3)) %>% 
  left_join(y = analyte_names, by = c("symbol" = "Analyte abbreviation")) %>% 
  dplyr::rename(Gain = Overall) %>% 
  mutate(Gain = round(Gain, 5)) %>% 
  select(2,3,1)
# write_excel_csv(analyte_importance, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/prot_only_analyte_importance_XGB-Gain.csv")

#Get predicted class probabilities so we can build ROC curve.
XGbProbs <- predict(XGbTune, segTest[,-1], type="prob")
head(XGbProbs)

#build ROC curve
# Setting levels: control = not_remitted, case = remitted --> same as in confusion matrix
# Setting direction: controls < cases
XGbROC <- roc(factor(segTest[,1], levels = c("not_remitted","remitted")), XGbProbs[,"remitted"])

# saveRDS(XGbROC, file = "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/roc_2ychronicity_proteomics.RDS")
auc(XGbROC)
pROC::ci.auc(XGbROC)

#plot ROC curve
closest <- function(ROC, optimal_threshold) {
  # function for determining coordinates of specificity and sensitivitywith optimal cutoff for accuracy
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
  geom_point(aes(x = cutoff$specificity, y = cutoff$sensitivity), size = 5, color = "#F87474") +
  geom_text(data = subset(plot_df, specificity==cutoff$specificity & sensitivity ==cutoff$sensitivity), 
            aes(label = paste0("balanced \naccuracy = ", as.character(round(cm$byClass[11],2))), 
                hjust = cutoff$specificity-0.85, 
                vjust = cutoff$sensitivity+0.5), 
            colour = "black", 
            size=6) +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
  theme_classic() +
  theme(text = element_text(size=18)) +
  labs(title = paste0("AUROC =", signif(auc(XGbROC), 2))) +
  theme(plot.title = element_text(size=16))
plot

# ggsave("/Users/philippehabets/Dropbox/STRESS_INDEX/manuscripts/prediction_2ychronicity/Figures/output/auc_proteomics.pdf",
#        plot,
#        width = 15, height = 15, units = "cm")



#####################################################################################################################  
# plots to check effect batch correction
#####################################################################################################################
library(Rtsne)
library(wesanderson)
library(umap)
library(gmodels)

train_uncorrected <- read_csv("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/output/prtm_train_uncorrected.csv")
train_batchCorrected <- read_csv("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/output/prtm_train_batchCorrected.csv")

tsne_uncorrected <- Rtsne(as.matrix(train_uncorrected[,-c(1:3)]), max_iter = 5000, check_duplicates = TRUE, perplexity = 15, theta = 0.0, pca = FALSE, normalize = F, verbose = TRUE) ##pca_scale (& -center) = TRUE if pca = TRUE
df_uncorrected <- cbind(pident = train_uncorrected$pident, batch = train_uncorrected$applate, as.data.frame(tsne_uncorrected$Y))

tsne_batchCorrected <- Rtsne(as.matrix(train_batchCorrected[,-c(1:3)]), max_iter = 5000, check_duplicates = TRUE, perplexity = 15, theta = 0.0, pca = FALSE, normalize = F, verbose = TRUE) ##pca_scale (& -center) = TRUE if pca = TRUE
df_batchCorrected <- cbind(pident = train_batchCorrected$pident, batch = train_batchCorrected$applate, as.data.frame(tsne_batchCorrected$Y))

## look at plate (batch) effect without 'correction':
pal <- wes_palette("Zissou1", 26, type = "continuous")
plot_uncorrected <- ggplot(df_uncorrected, aes(x = V1, y = V2, color = batch)) +
  geom_point() +
  scale_color_gradientn(colours = pal) +
  xlab("tSNE 1") +
  ylab("tSNE 2") +
  theme_minimal() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "white"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())

plot_uncorrected
# ggsave("/Users/philippehabets/Dropbox/STRESS_INDEX/manuscripts/prediction_2ychronicity/supplements/plots/proteomics_uncorrected_tsne.pdf",
#               plot_uncorrected,
#               width = 12,
#               height = 10,
#               units = "cm")

## look at plate (batch) effect with batch 'correction':
pal <- wes_palette("Zissou1", 26, type = "continuous")
plot_corrected <- ggplot(df_batchCorrected, aes(x = V1, y = V2, color = batch)) +
  geom_point() +
  scale_color_gradientn(colours = pal) +
  xlab("tSNE 1") +
  ylab("tSNE 2") +
  theme_minimal() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "white"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())

plot_corrected
# ggsave("/Users/philippehabets/Dropbox/STRESS_INDEX/manuscripts/prediction_2ychronicity/supplements/plots/proteomics_corrected_tsne.pdf",
#        plot_corrected,
#        width = 12,
#        height = 10,
#        units = "cm")


###### SHAP ########
analyte_names <- read_csv("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/analyte_names.csv")

# With SHAP for xgboost package:
dataX <- as.matrix(segTrain[,XGbTune$finalModel$feature_names])
shap_values <- shap.values(xgb_model = XGbTune$finalModel, X_train = dataX)
shap_values$mean_shap_score

shap.plot.summary.wrap1(XGbTune$finalModel,dataX, top_n=50, dilute = F)

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
    #if string does not start with "a", it is either 'Age' or "Sex" variable: keep all letters
    x <- x
  }
  return(x)
}

analyte_shap <- data.frame(symbol = sapply(names(shap_values$mean_shap_score), strip_ap), 
                           mean_shap_score =  shap_values$mean_shap_score) %>% 
  left_join(y = analyte_names, by = c("symbol" = "Analyte abbreviation")) %>% 
  select(1,3,2)
# write_excel_csv(analyte_shap, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/shap_xgb_proteomics.csv")

# save analytes with mean absolute Shapley value > 0 for inclusion in all-data-modalities model
shap_proteomics <-names(shap_values$mean_shap_score[shap_values$mean_shap_score>0])
# saveRDS(shap_proteomics, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/proteomics/data/shap_values_proteomics.RDS")























