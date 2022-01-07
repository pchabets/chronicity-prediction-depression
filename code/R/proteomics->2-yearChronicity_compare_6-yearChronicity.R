library(tidyverse)
library(plyr)
library(foreign)
library(caret)
library(xgboost)
library(pROC)

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
  text(10, 70, signif(as.numeric(cm$byClass[1]), 2), cex=1.2)
  text(30, 85, names(cm$byClass[2]), cex=1.2, font=2)
  text(30, 70, signif(as.numeric(cm$byClass[2]), 2), cex=1.2)
  text(50, 85, names(cm$byClass[5]), cex=1.2, font=2)
  text(50, 70, signif(as.numeric(cm$byClass[5]), 2), cex=1.2)
  text(70, 85, names(cm$byClass[6]), cex=1.2, font=2)
  text(70, 70, signif(as.numeric(cm$byClass[6]), 2), cex=1.2)
  text(90, 85, names(cm$byClass[7]), cex=1.2, font=2)
  text(90, 70, signif(as.numeric(cm$byClass[7]), 2), cex=1.2)
  
  # add in the accuracy information 
  text(10, 35, names(cm$overall[1]), cex=1.5, font=2)
  text(10, 20, round(as.numeric(cm$overall[1]), 3), cex=1.4)
  text(50, 35, names(cm$byClass[11]), cex=1.5, font=2)
  text(50, 20, round(as.numeric(cm$byClass[11]), 3), cex=1.4)
  text(90, 35, names(cm$overall[2]), cex=1.5, font=2)
  text(90, 20, round(as.numeric(cm$overall[2]), 3), cex=1.4)
}
###############

##create dataframe: 1170 individual (join of 1860 individuals with defined outcome and 1837 individuals with proteomics measurements)
w3Chronicity <- read.spss("/Users/philippehabets/Dropbox/STRESS_INDEX/data/outcome_groups/DSM_groups.sav", to.data.frame = TRUE, use.missings = FALSE)
w3Chronicity <- w3Chronicity[,c(1,3)]
prtm <- read_csv("/Users/philippehabets/Dropbox/STRESS_INDEX/data/blood_and_saliva_variables/W1/proteomics/output/proteomics_replaced_outliers.csv")
prtm$Pident <- as.factor(prtm$Pident) 

whole_set <- inner_join(w3Chronicity, prtm, by = c("pident" = "Pident"))
whole_set$Remitted_depression <- revalue(whole_set$Remitted_depression, c("remitted"="remitted", "non-remitted"="not_remitted"))
table(whole_set$Remitted_depression) # 312 (51%) remitted, 299 (49%) non-remitted

############################################################
## partition data in similar way as trained model.
############################################################
set.seed(42)

trainIndex <- createDataPartition(y=whole_set$Remitted_depression, times=1, p=0.8, list=F)

segTrain <- as.data.frame(whole_set[trainIndex,-c(1,3)])
segTrain[,1] <- as.factor(segTrain[,1])
segTrain <- segTrain %>% 
  #log10 transform as preprocessing step (not included as option in caret, so done manually here)
  mutate_at(vars(-1), log10)

segTest <- as.data.frame(whole_set[-trainIndex,-c(1,3)])
segTest[,1] <- as.factor(segTest[,1])
segTest <- segTest %>% 
  #log10 transform as preprocessing step (not included as option in caret, so done manually here)
  mutate_at(vars(-1), log10)

summary(segTrain$Remitted_depression) #remitted/non-remitted = 250/240
summary(segTest$Remitted_depression) #remitted/non-remitted = 62/59


##############################################################################################################
## Load in xgb model
##############################################################################################################

XGbTune <- readRDS("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/output/XGbTune.RDS")
XGbPred <- predict(XGbTune, segTest[,-1])
cm <- confusionMatrix(XGbPred, segTest[,1])

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
  labs(title = paste0("AUROC =", signif(auc(XGbROC), 3)))


############################################################################################################
# check for false negatives at 2 year chronicity, if 6-year chronicity has same FNR
############################################################################################################
cm$positive #"positive" class = "remitted" 

# we want to know hom many non-remitted were missclassified as "remitted", so false positives.
draw_confusion_matrix(cm, ratio = TRUE)
cm$byClass[2] #specificity = 0.69. So 69% of not_remitted is correctly classified. 
cm$byClass[1] #sensitivity = 0.71. So 71% of remitted is correctly classified. 

#Find out if specificity is higher for the not_remitted cases that were not_remitted for 6 consecutive years.
testSet <- as.data.frame(whole_set[-trainIndex,1:2]) #retrieve PID for matching with 6-y chronicity.
testSet$predicted <- XGbPred

##create dataframe: 
w5Chronicity <- read_csv("/Users/philippehabets/Dropbox/STRESS_INDEX/data/outcome_groups/6yChronicity.csv")
w5Chronicity <- w5Chronicity[,-2]
w5Chronicity <- w5Chronicity %>% 
  mutate_all(as.factor)
w5Chronicity <- w5Chronicity %>% 
  mutate_at(vars(year2Chronicity, year4Chronicity, year6Chronicity), ~revalue(.x, c("not_depressed"="remitted", "depressed"="not_remitted")))

chronicity <- inner_join(testSet, w5Chronicity, by = "pident") # only 75 cases left
sum(chronicity$year2Chronicity == chronicity$year6Chronicity)/nrow(chronicity) *100 # 64% of 2-year and 6-year status overlap.

#filter cases where 2-y, 4-y and 6-y chronicity is true (not remitted at all times)
chronicity_filtered <- chronicity %>% 
  filter(year2Chronicity == "not_remitted" & year4Chronicity == "not_remitted" & year6Chronicity == "not_remitted") # only 13 cases left

#filter cases where 2-y, 4-y and 6-y chronicity is not true (remitted at all times)
chronicity_filtered2 <- chronicity %>% 
  filter(year2Chronicity == "remitted" & year4Chronicity == "remitted" & year6Chronicity == "remitted") # only 26 controls left


sum(chronicity_filtered$predicted == chronicity_filtered$Remitted_depression)/nrow(chronicity_filtered) # specificity = 77%
sum(chronicity_filtered2$predicted == chronicity_filtered2$Remitted_depression)/nrow(chronicity_filtered2) # sensitivity = 77%



#predict on 6-year labels as if they were 2-year labels (overlap is only 64%)
testSetPrtm <- as.data.frame(whole_set[-trainIndex,-3])
testSetPrtm <- testSetPrtm %>% 
  #log10 transform as preprocessing step (not included as option in caret, so done manually here)
  mutate_at(vars(-c(1:2)), log10)

test6year <- inner_join(chronicity[,c(1,6)], testSetPrtm, by = "pident")
test6year$year6Chronicity <- factor(test6year$year6Chronicity, levels = c("remitted", "not_remitted"))

XGbPred6 <- predict(XGbTune, test6year[,-c(1:3)])
cm <- confusionMatrix(XGbPred6, test6year[,2])
draw_confusion_matrix(cm, ratio = TRUE)

#Get predicted class probabilities so we can build ROC curve.
XGbProbs6 <- predict(XGbTune, test6year[,-c(1:3)], type="prob")
head(XGbProbs6)

#build ROC curve
XGbROC6 <- roc(test6year[,2], XGbProbs6[,"remitted"])
auc(XGbROC6)

#plot ROC curve
ggplot(data.frame(specificity = XGbROC6$specificities, sensitivity = XGbROC6$sensitivities), aes(specificity, sensitivity)) +
  geom_path() +
  scale_x_reverse() +
  geom_abline(intercept = 1, slope = 1, colour='grey') +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
  theme_classic() +
  labs(title = paste0("AUROC =", signif(auc(XGbROC6), 3)))



















