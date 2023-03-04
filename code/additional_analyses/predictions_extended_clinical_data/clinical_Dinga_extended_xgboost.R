library(tidyverse)
library(plyr)
library(foreign)
library(caret)
library(xgboost)
library(SHAPforxgboost)
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
  text(10, 35, names(cm$overall[1]), cex=1.5, font=2)
  text(10, 20, round(as.numeric(cm$overall[1]), 3), cex=1.4)
  text(50, 35, names(cm$byClass[11]), cex=1.5, font=2)
  text(50, 20, round(as.numeric(cm$byClass[11]), 3), cex=1.4)
  text(90, 35, names(cm$overall[2]), cex=1.5, font=2)
  text(90, 20, round(as.numeric(cm$overall[2]), 3), cex=1.4)
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
  select(pident, aids) %>% 
  mutate(aids = as.character.factor(aids)) %>% 
  mutate(aids = sapply(aids, removeAddition)) %>% 
  mutate(aids = as.numeric(aids))

#CIDI Depression (n = 2981)
cidi_raw <- read.spss("psychiatry_and_somatic_variables/W1/CIDI/depression/raw.scores+diagnoses/N1_A256D.sav", to.data.frame = TRUE, use.missings = FALSE)
cidiDerived_raw <- read.spss("psychiatry_and_somatic_variables/W1/CIDI/depression/derived.diagnosis.variables/N1_A257D.sav", to.data.frame = TRUE, use.missings = FALSE)
as.numeric.factor <- function(x){
  #function for turning factor values (levels) into numeric values
  as.numeric(levels(x))[x]
}
cidi <- cidi_raw %>% 
  mutate_at(vars(AD2962xAO, AD2963xAO, AD3004AO), ~sapply(.x, as.numeric.factor)) %>% 
  rowwise() %>% 
  mutate(MDDxAO = if_else(is.na(AD2962xAO) & is.na(AD2963xAO), NaN, min(c(AD2962xAO,AD2963xAO), na.rm = TRUE))) %>% #if two Age of Onsets available: choose youngest age
  select(Pident, AD3004, ANDPBOXSX, AD3004AO, MDDxAO, ADPCORE1, ADPCORE2, ADPBOX1, ADPBOX2, ADPBOX3, ADPBOX4, ADPBOX5, ADPBOX7, ADPBOX8)
cidiDerived <- cidiDerived_raw %>% 
  select(pident, acidep02, acidep04, acidep06, acidep08, acidep03, acidep05, acidep07, acidep09, acidep13, acidep14) %>% 
  mutate(pident = as.factor(pident))
cidiAll <- inner_join(cidiDerived, cidi, by = c("pident"= "Pident"))
cidiSelected <- cidiAll %>% 
  select(pident, AD3004, acidep14, ANDPBOXSX) %>% 
  mutate_at(vars(ANDPBOXSX), ~sapply(.x, as.numeric.factor))
names(cidiSelected) <- c("pident", "dysthymia", "MDD_history", "CIDI_yes_boxes")
cidiSelected <- droplevels.data.frame(cidiSelected) %>% 
  mutate(MDD_history = ifelse(MDD_history == "MDD recurrent", "Has had depressive episode previously", ifelse(MDD_history == "MDD first episode", "This is the first episode", "No MDD")))

#CIDI Anxiety (n = 2981)
cidiA_raw <- read.spss("psychiatry_and_somatic_variables/W1/CIDI/anxiety/raw.scores+diagnoses/N1_A258D.sav", to.data.frame = TRUE, use.missings = FALSE)
cidiADerived_raw <- read.spss("psychiatry_and_somatic_variables/W1/CIDI/anxiety/derived.diagnosis.variables/N1_A259D.sav", to.data.frame = TRUE, use.missings = FALSE)
cidiA <- cidiA_raw %>% 
  select(Pident, AD30023, AD30022, AD30001, AD30021, AD30002)
cidiADerived <- cidiADerived_raw %>% 
  mutate(pident = as.factor(pident)) %>% 
  select(pident, aanxy21, aanxy01, aanxy06, aanxy11, aanxy16, aanxy04, 
         aanxy09, aanxy14, aanxy19, aanxy03, aanxy08, aanxy13, aanxy18, 
         aanxy02, aanxy07, aanxy12, aanxy17, aanxy05, aanxy10, aanxy15, aanxy20)
cidiA_all <- inner_join(cidiA, cidiADerived, by = c("Pident" = "pident"))
not_posORneg <- function(x){
  #function to see percentages of diagnosis values that are "Diag Indeterminate" or "Diag Pos Except Excl" in each AD300** column 
  sum(x=="Diag Pos Except Excl" | x == "Diag Indeterminate")/length(x)*100
}
sapply(cidiA_all, not_posORneg)
binarize_cidiA <- function(x){ 
  # function to binarize factors in cidiA
  ifelse(x == "Negative" | x == "Diag Negative", 0, 
         ifelse(x == "Positive" | x == "Diag Postive", 1, NA))
}
cidiA <- cidiA_all %>% 
  mutate_at(vars(-c(Pident, aanxy21)), ~as.character.factor(.x)) %>% 
  mutate_at(vars(-c(Pident, aanxy21)), ~binarize_cidiA(.x)) %>% 
  drop_na() %>% 
  dplyr::rename(pident=Pident)

# BAI data
BAI <- read_csv("psychological_and_lifestyle_variables/W1/BAI/BAI_processed.csv") %>% 
  mutate(pident = as.factor(pident)) %>% 
  select(pident, abaiscal)

# CT data Nemesis W1
CT_raw <- read.spss("psychiatry_and_somatic_variables/W1/NEMESIS(ChildhoodTrauma)/N1_A252D.sav", to.data.frame = TRUE, use.missings = FALSE)
as.character.factor <- function(x){
  #function for turning factor values (levels) into character values
  as.character(levels(x))[x]
}
turnNA_CT <- function(x){
  # turn missings into NA
  if(x %in% c("Questionnaire not used in interview","Too many missings")){
    x <- NA
  } 
  x
}
CT <- CT_raw %>% 
  mutate(Pident = as.factor(Pident)) %>% 
  select(Pident, ACTI_total) %>% 
  mutate(ACTI_total = as.numeric(sapply(as.character.factor(ACTI_total), turnNA_CT))) %>% 
  dplyr::rename(pident=Pident)

# Family history (taken at W6)
FH_raw <- read.spss("demographic_and_sociological_variables/W6/Family_inventory/N1_F121R.sav", to.data.frame = TRUE, use.missings = FALSE)
nrow(FH_raw) #2069 participants available
FH <- FH_raw %>% 
  select(pident, fftimo_g05, fftifa_g05, fftimo_d01, fftimo_a01, fftimo_a02, fftifa_d01, fftifa_a01, fftifa_a02) %>% 
  mutate_at(vars(-pident), ~as.character.factor(.x)) %>% 
  # one hot encode maternal and paternal, with "Asked, no answer" coded as 0,0,0
  mutate("maternal_pos" = ifelse(fftimo_g05 == "Ja", 1, 0),
         "maternal_unk" = ifelse(fftimo_g05 == "Weet niet", 1, 0),
         "maternal_neg" = ifelse(fftimo_g05 == "Nee", 1, 0),
         "paternal_pos" = ifelse(fftifa_g05 == "Ja", 1, 0),
         "paternal_unk" = ifelse(fftifa_g05 == "Weet niet", 1, 0),
         "paternal_neg" = ifelse(fftifa_g05 == "Nee", 1, 0)) %>% 
  select(pident, maternal_pos:paternal_neg)
         
#merge into dataframe
whole_set <- left_join(w3Chronicity, demographic, by = "pident") %>% 
  left_join(bigFive, by = "pident") %>%
  left_join(IDS, by = "pident") %>% 
  left_join(cidiSelected, by = "pident") %>% 
  left_join(cidiA, by = "pident") %>% 
  left_join(BAI, by = "pident") %>% 
  left_join(CT, by = "pident") %>% 
  left_join(FH, by = "pident") %>% 
  mutate_if(is.factor, droplevels) %>% 
  # feature engineer sex, dysthymia and mdd history
  mutate(female = ifelse(Sexe == "male", 0, 1),
         dysthymia = ifelse(dysthymia == "Diag Negative", 0 ,1),
         MDD_history = ifelse(MDD_history == "No MDD", 0, 
                              ifelse(MDD_history == "This is the first episode", 1, 0))) %>% 
  select(-Sexe)

whole_set$Remitted_depression <- revalue(whole_set$Remitted_depression, c("remitted"="remitted", "non-remitted"="not_remitted")) 
table(whole_set$Remitted_depression) # 407 (58%) remitted, 397 (42%) non-remitted

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

summary(segTrain$Remitted_depression) #remitted/non-remitted = 326/318
summary(segTest$Remitted_depression) #remitted/non-remitted = 81/79

## Fit model
## set seed list for reproduction
set.seed(42)
seeds <- vector(mode = "list", length = 101)
for(i in 1:100) seeds[[i]] <- sample.int(1000, 1000)
seeds[[101]] <- sample.int(100,1)

train_ctrl <- trainControl(method="repeatedcv",
                           number = 10,
                           repeats = 10,
                           summaryFunction = twoClassSummary,
                           search = "random",
                           classProbs = TRUE,
                           # seeds = seeds,
                           allowParallel = TRUE,
                           verboseIter = TRUE)

adaptControl <- trainControl(method = "adaptive_cv",
                             number = 10, repeats = 10,
                             adaptive = list(min = 10, alpha = 0.05, method = "gls", complete = TRUE), #using adaptive resampling to have more efficient computation, see http://topepo.github.io/caret/adaptive-resampling.html
                             classProbs = TRUE,
                             summaryFunction = twoClassSummary,
                             preProcOptions = list(knnSummary = median, k=5),
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
tic()
XGbTune <- train(x = segTrain[,-1],
                 y = segTrain[,1],
                 method = "xgbTree",
                 tuneLength = 1000,
                 metric = "ROC",
                 trControl = train_ctrl,
                 preProcess = "knnImpute",
                 verbose = TRUE)
toc_hour()
print("save model!")

# saveRDS(XGbTune, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/XGbTune_clinical_Dinga_extended.RDS")
XGbTune <- readRDS("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/XGbTune_clinical_Dinga_extended.RDS")

# saveRDS(XGbTune, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/XGbTune_clinical_Dinga_extended_noBigFive.RDS")
# XGbTune <- readRDS("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/XGbTune_clinical_Dinga_extended_noBigFive.RDS")

XGbTune$finalModel

## predictions on test set:
XGbPred <- predict(XGbTune, segTest[,-1])
cm <- confusionMatrix(XGbPred, segTest[,1])
draw_confusion_matrix(cm, ratio = TRUE)

accuracies <- list()
for(i in seq(0.4, 0.6, 0.001)){
  pred_ <- factor(ifelse(predict(XGbTune, segTest[,-1], type="prob")$remitted>i,"remitted","not_remitted"), levels = c("not_remitted", "remitted"))
  cm_ <- confusionMatrix(pred_, segTest %>% pull(Remitted_depression) %>% as.factor())
  accuracies <- append(accuracies, cm_$overall[1])
}
names(accuracies) <- seq(0.4, 0.6, 0.001)
opt_thresh <- as.numeric(names(accuracies[accuracies == max(sapply(accuracies, max))]))[1]
XGbPred <- factor(ifelse(predict(XGbTune, segTest[,-1], type="prob")$remitted>opt_thresh,"remitted","not_remitted"), levels = c("not_remitted", "remitted"))
cm <- confusionMatrix(XGbPred, as.factor(segTest$Remitted_depression))
draw_confusion_matrix(cm, ratio = TRUE)

var_imp <- varImp(XGbTune, scale = FALSE)
var_imp$importance
analyte_names <- read_csv("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/analyte_names.csv")

analyte_importance <- var_imp$importance %>% 
  mutate(symbol = str_sub(rownames(var_imp$importance), start = 3)) %>% 
  left_join(y = analyte_names, by = c("symbol" = "Analyte abbreviation")) %>% 
  dplyr::rename(Gain = Overall) %>% 
  mutate(Gain = round(Gain, 5)) %>% 
  select(2,3,1)
# write_excel_csv(analyte_importance, "/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/data/analyte_importance_XGB-Gain.csv")


#Get predicted class probabilities so we can build ROC curve.
XGbProbs <- predict(XGbTune, segTest[,-1], type="prob")
head(XGbProbs)

#build ROC curve
XGbROC <- roc(segTest[,1], XGbProbs[,"remitted"])
auc(XGbROC)
# save to file
saveRDS(XGbROC, "./../scripts/Predictions_explorative/R.scripts/data/roc_2ychr_dinga_extended_clinical.RDS")


#plot ROC curve
ggplot(data.frame(specificity = XGbROC$specificities, sensitivity = XGbROC$sensitivities), aes(specificity, sensitivity)) +
  geom_path() +
  scale_x_reverse() +
  geom_abline(intercept = 1, slope = 1, colour='grey') +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
  theme_classic() +
  labs(title = paste0("AUROC =", signif(auc(XGbROC), 2)))


