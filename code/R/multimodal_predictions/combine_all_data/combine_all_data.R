library(tidyverse)
library(plyr)
library(foreign)
library(caret)
library(doMC)

setwd("/Users/philippehabets/Dropbox/STRESS_INDEX/")

##create dataframe
w3Chronicity <- read.spss("data/outcome_groups/DSM_groups.sav", to.data.frame = TRUE, use.missings = FALSE)
w3Chronicity <- w3Chronicity[,c(1,3)]

prtm <- read_csv("data/blood_and_saliva_variables/W1/proteomics/output/proteomics_replaced_outliers.csv")
prtm <- prtm %>% 
  mutate(Pident = as.factor(Pident))

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

#transcriptomics
ge <- read_csv('data/blood_and_saliva_variables/W1/transcriptomics/transcriptomics_2ychronicity.csv') %>% 
  select(-c(plate, Remitted_depression)) %>% 
  mutate(pident = as.factor(pident))

# Merge to one frame
whole_set <- left_join(w3Chronicity, demographic, by = "pident") %>% 
  left_join(bigFive, by = "pident") %>% 
  left_join(IDS, by = "pident") %>% 
  inner_join(prtm, by = c("pident" = "Pident")) %>% #keeping only cases with proteomics measured
  mutate_if(is.factor, droplevels) %>% 
  mutate(applate = as.factor(applate)) %>% 
  mutate(aidssev = as.numeric(aidssev)) %>%  #turn ordinal data into numerical data
  inner_join(mtbl, by = 'pident') %>% 
  inner_join(prs, by = 'pident') %>% 
  inner_join(ge, by = 'pident')

options(na.action = "na.pass")
modelMatrix <- model.matrix(Remitted_depression ~ ., data = whole_set %>% select(Remitted_depression, Sexe))[,-1] #transform categorical nominal features into numeric values
female = as.data.frame(modelMatrix); colnames(female) <- "female"
whole_set <- bind_cols(pident = whole_set$pident, 
                       Remitted_depression = whole_set$Remitted_depression, 
                       applate = whole_set$applate, 
                       female = female, 
                       whole_set %>% select(-c(pident, Remitted_depression, applate, Sexe)))

whole_set$Remitted_depression <- revalue(whole_set$Remitted_depression, c("remitted"="remitted", "non-remitted"="not_remitted")) 
table(whole_set$Remitted_depression) 

##############################################################################################################
## Using heldout test set 
##############################################################################################################
## partition data
set.seed(42)
trainIndex <- createDataPartition(y=whole_set$Remitted_depression, times=1, p=0.8, list=F)

segTrain <- as.data.frame(whole_set[trainIndex,])
table(segTrain$Remitted_depression) # 319 vs 314

segTest <- as.data.frame(whole_set[-trainIndex,])
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

#impute clinical variables
pp <- preProcess(segTrain[, 4:13], method = c("medianImpute"))
transformedTrain <- predict(pp, newdata = segTrain[, 4:13])
segTrain[4:13] <- transformedTrain

#impute lipidomic variables
pp2 <- preProcess(segTrain[names(trainLipidMins)], method = c("zv", "bagImpute"))
transformedTrain2 <- predict(pp2, newdata = segTrain[names(trainLipidMins)])
segTrain[names(trainLipidMins)] <- transformedTrain2

## Test data:
# Replace -1 with minimum values of that analyte
segTest[names(trainLipidMins)] <- replaceMinCol(segTest[names(trainLipidMins)])
# Replace -2 with NA
segTest[names(trainLipidMins)] <- segTest[names(trainLipidMins)] %>% 
  mutate_if(is_any_NR, ~na_if(.x, -2)) 

#imputetest data based on traindata fit of preprocessing 
transformedTest <- predict(pp, newdata = segTest[, 4:13])
segTest[4:13] <- transformedTest #clinical variables

transformedTest2 <- predict(pp2, newdata = segTest[names(trainLipidMins)])
segTest[names(trainLipidMins)] <- transformedTest2 #lipidomics

#proteomic imputation and batch effect correction
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
  batch_median[[batch]] <- sapply(batch_list[[plates[batch]]][names(prtm)[-c(1:2)]], median.default, na.rm=T)
  batch_IQR[[batch]] <- sapply(batch_list[[plates[batch]]][names(prtm)[-c(1:2)]], IQR, na.rm=T)
}
names(batch_median) <- plates
names(batch_IQR) <- plates

batch_srs <- vector(mode = "list", length = length(plates))
for(batch in 1:length(plates)) {
  batch_srs[[batch]] <- setNames(data.frame(matrix(ncol = ncol(whole_set), nrow = nrow(batch_list[[batch]]))), colnames(whole_set))
  batch_srs[[batch]][1:13] <- batch_list[[batch]][1:13]
  for(col in 1:(ncol(prtm)-2)) { #-2 because first two columns are no analytes
    column <- batch_list[[batch]][col+13] # +13 because first thirteen columns are no analytes
    if(!all(is.na(column)) && batch_IQR[[batch]][col] !=0) { # keep constants
      column <- sapply(column, function(analyte) {
        analyte <- 1/(1+exp(-(analyte-(batch_median[[batch]][col]))/(batch_IQR[[batch]][col]/1.35)))
      })
    }
    batch_srs[[batch]][col+13] <- column
  }
  batch_srs[[batch]][(1+ncol(prtm)-2+13):ncol(batch_list[[batch]])] <- batch_list[[batch]][(1+ncol(prtm)-2+13):ncol(batch_list[[batch]])]
}
names(batch_srs) <- plates

# rescale to unit interval
batch_srs_min <- list()
batch_srs_max <- list()
for(batch in 1:length(plates)){
  batch_srs_min[[batch]] <- sapply(batch_srs[[plates[batch]]][names(prtm)[-c(1:2)]], min, na.rm=T)
  batch_srs_max[[batch]] <- sapply(batch_srs[[plates[batch]]][names(prtm)[-c(1:2)]], max, na.rm=T)
}
names(batch_srs_min) <- plates
names(batch_srs_max) <- plates

batch_srs_norm <- vector(mode = "list", length = length(plates))
for(batch in 1:length(plates)) {
  batch_srs_norm[[batch]] <- setNames(data.frame(matrix(ncol = ncol(whole_set), nrow = nrow(batch_srs[[batch]]))), colnames(whole_set))
  batch_srs_norm[[batch]][1:13] <- batch_srs[[batch]][1:13]
  for(col in 1:(ncol(prtm)-2)) { #-2 because first two columns are no analytes
    column <- batch_srs[[batch]][col+13] # +13 because first thirteen columns are no analytes
    if(!all(is.na(column)) && !is.na(batch_srs_max[[batch]][col]) && batch_srs_max[[batch]][col] != batch_srs_min[[batch]][col]) { # keep constants
      column <- sapply(column, function(analyte) {
        analyte <- (analyte-batch_srs_min[[batch]][col])/(batch_srs_max[[batch]][col]-batch_srs_min[[batch]][col])
      })
    }
    batch_srs_norm[[batch]][col+13] <- column
  }
  batch_srs_norm[[batch]][(1+ncol(prtm)-2+13):ncol(batch_srs[[batch]])] <- batch_srs[[batch]][(1+ncol(prtm)-2+13):ncol(batch_srs[[batch]])]
}
names(batch_srs_norm) <- plates

# transform back to single dataframe
segTrain <- bind_rows(batch_srs_norm) # or batch_srs if no rescaling is required
segTrain <- segTrain %>% mutate_at(vars(-c(1:13)), as.numeric)

# Because in some batches there was no variance and/or IQR, and min and max were equal, values were kept original here. Because either strategy to deal with this might result in preservation of batch information, 
# we drop these columns here.
no_IQR <- names(which(sapply(segTrain[names(prtm)[-c(1:2)]], function(x){ return(max(x, na.rm=T)>1)})))
segTrain <- segTrain %>% select(-no_IQR)

# store columns so we can select same columns in test set
trainCols <- colnames(segTrain) 

# make list of proteins kept for preProcessing pipeline
keptProt <- names(prtm)[-c(1,2)][!names(prtm)[-c(1,2)] %in% no_IQR]

## remove zero variance variables and impute missing variables
pp <- preProcess(segTrain[keptProt], method = c("zv", "bagImpute"))
transformedTrain <- predict(pp, newdata = segTrain[keptProt])
segTrain <- bind_cols(segTrain[,1:13], transformedTrain, segTrain %>% select(-c(1:13)) %>% select(-keptProt))

## shuffle rows back so no ordering on batch
shuffled <- as.data.frame(whole_set[trainIndex,])
segTrain <- inner_join(shuffled[,1:2], segTrain, by = c("pident","Remitted_depression"))

# Write to file
# write_csv(segTrain, "data/multimodal/wave1_omics_clinical_segTrain.csv")
# write_csv(segTrain, "data/multimodal/wave1_omics_clinical_norescale_segTrain.csv")

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
  batch_srs_y[[batch]][1:13] <- batch_list_y[[batch]][1:13]
  for(col in 1:(ncol(prtm)-2)) { #-2 because first two columns are no analytes
    column <- batch_list_y[[batch]][col+13] # +13 because first three columns are no analytes
    if(!all(is.na(column)) && batch_IQR[[batch]][col] !=0) { # keep constants
      column <- sapply(column, function(analyte) {
        analyte <- 1/(1+exp(-(analyte-(batch_median[[batch]][col]))/(batch_IQR[[batch]][col]/1.35)))
      })
    }
    batch_srs_y[[batch]][col+13] <- column
  }
  batch_srs_y[[batch]][(1+ncol(prtm)-2+13):ncol(batch_list_y[[batch]])] <- batch_list_y[[batch]][(1+ncol(prtm)-2+13):ncol(batch_list_y[[batch]])]
}
names(batch_srs_y) <- plates

# Rescale to unit interval using min and max of analyte from trainings data.
batch_srs_norm_y <- vector(mode = "list", length = length(plates))
for(batch in 1:length(plates)) {
  batch_srs_norm_y[[batch]] <- setNames(data.frame(matrix(ncol = ncol(whole_set), nrow = nrow(batch_srs_y[[batch]]))), colnames(whole_set))
  batch_srs_norm_y[[batch]][1:13] <- batch_srs_y[[batch]][1:13]
  for(col in 1:(ncol(prtm)-2)) { #-2 because first two columns are no analytes
    column <- batch_srs_y[[batch]][col+13] # +13 because first three columns are no analytes
    if(!all(is.na(column)) && !is.na(batch_srs_max[[batch]][col]) && batch_srs_max[[batch]][col] != batch_srs_min[[batch]][col]) { # keep constants
      column <- sapply(column, function(analyte) {
        analyte <- (analyte-batch_srs_min[[batch]][col])/(batch_srs_max[[batch]][col]-batch_srs_min[[batch]][col])
      })
    }
    batch_srs_norm_y[[batch]][col+13] <- column
  }
  batch_srs_norm_y[[batch]][(1+ncol(prtm)-2+13):ncol(batch_srs_y[[batch]])] <- batch_srs_y[[batch]][(1+ncol(prtm)-2+13):ncol(batch_srs_y[[batch]])]
}
names(batch_srs_norm_y) <- plates

# transform back to single dataframe
segTest <- bind_rows(batch_srs_norm_y) # or batch_srs_y if no rescaling required
segTest <- segTest %>% mutate_at(vars(-c(1:13)), as.numeric)

# remove columns with only NA values in trainset
segTest <- segTest %>% select(all_of(trainCols))

## remove zero variance variables and impute missing variables based on trainingsdata
transformedTest <- predict(pp, newdata = segTest[keptProt])
segTest <- bind_cols(segTest[,c(1:13)], transformedTest, segTest %>% select(-c(1:13)) %>% select(-keptProt))

# function to truncate x < 0 and x > 1 to 0 and 1 respectively if rescaled (this is due to scaling with minMax based on traindata)
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
segTest <- segTest %>% mutate_at(vars(keptProt), truncate) 

# Write to file
# write_csv(segTest, "data/multimodal/wave1_omics_clinical_segTest.csv")
# write_csv(segTest, "data/multimodal/wave1_omics_clinical_norescale_segTest.csv")



########## Alternatively: Feature pre-selection for transcriptomics before input in Boruta feature selection of all features#############
segTrain <- read_csv("data/multimodal/wave1_omics_clinical_segTrain.csv")
segTest <- read_csv("data/multimodal/wave1_omics_clinical_segTest.csv")
ge <- read_csv('data/blood_and_saliva_variables/W1/transcriptomics/transcriptomics_2ychronicity.csv') %>% 
  select(-c(plate, Remitted_depression)) %>% 
  mutate(pident = as.factor(pident))
probes <- colnames(ge)[-1]; ge <- NULL

# reduce with Mann-Whitney U test (cutoff p <= 0.05)
rem <- segTrain %>% 
  filter(Remitted_depression == "remitted") %>% 
  select(probes) %>% 
  as.data.frame()

n_rem <- segTrain %>% 
  filter(Remitted_depression == "not_remitted") %>% 
  select(probes) %>% 
  as.data.frame()

df <- data.frame("probe" = character(), "p_value" = numeric(), stringsAsFactors = FALSE)
for(probe in 1:ncol(rem)) {
  test <- wilcox.test(rem[,probe], n_rem[,probe])
  df <- rbind(df, c(colnames(rem)[probe], test$p.value), stringsAsFactors = FALSE)
}
colnames(df) <- c("probe", "p_value")

df <- df %>% 
  arrange(p_value)

top_probes <- df$probe[1:250]

# remove non-top250 probes
train_select <- segTrain %>% 
  select(-(probes[!probes %in% top_probes]))
test_select <- segTest %>% 
  select(-(probes[!probes %in% top_probes]))

# write to file
# write_csv(train_select, "data/multimodal/wave1_omics_clinical_segTrain_transcriptomics_preselected.csv")
# write_csv(test_select, "data/multimodal/wave1_omics_clinical_segTest_transcriptomics_preselected.csv")





