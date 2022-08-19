##### Make datasets to be used for human chronicity predictions ##### ##### ##### ##### 
# - 4 human raters
# - for each rater: 
#     -40 samples with clinical data only (10 predictors)
#     -40 samples with clinical + additional clinial data (comorbidities, disease history, psychiatric diagnoses)
# - 40 samples will overlap (20 clinical, 20 clinical + additional data) for rater agreement analysis
# - in total 40 + 4*40 = 200 samples rated by humans. These 200 samples are sampled from the test sets in ML evaluation.
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

library(foreign)
library(tidyverse)
library(hashids)
library(kappaSize)
library(xlsx)

setwd("/Users/philippehabets/Dropbox/STRESS_INDEX/data/")

# Import additional CIDI Depression data
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
  select(pident, AD3004, acidep14)
names(cidiSelected) <- c("pident", "dysthymia", "MDD")
cidiSelected <- droplevels.data.frame(cidiSelected)
cidiSelected <- cidiSelected %>% 
  mutate(MDD_history = ifelse(MDD=="MDD recurrent","Previous MDD episode(s)",ifelse(MDD=="MDD first episode", "First MDD episode", "no"))) %>%
  mutate(dysthymia = ifelse(dysthymia == "Diag Positive", "Dysthymia", "no")) %>% 
  rename(dysthymia_lifetime = dysthymia) %>% 
  select(-MDD)

# Import additional CIDI Anxiety data
cidiADerived <- read.spss("psychiatry_and_somatic_variables/W1/CIDI/anxiety/derived.diagnosis.variables/N1_A259D.sav", to.data.frame = TRUE, use.missings = FALSE)
cidiAnx <- cidiADerived %>% 
  mutate(pident = as.factor(pident)) %>% 
  select(pident, aanxy01:aanxy05, aanxy16:aanxy20) %>% 
  dplyr::rename("social_phobia" = aanxy16, "panic_with_agoraphobia" = aanxy17, "panic_without_agoraphobia" = aanxy18, "agoraphobia" = aanxy19, "GAD" = aanxy20) %>% 
  mutate(social_phobia = ifelse(social_phobia=="Positive", "Social phobia", "no")) %>% 
  mutate(panic_with_agoraphobia = ifelse(panic_with_agoraphobia =="Positive", "Panic with agoraphobia", "no")) %>% 
  mutate(panic_without_agoraphobia = ifelse(panic_without_agoraphobia == "Positive", "Panic without agoraphobia", "no")) %>% 
  mutate(agoraphobia = ifelse(agoraphobia == "Positive", "Agoraphobia", "no")) %>% 
  mutate(GAD = ifelse(GAD == "Positive", "GAD", "no")) %>% 
  mutate(aanxy01 = ifelse(aanxy01=="Positive", "Social phobia", "no")) %>% 
  mutate(aanxy02 = ifelse(aanxy02 =="Positive", "Panic with agoraphobia", "no")) %>% 
  mutate(aanxy03 = ifelse(aanxy03 == "Positive", "Panic without agoraphobia", "no")) %>% 
  mutate(aanxy04 = ifelse(aanxy04 == "Positive", "Agoraphobia", "no")) %>% 
  mutate(aanxy05 = ifelse(aanxy05 == "Positive", "GAD", "no"))
  
# Import additional CIDI Alcohol disorder data
alcohol <- read.spss("psychiatry_and_somatic_variables/W1/CIDI/alcohol/N1_A203D.sav", to.data.frame = TRUE, use.missings = FALSE)
as.character.factor <- function(x){
  #function for turning factor values (levels) into character values
  as.character(levels(x))[x]
}
alcohol <- alcohol %>% 
  dplyr::rename(pident = Pident, alcohol_diagnosis_status_lifetime = alcversl, alcohol_abuse_recency = AD30500RE, alcohol_dependent_recency = AD30390RE) %>% 
  mutate(`alcohol_abuse_or_dependent_recency` = ifelse(is.na(alcohol_abuse_recency), 
                                         ifelse(is.na(alcohol_dependent_recency), "no alcohol abuse or dependency", as.character.factor(alcohol_dependent_recency)), 
                                         as.character.factor(alcohol_abuse_recency))) %>%  
  select(pident, alcohol_diagnosis_status_lifetime, alcohol_abuse_or_dependent_recency)

# Import chronicity data
chronicity.raw <- read.spss("/Users/philippehabets/Dropbox/STRESS_INDEX/data/psychiatry_and_somatic_variables/W1/Chronical-diseases/N1_A250R.sav", to.data.frame = TRUE, use.missings = FALSE)
as.character.factor <- function(x){
  #function for turning factor values (levels) into character values
  as.character(levels(x))[x]
}
to_no <- function(x){
  # turn cases containing "no" into "
  if(str_detect(tolower(x),"no")){
    x <- "no"
  } 
  x
}
chronicity <- chronicity.raw %>% select(Pident, acnsld02, aheart06:aheart09, aheart12n, adiabe01, adiabe03, astroke1, aarthr01, aarthr02, areuma02:areuma04, acanc01:acanc03, acanc05, ahypert1, ahypert2, aulcer1, aulcer2, aintes02:aintes04, aintes07, aliver02, aliver03, aliver06, aepilep1, acfs1, aallerg2, aallerg3, aallerg5n, athygl02:athygl04, aneuro02:aneuro04, aother01, aother06, aother11) %>% 
  mutate_if(is.factor, ~sapply(.x, as.character.factor)) %>% 
  mutate_all(~sapply(.x, to_no)) %>% 
  mutate_all(~sapply(.x, tolower))
chronicity <- chronicity %>% 
  mutate(aheart06 = ifelse(aheart06 == "yes", "coronary heart disease", "no")) %>% 
  mutate(aheart07 = ifelse(aheart07 == "yes", "angina pectoris", "no")) %>% 
  mutate(aheart08 = ifelse(aheart08 == "yes", "heart failure", "no")) %>% 
  mutate(aheart09 = ifelse(aheart09 == "yes", "cardiac arrhythmia", "no")) %>% 
  mutate(diabetes = ifelse(adiabe01 == "yes" & adiabe03 == "no", "diabetes without medication use", 
                           ifelse(adiabe01 == "yes" & adiabe03 == "yes, only tablets", "diabetes (only tablets for medication)", 
                                  ifelse(adiabe01 == "yes" & adiabe03 == "yes, injections with insuline", "insuline dependent diabetes", "no")))) %>% 
  mutate(astroke1 = ifelse(astroke1 == "yes", "CVA", "no")) %>% 
  mutate(arthritis_arthrose = ifelse(aarthr01 == "yes", paste0("arthritis or arthrosis from age ", aarthr02), "no")) %>% 
  mutate(areuma02 = ifelse(areuma02 == "mentioned", "SLE", "no")) %>% 
  mutate(areuma03 = ifelse(areuma03 == "mentioned", "fibromyalgie", "no")) %>% 
  mutate(areuma04 = ifelse(areuma04 == "mentioned", "reumatoide arthritis", "no")) %>% 
  mutate(cancer = ifelse(acanc01 == "yes" & acanc05 == "yes", paste(acanc03, "malignancy with metastasis, ", "age onset", acanc02),
                         ifelse(acanc01 == "yes" & acanc05 == "no", paste(acanc03, "malignancy, ", "age onset", acanc02), "no"))) %>% 
  mutate(hypertension = ifelse(ahypert1 == "yes", paste("hypertension under treatment (age onset ", ahypert2, ")", sep = ""), "no")) %>% 
  mutate(ulcer = ifelse(aulcer1 == "yes", paste0("stomach/intestinal ulcer at age ", aulcer2), "no")) %>% 
  mutate(aintes02 = ifelse(aintes02 == "mentioned", "IBS", "no")) %>% 
  mutate(aintes03 = ifelse(aintes03 == "mentioned", "obstipation", "no")) %>% 
  mutate(aintes04 = ifelse(aintes04 == "mentioned", "IBD: Crohn's disease or Colitis Ulcerosa", "no")) %>% 
  mutate(aintes07 = ifelse(aintes07 == "                                                  ", "no", aintes07)) %>% 
  mutate(aliver02 = ifelse(aliver02 == "mentioned", "liver cirrhosis", "no")) %>% 
  mutate(aliver03 = ifelse(aliver03 == "mentioned", "Hepatitis", "no")) %>% 
  mutate(aliver06 = ifelse(aliver06 == "                                        ", "no", paste0("lever: ", aliver06))) %>% 
  mutate(aepilep1 = ifelse(aepilep1 == "yes", "epilepsy", "no")) %>% 
  mutate(acfs1 = ifelse(acfs1 == "yes", "chronisch vermoeidheidssyndroom (CVS)", "no")) %>% 
  mutate(aallerg2 = ifelse(aallerg2 == "mentioned", "hay fever", "no")) %>% 
  mutate(aallerg3 = ifelse(aallerg3 == "mentioned", "eczema", "no")) %>% 
  mutate(athygl02 = ifelse(athygl02 == "mentioned", "Graves disease", "no")) %>% 
  mutate(athygl03 = ifelse(athygl03 == "mentioned", "hypothyroidism", "no")) %>% 
  mutate(athygl04 = ifelse(athygl04 == "mentioned", "hyperthyroidism (not Graves)", "no")) %>% 
  mutate(aneuro02 = ifelse(aneuro02 == "mentioned", "migraine or headaches", "no")) %>% 
  mutate(aneuro03 = ifelse(aneuro03 == "mentioned", "MS", "no")) %>% 
  mutate(aneuro04 = ifelse(aneuro04 == "mentioned", "neuropathy", "no")) %>% 
  mutate_at(vars(aother01, aother06, aother11), ~sapply(.x, str_replace, "other", "other disease not specified")) %>% 
  select(-c(aulcer1, aulcer2, ahypert1, ahypert2, acanc01:acanc05, aarthr01, aarthr02, adiabe01,adiabe03 )) %>% 
  mutate(Pident = as.factor(Pident)) %>% 
  rename(pident = Pident)
           
disease_history <- apply(chronicity[,-1], 1, paste, collapse=", ")
disease_history <- str_remove_all(disease_history, "no, ") %>% 
  str_remove_all(", no") %>% 
  str_squish()
disease_history <- chronicity %>% 
  mutate(disease_history = ifelse(disease_history == "no", "", disease_history)) %>% 
  select(pident, disease_history)

# merge all data into one df "additional_data"
additional_data <- inner_join(cidiSelected, cidiAnx, by = "pident") %>% 
  inner_join(alcohol, by = "pident") %>% 
  inner_join(disease_history, by = "pident")

# import data with 10 clinical variables, remove missings, add the additional data
setwd("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/human_vs_machine/data")
whole_set_clinical <- read_csv("whole_set_clinical.csv") %>% 
  drop_na() %>% 
  mutate(pident = as.factor(pident))
whole_set <- inner_join(whole_set_clinical, additional_data, by = "pident")

# tidy df a bit
whole_set$aidssev <- recode(whole_set$aidssev, "1"="none to mild", "2"="mild", "3"="moderate", "4"="severe", "5"="very severe") 
whole_set <- whole_set %>%
  mutate(Sexefemale = ifelse(Sexefemale == 1, "female", "male")) %>% 
  mutate(aidssev = as.character(aidssev)) 

# combine anxiety related variables into 2 variables (lifetime and 1-month recency
anxiety_lifetime <- apply((whole_set %>% select(social_phobia:GAD)), 1, paste, collapse=", ")
anxiety_lifetime <- str_remove_all(anxiety_lifetime, "no, ") %>% 
  str_remove_all(", no") %>% 
  str_squish()
anxiety_month <- apply((whole_set %>% select(aanxy01:aanxy05)), 1, paste, collapse=", ")
anxiety_month <- str_remove_all(anxiety_month, "no, ") %>% 
  str_remove_all(", no") %>% 
  str_squish()

whole_set_combined <- whole_set %>% 
  mutate(`anxiety diagnosis lifetime` = ifelse(anxiety_lifetime == "no", "no anxiety diagnosis", anxiety_lifetime), .before = alcohol_diagnosis_status_lifetime) %>% 
  mutate(`anxiety 1-month recency` = ifelse(anxiety_month == "no", "no 1-month recency anxiety", anxiety_month), .after = `anxiety diagnosis lifetime`) %>% 
  rename(gender = Sexefemale, 
         `years eduction` = aedu, 
         age = Age,
         `neuroticism (range: 12-60)` = aneurot, 
         `extraversion (range: 12-60)` = aextrave, 
         `openness (range: 12-60)` = aopenes, 
         `agreeableness (range: 12-60)` = aagreeab, 
         `conscientiousness (range: 12-60)` = aconscie, 
         `IDS (range: 0-84)` = aids, 
         `IDS categorized (none to mild, mild, moderate, severe, very severe)` = aidssev,
         `dysthymia diagnosis (lifetime)` = dysthymia_lifetime,
         `alcohol diagnosis status (lifetime)` = alcohol_diagnosis_status_lifetime,
         `alcohol abuse or dependent: recency` = alcohol_abuse_or_dependent_recency,
         `MDD history` = MDD_history,
         `disease history` = disease_history) %>% 
  select(-c(aanxy01:aanxy05, social_phobia:GAD))


###################
# create sets for raters that at least overlap with test sets used in ML models (clinical based and +proteomics based), remove imputed rows:
test_clinical <- read_csv("clinical_test_imputed.csv") %>% 
  mutate(pident = as.factor(pident)) %>% 
  select(pident, Remitted_depression) %>% 
  inner_join(whole_set_combined, by = c("pident", "Remitted_depression"))
test_proteomics_clinical<- read_csv("clinical_proteomics_test.csv") %>% 
  mutate(pident = as.factor(pident)) %>% 
  select(pident, Remitted_depression) %>% 
  inner_join(whole_set_combined, by = c("pident", "Remitted_depression"))
  
# check proportions of outcome in combined test sample (imputed rows removed):
test_combined <- full_join(test_clinical, test_proteomics_clinical, by = c("pident", "Remitted_depression")) 
outcome_proportions <- prop.table(table(test_combined$Remitted_depression)) # 48% vs 52%

# evaluate required sample size for 40 shared ratings (assume kappa = 0.6, with precision 0.2 on both sides)
# result: minimum of 32 --> 40 is enough
CIBinary(kappa0 = 0.6, 
         kappaL = 0.4, 
         kappaU = 0.8, 
         props = c(outcome_proportions[1], outcome_proportions[2]), 
         raters = 4, 
         alpha = 0.05) 

#### create rater reliability set: 

# use 20 samples of test_clinical for rating set 1 (10vars)
#sample from 10 samples not_remitted
set.seed(1)
set_reliability_1_cases <- test_clinical %>% 
  filter(Remitted_depression =="not_remitted") %>% 
  slice_sample(n = 10)

set.seed(0)
set_reliability_1_controls <- test_clinical %>% 
  filter(Remitted_depression =="remitted") %>% 
  slice_sample(n = 10)

#merge both samples 
set_reliability_1 <- bind_rows(set_reliability_1_cases, set_reliability_1_controls)
prop.table(table(set_reliability_1$Remitted_depression)) # not_remitted/remitted = 50% vs 50%

# write these to file
# write_csv(set_reliability_1, "set_reliability_1(clinical).csv")

# use 20 non-overlapping samples of test_additional for rating set 2 (10vars + additional vars) (also use balanced outcomes)
index_additional <- which(!test_proteomics_clinical$pident %in% test_clinical$pident) # index of samples from test_proteomics_clinical that are not already included in test_clinical and therefore also not in set_reliability_1 

#sample from the non-overlapping samples 10 not_remitted
set.seed(3)
set_reliability_2_cases <- test_proteomics_clinical[index_additional,] %>% 
  filter(Remitted_depression =="not_remitted") %>% 
  slice_sample(n = 10) 

#sample from the non-overlapping samples 10 remitted
set.seed(2)
set_reliability_2_controls <- test_proteomics_clinical[index_additional,] %>% 
  filter(Remitted_depression =="remitted") %>% 
  slice_sample(n = 10)

#merge both samples 
set_reliability_2 <- bind_rows(set_reliability_2_cases, set_reliability_2_controls)
prop.table(table(set_reliability_2$Remitted_depression)) # not_remitted/remitted = 50% vs 50%
# write these to file
# write_csv(set_reliability_2, "set_reliability_2(additional).csv")

### create human rater set clinical: 20 samples from reliability set + 4*20 other clinical test samples (balanced)
# filter samples of test_clinical not in test_additional and not in set_reliability 1:
unique_clinical <- anti_join(test_clinical, test_proteomics_clinical, by = "pident") %>% 
  anti_join(set_reliability_1, by = "pident")
nrow(unique_clinical) # >80 (117)

#sample from these non-overlapping samples 40 not_remitted
set.seed(0)
set_clinical_cases <- unique_clinical %>% 
  filter(Remitted_depression =="not_remitted") %>% 
  slice_sample(n = 40) 

#sample from the non-overlapping samples 40 remitted
set.seed(1)
set_clinical_controls <- unique_clinical %>% 
  filter(Remitted_depression =="remitted") %>% 
  slice_sample(n = 40)

#merge 20 reliability + 80 unique cases and control samples:
set_clinical <- bind_rows(set_reliability_1, set_clinical_cases, set_clinical_controls) %>% 
  select(pident:`IDS categorized (none to mild, mild, moderate, severe, very severe)`) # remove additional columns
prop.table(table(set_clinical$Remitted_depression)) # not_remitted/remitted = 50% vs 50%


### create human rater set clinical + additional data: 20 samples from reliability set 2 + 4*20 other additional test samples (balanced)
# additional samples not in test_clinical and not in set_reliability_2:
unique_additional <- anti_join(test_proteomics_clinical, test_clinical, by = "pident") %>% 
  anti_join(set_reliability_2, by = "pident")
nrow(unique_additional) # 73, so additional samples should be added (sample from unused samples)
table(unique_additional$Remitted_depression) # 32 not_remitted, 41 remitted -> sample 8 not_remitted cases from unused whole_set samples

# filter samples in test_clinical not in set_clinical, and not in test_additional:
extra_samples <- anti_join(whole_set_combined, set_clinical, by = "pident") %>% #not in 100 samples with clinical data
  anti_join(set_reliability_2, by = "pident") %>% # not in 20 samples used for reliability scoring in clinical+additional data
  anti_join(unique_additional, by = "pident") # not in test_proteomics_clinical already used for sampling

# Sample 8 not_remitted cases from these non-overlapping extra samples 
extra_cases <- extra_samples %>% 
  filter(Remitted_depression =="not_remitted") %>% 
  slice_head(n = 8) 

#sample from non-overlapping test_additional samples 32 not_remitted (only 32 available) and add extra cases
set.seed(0)
set_additional_cases <- unique_additional %>% 
  filter(Remitted_depression =="not_remitted") %>% 
  slice_sample(n = 32) %>%  # will select 32, nrow=32
  bind_rows(extra_cases)

#sample from non-overlapping test_additional samples 40 remitted
set.seed(1)
set_additional_controls <- unique_additional %>% 
  filter(Remitted_depression =="remitted") %>% 
  slice_sample(n = 40) 

#merge 20 reliability + 80 unique cases and control samples:
set_additional <- bind_rows(set_reliability_2, set_additional_cases, set_additional_controls)
prop.table(table(set_additional$Remitted_depression)) # not_remitted/remitted = 50% vs 50%


################################################################################################################
# encode pident to hash id for distribution to human raters
################################################################################################################
settings_list <- hashid_settings(salt = "eval", min_length = 10, alphabet = DEFAULT_ALPHABET, sep = DEFAULT_SEPS)

as.numeric.factor <- function(x){
  #function for turning factor values (levels) into numeric values
  as.numeric(levels(x))[x]
}

set_clinical <- set_clinical %>% 
  mutate(pident = sapply(pident, as.numeric.factor)) %>%  
  mutate(hashID = sapply(pident, encode, settings_list)) %>% 
  relocate(hashID, .before = pident) %>% 
  mutate(pident = as.factor(pident))
           
set_additional <- set_additional %>% 
  mutate(pident = sapply(pident, as.numeric.factor)) %>% 
  mutate(hashID = sapply(pident, encode, settings_list)) %>% 
  relocate(hashID, .before = pident) %>% 
  mutate(pident = as.factor(pident))


### distribute clinical set in 4 sets for 4 raters (first 20 are the same samples)
raters <- c("rater_1","rater_2", "rater_3", "rater_4")

sets_clinical <- list()
for(i in 1:length(raters)){
  start_index <- 1 + 20*i
  stop_index <- 20 + 20*i
  sets_clinical[[i]] <- bind_rows(set_clinical[1:20,], set_clinical[start_index:stop_index,]) %>% 
    select(-c(pident, Remitted_depression)) %>%  #remove pident and outcome for raters
    slice_sample(n=40) %>%  # shuffle rows so first 20 samples are not the same for each rater
    mutate(predicted_outcome = "") # add column to fill in predictions
  names(sets_clinical)[i] <- raters[i]
}

# write to file
for(i in 1:length(sets_clinical)){
  write.xlsx(as.data.frame(sets_clinical[[i]]), paste("rater_files/", names(sets_clinical)[i], "_clinical_set.xlsx", sep = ""), row.names = FALSE)
}

### distribute clinical+additional set in 4 sets for 4 raters (first 20 are the same samples)
sets_additional <- list()
for(i in 1:length(raters)){
  start_index <- 1 + 20*i
  stop_index <- 20 + 20*i
  sets_additional[[i]] <- bind_rows(set_additional[1:20,], set_additional[start_index:stop_index,]) %>% 
    select(-c(pident, Remitted_depression)) %>%  #remove pident and outcome for raters
    slice_sample(n=40) %>%  # shuffle rows so first 20 samples are not the same for each rater
    mutate(predicted_outcome = "") # add column to fill in predictions
  names(sets_additional)[i] <- raters[i]
}

# write to file
for(i in 1:length(sets_additional)){
  write.xlsx(as.data.frame(sets_additional[[i]]), paste("rater_files/", names(sets_additional)[i], "_additional_data_set.xlsx", sep = ""), row.names = FALSE)
}
















