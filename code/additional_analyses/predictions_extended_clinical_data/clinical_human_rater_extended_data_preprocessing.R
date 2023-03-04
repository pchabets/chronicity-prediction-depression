library(foreign)
library(tidyverse)
library(hashids)
library(kappaSize)
library(readxl)

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

# combine anxiety related variables into 2 variables (lifetime and 1-month recency
anxiety_lifetime <- apply((whole_set %>% select(social_phobia:GAD)), 1, paste, collapse=", ")
anxiety_lifetime <- str_remove_all(anxiety_lifetime, "no, ") %>% 
  str_remove_all(", no") %>% 
  str_squish()
anxiety_month <- apply((whole_set %>% select(aanxy01:aanxy05)), 1, paste, collapse=", ")
anxiety_month <- str_remove_all(anxiety_month, "no, ") %>% 
  str_remove_all(", no") %>% 
  str_squish()

whole_set <- whole_set %>%
  rename(female = Sexefemale, 
         `years eduction` = aedu, 
         age = Age,
         `neuroticism` = aneurot, 
         `extraversion` = aextrave, 
         `openness` = aopenes, 
         `agreeableness` = aagreeab, 
         `conscientiousness` = aconscie, 
         `IDS` = aids, 
         `IDS categorized` = aidssev,
         `dysthymia diagnosis` = dysthymia_lifetime,
         `alcohol diagnosis status (lifetime)` = alcohol_diagnosis_status_lifetime,
         `alcohol abuse or dependent: recency` = alcohol_abuse_or_dependent_recency,
         `MDD history` = MDD_history,
         `disease history` = disease_history)



# encode categoricals for model training and tidy variables 
whole_set <- whole_set%>%
  mutate(dysthymia_lifetime = ifelse(dysthymia_lifetime == "Dysthymia", 1, 0), # dysthymia lifetime yes = 1, no = 0
         
         # engineer feature MDD_history: no history = 0, First episode = 1, Previous episodes = 2
         MDD_history = ifelse(MDD_history == "no", 0, ifelse(MDD_history == "First MDD episode", 1, 2)),
         
         # engineer feature social_phobia: no lifetime diagnosis = 0, lifetime diagnosis = 1, lifetime diagnosis & 1 month recency = 2
         social_phobia = ifelse(social_phobia == "no", 0, ifelse(aanxy01 == "no", 1, 2)),
         
         # engineer feature panic_with_agoraphobia: no lifetime diagnosis = 0, lifetime diagnosis = 1, lifetime diagnosis & 1 month recency = 2
         panic_with_agoraphobia = ifelse(panic_with_agoraphobia == "no", 0, ifelse(aanxy02 == "no", 1, 2)),
         
         # engineer feature panic_without_agoraphobia: no lifetime diagnosis = 0, lifetime diagnosis = 1, lifetime diagnosis & 1 month recency = 2
         panic_without_agoraphobia = ifelse(panic_without_agoraphobia == "no", 0, ifelse(aanxy03 == "no", 1, 2)),
         
         # engineer feature agoraphobia: no lifetime diagnosis = 0, lifetime diagnosis = 1, lifetime diagnosis & 1 month recency = 2
         agoraphobia = ifelse(agoraphobia == "no", 0, ifelse(aanxy04 == "no", 1, 2)),
         
         # engineer feature GAD: no lifetime diagnosis = 0, lifetime diagnosis = 1, lifetime diagnosis & 1 month recency = 2
         GAD = ifelse(GAD == "no", 0, ifelse(aanxy05 == "no", 1, 2)),
         
         # engineer alcohol lifetime feature: no positive diagnosis = 0, positive diagnosis = 1
         alcohol_diagnosis_status_lifetime = str_trim(as.character.factor(alcohol_diagnosis_status_lifetime)),
         alcohol_diagnosis_status_lifetime = ifelse(alcohol_diagnosis_status_lifetime == "Diagnose alcohol dependent or alcohol abuse", 1, 0),
         
         # engineer alcohol recency feature: no abuse/recency = 0, >1yr = 1, 1yr = 2, 1yr-6m =3, 6m-1m = 4, 1m-2wks = 5, <2wks = 6
         alcohol_abuse_or_dependent_recency = str_trim(alcohol_abuse_or_dependent_recency),
         alcohol_abuse_or_dependent_recency = ifelse(alcohol_abuse_or_dependent_recency == "no alcohol abuse or dependency", 0,
                                                     ifelse(alcohol_abuse_or_dependent_recency == "More than 1 yr ago", 1, 
                                                            ifelse(alcohol_abuse_or_dependent_recency == "In last 12 mths", 2,
                                                                   ifelse(alcohol_abuse_or_dependent_recency == "6 mths - 1 yr ago", 3,
                                                                          ifelse(alcohol_abuse_or_dependent_recency == "1mth - 6 mths ago", 4,
                                                                                 ifelse(alcohol_abuse_or_dependent_recency == "2 wks - 1 mth ago", 5,
                                                                                        ifelse(alcohol_abuse_or_dependent_recency == "Within last 2 weeks", 6, NA)))))))) %>% 
  # engineer disease_history: disease history is count of number of diseases 
  mutate(disease_history_length = str_length(whole_set$disease_history)) %>%  
  mutate(disease_history_number = str_count(str_remove_all(whole_set$disease_history, ", only symptoms mentioned"),",")+1) %>% 
  mutate(disease_history = ifelse(disease_history_length == 0, 0, disease_history_number)) %>% 
  select(-c(aanxy01:aanxy05, disease_history_length, disease_history_number)) %>% 
  rename(female = Sexefemale, 
         `years eduction` = aedu, 
         age = Age,
         `neuroticism` = aneurot, 
         `extraversion` = aextrave, 
         `openness` = aopenes, 
         `agreeableness` = aagreeab, 
         `conscientiousness` = aconscie, 
         `IDS` = aids, 
         `IDS categorized` = aidssev,
         `dysthymia diagnosis` = dysthymia_lifetime,
         `alcohol diagnosis status (lifetime)` = alcohol_diagnosis_status_lifetime,
         `alcohol abuse or dependent: recency` = alcohol_abuse_or_dependent_recency,
         `MDD history` = MDD_history,
         `social phobia` = social_phobia,
         `panic with agoraphobia` = panic_with_agoraphobia,
         `panic without agoraphobia` = panic_without_agoraphobia,
         `disease history` = disease_history)

# setwd("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/additional_analyses_revision/")
# write_csv(whole_set, "data/clinical_human_rater_extended_data.csv") 
###################














