library(tidyverse)
library(hashids)
library(xlsx)
library(irr)
library(reshape2)

setwd("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/human_vs_machine/data/")

################################################################################################################
# read in data, decode pident and check accuracy
################################################################################################################
# function to calculate accuracy:
accuracy <- function(df) {
  acc <- sum(df$predicted_outcome == df$Remitted_depression)/nrow(df)
  return(acc)
}
settings_list <- hashid_settings(salt = "eval", min_length = 10, alphabet = DEFAULT_ALPHABET, sep = DEFAULT_SEPS)

whole_set <- read_csv("whole_set_clinical.csv") %>% 
  mutate(pident = as.factor(pident))

# import data 
raters <- c("Christiaan","Joeri", "Jurjen", "Wouter")
ratings <- list()
for(i in 1:length(raters)){
  path_clin <- paste("rater_files/received/", raters[i], "_clinical_set.xlsx", sep = "")
  path_add <- paste("rater_files/received/", raters[i], "_additional_data_set.xlsx", sep = "")
  
  clinical <- read.xlsx(path_clin, sheetIndex = 1)
  additional <- read.xlsx(path_add, sheetIndex = 1)
  
  clinical_rating <- clinical %>% 
    mutate(pident = sapply(hashID, decode, settings_list)) %>% 
    relocate(pident, .before = hashID) %>% 
    mutate(pident = as.factor(pident)) %>% 
    inner_join(whole_set, by = "pident") %>% 
    select(pident, predicted_outcome, Remitted_depression)
  
  additional_rating <- additional %>% 
    mutate(pident = sapply(hashID, decode, settings_list)) %>% 
    relocate(pident, .before = hashID) %>% 
    mutate(pident = as.factor(pident)) %>% 
    inner_join(whole_set, by = "pident") %>% 
    select(pident, predicted_outcome, Remitted_depression)
  
  ratings[[i]] <- list(clinical_rating, additional_rating)
  names(ratings[[i]]) <- c("clinical_data", "additional_data")
  
  names(ratings)[i] <-raters[i]
}

# calculate accuracies
for(i in 1:length(raters)){
  accs <- lapply(ratings[[i]], accuracy)
  print(raters[i])
  print(accs)
}

# dataframe of accuracies:
accs <- data.frame(rater = raters,
                   `clinical data only` = c(accuracy(ratings$Christiaan$clinical_data), 
                                         accuracy(ratings$Joeri$clinical_data), 
                                         accuracy(ratings$Jurjen$clinical_data), 
                                         accuracy(ratings$Wouter$clinical_data)),
                   `with added data` = c(accuracy(ratings$Christiaan$additional_data), 
                                         accuracy(ratings$Joeri$additional_data), 
                                         accuracy(ratings$Jurjen$additional_data), 
                                         accuracy(ratings$Wouter$additional_data))
)

# plotting of results
plot_accs <- melt(accs, variable.name = "data", value.name = "accuracy") %>% 
  mutate(rater = recode(rater, Christiaan = "rater 1", Joeri = "rater 2", Jurjen = "rater 3", Wouter = "rater 4"))
ggplot(plot_accs, aes(x = accuracy, y = rater, fill = data)) +
  geom_dotplot(binaxis = "y", stackgroups = TRUE, binwidth = 0.1, method = "histodot")

# average score for every rater
apply(accs[,-1], 1, mean)

# average score per type of data
apply(accs[,-1], 2, mean)


################################################################################################################
# check interrater reliability
################################################################################################################
ratings$Christiaan$all <- bind_rows(ratings$Christiaan$clinical_data, ratings$Christiaan$additional_data) %>% 
  select(-Remitted_depression)
ratings$Joeri$all <- bind_rows(ratings$Joeri$clinical_data, ratings$Joeri$additional_data) %>% 
  select(-Remitted_depression)
ratings$Jurjen$all <- bind_rows(ratings$Jurjen$clinical_data, ratings$Jurjen$additional_data) %>% 
  select(-Remitted_depression)
ratings$Wouter$all <- bind_rows(ratings$Wouter$clinical_data, ratings$Wouter$additional_data) %>% 
  select(-Remitted_depression)

interRater <- inner_join(ratings$Christiaan$all, ratings$Joeri$all, by = "pident", suffix = c(".Christiaan", ".Joeri")) %>% 
  inner_join(ratings$Jurjen$all, by = "pident") %>% 
  inner_join(ratings$Wouter$all, by = "pident", suffix = c(".Jurjen", ".Wouter"))

kappa <- kappam.fleiss(interRater[,-1], exact = F, detail = T)


 
















