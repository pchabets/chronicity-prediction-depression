library(tidyverse)
library(readxl)
library(hashids)
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
raters <- c("rater_1","rater_2", "rater_3", "rater_4")
ratings <- list()
for(i in 1:length(raters)){
  path_clin <- paste("rater_files/received/", raters[i], "_clinical_set.xlsx", sep = "")
  path_add <- paste("rater_files/received/", raters[i], "_additional_data_set.xlsx", sep = "")
  
  clinical <- read_xlsx(path_clin)
  additional <- read_xlsx(path_add)
  
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
accs <- data.frame(rater = c("xgboost", raters),
                   `clinical data only` = c(.63,
                                         accuracy(ratings$rater_1$clinical_data), 
                                         accuracy(ratings$rater_2$clinical_data), 
                                         accuracy(ratings$rater_3$clinical_data), 
                                         accuracy(ratings$rater_4$clinical_data)),
                   `with added data` = c(.71,
                                         accuracy(ratings$rater_1$additional_data), 
                                         accuracy(ratings$rater_2$additional_data), 
                                         accuracy(ratings$rater_3$additional_data), 
                                         accuracy(ratings$rater_4$additional_data))
)

# plotting of results
plot_accs <- melt(accs, variable.name = "used data", value.name = "balanced accuracy") %>% 
  mutate(rater = recode(rater, xgboost = "xgboost", rater_1 = "rater 1", rater_2 = "rater 2", rater_3 = "rater 3", rater_4 = "rater 4")) %>% 
  mutate(`used data` = recode(`used data`, clinical.data.only = "10 clinical features", with.added.data = "with added data"))
plot <- ggplot(plot_accs, aes(x = `balanced accuracy`, y = rater, fill = `used data`)) +
  geom_point(aes(color=`used data`), size = 6) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10),limits = c(0.3, 0.8)) +
  geom_vline(xintercept=c(.63,.71), linetype="dotted") +
  annotate("text", x = c(0.63,0.715), y=5.2, label = c("0.63", "0.71"), size=6) +
  theme_gray(base_size = 24) +
  theme(legend.justification = 'top')
plot

# ggsave("/Users/philippehabets/Dropbox/STRESS_INDEX/manuscripts/prediction_2ychronicity/Figures/output/human_vs_machine.pdf",
#        plot,
#        width = 35, height = 15, units = "cm")

# average score for every rater
apply(accs[,-1], 1, mean)

# average score per type of data
apply(accs[,-1], 2, mean)


################################################################################################################
# check interrater reliability
################################################################################################################
ratings$rater_1$all <- bind_rows(ratings$rater_1$clinical_data, ratings$rater_1$additional_data) %>% 
  select(-Remitted_depression)
ratings$rater_2$all <- bind_rows(ratings$rater_2$clinical_data, ratings$rater_2$additional_data) %>% 
  select(-Remitted_depression)
ratings$rater_3$all <- bind_rows(ratings$rater_3$clinical_data, ratings$rater_3$additional_data) %>% 
  select(-Remitted_depression)
ratings$rater_4$all <- bind_rows(ratings$rater_4$clinical_data, ratings$rater_4$additional_data) %>% 
  select(-Remitted_depression)

interRater <- inner_join(ratings$rater_1$all, ratings$rater_2$all, by = "pident", suffix = c(".rater_1", ".rater_2")) %>% 
  inner_join(ratings$rater_3$all, by = "pident") %>% 
  inner_join(ratings$rater_4$all, by = "pident", suffix = c(".rater_3", ".rater_4"))

kappa <- kappam.fleiss(interRater[,-1], exact = F, detail = T)


 
















