library(tidyverse)
setwd("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/additional_analyses_revision/")

data <- read_csv("data/transcriptomics_kegg_selection.csv")

# STEP 1: select for each gene the probe that shows the highest variance across individuals
highest_var <- data %>%
  select(-c(PID,plate)) %>%
  gather() %>%
  mutate(gene = str_split(key, "\\.") %>%
           map_chr(1)) %>%
  group_by(gene) %>%
  mutate(variance = var(value)) %>%
  slice_tail(n=1) %>%
  select(-variance)

selected_probes <- data %>%
  select(PID, plate, any_of(highest_var$key))

write_csv(selected_probes, "data/gene_expression_var_filtered_probes.csv")





