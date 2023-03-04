library(tidyverse)
setwd("/Users/philippehabets/Dropbox/STRESS_INDEX/data/")

# load transcriptomics 
load("blood_and_saliva_variables/W1/transcriptomics/workspace_2_v2_philipe _select.Rdata")
load("blood_and_saliva_variables/W1/transcriptomics/U219_probeset_info_new_v2")

rownames(expr) <- sapply(ens_gene_names, function(x) ifelse(x=="", "no_annotation", x))
colnames(expr) <- PID
expr <- t(expr)
expr <- cbind(PID, plate, expr) 
expr <- data.frame(expr, check.names = TRUE)
expr$PID <- as.factor(expr$PID)

# load KEGG categories genes data and extract ensg ids
kegg_genes_df <- read_csv("./../scripts/Predictions_explorative/R.scripts/2-yearChronicity/additional_analyses_revision/data/KEGG_categories_genes.csv") 
kegg_genes <- unique(kegg_genes_df %>% 
                       pull(var = ENSEMBL))

# # load KEGG categories genes data and extract symbols and aliases 
# kegg_genes_df <- read_csv("./../scripts/Predictions_explorative/R.scripts/2-yearChronicity/additional_analyses_revision/data/KEGG_categories_genes.csv") 
# kegg_genes <- kegg_genes_df %>% 
#   pull(var = SYMBOL) %>% 
#   c(kegg_genes_df %>% pull(ALIAS))
# kegg_genes <- unique(kegg_genes)

# keep only KEGG selected genes and write to file
selection <- expr %>% 
  select(PID, plate, any_of(starts_with(kegg_genes)))
length(unique(str_remove_all(colnames(selection)[-c(1:2)], "[.].*"))) # 3738 unique genes, with a total of 9771 probes

write_csv(selection, "./../scripts/Predictions_explorative/R.scripts/2-yearChronicity/additional_analyses_revision/data/transcriptomics_kegg_selection.csv")











