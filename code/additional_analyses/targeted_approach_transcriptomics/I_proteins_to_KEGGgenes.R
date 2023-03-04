library(tidyverse)
library(enrichR)
library(KEGGREST)
library(org.Hs.eg.db)

setwd("/Users/philippehabets/Dropbox/STRESS_INDEX/scripts/Predictions_explorative/R.scripts/2-yearChronicity/additional_analyses_revision/")

# load in uniprotIDs of proteomics platform matched to genes
mapping <- read_csv("data/idmapping.csv")
mapped_genes <- mapping %>% 
  pull(var=gene_symbol)

# Get enriched KEGG pathways associated with the proteins, and then get all genes associated with these pathways
setEnrichrSite("Enrichr") # Human genes
websiteLive <- TRUE
dbs <- listEnrichrDbs()
# use KEGG 2021 human database (gene coverage = 8078)
dbs <- "KEGG_2021_Human"
if (websiteLive) {
  enriched <- enrichr(mapped_genes, dbs)
}

# only keep terms with q-value <= 0.05
enriched_terms <- lapply(enriched, function(dataframe){return(dataframe %>% filter(Adjusted.P.value <= 0.05))})

# Get KEGG pathways and associated genes per pathway
hsa_path_eg  <- keggLink("pathway", "hsa") %>% 
  tibble(pathway = ., entrez_id = sub("hsa:", "", names(.)))
hsa_pathways <- keggList("pathway", "hsa") %>% 
  tibble(pathway = names(.), description = .)
genes_pathways <- left_join(hsa_pathways, hsa_path_eg, by="pathway") %>% 
  # remove " - Homo sapiens (human)" substring to match with enrichR terms
  mutate(description = str_remove_all(str_trim(description, side="both"), " - Homo sapiens \\(human\\)"))

# Inner join to get only enriched pathways and their associated genes
selected_pathways <- inner_join(genes_pathways, enriched_terms[[1]], by = c("description"="Term")) 

# Get list of unique entrez_ids 
genes <- unique(selected_pathways$entrez_id) #n = 4241

# convert to symbols and aliases for matching with transcriptomics data
keytypes(org.Hs.eg.db)
selected_genes <- select(org.Hs.eg.db, keys=genes,
                       columns=c("ENSEMBL", "SYMBOL", "ALIAS", "GENENAME"), keytype="ENTREZID")

write_csv(selected_genes, "data/KEGG_categories_genes.csv")








