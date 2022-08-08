# Input list of genes for the search bar

library(tidyverse)

load("all_data_prioritization_report.RData")

biomarker_list %>% 
  select(ENTRY, `ENTRY NAME`, `PROTEIN NAMES`, `GENE NAMES`, `GENE SYNONYM`) %>% 
  write_tsv("uniprot_mapping.tsv")


dat <- read_tsv("../Biomarker_prioritisation_matrices/EFO_0000400_Master_Matrix.tsv")

dat %>% 
  select(`UniProt AC`, Ensembl, Gene, `Gene synonym`, `name (PHAROS)`, `Gene description (HPA)`) %>% 
  mutate(`Gene synonym` = str_split(`Gene synonym`, ", ")) %>% 
  unnest(`Gene synonym`)


dat %>% 
  select(`UniProt AC`, Gene, `Gene synonym`) %>% 
  mutate(`Gene synonym` = str_split(`Gene synonym`, ", ")) %>% 
  unnest(`Gene synonym`) %>% 
  pivot_longer(cols = -`UniProt AC`,
               names_to = "source",
               values_to = "gene_names") %>% 
  filter(!is.na(gene_names)) %>% 
  select(-source) %>% 
  distinct_all() %>% 
  group_by(`UniProt AC`) %>% 
  summarise(gene_names = paste(gene_names, collapse = " ")) %>% 
  write_tsv("gene_list.tsv")

