# Prep idmapping file
library(tidyverse)
library(vroom)
library(janitor)

# directory
data_dir <- here::here("data/")


# UniProt ID and ensembl mapping

# Human data
human_uniprot_data <- read_tsv(str_c(data_dir, "UniProt_mapping_files/", "HUMAN_9606_idmapping.dat"), col_names = c("uniprot", "id_type", "id"))

human_uniprotkb <- human_uniprot_data %>% 
  filter(id_type %in% c("UniProtKB-ID")) %>% 
  select(uniprot, uniprotkb_id = id)

human_gene_name <- human_uniprot_data %>% 
  filter(id_type %in% c("Gene_Name")) %>% 
  select(uniprot, gene_name = id)

human_ensembl <- human_uniprot_data %>% 
  filter(id_type %in% c("Ensembl"))  %>%
  select(uniprot, ensembl = id) %>% 
  distinct_all()

human_string <- human_uniprot_data %>% 
  filter(id_type %in% c("STRING"))  %>%
  select(uniprot, string = id) %>% 
  distinct_all()

human_full <- full_join(human_uniprotkb, human_gene_name) %>% 
  full_join(., human_ensembl) %>% 
  full_join(., human_string)

human_full %>% 
  write_tsv(str_c(data_dir, "UniProt_mapping_files/", "uniprot_mapping_human.tsv"))


# mouse data
mouse_uniprot_data <-  read_tsv(str_c(data_dir, "UniProt_mapping_files/", "MOUSE_10090_idmapping.dat"), col_names = c("uniprot", "id_type", "id"))

mouse_uniprotkb <- mouse_uniprot_data %>% 
  filter(id_type %in% c("UniProtKB-ID")) %>% 
  select(uniprot, uniprotkb_id = id)

mouse_gene_name <- mouse_uniprot_data %>% 
  filter(id_type %in% c("Gene_Name")) %>% 
  select(uniprot, gene_name = id)

mouse_ensembl <-  mouse_uniprot_data %>% 
  filter(id_type %in% c("Ensembl"))  %>%
  select(uniprot, ensembl = id) %>% 
  distinct_all()

mouse_full <- full_join(mouse_uniprotkb, mouse_gene_name) %>% 
  full_join(., mouse_ensembl)

mouse_full %>% 
  write_tsv(str_c(data_dir, "UniProt_mapping_files/", "uniprot_mapping_mouse.tsv"))

# ortholog information
human_prot_names <- read_tsv(str_c(data_dir, "Phylomedb_data/", "phylomedb_human_proteinname.txt")) %>% clean_names() %>% rename(human_uniprot = protein_name, number_seqid = number_number_phylome_db_id)
mouse_prot_names <- read_tsv(str_c(data_dir, "Phylomedb_data/", "phylomedb_mouse_proteinname.txt")) %>% clean_names() %>% rename(mouse_uniprot = protein_name, orthologid = number_number_phylome_db_id)
human_ortho <- read_tsv(str_c(data_dir, "Phylomedb_data/", "phylomedb_human_ortholog_info.txt")) %>% clean_names()

mapped_to_ortho <- human_ortho %>% 
  inner_join(., human_prot_names) %>% 
  left_join(., mouse_prot_names) %>% 
  select(ends_with("uniprot")) %>% 
  distinct_all() %>% 
  write_tsv(str_c(data_dir, "Phylomedb_data/", "human_mouse_uniprot_orthologs.tsv"))

# uniprot history
uniprot_history <- read_delim(str_c(data_dir, "UniProt_mapping_files/", "sec_ac.txt"), skip = 31, col_names = c("uniprot_all", "uniprot_primary"), delim = " ") %>% 
  mutate(across(everything(), str_trim)) %>% 
  nest(data = everything()) %>% 
  mutate(data = map(data, ~.x %>% 
                      bind_rows(., .x %>% select(uniprot_primary) %>% mutate(uniprot_all = uniprot_primary)) %>% 
                      arrange(uniprot_primary)
  )
  ) %>% 
  unnest(data) %>% 
  distinct_all()

uniprot_history %>% 
  write_tsv(str_c(data_dir, "UniProt_mapping_files/", "uniprot_id_history.tsv"))

# Ensembl mapping 
human_ensembl_data <- read_tsv(str_c(data_dir, "ensembl/", "human_stable_id_event.txt"), 
                               col_names = c("old_stable_id", "old_version", "new_stable_id", "new_version", "mapping_session_id", "type", "score"), 
                               col_types = cols(.default = "c")) %>% 
  mutate(across(everything(), ~recode(., "\\N" = NA_character_))) %>% 
  transmute(new_stable_id = new_stable_id, 
            old_stable_id = ifelse(is.na(old_stable_id), new_stable_id, old_stable_id)) %>% 
  filter(!is.na(new_stable_id)) %>% 
  distinct_all()
            
mouse_ensembl_data <- read_tsv(str_c(data_dir, "ensembl/", "mouse_stable_id_event.txt"), 
                               col_names = c("old_stable_id", "old_version", "new_stable_id", "new_version", "mapping_session_id", "type", "score"), 
                               col_types = cols(.default = "c")) %>% 
  mutate(across(everything(), ~recode(., "\\N" = NA_character_))) %>% 
  transmute(new_stable_id = new_stable_id, 
            old_stable_id = ifelse(is.na(old_stable_id), new_stable_id, old_stable_id)
  ) %>% 
  filter(!is.na(new_stable_id)) %>%
  distinct_all()

bind_rows(human_ensembl_data, mouse_ensembl_data) %>% 
  distinct_all() %>% 
  write_tsv(str_c(data_dir, "ensembl/", "human_mouse_ensembl_mapping.tsv"))
