# Download and prepare blood expression and excretion data 
library(tidyverse)
library(janitor)
# library(R.utils)

data_dir <- here::here("data/HPA_blood_proteins/")
uniprot_map <- read_tsv(here::here("data/prepped_data/uniprot_mapping_file_mousehuman.tsv")) %>%
  distinct(uniprot, ensembl)

# download
hpa_link <- "https://www.proteinatlas.org/download/proteinatlas.tsv.zip" 
hpa_dest <- str_c(data_dir, "proteinatlas.tsv.zip")

download.file(hpa_link, destfile = hpa_dest)
unzip(hpa_dest, exdir = data_dir)


# Prepare data
hpa_file <- read_tsv(str_c(data_dir, "proteinatlas.tsv"))

hpa_blood <- hpa_file %>% 
  clean_names() %>% 
  select(uniprot, gene, ensembl, 
         secretome_location, single_cell_type_rna_monocytes_n_tpm, single_cell_type_rna_plasma_cells_n_tpm,
         starts_with("blood_concentration"), starts_with("rna_blood_lineage"), starts_with("rna_blood_cell"), starts_with("blood_rna")
         ) %>% 
  left_join(., uniprot_map, by = "ensembl") %>% 
  mutate(uniprot = ifelse(is.na(uniprot.x), uniprot.y, uniprot.x)) %>% 
  select(-c(uniprot.x, uniprot.y)) %>% 
  distinct_all() %>% 
  filter(!is.na(uniprot)) %>% 
  nest(data = everything()) %>% 
  mutate(found_in_blood = map(data, ~.x %>% filter(secretome_location == "Secreted to blood" | 
                                             across(starts_with("single_cell_type_rna_"), ~. != 0) |
                                               across(starts_with("blood_concentration"), ~!is.na(.)) | 
                                               across(starts_with("rna_blood_"), ~!is.na(.) & . != "Not detected") | 
                                               across(starts_with("blood_rna_"), ~. != 0)) %>% 
                                mutate(found_in_blood = TRUE) %>% 
                                select(uniprot, found_in_blood)
                              ),
         not_found_in_blood = map2(data, found_in_blood, ~ anti_join(.x, .y) %>% 
                                     mutate(found_in_blood = FALSE) %>% 
                                     select(uniprot, found_in_blood)
                                   )
         ) %>% 
  transmute(hpa_blood = map2(found_in_blood, not_found_in_blood, full_join)) %>%
  unnest(hpa_blood)

# Uniprot IDs are up to date
# hpa_blood %>%
#   distinct(uniprot) %>% 
#   inner_join(., uniprot_map) %>% 
#   distinct(uniprot)

hpa_blood %>% 
  write_tsv(str_c(data_dir, "hpa_proteins_in_blood.tsv"))

            