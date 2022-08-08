# targets for report example
library(tidyverse)

data_dir <- here::here("data/")

abbasi_2016 <- read_tsv(str_c(data_dir, "prepped_data/", "known_biomarkers_list_prepped.tsv")) #%>% filter(source == "Abbasi 2016")
huang_deng_2019 <- read_tsv(str_c(data_dir, "example_data/", "breastcancer_targets_HuangDeng2019.txt"))
zakynthinos_pappa_2009 <- read_tsv(str_c(data_dir, "example_data/", "cad_targets_ZakynthinosPappa2009.txt"))
uniprot_mapping <- read_tsv(str_c(data_dir, "prepped_data/", "uniprot_mapping_file_mousehuman.tsv"))

n_targets <- 3

seed <- 999 #111

# diabetes targets 
set.seed(seed)

dm_targets <- abbasi_2016 %>% 
  filter(source == "Abbasi 2016") %>% 
  distinct(uniprot) %>% 
  sample_n(size = n_targets, replace = FALSE) %>% 
  mutate(disease = "diabetes")
  
# breast cancer targets
set.seed(seed)
bc_targets <- huang_deng_2019 %>% 
  distinct(uniprot) %>% 
  filter(!uniprot %in% abbasi_2016$uniprot) %>% 
  sample_n(size = n_targets, replace = FALSE) %>% 
  mutate(disease = "breast cancer")

# coronary artery disease targets
set.seed(seed)
cad_targets <- zakynthinos_pappa_2009 %>% 
  filter(!is.na(uniprot)) %>% 
  distinct(uniprot) %>% 
  filter(!uniprot %in% abbasi_2016$uniprot) %>% 
  sample_n(size = n_targets, replace = FALSE) %>% 
  mutate(disease = "CAD")

# all targets
all_targets <- bind_rows(dm_targets, bc_targets) %>% 
  bind_rows(., cad_targets) %>% 
  left_join(., uniprot_mapping %>% distinct(uniprot, uniprotkb_id, gene_name))

## save as dataframe
all_targets %>% 
  write_tsv(str_c(data_dir, "example_data/", "manuscript_targets_info.tsv"))

## save as .csv vector
all_targets %>% 
  pull(uniprot) %>%
  str_c(collapse = ",") %>% 
  write_lines(str_c(data_dir, "example_data/", "manuscript_targets_uniprot.csv"))



# Save Abbasi targets only
abbasi_2016 %>% 
  filter(source == "Abbasi 2016") %>% 
  distinct(uniprot) %>% 
  sample_n(size = 10, replace = FALSE) %>% 
  pull(uniprot) %>%
  str_c(collapse = ",") %>% 
  write_lines(str_c(data_dir, "example_data/", "manuscript_abbasi_targets_uniprot.csv"))


