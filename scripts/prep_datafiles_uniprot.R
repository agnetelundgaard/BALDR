# map all files to reference UniProt mapping

# libraries
library(tidyverse)
library(readxl)
library(janitor)

data_dir <- here::here("data/")
# data_dir <- "OneDrive - University of Copenhagen/Documents/RHAPSODY/customized_report/data/"

# Mapping files from UniProt ----

# ~ Primary and secondary UniProt IDs ----
# contains UniProt IDs where there have been made changes (exists secondary accession numbers)
uniprot_history <- read_tsv(str_c(data_dir, "UniProt_mapping_files/", "uniprot_id_history.tsv"))
ensembl_history <- read_tsv(str_c(data_dir, "ensembl/", "human_mouse_ensembl_mapping.tsv")) 


# ~ Human and mouse mapping IDs ----
# Human UniProt + ENSG mapping
human_uniprot_data <- read_tsv(str_c(data_dir, "UniProt_mapping_files/", "uniprot_mapping_human.tsv")) %>% 
  left_join(., uniprot_history, by = c("uniprot" = "uniprot_primary")) %>% 
  mutate(uniprot_all = ifelse(is.na(uniprot_all), uniprot, uniprot_all)) 

human_ensembl_data <- ensembl_history %>% 
  right_join(., human_uniprot_data %>% distinct(uniprot, ensembl),
            by = c("old_stable_id" = "ensembl")) %>%
  left_join(., ensembl_history, by = "new_stable_id") %>% 
  select(uniprot, 
         ensembl = new_stable_id, 
         ensembl_all = old_stable_id.y) 

human_mapping <- human_uniprot_data %>% 
  left_join(., human_ensembl_data %>% select(uniprot, ensembl_all)) %>% 
  select(uniprot:gene_name, ensembl = ensembl_all, string, uniprot_all) %>% 
  distinct_all()


# Mouse UniProt + ENSMUS mapping w/ ortholog
mouse_uniprot <- read_tsv(str_c(data_dir, "UniProt_mapping_files/", "uniprot_mapping_mouse.tsv")) %>% 
  left_join(., uniprot_history, by = c("uniprot" = "uniprot_primary")) %>% 
  mutate(uniprot_all = ifelse(is.na(uniprot_all), uniprot, uniprot_all)) %>% 
  distinct_all() 

ortholog_file <- read_tsv(str_c(data_dir, "Phylomedb_data/", "human_mouse_uniprot_orthologs.tsv"))
orthologs <- ortholog_file %>% 
  left_join(., human_mapping, by = c("human_uniprot" = "uniprot_all")) %>% 
  distinct_all() %>% 
  select(-human_uniprot) %>% 
  rename_with(~str_c("human_", .), .cols = c(uniprot:ensembl)) %>% 
  left_join(., mouse_uniprot %>% select(uniprot, uniprot_all), by = c("mouse_uniprot" = "uniprot_all")) %>% 
  filter(!is.na(uniprot) & !is.na(human_uniprot)) %>% 
  select(mouse_uniprot = uniprot, everything(), -mouse_uniprot)

mouse_ensembl_data <- ensembl_history %>% 
  right_join(., mouse_uniprot %>% distinct(uniprot, ensembl),
             by = c("old_stable_id" = "ensembl")) %>% 
  left_join(., ensembl_history, by = "new_stable_id") %>%
  select(uniprot, ensembl = old_stable_id.y) %>% 
  left_join(., orthologs, by = c("uniprot" = "mouse_uniprot")) %>% 
  filter(!is.na(human_uniprot)) %>% 
  select(uniprot = human_uniprot, 
         ensembl)

mouse_mapping <- mouse_uniprot %>% 
  left_join(., orthologs, by = c("uniprot" = "mouse_uniprot")) %>% 
  distinct_all() %>% 
  mutate(uniprot = ifelse(!is.na(human_uniprot), human_uniprot, uniprot),
         uniprotkb_id = ifelse(!is.na(human_uniprotkb_id), human_uniprotkb_id, uniprotkb_id),
         gene_name = ifelse(!is.na(human_gene_name), human_gene_name, gene_name),
         ensembl = ifelse(!is.na(human_ensembl), human_ensembl, ensembl)) %>% 
  select(-starts_with("human")) %>% 
  left_join(., mouse_ensembl_data, by = "uniprot")%>% 
  select(uniprot:gene_name, ensembl = ensembl.y, uniprot_all) %>% 
  distinct_all() 

# ~ Shared mouse and human mapping ----

humanmouse_mapping <- bind_rows(human_mapping, mouse_mapping) 

humanmouse_mapping %>% 
  write_tsv(str_c(data_dir, "prepped_data/", "uniprot_mapping_file_mousehuman.tsv"))

# Mapping of data files ----

# ~ RHAPSODY experimental data  ----

# ~~ Reference data ----
# 10,946 commonIDs not mapping. ENSMUS and ENSG IDs 
ref_exp_file <- read_csv(str_c(data_dir, "/data_BCDP_rhapsody.csv"), col_types = cols(More.info.on.statistical.test = col_character(),
                                                                                      NumberAnalytes = col_double())) 
common_id_mapping <- ref_exp_file %>%
  clean_names() %>% 
  filter(omics %in% c("Transcriptomics", "Proteomics")) %>% 
  distinct(common_id) %>% 
  left_join(., humanmouse_mapping %>% 
              pivot_longer(cols = c(ensembl, uniprot_all), 
                           names_to = "source_type", 
                           values_to = "mapping") %>% 
              select(-source_type), 
            by = c("common_id" = "mapping")) %>% 
  filter(!is.na(uniprot)) %>% 
  distinct_all() 

ref_exp_file_mapped <- ref_exp_file %>%
  clean_names() %>% 
  filter(omics %in% c("Transcriptomics", "Proteomics")) %>% 
  left_join(., common_id_mapping) %>% 
  filter(!is.na(uniprot)) %>% 
  distinct_all()

ref_exp_file_mapped %>% 
  write_tsv(str_c(data_dir, "prepped_data/", "/data_BCDP_rhapsody_prepped.tsv"))

# ~~ Experimental properties
load(str_c(data_dir, "experimentPropertiesDFwide.RData"))

experimentPropertiesDFwide_prepped <- experimentPropertiesDFwide %>% 
  as_tibble() %>% 
  clean_names()

experimentPropertiesDFwide_prepped %>% 
  write_tsv(str_c(data_dir, "prepped_data/", "experimentPropertiesDFwide_prepped.tsv"))

# ~ PPI data ----
# ~~ Intomics ----
# 1 Uniprot NA removed. 10 UniProt IDs without gene names kept.
ppi_data_file <- read_tsv(str_c(data_dir, "core.psimitab"),
                          col_names = c("A", "B", "altA", "altB", "aliasA", "aliasB", "detectionMethod", "firstAuthor", "pubID", "taxonA", "taxonB", "interactionType", "sourceDatabases", "interactionID", "confidence_score", "provider"))

ppi_data_prep <- ppi_data_file %>% 
  nest(ppi_data = everything()) %>% 
  bind_cols(., human_mapping %>% nest(mapping_data = everything())) %>% 
  mutate(ppi_data_prep = map(ppi_data, ~.x %>% 
                               mutate(confidence_score = str_remove(confidence_score, "\\|[[:graph:]]{1,6}$") %>% as.numeric()) %>% 
                               mutate(across(c(A, B), ~str_remove(., "^uniprotkb:")))
                               ),
         ppi_data_mapped = map2(ppi_data_prep, mapping_data, ~.x %>% 
                                  left_join(., .y %>% select(uniprot_all, A_new = uniprot, A_name = gene_name), by = c("A" = "uniprot_all")) %>% 
                                  left_join(., .y %>% select(uniprot_all, B_new = uniprot, B_name = gene_name), by = c("B" = "uniprot_all")) %>% 
                                  select(A = A_new, B = B_new, A_name, B_name, confidence_score) %>% 
                                  filter(across(A:B_name, ~!is.na(.))) %>% 
                                  distinct_all()
                                )) %>% 
  transmute(ppi_AB_all = map(ppi_data_mapped, ~ bind_rows(.x %>% select(A, B, confidence_score),
                                                          .x %>% select(A = B, B = A, confidence_score)) %>% 
                               distinct_all()
                             ),
            ppi_names = map(ppi_data_mapped, ~ bind_rows(.x %>% select(uniprot = A, gene_name = A_name),
                                                         .x %>% select(uniprot = B, gene_name = B_name)) %>% 
                              distinct_all()
           
         )) 

ppi_data_prep %>% 
  saveRDS(str_c(data_dir, "prepped_data/", "ppi_data_prepped.Rds"))

# ~~ STRING DB ----
string_db <- read_delim(str_c(data_dir, "9606.protein.physical.links.detailed.v11.5.txt"), delim = " ") %>% 
  nest(ppi_data = everything()) %>% 
  bind_cols(.,  human_mapping %>% nest(mapping_data = everything())) %>% 
  mutate(ppi_data_mapped = map2(ppi_data, mapping_data, ~.x %>% 
                                  rename(A_string = protein1, B_string = protein2) %>% 
                                  left_join(., .y %>% select(string, A = uniprot, A_name = gene_name), by = c("A_string" = "string")) %>% 
                                  left_join(., .y %>% select(string, B = uniprot, B_name = gene_name), by = c("B_string" = "string")) %>% 
                                  mutate(combined_score = combined_score / 1000,
                                         experimental_score = experimental / 1000, 
                                         database_score = database / 1000, 
                                         textmining_score = textmining / 1000
                                         ) %>% 
                                  select(A, B, A_name, B_name, experimental_score, database_score, textmining_score, combined_score) %>% 
                                  filter(across(A:B_name, ~!is.na(.))) %>% 
                                  distinct_all()
                                  )
         ) %>% 
  transmute(ppi_AB_all = map(ppi_data_mapped, ~ bind_rows(.x %>% select(A, B, combined_score),
                                                          .x %>% select(A = B, B = A, combined_score)) %>% 
                                distinct_all()),
            ppi_names = map(ppi_data_mapped, ~ bind_rows(.x %>% select(uniprot = A, gene_name = A_name),
                                                         .x %>% select(uniprot = B, gene_name = B_name)) %>% 
                              distinct_all()),
            ppi_detailed_scores = map(ppi_data_mapped, ~ bind_rows(.x %>% select(A, B, experimental_score, database_score, textmining_score, combined_score),
                                                                 .x %>% select(A = B, B = A, experimental_score, database_score, textmining_score, combined_score)) %>% 
                                        distinct_all())
            
            )

string_db %>% 
  saveRDS(str_c(data_dir, "prepped_data/", "string_ppi_data_prepped.Rds"))



# # ~ PPI enriched gene set data ----
# ppi_enriched_gene_sets <- read_tsv(str_c(data_dir, "/geneenrich_all_databases_only_enriched.tsv"))
# 
# ppi_enriched_gene_sets_mapped <- ppi_enriched_gene_sets %>% 
#   left_join(., human_mapping, by = c("uniprot" = "uniprot_all")) 


# ~ Biomarker matrix and known biomarkers ----

# ~~ Text mining data from old matrix ----
load(str_c(data_dir, "all_data_prioritization_report.RData"))

text_mining_prepped <- biomarker_list %>% 
  clean_names() %>% 
  left_join(., human_mapping, by = c("entry" = "uniprot_all")) %>%
  select(uniprot,
         uniprotkb_id,
         uniprot_gene = gene_name,
         everything(),
         -c(entry, entry_name, starts_with("ensembl"))
         ) %>%
  filter(!is.na(uniprot) & !(is.na(gene_synonym) & is.na(uniprot_gene))) %>%
  distinct_all() %>% 
  select(uniprot, uniprotkb_id, uniprot_gene, pubmed_diabetes, pubmed_all_other, fulltext_diabetes, fulltext_all_other)

text_mining_prepped %>% 
  write_tsv(str_c(data_dir, "prepped_data/", "textmining_results_prepped.tsv"))

# ~~ Biomarker matrix ---- 
biomarker_list <- read_tsv(str_c(data_dir, "2022-07-11-baldr_matrix.tsv")) 

biomarker_list_prepped <- biomarker_list %>% 
  rename(uniprotkb_id = uniprot_entry_name, 
         uniprot = uniprot_id) %>% 
  mutate(.,
         opentargets_assoc_label = str_replace_all(opentargets_assoc_label, "Type 2", "T2DM"),
         opentargets_assoc_label = str_replace_all(opentargets_assoc_label, "Type 1", "T1DM"),
         opentargets_assoc_which_others = str_to_title(opentargets_assoc_which_others),
         # opentargets_assoc_which_others = str_replace_all(opentargets_assoc_which_others, "^-$", NA_character_)
         )
  

biomarker_list_prepped %>% 
    write_tsv(str_c(data_dir, "prepped_data/", "biomarker_list_prepped.tsv"))

# ~~ OLD VERSION OF BIOMARKER MATRIX
# load(str_c(data_dir, "all_data_prioritization_report.RData"))
# 51 UniProt IDs not mapping - seem to obsolete
# biomarker_list_prepped <- biomarker_list %>% 
#   clean_names() %>% 
#   left_join(., human_mapping, by = c("entry" = "uniprot_all")) %>%
#   select(uniprot,
#          uniprotkb_id, 
#          gene_name, 
#          everything(), 
#          -c(entry, entry_name, starts_with("ensembl"))
#          ) %>% 
#   filter(!is.na(uniprot) & !(is.na(gene_synonym) & is.na(gene_name))) %>% 
#   distinct_all()
# 
# biomarker_list_prepped %>% 
#   write_tsv(str_c(data_dir,"biomarker_list_prepped_Danai.tsv"))




# ~~ Known biomarkers ----
# The list is using updated versions of UniProt IDs
known_biomarkers_list <- read_xlsx(str_c(data_dir, "BMKTF_Rhapsody_full_table_for_Karla_080317.xlsx"))

known_biomarkers_list_prepped <- known_biomarkers_list %>% 
  clean_names() %>% 
  rename(ensembl = ensembl_gene_id) %>%
  left_join(., human_mapping %>% select(-c(uniprot_all, ensembl)) %>% distinct_all()) %>% 
  mutate(source = ifelse(uniprot %in% c("Q15848", "P02647", "P04114", "P01185", "P01034", "P02765", "P01160", "P05164", "Q9BZZ2", "P04004"),
                         "Abbasi 2016",
                         source
                         ),
         source = recode(source, 
                         "Rhapsody partner suggestions" = "RHAPSODY partner suggestion",
                         "Abbasi 2016" = "Abbasi et al, 2016"
                         ),
         relevance_score = as.numeric(relevance_score)
         )

known_biomarkers_list_prepped %>% 
  write_tsv(str_c(data_dir, "prepped_data/", "known_biomarkers_list_prepped.tsv"))


# ~~ DrugBank ----
drugbank_data <- read_tsv(str_c(data_dir, "2022-08-01-drugbank_protein_interactions.tsv"))

drugbank_prepped <- drugbank_data %>% 
  select(uniprot = target_uniprot, drug_name, drug_drugbank_id, drug_groups, target_type, pharmacological_action, drug_url, interaction_url)

drugbank_prepped %>% 
  write_tsv(str_c(data_dir, "prepped_data/", "drugbank_data_prepped.tsv"))
