library(tidyverse)
library(DBI)

# Load uniprot alias history from Pharos data -----------------------
con <- dbConnect(RMariaDB::MariaDB(), group = "tcrd") # connection defined in the file ~/.my.cnf (on Mac)

pharos_protein_raw <- collect(tbl(con, "protein"))
pharos_alias_raw <- collect(tbl(con, "alias")) 

pharos_uniprot_alias <- pharos_alias_raw %>%
  filter(type == "uniprot") %>% 
  select(protein_id, value) %>% 
  mutate(version = "historic") %>% 
  bind_rows(select(pharos_protein_raw, protein_id = id, value = uniprot)) %>% 
  mutate(version = if_else(is.na(version), "current", "historic")) %>% 
  arrange(protein_id)

# Uniprot ----------------------
uniprot_raw <- read_delim("data/2022_03_uniprot_human_reviewed.tsv", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)

uniprot_df <- uniprot_raw %>% 
  select(Entry, `Entry Name`, `Protein names`, `Gene Names`) %>% 
  mutate(gene = str_replace(`Gene Names`, "\\s+.+", ""))

# HPA --------------------------
hpa_raw <- read_delim("data/v21.1-proteinatlas.tsv", 
                      "\t", escape_double = FALSE, trim_ws = TRUE)

hpa_temp <- hpa_raw %>% 
  select(uniprot_id = Uniprot, hpa_gene = Gene, hpa_gene_synonyms = `Gene synonym`,
         hpa_protein_class = `Protein class`, hpa_subcellular_location = `Subcellular location`) %>% 
  filter(! is.na(uniprot_id))

hpa_duplicate_ids <- hpa_temp %>% 
  group_by(uniprot_id) %>% 
  filter(n() > 1) %>% 
  pull(uniprot_id) %>% 
  unique()

custom_str_split <- function(str_in, split_on = ", ", return_i = "all") {
  str1 <- str_split(string = str_in, pattern = split_on)
  str2 <- str1[[1]][str1[[1]] != "NA"]
  str3 <- unique(str2)
  if (return_i == "all") return(paste(str3, collapse = split_on))
  if (return_i == "first") return(str3[1])
  if (return_i == "rest") return(paste(str3[-1], collapse = split_on))
}

hpa_fixed_duplicates <- hpa_temp %>% 
  filter(uniprot_id %in% hpa_duplicate_ids) %>% 
  group_by(uniprot_id) %>% 
  summarise(all_genes = paste(c(hpa_gene, hpa_gene_synonyms), collapse = ", "),
            all_protein_class = paste(hpa_protein_class, collapse = ", "),
            all_locations = paste(hpa_subcellular_location, collapse = ",")) %>% 
  ungroup() %>% 
  rowwise() %>% 
  mutate(hpa_gene = custom_str_split(all_genes, return_i = "first"),
         hpa_gene_synonyms = custom_str_split(all_genes, return_i = "rest"),
         hpa_protein_class = custom_str_split(all_protein_class),
         hpa_subcellular_location = na_if(custom_str_split(all_locations, split_on = ","), "")) %>% 
  select(uniprot_id, hpa_gene, hpa_gene_synonyms,
         hpa_protein_class, hpa_subcellular_location)

hpa_df <- hpa_temp %>% 
  filter(! uniprot_id %in% hpa_duplicate_ids) %>% 
  bind_rows(hpa_fixed_duplicates)
  
# Opentargets ------------------------
optarg_con <- file("./data/2022-08-04-opentargets_alltargets/combined_targets.json", "r")
optarg_all_targets <- jsonlite::stream_in(optarg_con)

optarg_proteins <- optarg_all_targets %>% 
  filter(biotype == "protein_coding") %>% 
  select(optarg_id = id, approvedSymbol, approvedName, proteinIds) %>% 
  unnest(cols = c(proteinIds)) #%>% 
# filter(source != "ensembl_PRO", source != "uniprot_trembl") %>%
# group_by(optarg_id) %>% 
# mutate(running_number = row_number()) %>% 
# ungroup() %>% 
# filter(running_number == 1)

optarg_current_uniprot <- optarg_proteins %>% 
  filter(source == "uniprot_swissprot") %>% 
  select(optarg_symbol = approvedSymbol) %>% 
  distinct(optarg_symbol)

read_tsv_with_efo <- function(path) {
  df <- read_tsv(path, 
                 trim_ws = TRUE,
                 col_types = cols(.default = col_double(),
                                  symbol = col_character(),
                                  targetName = col_character())) 
  
  # downloads of release 22.06 have missing column names, so we replace them here
  colnames(df) <- c("symbol", "overallAssociationScore", "geneticAssociations", 
                    "somaticMutations", "drugs", "pathwaysSystemsBiology", "textMining", 
                    "rnaExpression", "animalModels", "targetName")
  
  df$efo_id <- str_extract(path, "(?<=tsv/).+(?=-associated)")
  
  df <- df %>% 
    group_by(symbol) %>% 
    mutate(highest_overall = min_rank(desc(overallAssociationScore))) %>% 
    ungroup() %>% 
    filter(highest_overall == 1) %>% 
    select(- highest_overall)
  
  return(distinct(df))
}

# another renaming for easier compatibility with existing downstream scripts
optarg_out_column_names <- c("overall", "genetic_association", "somatic_mutation", "known_drug", "affected_pathway", "literature", "rna_expression", "animal_model")

# Manual downloads from the landing page when looking up the EFO disease codes for Diabetes mellitus and its EFO ontology descendants (of which there are 12 as of 2022-03-17)
optarg_descendant_info <- tibble::tribble(
  ~efo_id, ~disease, ~diabetes_subtype,
  "Orphanet_183625", "Rare genetic diabetes mellitus", "other",
  "EFO_1001511", "monogenic diabetes", "other",
  "EFO_1001121", "prediabetes syndrome", "other",
  "EFO_0010164", "insulin-resistant diabetes mellitus", "other",
  "EFO_0007498", "Stiff-Person syndrome", "other",
  "EFO_0004593", "gestational diabetes", "gdm",
  "MONDO_0019193", "acquired generalized lipodystrophy", "other",
  "MONDO_0018575", "microcephalic primordial dwarfism-insulin resistance syndrome", "other",
  "MONDO_0014497", "polyendocrine-polyneuropathy syndrome", "other",
  "MONDO_0005148", "type 2 diabetes mellitus", "type_2",
  "MONDO_0005147", "type 1 diabetes mellitus", "type_1",
  "HP_0000857", "Neonatal insulin-dependent diabetes mellitus", "other") %>% 
  mutate(url = paste0("https://platform.opentargets.org/disease/", efo_id),
         assoc_file = paste0("./data/2022-08-04-opentargets_tsv/", efo_id, "-associated-diseases.tsv"))

# This will throw warnings for some files, but all the problematic fields just turn into NAs (which is what we want)
optarg_all_subtypes <- map(optarg_descendant_info$assoc_file, 
                           read_tsv_with_efo) %>% 
  bind_rows() %>% 
  left_join(select(optarg_descendant_info, efo_id, disease, diabetes_subtype)) %>% 
  filter(overallAssociationScore > 0) 

optarg_other_assoc_details <- optarg_all_subtypes %>% 
  filter(diabetes_subtype == "other") %>% 
  group_by(symbol) %>% 
  summarise(assoc_which_others = paste(disease, collapse = ", "))

optarg_assoc_yesno <- optarg_all_subtypes %>% 
  distinct(symbol, diabetes_subtype) %>% 
  mutate(present = T) %>% 
  pivot_wider(names_from = diabetes_subtype, values_from = present, values_fill = F) %>% 
  mutate(assoc_concat = paste0(if_else(type_1, "Type 1, ", ""), 
                               if_else(type_2, "Type 2, ", ""),
                               if_else(gdm, "GDM, ", ""),
                               if_else(other, "Other", "")),
         assoc_label = str_replace(assoc_concat, ", $", "")) %>% # remove comma at the end of the label
  left_join(optarg_other_assoc_details)

optarg_supergroup <- read_tsv_with_efo("./data/2022-08-04-opentargets_tsv/EFO_0000400-associated-diseases.tsv") %>% 
  filter(overallAssociationScore > 0) %>% 
  replace_na(list(overallAssociationScore = 0,
                  geneticAssociations = 0,
                  somaticMutations = 0,
                  drugs = 0,
                  pathwaysSystemsBiology = 0,
                  textMining = 0,
                  rnaExpression = 0,
                  animalModels = 0)) %>% 
  group_by(symbol) %>% 
  mutate(highest_overall = min_rank(desc(overallAssociationScore))) %>% 
  ungroup() %>% 
  filter(highest_overall == 1) %>% 
  select(- highest_overall)

# Supergroup scores
optarg_scores_supergroup <- optarg_supergroup %>% 
  select(symbol, where(is.numeric))

colnames(optarg_scores_supergroup) <- c("optarg_symbol", 
                                        paste0("opentargets_supergroup_association_", optarg_out_column_names))

# T1d
optarg_scores_t1d <- optarg_all_subtypes %>% 
  filter(diabetes_subtype == "type_1") %>% 
  select(symbol, where(is.numeric)) %>% 
  replace_na(list(overallAssociationScore = 0,
                  geneticAssociations = 0,
                  somaticMutations = 0,
                  drugs = 0,
                  pathwaysSystemsBiology = 0,
                  textMining = 0,
                  rnaExpression = 0,
                  animalModels = 0))

colnames(optarg_scores_t1d) <- c("optarg_symbol", 
                                 paste0("opentargets_t1d_association_", optarg_out_column_names))

# T2d
optarg_scores_t2d <- optarg_all_subtypes %>% 
  filter(diabetes_subtype == "type_2") %>% 
  select(symbol, where(is.numeric)) %>% 
  replace_na(list(overallAssociationScore = 0,
                  geneticAssociations = 0,
                  somaticMutations = 0,
                  drugs = 0,
                  pathwaysSystemsBiology = 0,
                  textMining = 0,
                  rnaExpression = 0,
                  animalModels = 0))

colnames(optarg_scores_t2d) <- c("optarg_symbol", 
                                 paste0("opentargets_t2d_association_", optarg_out_column_names))

# gdm
optarg_scores_gdm <- optarg_all_subtypes %>% 
  filter(diabetes_subtype == "gdm") %>% 
  select(symbol, where(is.numeric)) %>% 
  replace_na(list(overallAssociationScore = 0,
                  geneticAssociations = 0,
                  somaticMutations = 0,
                  drugs = 0,
                  pathwaysSystemsBiology = 0,
                  textMining = 0,
                  rnaExpression = 0,
                  animalModels = 0))

colnames(optarg_scores_gdm) <- c("optarg_symbol", 
                                 paste0("opentargets_gdm_association_", optarg_out_column_names))

# other
optarg_scores_other <- optarg_all_subtypes %>% 
  filter(diabetes_subtype == "other") %>% 
  group_by(symbol) %>% 
  summarise(overallAssociationScore = mean(overallAssociationScore, na.rm = T),
            geneticAssociations = mean(geneticAssociations, na.rm = T),
            somaticMutations = mean(somaticMutations, na.rm = T),
            drugs = mean(drugs, na.rm = T),
            pathwaysSystemsBiology = mean(pathwaysSystemsBiology, na.rm = T),
            textMining = mean(textMining, na.rm = T),
            rnaExpression = mean(rnaExpression, na.rm = T),
            animalModels = mean(animalModels, na.rm = T)) %>% 
  ungroup() %>% 
  replace_na(list(overallAssociationScore = 0,
                  geneticAssociations = 0,
                  somaticMutations = 0,
                  drugs = 0,
                  pathwaysSystemsBiology = 0,
                  textMining = 0,
                  rnaExpression = 0,
                  animalModels = 0))

colnames(optarg_scores_other) <- c("optarg_symbol", 
                                   paste0("opentargets_other_association_", optarg_out_column_names))

optarg_scores_final <- optarg_current_uniprot %>% 
  left_join(optarg_scores_supergroup) %>%
  left_join(optarg_scores_t1d) %>% 
  left_join(optarg_scores_t2d) %>%
  left_join(optarg_scores_gdm) %>%
  left_join(optarg_scores_other) 

optarg_scores_final[is.na(optarg_scores_final)] <- 0

optarg_master <- optarg_scores_final %>% 
  left_join(select(optarg_assoc_yesno, optarg_symbol = symbol, 
                   opentargets_assoc_label = assoc_label, 
                   opentargets_assoc_which_others = assoc_which_others)) %>% 
  replace_na(list(opentargets_assoc_label = "-",
                  opentargets_assoc_which_others = "-")) %>% 
  select(optarg_symbol, opentargets_assoc_label, opentargets_assoc_which_others, everything())

# Pharos data ------------------
pharos_target_raw <- collect(tbl(con, "target"))
pharos_idg_df <- pharos_target_raw %>% 
  mutate(protein_family = case_when(is.na(fam) ~ "Non-IDG",
                                    fam == "IC" ~ "Ion Channel",
                                    fam == "NR" ~ "Nuclear Receptor",
                                    fam == "TF" ~ "Transcription Factor",
                                    fam == "TF; Epigenetic" ~ "TF/Epigenetic",
                                    T ~ fam)) %>% 
  select(protein_id = id, pharos_tdl = tdl, pharos_protein_family = protein_family)

pharos_novelty <- collect(tbl(con, "tinx_novelty")) %>% 
  mutate(log10novelty = log10(score)) %>% 
  select(protein_id, pharos_novelty = score, pharos_log10novelty = log10novelty)

info_type_tbl <- collect(tbl(con, "info_type")) # describes features in tdl_info
pharos_tdl_info <- collect(tbl(con, "tdl_info"))

pharos_meta_integer <- pharos_tdl_info %>% 
  filter(itype %in% c("Ab Count", "EBI Total Patent Count", "MAb Count")) %>% 
  mutate(clean_category = str_to_lower(str_replace_all(itype, " ", "_"))) %>% 
  pivot_wider(id_cols = protein_id, names_from = clean_category, values_from = integer_value, names_prefix = "pharos_") %>% 
  mutate(pharos_derived_pab_count = pharos_ab_count - pharos_mab_count)

pharos_meta_numeric <- pharos_tdl_info %>% 
  filter(itype %in% c("JensenLab PubMed Score", "PubTator Score")) %>% 
  mutate(clean_category = str_to_lower(str_replace_all(itype, " ", "_"))) %>% 
  pivot_wider(id_cols = protein_id, names_from = clean_category, values_from = number_value, names_prefix = "pharos_")

pharos_meta_character <- pharos_tdl_info %>% 
  filter(itype == "UniProt Function") %>% 
  mutate(clean_category = str_to_lower(str_replace_all(itype, " ", "_"))) %>%
  # some proteins have more than one entry, we just grab the first one mentioned (arbitrary decision)
  group_by(protein_id) %>% 
  mutate(field_order = row_number()) %>% 
  ungroup() %>% 
  filter(field_order == 1) %>% 
  pivot_wider(id_cols = protein_id, names_from = clean_category, values_from = string_value) %>% 
  rename(pharos_protein_description_long = uniprot_function)

# Merge data sources
pharos_all_features <- pharos_protein_raw %>% 
  select(protein_id = id, pharos_protein_name = description) %>% 
  left_join(pharos_idg_df) %>% 
  left_join(pharos_meta_integer) %>% 
  left_join(pharos_meta_numeric) %>% 
  left_join(pharos_meta_character) %>% 
  left_join(pharos_novelty)

# Map IDs to Uniprot alias master table -------
lookup_pharos_x_uniprot <- uniprot_df %>% 
  select(Entry, gene) %>% 
  left_join(select(pharos_protein_raw, pharos_protein_id = id, uniprot), by = c("Entry" = "uniprot"))

lookup_pharos_x_symbol <- lookup_pharos_x_uniprot %>% 
  filter(is.na(pharos_protein_id)) %>%
  filter(Entry != "C0HLV8") %>% # Uniprot ID with no corresponding gene (as of 2022-05-10 https://www.uniprot.org/uniprot/C0HLV8), creates problems during left_join on gene names 
  select(-pharos_protein_id) %>%
  left_join(select(pharos_protein_raw, pharos_protein_id = id, sym), by = c("gene" = "sym")) 

final_lookup_pharos <- lookup_pharos_x_uniprot %>% 
  filter(!is.na(pharos_protein_id)) %>%
  bind_rows(lookup_pharos_x_symbol) %>% 
  select(uniprot_id = Entry, pharos_protein_id) %>% 
  bind_rows(tibble(uniprot_id = "C0HLV8", pharos_protein_id = as.integer(NA)))

# Collect all measures ----
master_matrix_wide <- uniprot_df %>% 
  select(uniprot_id = Entry, uniprot_entry_name = `Entry Name`, uniprot_protein_name = `Protein names`, 
         uniprot_gene = gene, uniprot_all_gene_names = `Gene Names`) %>% 
  # HPA
  left_join(hpa_df) %>% 
  # pharos
  left_join(final_lookup_pharos) %>% 
  left_join(pharos_all_features, by = c("pharos_protein_id" = "protein_id")) %>% 
  select(- pharos_protein_id) %>% 
  # opentargets
  left_join(optarg_master, by = c("uniprot_gene" = "optarg_symbol")) %>% 
  select(uniprot_id,
         uniprot_entry_name,
         uniprot_protein_name,
         pharos_protein_name, everything()) %>% 
  mutate(uniprot_version = "2022_03",
         uniprot_version_release = "2022-08-03",
         uniprot_version_download = "2022-08-04",
         hpa_version = "21.1",
         hpa_version_release = "2022-05-31",
         hpa_version_download = "2022-08-04",
         pharos_version = "TCRD v6.12.4",
         pharos_version_release = "2021-10-29",
         pharos_version_download = "2022-01-28",
         opentargets_version = "22.06",
         opentargets_version_release = "2022-06-24",
         opentargets_version_download = "2022-08-04")

write_tsv(master_matrix_wide, 
          file = "2022-08-05-corrected_baldr_matrix.tsv.gz")
