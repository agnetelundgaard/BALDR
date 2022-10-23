### Packages and Import -----
library(tidyverse)
library(xml2)

drugbank_raw <- read_xml("data/2022-07-18-namespace_free_drugbank.xml")

### Drug Metadata ------
drug_nodeset <- xml_find_all(drugbank_raw, "drug")

meta_drugs <- tibble(node_index = seq_along(drug_nodeset),
                     drug_id = map_chr(drug_nodeset, ~xml_text(xml_find_first(., xpath = "./drugbank-id"))),
                     drug_name = map_chr(drug_nodeset, ~xml_text(xml_find_first(., xpath = "./name"))),
                     drug_groups = map_chr(drug_nodeset, ~paste(xml_text(xml_find_all(., xpath = "./groups/group")), collapse = ","))
)

### Interactions -------
carrier_nodeset <- map(drug_nodeset, xml_find_all, xpath = "./carriers/carrier")
enzyme_nodeset <- map(drug_nodeset, xml_find_all, xpath = "./enzymes/enzyme")
target_nodeset <- map(drug_nodeset, xml_find_all, xpath = "./targets/target")
transporter_nodeset <- map(drug_nodeset, xml_find_all, xpath = "./transporters/transporter")

# in_nodeset <- target_nodeset
# category <- "standard"

extract_interaction_info <- function(in_nodeset, category) {
  first_df <- tibble(node_index = seq_along(in_nodeset),
                     target_type = category,
                     nodeset = in_nodeset)
  
  second_df <- first_df %>% 
    mutate(i = map(nodeset, ~ seq_along(.)),
           target_id = map(nodeset, ~xml_text(xml_find_first(., xpath = "./id"))),
           target_name = map(nodeset, ~xml_text(xml_find_first(., xpath = "./name"))),
           target_organism = map(nodeset, ~xml_text(xml_find_first(., xpath = "./organism"))),
           known_action = map(nodeset, ~xml_text(xml_find_first(., xpath = "./known-action"))),
           target_uniprot = map(nodeset, ~xml_text(xml_find_first(., xpath = "./polypeptide/external-identifiers/external-identifier[resource='UniProtKB']/identifier")))
    ) %>% 
    select(-nodeset) %>%
    unnest(cols = c(i, target_id, target_name, target_organism, known_action, target_uniprot))
  
  return(second_df)
}

extracted_carrier_info <- extract_interaction_info(carrier_nodeset, "carrier")
extracted_enzyme_info <- extract_interaction_info(enzyme_nodeset, "enzyme")
extracted_target_info <- extract_interaction_info(target_nodeset, "target")
extracted_transporter_info <- extract_interaction_info(transporter_nodeset, "transporter")

all_interactions <- extracted_carrier_info %>%
  bind_rows(extracted_enzyme_info) %>%
  bind_rows(extracted_target_info) %>%
  bind_rows(extracted_transporter_info)

# ATC Codes
atc_meta <- tibble(
  node_index = seq_along(drug_nodeset),
  drug_atc = map_chr(drug_nodeset, ~paste(xml_attr(xml_find_all(., xpath = "./atc-codes/atc-code"), attr = "code"), collapse = ","))
)

initial_out <- all_interactions %>% 
  filter(target_organism == "Humans") %>% 
  filter(! is.na(target_uniprot)) %>% 
  left_join(meta_drugs) %>% 
  left_join(atc_meta) %>% 
  mutate(drug_atc = if_else(drug_atc == "", as.character(NA), drug_atc),
         drug_url = paste0("https://go.drugbank.com/drugs/", drug_id),
         interaction_url = paste0(drug_url, "#", target_id), 
         drugbank_version = "5.1.9",
         drugbank_version_release = "2022-01-04",
         drugbank_version_download = "2022-07-18") %>% 
  select(target_uniprot, target_name, drug_name, drug_groups, drug_atc, target_type, pharmacological_action = known_action, drug_drugbank_id = drug_id, target_drugbank_id = target_id, drug_url, interaction_url, drugbank_version, drugbank_version_release, drugbank_version_download, node_index, i)

write_tsv(initial_out, file = "2022-08-02-drugbank_protein_interactions.tsv.gz")
