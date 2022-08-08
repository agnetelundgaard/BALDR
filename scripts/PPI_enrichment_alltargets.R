# PPI networks and gene set enrichment ----

# ~ Libraries ----
library(tidyverse)
library(WebGestaltR)
library(furrr)
library(tictoc)

# ~ Functions ---- 
  
  # Run WGR
WGR_run_ppi <- function(data, database){
    ppi_enrichResult <- NULL
    
    ppi_enrichResult <- 
      WebGestaltR(
        enrichMethod="ORA", 
        organism="hsapiens",
        enrichDatabase = database,
        interestGene = data$uniprot, interestGeneType="uniprotswissprot",
        referenceGene = ref_exp_uniprot$uniprot, referenceGeneType="uniprotswissprot", 
        fdrThr = 1, minNum = 10, fdrMethod = "BH",
        isOutput=FALSE, 
        projectName=NULL)
    
    return(ppi_enrichResult)
  }
  
  WGR_run_ppi_safe <- safely(WGR_run_ppi)
  
  # Prep WGR
  WGR_prep_ppi <- function(WGR_result, input_uniprot){
    if(nrow(WGR_result) > 0){
      prepped <- WGR_result %>%
        filter(str_detect(userId, input_uniprot)) %>%
        mutate(p_threshold = ifelse(FDR < 0.05, TRUE, FALSE)) #%>% 
      # ungroup() %>% 
      # mutate(description = fct_reorder(description, desc(enrichmentRatio)))
      
      nrow_prepped <- nrow(prepped)
      
      if(nrow_prepped > 0){
        ppi_enrichResult_prep <- prepped
      } else {
        ppi_enrichResult_prep <- NULL
      }
      
    } else {
      ppi_enrichResult_prep <- NULL
    }
    return(ppi_enrichResult_prep)
  }  
  


# ~ Data ---- 

data_dir <- here::here("data", "prepped_data/")
# data_dir <- "/Users/xgh849/OneDrive - University of Copenhagen/Documents/RHAPSODY/customized_report/data/prepped_data/"

output_dir <- here::here("data", "ppi_geneset_precalc_STRING/")
# output_dir <- "Users/xgh849/OneDrive - University of Copenhagen/Documents/RHAPSODY/customized_report/data/ppi_geneset_precalc/"
# dir.create(output_dir)

ppi_data <- readRDS(str_c(data_dir, "string_ppi_data_prepped.Rds"))

ref_exp <- read_tsv(str_c(data_dir, "data_BCDP_rhapsody_prepped.tsv"), col_types = cols(.default = "c")) %>% 
  distinct(uniprot) 


# ~ Gene sets ----
ORA_features <- tibble(database = c("geneontology_Biological_Process_noRedundant",
                                    "geneontology_Cellular_Component_noRedundant",
                                    "geneontology_Molecular_Function_noRedundant",
                                    "pathway_KEGG",
                                    "pathway_Reactome",
                                    "pathway_Wikipathway",
                                    "disease_GLAD4U",
                                    "drug_GLAD4U"),
                       database_formatted_name = c("Biological process GO-terms",
                                                   "Cellular component GO-terms",
                                                   "Molecular function GO-terms",
                                                   "KEGG pathways",
                                                   "Reactome pathways",
                                                   "Wikipathways",
                                                   "GLAD4U diseases",
                                                   "GLAD4U drugs"
                       )) %>% 
  mutate(database = factor(database))



# PPI ----

# ~ Interactions ---- 
ppi_names <- ppi_data %>% pluck("ppi_names", 1)
ppi_AB_all <- ppi_data %>% pluck("ppi_AB_all", 1)

tic()
interactions_alltargets <- ppi_names %>% 
  inner_join(ppi_AB_all, ., by = c("A" = "uniprot")) %>%
  group_by(gene_name) %>% 
  nest(hit_data = c(A, B, combined_score)) %>%
  left_join(ppi_names, ., by = "gene_name") %>% 
  bind_cols(., ppi_AB_all %>% nest(full_data = everything())) %>% 
  mutate(full_data = map2(full_data, hit_data, ~.x %>% filter(A %in% .y$B & B %in% .y$B))) %>% 
  rename(int_prim = hit_data, 
         int_sec = full_data) %>% 
  filter(!is.na(uniprot))
toc()


# ~ Nodes and edges ----
tic() 

nodes_edges_alltargets <- interactions_alltargets %>%
  mutate(.,
         edges_prim = map(int_prim, ~.x %>%
                            select(A, B, combined_score)),
         edges_sec = int_sec,
         edges_all = map2(edges_prim, edges_sec, ~bind_rows(.x, .y) %>% 
                            distinct_all() %>% 
                            mutate(confint_cutoff = if_else(combined_score > 0.7, TRUE, FALSE)) %>% 
                            rename(from = A, to = B)
         ),
  ) %>% 
  select(uniprot, 
         gene_name, 
         edges_all, 
  )

toc()

nodes_edges_alltargets %>% 
  write_rds(str_c(output_dir, "nodes_edges_alltargets.rds"))


# nodes_edges_alltargets <- read_rds(str_c(output_dir, "nodes_edges_alltargets.rds"))

# Gene set enrichment ---- 

# ~ Data set ----
ref_exp_uniprot <- ref_exp %>% distinct(uniprot) %>% filter(!is.na(uniprot)) #%>% nest(refFile = everything())

tic()
alltargets_oradata <- ORA_features %>%
  tidyr::crossing(., nodes_edges_alltargets) %>%  
  group_by(database) %>%
  mutate(edges_all = map2(edges_all, uniprot, ~.x %>%
                            filter(from == .y) %>%
                            select(-c(combined_score, confint_cutoff)) %>%
                            pivot_longer(names_to = "status", values_to = "uniprot", cols = everything()) %>%
                            distinct(uniprot)
  )
  )
toc()

alltargets_oradata %>% 
  write_rds(str_c(output_dir, "alltargets_oradata.rds"))

# alltargets_oradata <- read_rds(str_c(output_dir, "alltargets_oradata.rds"))

alltargets_oradata_refmatch <- alltargets_oradata %>% 
  inner_join(., ref_exp_uniprot)


# ~ PPI interactions ----

# ~~ Run WebGestalt ----


tic()
for(i in ORA_features$database){
  plan(multisession, workers = 7)
  
  ppi_interactions_data_func <- alltargets_oradata_refmatch %>%
    filter(database == i) %>%
    mutate(enrichAnalysis = future_map2(edges_all, database,  ~WGR_run_ppi_safe(data = .x, database = .y)))
  
  plan(sequential)
  ppi_interactions_data_func <- ppi_interactions_data_func %>% 
    mutate(enrichResult = map(enrichAnalysis, ~.x$result %>% as_tibble()),
           enrichResult_Prep = map2(enrichResult, uniprot, ~WGR_prep_ppi(.x, .y)))
  
  ppi_interactions_data_func %>%
    select(-c(edges_all, enrichResult)) %>%
    write_rds(str_c(output_dir, i, ".Rds"))
  
}
toc()




for(i in ORA_features$database){
  read_rds(str_c(output_dir, i, ".Rds")) %>% 
    select(-enrichAnalysis) %>% 
    unnest(enrichResult_Prep) %>%
    write_tsv(str_c(output_dir, i, "_only_enriched", ".tsv"))
  }


full_data <- NULL

for (i in ORA_features$database){
  data <- read_tsv(str_c(output_dir, i, "_only_enriched", ".tsv"))
  full_data <- full_data %>% 
    bind_rows(., data)
  
}

full_data %>% 
  write_tsv(str_c(output_dir, "all_databases_only_enriched", ".tsv"))


