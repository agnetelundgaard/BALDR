# Get data files
library(R.utils)

data_dir <- here::here("data/")

# Primary and secondary uniprot accession numbers
download.file("https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/docs/sec_ac.txt", 
              destfile = str_c(data_dir, "UniProt_mapping_files/", "sec_ac.txt"))

# Uniprot and ENSG mapping
human_idmapping <- str_c(data_dir, "UniProt_mapping_files/", "HUMAN_9606_idmapping.dat.gz")
download.file("https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz",  destfile = human_idmapping)
gunzip(human_idmapping)

file.info("data/UniProt_mapping_files/HUMAN_9606_idmapping.dat")

mouse_idmapping <- str_c(data_dir, "UniProt_mapping_files/", "MOUSE_10090_idmapping.dat.gz")
download.file("https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/idmapping/by_organism/MOUSE_10090_idmapping.dat.gz", destfile = mouse_idmapping)
gunzip(mouse_idmapping)


# Ortholog information from Phylomedb
human_protein_name <- str_c(data_dir, "Phylomedb_data/", "phylomedb_human_proteinname.txt.gz")
download.file("ftp://phylomedb.org/phylomedb//phylomes/phylome_0514/all_protein_names.txt.gz", destfile = human_protein_name)
gunzip(human_protein_name)

human_ortho_info <- str_c(data_dir, "Phylomedb_data/", "phylomedb_human_ortholog_info.txt.gz")
download.file("ftp://phylomedb.org/phylomedb//phylomes/phylome_0514/orthologs.txt.gz", destfile = human_ortho_info)
gunzip(human_ortho_info)

mouse_protein_name <- str_c(data_dir, "Phylomedb_data/", "phylomedb_mouse_proteinname.txt.gz")
download.file("ftp://phylomedb.org/phylomedb//phylomes/phylome_0518/all_protein_names.txt.gz", destfile = mouse_protein_name)
gunzip(mouse_protein_name)


# ensembl historic mapping files
human_ensembl <- str_c(data_dir, "ensembl/", "human_stable_id_event.txt.gz")
download.file("http://ftp.ensembl.org/pub/release-104/mysql/homo_sapiens_core_104_38/stable_id_event.txt.gz", destfile = human_ensembl)
gunzip(human_ensembl)

mouse_ensembl <- str_c(data_dir, "ensembl/", "mouse_stable_id_event.txt.gz")
download.file("http://ftp.ensembl.org/pub/release-104/mysql/mus_musculus_core_104_39/stable_id_event.txt.gz", destfile = mouse_ensembl)
gunzip(mouse_ensembl)


# STRING database
string_data <- str_c(data_dir, "ppi_data/", "9606.protein.physical.links.detailed.v11.5.txt.gz")

download.file("https://stringdb-static.org/download/protein.physical.links.detailed.v11.5/9606.protein.physical.links.detailed.v11.5.txt.gz" , destfile = string_data)
gunzip(string_data)

# HPA data
hpa_data <- str_c(data_dir, "HPA_blood_proteins/", "proteinatlas.tsv.zip")

download.file("https://www.proteinatlas.org/download/proteinatlas.tsv.zip", destfile = hpa_data)
unzip(hpa_data, exdir = str_c(data_dir, "HPA_blood_proteins/"))


