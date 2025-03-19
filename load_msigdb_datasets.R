# Load necessary libraries
library(msigdbr)  # MSigDB gene set retrieval

# Retrieve MSigDB gene sets for Homo sapiens
# Hallmark gene sets
an_df_hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol) %>% 
  as.data.frame()

# BIOCARTA pathway gene sets
an_df_biocarta <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:BIOCARTA") %>% 
  dplyr::select(gs_name, gene_symbol) %>% 
  as.data.frame()

# KEGG pathway gene sets
an_df_kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>% 
  dplyr::select(gs_name, gene_symbol) %>% 
  as.data.frame()

# Reactome pathway gene sets
an_df_reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>% 
  dplyr::select(gs_name, gene_symbol) %>% 
  as.data.frame()

# Gene Ontology (GO) gene sets
an_df_all_GO <- msigdbr(species = "Homo sapiens", category = "C5") %>% 
  dplyr::select(gs_name, gene_symbol) %>% 
  as.data.frame()

# Transcription Factor Targets (C3 TFT:GTRD)
an_df_C3_TFT <- msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:GTRD") %>% 
  dplyr::select(gs_name, gene_symbol) %>% 
  as.data.frame()
an_df_C3_TFT$gs_name <- paste0("TF_", as.character(an_df_C3_TFT$gs_name))

# Gene Ontology - Biological Process (GO:BP)
an_df_GOBP <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>% 
  dplyr::select(gs_name, gene_symbol) %>% 
  as.data.frame()

# Gene Ontology - Cellular Component (GO:CC)
an_df_GOCC <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC") %>% 
  dplyr::select(gs_name, gene_symbol) %>% 
  as.data.frame()

# Gene Ontology - Molecular Function (GO:MF)
an_df_GOMF <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF") %>% 
  dplyr::select(gs_name, gene_symbol) %>% 
  as.data.frame()

# Cell Type gene sets (C8)
an_df_Celltype <- msigdbr(species = "Homo sapiens", category = "C8") %>% 
  dplyr::select(gs_name, gene_symbol) %>% 
  as.data.frame()

# Human Phenotype Ontology (HPO) gene sets
an_df_HPO <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "HPO") %>% 
  dplyr::select(gs_name, gene_symbol) %>% 
  as.data.frame()


gene_annodf_list <- list(
  hallmark = an_df_hallmark,
  biocarta = an_df_biocarta,
  kegg = an_df_kegg,
  reactome = an_df_reactome,
  GO = an_df_all_GO,
  GO_BP = an_df_GOBP,
  GO_CC = an_df_GOCC,
  GO_MF = an_df_GOMF,
  celltype = an_df_Celltype,
  HPO = an_df_HPO,
  TFT = an_df_C3_TFT
)
# Combine selected gene annotation data frames into a single dataset
# Uncomment additional gene sets if needed
db_of_choice = "GO"
gene_annodf <- gene_annodf_list[[db_of_choice]]

gene_annodf <- gene_annodf[, c("gs_name", "gene_symbol")]
colnames(gene_annodf) <- c("TERM", "GENE")

# Display a message indicating how to use the generated gene annotation data
cat("Gene annotation data frame created as 'gene_annodf'.\n")
cat("To include GO:BP, GO:CC, or GO:MF, uncomment relevant lines in the script.\n")



##### Perform GSEA or enricher with Clusterprofiler with TERM2GENE dataframe from the list#######