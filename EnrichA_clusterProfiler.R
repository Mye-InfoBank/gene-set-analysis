# Gen Set Enrichment Anlysis using clusterProfiler

# NOTE: You need run load_msigdb_datasets_updated.R to download
# gene annotation databases need fo enrichment analysis 


# BiocManager::install("clusterProfiler")
# install.packages("msigdbr")
library("clusterProfiler")
source("path/load_msigdb_datasets.R")

res_CP <- enricher(genes, # gene vector ID
            pvalueCutoff = 0.05, # adjusted pvalue cutoff on enrichment tests to report
            pAdjustMethod = "BH", # one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
            universe = NULL, # If missing, the all genes listed in the database (eg TERM2GENE table) will be used as bbackground genes. # 
            minGSSize = 10, # minimal size of genes annotated for testing
            maxGSSize = 500, # maximal size of genes annotated for testing
            qvalueCutoff = 0.2, # qvalue cutoff on enrichment tests to report as significant
            TERM2GENE = an_df_all_GO) # user input annotation of TERM TO GENE mapping, a data.frame of 2 column with term and gene. Example: GO annotations

# Show results
colnames(res_CP@result)
# [1] "ID"          "Description" "GeneRatio"   "BgRatio"    
# [5] "pvalue"      "p.adjust"    "qvalue"      "geneID"     
# [9] "Count"