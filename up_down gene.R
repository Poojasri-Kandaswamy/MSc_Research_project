library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(msigdbr)
library(dplyr)

#accessing the hallmark gene set
m_df <- msigdbr(species = "Homo sapiens")
head(m_df)

msigdbr_collections <- as.data.frame(msigdbr_collections())
View(msigdbr_collections)

hallmark_gene_set <- msigdbr(species = "Homo sapiens", collection = "H")
View(hallmark_gene_set)

hallmark_gene_set_2 <- hallmark_gene_set %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
View(hallmark_gene_set_2)



godata <- Treatment_Hpx_vs_Nx_down
View(godata)
godata <- godata[c("log2FoldChange","external_gene_name")] 
View(godata)

data_vector <- as.character(godata$external_gene_name)
View(data_vector)

#run ORA analysis with cluster profiler using enricher function
enricher_test = enricher(
  data_vector, #gene list of interest
  pAdjustMethod = "BH",
  universe = NULL,
  minGSSize = 10,
  maxGSSize = 500,
  gson = NULL,
  TERM2GENE = hallmark_gene_set_2, #geneset
  TERM2NAME = NA) #background list/gene universe if you have one) 

fit <- plot(barplot(enricher_test, showCategory = 10))


# Get gene IDs for RNA splicing GO term (GO:0008380)
genes <-AnnotationDbi::select(org.Hs.eg.db, keys=keys(org.Hs.eg.db),
                              columns=c("SYMBOL","GO"), keytype="ENTREZID")

# Filter genes by RNA splicing GO term
splicing_genes <- genes[grep("GO:0008380", genes$GO), ]

# Print the list of genes
head(splicing_genes)
View(splicing_genes)
write.csv(splicing_genes,"u87_splicing_genes_down.csv",row.names = FALSE)
