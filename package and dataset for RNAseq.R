#install R bioconductor core packages 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")



#install specific R bioconductor packages
BiocManager::install(c("DESeq2","ggplot2","pheatmap","biomaRt", "GenomicFeatures","EnhancedVolcano","GenomicAlignments", "GenomicRanges", "IRanges","GeneOverlap","dp
lyr","txdbmaker","apeglm","clusterProfiler","org.Hs.eg.db","AnnotationDbi"))




#upgrade packages
BiocManager::install()


#load all packages each time a new R session is started
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(biomaRt)
library(GenomicFeatures)
library(EnhancedVolcano)
library(GenomicAlignments)
library(GenomicRanges)
library(IRanges)
library(GeneOverlap)
library(dplyr)
library(txdbmaker)
library(apeglm)
library(clusterProfiler)
library(AnnotationDbi)



#get ensembl databases from biomart for gene information 
mart38 <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl", host="https://www.ensembl.org")

#if the above doesnt work try one of the below
mart38 <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl", host="uswest.ensembl.org")
mart38 <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl", host="useast.ensembl.org")
mart38 <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl", host="asia.ensembl.org")


#get ensembl gene info from biomart for all ensembl genes 
mart38_data<-getBM(attributes=c("external_gene_name",
                                "ensembl_gene_id",
                                "gene_biotype",
                                "chromosome_name",
                                "start_position",
                                "end_position",
                                "strand"),
                   values="ensembl_gene_id",
                   mart=mart38)

#output ensembl gene database 
write.table(mart38_data, "/users/poojasrikandaswamy/Downloads/mart38_data.txt", sep="\t", row.names=FALSE)

