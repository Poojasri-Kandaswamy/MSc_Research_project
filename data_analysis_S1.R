test <- read.delim("~/Documents/feature_counts/501featurecount.txt",skip=1,head=T)
test <- test[-c(2,3,4,5,6)]
colnames(test) <- c("ensembl_gene_id","Nx_N1","Nx_N2","Nx_N3","Hpx_N1","Hpx_N2","Hpx_N3")
write.csv(test,"501_countdata.csv", row.names=FALSE)
countdata <- as.matrix(read.csv("501_countdata.csv",row.names="ensembl_gene_id"))
coldata <- as.matrix(read.csv("~/Documents/coldata.csv",row.names=1))
coldata

#setup for differential expression analysis using the DESeqDataSetFromMatrix function of deseq2
dds501 <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design= ~ Treatment)


#make heatmap 
library(pheatmap)
library(RColorBrewer)
rld1 <- rlog(dds501, blind=FALSE)
sampleDists <- dist(t(assay(rld1)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld1$Treatment,rld1$Names, sep="-")
colnames(sampleDistMatrix) <- paste(rld1$Treatment,rld1$Names, sep="-")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#make PCA
rld1 <- rlog(dds501, blind=FALSE)
pcaData <- plotPCA(rld1, intgroup = c( "Treatment", "Names"), returnData = TRUE)  
percentVar <- round(attr(pcaData, "percentVar")) 
ggplot(pcaData, aes(x = PC1, y = PC2, color = factor(Treatment), label = factor(Names))) + 
  geom_point(size =4) + 
  scale_color_manual(values=c("grey", "blue")) +
  theme_bw() + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  ggtitle("Title") + 
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())


#Get normalised gene counts for each gene from each replicate
dds501<- estimateSizeFactors(dds501)
count=as.data.frame(counts(dds501, normalized=TRUE))
count$ensembl_gene_id <- row.names(count)
row.names(count) <- NULL
write.table(count, "~/Documents/501_DESEQ2_counts.txt", sep="\t", row.names=TRUE)

#specify which is control (reference group) using relevel function of deseq2
dds501$Treatment <- relevel(dds501$Treatment, ref = "Nx")



#run differential expression analysis using DESeq function of deseq2
deseq2501<- DESeq(dds501)

#check coefs and use the in next step
resultsNames(deseq2501)

#generate normalized results using lfcShrink function of deseq2
#you can see the coef name from the above lines result 
Treatment_Hpx_vs_Nx_diff <-  lfcShrink(deseq2501, coef="Treatment_Hpx_vs_Nx", type="apeglm")


#convert this file to a dataframe
Treatment_Hpx_vs_Nx_raw <- as.data.frame(Treatment_Hpx_vs_Nx_diff)


#we want the row names (ensembl_gene_id) to be a column in our data frame 
Treatment_Hpx_vs_Nx_raw$ensembl_gene_id <- row.names(Treatment_Hpx_vs_Nx_raw)


#now we can delete row names
row.names(Treatment_Hpx_vs_Nx_raw) <- NULL
write.table(Treatment_Hpx_vs_Nx_raw, "Treatment_Hpx_vs_Nx_raw.txt", sep="\t", row.names=FALSE)

# we only care about the log2FoldChange column and the padj (also known as FDR) column and the ensembl_gene_id column
#so we can delete the other columns (columns 1,3 and 4 )
Treatment_Hpx_vs_Nx_raw <- Treatment_Hpx_vs_Nx_raw[-c(1,3,4) ]


#we can rename padj column name to FDR (i always call it FDR)
colnames(Treatment_Hpx_vs_Nx_raw) <- c("log2FoldChange","FDR","ensembl_gene_id")


#an NA value for the FDR  means the seqeuncing experiment has captured enough information (sequencing reads) for those genes
#so we can delete those genes which have a NA value for FDR
test <-Treatment_Hpx_vs_Nx_raw[!is.na(Treatment_Hpx_vs_Nx_raw$FDR), ]


#this is now final full differential expression list 
#so can rename 
Treatment_Hpx_vs_Nx_full <- as.data.frame(test)

#we also want make upregulated and downregulated lists applying cut offs
#we typcally use <0.05 FDR and > -/+ 0.58 log2FoldChange
#Filter to make upregulated and downregulated lists
testfilter <- filter(Treatment_Hpx_vs_Nx_full, FDR < 0.0501)
Treatment_Hpx_vs_Nx_up <- filter(testfilter, log2FoldChange > 0.5799)
Treatment_Hpx_vs_Nx_down <-filter(testfilter, log2FoldChange < - 0.5799)

mart38_data <- read.delim("~/Documents/mart38_data.txt",head=T)


#merge ensemble gene info list and our differential expression lists
Treatment_Hpx_vs_Nx_full<-merge(Treatment_Hpx_vs_Nx_full, mart38_data, by.x="ensembl_gene_id", by.y="ensembl_gene_id")
Treatment_Hpx_vs_Nx_up<-merge(Treatment_Hpx_vs_Nx_up, mart38_data, by.x="ensembl_gene_id", by.y="ensembl_gene_id")
Treatment_Hpx_vs_Nx_down<-merge(Treatment_Hpx_vs_Nx_down, mart38_data, by.x="ensembl_gene_id", by.y="ensembl_gene_id")

#export data
write.table(Treatment_Hpx_vs_Nx_full, "~/Documents/501_Treatment_Hpx_vs_Nx_full.txt", sep="\t", row.names=FALSE)
write.table(Treatment_Hpx_vs_Nx_up, "~/Documents/501_Treatment_Hpx_vs_Nx_up.txt", sep="\t", row.names=FALSE)
write.table(Treatment_Hpx_vs_Nx_down, "~/Documents/501_Treatment_Hpx_vs_Nx_down.txt", sep="\t", row.names=FALSE)

