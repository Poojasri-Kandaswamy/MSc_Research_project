#alternatively get this info by just converting your _full dataframe in your R environmnet that was made when using deseq2
# in my case its called heladiff_full
volcanodf <- Treatment_Hpx_vs_Nx_full
View(volcanodf)
#remove columns to be left with external_gene_name
volcanodf <- volcanodf[c("log2FoldChange","FDR","external_gene_name")] 
View(volcanodf)
colnames(volcanodf) <- c("log2FoldChange","FDR","ID")

#make volcano with point labels
#note ylim values are the scale of the y axis (log10 fdr), change to match your data
#note xlim values are the scale of the x axis (log2 FC), change to match your data
#SOME OTHER USEFUL CUTOMISABLE OPTIONS
#colAlpha is transparency of the points, 1 low, 0.5 medium, 0 = high
#pointSize is the size of the points, I think higher number = higher point size (package info will tell you more)
EnhancedVolcano(volcanodf,lab = volcanodf$ID,x = "log2FoldChange",y = "FDR",pCutoff = 0.049,FCcutoff = 0.579,
                pointSize = 3, col=c("grey", "grey", "grey", "blue"),
                colAlpha = 1, border = "full",borderWidth = 1.5, labSize = 5,
                selectLab = c("X"),
                borderColour = "black",gridlines.major = FALSE,gridlines.minor = FALSE, xlab = bquote(~Log[2]~ "fold change"),
                ylab = bquote(~-Log[10]~ "FDR"),legendLabels=c("NS","NS","NS","S"), xlim = c(-7, 13), ylim = c(0, 350),legendPosition = 'right')




#make volcano with point labels (in this example CA9 and BNIP3L are labelled)
EnhancedVolcano(volcanodf,lab = volcanodf$ID,x = "log2FoldChange",y = "FDR",pCutoff = 0.049,FCcutoff = 0.579,
                pointSize = 3, col=c("grey", "grey", "grey", "blue"),
                colAlpha = 1, border = "full",borderWidth = 1.5, labSize = 5,
                selectLab = c('CA9','BNIP3L'),
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                borderColour = "black",gridlines.major = FALSE,gridlines.minor = FALSE, xlab = bquote(~Log[2]~ "fold change"),
                ylab = bquote(~-Log[10]~ "FDR"),legendLabels=c("NS","NS","NS","S"), xlim = c(-7, 13), ylim = c(0, 350),legendPosition = 'right')

#labsize is label size
#draw connectors is having a line for label connecing to point turning to FALSE removes it
#widthConnectors is the width of the connectors
library(ggplot2)
library(EnhancedVolcano)

#make volcano plots with up and down genes as different colors 
keyvals <- ifelse(volcanodf$log2FoldChange <= -0.579 & volcanodf$FDR<0.049, 'royalblue',  ifelse(volcanodf$log2FoldChange >=0.579 & volcanodf$FDR<0.049, 'red2', 'grey'))
names(keyvals)[keyvals == 'red2'] <- 'Up'
names(keyvals)[keyvals == 'grey'] <- 'NS'
names(keyvals)[keyvals == 'royalblue'] <- 'Down'
up_genes <- sum(volcanodf$log2FoldChange >=0.579 & volcanodf$FDR<0.049)
down_genes <- sum(volcanodf$log2FoldChange <= -0.579 & volcanodf$FDR<0.049)

p <- EnhancedVolcano(volcanodf,lab = volcanodf$ID,x = "log2FoldChange",y = "FDR",pCutoff = 0.049,FCcutoff = 0.579,
                pointSize = 3, colCustom = keyvals,
                colAlpha = 1, border = "full",borderWidth = 1.5, labSize = 5,
                selectLab = c("X"),
                borderColour = "black",gridlines.major = FALSE,gridlines.minor = FALSE, xlab = bquote(~Log[2]~ "fold change"),
                ylab = bquote(~-Log[10]~ "FDR"),legendLabels=c("Up","NS","Down"), xlim = c(-7, 13), ylim = c(0, 350),legendPosition = 'right', title = "u87")
p + 
  annotate("text", x = 8, y = 280, label = paste(up_genes), size = 5) +
  annotate("text", x = -5, y = 280, label = paste(down_genes), size = 5)


#make volcano plot with selective gene labels 
EnhancedVolcano(volcanodf,lab = volcanodf$ID,x = "log2FoldChange",y = "FDR",pCutoff = 0.049,FCcutoff = 0.579,
                pointSize = 3, colCustom = keyvals,,
                colAlpha = 1, border = "full",borderWidth = 1.5, labSize = 5,
                selectLab = c('CA9','BNIP3L'),
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                borderColour = "black",gridlines.major = FALSE,gridlines.minor = FALSE, xlab = bquote(~Log[2]~ "fold change"),
                ylab = bquote(~-Log[10]~ "FDR"),legendLabels=c("Up","NS","Down"), xlim = c(-7, 13), ylim = c(0, 350),legendPosition = 'right')

