### https://github.com/bli25/RSEM_tutorial
rm(list = ls())

# Load packages
library(edgeR)
library(dplyr)
library(gplots)
library(RColorBrewer)
library(plyr)
library(pheatmap)
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(plotly)
library(limma)
library(ggfortify)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)

# A short function for outputting the tables
knitr_table <- function(x) {
  x %>% 
    knitr::kable(format = "html", digits = Inf, 
                 format.args = list(big.mark = ",")) %>%
    kableExtra::kable_styling(font_size = 15)
}

# 
# meta_table <- read.delim("files/samples_meta_table.csv", sep = ",")
# rownames(meta_table) <- meta_table$Sample.ID
# x = read.delim("data/12_samples.RSEM.isoforms.counts.matrix")
# rownames(x) <- x$X
# 
# colnames(x) <- gsub(".isoforms.results", "", colnames(x))
# x = x[,-1]
# head(x)

# y <- DGEList(counts=x)
# 
# meta_table = meta_table[rownames(y$samples),]
# y$samples$group = meta_table$Treatment
# 
# keep <- filterByExpr(y)
# y <- y[keep, , keep.lib.sizes=FALSE]
# dim(y)
# head(y)
# 
# y <- calcNormFactors(y)
# 
# plotMDS(y)
# 
# design <- model.matrix(~meta_table$Treatment+meta_table$Timepoint)
# 
# y <- estimateDisp(y, design)
# 
# 
# 
# fit <- glmQLFit(y, design)
# 
# qlf <- glmQLFTest(fit, coef=2:3)
# topTags(qlf)
# 
# Group <- factor(paste(targets$Treat,targets$Time,sep="."))
# cbind(targets,Group=Group)


###############
# gtf = read.delim("/media/data01/common/ReferenceDB/GRCm39/Annotations/GRCm39_GENCODE_VM27_refseq.bed", header = F)

meta_table <- read.delim("files/84_samples_meta_data.csv", sep = ",")
#meta_table <- read.delim("files/30_samples_Psilo_saline_meta_data.csv", sep = ",")
meta_table$Sample.ID = gsub("\\-", "", meta_table$Sample.ID)
rownames(meta_table) <- meta_table$Sample.ID
meta_table = meta_table[order(meta_table$CellType, meta_table$Treatment, meta_table$Timepoint, decreasing = T),]
meta_table$Treatment[meta_table$Treatment == "Saline"] = "Ctrl"
Group <- factor(paste(meta_table$Treatment,meta_table$Timepoint,sep="."))
meta_table = cbind(meta_table,Group=Group)
Group2 <- factor(paste(meta_table$CellType,meta_table$Treatment,meta_table$Timepoint,sep="."))
meta_table = cbind(meta_table,Group2=Group2)
meta_table$Batch = paste0("B", meta_table$Batch)

### select cell type
#meta_table = meta_table[meta_table$CellType == "PV",]
#meta_table = meta_table[meta_table$CellType == "PV" & meta_table$Timepoint == '48hrs',]
meta_table = meta_table[meta_table$CellType == "PV" & meta_table$Timepoint == '1week',]
#meta_table = meta_table[meta_table$CellType == "PV" & meta_table$Timepoint == '1month',]

#meta_table = meta_table[meta_table$CellType == "Rbp4",]
#meta_table = meta_table[meta_table$CellType == "Rbp4" & meta_table$Timepoint == '48hrs',]
#meta_table = meta_table[meta_table$CellType == "Rbp4" & meta_table$Timepoint == '1week',]
#meta_table = meta_table[meta_table$CellType == "Rbp4" & meta_table$Timepoint == '1month',]
#meta_table$Batch = paste0("B", meta_table$Batch)


#meta_table = meta_table[meta_table$CellType == "Cux2ERT" & meta_table$Timepoint == '48hrs',]
#meta_table = meta_table[meta_table$CellType == "Cux2ERT" & meta_table$Timepoint == '1week',]
#meta_table = meta_table[meta_table$CellType == "Cux2ERT" & meta_table$Timepoint == '1month',]

meta_table = droplevels(meta_table)
grid::grid.newpage()
gridExtra::grid.table(meta_table, rows = NULL)


Count_data <- read.delim("data/84_samples_merged_RSEM_bowtie2.genes.counts.matrix")
#Count_data <- read.delim("data/30_samples_48h_psilo_merged_RSEM_bowtie2.genes.counts.matrix")
#colnames(Count_data) <- gsub(".genes.results", "", colnames(Count_data))
colnames(Count_data) <- gsub("\\_.*", "", colnames(Count_data))
colnames(Count_data) <- gsub("\\..*", "", colnames(Count_data))

Counts <- Count_data
rownames(Counts) <- Count_data$X
counts_IDs <- Count_data
Counts_only <- Count_data[, -1]  # create the table with only counts here
rownames(Counts_only) <- counts_IDs$X
Counts_only = Counts_only[,colnames(Counts_only) %in% meta_table$Sample.ID]
Counts_only = Counts_only[, meta_table$Sample.ID]


y <- DGEList(counts=Counts_only, genes = counts_IDs$X, group=meta_table$Group)

#HEATMAPS:
#design <- model.matrix(~meta_table$Group)
keep <- filterByExpr(y, min.count=10, large.n=10, min.prop=0.7)
#keep <- filterByExpr(y, large.n=4)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]

# keep <- rowSums(y$counts) > 60
# samples_g1 = meta_table$Sample.ID[meta_table$Treatment == 'Ctrl']
# samples_g2 = meta_table$Sample.ID[meta_table$Treatment == 'DOI']
# apply(y$counts[, samples_g2], 1, function(i) sum(i > 0))
#   
# y$counts['E330020D12Rik',]

# y$samples
# y <- normLibSizes(y)
# y$samples

logcounts <- cpm(y, normalized.lib.sizes = TRUE, log=TRUE)
head(logcounts)

logcounts['Dot1l',]

ZlogCPM <- t(scale(t(logcounts)))
head(ZlogCPM)

boxplot( t(logcounts[1:30,])   , ylim=c(0,12), las=2 , col="grey" , main="logCPM raw_data" ,cex=.2)
boxplot( t(ZlogCPM[1:30,]) , ylim=c(-4,4), las=2 , col="grey" , main="Z-score data" ,cex=.2)


####




####
counts_vst = vst(round(y$counts))

# Perform the PCA
pca_res <- prcomp(t(counts_vst), scale=TRUE)
pca_res <- prcomp(t(logcounts), scale=TRUE)

# pca_res_data = as.data.frame(pca_res$x)
# 
# pca_res_data = cbind(pca_res_data, treatment=meta_table[rownames(pca_res_data),]$Group2)
# pca_res_data = cbind(pca_res_data, Batch=factor(meta_table[rownames(pca_res_data),]$Batch))
# pca_res_data = cbind(pca_res_data, Sex=factor(meta_table[rownames(pca_res_data),]$Sex))
# 
# var_explained <- pca_res$sdev^2/sum(pca_res$sdev^2)
# var_explained[1:5]
# 
# ggplot(pca_res_data, aes(x=PC1,y=PC2)) +
#   geom_point(aes(color=treatment),size=4) +
#   theme_classic(base_size=17) +
#   labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
#        y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
#   theme(legend.position="top")
#   

autoplot(pca_res,
         data = meta_table, 
         colour="Treatment", 
         shape="Sex",
         size=5) +
  theme_classic2(base_size = 17) +
  geom_text_repel(aes(x=PC1, y=PC2, label=Sample.ID), box.padding = 0.8)


# autoplot(pca_res,
#          data = meta_table, 
#          colour="Timepoint", 
#          shape="Treatment",
#          size=5) +
#   theme_classic2(base_size = 17) +
#   geom_text_repel(aes(x=PC1, y=PC2, label=Sample.ID), box.padding = 0.8)

######## Gene selection
#Compute the mean, variance and cv for each gene and sort them in decreasing order
gene_stats <- data.frame(row.names = rownames(Counts_only))
means <- rowMeans(logcounts)
vars <- apply(logcounts,1,var)

#Plot the row means versus row means
#mypar(1,2,mar = c(5,3,2,1))
plot(x=means, y=vars, cex=.1 )

#Sort and select the top 500 highly variable genes from the data
vars <- sort( vars , decreasing=T)
top_var <- names( vars ) [1:500]
boxplot( t(logcounts[ top_var[1:15],])   , ylim=c(0,12), las=2 , col="grey" , main="logCPM (top var genes)" ,cex=.2)


####

#PC <-  prcomp( t( ZlogCPM[ top_var, ]) ) #Method1
#PC <-  prcomp( t( logcounts[ top_var, ]), center = TRUE, scale. = TRUE) #Method2
PC <-  prcomp( t( counts_vst[ top_var, ])) #Method2

plot(PC$x[,1] , PC$x[,2], cex=2, col=factor(meta_table$Group), xlab="PC1", ylab="PC2", pch=16, main="PCA", las=1)
text(PC$x[,1] , PC$x[,2], cex=.7, labels = paste0(meta_table$Sample.ID), pos = 3)

plot(PC$x[,3] , PC$x[,4], cex=2, col=factor(meta_table$Group), xlab="PC3", ylab="PC4", pch=16, main="PCA", las=1)
text(PC$x[,3] , PC$x[,4], cex=.7, labels = paste0(meta_table$Sample.ID), pos = 3)


autoplot(PC,
         data = meta_table, 
         colour="Treatment", 
         shape="Sex",
         size=5) +
  theme_classic2(base_size = 17) +
  geom_text_repel(aes(x=PC1, y=PC2, label=Sample.ID), box.padding = 0.8)


PC_sd <- setNames(PC$sdev , paste0("PC",1:length(PC$sdev)))
PC_var_expl <- ( PC_sd^2 ) / sum(PC_sd^2) * 100

# Which PCs explain at least 2% of variance?
barplot(PC_var_expl, las=2, ylab="% variance explained")
abline(h=2, lty=2)

#Compute sample correlations
sample_cor <- cor( logcounts )
round(sample_cor,4)

#Transform the scale from correlations
cor_distance <- -(sample_cor - 1)/2
round(cor_distance,4)

#Convert it to a distance object
d2 <- as.dist(cor_distance)
d2
h2 <- hclust(d2, method="complete")
plot( as.dendrogram(h2) , las=1, main="d=correlation\nh=complete")
points(1:ncol(Counts_only) ,rep(0,ncol(Counts_only)), pch= 16, cex=2, col=as.factor(meta_table$Batch[h2$order]))

###############
mod = model.matrix(~meta_table$Treatment)
mod = mod[, 1:2]
combat_alln_mean = ComBat(dat=y, batch = meta_table$Sex, mod = mod)
#combat_alln_mean = ComBat(dat=counts_rlog, batch = meta_table$Batch, mod = mod ,  par.prior = T, mean.only = T)
# 
pca_combat = prcomp(t(combat_alln_mean))

# autoplot(PC,
#          data = meta_table, 
#          colour="Timepoint", 
#          shape="Treatment",
#          size=5) +
#   theme_classic2(base_size = 17) +
#   geom_text_repel(aes(x=PC1, y=PC2, label=Sample.ID), box.padding = 0.8)

# 
autoplot(pca_combat,
         data = meta_table,
         colour="Treatment",
         shape="Sex",
         size=5) +
  theme_classic2(base_size = 17) +
  geom_text_repel(aes(x=PC1, y=PC2, label=Sample.ID), box.padding = 0.8)
# 
# 
# PC <-  prcomp( t( combat_alln_mean[ top_var, ])) #Method2
# 
# autoplot(PC,
#          data = meta_table, 
#          colour="Sex", 
#          shape="Treatment",
#          size=5) +
#   theme_classic2(base_size = 17) +
#   geom_text_repel(aes(x=PC1, y=PC2, label=Sample.ID), box.padding = 0.8)
# 


# pca_res_data = as.data.frame(pca_combat$x)
# 
# #pca_res_data = cbind(pca_res_data, treatment=meta_table[rownames(pca_res_data),]$Group2)
# pca_res_data = cbind(pca_res_data, Batch=factor(meta_table[rownames(pca_res_data),]$Batch))
# #pca_res_data = cbind(pca_res_data, Sex=factor(meta_table[rownames(pca_res_data),]$Sex))
# 
# var_explained <- pca_combat$sdev^2/sum(pca_combat$sdev^2)
# var_explained[1:5]
# 
# ggplot(pca_res_data, aes(x=PC1,y=PC2)) +
#   geom_point(aes(color=Batch),size=4) +
#   theme_classic(base_size=17) +
#   labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
#        y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
#   theme(legend.position="top")


#biasPlot(, "GC", log = T, ylim = c(-5,5), main = "Bias plot")

# autoplot(pca_combat,
#          data = meta_table, 
#          colour="CellType", 
#          size = 4,
#          shape="Timepoint") +
#   theme_classic2(base_size = 17) +
#   geom_text_repel(aes(x=PC1, y=PC2, label=Sample.ID), box.padding = 0.8)
# 
# y.corrected <- removeBatchEffect(y, batch=meta_table$Sex, group=meta_table$Group)
# 
# plotMDS(y,main="Original")
# plotMDS(y.corrected,main="Batch corrected")
# 
# y.corrected
# 
# 
# PC <-  prcomp( t( y.corrected)) #Method2
# 
# autoplot(PC,
#          data = meta_table, 
#          colour="Timepoint", 
#          shape="Treatment",
#          size=5) +
#   theme_classic2(base_size = 17) +
#   geom_text_repel(aes(x=PC1, y=PC2, label=Sample.ID), box.padding = 0.8)
# 
# autoplot(PC,
#          data = meta_table, 
#          colour="Sex", 
#          shape="Treatment",
#          size=5) +
#   theme_classic2(base_size = 17) +
#   geom_text_repel(aes(x=PC1, y=PC2, label=Sample.ID), box.padding = 0.8)
# 


#####

# len <- read.csv("/media/data01/common/ReferenceDB/GRCm39/Annotations/GRCm39_GENCODE_VM27_refseq_lengths.cvs")
# head(len)
# 
# # Subset data for transcripts that are in the transcript length reference file
# Counts_only <- Counts_only[rownames(Counts_only) %in% rownames(len),]
# # Sort data by genes in hnsc
# len <- len[rownames(Counts_only),]
# 
# # Divide each length by 1000 to get number of kilobases
# len$KB <- len$Length/1000
# 
# # Create new hnsc object, don't overwrite the previous
# countsTPM <- Counts_only
# 
# # Divide each gene by transcript length
# expPKB <- apply( countsTPM, 2, function(x){ x / len$length } )
# 
# # Divide by the transcript length
# exprs(countsTPM) <- apply( expPKB, 2, function(x) { x / sum(x) * 1E6} )
# exprs(countsTPM)[1:5, 1:5]

genes = c("Rab5a", "Syp", "Eif5", "Scn1b", "Bdnf", "Syngap1", "Nf1", "Dlg4", "Grin2a", "Kmt2d", "Grin1")

logcounts_sel = ZlogCPM[rownames(logcounts) %in% genes,]
logcounts_sel = logcounts_sel[genes,]
dim(logcounts_sel)

logcounts_sel = logcounts_sel[, meta_table$Sample.ID]

# log2matrix <- log2(assay(vsd))
# 
# x = t(scale(t(logcounts_GSE161626[-which(colnames(logcounts_GSE161626) %in% c("Row.names", "external_gene_name"))])))
# log2matrix_scaled = range01(x)

# 
# pheatmap(logcounts_sel,
#          cluster_cols = F,
#          cluster_rows = F,
#          annotation_col = meta_table[,c(2, 3)],
#          #clustering_distance_rows=sampleDists,
#          #clustering_distance_cols=sampleDists,
#          main=paste0(nrow(logcounts_sel), " genes across conditions"),
#          #color = colorRampPalette(c("white", "#af8dc3", "#762a83"))(50),
#          show_rownames = T)

breaksList = seq(-2, 2, by = 0.1)

my_colour = list(
  Treatment = c(DOI = "#5977ff", Ctrl = "#f74747"),
  CellType = c(Rbp4 = "#f6e8c3", PV = "#9e82ed", Cux2ERT = "#5ab4ac"),
  Sex = c(`F` = "#d73027", `M` = "#4575b4"),
  Batch = c(`B1` = "#1b7837", `B2` = "#af8dc3"),
  Timepoint = c(`1week` = "#e89829", `1month` = "#82ed82", `48hrs` = "brown")
)


pheatmap(logcounts_sel,
         cluster_cols = F,
         cluster_rows = T,
         annotation_col = meta_table[,c(2:6)],
         annotation_colors = my_colour, 
         show_rownames = T, 
         color=colorRampPalette(c("#4575b4", "white", "#d73027"))(length(breaksList)),
         breaks = breaksList,
         main=paste0(nrow(logcounts_sel), " genes across conditions (z-score)"))



####



# Perform the PCA
# pca_res <- prcomp(t(logcounts), scale=TRUE)
# pca_res_data = as.data.frame(pca_res$x)
# 
# pca_res_data = cbind(pca_res_data, treatment=meta_table[rownames(pca_res_data),]$Group)
# 
# var_explained <- pca_res$sdev^2/sum(pca_res$sdev^2)
# var_explained[1:5]
# 
# ggplot(pca_res_data, aes(x=PC1,y=PC2)) +
#   geom_point(aes(color=treatment),size=4) +
#   theme_classic(base_size=17) +
#   labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
#        y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
#   theme(legend.position="top")

#####
# 
# prin_comp = prcomp(t(logcounts), rank. = 3)
# prin_comp = prcomp(t(log2(y$counts + 1)), rank. = 3)
# 
# 
# components <- prin_comp[["x"]]
# components <- data.frame(components)
# components$PC2 <- -components$PC2
# components$PC3 <- -components$PC3
# components = cbind(components, treatment=meta_table[rownames(components),]$Group)
# 
# tot_explained_variance_ratio <- summary(prin_comp)[["importance"]]['Proportion of Variance',]
# tot_explained_variance_ratio <- 100 * sum(tot_explained_variance_ratio)
# 
# #tit = 'Total Explained Variance = 99.48'
# 
# fig <- plot_ly(components, 
#                x = ~PC1, 
#                y = ~PC2, 
#                z = ~PC3, 
#                color = ~meta_table$Treatment, 
#                text = rownames(components), 
#                colors = c('#636EFA','#EF553B','#00CC96', 'orange') ) %>%
#   add_markers(size = 17)
# 
# 
# fig <- fig %>%
#   layout(
#     title = "PCA genes",
#     scene = list(bgcolor = "#e5ecf6")
#   )
# 
# fig
# 
# 
# ###############
# 
# 
# pca_combat = prcomp(t(logcounts))
# 
# #pdf("results/PCA_Michael_samples.pdf", width = 5.5, height = 4)
# ggplot(data= data.frame(pca_combat$x[,1:2]),aes(x = PC1, y = PC2, shape=as.character(meta_table$Batch), color= as.factor(meta_table$Group))) +
#   ggtitle("log CPM counts") +
#   geom_point(size = 3) +
#   theme_classic(base_size = 17) +
#   scale_colour_manual(name = "Group", labels = c("Ctrl.1month", "Ctrl.1week", "Ctrl.48hrs", "DOI.1month", "DOI.1week", "DOI.48hrs"), values = c("blue", "#377EB8","darkgreen", "red", "orange", "brown")) +
#   scale_shape_manual(name = "Batch", labels = c("1", "2"), values = c(19, 17)) +
# 
#   xlab(paste0("PC1 ", prettyNum(summary(pca_combat)$importance[2,1]*100,
#                                 digits = 2, decimal.mark = "."), "%")) +
#   ylab(paste0("PC2 ", prettyNum(summary(pca_combat)$importance[2,2]*100,
#                                 digits = 2, decimal.mark = "."), "%"))+
#   theme(plot.title = element_text(hjust = 0.5))
# #dev.off()
# 
# 
# pca_combat = prcomp(t(Counts_only))
# 
# #pdf("results/PCA_Michael_samples.pdf", width = 5.5, height = 4)
# ggplot(data= data.frame(pca_combat$x[,1:2]),aes(x = PC1, y = PC2, shape=as.character(meta_table$Batch), color= as.factor(meta_table$Group))) +
#   ggtitle("raw counts") +
#   geom_point(size = 3) +
#   theme_classic(base_size = 17) +
#   scale_colour_manual(name = "Group", labels = c("Ctrl.1month", "Ctrl.1week", "Ctrl.48hrs", "DOI.1month", "DOI.1week", "DOI.48hrs"), values = c("blue", "#377EB8","darkgreen", "red", "orange", "brown")) +
#   scale_shape_manual(name = "Batch", labels = c("1", "2"), values = c(19, 17)) +
#   
#   xlab(paste0("PC1 ", prettyNum(summary(pca_combat)$importance[2,1]*100,
#                                 digits = 2, decimal.mark = "."), "%")) +
#   ylab(paste0("PC2 ", prettyNum(summary(pca_combat)$importance[2,2]*100,
#                                 digits = 2, decimal.mark = "."), "%"))+
#   theme(plot.title = element_text(hjust = 0.5))
# 
# 
# #####################
# 
# 
# # # Load the Bioconductor package "TxDb.Mmus.UCSC.mm10.knownGene"
# # library(TxDb.Mmusculus.UCSC.mm10.knownGene)
# # 
# # # Get the list of genes
# # genes <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
# # 
# # # Get the lengths of the genes
# # gene_lengths <- width(genes)
# # 
# # # Create a table of gene lengths
# # gene_lengths_table <- data.frame(Gene = genes$gene_id, Length = gene_lengths)
# # 
# # # Print the table
# # print(gene_lengths_table)
# 
# ####
# 
# 
# pc_scores <- prin_comp[["x"]]
# 
# pc_scores <- pc_scores %>% 
#   # convert to a tibble retaining the sample names as a new column
#   as_tibble(rownames = "sample")
# 
# # print the result  
# head(pc_scores)
# 
# pc_scores %>% 
#   # create the plot
#   ggplot(aes(x = PC1, y = PC2, label=sample)) +
#   geom_point() +
#   geom_text(hjust=0, vjust=0) +
#   theme_classic(base_size = 17)
# 
# cor_matrix = cor(logcounts)
# 
my_colour = list(
  #Treatment = c(Psilocybin = "#5977ff", Ctrl = "#f74747"),
  Treatment = c(DOI = "#5977ff", Ctrl = "#f74747"),
  CellType = c(Rbp4 = "#f6e8c3", PV = "#9e82ed", Cux2ERT = "#5ab4ac"),
  Sex = c(`F` = "#d73027", `M` = "#4575b4"),
  #Batch = c(`1` = "#1b7837", `2` = "#af8dc3"),
  Batch = c('B2' = "#1b7837", 'B3' = "#af8dc3"),
  Timepoint = c(`1week` = "#e89829", `1month` = "#82ed82", `48hrs` = "brown")
)
# 
# 
# pheatmap(cor_matrix,
#          annotation_col = meta_table[, c(2:6)],
#          annotation_colors = my_colour, 
#          show_rownames = F,
#          color=colorRampPalette(c("#4575b4", "white", "#d73027"))(50),
#          clustering_distance_cols = "euclidean",
#          clustering_distance_rows = "euclidean",
#          clustering_method = "ward.D2")


####

gene_stats <- data.frame(row.names = rownames(Counts_only))
means <- rowMeans(logcounts)
vars <- apply(logcounts,1,var)

#Plot the row means versus row means
#mypar(1,2,mar = c(5,3,2,1))
plot(x=means, y=vars, cex=.1 )

#Sort and select the top 500 highly variable genes from the data
vars <- sort( vars , decreasing=T)
top_var <- names( vars ) [1:300]
boxplot( t(logcounts[ top_var[1:15],])   , ylim=c(0,12), las=2 , col="grey" , main="logCPM (top var genes)" ,cex=.2)

# var_genes <- apply(logcounts, 1, var)
# head(var_genes)
# length(var_genes)
# # Get the gene names for the top 500 most variable genes
# select_var <- names(sort(var_genes, decreasing=TRUE))[1:50]
# head(select_var)

highly_variable_lcpm <- ZlogCPM[top_var, ]
#highly_variable_lcpm <- y.corrected[top_var, ]
dim(highly_variable_lcpm)


## Get some nicer colours

mypalette <- brewer.pal(7,"RdYlBu")
#morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
#col.cell <- revalue(meta_table$Treatment, c("DOI" = "purple", "Saline" = "orange"))
#col.cell <- as.matrix(cbind(col.cell,revalue(meta_table$Timepoint, c("1week" = "blue", "1month" = "red"))))

pheatmap(highly_variable_lcpm,
         #color = rev(mypalette),
         cluster_rows = T, 
         cluster_cols = T,
         annotation_col = meta_table[, c(2:6)],
         annotation_colors = my_colour,
         gaps_col =  3,
         color=colorRampPalette(c("#4575b4", "white", "#d73027"))(50),
         main="Top 50 variable genes across conditions",
         border_color = F)


####

# pseudoCounts <- log2(y$counts+1)
# head(pseudoCounts)
# 
# boxplot(logcounts, col="gray")
# plotMDS(logcounts)
# 
# sampleDists <- as.matrix(dist(t(pseudoCounts)))
# sampleDists
# 
# cimColor <- colorRampPalette(rev(brewer.pal(9, "Blues")))(16)
# pheatmap(sampleDists, col=cimColor, symkey=FALSE)


# Plot the heatmap
# heatmap.2(highly_variable_lcpm,
#           col=rev(morecols(50)),
#           trace="none", 
#           main="Top variable genes across conditions",
#           ColSideColors=col.cell,
#           scale="row",
#           margins=c(15,15),
#           cexRow=1.5)

# 
# go <- goana(qlf, species="Mm")
#  topGO(go, sort="up")


# 
# Differential isoform expression analysis
# Transcript counts were first obtained using Salmon software v.0.11.2 [18]. Differential isoform usage was then analyzed using edgeR [66, 67], 
# considering an isoform as differentially expressed when the adjusted p-value in the comparison between samples was below 0.01.
# 
# Differential splicing analysis
# Differential splicing was assessed using the vast-tools software package v2.0.1 (https://github.com/vastgroup/vast-tools), 
# normalizing the distribution of each event by the frequency of that event in the mouse transcriptome (Mmu10). Differentially spliced events in males vs females (at E11 and E12) 
# were determined by calculating the difference in their average inclusion levels (ΔPSI), considering only those ΔPSI that were higher than 0.1. Events were classified as: 
#   exon skipping (ES), alternative 3’splice site (3’ss), alternative 5′ splice site (5’ss), intron retention (IR) and microexon (MIC). Intron retention events, which correspond to U12 events, 
# were then identified comparing our results with mouse U12 events obtained from U12db (http://genome.crg.es/cgi-bin/u12db/u12db.cgi) [39].

# Group <- factor(paste(meta_table$Treatment,meta_table$Timepoint,sep="."))
# meta_table = cbind(meta_table,Group=Group)

# https://support.bioconductor.org/p/45640/


# time <- factor(time, levels=c("0h","6h","12h","18h","24h","36h","48h"))
# #design <- model.matrix(~Timepoint+Treatment:Timepoint, meta_table)
# design <- model.matrix(~Group, meta_table)
# 
# y <- estimateGLMCommonDisp(y, design)
# y <- estimateGLMTrendedDisp(y, design)
# y <- estimateGLMTagwiseDisp(y, design)
# fit <- glmQLFit(y, design)
# 
# qlf <- glmQLFTest(fit, coef=1)
# topTags(qlf)
# 
# plotQLDisp(fit)
# #To test for treatment effect at time 0:
#   
# lrt0 <- glmLRT(y,fit,coef=8)
# 
# topTags(lrt0)
# 
# # To test for treatment effect at time 1:
#   
# lrt1 <- glmLRT(y,fit,coef=9)
# topTags(lrt1)

# ##############################
#  my.contrasts <- makeContrasts(BvsA=B-A, CvsB=C-B, CvsA=C-A, levels=design)
# qlf.BvsA <- glmQLFTest(fit, contrast=my.contrasts[,"BvsA"])
# topTags(qlf.BvsA)
# qlf.CvsB <- glmQLFTest(fit, contrast=my.contrasts[,"CvsB"])
# topTags(qlf.CvsB)
# qlf.CvsA <- glmQLFTest(fit, contrast=my.contrasts[,"CvsA"])
# topTags(qlf.CvsA)

design <- model.matrix(~meta_table$Sex+meta_table$Treatment)
design <- model.matrix(~meta_table$Group)
colnames(design) <- levels(meta_table$Group)
rownames(design) <- colnames(y)

## https://www.biostars.org/p/400408/

y <- estimateDisp(y, design)


condition = "DOI vs Saline - 48 hours"
condition = "DOI vs Saline - 1 week"
condition = "DOI vs Saline - 1 month"
et <- exactTest(y)
topTags(et)

data = et$table
colnames(data) <- c("logFC", "logCPM", "PValue")
data$Genes = rownames(data)

p_fc = 0
p_pval = 0.05

data <- data %>% 
  mutate(
    Expression = case_when(logFC >= p_fc & PValue <= p_pval ~ "Up-regulated",
                           logFC <= -p_fc & PValue <= p_pval ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )
head(data) 
table(data$Expression)




##############
# 
# fit <- glmQLFit(y, design)
# 
# #### edgeR offers a rigorous statistical test for thresholded hypotheses under the GLM framework
# tr <- glmTreat(fit, coef=2, lfc=1)
# topTags(tr)
# 
# data = tr$table[, c(1, 3, 4)]
# colnames(data) <- c("logFC", "logCPM", "PValue")
# data$Genes = rownames(data)
# 
# p_fc = 0
# p_pval = 0.05
# 
# data <- data %>% 
#   mutate(
#     Expression = case_when(logFC >= p_fc & PValue <= p_pval ~ "Up-regulated",
#                            logFC <= -p_fc & PValue <= p_pval ~ "Down-regulated",
#                            TRUE ~ "Unchanged")
#   )
# head(data) 
# table(data$Expression)
# 
# 
# ##
# qlf <- glmQLFTest(fit)
# topTags(qlf)
# 
# 
# # my.contrasts <- makeContrasts(
# #   DOIvsSaline.48h = DOI.48hrs-Ctrl.48hrs,
# #   DOIvsSaline.1W = (DOI.1week-DOI.48hrs)-(Ctrl.1week-Ctrl.48hrs),
# #   DOIvsSaline.1M = (DOI.1month-DOI.1week-DOI.48hrs)-(Ctrl.1month-Ctrl.1week-Ctrl.48hrs),
# #   levels=design)
# 
# my.contrasts <- makeContrasts(
#   DOIvsSaline.48h = DOI.48hrs-Ctrl.48hrs,
#   levels=design)
# 
# 
# my.contrasts <- makeContrasts(
#   DOIvsSaline.1W = DOI.1week-Ctrl.1week,
#   levels=design)
# 
# 
# my.contrasts <- makeContrasts(
#   DOIvsSaline.1M = DOI.1month-Ctrl.1month,
#   levels=design)
# 
# 
# plotMDS(y)
# # 
# plotBCV(y)
# # 
# # # ###
# condition = "DOI vs Saline - 48 hours"
# qlf <- glmQLFTest(fit, contrast=my.contrasts[,"DOIvsSaline.48h"])
# topTags(qlf)
# 
# 
# condition = "DOI vs Saline - 1 week"
# #qlf <- glmQLFTest(fit, contrast=my.contrasts[,"DOIvsSaline.1W"])
# qlf <- glmQLFTest(fit)
# topTags(qlf)
# 
# # ###
# # qlf <- glmQLFTest(fit, contrast=my.contrasts[,"DOIvsSaline.1W"])
# # topTags(qlf)
# 
# # ###
# condition = "DOI vs Saline - 1 month"
# qlf <- glmQLFTest(fit, contrast=my.contrasts[,"DOIvsSaline.1M"])
# topTags(qlf)
# # 
# # qlf <- glmQLFTest(fit, contrast=my.contrasts[,"DOIvsSaline.1M.2"])
# # topTags(qlf)
# 
# #summary(decideTests(qlf))
# 
# #plotMD(qlf)
# #abline(h=c(-1, 1), col="blue")
# 
# # volcanoData <- cbind(qlf$table$logFC, -log10(qlf$table$PValue))
# # colnames(volcanoData) <- c("logFC", "negLogPval")
# # head(volcanoData)
# # 
# # plot(volcanoData, pch=19)
# 
# 
# data = qlf$table[, c(1, 2, 3, 4)]
# colnames(data) <- c("logFC", "logCPM", "FDR", "PValue")
# data$Genes = rownames(data)
# 
# # data <- data %>% 
# #   mutate(
# #     Significance = case_when(
# #       abs(logFC) >= log(2) & PValue <= 0.05 & PValue > 0.01 ~ "PValue 0.05", 
# #       abs(logFC) >= log(2) & PValue <= 0.01 & PValue > 0.001 ~ "PValue 0.01",
# #       abs(logFC) >= log(2) & PValue <= 0.001 ~ "PValue 0.001", 
# #       TRUE ~ "Unchanged")
# #   )
# # head(data)
# # 
# # 
# # p3 <- ggplot(data, aes(logFC, -log10(PValue))) +
# #   geom_point(aes(color = Significance), size = 1) +
# #   xlab(expression("log"[2]*"FC")) + 
# #   ylab(expression("-log"[10]*"PValue")) +
# #   scale_color_viridis_d() +
# #   guides(colour = guide_legend(override.aes = list(size=1.5))) +
# #   theme_classic(base_size = 17)
# # p3
# 
# p_fc = 0
# p_pval = 0.05
# 
# 
# 
# data <- data %>% 
#   mutate(
#     Expression = case_when(logFC >= p_fc & PValue <= p_pval ~ "Up-regulated",
#                            logFC <= -p_fc & PValue <= p_pval ~ "Down-regulated",
#                            TRUE ~ "Unchanged")
#   )
# head(data) 
# table(data$Expression)
###########################################################################################

#data = data %>% arrange(PValue, desc(abs(logFC)))


write.table(data, "results/DEG/Episode III/Timepoint/genes/PV/DEG_DOIvsSaline_48hrs_only_EdgeR.txt", sep = "\t", quote = F, row.names = F)
write.table(data, "results/DEG/Episode III/Timepoint/genes/1week/PV/pipe0/DEG_DOIvsSaline_1week_only_EdgeR.txt", sep = "\t", quote = F, row.names = F)
write.table(data, "results/DEG/Episode III/Timepoint/genes/PV/DEG_DOIvsSaline_1month_only_EdgeR.txt", sep = "\t", quote = F, row.names = F)


write.table(data, "results/DEG/Episode III/Timepoint/genes/Rbp4/DEG_DOIvsSaline_48hrs_only_EdgeR.txt", sep = "\t", quote = F, row.names = F)
write.table(data, "results/DEG/Episode III/Timepoint/genes/Rbp4/DEG_DOIvsSaline_1week_only_EdgeR.txt", sep = "\t", quote = F, row.names = F)
write.table(data, "results/DEG/Episode III/Timepoint/genes/Rbp4/DEG_DOIvsSaline_1month_only_EdgeR.txt", sep = "\t", quote = F, row.names = F)

write.table(data, "results/DEG/Episode III/Timepoint/genes/Cux2/DEG_DOIvsSaline_48hrs_only_EdgeR.txt", sep = "\t", quote = F, row.names = F)
write.table(data, "results/DEG/Episode III/Timepoint/genes/Cux2/DEG_DOIvsSaline_1week_only_EdgeR.txt", sep = "\t", quote = F, row.names = F)

#write.table(data, "results/DEG/Cux2/Timepoint/DEG_DOIvsSaline_48hrs_only_EdgeR.txt", sep = "\t", quote = F, row.names = F)
#write.table(data, "results/DEG/Cux2/Timepoint/DEG_DOIvsSaline_1week_only_EdgeR.txt", sep = "\t", quote = F, row.names = F)
#write.table(data, "results/Rbp4/Timepoint/DEG_DOIvsSaline_1month_only_EdgeR.txt", sep = "\t", quote = F, row.names = F)


top <- 15
top_genes <- bind_rows(
  data %>% 
    filter(Expression == 'Up-regulated') %>% 
    arrange(PValue, desc(abs(logFC))) %>% 
    head(top),
  data %>% 
    filter(Expression == 'Down-regulated') %>% 
    arrange(PValue, desc(abs(logFC))) %>% 
    head(top)
)
top_genes


p2 <- ggplot(data, aes(logFC, -log(PValue,10))) +
  geom_point(aes(color = Expression), size = 2, alpha = 0.4,) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"PValue")) +
  scale_color_manual(values = c("#7fbf7b", "#f7f7f7", "#af8dc3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  labs(title = paste0("Gene expression changes in ", condition)) +
  theme_classic(base_size = 22) +
  # Change line type and color
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "black") +
  geom_vline(xintercept=-p_fc, linetype="dashed", color = "black") +
  geom_vline(xintercept=p_fc, linetype="dashed", color = "black")+
  theme(plot.title = element_text(hjust = 0.3)) +
  #scale_y_continuous(limits = c(0, 7), breaks = c(0, 2, 4, 6)) +
  #scale_x_continuous(limits = c(-15, 15), breaks = c(-15, -10, -5, 0, 5, 10))
  scale_x_continuous(limits = c(-5, 5), breaks = c(-7, -5, -2, 0, 2, 5, 7))
p2


p2 <- p2 +   
  geom_text_repel(
  data          = subset(top_genes, Expression == "Up-regulated"),
  aes(label = top_genes$Genes[which(top_genes$Expression == "Up-regulated")]),
  segment.size  = 0.4,
  direction     = "y",
  size = 6,
  segment.color = "grey50",
  hjust         = 0
  ) +
  geom_text_repel(
    data          = subset(top_genes, Expression == "Down-regulated"),
    aes(label = top_genes$Genes[which(top_genes$Expression == "Down-regulated")]),
    segment.size  = 0.4,
    direction     = "y",
    size = 6,
    segment.color = "grey50",
    hjust         = 1
  )
p2

# pdf(paste0("results/volcano_genelist_Michael_samples_", condition, "_EdgeR_low_FC.pdf"), height = 9.5, width = 10.5)
# p2
# dev.off()


# 
# p2 <-  p2 +
#   geom_label_repel(data = top_genes,
#                    mapping = aes(logFC, -log(PValue,10), label = Genes),
#                    #direction = "y",
#                    nudge_y = 2,
#                    size = 6)
# p2
# 
# 
# 
# p3 <-  p3 +
#   geom_label_repel(data = top_genes,
#                    mapping = aes(logFC, -log(PValue,10), label = Genes),
#                    position = dodge,
#                    size = 4)
# p3


# grid::grid.newpage()
# gridExtra::grid.table(table(data$Expression, data$Significance))


# top <- 15
# top_genes <- bind_rows(
#   data %>%
#     filter(Expression == 'Up-regulated') %>%
#     arrange(logFC, PValue) %>%
#     head(top),
#   data %>%
#     filter(Expression == 'Down-regulated') %>%
#     arrange(desc(logFC), PValue) %>%
#     head(top)
# )
# top_genes

#top_genes = head(data[order(data$logFC, decreasing = T),], 15)


####
top_genes

# sel_col = meta_table$Sample.ID[meta_table$Timepoint == "1week"]
# sel_col = meta_table$Sample.ID[meta_table$Timepoint == "1month"]
# sel_col = meta_table$Sample.ID[meta_table$Timepoint == "48hrs"]

#colnames(logcounts)

#highly_variable_lcpm <- logcounts[top_genes$Genes[which(top_genes$Expression == 'Down-regulated')], colnames(logcounts) %in% sel_col]
# highly_variable_lcpm <- ZlogCPM[top_genes$Genes[which(top_genes$Expression == 'Down-regulated')],]
# title_txt = paste0("Top FC", nrow(highly_variable_lcpm),  " Down-regulated DEG - ", condition)


#highly_variable_lcpm <- ZlogCPM[top_genes$Genes[which(top_genes$Expression == 'Up-regulated')],]

highly_variable_lcpm <- ZlogCPM[top_genes$Genes,]
title_txt = paste0("Top ", nrow(highly_variable_lcpm)/2,  " Up and Down-regulated DEG - \n", condition, " (z-score expression)")

# highly_variable_lcpm <- ZlogCPM[data$Genes[which(data$Expression == 'Up-regulated')],]
# title_txt = paste0("Top FC ", nrow(highly_variable_lcpm),  " Up-regulated DEG - ", condition)
# title_txt = paste0("Top ", nrow(highly_variable_lcpm),  " Up-regulated DEG - ", condition)


#highly_variable_lcpm <- counts_vst[top_genes$Genes[which(top_genes$Expression == 'Up-regulated')],]
#highly_variable_lcpm = highly_variable_lcpm[,sel_col]
#dim(highly_variable_lcpm)


# ####
# y <- calcNormFactors(y)
# 
# logcounts <- cpm(y, log=TRUE)
# head(logcounts)
# 
# 
# title_txt = paste0("Collagen formation genes acroos timepoints")
# genes = c("Bmp1", "Col11a1", "Col11a2", "Col12a1", "Col15a1", "Col16a1", "Col19a1", "Col1a2", "Col23a1", "Col27a1", "Col4a1", "Col4a2", "Col4a6", "Col5a1",
#   "Col5a2","Col5a3","Col6a1","Col6a2","Col9a3","P3h3","Plod3","Pxdn")
# sel_col = meta_table$Sample.ID[meta_table$Timepoint %in% c("48hrs", "1week", "1month")]
# highly_variable_lcpm = ZlogCPM[genes, sel_col]
# 
# highly_variable_lcpm = logcounts[genes, sel_col]
# nrow(highly_variable_lcpm) == length(genes)
# 
# title_txt = paste0("Selected genes ", condition)
# genes = c("Rab5a", "Syp", "Eif5", "Scn1b", "Bdnf", "Syngap1", "Nf1", "Dlg4", "Grin2a", "Kmt2d", "Grin1")
# highly_variable_lcpm = ZlogCPM[genes, sel_col]
# 
# highly_variable_lcpm = logcounts[genes, sel_col]
# nrow(highly_variable_lcpm) == length(genes)
# 
# 
# mypalette <- brewer.pal(7,"RdYlBu")

# my_colour = list(
#   Treatment = c(DOI = "#5977ff", Ctrl = "#f74747"),
#   Timepoint = c(`1week` = "#e89829", `1month` = "#82ed82")
# )

breaksList = seq(-2, 2, by = 0.1)

my_colour = list(
  Treatment = c(DOI = "#5977ff", Ctrl = "#f74747"),
  CellType = c(Rbp4 = "#f6e8c3", PV = "#9e82ed", Cux2ERT = "#5ab4ac"),
  Sex = c(`F` = "#d73027", `M` = "#4575b4"),
  Batch = c(`B1` = "#1b7837", `B2` = "#af8dc3"),
  Timepoint = c(`1week` = "#e89829", `1month` = "#82ed82", `48hrs` = "brown")
)

#highly_variable_lcpm = highly_variable_lcpm[, meta_table$Sample.ID]


pheatmap(highly_variable_lcpm,
       cluster_rows = F, 
       cluster_cols = F,
       annotation_col = meta_table[, c(2:6)],
       annotation_colors = my_colour,
       color=colorRampPalette(c("#4575b4", "white", "#d73027"))(length(breaksList)),
       breaks = breaksList,
       gaps_col = 5,
       #gaps_col =  c(5, 10, 15, 20, 25),
       main=title_txt,
       border_color = F)

library('org.Mm.eg.db')
library(ReactomePA)
library(clusterProfiler)
library(enrichplot)
# 
# top <- 10
# top_genes <- bind_rows(
#   data %>% 
#     filter(Expression == 'Up-regulated') %>% 
#     arrange(PValue, desc(abs(logFC))) %>% 
#     head(top),
#   data %>% 
#     filter(Expression == 'Down-regulated') %>% 
#     arrange(PValue, desc(abs(logFC))) %>% 
#     head(top)
# )
# top_genes


# gs.ephi <- lapply(colnames(data.m.prob.o), function(x) {
#   print(x)
#   
#   clu <- ephi.markers.extra[ephi.markers.extra$cluster == x,]
#   clu <- clu[clu$p_val_adj < 0.05 & clu$avg_logFC > 0,]
#   print(head(clu))
#   
#   symbols = clu$gene
#   IDs = mapIds(org.Mm.eg.db, symbols, 'ENTREZID', 'SYMBOL')
#   
#   return(IDs)
#   
# })


top <- 500
top_genes <- bind_rows(
  data %>%
    filter(Expression == 'Up-regulated') %>%
    arrange(PValue, desc(abs(logFC))) %>%
    head(top),
  data %>%
    filter(Expression == 'Down-regulated') %>%
    arrange(PValue, desc(abs(logFC))) %>%
    head(top)
)
top_genes
# 

data = read.delim("results/DEG/Episode II/Cux2/Timepoint/genes/Cux2_table_all_genes_mean_SD.csv", sep = ",")

timep = "48hrs"
timep = "1week"
timep = "1month"
status = "Up"
symbols = data$Genes[data$Expression == 'Up-regulated']
#symbols = top_genes$Genes[top_genes$Expression == 'Up-regulated']

#timep = "48hrs"
status = "Down"
symbols = data$Genes[data$Expression == 'Down-regulated' & data$id == "Cux2_1week"]
#symbols = top_genes$Genes[top_genes$Expression == 'Down-regulated']

gene <- mapIds(org.Mm.eg.db, symbols, 'ENTREZID', 'SYMBOL')
gene = gene[!is.na(gene)]
#   
yy = enrichPathway(gene, pvalueCutoff=0.05, organism = "mouse", readable = TRUE)
head(summary(yy))
clusterProfiler::dotplot(yy, showCategory=10)
write.table(as.data.frame(yy) ,paste0("results/DEG/Episode III/Timepoint/genes/Cux2/DEG_DOIvsSaline_", timep, "_only_EdgeR_Reactome_analysis_", status, "_genes.txt"), sep = "\t", quote = F, row.names = F)
#write.table(as.data.frame(yy) ,"results/Rbp4/Timepoint/DEG_DOIvsSaline_1month_only_EdgeR_Reactome_analysis_Down_genes.txt", sep = "\t", quote = F, row.names = F)

# 
ego <- enrichGO(gene          = gene,
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)

pdf("figures/GO_CUX2_1week_down.pdf", height = 10, width = 8)
barplot(ego, showCategory=20) 
dev.off()
write.table(as.data.frame(ego) ,paste0("results/DEG/Episode III/Timepoint/genes/Cux2/DEG_DOIvsSaline_", timep, "_only_EdgeR_GO_analysis_", status, "_genes.txt"), sep = "\t", quote = F, row.names = F)
#write.table(as.data.frame(ego) ,"results/PV/Timepoint/DEG_DOIvsSaline_1month_only_EdgeR_GO_analysis_Down_genes.txt", sep = "\t", quote = F, row.names = F)
 
# http://datatripping:8787/graphics/plot_zoom_png?width=687&height=340
# y <- gsePathway(gene, 
#                 organism = "mouse",
#                 pvalueCutoff = 0.2,
#                 pAdjustMethod = "BH", 
#                 verbose = FALSE)
# head(y)
# 
# 
# 




