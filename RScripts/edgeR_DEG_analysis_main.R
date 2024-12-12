### https://github.com/bli25/RSEM_tutorial
rm(list = ls())

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
library(ggpubr)
library(sva)
library(ggpmisc)
library(rstatix)
library(ComplexUpset)
library(ComplexHeatmap)
library(MuSiC)
library(biomaRt)

#set seed for reproductive results
set.seed(123) 

# Load scripts
source("src/edgeR_DEG_analysis_IO.R")
source("src/edgeR_DEG_analysis_visualization.R")
source("src/edgeR_DEG_analysis_processing.R")

figure_colors = list (
  Treatment = c(DOI = "#e54291", Ctrl = "black", Psilocybin="#1f87c5"),
  CellType = c(Rbp4 = "#0f7733", PV = "#a94498", Cux2ERT = "#45aa9a"),
  Sex = c('F' = "#d9d9d9", 'M' = "#525252"),
  Batch = c('B1' = "#dfc27d", 'B2' = "#c7eae5", 'B3' = "#bac4e3", 'B4' = "#b3be62"),
  Timepoint = c("48hrs" = "#fdae6b", "1week" = "#f16913", "1month" = "#8c2d04"),
  Status = c("Down-regulated" = "orange", "Unchanged" = "#f7f7f7", "Up-regulated" = "blue")
)
  
################################################################################
#                          all cell types PCA visualization

ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
bm <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'gene_biotype', 'ensembl_transcript_id'), mart = ensembl)
#save(bm, file = "data/rdata/bm.rds")
head(bm)

meta_table = get_metatable(all = T)
plot_table(meta_table)

y = get_DGEList_counts(meta_table)

y = filter_genes(y)

logcounts <- cpm(y, normalized.lib.sizes = TRUE, log=TRUE)
head(logcounts)
p1 = plot_PCA(logcounts, meta_table, "Batch", "Treatment", TRUE)
p1 = plot_PCA(logcounts, meta_table, "Batch", "CellType", TRUE)
p1 = p1 + scale_color_manual(values = figure_colors$Batch)

export_svg(p1, dir = "figures/drafts/", file_name = "PCA_raw_counts", h = 3, w = 5)

dds <- DESeqDataSetFromMatrix(countData=round(y$counts,0),colData=meta_table,design = ~ Batch + CellType)
rlog <- assay(vst(dds,blind=T))
mod <- model.matrix(~ CellType, meta_table)
cbat <- ComBat(rlog, batch=meta_table$Batch, mod=mod)

pca_combat = prcomp(t(cbat))

p1 = autoplot(pca_combat,
              data = meta_table,
              colour="Treatment",
              shape="Sex",
              size=5) +
  theme_classic2(base_size = 17) +
  geom_text_repel(aes(x=PC1, y=PC2, label=Sample.ID), box.padding = 0.8) +
  scale_color_manual(values = figure_colors$Treatment)
p1
p1 = autoplot(pca_combat,
              data = meta_table,
              colour="CellType",
              shape="Treatment",
              size=5) +
  theme_classic2(base_size = 17) +
  geom_text_repel(aes(x=PC1, y=PC2, label=Sample.ID), box.padding = 0.8) +
  scale_color_manual(values = figure_colors$CellType)
p1

export_svg(p1, dir = "figures/drafts/", file_name = "PCA_combat", h = 3, w = 5)

p1 = autoplot(pca_combat,
              data = meta_table,
              colour="Batch",
              shape="CellType",
              size=5) +
  theme_classic2(base_size = 17) +
  geom_text_repel(aes(x=PC1, y=PC2, label=Sample.ID), box.padding = 0.8) +
  scale_color_manual(values = figure_colors$Batch)
p1

export_svg(p1, dir = "figures/drafts/", file_name = "PCA_batch_combat", h = 3, w = 5)

#########################################################################################
###                           Bowtie2 + RSEM + EdgeR DEG Analysis
#########################################################################################

##### Pipe0 EdgeR ####
pipeline = "pipe0"
cell_type = c("PV", "Cux2ERT", "Rbp4") # [PV Cux2ERT Rbp4]
timepoint = c("48hrs", "1week", "1month") # [48hrs 1week 1month]
treat = "DOI" # [DOI Psilocybin Ctrl]
title_condition = paste(treat,"vs Saline -", cell_type, timepoint, pipeline, sep = " ")
title_condition

# for getting all Ctrl samples, and order by PV, L2/3, L5
meta_table <- meta_table[meta_table$Treatment == "Ctrl",]
my_order <- c("PV", "Cux2ERT", "Rbp4")
meta_table = meta_table %>%
  arrange(match(CellType, my_order))

y = get_DGEList_counts(meta_table, treat = treat)
dim(y$counts)
all(colnames(y$counts %in% meta_table$Sample.ID))
y = filter_genes(y, large_n = 7, min_count = 50)

logcounts <- cpm(y, normalized.lib.sizes = TRUE, log=TRUE)
head(logcounts)
plot_boxplot_topN_genes(logcounts, 30, ylimits=c(-5, 12), "logCPM raw_data")

ZlogCPM <- t(scale(t(logcounts)))
head(ZlogCPM)
plot_boxplot_topN_genes(ZlogCPM, 30, ylimits=c(-4, 4), "Z-score data")

##
plot_PCA(logcounts, meta_table, "Treatment", "Sex", TRUE)
counts_vst = vst(round(y$counts))
plot_PCA(counts_vst, meta_table, "Treatment", "Sex", TRUE)

top_var <- get_topn_var_genes(logcounts, 500)
plot_PCA(counts_vst[ top_var, ], meta_table, "Treatment", "Sex", TRUE)

plot_heatmap(ZlogCPM, meta_table, top_var = top_var, title="Top 500 variable genes across conditions")

## add batch correction TRUE, FALSE
data = get_DEG(y, batch = FALSE)

data <- data %>% 
  mutate(
    Expression = case_when(logFC >= 0 & FDR <= 0.1 ~ "Up-regulated",
                           logFC <= -0 & FDR <= 0.1 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )
colnames(data)[1] <- "Genes"

head(data, 20) 
table(data$Expression)
top_genes = head(data$Genes, 20)

top_expr = ZlogCPM[rownames(ZlogCPM) %in% top_genes, meta_table$Sample.ID]

p1 = pheatmap(top_expr,
              cluster_rows = T, 
              cluster_cols = F,
              show_rownames = T,
              color=colorRampPalette(c("#2e4a95", "black", "#f5ec03"))(50),
              main=paste0(nrow(top_expr), " top expressed genes"
              )
)

save_file(data, 
          dir = paste0("results/DEG/Episode VI/", pipeline, "/", cell_type, "/", timepoint, "/"),
          file_name = paste0("DEG_", title_condition, ".txt"),
          row = FALSE,
          sep = "\t")

#########################################################################################
###                           Bowtie2 + RSEM + DESeq2 DEG Analysis
#########################################################################################

##### Pipe6 DESeq2 ####
pipeline = "pipe6"
cell_type = c("PV","Cux2ERT","Rbp4") # [PV Cux2ERT Rbp4]
timepoint = c("48hrs", "1week", "1month") # [48hrs 1week 1month]
treat = c("DOI") # [DOI Psilocybin Ctrl]
title_condition = paste(treat,"vs Saline -", cell_type, timepoint, pipeline, sep = " ")
title_condition

meta_table <- get_metatable(cell_type, timepoint, treat)
#meta_table <- meta_table[1:nrow(meta_table)-1,] # PV 1 month
meta_table

## for plotting all Ctrl samples
meta_table <- meta_table[meta_table$Treatment == "Ctrl",]
my_order <- c("PV", "Cux2ERT", "Rbp4")
meta_table = meta_table %>%
  arrange(match(CellType, my_order))

files = get_genes_results_files(meta_table, treat)
files

txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
head(txi.rsem$counts)
dim(txi.rsem$counts)

zero_length_and_unexpressed = (apply(txi.rsem$counts, 1, max) == 0) |
  (apply(txi.rsem$length, 1, min) == 0)

txi.rsem$length = txi.rsem$length[!zero_length_and_unexpressed,]
txi.rsem$abundance = txi.rsem$abundance[!zero_length_and_unexpressed,]
txi.rsem$counts = txi.rsem$counts[!zero_length_and_unexpressed,]
dim(txi.rsem$counts)
head(txi.rsem$counts)

dds <- DESeqDataSetFromTximport(txi.rsem, meta_table, ~Treatment)
all(rownames(meta_table) == colnames(txi.rsem))

smallestGroupSize <- 8
keep <- rowSums(counts(dds) >= 50) >= smallestGroupSize
table(keep)
out <- bm[match(names(keep), bm$external_gene_name),]
keep = keep & (out$gene_biotype=="protein_coding" & !is.na(out$gene_biotype))
table(keep)
dds <- dds[keep,]

dds$condition <- dds$Treatment
dds$condition <- relevel(dds$condition, ref = "Ctrl")
dds$condition <- droplevels(dds$condition)

dds <- DESeq(dds)
res <- results(dds)
res
summary(res)

# 1 month PV
#res.ihw <- results(dds, filterFun=ihw)
#summary(res.ihw)
#res <- res.ihw
resOrdered <- res[order(res$padj),]

resOrdered_df = as.data.frame(resOrdered)
resOrdered_df$Genes = rownames(resOrdered_df)

resOrdered_df <- resOrdered_df %>% 
  mutate(
    Expression = case_when(log2FoldChange >= 0 & padj <= 0.1 ~ "Up-regulated",
                           log2FoldChange <= -0 & padj <= 0.1 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )
head(resOrdered_df,20)
save_file(resOrdered_df, 
          dir = paste0("results/DEG/Episode VI/", pipeline, "/", cell_type, "/", timepoint, "/"),
          file_name = paste0("DEG_", title_condition, ".txt"),
          row = FALSE,
          sep = "\t")

topVarGenes <- head(rownames(resOrdered), 500)
plotCounts(dds, gene=which.min(res$pvalue), intgroup="condition")
vsd <- vst(dds, blind = FALSE)
Z <- t(scale(t(assay(vsd))))
mat <- Z[topVarGenes, ]

breaksList = c(seq(-3, 3, by = 0.1))
breaksList = unique(breaksList)
my_ramp = colorRampPalette(c("#2e4a95", "black", "#f5ec03"))(length(breaksList))

p1 = pheatmap(mat,
              cluster_rows = T, 
              cluster_cols = F,
              show_rownames = F,
              show_colnames = F,
              border_color = F,
              annotation_col = meta_table[, c(2:6)],
              breaks = breaksList,
              color=my_ramp,
              annotation_colors = my_colour,
              main=paste0(nrow(mat), " top differentially expressed genes"
              )
)
p1

#########################################################################################
##### supp fig 2 scatter plot comparing pipelines #####

cell_type = c("PV") # [PV Cux2ERT Rbp4]
timepoint = c("48hrs") # [48hrs 1week 1month]
treat = "DOI" # [DOI Psilocybin Ctrl]

res_pipe6 <- read_tsv(paste0("results/DEG/Episode VI/pipe6/", cell_type, "/", timepoint, "/DEG_", treat, " vs Saline - ", cell_type, " ", timepoint, " pipe6.txt"), show_col_types = FALSE)
#head(res_pipe6)
#nrow(res_pipe6)
table(res_pipe6$padj < .05)
table(res_pipe6$padj < .10)

res_pipe0 <- read_tsv(paste0("results/DEG/Episode VI/pipe0/", cell_type, "/", timepoint, "/DEG_", treat, " vs Saline - ", cell_type, " ", timepoint, " pipe0.txt"), show_col_types = FALSE)
#head(res_pipe0)
#nrow(res_pipe0)
table(res_pipe0$FDR < .05)
table(res_pipe6$padj < .10)

aux1 <- res_pipe6 %>% 
  dplyr::select(Genes, log2FoldChange) %>%
  na.omit()
#head(aux)
#nrow(aux)

aux2 <- res_pipe0 %>% 
  dplyr::select(Genes, logFC) %>%
  na.omit()
#head(aux2)
#nrow(aux2)

aux1$ranks <- rownames(aux1)
aux2$ranks <- rownames(aux2)

aux3 <- merge(aux1, aux2, by = "Genes", all = TRUE) %>% na.omit()
#head(aux3)
aux3[which(aux3$Genes=="Htr1a"),]

p <- ggplot(aux3, aes(x=log2FoldChange, y=logFC)) + 
  geom_point() +
  theme(panel.background = element_blank()) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "R2"))) +
  xlim(-2.5,2.5) +
  ylim(-2.5,2.5) +
  labs(x = "DESeq2 log2FC") +
  labs(y = "EdgeR log2FC") +
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  ggtitle(paste(cell_type, timepoint, treat, sep = " "))
p
export_svg(p, dir = "results/SuppFig2NoBG/", file_name = paste(cell_type, timepoint, treat, sep = " "), h = 6, w = 6)

p <- ggplot(head(aux3,500), aes(x=ranks.x, y=ranks.y)) + 
  geom_point() +
  theme(panel.background = element_blank(), axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "R2"))) +
  labs(x = "DESeq2 ranks") +
  labs(y = "EdgeR ranks") +
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  ggtitle(paste(cell_type, timepoint, treat, sep = " "))
p
export_svg(p, dir = "results/SuppFig2NoBG/", file_name = paste(cell_type, timepoint, treat, "ranks", sep = " "), h = 6, w = 6)


########################################################################################## Marker Genes Expression

## Hrvatin gene markers
Hrvatin_cell_markers = get_Hrvatin_cell_marker_genes()
Hrvatin_cell_markers_filt = Hrvatin_cell_markers[Hrvatin_cell_markers$Cell_type %in% c('Int_Pv', 'ExcL23', 'ExcL5_1', 'ExcL5_2', 'ExcL5_3'),]
order_vector <- c('Int_Pv', 'ExcL23', 'ExcL5_1', 'ExcL5_2', 'ExcL5_3')

# Reorder the rows based on the group and the order vector
Hrvatin_cell_markers_filt <- Hrvatin_cell_markers_filt %>%
  group_by(Cell_type) %>%
  mutate(order = factor(Cell_type, levels = order_vector)) %>%
  arrange(order) %>%
  ungroup() %>%
  select(-order)

Hrvatin <- as.data.frame(Hrvatin_cell_markers_filt[which(Hrvatin_cell_markers_filt$Gene_name %in% resOrdered_df$Genes),])
dim(Hrvatin)
rownames(Hrvatin) <- Hrvatin$Gene_name

vsd <- vst(dds, blind = FALSE)
Z <- t(scale(t(assay(vsd))))
mat <- Z[topVarGenes, ]

breaksList = c(seq(-3, 3, by = 0.1))
breaksList = unique(breaksList)
my_ramp = colorRampPalette(c("#2e4a95", "black", "#f5ec03"))(length(breaksList))

my_colour = list(
  Treatment = c(DOI = "#e54291", Ctrl = "black", Psilocybin="#1f87c5"),
  CellType = c(Rbp4 = "#0f7733", PV = "#a94498", Cux2ERT = "#45aa9a"),
  Sex = c('F' = "#d9d9d9", 'M' = "#525252"),
  Batch = c('B1' = "#dfc27d", 'B2' = "#c7eae5", 'B3' = "#bac4e3", 'B4' = "#b3be62"),
  Timepoint = c("48hrs" = "#fdae6b", "1week" = "#f16913", "1month" = "#8c2d04"),
  Cell_type = c('Int_Pv'= "#a94498", 'ExcL23' = "lightblue", 'ExcL5_1'= "#e6f5d0", 'ExcL5_2'= "#a1d76a", 'ExcL5_3' = "#4d9221")
)

p1 = pheatmap(mat,
              cluster_rows = T, 
              cluster_cols = F,
              show_rownames = T,
              show_colnames = F,
              border_color = F,
              annotation_col = meta_table[, c(2:6)],
              annotation_row = Hrvatin[, 3, drop=FALSE],
              breaks = breaksList,
              color=my_ramp,
              annotation_colors = my_colour,
              main=paste0(nrow(mat), " top differentially expressed genes"
              )
)
p1
export_pdf(p1, dir = "figures/drafts/", file_name = "Hrvatin_cell_markers_heatmap", h = 9, w = 10)

############################################################################## Plot FPKM TPM gene

cell_type = c("Cux2ERT") # [PV Rbp4 Cux2ERT]
timepoint = c("1month") # [48hrs 1week 1month]
treat = "DOI"

meta_table = get_metatable(cell_type, timepoint, all = T)
meta_table = meta_table[meta_table$Treatment == "Ctrl",]
plot_table(meta_table)

typev = "FPKM"
fpkm_matrix = get_counts_RNA(meta_table = meta_table, treat = "DOI", type = typev) # [DOI or Psilocybin samples ]
head(fpkm_matrix)

rt_genes = get_Receptors_T_genes()
gene_list = rt_genes$gene_name[which(rt_genes$description == " serotonin receptor")]
gene_list = gene_list[which(!gene_list %in% c("Htr2b", "Htr3a", "Htr3b"))]

fpkm_matrix_filt = fpkm_matrix[gene_list,]

fpkm_matrix_filt$gene <- rownames(fpkm_matrix_filt)
fpkm_matrix_m = melt(fpkm_matrix_filt)
fpkm_matrix_m$Timepoint = NA
fpkm_matrix_m$CellType = NA
fpkm_matrix_m$Treatment = NA
fpkm_matrix_m$CellType[which(fpkm_matrix_m$variable %in% meta_table$Sample.ID[meta_table$CellType == "PV"])] = "PV"
fpkm_matrix_m$CellType[which(fpkm_matrix_m$variable %in% meta_table$Sample.ID[meta_table$CellType == "Cux2ERT"])] = "Cux2ERT"
fpkm_matrix_m$CellType[which(fpkm_matrix_m$variable %in% meta_table$Sample.ID[meta_table$CellType == "Rbp4"])] = "Rbp4"

fpkm_matrix_m$Treatment[which(fpkm_matrix_m$variable %in% meta_table$Sample.ID[meta_table$Treatment == "DOI"])] = "DOI"
fpkm_matrix_m$Treatment[which(fpkm_matrix_m$variable %in% meta_table$Sample.ID[meta_table$Treatment == "Ctrl"])] = "Ctrl"

fpkm_matrix_m$Timepoint[which(fpkm_matrix_m$variable %in% meta_table$Sample.ID[meta_table$Timepoint == "48hrs"])] = "48hrs"
fpkm_matrix_m$Timepoint[which(fpkm_matrix_m$variable %in% meta_table$Sample.ID[meta_table$Timepoint == "1week"])] = "1week"
fpkm_matrix_m$Timepoint[which(fpkm_matrix_m$variable %in% meta_table$Sample.ID[meta_table$Timepoint == "1month"])] = "1month"
head(fpkm_matrix_m)
table(fpkm_matrix_m$CellType)


# https://www.datanovia.com/en/blog/how-to-add-p-values-to-ggplot-facets/
table(fpkm_matrix_m$CellType)
cell_name = "Rbp4"
fpkm_matrix_m_ = fpkm_matrix_m[fpkm_matrix_m$CellType == cell_name,]

head(fpkm_matrix_m_)
stat.test <- fpkm_matrix_m_ %>%
  group_by(gene, Timepoint) %>%
  wilcox_test(value ~ Treatment) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test

txt_title = paste(cell_name, " samples - wilcox_test BH adjustment")
# Box plots
bxp <- ggboxplot(
  fpkm_matrix_m_, x = "Timepoint", y = "value", color = "darkgrey", fill = "Treatment", add = c("dotplot"), size = 1, scales = "free",
  palette = figure_colors$Treatment, facet.by = "gene", title = txt_title
)
bxp 

# Box plots with p-values
# Hide ns (non-significant)
stat.test <- stat.test %>%
  add_xy_position(x = "Timepoint", dodge = 0.8)
stat.test$y.position = 2

stat.test$round_p.adj = round(stat.test$p.adj, 2)

p1 <- bxp + 
  stat_pvalue_manual(
    stat.test, tip.length = 0.01,
    label = "{round_p.adj} {p.adj.signif}",
    hide.ns = FALSE
  ) +
  scale_y_continuous(limits=c(0,9), breaks=seq(0,9,2), expand = expansion(mult = c(0.01, 0.1))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5))
p1

export_pdf(p1, dir = "figures/drafts/", file_name = paste0(cell_name, "_serotonin_FPKM"), h = 10, w = 10)

#################################


cell_name = "Saline"
fpkm_matrix_m_ = fpkm_matrix_m[fpkm_matrix_m$Treatment == "Ctrl",]
fpkm_matrix_m_$Timepoint <- factor(fpkm_matrix_m_$Timepoint, levels = c("48hrs", "1week", "1month"))
fpkm_matrix_m_$CellType <- factor(fpkm_matrix_m_$CellType, levels = c("PV", "Cux2ERT", "Rbp4"))

head(fpkm_matrix_m_)
stat.test <- fpkm_matrix_m_ %>%
  group_by(gene, CellType) %>%
  t_test(value ~ Timepoint)  %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
stat.test

txt_title = paste(cell_name, " samples - t-test test BH adjusted p-value")

# Box plots
bxp <- ggboxplot(
  fpkm_matrix_m_, x = "CellType", y = "value", color = "black", fill = "Timepoint", size = 1, scales = "free", outlier.shap=NA, add=c('dotplot', 'mean_sd'),
  palette = figure_colors$Timepoint, facet.by = "gene", title = txt_title
)
bxp 

stat.test <- stat.test %>%
  add_xy_position(x = "CellType", dodge = 0.8)

stat.test$y.position = ifelse(stat.test$y.position > 8, 9, stat.test$y.position)

stat.test$round_p.adj = round(stat.test$p.adj, 2)

p1 <- bxp + 
  stat_pvalue_manual(
    stat.test, tip.length = 0.01,
    label = "{round_p.adj} {p.adj.signif}",
    hide.ns = FALSE
  ) +
  scale_y_continuous(limits=c(0,9), breaks=seq(0,9,2)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.9)) 
p1


export_pdf(p1, dir = "figures/drafts/", file_name = paste0(cell_name, "_by_timepoint_serotonin_FPKM"), h = 10, w = 10)

########################################################################################## Splicing factors expression
#### Splicing factors expression

sf_genes = get_splicing_factor_genes()
length(sf_genes)
sf_genes <- sf_genes[which(sf_genes %in% resOrdered_df$Genes)]
length(sf_genes)

vsd <- vst(dds, blind = FALSE)
Z <- t(scale(t(assay(vsd))))
mat <- Z[sf_genes, ]

breaksList = c(seq(-3, 3, by = 0.1))
breaksList = unique(breaksList)
my_ramp = colorRampPalette(c("#2e4a95", "black", "#f5ec03"))(length(breaksList))

p1 = pheatmap(mat,
              cluster_rows = T, 
              cluster_cols = F,
              show_rownames = T,
              show_colnames = F,
              border_color = F,
              annotation_col = meta_table[, c(2:6)],
              breaks = breaksList,
              color=my_ramp,
              annotation_colors = my_colour,
              main=paste0(nrow(mat), " top differentially expressed genes"
              )
)
p1


export_pdf(p1, dir = "figures/drafts/", file_name = "heatmap_splicing_factors", h = 9, w = 9)

########################################################################################## Receptors/transporter expression
#### 

txt_title = "Receptors, transporter \n z-score expression"
rt_genes = get_Receptors_T_genes()

cell_type = c("PV", "Rbp4", "Cux2ERT") # [PV Rbp4 Cux2ERT]
timepoint = c("48hrs", "1week", "1month") # [48hrs 1week 1month]

meta_table = get_metatable(cell_type, timepoint)
meta_table = meta_table[meta_table$Treatment == "Ctrl",]
plot_table(meta_table)

order_vector = c("PV", "Cux2ERT", "Rbp4")

meta_table <- meta_table %>%
  group_by(CellType) %>%
  mutate(order = factor(CellType, levels = order_vector)) %>%
  arrange(order, desc(Timepoint)) %>%
  ungroup() %>%
  select(-order)
meta_table = as.data.frame(meta_table)
rownames(meta_table) <- meta_table$Sample.ID

##
y = get_DGEList_counts(meta_table)
y = filter_genes(y)

logcounts <- cpm(y, normalized.lib.sizes = TRUE, log=TRUE)
head(logcounts)
plot_boxplot_topN_genes(logcounts, 30, ylimits=c(-5, 12), "logCPM raw_data")

ZlogCPM <- t(scale(t(logcounts)))
head(ZlogCPM)
plot_boxplot_topN_genes(ZlogCPM, 30, ylimits=c(-4, 4), "Z-score data")

expr = ZlogCPM[rownames(ZlogCPM) %in% rt_genes$gene_name,]
r_order = match(rt_genes$gene_name, rownames(expr))
r_order = r_order[!is.na(r_order)]

expr = expr[rev(r_order), meta_table$Sample.ID]
rownames(rt_genes) <- rt_genes$gene_name

p1 = pheatmap(expr,
              cluster_rows = T, 
              cluster_cols = F,
              show_rownames = T,
              annotation_col = meta_table[, c(2:6)],
              annotation_colors = my_colour,
              color=colorRampPalette(c("#2e4a95", "white", "#f5ec03"))(50),
              main=txt_title,
              border_color = F)
p1
export_svg(p1, dir = "figures/drafts/", file_name = "Heatmap_Receptors_transporters", h = 8, w = 9)

export_pdf(p1, dir = "figures/drafts/", file_name = "Heatmap_Receptors_transporters", h = 8, w = 9)

########################################################################################################################
##### GSEA ####
library(fgsea)
set.seed(12345)

cell_type = c("PV") # [PV Cux2ERT Rbp4]
timepoint = c("48hrs") # [48hrs 1week 1month]
treat = "DOI" # [DOI Psilocybin Ctrl]

res <- read_tsv(paste0("results/DEG/Episode VI/pipe6/", cell_type, "/", timepoint, "/DEG_", treat, " vs Saline - ", cell_type, " ", timepoint, " pipe6.txt"), show_col_types = FALSE)
res2 <- res %>% 
  dplyr::select(Genes, stat) %>%
  na.omit()
ranks <- deframe(res2)

#uncomment depending on which gene set you want to use as reference
#reactome.pathways <- gmtPathways("data/GO_annot_gsea-msigdb/m2.cp.reactome.v2024.1.Mm.symbols.gmt")
#GO.pathways <- gmtPathways("data/GO_annot_gsea-msigdb/m5.go.bp.v2024.1.Mm.symbols.gmt")

#changing background set of genes https://biostatsquid.com/fgsea-tutorial-gsea/
my_genes <- res2$Genes

### run this once ###
matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}

gmt <- gmtPathways("data/GO_annot_gsea-msigdb/m5.go.bp.v2024.1.Mm.symbols.gmt")
length(gmt) #number of GO terms
hidden <- unique(unlist(gmt))
length(hidden) #number of genes in all go terms

# Convert gmt file to a matrix with the genes as rows and for each go annotation (columns) the values are 0 or 1
mat <- matrix(NA, dimnames = list(hidden, names(gmt)),
              nrow = length(hidden), ncol = length(gmt))
for (i in 1:dim(mat)[2]){
  mat[,i] <- as.numeric(hidden %in% gmt[[i]])
}
### run this once ###

#Subset to the genes that are present in our data to avoid bias
hidden1 <- intersect(my_genes, hidden)
mat2 <- mat[hidden1, colnames(mat)[which(colSums(mat[hidden1,])>5)]] # filter for gene sets with more than 5 genes annotated

# And get the list again using the function we previously defined
final_list <- matrix_to_list(mat2)
#length(final_list) #number of GO terms left
#length(unique(unlist(final_list))) #number of genes in leftover go terms

fgseaRes <- fgsea(pathways=final_list, stats=ranks, minSize = 15, maxSize = 200)
sig_pathways <- fgseaRes[padj < .01,]

dim(sig_pathways) #number of significant pathways w/ cutoff of .01
length(unique(unlist(sig_pathways$leadingEdge))) #number of unique genes in leading edge of significant pathways

current_genes <- unique(unlist(sig_pathways$leadingEdge))
current_genes <- as.data.frame(current_genes)
current_genes$log2FC <- res$log2FoldChange[match(current_genes$current_genes, res$Genes)]
current_genes <- current_genes %>%
  arrange(log2FC)
dim(current_genes)
genes_current2 <- current_genes[abs(current_genes$log2FC) > .2,]
dim(genes_current2)
#genes_current = unique(unlist(res$leadingEdge))
#genes_current2 = unique(unlist(res2$leadingEdge))
#genes_all <- genes_current2$current_genes
#genes_all <- append(genes_all, genes_current2$current_genes)
#genes_all <- unique(genes_all)
#length(genes_all)

save_file(as.data.frame(genes_current2), 
          dir = "results/GSEA/",
          file_name = paste0(cell_type, timepoint, treat, " FC Genes.txt"),
          row = FALSE,
          sep = "\t")

sig_pathwaysTidy <- sig_pathways %>%
  as_tibble() %>%
  arrange(padj)

# fgseaResTidy <- fgseaRes %>%
#   as_tibble() %>%
#   arrange(padj)

#table output
sig_pathwaysTidy %>% 
  DT::datatable()

## for visualization  
#plot table of top 10 up and down-regulated pathways
topPathwaysUp <- fgseaRes[ES > 0][head(order(padj), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(padj), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(final_list[topPathways], ranks, fgseaRes, 
              gseaParam=0.5)

#plot top 10 up and down-regulated pathways w/ NES
top_fgseaRes_up <- fgseaRes[ES > 0][head(order(padj), n=10), ]
top_fgseaRes_down <- fgseaRes[ES < 0][head(order(padj), n=10), ]
top_fgseaRes <- bind_rows(top_fgseaRes_up, top_fgseaRes_down)
p <- ggplot(top_fgseaRes, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.01)) +
  scale_fill_manual(values=c("#56BBC1", "#E77D72")) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GO Normalized Enrichment Score from GSEA") + 
  theme_minimal()
p
export_svg(p, dir = "results/GSEA/Figures/", file_name = paste0(cell_type, timepoint, treat, " Top 20 Enriched GO terms"), h = 5, w = 10)

#plot individual pathways
pathway_interest = "GOBP_CELL_CELL_ADHESION_VIA_PLASMA_MEMBRANE_ADHESION_MOLECULES"
plotEnrichment(final_list[[pathway_interest]],
               ranks) + labs(title=pathway_interest)

#OR

library(easybio)
p <- plotGSEA(fgseaRes, pathways = final_list, pwayname = pathway_interest, stats = ranks, save = FALSE)
p
export_svg(p, dir = "results/GSEA/Figures/", file_name = paste0(cell_type, timepoint, treat, "plotGSEA_", pathway_interest), h = 5, w = 6)

## heatmap if dds is saved from DESeq2 
pathway_genes <- unlist(fgseaRes[pathway == pathway_interest,leadingEdge])
vsd <- vst(dds, blind = FALSE)
Z <- t(scale(t(assay(vsd))))
mat <- Z[head(pathway_genes,30), ]

p1 = pheatmap(mat,
              cluster_rows = T, 
              cluster_cols = F,
              show_rownames = T,
              show_colnames = F,
              border_color = NA,
              annotation_col = meta_table[, c(2:6)],
              annotation_colors = my_colour,
              color=colorRampPalette(c("#2e4a95", "black", "#f5ec03"))(50),
              main=paste0(nrow(mat), " genes in ", pathway_interest)
)
p1

export_svg(p1, dir = "results/GSEA/Figures/", file_name = paste0(cell_type, timepoint, treat, "Heatmap ", pathway_interest), h = 6, w = 6)

library(data.table)
file_name = paste0("FGSEA_", cell_type, "_", timepoint, "_", treat, "_all.txt")
file_name
fwrite(sig_pathwaysTidy, file=file_name, sep="\t", sep2=c("", " ", ""))

########################################################################################################
##### fig 2D log2FC heatmap all FC genes #####

aux <- read_table("results/GSEA/PV48hrsPsilocybin FC Genes.txt", show_col_types = F)
union_set <- as.data.frame(aux$current_genes)
aux <- read_table("results/GSEA/PV48hrsDOI FC Genes.txt", show_col_types = F)
union_set <- rbind(union_set,as.data.frame(aux$current_genes))
aux <- read_table("results/GSEA/PV1weekDOI FC Genes.txt", show_col_types = F)
union_set <- rbind(union_set,as.data.frame(aux$current_genes))
aux <- read_table("results/GSEA/Cux2ERT48hrsPsilocybin FC Genes.txt", show_col_types = F)
union_set <- rbind(union_set,as.data.frame(aux$current_genes))
aux <- read_table("results/GSEA/Cux2ERT48hrsDOI FC Genes.txt", show_col_types = F)
union_set <- rbind(union_set,as.data.frame(aux$current_genes))
aux <- read_table("results/GSEA/Cux2ERT1weekDOI FC Genes.txt", show_col_types = F)
union_set <- rbind(union_set,as.data.frame(aux$current_genes))
aux <- read_table("results/GSEA/Cux2ERT1monthDOI FC Genes.txt", show_col_types = F)
union_set <- rbind(union_set,as.data.frame(aux$current_genes))
aux <- read_table("results/GSEA/Rbp448hrsPsilocybin FC Genes.txt", show_col_types = F)
union_set <- rbind(union_set,as.data.frame(aux$current_genes))
aux <- read_table("results/GSEA/Rbp448hrsDOI FC Genes.txt", show_col_types = F)
union_set <- rbind(union_set,as.data.frame(aux$current_genes))
aux <- read_table("results/GSEA/Rbp41monthDOI FC Genes.txt", show_col_types = F)
union_set <- rbind(union_set,as.data.frame(aux$current_genes))
dim(union_set)
union_set <- unique(union_set)
dim(union_set)

data_FC = data.frame(row.names = unique(union_set$`aux$current_genes`))
dim(data_FC)

data_FC = get_log2FC_values(data_FC, loc_res)
head(data_FC)

meta_table = data.frame(row.names = c('PV_48hrs_DOI', 'PV_48hrs_Psilocybin', 'PV_1week_DOI', 'PV_1month_DOI', 'Cux2ERT_48hrs_DOI', 'Cux2ERT_48hrs_Psilocybin', 'Cux2ERT_1week_DOI', 'Cux2ERT_1month_DOI', 'Rbp4_48hrs_DOI', 'Rbp4_48hrs_Psilocybin', 'Rbp4_1week_DOI', 'Rbp4_1month_DOI'))
meta_table = cbind(meta_table, Timepoint=c('48hrs', '48hrs', '1week', '1month', '48hrs', '48hrs','1week', '1month', '48hrs','48hrs',  '1week', '1month'))
meta_table = cbind(meta_table, CellType=c('PV', 'PV', 'PV', 'PV', 'Cux2ERT', 'Cux2ERT', 'Cux2ERT', 'Cux2ERT','Rbp4', 'Rbp4','Rbp4', 'Rbp4'))
meta_table = cbind(meta_table, Treatment=c('DOI', 'Psilocybin', 'DOI', 'DOI', 'DOI', 'Psilocybin','DOI', 'DOI', 'DOI','Psilocybin', 'DOI', 'DOI'))

meta_table = meta_table[rownames(meta_table) %in% colnames(data_FC),]
meta_table = meta_table[colnames(data_FC),]

breaksList = c(seq(-2, -1, by = 0.5), seq(-1, 1, 0.1), seq(1, 2, 0.5))
breaksList = unique(breaksList)

my_ramp = colorRampPalette(c("blue", "white", "orange"))(length(breaksList))
#my_ramp[4] = "#FFFFFF"

p1 = pheatmap(data_FC,
              cluster_rows = F, 
              cluster_cols = F,
              show_rownames = F,
              show_colnames = F,
              annotation_col = figure_colors,
              annotation_colors = my_colour,
              breaks = breaksList,
              color=my_ramp,
              border_color = NA
)
p1

export_svg(p1, dir = "results/GSEA/", file_name = "log2FC_Heatmap_GSEA_Genes_reordered", h = 8, w = 8)


####################################################################################################
## supp figure 5-G volcano plots
#### PV
data <- read_table("results/DEG/Episode VI/pipe6/PV/48hrs/DEG_Psilocybin vs Saline - PV 48hrs pipe6.txt", show_col_types = F)
table(data$Expression)
head(data)
title_condition = 'PsilovsSaline PV 48hrs pipe6'

data <- read_table("results/DEG/Episode VI/pipe6/PV/48hrs/DEG_DOI vs Saline - PV 48hrs pipe6.txt", show_col_types = F)
table(data$Expression)
title_condition = 'DOIvsSaline PV 48hrs pipe6'

data <- read_table("results/DEG/Episode VI/pipe6/PV/1week/DEG_DOI vs Saline - PV 1week pipe6.txt", show_col_types = F)
table(data$Expression)
title_condition = 'DOIvsSaline PV 1week pipe6'

data <- read_table("results/DEG/Episode VI/pipe6/PV/1month/DEG_DOI vs Saline - PV 1month pipe6.txt", show_col_types = F)
table(data$Expression)
title_condition = 'DOIvsSaline PV 1month pipe6'

#### Cux2
data <- read_table("results/DEG/Episode VI/pipe6/Cux2ERT/48hrs/DEG_Psilocybin vs Saline - Cux2ERT 48hrs pipe6.txt", show_col_types = F)
table(data$Expression)
title_condition = 'PsilovsSaline Cux2ERT 48hrs pipe6'

data <- read_table("results/DEG/Episode VI/pipe6/Cux2ERT/48hrs/DEG_DOI vs Saline - Cux2ERT 48hrs pipe6.txt", show_col_types = F)
table(data$Expression)
title_condition = 'DOIvsSaline Cux2ERT 48hrs pipe6'

data <- read_table("results/DEG/Episode VI/pipe6/Cux2ERT/1week/DEG_DOI vs Saline - Cux2ERT 1week pipe6.txt", show_col_types = F)
table(data$Expression)
title_condition = 'DOIvsSaline Cux2ERT 1week pipe6'

data <- read_table("results/DEG/Episode VI/pipe6/Cux2ERT/1month/DEG_DOI vs Saline - Cux2ERT 1month pipe6.txt", show_col_types = F)
table(data$Expression)
title_condition = 'DOIvsSaline Cux2ERT 1month pipe6'

#### Rbp4
data <- read_table("results/DEG/Episode VI/pipe6/Rbp4/48hrs/DEG_Psilocybin vs Saline - Rbp4 48hrs pipe6.txt", show_col_types = F)
table(data$Expression)
title_condition = 'PsilovsSaline Rbp4 48hrs pipe6'

data <- read_table("results/DEG/Episode VI/pipe6/Rbp4/48hrs/DEG_DOI vs Saline - Rbp4 48hrs pipe6.txt", show_col_types = F)
table(data$Expression)
title_condition = 'DOIvsSaline Rbp4 48hrs pipe6'

data <- read_table("results/DEG/Episode VI/pipe6/Rbp4/1week/DEG_DOI vs Saline - Rbp4 1week pipe6.txt", show_col_types = F)
table(data$Expression)
title_condition = 'DOIvsSaline Rbp4 1week pipe6'

data <- read_table("results/DEG/Episode VI/pipe6/Rbp4/1month/DEG_DOI vs Saline - Rbp4 1month pipe6.txt", show_col_types = F)
table(data$Expression)
title_condition = 'DOIvsSaline Rbp4 1month pipe6'

top_genes = get_topN_DEG_v2(data, 10)

obj = plot_volcano_v2(data, top_genes, title_condition, p_pv = 0.1, p_fc = 0)
obj

export_svg(obj,
           dir = 'figures/drafts/',
           file_name = paste(title_condition),
           h = 11,
           w = 12
)


####################################################################################################
##                    Deconvolution analysis

mPFC_saline = get_scRNA_addition_dataset()
head(mPFC_saline@meta.data)
table(mPFC_saline@meta.data$CellType)

prop_cells = data.frame(table(mPFC_saline@meta.data$CellType))
plot_cell_prop(prop_cells)

prop_cells = data.frame(table(mPFC_saline@meta.data$L1_clusters))
plot_cell_prop(prop_cells)


Idents(mPFC_saline) <- mPFC_saline@meta.data$CellType
single.cell.expression.set <- as.SingleCellExperiment(mPFC_saline)

##
v1 = get_V1_dataset()
head(v1@meta.data)
table(v1@meta.data$maintype)
table(v1@meta.data$stim)

prop_cells = data.frame(table(v1@meta.data$maintype))
p1 = plot_cell_prop(prop_cells)
export_pdf(obj = p1, dir = "figures/drafts/", file_name = "V1_number_of_cells", h = 5.5, w = 5.5)

prop_cells = data.frame(table(v1@meta.data$celltype))
plot_cell_prop(prop_cells)

table(v1@meta.data$sample)

Idents(v1) <- v1@meta.data$maintype
single.cell.expression.set <- as.SingleCellExperiment(v1)


######################################### Ribo-seq

cell_type = c("PV", "Rbp4", "Cux2ERT") # [PV Rbp4 Cux2ERT]
timepoint = c("48hrs", "1week", "1month") # [48hrs 1week 1month]
treat = "DOI" # [DOI Psilocybin]

meta_table = get_metatable(cell_type, timepoint, treat)
plot_table(meta_table)
head(meta_table)

matrix = get_counts_RNA(meta_table)

sd <- apply(X = matrix, MARGIN = 1, FUN = sd)
rank <- (length(sd) - rank(sd) + 1)
select <- rank <= 500

pheatmap(as.matrix(log2(matrix[select,] +1)), cluster_rows=T, show_rownames=F, cluster_cols=T, annotation_col = meta_table[, c(2:6)], 
         main = paste0('Hisao samples Top ',sum(select),' protein coding genes (n=', nrow(matrix), ")"))

pData = new("AnnotatedDataFrame",data=meta_table[, c(2:6)])
rownames(pData)
exprs = as.matrix(matrix)
all(rownames(pData)==colnames(exprs))

matrix.eset = ExpressionSet(assayData=exprs, phenoData = pData)

###

scdata = "mPFC_addiction_study"
cluster = "CellType"
cluster = "L1_clusters"
subset = "all_cells"
bulk = "Hisao"

scdata = "V1_study_all_stim"
cluster = "maintype"
cluster = "celltype"
subset = "all_cells"
bulk = "Hisao"
bulk = "de_la_fuente_revenga"


Est.prop = music_prop(bulk.mtx = exprs(matrix.eset),
                      sc.sce = single.cell.expression.set, 
                      clusters = cluster,
                      samples = 'sample',
                      verbose = T)



est.prop = Est.prop$Est.prop.weighted
est.prop[order(row.names(est.prop)), ]
est.prop

colnames(est.prop)[which(is.na(colnames(est.prop)))] <- "Unknown"

## L1 subclusters
est.prop = est.prop[,colSums(est.prop) > 0]

bp.esti.pro = merge(est.prop, meta_table[, c(2:6)], by=0)
bp.esti.pro = bp.esti.pro[match(rownames(meta_table), bp.esti.pro$Row.names), ]

ext = "tsv"
res = "Est.prop.weighted"
fileout = paste(scdata, bulk, cluster, subset, res, ext, sep = ".")
fileout
write.table(bp.esti.pro, paste0("results/Decon_Music/", fileout), quote = F, sep = "\t", row.names = F)
res = "Weight.gene"
fileout = paste(scdata, bulk, cluster, subset, res, ext, sep = ".")
fileout
write.table(Est.prop$Weight.gene, paste0("results/Decon_Music/", fileout), quote = F, sep = "\t", row.names = T)

# 
# collorder = c("Astro" = "#8dd3c7", "Endo" = "#ffffb3", "Excitatory" = "#bebada", "Inhibitory" = "#fb8072", "Microglia" = "#80b1d3", 
#               "NF Oligo" = "#fdb462", "Oligo"= "#b3de69", "OPC" = "#fccde5", my_colour$CellType)

collorder = c("Astrocytes" = "#8dd3c7", "Endothelial_SmoothMuscle" = "#ffffb3", "Excitatory" = "#bebada", "Interneurons" = "#fb8072", "Microglia" = "#80b1d3", 
              "Mural" = "#fdb462", "Oligodendrocytes"= "#b3de69", "Macrophages" = "#fccde5", "Unkmown" = "white")

anno_width = unit(10, "cm")

my_colour = list(
  Treatment = c(DOI = "#5977ff", Ctrl = "#f74747", Psilocybin="#f5945c"),
  CellType = c(Rbp4 = "#0f7733", PV = "#a94498", Cux2ERT = "#45aa9a"),
  Sex = c('F' = "#f1b6da", 'M' = "#8073ac"),
  Batch = c('B1' = "#dfc27d", 'B2' = "#c7eae5", 'B3' = "#bac4e3", 'B4' = "#b3be62"),
  Timepoint = c("48hrs" = "#80afd2", "1week" = "#0072b3", "1month" = "#313695"),
  Cell_type = c('Int_Pv'= "#a94498", 'ExcL23' = "lightblue", 'ExcL5_1'= "darkgreen", 'ExcL5_2'= "darkgreen", 'ExcL5_3' = "darkgreen")
)

head(bp.esti.pro)
n_meta_columns = 5
.cols = c(2: (ncol(bp.esti.pro)-n_meta_columns))
## minor classes
collorder = rainbow(length(.cols), alpha = 0.6)

ht_list = rowAnnotation(text = anno_text(bp.esti.pro$Row.names, location = unit(1, "npc"), just = "right", 
                                         gp = gpar(fontsize = 12)))
ht_list = ht_list + rowAnnotation("dist_tss" = anno_barplot(bp.esti.pro[, .cols], bar_width = 1, gp = gpar(fill = collorder, fontsize = 14), 
                                                            width = anno_width), show_annotation_name = FALSE)
ht_list = ht_list + rowAnnotation(CellType = anno_simple(bp.esti.pro$CellType, col = my_colour$CellType))
ht_list = ht_list + rowAnnotation(Treatment = anno_simple(bp.esti.pro$Treatment, col = my_colour$Treatment))
ht_list = ht_list + rowAnnotation(Timepoint = anno_simple(bp.esti.pro$Timepoint, col = my_colour$Timepoint))
ht_list = ht_list + rowAnnotation(Batch = anno_simple(bp.esti.pro$Batch, col = my_colour$Batch))

p1 <- draw(ht_list, 
     heatmap_legend_list = Legend(title = "mPFC major cell types", labels = colnames(bp.esti.pro[, .cols]), legend_gp = gpar(fill = collorder)))

p1
export_svg(p1, dir = "figures/drafts/", file_name = "V1_all_stim_decon_Hisao_major", h = 13, w = 8)

rownames(est.prop) <- paste(rownames(est.prop), bp.esti.pro$CellType, bp.esti.pro$Treatment, bp.esti.pro$Timepoint, sep = "_")

prop_all = cbind('proportion'=c(est.prop), 'sampleID'=rep(rownames(est.prop),times=ncol(est.prop)), 'celltype'=rep(colnames(est.prop), each=nrow(est.prop)))
prop_all = as.data.frame(prop_all)
prop_all$proportion = as.numeric(as.character(prop_all$proportion))

prop_all$group = ifelse(grepl("Ctrl", prop_all$sampleID), 'Ctrl', 'DOI')
prop_all$timepoint[grepl("48hrs", prop_all$sampleID)] = '48hrs'
prop_all$timepoint[grepl("1month", prop_all$sampleID)] = '1month'
prop_all$timepoint[grepl("1week", prop_all$sampleID)] = '1week'
prop_all$CellType[grepl("PV", prop_all$sampleID)] = 'PV'
prop_all$CellType[grepl("Rbp4", prop_all$sampleID)] = 'Rbp4'
prop_all$CellType[grepl("Cux2", prop_all$sampleID)] = 'Cux2'

ggplot(prop_all, aes(x=celltype, y=proportion, color=celltype)) + xlab('')+
  geom_jitter(width=0.25,alpha=0.8)+ylab('Cell Type Proportions')+theme_bw()+
  stat_summary(fun = median,
               geom = "crossbar", width = 0.5,size=0.5,color='gray36')+
  facet_wrap(~ CellType, ncol = 1) +
  #facet_grid(.~group)+
  theme(plot.title = element_text(hjust = 0.5, size=12),
        axis.text.x = element_text(size=12,angle = 45,hjust=1),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.line = element_line(colour = "black"),
        strip.text.x = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')
####### de la fuente revenga

bulk = "de_la_fuente_revenga"

ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
bm <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'gene_biotype', 'ensembl_transcript_id'), mart = ensembl)
head(bm)

meta_table = read.delim("/media/data01/processed/public_data/RNA-Seq/GSE161626_de_la_Fuente_Revenga/meta_table.csv", sep = ",")
meta_table_RNA = meta_table[meta_table$Experiment == "RNA",]
meta_table_RNA$Timepoint = gsub(" ", "_", meta_table_RNA$Timepoint)
Group_GSE161626 = factor(paste(meta_table_RNA$Treatment,meta_table_RNA$Timepoint,sep="."))
meta_table_RNA = cbind(meta_table_RNA,Group=Group_GSE161626)
meta_table_RNA$Name[9] <- "24h6.1.RNA"
meta_table_RNA$Name[10] <- "24h6.2.RNA"
meta_table_RNA$Name[19] <- "48h6.1.RNA"
meta_table_RNA$Name[20] <- "48h6.2.RNA"
meta_table_RNA$Name[28] <- "7d6-2-RNA"
meta_table_RNA$Name[27] <- "7d6-1-RNA"
meta_table_RNA$Name <- gsub("-", ".", meta_table_RNA$Name)

path_data = "/media/data01/processed/public_data/RNA-Seq/GSE161626_de_la_Fuente_Revenga/GSE161626_RAW/output.txt"
counts_GSE161626 = read.delim(path_data, header = T, sep = "\t", comment.char = "#")

head(counts_GSE161626)
colnames(counts_GSE161626) = gsub("X.work.cascades.bz10.RNA_seq.data.VCU_psychedelics.VEH.Aligned_BAM.", "", colnames(counts_GSE161626))
colnames(counts_GSE161626) = gsub("X.work.cascades.bz10.RNA_seq.data.VCU_psychedelics.24h.Aligned_BAM.", "", colnames(counts_GSE161626))
colnames(counts_GSE161626) = gsub("X.work.cascades.bz10.RNA_seq.data.VCU_psychedelics.48h.Aligned_BAM.", "", colnames(counts_GSE161626))
colnames(counts_GSE161626) = gsub("X.work.cascades.bz10.RNA_seq.data.VCU_psychedelics.7d.Aligned_BAM.", "", colnames(counts_GSE161626))
colnames(counts_GSE161626) = gsub(".bam", "", colnames(counts_GSE161626))
head(counts_GSE161626)


rownames(counts_GSE161626) <- counts_GSE161626$Geneid
counts_GSE161626 <- counts_GSE161626[, -1]  # create the table with only counts here

matrix = counts_GSE161626
matrix = matrix[!duplicated(rownames(matrix)),]

bm_genes = bm[bm$ensembl_gene_id[!is.na(bm$gene_biotype)] %in% rownames(matrix), ]
rownames(matrix) = make.names(bm_genes$external_gene_name[match(rownames(matrix), bm_genes$ensembl_gene_id)], unique = TRUE)

meta_table_RNA$Name
colnames(counts_GSE161626)
rownames(meta_table_RNA) <- colnames(counts_GSE161626)

df = meta_table_RNA
df = df[df$Name %in% colnames(matrix),]
rownames(df) <- df$Name

sd <- apply(X = matrix, MARGIN = 1, FUN = sd)
rank <- (length(sd) - rank(sd) + 1)
select <- rank <= 500

pheatmap(log2(matrix[select,] +1),
         cluster_rows=T, 
         show_rownames=F, 
         cluster_cols=T, 
         annotation_col = df[, -c(1, 2)], 
         main = paste0('de La Fuente Revenga ',sum(select),' protein coding genes (n=', nrow(matrix), ")"))


pData = new("AnnotatedDataFrame",data=df)
rownames(pData)
exprs = as.matrix(matrix)
all(rownames(pData)==colnames(exprs))

matrix.eset = ExpressionSet(assayData=exprs, phenoData = pData)

##

scdata = "V1_study_all_stim"
cluster = "maintype"
cluster = "celltype"
subset = "all_cells"
bulk = "Hisao"
bulk = "de_la_fuente_revenga"

Est.prop = music_prop(bulk.mtx = exprs(matrix.eset),
                      sc.sce = single.cell.expression.set, 
                      clusters = cluster,
                      samples = 'sample',
                      verbose = T)

est.prop = Est.prop$Est.prop.weighted
est.prop[order(row.names(est.prop)), ]
est.prop
## L1 subclusters
est.prop = est.prop[,colSums(est.prop) > 0]
est.prop

colnames(est.prop)[which(is.na(colnames(est.prop)))] <- "Unknown"

bp.esti.pro = merge(est.prop, df[, c(2:4)], by=0)
bp.esti.pro = bp.esti.pro[match(rownames(df), bp.esti.pro$Row.names), ]

ext = "tsv"
res = "Est.prop.weighted"
fileout = paste(scdata, bulk, cluster, subset, res, ext, sep = ".")
fileout
write.table(bp.esti.pro, paste0("results/Decon_Music/", fileout), quote = F, sep = "\t", row.names = F)
res = "Weight.gene"
fileout = paste(scdata, bulk, cluster, subset, res, ext, sep = ".")
fileout
write.table(Est.prop$Weight.gene, paste0("results/Decon_Music/", fileout), quote = F, sep = "\t", row.names = T)


collorder = c("Astro" = "#8dd3c7", "Endo" = "#ffffb3", "Excitatory" = "#bebada", "Inhibitory" = "#fb8072", "Microglia" = "#80b1d3", 
              "NF Oligo" = "#fdb462", "Oligo"= "#b3de69", "OPC" = "#fccde5", my_colour$CellType)

collorder = c("Astrocytes" = "#8dd3c7", "Endothelial_SmoothMuscle" = "#ffffb3", "Excitatory" = "#bebada", "Interneurons" = "#fb8072", "Microglia" = "#80b1d3", 
              "Mural" = "#fdb462", "Oligodendrocytes"= "#b3de69", "Macrophages" = "#fccde5", "Unkmown" = "white")

anno_width = unit(10, "cm")


my_colour_df = list(
  Treatment = c(DOI = "#e54291", VEH = "#525252"),
  Timepoint = c("0" = "grey", "24_hours" = "#80afd2", "48_hours" = "#0072b3", "7_days" = "#313695")
)
head(bp.esti.pro)
n_meta_columns = 3
.cols = c(2: (ncol(bp.esti.pro)-n_meta_columns))
## minor 
collorder = rainbow(length(.cols), alpha = 0.6)

ht_list = rowAnnotation(text = anno_text(bp.esti.pro$Row.names, location = unit(1, "npc"), just = "right", 
                                         gp = gpar(fontsize = 12)))
ht_list = ht_list + rowAnnotation("dist_tss" = anno_barplot(bp.esti.pro[, .cols], bar_width = 1, gp = gpar(fill = collorder, fontsize = 14), 
                                                            width = anno_width), show_annotation_name = FALSE)
ht_list = ht_list + rowAnnotation(Treatment = anno_simple(bp.esti.pro$Treatment, col = my_colour_df$Treatment))
ht_list = ht_list + rowAnnotation(Timepoint = anno_simple(bp.esti.pro$Timepoint, col = my_colour_df$Timepoint))


p1 <- draw(ht_list, 
     heatmap_legend_list = Legend(title = "V1 major cell types", labels = colnames(bp.esti.pro[, .cols]), legend_gp = gpar(fill = collorder)))

export_svg(p1, dir = "figures/drafts/", file_name = "V1_all_stim_decon_delafuente_major", h = 6, w = 8)

#############################################################################################################

##### fig 2C dot plot #####
data1 = read_tsv('results/GSEA/FGSEA_PV_48hrs_DOI.txt', show_col_types = FALSE)
data1$leadingEdge <- strsplit(data1$leadingEdge," ")
for(i in 1:dim(data1)[1]){
  data1$size[i] = length(unlist(data1$leadingEdge[i]))
}
data1$Condition = "PV_48hrs_DOI"
dim(data1)
length(unique(unlist(data1$leadingEdge)))
data2 = read_tsv('results/GSEA/FGSEA_PV_48hrs_Psilocybin.txt', show_col_types = FALSE)
data2$leadingEdge <- strsplit(data2$leadingEdge," ")
for(i in 1:dim(data2)[1]){
  data2$size[i] = length(unlist(data2$leadingEdge[i]))
}
data2$Condition = "PV_48hrs_Psilo"
dim(data2)
length(unique(unlist(data2$leadingEdge)))
data3 = read_tsv('results/GSEA/FGSEA_PV_1week_DOI.txt', show_col_types = FALSE)
data3$leadingEdge <- strsplit(data3$leadingEdge," ")
for(i in 1:dim(data3)[1]){
  data3$size[i] = length(unlist(data3$leadingEdge[i]))
}
data3$Condition = "PV_1week_DOI"
dim(data3)
length(unique(unlist(data3$leadingEdge)))
data4 = read_tsv('results/GSEA/FGSEA_Cux2ERT_48hrs_DOI.txt', show_col_types = FALSE)
data4$leadingEdge <- strsplit(data4$leadingEdge," ")
for(i in 1:dim(data4)[1]){
  data4$size[i] = length(unlist(data4$leadingEdge[i]))
}
data4$Condition = "Cux2ERT_48hrs_DOI"
dim(data4)
length(unique(unlist(data4$leadingEdge)))
data5 = read_tsv('results/GSEA/FGSEA_Cux2ERT_48hrs_Psilocybin.txt', show_col_types = FALSE)
data5$leadingEdge <- strsplit(data5$leadingEdge," ")
for(i in 1:dim(data5)[1]){
  data5$size[i] = length(unlist(data5$leadingEdge[i]))
}
data5$Condition = "Cux2ERT_48hrs_Psilo"
dim(data5)
length(unique(unlist(data5$leadingEdge)))
data6 = read_tsv('results/GSEA/FGSEA_Cux2ERT_1week_DOI.txt', show_col_types = FALSE)
data6$leadingEdge <- strsplit(data6$leadingEdge," ")
for(i in 1:dim(data6)[1]){
  data6$size[i] = length(unlist(data6$leadingEdge[i]))
}
data6$Condition = "Cux2ERT_1week_DOI"
dim(data6)
length(unique(unlist(data6$leadingEdge)))
data7 = read_tsv('results/GSEA/FGSEA_Cux2ERT_1month_DOI.txt', show_col_types = FALSE)
data7$leadingEdge <- strsplit(data7$leadingEdge," ")
for(i in 1:dim(data7)[1]){
  data7$size[i] = length(unlist(data7$leadingEdge[i]))
}
data7$Condition = "Cux2ERT_1month_DOI"
dim(data7)
length(unique(unlist(data7$leadingEdge)))
data8 = read_tsv('results/GSEA/FGSEA_Rbp4_48hrs_DOI.txt', show_col_types = FALSE)
data8$leadingEdge <- strsplit(data8$leadingEdge," ")
for(i in 1:dim(data8)[1]){
  data8$size[i] = length(unlist(data8$leadingEdge[i]))
}
data8$Condition = "Rbp4_48hrs_DOI"
dim(data8)
length(unique(unlist(data8$leadingEdge)))
data9 = read_tsv('results/GSEA/FGSEA_Rbp4_48hrs_Psilocybin.txt', show_col_types = FALSE)
data9$leadingEdge <- strsplit(data9$leadingEdge," ")
for(i in 1:dim(data9)[1]){
  data9$size[i] = length(unlist(data9$leadingEdge[i]))
}
data9$Condition = "Rbp4_48hrs_Psilo"
dim(data9)
length(unique(unlist(data9$leadingEdge)))
data10 = read_tsv('results/GSEA/FGSEA_Rbp4_1month_DOI.txt', show_col_types = FALSE)
data10$leadingEdge <- strsplit(data10$leadingEdge," ")
for(i in 1:dim(data10)[1]){
  data10$size[i] = length(unlist(data10$leadingEdge[i]))
}
data10$Condition = "Rbp4_1month_DOI"
dim(data10)
length(unique(unlist(data10$leadingEdge)))

allpathways <- rbind(data1,data2,data3,data4,data5,data6,data7,data8,data9,data10)

# pathways_to_plot <- allpathways %>% 
#   filter(-log10(allpathways$padj) > 10)
# pathways_to_plot <- unique(pathways_to_plot$pathway)
# length(pathways_to_plot)
# pathways_to_plot <- allpathways[allpathways$pathway %in% pathways_to_plot,]

pathways_to_plot <- read_xlsx('files/Fig2C_GO_terms.xlsx',sheet = 3,col_names = FALSE)
dim(pathways_to_plot)
pathways_to_plot <- allpathways[which(allpathways$pathway %in% pathways_to_plot$...1),]
dim(pathways_to_plot)

#hypothetical data, just to include the time points with no pathways #
pathway <- c('GOBP_CYTOPLASMIC_TRANSLATION','GOBP_CYTOPLASMIC_TRANSLATION','GOBP_CYTOPLASMIC_TRANSLATION')
pval <- c(0,0,0)
padj <- c(0,0,0)
log2err <- c(.5,.5,.5)
ES <- c(.5,.5,.5)
NES <- c(.5,.5,.5)
size <- c(15,15,15)
leadingEdge <- c('Egr1','Egr1','Egr1')
Condition <- c('PV_1month_DOI','Rbp4_1week_DOI','Rbp4_1month_DOI')
made_up <- data.frame(pathway, pval, padj, log2err, ES, NES, size, leadingEdge, Condition)
print(made_up)

pathways_to_plot2 <- rbind(pathways_to_plot, made_up)
order <- as.data.frame(table(pathways_to_plot$pathway))
order <- arrange(order, Freq)

p1 <- ggplot(pathways_to_plot2, aes(x = factor(Condition, level = c('PV_48hrs_Psilo','PV_48hrs_DOI', 'PV_1week_DOI','PV_1month_DOI','Cux2ERT_48hrs_Psilo','Cux2ERT_48hrs_DOI','Cux2ERT_1week_DOI','Cux2ERT_1month_DOI','Rbp4_48hrs_Psilo','Rbp4_48hrs_DOI','Rbp4_1week_DOI','Rbp4_1month_DOI')), 
                                    y = factor(pathway, level = order$Var1),
                                    fill = sign(NES), 
                                    #alpha = -log10(padj),
                                    size = -log10(padj)
)) +
  geom_point(shape = 21, stroke = 0.5) + # change the thickness of the boarder with stroke
  scale_fill_gradient(low = "#0000f4", 
                      high = "#f2af42",) +
  scale_alpha(range = c(.3, 1)) +
  ylab("") +
  xlab("") + 
  guides(color = guide_colorbar(reverse = TRUE)) +
  theme_bw(base_size = 15) +
  theme(axis.text.y=element_text(colour="black"),
        axis.text.x=element_text(colour="black",angle = 90)) +
  labs(
    x = "Treatment",
    y = "",
    size = "-log10(padj)"
    #alpha = "-log10(padj)"
  )
p1

export_svg(p1, dir = "results/GSEA/", file_name = "dot plot GSEA pathways NEW", h = 6, w = 12)


#############################################################################################################
##### fig 2D log2FC heatmap all FC genes #####

aux <- read_table("results/GSEA/PV48hrsPsilocybin FC Genes.txt", show_col_types = F)
union_set <- as.data.frame(aux$current_genes)
aux <- read_table("results/GSEA/PV48hrsDOI FC Genes.txt", show_col_types = F)
union_set <- rbind(union_set,as.data.frame(aux$current_genes))
aux <- read_table("results/GSEA/PV1weekDOI FC Genes.txt", show_col_types = F)
union_set <- rbind(union_set,as.data.frame(aux$current_genes))
aux <- read_table("results/GSEA/Cux2ERT48hrsPsilocybin FC Genes.txt", show_col_types = F)
union_set <- rbind(union_set,as.data.frame(aux$current_genes))
aux <- read_table("results/GSEA/Cux2ERT48hrsDOI FC Genes.txt", show_col_types = F)
union_set <- rbind(union_set,as.data.frame(aux$current_genes))
aux <- read_table("results/GSEA/Cux2ERT1weekDOI FC Genes.txt", show_col_types = F)
union_set <- rbind(union_set,as.data.frame(aux$current_genes))
aux <- read_table("results/GSEA/Cux2ERT1monthDOI FC Genes.txt", show_col_types = F)
union_set <- rbind(union_set,as.data.frame(aux$current_genes))
aux <- read_table("results/GSEA/Rbp448hrsPsilocybin FC Genes.txt", show_col_types = F)
union_set <- rbind(union_set,as.data.frame(aux$current_genes))
aux <- read_table("results/GSEA/Rbp448hrsDOI FC Genes.txt", show_col_types = F)
union_set <- rbind(union_set,as.data.frame(aux$current_genes))
aux <- read_table("results/GSEA/Rbp41monthDOI FC Genes.txt", show_col_types = F)
union_set <- rbind(union_set,as.data.frame(aux$current_genes))
dim(union_set)
union_set <- unique(union_set)
dim(union_set)

data_FC = data.frame(row.names = unique(union_set$`aux$current_genes`))
dim(data_FC)

data_FC = get_log2FC_values(data_FC, loc_res)
head(data_FC)

meta_table = data.frame(row.names = c('PV_48hrs_DOI', 'PV_48hrs_Psilocybin', 'PV_1week_DOI', 'PV_1month_DOI', 'Cux2ERT_48hrs_DOI', 'Cux2ERT_48hrs_Psilocybin', 'Cux2ERT_1week_DOI', 'Cux2ERT_1month_DOI', 'Rbp4_48hrs_DOI', 'Rbp4_48hrs_Psilocybin', 'Rbp4_1week_DOI', 'Rbp4_1month_DOI'))
meta_table = cbind(meta_table, Timepoint=c('48hrs', '48hrs', '1week', '1month', '48hrs', '48hrs','1week', '1month', '48hrs','48hrs',  '1week', '1month'))
meta_table = cbind(meta_table, CellType=c('PV', 'PV', 'PV', 'PV', 'Cux2ERT', 'Cux2ERT', 'Cux2ERT', 'Cux2ERT','Rbp4', 'Rbp4','Rbp4', 'Rbp4'))
meta_table = cbind(meta_table, Treatment=c('DOI', 'Psilocybin', 'DOI', 'DOI', 'DOI', 'Psilocybin','DOI', 'DOI', 'DOI','Psilocybin', 'DOI', 'DOI'))

meta_table = meta_table[rownames(meta_table) %in% colnames(data_FC),]
meta_table = meta_table[colnames(data_FC),]


breaksList = c(seq(-2, -1, by = 0.5), seq(-1, 1, 0.1), seq(1, 2, 0.5))
breaksList = unique(breaksList)

my_ramp = colorRampPalette(c("blue", "white", "orange"))(length(breaksList))
#my_ramp[4] = "#FFFFFF"

p1 = pheatmap(data_FC,
              cluster_rows = F, 
              cluster_cols = F,
              show_rownames = F,
              show_colnames = F,
              annotation_col = meta_table,
              annotation_colors = figure_colors,
              breaks = breaksList,
              color=my_ramp,
              border_color = NA
)
p1

export_svg(p1, dir = "results/GSEA/", file_name = "log2FC_Heatmap_GSEA_Genes_reordered", h = 8, w = 8)

##############################################################################################
##### fig 2E UpSet plot ####

#PV
data1 = read_tsv('results/GSEA/FGSEA_PV_48hrs_DOI.txt', show_col_types = FALSE)
data1$leadingEdge <- strsplit(data1$leadingEdge," ")
data2 = read_tsv('results/GSEA/FGSEA_PV_48hrs_Psilocybin.txt', show_col_types = FALSE)
data2$leadingEdge <- strsplit(data2$leadingEdge," ")
data3 = read_tsv('results/GSEA/FGSEA_PV_1week_DOI.txt', show_col_types = FALSE)
data3$leadingEdge <- strsplit(data3$leadingEdge," ")
data4 = read_tsv('results/GSEA/FGSEA_PV_1month_DOI.txt', show_col_types = FALSE)
data4$leadingEdge <- strsplit(data4$leadingEdge," ")

head(data1)

x <- list(
  PV_48h_DOI = unique(unlist(data1$leadingEdge)),
  PV_48_Psilo = unique(unlist(data2$leadingEdge)),
  PV_1week_DOI = unique(unlist(data3$leadingEdge)),
  PV_1month_DOI = unique(unlist(data4$leadingEdge))
)

anno_table = data.frame(Timepoint = c("48hrs", "48hrs","1week","1month"),
                        Treatment = c("DOI", "Psilocybin", "DOI", "DOI"))

row_a = rowAnnotation(
  df = anno_table,
  col = figure_colors
)

### compare gene lists
x = lapply(x, function(z){ z[!is.na(z) & z != ""]})

m1 = make_comb_mat(x)
m1

#rev(order(comb_degree(m1)))

p1 <- UpSet(m1,top_annotation = upset_top_annotation(m1, add_numbers = TRUE),
            left_annotation = upset_left_annotation(m1, add_numbers = TRUE),  
            comb_order = rev(order(comb_size(m1))),
            right_annotation = row_a,
            comb_col = c("red", "blue", "black", "#525252")[comb_degree(m1)],
)
p1

export_svg(p1, dir = "results/GSEA/", file_name = "upset_plot_PV_GSEA", h = 3, w = 7)

#L2_3
data1 = read_tsv('results/GSEA/FGSEA_Cux2ERT_48hrs_DOI.txt', show_col_types = FALSE)
data1$leadingEdge <- strsplit(data1$leadingEdge," ")
data2 = read_tsv('results/GSEA/FGSEA_Cux2ERT_48hrs_Psilocybin.txt', show_col_types = FALSE)
data2$leadingEdge <- strsplit(data2$leadingEdge," ")
data3 = read_tsv('results/GSEA/FGSEA_Cux2ERT_1week_DOI.txt', show_col_types = FALSE)
data3$leadingEdge <- strsplit(data3$leadingEdge," ")
data4 = read_tsv('results/GSEA/FGSEA_Cux2ERT_1month_DOI.txt', show_col_types = FALSE)
data4$leadingEdge <- strsplit(data4$leadingEdge," ")

head(data1)

x <- list(
  Cux2ERT_48h_DOI = unique(unlist(data1$leadingEdge)),
  Cux2ERT_48_Psilo = unique(unlist(data2$leadingEdge)),
  Cux2ERT_1week_DOI = unique(unlist(data3$leadingEdge)),
  Cux2ERT_1month_DOI = unique(unlist(data4$leadingEdge))
)

anno_table = data.frame(Timepoint = c("48hrs", "48hrs","1week","1month"),
                        Treatment = c("DOI", "Psilocybin", "DOI", "DOI"))

row_a = rowAnnotation(
  df = anno_table,
  col = figure_colors
)

### compare gene lists
x = lapply(x, function(z){ z[!is.na(z) & z != ""]})

m1 = make_comb_mat(x)
m1


p1 <- UpSet(m1,top_annotation = upset_top_annotation(m1, add_numbers = TRUE),
            left_annotation = upset_left_annotation(m1, add_numbers = TRUE),  
            comb_order = rev(order(comb_size(m1))),
            right_annotation = row_a,
            comb_col = c("#a4bade", "#7a1cac", "#2e037f", "black")[comb_degree(m1)],
)
p1

export_svg(p1, dir = "results/GSEA/Figures/", file_name = "upset_plot_Cux2ERT_GSEA", h = 3, w = 7)
#L5
data1 = read_tsv('results/GSEA/FGSEA_Rbp4_48hrs_DOI.txt', show_col_types = FALSE)
data1$leadingEdge <- strsplit(data1$leadingEdge," ")
data2 = read_tsv('results/GSEA/FGSEA_Rbp4_48hrs_Psilocybin.txt', show_col_types = FALSE)
data2$leadingEdge <- strsplit(data2$leadingEdge," ")
data3 = read_tsv('results/GSEA/FGSEA_Rbp4_1week_DOI.txt', show_col_types = FALSE)
data3$leadingEdge <- strsplit(data3$leadingEdge," ")
data4 = read_tsv('results/GSEA/FGSEA_Rbp4_1month_DOI.txt', show_col_types = FALSE)
data4$leadingEdge <- strsplit(data4$leadingEdge," ")

head(data1)

x <- list(
  Rbp4_48h_DOI = unique(unlist(data1$leadingEdge)),
  Rbp4_48_Psilo = unique(unlist(data2$leadingEdge)),
  Rbp4_1week_DOI = unique(unlist(data3$leadingEdge)),
  Rbp4_1month_DOI = unique(unlist(data4$leadingEdge))
)

anno_table = data.frame(Timepoint = c("48hrs", "48hrs","1week","1month"),
                        Treatment = c("DOI", "Psilocybin", "DOI", "DOI"))

row_a = rowAnnotation(
  df = anno_table,
  col = figure_colors
)

### compare gene lists
x = lapply(x, function(z){ z[!is.na(z) & z != ""]})

m1 = make_comb_mat(x)
m1

p1 <- UpSet(m1,top_annotation = upset_top_annotation(m1, add_numbers = TRUE),
            left_annotation = upset_left_annotation(m1, add_numbers = TRUE),  
            comb_order = rev(order(comb_size(m1))),
            right_annotation = row_a,
            comb_col = c("#a4bade", "#7a1cac", "#2e037f", "black")[comb_degree(m1)],
)
p1

export_svg(p1, dir = "results/GSEA/Figures/", file_name = "upset_plot_Rbp4_GSEA", h = 3, w = 7)

##############################################################################################
##### fig 3A rMATs filtering #####

loc_res = read.delim("files/AS_results_location-copy.csv", sep = ",", comment.char = "#")
head(loc_res)

bar_data <- get_summary_pie_values(loc_res, celltype = "PV")
bar_data

euler_data = aggregate(bar_data$value, by=list(group=bar_data$group), FUN=sum)
colnames(euler_data)[2] <- 'value'

plot_euler <- function(pie_data) {
  fit = euler(c('total' = pie_data$value[4],
                'total&events_passed_lc' = pie_data$value[1],
                'total&events_passed_lc&Pval.adjust' = pie_data$value[3],
                'total&events_passed_lc&Pval.adjust&PSI.filter' = pie_data$value[2]))
  plot(fit,
       fills=list(fill = c("#f7f7f7", "#fde0ef", "#f1b6da","#de77ae")),
       labels=list(col="black", font=4), quantities=T)
}

p1 = plot_euler(euler_data)
p1
export_svg(p1, dir = "figures/drafts/rMATS/", file_name = "PV48DOI_AS_events_euler_plot", h = 5, w = 5)

##############################################################################################
##### supp fig 3 GSEA genes/GO term overlaps ####
### PV
# all genes
PV <- read_table("results/GSEA/All PV Genes.txt", show_col_types = F)
PV <- unique(PV)
dim(PV)
# FC genes only 
PV <- read_table("results/GSEA/PV48hrsPsilocybin FC Genes.txt", show_col_types = F)
aux <- read_table("results/GSEA/PV48hrsDOI FC Genes.txt", show_col_types = F)
PV <- rbind(PV,aux)
aux <- read_table("results/GSEA/PV1weekDOI FC Genes.txt", show_col_types = F)
PV <- rbind(PV,aux)
dim(PV)
PV <- PV$current_genes
PV <- unique(PV)
length(PV)

# all GO terms
PV_GO <- read_tsv('results/GSEA/FGSEA_PV_48hrs_DOI.txt', show_col_types = FALSE)
aux <- read_tsv('results/GSEA/FGSEA_PV_48hrs_Psilocybin.txt', show_col_types = FALSE)
PV_GO <- rbind(PV_GO, aux)
aux <- read_tsv('results/GSEA/FGSEA_PV_1week_DOI.txt', show_col_types = FALSE)
PV_GO <- rbind(PV_GO, aux)
aux <- read_tsv('results/GSEA/FGSEA_PV_1month_DOI.txt', show_col_types = FALSE)
PV_GO <- rbind(PV_GO, aux)
PV_GO <- PV_GO$pathway
PV_GO <- unique(PV_GO)
length(PV_GO)

### Cux2
Cux2 <- read_table("results/GSEA/All Cux2ERT Genes.txt", show_col_types = F)
Cux2 <- unique(Cux2)
dim(Cux2)

Cux2ERT <- read_table("results/GSEA/Cux2ERT48hrsPsilocybin FC Genes.txt", show_col_types = F)
aux <- read_table("results/GSEA/Cux2ERT48hrsDOI FC Genes.txt", show_col_types = F)
Cux2ERT <- rbind(Cux2ERT,aux)
aux <- read_table("results/GSEA/Cux2ERT1weekDOI FC Genes.txt", show_col_types = F)
Cux2ERT <- rbind(Cux2ERT,aux)
aux <- read_table("results/GSEA/Cux2ERT1monthDOI FC Genes.txt", show_col_types = F)
Cux2ERT <- rbind(Cux2ERT,aux)
dim(Cux2ERT)
Cux2ERT <- Cux2ERT$current_genes
Cux2ERT <- unique(Cux2ERT)
length(Cux2ERT)

Cux2_GO <- read_tsv('results/GSEA/FGSEA_Cux2ERT_48hrs_DOI.txt', show_col_types = FALSE)
aux <- read_tsv('results/GSEA/FGSEA_Cux2ERT_48hrs_Psilocybin.txt', show_col_types = FALSE)
Cux2_GO <- rbind(Cux2_GO, aux)
aux <- read_tsv('results/GSEA/FGSEA_Cux2ERT_1week_DOI.txt', show_col_types = FALSE)
Cux2_GO <- rbind(Cux2_GO, aux)
aux <- read_tsv('results/GSEA/FGSEA_Cux2ERT_1month_DOI.txt', show_col_types = FALSE)
Cux2_GO <- rbind(Cux2_GO, aux)
Cux2_GO <- Cux2_GO$pathway
Cux2_GO <- unique(Cux2_GO)
length(Cux2_GO)

### Rbp4
Rbp4 <- read_table("results/GSEA/All Rbp4 Genes.txt", show_col_types = F)
Rbp4 <- unique(Rbp4)
dim(Rbp4)

Rbp4 <- read_table("results/GSEA/Rbp448hrsPsilocybin FC Genes.txt", show_col_types = F)
aux <- read_table("results/GSEA/Rbp448hrsDOI FC Genes.txt", show_col_types = F)
Rbp4 <- rbind(Rbp4,aux)
aux <- read_table("results/GSEA/Rbp41monthDOI FC Genes.txt", show_col_types = F)
Rbp4 <- rbind(Rbp4,aux)
dim(Rbp4)
Rbp4 <- Rbp4$current_genes
Rbp4 <- unique(Rbp4)
length(Rbp4)

Rbp4_GO <- read_tsv('results/GSEA/FGSEA_Rbp4_48hrs_DOI.txt', show_col_types = FALSE)
aux <- read_tsv('results/GSEA/FGSEA_Rbp4_48hrs_Psilocybin.txt', show_col_types = FALSE)
Rbp4_GO <- rbind(Rbp4_GO, aux)
aux <- read_tsv('results/GSEA/FGSEA_Rbp4_1week_DOI.txt', show_col_types = FALSE)
Rbp4_GO <- rbind(Rbp4_GO, aux)
aux <- read_tsv('results/GSEA/FGSEA_Rbp4_1month_DOI.txt', show_col_types = FALSE)
Rbp4_GO <- rbind(Rbp4_GO, aux)
Rbp4_GO <- Rbp4_GO$pathway
Rbp4_GO <- unique(Rbp4_GO)
length(Rbp4_GO)

Hsiao <- rbind(PV,Cux2,Rbp4)
Hsiao <- unique(Hsiao)

x <- list(
  PV = PV$genes_all,
  Cux2 = Cux2$genes_all,
  Rbp4 = Rbp4$genes_all
)

x <- list(
  PV = PV,
  Cux2 = Cux2ERT,
  Rbp4 = Rbp4
)

plot_venn(x, title = "PV vs Cux2 vs Rbp4 gene overlap")
plot_venn_v2(x, "PV vs Cux2 vs Rbp4 FC > 0.2 genes overlap")
VennDiag <- euler(x, shape = "ellipse")
p1 <- plot(VennDiag, 
           counts = TRUE, 
           font=1, 
           cex=1, 
           alpha=0.6,
           main=title,
           fill=c("#a94498", "#45aa9a","#0f7733"),
           quantities = TRUE)
p1
export_svg(p1, dir = "results/GSEA/", file_name = "VennDiagram Genes", h = 8, w = 8)

x3 <- list(
  PV = PV_GO,
  Cux2 = Cux2_GO,
  Rbp4 = Rbp4_GO
)

intersect(intersect(PV_GO,Cux2_GO),Rbp4_GO)
plot_venn(x3, title = "PV vs Cux2 vs Rbp4 GO terms overlap")
plot_venn_v2(x3, "PV vs Cux2 vs Rbp4 GO terms overlap")
VennDiag <- euler(x3, shape = "ellipse")
p1 <- plot(VennDiag, 
           counts = TRUE, 
           font=1, 
           cex=1, 
           alpha=0.6,
           main=title,
           fill=c("#a94498", "#45aa9a","#0f7733"),
           quantities = TRUE)
p1
export_svg(p1, dir = "results/DEG/", file_name = "VennDiagram GO terms", h = 8, w = 8)

##############################################################################################
##### supp fig 3 comparing GSEA genes w/ AS genes ####

cell_type = c("Rbp4") # [PV Cux2ERT Rbp4]
timepoint = c("48hrs") # [48hrs 1week 1month]
treat = "Psilocybin" # [DOI Psilocybin Ctrl]

#get GSEA genes
aux = read_tsv(paste0("results/GSEA/FGSEA_", cell_type, "_", timepoint, "_", treat, ".txt"), show_col_types = FALSE)
aux$leadingEdge <- strsplit(aux$leadingEdge," ")
GSEA_genes = unique(unlist(aux$leadingEdge))
length(GSEA_genes)

#or get GSEA FC genes
PV <- read_table(paste0("results/GSEA/", cell_type, timepoint, treat, " FC Genes.txt"), show_col_types = F)
dim(PV)
PV <- unique(unlist(PV$current_genes))
length(PV)

#get AS genes
loc_res = read.delim("files/AS_results_location.csv", sep = ",", comment.char = "#")
head(loc_res)
table(loc_res$CellType)

AS_genes = get_combine_gene_events(loc_res, celltype = cell_type, treat = treat, time=timepoint)

x <- list(
  AS = AS_genes,
  GSEA = PV
)

p1 <- plot_venn_v2(x, paste0(cell_type, "_", timepoint, "_", treat))
p1
export_svg(p1, dir = "figures/drafts/rMATS/", file_name = paste0(cell_type, "_", timepoint, "_", treat, "_ASvsGSEA_genes"), h = 4, w = 4)

##############################################################################################
##### fig 3C tile plot of AS GO terms ####

AS_GO1 <- read_tsv(paste0("results/rmats/Episode_VI/GO/GO_Hsiao_PV_48hrs_DOI_AS_genes_combined.txt"), show_col_types = F)
AS_GO1$Condition <- "PV_48hrs_DOI"
AS_GO2 <- read_tsv(paste0("results/rmats/Episode_VI/GO/GO_Hsiao_PV_48hrs_Psilocybin_AS_genes_combined.txt"), show_col_types = F)
AS_GO2$Condition <- "PV_48hrs_Psilo"
AS_GO3 <- read_tsv(paste0("results/rmats/Episode_VI/GO/GO_Hsiao_PV_1week_DOI_AS_genes_combined.txt"), show_col_types = F)
AS_GO3$Condition <- "PV_1week_DOI"
AS_GO4 <- read_tsv(paste0("results/rmats/Episode_VI/GO/GO_Hsiao_PV_1month_DOI_AS_genes_combined.txt"), show_col_types = F)
AS_GO4$Condition <- "PV_1month_DOI"
AS_GO5 <- read_tsv(paste0("results/rmats/Episode_VI/GO/GO_Hsiao_Cux2ERT_48hrs_DOI_AS_genes_combined.txt"), show_col_types = F)
AS_GO5$Condition <- "Cux2ERT_48hrs_DOI"
AS_GO6 <- read_tsv(paste0("results/rmats/Episode_VI/GO/GO_Hsiao_Cux2ERT_48hrs_Psilocybin_AS_genes_combined.txt"), show_col_types = F)
AS_GO6$Condition <- "Cux2ERT_48hrs_Psilo"
AS_GO7 <- read_tsv(paste0("results/rmats/Episode_VI/GO/GO_Hsiao_Cux2ERT_1week_DOI_AS_genes_combined.txt"), show_col_types = F)
AS_GO7$Condition <- "Cux2ERT_1week_DOI"
AS_GO8 <- read_tsv(paste0("results/rmats/Episode_VI/GO/GO_Hsiao_Cux2ERT_1month_combat_combined_DOI_AS_genes_combined.txt"), show_col_types = F)
AS_GO8$Condition <- "Cux2ERT_1month_DOI"
AS_GO9 <- read_tsv(paste0("results/rmats/Episode_VI/GO/GO_Hsiao_Rbp4_48hrs_Psilocybin_AS_genes_combined.txt"), show_col_types = F)
AS_GO9 <- tail(AS_GO9,1)
AS_GO9$Description <- "regulation of synapse structure or activity"
AS_GO9$FoldEnrichment <- 0
AS_GO9$Condition <- "Rbp4_48hrs_DOI"
AS_GO10 <- read_tsv(paste0("results/rmats/Episode_VI/GO/GO_Hsiao_Rbp4_48hrs_Psilocybin_AS_genes_combined.txt"), show_col_types = F)
AS_GO10$Condition <- "Rbp4_48hrs_Psilo"
AS_GO11 <- read_tsv(paste0("results/rmats/Episode_VI/GO/GO_Hsiao_Rbp4_1week_combat_combined_DOI_AS_genes_combined.txt"), show_col_types = F)
AS_GO11$Condition <- "Rbp4_1week_DOI"
AS_GO12 <- read_tsv(paste0("results/rmats/Episode_VI/GO/GO_Hsiao_Rbp4_1month_combat_combined_DOI_AS_genes_combined.txt"), show_col_types = F)
AS_GO12$Condition <- "Rbp4_1month_DOI"
allGO <- rbind(AS_GO1,AS_GO2,AS_GO3,AS_GO4,AS_GO5,AS_GO6,AS_GO7,AS_GO8,AS_GO9,AS_GO10,AS_GO11,AS_GO12)
dim(allGO)

GO_to_plot <- allGO %>% 
  filter(-log10(allGO$p.adjust) > 6)
dim(GO_to_plot)
GO_to_plot <- unique(GO_to_plot$Description)
length(GO_to_plot)
GO_to_plot <- allGO[allGO$Description %in% GO_to_plot,]
dim(GO_to_plot)

order <- as.data.frame(table(GO_to_plot$Description))
order <- arrange(order, Freq)

p2 <- ggplot(GO_to_plot, aes(x = factor(Condition, level = c('PV_48hrs_Psilo', 'PV_48hrs_DOI', 'PV_1week_DOI','PV_1month_DOI','Cux2ERT_48hrs_Psilo','Cux2ERT_48hrs_DOI','Cux2ERT_1week_DOI','Cux2ERT_1month_DOI','Rbp4_48hrs_Psilo','Rbp4_48hrs_DOI','Rbp4_1week_DOI','Rbp4_1month_DOI')),
                             y = factor(Description, level = order$Var1), fill = p.adjust)) +
  geom_tile(color = "black",
            linetype = 1) +
  scale_fill_gradient(low = "white", 
                      high = "#757575",
                      trans = 'reverse') +
  theme_classic2() +
  theme(axis.text.y=element_text(colour="black"),
        axis.text.x=element_text(colour="black",angle = 90)) +
  labs(
    x = "Treatment",
    y = "",
    fill = "padj"
  )
p2

export_svg(p2, dir = "figures/drafts/rMATS/", file_name = "tile plot AS GO terms", h = 8, w = 8)

##############################################################################################
##### fig 3D Venn Diagrams #####
loc_res = read.delim("files/AS_results_location.csv", sep = ",", comment.char = "#")

x <- list(
  PV = get_combine_gene_events(loc_res, celltype = "PV", treat = "DOI",time="48hrs"),
  L2_3 = get_combine_gene_events(loc_res, celltype = "Cux2ERT", treat = "DOI",time="48hrs"),
  L5 = get_combine_gene_events(loc_res, celltype = "Rbp4", treat = "DOI",time="48hrs")
)

x <- list(
  PV = get_combine_gene_events(loc_res, celltype = "PV", treat = "Psilocybin",time="48hrs"),
  L2_3 = get_combine_gene_events(loc_res, celltype = "Cux2ERT", treat = "Psilocybin",time="48hrs"),
  L5 = get_combine_gene_events(loc_res, celltype = "Rbp4", treat = "Psilocybin",time="48hrs")
)

x <- list(
  PV = get_combine_gene_events(loc_res, celltype = "PV", treat = "DOI",time="1week"),
  L2_3 = get_combine_gene_events(loc_res, celltype = "Cux2ERT", treat = "DOI",time="1week"),
  L5 = get_combine_gene_events(loc_res, celltype = "Rbp4", treat = "DOI",time="1week", batch = TRUE)
)

x <- list(
  PV = get_combine_gene_events(loc_res, celltype = "PV", treat = "DOI",time="1month"),
  L2_3 = get_combine_gene_events(loc_res, celltype = "Cux2ERT", treat = "DOI",time="1month", batch = TRUE),
  L5 = get_combine_gene_events(loc_res, celltype = "Rbp4", treat = "DOI",time="1month", batch = TRUE)
)

CellType = c(PV = "#a94498", L2_3 = "#45aa9a", L5 = "#0f7733")
plot_venn_v2 <- function(x, title) {
  VennDiag <- euler(x, shape = "ellipse")
  p1 <- plot(VennDiag, 
             counts = TRUE, 
             font=1, 
             cex=2, 
             alpha=0.3,
             main=title,
             fill=CellType,
             quantities = TRUE)
  return(p1)
}
p1 <- plot_venn_v2(x, title = "DOI AS genes overlap")
p1
export_svg(p1, dir = "figures/drafts/rMATS/", file_name = "Venn Diagram AS genes", h = 8, w = 8)

##############################################################################################
##### fig 4A ECM genes ####
loc_res = read.delim("files/DESeq2_results_location.csv", sep = ",", comment.char = "#")
loc_res

ECM_genes = read.delim("files/ECM GENE LIST ORDER.txt", sep = " ", header = T)
head(ECM_genes)
dim(ECM_genes)
order_vector <- c('CSPG','link_protein','ECM_glycoprotein','hyaluronan_ligands','CS-GAG_chains_regulators','CS-GAG_chains_ligands','Regulators_of_CSPG_signaling','perisynaptic_matrix_regulation','ECM_components','ECM_regulators')

ECM_genes <- ECM_genes %>%
  group_by(PNN_components) %>%
  mutate(order = factor(PNN_components, levels = order_vector)) %>%
  arrange(order) %>%
  ungroup() %>%
  select(-order)

data_FC = data.frame(row.names = ECM_genes$Gene_name)
data_FC = get_log2FC_values(data_FC, loc_res)
head(data_FC)

table(complete.cases(data_FC))
data_FC_ = data_FC[complete.cases(data_FC),]
ECM_genes <- ECM_genes[complete.cases(data_FC),]

meta_table = data.frame(row.names = c('PV_48hrs_DOI', 'PV_48hrs_Psilocybin', 'PV_1week_DOI', 'PV_1month_DOI', 'Cux2ERT_48hrs_DOI', 'Cux2ERT_48hrs_Psilocybin', 'Cux2ERT_1week_DOI', 'Cux2ERT_1month_DOI', 'Rbp4_48hrs_DOI', 'Rbp4_48hrs_Psilocybin', 'Rbp4_1week_DOI', 'Rbp4_1month_DOI'))
meta_table = cbind(meta_table, Timepoint=c('48hrs', '48hrs', '1week', '1month', '48hrs', '48hrs','1week', '1month', '48hrs','48hrs',  '1week', '1month'))
meta_table = cbind(meta_table, CellType=c('PV', 'PV', 'PV', 'PV', 'Cux2ERT', 'Cux2ERT', 'Cux2ERT', 'Cux2ERT','Rbp4', 'Rbp4','Rbp4', 'Rbp4'))
meta_table = cbind(meta_table, Treatment=c('DOI', 'Psilocybin', 'DOI', 'DOI', 'DOI', 'Psilocybin','DOI', 'DOI', 'DOI','Psilocybin', 'DOI', 'DOI'))

meta_table = meta_table[rownames(meta_table) %in% colnames(data_FC),]
meta_table = meta_table[colnames(data_FC),]

breaksList = c(seq(-2, -1, by = 0.5), seq(-1, 1, 0.1), seq(1, 2, 0.5))
breaksList = unique(breaksList)

my_ramp = colorRampPalette(c("blue", "white", "orange"))(length(breaksList))
#my_ramp[4] = "#FFFFFF"

p1 = pheatmap(data_FC_,
              cluster_rows = F, 
              cluster_cols = F,
              show_rownames = T,
              show_colnames = F,
              annotation_col = meta_table,
              annotation_colors = figure_colors,
              annotation_row = ECM_genes[, 2, drop=FALSE],
              breaks = breaksList,
              color=my_ramp,
              border_color = NA
)
p1

export_svg(p1, dir = "results/GSEA/", file_name = "log2FC_Heatmap_ECM_Genes", h = 8, w = 8)

######################################################################################################
##### fig 5A logFC heatmap of pathway specific genes ####
loc_res = read.delim("files/DESeq2_results_location.csv", sep = ",", comment.char = "#")
loc_res

res <- read_tsv(paste0("results/GSEA/FGSEA_PV_48hrs_Psilocybin.txt"), show_col_types = FALSE)
res$leadingEdge <- strsplit(res$leadingEdge," ")
pathway_interest <- res[which(res$pathway=='GOBP_POTASSIUM_ION_TRANSPORT'),]
genes <- as.data.frame(unlist(pathway_interest$leadingEdge))
genes$GO <- 'GOBP_POTASSIUM_ION_TRANSPORT'
pathway_genes <- genes
pathway_interest <- res[which(res$pathway=='GOBP_SODIUM_ION_TRANSPORT'),]
genes <- as.data.frame(unlist(pathway_interest$leadingEdge))
genes$GO <- 'GOBP_SODIUM_ION_TRANSPORT'
pathway_genes <- rbind(pathway_genes,genes)
pathway_interest <- res[which(res$pathway=='GOBP_REGULATION_OF_POSTSYNAPSE_ORGANIZATION'),]
genes <- as.data.frame(unlist(pathway_interest$leadingEdge))
genes$GO <- 'GOBP_REGULATION_OF_POSTSYNAPSE_ORGANIZATION'
pathway_genes <- rbind(pathway_genes,genes)
pathway_interest <- res[which(res$pathway=='GOBP_SYNAPTIC_TRANSMISSION_GLUTAMATERGIC'),]
genes <- as.data.frame(unlist(pathway_interest$leadingEdge))
genes$GO <- 'GOBP_SYNAPTIC_TRANSMISSION_GLUTAMATERGIC'
pathway_genes <- rbind(pathway_genes,genes)
colnames(pathway_genes) <- c('gene_name','GO')
dim(pathway_genes)

data <- read_table("results/GSEA/PV48hrsPsilocybin FC Genes.txt", show_col_types = F)
pathway_genes_filt <- intersect(pathway_genes$gene_name, data$current_genes)
pathway_genes = pathway_genes[match(pathway_genes_filt, pathway_genes$gene_name),]

data_FC = data.frame(row.names = pathway_genes_filt)
data_FC = get_log2FC_values(data_FC, loc_res)
head(data_FC)

table(complete.cases(data_FC))
data_FC_ = data_FC[complete.cases(data_FC),]
pathway_genes = pathway_genes[complete.cases(data_FC),]

meta_table = data.frame(row.names = c('PV_48hrs_DOI', 'PV_48hrs_Psilocybin', 'PV_1week_DOI', 'PV_1month_DOI', 'Cux2ERT_48hrs_DOI', 'Cux2ERT_48hrs_Psilocybin', 'Cux2ERT_1week_DOI', 'Cux2ERT_1month_DOI', 'Rbp4_48hrs_DOI', 'Rbp4_48hrs_Psilocybin', 'Rbp4_1week_DOI', 'Rbp4_1month_DOI'))
meta_table = cbind(meta_table, Timepoint=c('48hrs', '48hrs', '1week', '1month', '48hrs', '48hrs','1week', '1month', '48hrs','48hrs',  '1week', '1month'))
meta_table = cbind(meta_table, CellType=c('PV', 'PV', 'PV', 'PV', 'Cux2ERT', 'Cux2ERT', 'Cux2ERT', 'Cux2ERT','Rbp4', 'Rbp4','Rbp4', 'Rbp4'))
meta_table = cbind(meta_table, Treatment=c('DOI', 'Psilocybin', 'DOI', 'DOI', 'DOI', 'Psilocybin','DOI', 'DOI', 'DOI','Psilocybin', 'DOI', 'DOI'))

meta_table = meta_table[rownames(meta_table) %in% colnames(data_FC),]
meta_table = meta_table[colnames(data_FC),]

breaksList = c(seq(-2, -1, by = 0.5), seq(-1, 1, 0.1), seq(1, 2, 0.5))
breaksList = unique(breaksList)

my_ramp = colorRampPalette(c("blue", "white", "orange"))(length(breaksList))
#my_ramp[4] = "#FFFFFF"

p1 = pheatmap(data_FC_,
              cluster_rows = F, 
              cluster_cols = F,
              show_rownames = T,
              show_colnames = F,
              annotation_col = meta_table,
              annotation_colors = figure_colors,
              breaks = breaksList,
              color=my_ramp,
              border_color = NA
)
p1

export_svg(p1, dir = "results/GSEA/", file_name = "log2FC_Heatmap_pathway_Genes", h = 8, w = 8)
##### fig 5A specific GO terms for GSEA and AS ####
data1 = read_tsv('results/GSEA/FGSEA_PV_48hrs_Psilocybin_all.txt', show_col_types = FALSE)
GO_terms = read_xlsx('files/Fig5_GO_reordered2.xlsx',sheet = 1,col_names = FALSE)
dim(GO_terms)
GO_terms_GSEA <- data1[which(data1$pathway %in% GO_terms$...1),]
dim(GO_terms_GSEA)

order_vector <- GO_terms$...1

# Reorder the rows based on the group and the order vector
GO_terms_GSEA <- GO_terms_GSEA %>%
  group_by(pathway) %>%
  mutate(order = factor(pathway, levels = order_vector)) %>%
  arrange(order) %>%
  ungroup() %>%
  select(-order)

p <- ggplot(GO_terms_GSEA, aes(factor(pathway, levels = rev(order_vector)), NES)) +
  geom_col(aes(fill=padj < .01),color = "black") +
  scale_fill_manual(values=c("#f1eef6","#a8ddb5")) +
  geom_text(aes(label=round(padj,digits = 3)), position=position_dodge(width=0.9), vjust=.5, hjust =-2) +
  ylim(-2.5,2.5) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GO Normalized Enrichment Score from GSEA") + 
  theme_minimal()
p

export_svg(p, dir = "results/GSEA/", file_name = "GSEA Fig5 GO terms REORDERED2", h = 8, w = 10)

data2 = read_tsv(paste0("results/rmats/Episode_VI/GO/PV_48hrs_Psilo_qval_1.txt"), show_col_types = F)
GO_terms = read_xlsx('files/Fig5_GO_reordered2.xlsx',sheet = 2,col_names = FALSE)
dim(GO_terms)
GO_terms_AS <- data2[which(data2$Description %in% GO_terms$...1),]
dim(GO_terms_AS)

order_vector <- GO_terms$...1

GO_terms_AS <- GO_terms_AS %>%
  group_by(Description) %>%
  mutate(order = factor(Description, levels = order_vector)) %>%
  arrange(order) %>%
  ungroup() %>%
  select(-order)

p1 <- ggplot(GO_terms_AS, aes(factor(Description, levels = rev(order_vector)), x = Count, fill = p.adjust < .01)) +
  geom_bar(stat="identity", width=0.9, color = "black") +
  scale_fill_manual(values=c("#f1eef6","#a8ddb5")) +
  guides(color = guide_colorbar(reverse = TRUE)) +
  geom_text(aes(label=round(p.adjust,digits = 3)), position=position_dodge(width=0.9), vjust=.5, hjust = -.5) +
  theme_bw(base_size = 15) +
  theme(axis.text.y=element_text(colour="black"),
        axis.text.x=element_text(colour="black")) +
  labs(
    x = "Count",
    y = ""
    #size = "Enrichment"
  )
p1

export_svg(p1, dir = "figures/drafts/rMATS/", file_name = "PV48hrsPsilo Fig5 GO terms AS REORDERED2", h = 8, w = 10)

########################################################################################## ECM Genes Expression
#### FC
# 
# ECM_genes = read.delim("files/ECM GENE LIST ORDER.txt", sep = " ", header = T)
# head(ECM_genes)
# 
# ECM_genes$Gene_name <- str_to_title(ECM_genes$Gene_name) 
# rownames(ECM_genes) <- ECM_genes$Gene_name
# 
# ### check which are DE and add annotation
# data_FC = data.frame(row.names = ECM_genes$Gene_name)
# data = read_file('results/DEG/Episode III/Timepoint/genes/48hrs/PV/pipe0/DEG_DOIvsSaline_48hrs_only_EdgeR.txt')
# head(data)
# table(data$Expression)
# data_FC = cbind(data_FC, pv_48h=data$logFC[match(ECM_genes$Gene_name, data$Genes)])
# 
# data = read_file('results/DEG/Episode III/Timepoint/genes/48hrs/PV/pipe0/DEG_PsilovsSaline_48hrs_only_EdgeR.txt')
# head(data)
# table(data$Expression)
# data_FC = cbind(data_FC, pv_48h_psilo=data$logFC[match(ECM_genes$Gene_name, data$Genes)])
# 
# 
# data = read_file('results/DEG/Episode III/Timepoint/genes/1week/PV/pipe0/DEG_DOIvsSaline_1week_only_EdgeR.txt')
# head(data)
# table(data$Expression)
# data_FC = cbind(data_FC, pv_1week=data$logFC[match(ECM_genes$Gene_name, data$Genes)])
# 
# data = read_file('results/DEG/Episode III/Timepoint/genes/1month//PV/pipe0/DEG_DOIvsSaline_1month_only_EdgeR.txt')
# head(data)
# table(data$Expression)
# data_FC = cbind(data_FC, pv_1month=data$logFC[match(ECM_genes$Gene_name, data$Genes)])
# 
# ## Cux2 
# data = read_file('results/DEG/Episode III/Timepoint/genes/48hrs/Cux2ERT/pipe0/DEG_DOIvsSaline_48hrs_only_EdgeR.txt')
# head(data)
# table(data$Expression)
# data_FC = cbind(data_FC, cux2_48h=data$logFC[match(ECM_genes$Gene_name, data$Genes)])
# 
# data = read_file('results/DEG/Episode III/Timepoint/genes/1week/Cux2/pipe0/DEG_DOIvsSaline_1week_only_EdgeR.txt')
# head(data)
# table(data$Expression)
# data_FC = cbind(data_FC, cux2_1week=data$logFC[match(ECM_genes$Gene_name, data$Genes)])
# 
# ## Rbp4
# data = read_file('results/DEG/Episode III/Timepoint/genes/48hrs/Rbp4/pipe0/DEG_DOIvsSaline_48hrs_only_EdgeR.txt')
# head(data)
# table(data$Expression)
# data_FC = cbind(data_FC, rbp4_48h=data$logFC[match(ECM_genes$Gene_name, data$Genes)])
# 
# data = read_file('results/DEG/Episode III/Timepoint/genes/1week/Rbp4/pipe0/DEG_DOIvsSaline_1week_only_EdgeR.txt')
# head(data)
# table(data$Expression)
# data_FC = cbind(data_FC, rbp4_1week=data$logFC[match(ECM_genes$Gene_name, data$Genes)])
# 
# data = read_file('results/DEG/Episode III/Timepoint/genes/1month/Rbp4/pipe0/DEG_DOIvsSaline_1month_only_EdgeR.txt')
# head(data)
# table(data$Expression)
# data_FC = cbind(data_FC, rbp4_1month=data$logFC[match(ECM_genes$Gene_name, data$Genes)])
# head(data_FC)
# 
# complete.cases(data_FC)
# 
# breaksList = c(seq(-2, -1, by = 0.5), seq(-1, 1, 0.5), seq(1, 2, 0.5))
# #breaksList = seq(-2, 2, by = 0.5)
# breaksList = unique(breaksList)
# 
# 
# meta_table = data.frame(row.names = c('pv_48h', 'pv_48', 'pv_1week', 'pv_1month', 'cux2_48h', 'cux2_1week', 'rbp4_48h', 'rbp4_1week', 'rbp4_1month'))
# meta_table = cbind(meta_table, Timepoint=c('48hrs', '48hrs', '1week', '1month', '48hrs', '1week', '48hrs', '1week', '1month'))
# meta_table = cbind(meta_table, CellType=c('PV', 'PV', 'PV', 'PV', 'Cux2ERT', 'Cux2ERT', 'Rbp4', 'Rbp4', 'Rbp4'))
# meta_table = cbind(meta_table, Treatment=c('DOI', 'Psilocybin', 'DOI', 'DOI', 'DOI', 'DOI', 'DOI', 'DOI', 'DOI'))
# 
# my_ramp = colorRampPalette(c("#4575b4", "white", "#d73027"))(length(breaksList))
# my_ramp[4] = "#FFFFFF"
# 
# p1 = pheatmap(data_FC,
#               cluster_rows = F, 
#               cluster_cols = F,
#               show_rownames = T,
#               show_colnames = F,
#               annotation_col = meta_table,
#               annotation_colors = my_colour,
#               annotation_row = ECM_genes[, 2, drop=FALSE],
#               #cutree_rows = 3,
#               breaks = breaksList,
#               color=my_ramp,
#               #main=title,
#               border_color = F)
# p1
# export_svg(p1, dir = "figures/drafts/", file_name = "pnn_genes_heatmap", h = 9, w = 5)
# 
# ####
# 
# 
# ZlogCPM_filt
# samples_g1 = meta_table$Sample.ID[meta_table$Treatment == 'Ctrl']
# samples_g2 = meta_table$Sample.ID[meta_table$Treatment == 'DOI']
# avg_Zscore = cbind(ctrl=apply(ZlogCPM_filt[,samples_g1], 1, mean), doi=apply(ZlogCPM_filt[,samples_g2], 1, mean))
# head(avg_Zscore)
# 
# p1 = pheatmap(avg_Zscore,
#               cluster_rows = F, 
#               cluster_cols = F,
#               show_rownames = T,
#               #annotation_col = meta_table[, c(2:6)],
#               annotation_colors = figure_colors,
#               annotation_row = ECM_genes[, c(2, 3)],
#               cutree_rows = 3,
#               color=colorRampPalette(c("#4575b4", "white", "#d73027"))(50),
#               #main=title,
#               border_color = F)
# p1
# 
