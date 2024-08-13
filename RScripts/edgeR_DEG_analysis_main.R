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

# Load scripts
source("src/edgeR_DEG_analysis_IO.R")
source("src/edgeR_DEG_analysis_visualization.R")
source("src/edgeR_DEG_analysis_processing.R")

### select cell type
cell_type = "PV" # [PV Rbp4 Cux2ERT]
timepoint = "48hrs" # [48hrs 1week 1month]

meta_table = get_metatable(cell_type, timepoint)
plot_table(meta_table)

##
y = get_DGEList_counts(meta_table)
y = filter_genes(y)

logcounts <- cpm(y, normalized.lib.sizes = TRUE, log=TRUE)
head(logcounts)
plot_boxplot_topN_genes(logcounts, 30, ylimits=c(-5, 12), "logCPM raw_data")

ZlogCPM <- t(scale(t(logcounts)))
head(ZlogCPM)
plot_boxplot_topN_genes(ZlogCPM, 30, ylimits=c(-4, 4), "Z-score data")

##
counts_vst = vst(round(y$counts))
plot_PCA(counts_vst, meta_table, "Treatment", "Sex", TRUE)
plot_PCA(logcounts, meta_table, "Treatment", "Sex", TRUE)

top_var <- get_topn_var_genes(logcounts, 500)
plot_PCA(counts_vst[ top_var, ], meta_table, "Treatment", "Sex", TRUE)

plot_heatmap(ZlogCPM, meta_table, top_var, title="Top 500 variable genes across conditions")

title_condition = "DOI vs Saline - 48 hours"

data = get_DEG(y)

data = set_DEG_status(data, p_fc = 0, p_pval = 0.05)
head(data) 
table(data$Expression)

save_file(data, 
          dir = "results/DEG/Episode III/Timepoint/genes/PV/", 
          file_name = "DEG_DOIvsSaline_48hrs_only_EdgeR.txt",
          row = FALSE,
          sep = "\t")

top_genes = get_topN_DEG(data, 15)

plot_volcano(data, top_genes)

plot_heatmap(ZlogCPM, 
             meta_table, 
             top_genes$Genes, 
             title=paste0("Top 15 Up and Down-regulated DEGs \n", condition, " (z-score expression)")
)

## UP Genes
symbols = data$Genes[data$Expression == 'Up-regulated']
# number of terms to be plot ncat
resp = get_enrichPathway(symbols, pcutoff=0.05, ncat=10)
save_file(resp, 
          dir = "results/DEG/Episode III/Timepoint/genes/PV/",
          file_name = paste("DEG_DOIvsSaline", cell_type, timepoint, "only_EdgeR_Reactome_analysis_UP_genes.txt", sep = "_"),
          row = FALSE,
          sep = "\t")
# number of categories to plot (ncat), Biological Process GO terms (ont)
ego = get_enrichGO(symbols, pcutoff = 0.01, qcutoff = 0.05, ncat = 20, ont = "BP")
save_file(ego, 
          dir = "results/DEG/Episode III/Timepoint/genes/PV/",
          file_name = paste("DEG_DOIvsSaline", cell_type, timepoint, "only_EdgeR_GO_analysis_UP_genes.txt", sep = "_"),
          row = FALSE,
          sep = "\t")

## DOWN Genes
symbols = data$Genes[data$Expression == 'Down-regulated']
resp = get_enrichPathway(symbols, pcutoff=0.05, ncat=10)
save_file(resp, 
          dir = "results/DEG/Episode III/Timepoint/genes/PV/",
          file_name = paste("DEG_DOIvsSaline", cell_type, timepoint, "only_EdgeR_Reactome_analysis_DOWN_genes.txt", sep = "_"),
          row = FALSE,
          sep = "\t")

ego = get_enrichGO(symbols, pcutoff = 0.01, qcutoff = 0.05, ncat = 20, ont = "BP")
save_file(ego, 
          dir = "results/DEG/Episode III/Timepoint/genes/PV/",
          file_name = paste("DEG_DOIvsSaline", cell_type, timepoint, "only_EdgeR_GO_analysis_DOWN_genes.txt", sep = "_"),
          row = FALSE,
          sep = "\t")


