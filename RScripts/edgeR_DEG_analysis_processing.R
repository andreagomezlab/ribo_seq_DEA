library(biomaRt)
library(edgeR)
library(DESeq2)
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
library("IHW")
library(tximportData)
library(tximport)


ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
bm <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'gene_biotype', 'ensembl_transcript_id'), mart = ensembl)
#save(bm, file = "data/rdata/bm.rds")
head(bm)


filter_genes <- function(y, min_count=10, large_n=10, min_prop=0.7) {
  keep <- filterByExpr(y, min.count=min_count, large.n=large_n, min.prop=min_prop)
  out <- bm[match(names(keep), bm$external_gene_name),]
  keep = keep & (out$gene_biotype=="protein_coding" & !is.na(out$gene_biotype))
  print(table(keep))
  y <- y[keep, ,keep.lib.sizes=FALSE]
  return(y)
}

get_DEG <- function(y, batch=FALSE) {
  if (batch){
    design <- model.matrix(~meta_table$Batch+meta_table$Treatment)
  } else {
    design <- model.matrix(~meta_table$Treatment)
  }
  
  rownames(design) <- colnames(y)
  y <- normLibSizes(y)
  y <- estimateDisp(y, design, robust = TRUE)
  #y <- estimateTagwiseDisp(y)
  
  print(plotBCV(y))
  
  et <- exactTest(y)
  # The classic edgeR pipeline: pairwise comparisons between two or more groups
  fit <- glmQLFit(y, design, robust = TRUE)
  #fit <- glmFit(y, design, robust = TRUE)
  
  #plotQLDisp(fit)
  
  res <- glmQLFTest(fit)
  #res <- glmLRT(fit)
  #topTags(res)
  
  print(topTags(et, adjust.method = "BH", n = 20))
  #topTags(et, adjust.method = "BH", n=nrow(et$table))
  data = topTags(et, adjust.method = "BH", n=nrow(et$table), p.value = 1) #bonferroni
  data = as.data.frame(data)
  #colnames(data)[1] <- c("Genes")
  #data$Genes = rownames(data)
  return(as.data.frame(data))
}

set_DEG_status <- function(data, p_fc, p_pval) {
  data <- data %>% 
    mutate(
      Expression = case_when(logFC >= p_fc & PValue <= p_pval ~ "Up-regulated",
                             logFC <= -p_fc & PValue <= p_pval ~ "Down-regulated",
                             TRUE ~ "Unchanged")
    )
  return(data)
}

get_topN_DEG <- function(data, top=15) {
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
  return(top_genes)
}

get_topN_DEG_v2 <- function(data, top=15) {
  top_genes <- bind_rows(
    data %>% 
      filter(Expression == 'Up-regulated') %>% 
      arrange(padj, desc(abs(log2FoldChange))) %>% 
      head(top),
    data %>% 
      filter(Expression == 'Down-regulated') %>% 
      arrange(padj, desc(abs(log2FoldChange))) %>% 
      head(top)
  )
  return(top_genes)
}

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

overlap_coef <- function(a, b){
  intersection = length(intersect(a, b))
  den = min(length(a), length(b))
  return (intersection/den)
}
