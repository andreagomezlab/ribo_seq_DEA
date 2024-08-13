library('org.Mm.eg.db')
library(ReactomePA)
library(clusterProfiler)
library(enrichplot)


plot_table <- function(table) {
  grid.newpage()
  return(gridExtra::grid.table(table, rows = NULL))
}

plot_PCA = function(matrix, meta_table, var1="Treatment", var2="Sex", scale=FALSE) {
  pca_res <- prcomp(t(matrix), scale=scale)
  p1 <- autoplot(pca_res,
                 data = meta_table, 
                 colour=var1, 
                 shape=var2,
                 size=5) +
    theme_classic2(base_size = 17) +
    geom_text_repel(aes(x=PC1, y=PC2, label=Sample.ID), box.padding = 0.8)
  
  plot_var_explained(pca_res)
  
  return(p1)
}


plot_var_explained <- function(PC) {
  PC_sd <- setNames(PC$sdev , paste0("PC",1:length(PC$sdev)))
  PC_var_expl <- ( PC_sd^2 ) / sum(PC_sd^2) * 100
  
  # Which PCs explain at least 2% of variance?
  p1 = barplot(PC_var_expl, las=2, ylab="% variance explained") +
    abline(h=2, lty=2)
  
  return(p1)
}

plot_boxplot_topN_genes <- function (matrix, topn = 30, ylimits=c(0, 10), title="") {
  p1 = boxplot( t(matrix[1:topn,]),
                ylim=ylimits,
                las=2 ,
                col="grey" ,
                main=title,
                cex=.2)
  return(p1)
}

get_topn_var_genes <- function(matrix, topn=500) {
  vars <- apply(matrix,1,var)
  vars <- sort( vars , decreasing=T)
  top_var <- names( vars ) [1:topn]
  return(top_var)
}

plot_boxplot_topN_var_genes <- function(matrix, topn) {
  means <- rowMeans(matrix)
  vars <- apply(matrix,1,var)
  plot(x=means, y=vars, cex=.1 )
  
  vars <- sort( vars , decreasing=T)
  top_var <- names( vars ) [1:topn]
  #boxplot( t(matrix[ top_var[1:15],])   , ylim=c(0,12), las=2 , col="grey" , main="logCPM (top var genes)" ,cex=.2)
}

plot_sample_correlation_dend <- function(logcounts, meta_table) {
  #Compute sample correlations
  sample_cor <- cor( logcounts )
  round(sample_cor,4)
  
  #Transform the scale from correlations
  cor_distance <- -(sample_cor - 1)/2
  round(cor_distance,4)
  
  #Convert it to a distance object
  d2 <- as.dist(cor_distance)
  h2 <- hclust(d2, method="complete")
  p1 = plot( as.dendrogram(h2) , las=1, main="d=correlation\nh=complete")
  p1 + points(1:ncol(Counts_only) ,rep(0,ncol(Counts_only)), pch= 16, cex=2, col=as.factor(meta_table$Batch[h2$order]))
  return(p1)
}


plot_heatmap <- function(matrix, meta_table, top_var, title="") {
  my_colour = list(
    Treatment = c(DOI = "#5977ff", Ctrl = "#f74747"),
    CellType = c(Rbp4 = "#f6e8c3", PV = "#9e82ed", Cux2ERT = "#5ab4ac"),
    Sex = c(`F` = "#d73027", `M` = "#4575b4"),
    Batch = c('B2' = "#1b7837", 'B3' = "#af8dc3"),
    Timepoint = c(`1week` = "#e89829", `1month` = "#82ed82", `48hrs` = "brown")
  )
  
  highly_variable_lcpm <- matrix[top_var, ]
  
  p1 = pheatmap(highly_variable_lcpm,
           cluster_rows = T, 
           cluster_cols = T,
           annotation_col = meta_table[, c(2:6)],
           annotation_colors = my_colour,
           gaps_col =  3,
           color=colorRampPalette(c("#4575b4", "white", "#d73027"))(50),
           main=title,
           border_color = F)
  return(p1)
}


plot_volcano <- function(data, top_genes) {
  
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
    scale_x_continuous(limits = c(-5, 5), breaks = c(-7, -5, -2, 0, 2, 5, 7))

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
  return(p2)
}

get_enrichPathway <- function(symbols, pcutoff=0.05, ncat=10) {
  
  gene <- mapIds(org.Mm.eg.db, symbols, 'ENTREZID', 'SYMBOL')
  gene = gene[!is.na(gene)]
  #print(paste0(length(gene), ' mapped genes'))
  yy = enrichPathway(gene, pvalueCutoff=pcutoff, organism = "mouse", readable = TRUE)
  print(head(as.data.frame(yy)))
  if (length(yy) > 0) {
    print(dotplot(yy, showCategory=ncat))
    return(as.data.frame(yy))
  }
  return(NULL)
}

get_enrichGO <- function(symbols, pcutoff=0.01, qcutoff=0.05, ncat=20, ont="BP") {
  
  gene <- mapIds(org.Mm.eg.db, symbols, 'ENTREZID', 'SYMBOL')
  gene = gene[!is.na(gene)]
  print(paste0(length(gene), ' mapped genes'))
  
  ego <- enrichGO(gene          = gene,
                  OrgDb         = org.Mm.eg.db,
                  ont           = ont,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = pcutoff,
                  qvalueCutoff  = qcutoff,
                  readable      = TRUE)
  head(ego)
  if (length(ego) > 0) {
    print(barplot(ego, showCategory=ncat))
    return(as.data.frame(ego))
  }
  return(NULL)
}
