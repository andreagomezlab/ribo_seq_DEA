library('org.Mm.eg.db')
library(ReactomePA)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(grid)
library(ggVennDiagram)
library("VennDiagram")
library(eulerr)

figure_colors = list (
  Treatment = c(DOI = "#e54291", Ctrl = "#525252", Psilocybin="#1f87c5"),
  CellType = c(Rbp4 = "#0f7733", PV = "#a94498", Cux2ERT = "#45aa9a"),
  Sex = c('F' = "#d9d9d9", 'M' = "#525252"),
  Batch = c('B1' = "#dfc27d", 'B2' = "#c7eae5", 'B3' = "#bac4e3", 'B4' = "#b3be62"),
  Timepoint = c("48hrs" = "#80afd2", "1week" = "#0072b3", "1month" = "#313695"),
  Status = c("Down-regulated" = "orange", "Unchanged" = "#f7f7f7", "Up-regulated" = "blue")
)

plot_cell_prop <- function(prop.cells) {
  prop.cells <- prop.cells[order(prop.cells$Freq),]
  prop.cells$Var1 <- factor(prop.cells$Var1, levels = prop.cells$Var1)
  
  p1 <- ggplot(prop.cells, aes(Var1, Freq)) +
    geom_col() +
    theme_minimal(base_size = 15) +
    geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=-0.25) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title="Number of cells in each cell type", x ="", y = "Number of Cells")
  return(p1)
}


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


plot_heatmap <- function(matrix, meta_table, row_meta=NULL, top_var, title="") {
  my_colour = list(
    Treatment = c(DOI = "#5977ff", Ctrl = "#f74747", Psilocybin="#f5945c"),
    CellType = c(Rbp4 = "#0f7733", PV = "#a94498", Cux2ERT = "#45aa9a"),
    Sex = c('F' = "#f1b6da", 'M' = "#8073ac"),
    Batch = c('B1' = "#dfc27d", 'B2' = "#c7eae5", 'B3' = "#bac4e3", 'B4' = "#b3be62"),
    Timepoint = c("48hrs" = "#80afd2", "1week" = "#0072b3", "1month" = "#313695"),
    ### row annotation colors
    ## Functions
    Essential.Genes = c('Yes' = "black", "No"="white"),
    Splicing.regulation = c('Yes' = "black", "No"="white"),
    Spliceosome = c('Yes' = "black", "No"="white"),
    RNA.modification = c('Yes' = "black", "No"="white"),
    X3..end.processing = c('Yes' = "black", "No"="white"),
    rRNA.processing = c('Yes' = "black", "No"="white"),
    Ribosome...basic.translation = c('Yes' = "black", "No"="white"),
    RNA.stability...decay = c('Yes' = "black", "No"="white"),
    microRNA.processing = c('Yes' = "black", "No"="white"),
    RNA.localization = c('Yes' = "black", "No"="white"),
    RNA.export = c('Yes' = "black", "No"="white"),
    Translation.regulation = c('Yes' = "black", "No"="white"),
    tRNA.regulation = c('Yes' = "black", "No"="white"),
    mitochondrial.RNA.regulation = c('Yes' = "black", "No"="white"),
    Viral.RNA.regulation = c('Yes' = "black", "No"="white"),
    snoRNA...snRNA...telomerase = c('Yes' = "black", "No"="white"),
    P.body...stress.granules = c('Yes' = "black", "No"="white"),
    Exon.Junction.Complex = c('Yes' = "black", "No"="white"),
    Novel.RBP = c('Yes' = "black", "No"="white"),
    Other = c('Yes' = "black", "No"="white"),
    ## Domains
    RRM = c('Yes' = "black", "No"="white"),
    ZNF = c('Yes' = "black", "No"="white"),
    KH = c('Yes' = "black", "No"="white"),
    Helicase = c('Yes' = "black", "No"="white"),
    Nuclease = c('Yes' = "black", "No"="white"),
    dRBM = c('Yes' = "black", "No"="white"),
    PUM_HD = c('Yes' = "black", "No"="white")
  )
  
  highly_variable_lcpm <- matrix[rownames(matrix) %in% top_var, ]
  highly_variable_lcpm <- highly_variable_lcpm[top_var, ]
  
  if (!is.null(row_meta)) {
    p1 = pheatmap(highly_variable_lcpm,
                  cluster_rows = T, 
                  cluster_cols = T,
                  show_rownames = T,
                  annotation_col = meta_table[, c(2:6)],
                  annotation_colors = my_colour,
                  annotation_row = row_meta,
                  #gaps_col =  3,
                  cutree_rows = 3,
                  #cutree_cols = 3,
                  color=colorRampPalette(c("#4575b4", "white", "#d73027"))(50),
                  main=title,
                  border_color = F)
  } else {
    p1 = pheatmap(highly_variable_lcpm,
                  cluster_rows = F, 
                  cluster_cols = F,
                  show_rownames = T,
                  annotation_col = meta_table[, c(2:6)],
                  annotation_colors = my_colour,
                  #gaps_col =  3,
                  #cutree_rows = 3,
                  #cutree_cols = 3,
                  color=colorRampPalette(c("#4575b4", "white", "#d73027"))(50),
                  main=title,
                  border_color = F)
  }
  
  return(p1)
}

plot_heatmap_simple <- function(matrix, meta_table, row_meta=NULL, top_var, title="") {
  my_colour = list(
    Treatment = c(DOI = "#e54291", Ctrl = "#525252", Psilocybin="#1f87c5"),
    CellType = c(Rbp4 = "#0f7733", PV = "#a94498", Cux2ERT = "#45aa9a"),
    Sex = c('F' = "#d9d9d9", 'M' = "#525252"),
    Batch = c('B1' = "#dfc27d", 'B2' = "#c7eae5", 'B3' = "#bac4e3", 'B4' = "#b3be62"),
    Timepoint = c("48hrs" = "#80afd2", "1week" = "#0072b3", "1month" = "#313695")
  )
  
  highly_variable_lcpm <- matrix[rownames(matrix) %in% top_var, ]
  head(highly_variable_lcpm)
  
  if (!is.null(row_meta)) {
    p1 = pheatmap(highly_variable_lcpm,
                  cluster_rows = F, 
                  cluster_cols = F,
                  show_rownames = T,
                  annotation_col = meta_table[, c(2:6)],
                  annotation_colors = my_colour,
                  annotation_row = row_meta,
                  cutree_rows = 3,
                  color=colorRampPalette(c("#4575b4", "white", "#d73027"))(50),
                  main=title,
                  border_color = F)
  } else {
    p1 = pheatmap(highly_variable_lcpm,
                  cluster_rows = T, 
                  cluster_cols = T,
                  show_rownames = F,
                  annotation_col = meta_table[, c(2:6)],
                  annotation_colors = my_colour,
                  color=colorRampPalette(c("#4575b4", "white", "#d73027"))(50),
                  main=title,
                  border_color = F)
  }
  
  return(p1)
}


plot_volcano <- function(data, top_genes, condition, p_pv, p_fc) {
  
  p2 <- ggplot(data, aes(logFC, -log(FDR,10))) +
    geom_point(aes(color = Expression), size = 3, alpha = 0.4,) +
    xlab(expression("log"[2]*"FC")) + 
    ylab(expression("-log"[10]*"FDR")) +
    scale_color_manual(values = c("orange", "#e0e0e0", "blue")) +
    guides(colour = guide_legend(override.aes = list(size=1.5))) +
    labs(title = paste0("Gene expression changes in ", condition)) +
    theme_classic(base_size = 22) +
    # Change line type and color
    geom_hline(yintercept=-log10(p_pv), linetype="dashed", color = "black") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red") +
    geom_vline(xintercept=-p_fc, linetype="dashed", color = "black") +
    geom_vline(xintercept=p_fc, linetype="dashed", color = "black")+
    theme(plot.title = element_text(hjust = 0.3)) +
    scale_x_continuous(limits = c(-4, 4), breaks = c(-4, -2, 0, 2, 4)) +
    scale_y_continuous(limits = c(0, 2))

  p2 <- p2 +   
    geom_label_repel(
      data          = subset(top_genes, Expression == "Down-regulated"),
      aes(label = top_genes$Genes[which(top_genes$Expression == "Down-regulated")]),
      segment.size  = 0.3,
      box.padding = 0.25,
      point.padding = 0.3,
      #force        = 0.5,
      #nudge_x      = -0.25,
      direction     = "y",
      size = 5,
      segment.color = "grey50",
      hjust         = 0.5
    ) +
    geom_label_repel(
      data          = subset(top_genes, Expression == "Up-regulated"),
      aes(label = top_genes$Genes[which(top_genes$Expression == "Up-regulated")]),
      segment.size  = 0.3,
      box.padding = 0.25,
      point.padding = 0.3,
      #nudge_y      =  6,
      #fill = "white",
      direction     = "y",
      size = 5,
      segment.color = "grey50",
      hjust         = -0.5
    )
  return(p2)
}



plot_volcano_v2 <- function(data, top_genes, condition, p_pv, p_fc) {
  
  p2 <- ggplot(data, aes(log2FoldChange, -log(padj,10))) +
    geom_point(aes(color = Expression), size = 3, alpha = 0.4,) +
    xlab(expression("log"[2]*"FC")) + 
    ylab(expression("-log"[10]*"padj")) +
    scale_color_manual(values = c('Down-regulated' = "orange", 'Unchanged' = "#e0e0e0", 'Up-regulated' = "blue")) +
    guides(colour = guide_legend(override.aes = list(size=1.5))) +
    labs(title = paste0("Gene expression changes in ", condition)) +
    theme_classic(base_size = 22) +
    # Change line type and color
    geom_hline(yintercept=-log10(p_pv), linetype="dashed", color = "black") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red") +
    geom_vline(xintercept=-p_fc, linetype="dashed", color = "black") +
    geom_vline(xintercept=p_fc, linetype="dashed", color = "black")+
    theme(plot.title = element_text(hjust = 0.3)) +
    scale_x_continuous(limits = c(-3, 3), breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
    scale_y_continuous(limits = c(0, 4.5))
  
  p2 <- p2 +   
    geom_label_repel(
      data          = subset(top_genes, Expression == "Down-regulated"),
      aes(label = top_genes$Genes[which(top_genes$Expression == "Down-regulated")]),
      segment.size  = 0.5,
      nudge_x = -0.35,
      nudge_y = 0.1,
      #segment.alpha = 0.5,
      box.padding = 0.25,
      point.padding = 0.3,
      direction     = "y",
      size = 5,
      segment.color = "grey50"
      #hjust         = 0.5
    ) +
    geom_label_repel(
      data          = subset(top_genes, Expression == "Up-regulated"),
      aes(label = top_genes$Genes[which(top_genes$Expression == "Up-regulated")]),
      segment.size  = 0.5,
      nudge_x = 0.35,
      nudge_y = 0.1,
      #segment.alpha = 0.5,
      box.padding = 0.25,
      point.padding = 0.3,
      direction     = "y",
      size = 5,
      segment.color = "grey50"
    )
  return(p2)
}


get_enrichPathway <- function(symbols, pcutoff=0.05, ncat=10) {
  
  gene <- mapIds(org.Mm.eg.db, symbols, 'ENTREZID', 'SYMBOL')
  gene = gene[!is.na(gene)]
  print(paste0(length(gene), ' mapped genes'))
  yy = enrichPathway(gene, pvalueCutoff=pcutoff, organism = "mouse", readable = TRUE)
  print(head(as.data.frame(yy)))
  if (length(yy) > 0) {
    print(dotplot(yy, showCategory=ncat))
    return(yy)
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
    return(ego)
  }
  return(NULL)
}

make_ratio_numeric = function(ratio) {
  ratio = as.numeric(unlist(str_split(ratio, "/"))[1])/as.numeric(unlist(str_split(ratio, "/"))[2])
  ratio
}

plot_enrich_dotplot <- function(enrich.df, ncat=10) {
  enrich.df = enrich.df[order(enrich.df$Count),]
  enrich.df <- cbind(enrich.df,order=1:nrow(enrich.df))
  enrich.df$GeneRatio_num =  sapply(enrich.df$GeneRatio, make_ratio_numeric)
  print(head(enrich.df))
  p1 <- ggplot(head(enrich.df, ncat), aes(x = GeneRatio_num, y = order, fill = p.adjust, size = Count)) +
    geom_point(shape = 21, stroke = 0.5) + # change the thickness of the boarder with stroke
    scale_fill_gradient(low = "white", 
                        high = "#317ebb",
                        trans = 'reverse') +
    ylab("") +
    xlab("") + 
    scale_y_continuous(breaks=1:nrow(enrich.df),labels=enrich.df$Description) +
    guides(color = guide_colorbar(reverse = TRUE)) +
    theme_bw(base_size = 15) +
    theme(axis.text.y=element_text(colour="black"),
          axis.text.x=element_text(colour="black")) +
    labs(
      x = "Gene Ratio",
      y = "",
      fill = "P-value adjusted",
      size = "Enrichment"
    )
  return(p1)
}

plot_enrich_barplot <- function(enrich.df, ncat=20) {
  
  top = head(enrich.df, ncat)
  top$Description <- factor(top$Description, levels = rev(top$Description))
  
  p1 <- ggplot(top, aes(x = Count, y = Description, fill = p.adjust)) +
    geom_bar(stat="identity", width=0.9) +
    scale_fill_gradient(low = "white", 
                        high = "#317ebb",
                        trans = 'reverse') +
    guides(color = guide_colorbar(reverse = TRUE)) +
    theme_bw(base_size = 15) +
    theme(axis.text.y=element_text(colour="black"),
          axis.text.x=element_text(colour="black")) +
    labs(
      x = "Enrichment",
      y = "",
      fill = "P-value adjusted",
      #size = "Enrichment"
    )
  
  return(p1)
}

plot_venn <- function(x, title) {
  p1 = ggVennDiagram(x) +
    #scale_color_brewer(palette = "Paired") +
    labs(title = title)
  
  return(p1)
}

plot_venn_v2 <- function(x, title) {
  VennDiag <- euler(x, shape = "ellipse")
  p1 <- plot(VennDiag, 
             counts = TRUE, 
             font=1, 
             cex=1, 
             alpha=0.3,
             main=title,
             fill=c("red", "green"),
             quantities = TRUE)
  return(p1)
}

