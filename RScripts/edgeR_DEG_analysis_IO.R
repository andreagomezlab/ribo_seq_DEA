library(edgeR)
library(readxl)
library(svglite)
library(sjPlot)
library(Seurat)


get_genes_results_files <- function(meta_table, treat="DOI") {
  if (treat == "DOI") {
    dir_ = "/media/data01/processed/private_data/RNA-Seq/H202SC23053594/counts/"
    
  } else if (treat == "Psilocybin") {
    dir_ = "/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/counts/"
  }
  files = paste0(dir_, meta_table$File.ID, ".genes.results")
  names(files) <- paste0(meta_table$Sample.ID)
  return(files)
}


get_cluster_genes <- function() {
  
  clu_data = read.delim(paste0("data/de_la_fuente_revenga_RNA_clusters/cluster_1_RNA_noorder.csv"), sep = ",")
  df = clu_data[, c("X", "cluster")]
  df$cluster = paste0("Cluster1")
  
  for (i in 2:7){
    clu_data = read.delim(paste0("data/de_la_fuente_revenga_RNA_clusters/cluster_", i,"_RNA_noorder.csv"), sep = ",")
    aux = clu_data[, c("X", "cluster")]
    aux$cluster = paste0("Cluster", i)
    df = rbind(df, aux)
  }
  
  return(df)
}


get_scRNA_addition_dataset <- function() {
  load("/media/data01/processed/public_data/scRNA-Seq/GSE124952_mPFC_cocaine/all_cell_inDD_sobj.RData")
  sc_dataset <- UpdateSeuratObject(object = all_cell_inDD_sobj)
  head(sc_dataset@meta.data)
  print(table(sc_dataset$treatment))
  sc_dataset_ = subset(sc_dataset, subset = treatment == "Saline")
  return(sc_dataset_)
}


get_V1_dataset <- function() {
  neurons = readRDS("../FDLM/data/GSE102817_mouse_vc/seurat.rds")
  head(neurons@meta.data)
  #neurons = subset(neurons, subset = stim == "0h")
  table(neurons@meta.data$sample)
  return(neurons)
}

get_metatable <- function(cell_type='PV', timepoint='48hrs', treat="DOI", all=F) {
  if (treat == "DOI") {
    meta_table <- read.delim("files/90_samples_meta_data_DOI_fileID.csv", sep = ",")
  }else if (treat == "Psilocybin"){
    meta_table <- read.delim("files/28_samples_meta_data_Psilo_fileID.csv", sep = ",")
  }
  
  meta_table$Sample.ID = gsub("\\-", "", meta_table$Sample.ID)
  rownames(meta_table) <- meta_table$Sample.ID
  meta_table = meta_table[order(meta_table$CellType, meta_table$Treatment, meta_table$Timepoint, decreasing = T),]
  meta_table$Treatment[meta_table$Treatment == "Saline"] = "Ctrl"
  meta_table$Batch = paste0("B", meta_table$Batch)
  meta_table$Group = factor(paste(meta_table$Treatment,meta_table$Timepoint,sep="."))
  
  if (all) {
    return(meta_table)
  } else {
    meta_table = meta_table[meta_table$CellType %in% cell_type & meta_table$Timepoint %in% timepoint,]
    meta_table = droplevels(meta_table)
  }
  
  return(meta_table)
}

get_counts_RNA <- function(meta_table, treat="DOI", type="counts") {
  if (treat == "DOI") {
    if (type == "counts") {
      Count_data <- read.delim("data/90_samples_merged_RSEM_bowtie2.genes.counts.matrix")
    } else if (type == "FPKM") {
      Count_data <- read.delim("data/90_samples_merged_RSEM_bowtie2.genes.FPKM.matrix")
    } else if (type == "TPM") {
      Count_data <- read.delim("data/90_samples_merged_RSEM_bowtie2.genes.TPM.matrix")
    }
    colnames(Count_data) <- gsub("\\_.*", "", colnames(Count_data))
    colnames(Count_data) <- gsub("\\..*", "", colnames(Count_data))
    rownames(Count_data) <- Count_data$X
    Count_data <- Count_data[, -1] 
  } else if (treat == "Psilocybin"){
    if (type == "counts") {
      Count_data <- read.delim("data/28_samples_psilo_RSEM_bowtie2.genes.counts.matrix")
    } else if (type == "fpkm") {
      Count_data <- read.delim("data/28_samples_psilo_RSEM_bowtie2.genes.FPKM.matrix")
    } else if (type == "tpm") {
      Count_data <- read.delim("data/28_samples_psilo_RSEM_bowtie2.genes.TPM.matrix")
    }
    colnames(Count_data) <- gsub("\\_.*", "", colnames(Count_data))
    Counts_only = Count_data
  }
  Count_data = Count_data[, colnames(Count_data) %in% meta_table$Sample.ID]
  Count_data = Count_data[, match(meta_table$Sample.ID, colnames(Count_data))]
  return(Count_data)
}

get_DGEList_counts <- function(meta_table, treat="DOI") {
  
  if (treat == "DOI") {
    Count_data <- read.delim("data/90_samples_merged_RSEM_bowtie2.genes.counts.matrix")
    colnames(Count_data) <- gsub("\\_.*", "", colnames(Count_data))
    colnames(Count_data) <- gsub("\\..*", "", colnames(Count_data))
    rownames(Count_data) <- Count_data$X
    Counts_only <- Count_data[, -1] 
  } else if (treat == "Psilocybin") {
    Count_data <- read.delim("data/28_samples_psilo_RSEM_bowtie2.genes.counts.matrix")
    colnames(Count_data) <- gsub("\\_.*", "", colnames(Count_data))
    rownames(Count_data) <- Count_data$X
    Counts_only <- Count_data[, -1] 
  }
  
  Counts_only = Counts_only[, colnames(Counts_only) %in% meta_table$Sample.ID]
  Counts_only = Counts_only[, match(meta_table$Sample.ID, colnames(Counts_only))]
  
  y <- DGEList(counts=Counts_only, genes = Count_data$X, group=meta_table$Group)
  return(y)
}

save_file <- function(obj, dir, file_name, row=F, sep="\t") {
  write.table(obj, 
              file = paste0(dir, file_name), 
              sep = sep, 
              row.names = row, 
              quote = FALSE)
}

export_pdf <- function(obj, dir, file_name, h, w) {
  print("Exporting pdf file")
  filename = paste0(dir, file_name, '.pdf')
  print(filename)
  pdf(file = filename, width = w, height = h)
  print(obj)
  dev.off()
}

export_svg <- function(obj, dir, file_name, h=8, w=11){
  print("Exporting svg file")
  filename = paste0(dir, file_name, '.svg')
  print(filename)
  svglite(filename, width = w, height = h)
  if (!is.null(obj))
    print(obj)
  dev.off()
}

get_ECM_genes <- function() {
  ECM_genes = read.delim("files/ECM GENE LIST.txt", sep = " ", header = T)
  head(ECM_genes)
  
  ECM_genes$Gene_name <- str_to_title(ECM_genes$Gene_name) 
  rownames(ECM_genes) <- ECM_genes$Gene_name
  return(ECM_genes)
}

get_Receptors_T_genes <- function() {
  receptors = read.delim("files/receptor_transporter_list.csv", sep = ",", header = T)
  return(receptors)
}

get_Hrvatin_cell_marker_genes <- function() {
  Hrvatin_cell_markers = read.delim("files/Hrvatin_markers.csv", header = F, sep = ",")
  colnames(Hrvatin_cell_markers) <- c("ID", "Gene_name", "Cell_type")
  return(Hrvatin_cell_markers)
}


get_splicing_factor_genes <- function(id="Symbol") {
    table_file = read.delim("files/41593_2019_465_MOESM8_ESM.csv", sep = ",")
    if (id == "Symbol")
      return(table_file$Gene.Symbol)
    return(table_file$Ensembl_ID)
}

get_RBP_genes <- function() {
   
  table_file <- read_excel("files/41586_2020_2077_MOESM3_ESM.xlsx")
  rbps = table_file$`RBP name`
  rbps = rbps[!is.na(rbps)]
  rbps = str_to_title(rbps)
  
  return(rbps)
}

get_RBP_annotations <- function(type=NULL) {
  
  table_file <- read_excel("files/41586_2020_2077_MOESM3_ESM.xlsx", skip = 1, )
  if (type == 'functions') {
    cols = c(1, 3:22)
  } else if (type == 'domains') {
    cols = c(1, 38:44)
  } else if (type == 'localization') {
    cols = c(1, 23:37)
  } else if (type == 'experiments') {
    cols = c(1, 45:52)
  }
  
  table_file_filt = table_file[, cols]
  colnames(table_file_filt)[1] <- 'gene'
  table_file_filt$gene = str_to_title(table_file_filt$gene)
  table_file_filt = as.data.frame(table_file_filt)
  rownames(table_file_filt) <- table_file_filt$gene
  colnames(table_file_filt) <- make.names(colnames(table_file_filt), unique=TRUE)
  table_file_filt[table_file_filt == 0] = 'No'
  table_file_filt[table_file_filt == 1] = 'Yes'
  return(table_file_filt[,-1])
  
}

read_file <- function(dir, h=TRUE, s="\t"){
  file = read.delim(dir, header = h, sep = s)
  return(file)
}

get_GO_gene_terms <- function(go_term="GO:0050767", ct="PV", tp="48hrs") {
  table_genes = read.delim("files/GO_neurogenesis_terms.csv", sep = ",", header = T)
  head(table_genes)
  table_genes_filt = table_genes[table_genes$go == go_term & table_genes$celltype == ct & table_genes$timepoint == tp,]
  genes = unique(unlist(strsplit(table_genes_filt$genes, "/")))
  return(genes)
}

get_DEG_genes <- function(loc_res) {
  DEG_union = NULL
  for (row in 1:nrow(loc_res)) {
    if (!is.na(loc_res$Location[row])) {
      data = read_file(loc_res$Location[row])
      data_ = data[data$Expression != 'Unchanged',]
      print(paste(loc_res$CellType[row], loc_res$Timepoint[row], loc_res$Treatment[row], sep = " : "))
      print(table(data_$Expression))
      data_$CellType = loc_res$CellType[row]
      data_$Timepoint = loc_res$Timepoint[row]
      data_$Treatment = loc_res$Treatment[row]
      DEG_union = rbind(DEG_union, data_[, c('CellType', 'Timepoint', 'Treatment', 'Genes', 'Expression')])
    }
  }
  return(DEG_union)
}

get_FC_values <- function(data_FC, loc_res) {
  for (row in 1:nrow(loc_res)) {
    if (!is.na(loc_res$Location[row])) {
      data = read_file(loc_res$Location[row])
      colname = paste(loc_res$CellType[row], loc_res$Timepoint[row], loc_res$Treatment[row], sep = "_")
      print(colname)
      data_FC[, colname] = data$log2FoldChange[match(rownames(data_FC), data$Genes)]
    }
  }
  return(data_FC)
}


get_log2FC_values <- function(data_FC, loc_res) {
  for (row in 1:nrow(loc_res)) {
    if (!is.na(loc_res$Location[row])) {
      data = read_file(loc_res$Location[row])
      colname = paste(loc_res$CellType[row], loc_res$Timepoint[row], loc_res$Treatment[row], sep = "_")
      print(colname)
      data_FC[, colname] = data$log2FoldChange[match(rownames(data_FC), data$Genes)]
    }
  }
  return(data_FC)
}
