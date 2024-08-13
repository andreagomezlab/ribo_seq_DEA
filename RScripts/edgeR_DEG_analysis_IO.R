library(edgeR)


get_metatable <- function(cell_type='PV', timepoint='48hrs') {
  meta_table <- read.delim("files/84_samples_meta_data.csv", sep = ",")
  meta_table$Sample.ID = gsub("\\-", "", meta_table$Sample.ID)
  rownames(meta_table) <- meta_table$Sample.ID
  meta_table = meta_table[order(meta_table$CellType, meta_table$Treatment, meta_table$Timepoint, decreasing = T),]
  meta_table$Treatment[meta_table$Treatment == "Saline"] = "Ctrl"
  meta_table$Batch = paste0("B", meta_table$Batch)
  meta_table$Group = factor(paste(meta_table$Treatment,meta_table$Timepoint,sep="."))
  
  meta_table = meta_table[meta_table$CellType == cell_type & meta_table$Timepoint == timepoint,]
  meta_table = droplevels(meta_table)
  return(meta_table)
}

get_DGEList_counts <- function(meta_table) {
  Count_data <- read.delim("data/84_samples_merged_RSEM_bowtie2.genes.counts.matrix")
  colnames(Count_data) <- gsub("\\_.*", "", colnames(Count_data))
  colnames(Count_data) <- gsub("\\..*", "", colnames(Count_data))
  
  Counts_only <- Count_data[, -1] 
  rownames(Counts_only) <- Count_data$X
  Counts_only = Counts_only[,colnames(Counts_only) %in% meta_table$Sample.ID]
  Counts_only = Counts_only[, meta_table$Sample.ID]
  
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


