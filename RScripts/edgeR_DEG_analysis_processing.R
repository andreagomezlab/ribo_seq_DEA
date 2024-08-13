

filter_genes <- function(y, min_count=10, large_n=10, min_prop=0.7) {
  keep <- filterByExpr(y, min.count=min_count, large.n=large_n, min.prop=min_prop)
  print(table(keep))
  y <- y[keep, , keep.lib.sizes=FALSE]
  return(y)
}

get_DEG <- function(y) {
  design <- model.matrix(~meta_table$Sex+meta_table$Treatment)
  rownames(design) <- colnames(y)
  y <- estimateDisp(y, design)
  et <- exactTest(y)
  topTags(et)
  data = et$table
  colnames(data) <- c("logFC", "logCPM", "PValue")
  data$Genes = rownames(data)
  return(data)
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