#
#
### TUMOR TYPE DIFFERENTIATION


## Two types, using classical binary ROC
##
analysis.diff <- function(value.matrix, items, clin.data, output.file,
                          mode = c('pal', 'expr')) {
  # items = list of pathways or genes to be analyzed (the rest are discarded)
  # value = PAL or expression, respectively
  
  library('pROC')
  # Initializing the metrics
  Ws <- numeric()
  p.vals <- numeric()
  AUCs <- numeric()
  directions <- character()
  # Selecting items to be analyzed
  small.matrix <- nrow(value.matrix) == 1
  if(!small.matrix) {
    items <- items[items %in% rownames(value.matrix)]
  }
  #
  # Looping over all items
  for(i in items) {
    if(small.matrix) {
      values <- value.matrix[1,]
    } else {
      values <- value.matrix[i,]
    }
    values <- values[clin.data$case_submitter_id]
    clin.value.df <- clin.data
    clin.value.df$value <- values
    #clin.value.df <- subset(clin.value.df, value != 0.0) # <--- NOPE
    # but we still need to check it:
    if(all(clin.value.df$value == 0.0)) {
      Ws <- c(Ws, NA)
      p.vals <- c(p.vals, NA)
      AUCs <- c(AUCs, NA)
      directions <- c(directions, NA)
      next
    }
    ## Wilcoxon's test
    wtest <- wilcox.test(value ~ tumor_type, data = clin.value.df)
    Ws <- c(Ws, wtest$statistic)
    p.vals <- c(p.vals, wtest$p.value)
    ## ROC AUC
    # direction = 'auto' MOSTLY gets it right. But not always.
    direction.i <- '<'
    auc.i <- auc(roc(clin.value.df$tumor_type, clin.value.df$value,
                     direction = direction.i, quiet = TRUE))
    if(auc.i < 0.5) {
      auc.i <- 1.0-auc.i
      direction.i <- '>'
    }
    AUCs <- c(AUCs, auc.i)
    directions <- c(directions, direction.i)
  }
  # FDR-adjusting the p-values
  p.vals.adj <- p.adjust(p.vals, method = 'fdr')
  # Putting the results together, writing to file
  stat.results <- data.frame(item = items,
                             w = Ws, w.pval = p.vals, w.pval.adj = p.vals.adj,
                             auc = AUCs, direction = directions)
  if(mode == 'pal') {
    colnames(stat.results)[1] <- 'pathway'
  } else if(mode == 'expr') {
    colnames(stat.results)[1] <- 'gene'
  }
  write.csv(stat.results, file = output.file,
            quote = FALSE, row.names = FALSE)
}


## More than two types, using multiclass ROC
##
analysis.diff.multi <- function(value.matrix, items, clin.data, output.file,
                                mode = c('pal', 'expr')) {
  # items = list of pathways or genes to be analyzed (the rest are discarded)
  # value = PAL or expression, respectively
  
  library('pROC')
  # Initializing the metrics
  AUCs <- numeric()
  # Selecting items to be analyzed
  small.matrix <- nrow(value.matrix) == 1
  if(!small.matrix) {
    items <- items[items %in% rownames(value.matrix)]
  }
  #
  # Looping over all items
  for(i in items) {
    if(small.matrix) {
      values <- value.matrix[1,]
    } else {
      values <- value.matrix[i,]
    }
    values <- values[clin.data$case_submitter_id]
    clin.value.df <- clin.data
    clin.value.df$value <- values
    #clin.value.df <- subset(clin.value.df, value != 0.0) # <--- NOPE
    # but we still need to check it:
    if(all(clin.value.df$value == 0.0)) {
      AUCs <- c(AUCs, NA)
      next
    }
    ## ROC AUC
    auc.i <- multiclass.roc(clin.value.df$tumor_type, clin.value.df$value,
                            quiet = TRUE)$auc
    if(auc.i < 0.5) {
      auc.i <- 1.0-auc.i
    }
    AUCs <- c(AUCs, auc.i)
  }
  # Putting the results together, writing to file
  stat.results <- data.frame(item = items,
                             auc = AUCs)
  if(mode == 'pal') {
    colnames(stat.results)[1] <- 'pathway'
  } else if(mode == 'expr') {
    colnames(stat.results)[1] <- 'gene'
  }
  write.csv(stat.results, file = output.file,
            quote = FALSE, row.names = FALSE)
  
}
