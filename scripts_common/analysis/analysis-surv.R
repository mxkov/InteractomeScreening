replace.null <- function(x) {
  if(length(x) == 0) {
    return(NA)
  }
  if(is.null(x) | is.na(x) | is.nan(x) | is.infinite(x)) {
    return(NA)
  } else {
    return(x)
  }
}


#
#
### SURVIVAL ANALYSIS
### (PATHWAYS)
#
#
analysis.surv <- function(PAL.matrix, pathways, clin.data, do.pfs,
                          doubt.patients, output.file,
                          debug = F, mode = c('median', 'quantiles'),
                          quant.probs = c(0.25, 0.75),
                          days.cutoff = NULL) {
  library('survival')
  library('survminer')
  library('survcomp')
  # os = overall survival, pfs = progression-free survival
  if(do.pfs) {
    surv.types <- c('os', 'pfs')
  } else {
    surv.types <- c('os')
  }
  days.col <- list('os' = 'days_os', 'pfs' = 'days_pfs')
  event.col <- list('os' = 'event_os', 'pfs' = 'event_pfs')
  for(st in surv.types) {
    if(!is.null(days.cutoff[[st]])) {
      message('Censoring events after ', days.col[[st]], ' = ',
              days.cutoff[[st]])
    }
  }
  # Initializing survival metrics
  pvals <- list()
  pvals.adj <- list()
  hr <- list()
  hr.lower <- list()
  hr.upper <- list()
  hr.pvals <- list()
  hr.pvals.adj <- list()
  for(st in surv.types) {
    pvals[[st]] <- numeric()
    pvals.adj[[st]] <- numeric()
    hr[[st]] <- numeric()
    hr.lower[[st]] <- numeric()
    hr.upper[[st]] <- numeric()
    hr.pvals[[st]] <- numeric()
    hr.pvals.adj[[st]] <- numeric()
  }
  #
  # Looping over all pathways
  #pathways <- rownames(PAL.matrix)
  pathways <- pathways[pathways %in% rownames(PAL.matrix)]
  for(pw in pathways) {
    pal.data <- PAL.matrix[pw,]
    pal.data <- pal.data[clin.data$case_submitter_id]
    clin.pal.df <- clin.data
    clin.pal.df$pal <- pal.data
    #clin.pal.df <- subset(clin.pal.df, pal != 0.0) # <--- NOPE
    # but we still need to check it!!!
    #clin.pal.df <- subset(clin.pal.df, !is.infinite(pal))  # no inf PALs
    #
    for(st in surv.types) {
      # Cleaning up the data
      clin.pal.df.st <- clin.pal.df[, c('case_submitter_id', 'pal',
                                        days.col[[st]],
                                        event.col[[st]])]
      cond1 <- !is.na(clin.pal.df.st[, days.col[[st]]])
      cond2 <- !is.na(clin.pal.df.st[, event.col[[st]]])
      cond3 <- clin.pal.df.st[, days.col[[st]]] >= 0
      cond4 <- !(clin.pal.df.st[, 'case_submitter_id'] %in% doubt.patients)
      clin.pal.df.st <- clin.pal.df.st[cond1 & cond2 & cond3 & cond4, ]
      #
      # that condition needs to be checked NOW, after the filtering:
      if(all(clin.pal.df.st$pal == 0.0)) {
        pvals[[st]] <- c(pvals[[st]], NA)
        hr[[st]] <- c(hr[[st]], NA)
        hr.lower[[st]] <- c(hr.lower[[st]], NA)
        hr.upper[[st]] <- c(hr.upper[[st]], NA)
        hr.pvals[[st]] <- c(hr.pvals[[st]], NA)
        next
      }
      # Censoring events after cutoff date
      if(!is.null(days.cutoff[[st]])) {
        filter <- (clin.pal.df.st[,days.col[[st]]] > days.cutoff[[st]])
        clin.pal.df.st[filter, event.col[[st]]] <- 0
        clin.pal.df.st[filter, days.col[[st]]] <- days.cutoff[[st]]
      }
      # Splitting the cases into groups
      if(mode == 'median') {
        value <- (clin.pal.df.st$pal >= median(clin.pal.df.st$pal))
        clin.pal.df.st$high_pal <- value
      }
      else if(mode == 'quantiles') {
        quant.probs <- quant.probs[order(quant.probs)]
        qs <- quantile(clin.pal.df.st$pal, probs = quant.probs)
        filter1 <- (clin.pal.df.st$pal <= qs[1])
        filter2 <- (clin.pal.df.st$pal >= qs[2])
        clin.pal.df.st <- clin.pal.df.st[filter1 | filter2,]
        clin.pal.df.st$high_pal <- (clin.pal.df.st$pal >= qs[2])
      }
      # Survival analysis
      surv.obj <- Surv(clin.pal.df.st[,days.col[[st]]],
                       clin.pal.df.st[,event.col[[st]]])
      surv.fit <- surv_fit(surv.obj ~ clin.pal.df.st$high_pal,
                           data = clin.pal.df.st)
      surv.fit.pval <- surv_pvalue(surv.fit)$pval
      hratio <- hazard.ratio(x = clin.pal.df.st$high_pal,
                             surv.time = clin.pal.df.st[,days.col[[st]]],
                             surv.event = clin.pal.df.st[,event.col[[st]]],
                             method.test = 'wald', na.rm = FALSE)
      # Record this
      pvals[[st]] <- c(pvals[[st]], replace.null(surv.fit.pval))
      hr[[st]] <- c(hr[[st]], replace.null(hratio$hazard.ratio))
      hr.lower[[st]] <- c(hr.lower[[st]], replace.null(hratio$lower))
      hr.upper[[st]] <- c(hr.upper[[st]], replace.null(hratio$upper))
      hr.pvals[[st]] <- c(hr.pvals[[st]], replace.null(hratio$p.value))
      
    }
  }
  
  # FDR-adjusting all p-values
  for(st in surv.types) {
    pvals.adj[[st]] <- p.adjust(pvals[[st]], method = 'fdr')
    hr.pvals.adj[[st]] <- p.adjust(hr.pvals[[st]], method = 'fdr')
  }
  # Putting the results together
  if(debug) {
    print('LENGTHS:')
    print(length(pathways))
    for(st in surv.types) {
      print(paste0('____', st))
      print(length(pvals[[st]]))
      print(length(pvals.adj[[st]]))
      print(length(hr[[st]]))
      print(length(hr.lower[[st]]))
      print(length(hr.upper[[st]]))
      print(length(hr.pvals[[st]]))
      print(length(hr.pvals.adj[[st]]))
    }
  }
  stat.results <- data.frame(pathway = pathways)
  for(st in surv.types) {
    stat.results[[paste0(st, '_pval')]] <- pvals[[st]]
    stat.results[[paste0(st, '_pval.adj')]] <- pvals.adj[[st]]
    stat.results[[paste0(st, '_hr')]] <- hr[[st]]
    stat.results[[paste0(st, '_hr.ci.lower')]] <- hr.lower[[st]]
    stat.results[[paste0(st, '_hr.ci.upper')]] <- hr.upper[[st]]
    stat.results[[paste0(st, '_hr.pval')]] <- hr.pvals[[st]]
    stat.results[[paste0(st, '_hr.pval.adj')]] <- hr.pvals.adj[[st]]
  }
  write.csv(stat.results, file = output.file,
            quote = FALSE, row.names = FALSE)
}


#
#
### SURVIVAL ANALYSIS
### (GENES)
#
#
analysis.surv.genes <- function(expr.matrix, genes, clin.data, surv.types,
                                doubt.patients, output.file,
                                debug = F, mode = c('median', 'quantiles'),
                                quant.probs = c(0.25, 0.75),
                                days.cutoff = NULL) {
  library('survival')
  library('survminer')
  library('survcomp')
  small.matrix <- nrow(expr.matrix) == 1
  # os = overall survival, pfs = progression-free survival
  days.col <- list('os' = 'days_os', 'pfs' = 'days_pfs')
  event.col <- list('os' = 'event_os', 'pfs' = 'event_pfs')
  for(st in surv.types) {
    if(!is.null(days.cutoff[[st]])) {
      message('Censoring events after ', days.col[[st]], ' = ',
              days.cutoff[[st]])
    }
  }
  # Initializing survival metrics
  pvals <- list()
  pvals.adj <- list()
  hr <- list()
  hr.lower <- list()
  hr.upper <- list()
  hr.pvals <- list()
  hr.pvals.adj <- list()
  for(st in surv.types) {
    pvals[[st]] <- numeric()
    pvals.adj[[st]] <- numeric()
    hr[[st]] <- numeric()
    hr.lower[[st]] <- numeric()
    hr.upper[[st]] <- numeric()
    hr.pvals[[st]] <- numeric()
    hr.pvals.adj[[st]] <- numeric()
  }
  #
  # Looping over all genes
  if(!small.matrix) {
    genes <- genes[genes %in% rownames(expr.matrix)]
  }
  for(g in genes) {
    if(small.matrix) {
      expr.data <- expr.matrix[1,]
    } else {
      expr.data <- expr.matrix[g,]
    }
    expr.data <- expr.data[clin.data$case_submitter_id]
    clin.expr.df <- clin.data
    clin.expr.df$expr <- expr.data
    #clin.expr.df <- subset(clin.expr.df, expr != 0.0) # <--- NOPE
    # but we still need to check it!!!
    #
    for(st in surv.types) {
      clin.expr.df.st <- clin.expr.df[, c('case_submitter_id', 'expr',
                                        days.col[[st]],
                                        event.col[[st]])]
      cond1 <- !is.na(clin.expr.df.st[, days.col[[st]]])
      cond2 <- !is.na(clin.expr.df.st[, event.col[[st]]])
      cond3 <- clin.expr.df.st[, days.col[[st]]] >= 0
      cond4 <- !(clin.expr.df.st[, 'case_submitter_id'] %in% doubt.patients)
      clin.expr.df.st <- clin.expr.df.st[cond1 & cond2 & cond3 & cond4, ]
      #
      # that condition needs to be checked NOW, after the filtering:
      if(all(clin.expr.df.st$expr == 0.0)) {
        pvals[[st]] <- c(pvals[[st]], NA)
        hr[[st]] <- c(hr[[st]], NA)
        hr.lower[[st]] <- c(hr.lower[[st]], NA)
        hr.upper[[st]] <- c(hr.upper[[st]], NA)
        hr.pvals[[st]] <- c(hr.pvals[[st]], NA)
        next
      }
      # Censoring events after cutoff date
      if(!is.null(days.cutoff[[st]])) {
        filter <- (clin.expr.df.st[,days.col[[st]]] > days.cutoff[[st]])
        clin.expr.df.st[filter, event.col[[st]]] <- 0
        clin.expr.df.st[filter, days.col[[st]]] <- days.cutoff[[st]]
      }
      # Splitting the cases into groups
      if(mode == 'median') {
        value <- (clin.expr.df.st$expr >= median(clin.expr.df.st$expr))
        clin.expr.df.st$high_expr <- value
      } else if(mode == 'quantiles') {
        quant.probs <- quant.probs[order(quant.probs)]
        qs <- quantile(clin.expr.df.st$expr, probs = quant.probs)
        filter1 <- (clin.expr.df.st$expr <= qs[1])
        filter2 <- (clin.expr.df.st$expr >= qs[2])
        clin.expr.df.st <- clin.expr.df.st[filter1 | filter2,]
        clin.expr.df.st$high_expr <- (clin.expr.df.st$expr >= qs[2])
      }
      # Survival analysis
      surv.obj <- Surv(clin.expr.df.st[,days.col[[st]]],
                       clin.expr.df.st[,event.col[[st]]])
      surv.fit <- surv_fit(surv.obj ~ clin.expr.df.st$high_expr,
                           data = clin.expr.df.st)
      surv.fit.pval <- surv_pvalue(surv.fit)$pval
      hratio <- hazard.ratio(x = clin.expr.df.st$high_expr,
                             surv.time = clin.expr.df.st[,days.col[[st]]],
                             surv.event = clin.expr.df.st[,event.col[[st]]],
                             method.test = 'wald', na.rm = FALSE)
      # Record this
      pvals[[st]] <- c(pvals[[st]], replace.null(surv.fit.pval))
      hr[[st]] <- c(hr[[st]], replace.null(hratio$hazard.ratio))
      hr.lower[[st]] <- c(hr.lower[[st]], replace.null(hratio$lower))
      hr.upper[[st]] <- c(hr.upper[[st]], replace.null(hratio$upper))
      hr.pvals[[st]] <- c(hr.pvals[[st]], replace.null(hratio$p.value))
      
    }
  }
  
  # FDR-adjusting all p-values
  for(st in surv.types) {
    pvals.adj[[st]] <- p.adjust(pvals[[st]], method = 'fdr')
    hr.pvals.adj[[st]] <- p.adjust(hr.pvals[[st]], method = 'fdr')
  }
  # Putting the results together
  if(debug) {
    print('LENGTHS:')
    print(length(genes))
    for(st in surv.types) {
      print(paste0('____', st))
      print(length(pvals[[st]]))
      print(length(pvals.adj[[st]]))
      print(length(hr[[st]]))
      print(length(hr.lower[[st]]))
      print(length(hr.upper[[st]]))
      print(length(hr.pvals[[st]]))
      print(length(hr.pvals.adj[[st]]))
    }
  }
  stat.results <- data.frame(gene = genes)
  for(st in surv.types) {
    stat.results[[paste0(st, '_pval')]] <- pvals[[st]]
    stat.results[[paste0(st, '_pval.adj')]] <- pvals.adj[[st]]
    stat.results[[paste0(st, '_hr')]] <- hr[[st]]
    stat.results[[paste0(st, '_hr.ci.lower')]] <- hr.lower[[st]]
    stat.results[[paste0(st, '_hr.ci.upper')]] <- hr.upper[[st]]
    stat.results[[paste0(st, '_hr.pval')]] <- hr.pvals[[st]]
    stat.results[[paste0(st, '_hr.pval.adj')]] <- hr.pvals.adj[[st]]
  }
  write.csv(stat.results, file = output.file,
            quote = FALSE, row.names = FALSE)
}
